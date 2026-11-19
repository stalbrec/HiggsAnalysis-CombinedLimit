// HS3 JSON export and import for CMSHistSum and CMSHistErrorPropagator (RooFit::JSONIO).
#include <algorithm>
#include "../interface/CMSHistSum.h"
#include "../interface/CMSHistFuncWrapper.h"
#include "../interface/ProcessNormalization.h"
#include "../interface/AsymPow.h"
#include <stdexcept>
#include <vector>
#include <ostream>
#include <memory>
#include <functional>
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "vectorized.h"
#include <charconv>
#include <string>
#include <set>
#include <map>
#include <sstream>

#include "RooFitHS3/RooJSONFactoryWSTool.h"
#include "RooFitHS3/JSONIO.h"
#include "../interface/CombineHS3BBLite.h"
#include "../interface/CombineHS3Channel.h"
#include <RooWorkspace.h>
#include "../interface/SimpleGaussianConstraint.h"
#include <RooAbsRealLValue.h>

using RooFit::Detail::JSONNode;

// Helper: try to pick a reasonable process name.
static std::string cms_hist_sum_process_name(const CMSHistSum &hs, unsigned ip)
{
  auto *coeff = hs.coeffAt(ip);   
  
  if (coeff && coeff->stringAttributes().count("combine.process")) {
    return coeff->getStringAttribute("combine.process");
  }
  
  if(coeff){
    return coeff->GetName();
  }
  
  return "proc_" + std::to_string(ip);
}

// Constraint names on the modifiers of a CMSHistSum channel (CombineHS3Channel.h, as the
// FastVerticalInterpHistPdf2 exporter does): set by tryExportHistFactory for the duration of a
// channel's export, nullptr otherwise (CMSHistErrorPropagator channels name none).
struct CMSHistSumConstraintCtx {
  RooJSONFactoryWSTool *tool;
  RooArgSet observables;
  RooAbsArg const *channel;
  std::map<std::string, std::string> cache; // parameter name -> constraint name ("" = none)
};
static CMSHistSumConstraintCtx *g_hsCtx = nullptr;

/// The constraint pdf of a parameter: "<par>_Pdf" (ShapeTools.py), if it depends on the parameter and
/// on no observable; for the y_j of a type-1 bin (named <x>_y, CombineHS3BBLite.h) the one of x_j.
/// By name only: a structural search (as in VerticalInterpHistPdf.cc) would also find, in an imported
/// model, the type-1 constraints, whose mean y_j / E_j depends on every rate parameter.
static RooAbsPdf *findConstraintPdf(RooWorkspace &ws, RooAbsArg const &par, RooArgSet const &observables)
{
  std::string name = par.GetName();
  RooAbsArg const *dep = &par;
  auto *byName = ws.pdf((name + "_Pdf").c_str());
  if (!byName && name.size() > 2 && name.compare(name.size() - 2, 2, "_y") == 0)
    byName = ws.pdf((name.substr(0, name.size() - 2) + "_Pdf").c_str());
  if (byName && byName->dependsOn(*dep) && !byName->dependsOn(observables)) return byName;
  return nullptr;
}

/// Name `constraint` (if any) on modifier m, unless an enclosing product applies it already.
static void nameConstraint(JSONNode &m, std::string const &par, RooAbsPdf const *constraint)
{
  if (!g_hsCtx || !constraint) return;
  auto it = g_hsCtx->cache.find(par);
  if (it == g_hsCtx->cache.end()) {
    std::string name;
    if (!CombineHS3::constraintAppliedOutside(*g_hsCtx->channel, *constraint)) name = constraint->GetName();
    g_hsCtx->tool->queueExport(*constraint);
    it = g_hsCtx->cache.emplace(par, name).first;
  }
  if (!it->second.empty()) m["constraint"] << it->second;
}

static void nameConstraintOf(JSONNode &m, std::string const &par)
{
  if (!g_hsCtx) return;
  if (auto it = g_hsCtx->cache.find(par); it != g_hsCtx->cache.end()) {
    if (!it->second.empty()) m["constraint"] << it->second;
    return;
  }
  RooAbsArg const *arg = g_hsCtx->tool->workspace()->arg(par.c_str());
  RooAbsPdf const *c = arg ? findConstraintPdf(*g_hsCtx->tool->workspace(), *arg, g_hsCtx->observables) : nullptr;
  if (c) nameConstraint(m, par, c);
  else g_hsCtx->cache.emplace(par, "");
}

void exportOtherFactorModifiers(RooJSONFactoryWSTool* tool, JSONNode& mods, RooAbsReal* factor) {
  if (auto* rcv = dynamic_cast<RooConstVar*>(factor)) {
    auto& m = mods.append_child().set_map();
    m["name"] << rcv->GetName();
    m["parameter"] << rcv->getVal();
    m["type"] << "normfactor";
  } else if (auto* rrv = dynamic_cast<RooRealVar*>(factor)) {
    auto& m = mods.append_child().set_map();
    m["name"] << rrv->GetName();
    m["parameter"] << rrv->GetName();
    m["type"] << "normfactor";
    nameConstraintOf(m, rrv->GetName());
    tool->queueExport(*rrv);
  } else if (auto* prod = dynamic_cast<RooProduct*>(factor)) {
    for (auto* comp : prod->components()) {
      exportOtherFactorModifiers(tool, mods, static_cast<RooAbsReal*>(comp));
    }
  } else {
    // handling other (e.g. RooAddition)
    auto& m = mods.append_child().set_map();
    m["name"] << factor->GetName();
    m["type"] << "custom";
    tool->queueExport(*factor);
  }
}

// Shared by tryExportHistFactory (CMSHistSum) and tryExportHistFactoryFromPropagator
// (CMSHistErrorPropagator): a sample's coefficient is often a RooProduct combining a
// ProcessNormalization (the lnN / plain-rate part) with one AsymPow per `shape` systematic
// (ShapeTools' systeff_* rate term paired with that systematic's histosys). Recurse into
// RooProduct factors instead of requiring the whole coefficient to match one known type: each
// AsymPow is exactly the normsys it represents, so there is nothing to preserve by serializing
// the AsymPow node itself.
// DECISIONS.md 2026-11-11 item 4a: an HS3 histosys/normsys "parameter" is the name of a scalar free
// parameter (/spec/docs/chapters/2.1.1_fundamental_distributions.md, histfactory_dist). A coefficient
// that is not a plain RooRealVar (e.g. RooProduct(alpha, c)) has no exact representation there, so
// the export stops instead of writing the compound node's name. Mirrored in VerticalInterpHistPdf.cc,
// CMSHistSum.cc and JSONFactories_Combine.cxx.
static void requirePlainModifierParameter(RooAbsArg const &par, std::string const &pdfName,
                                          std::string const &sampleName, const char *modType)
{
   if (dynamic_cast<RooRealVar const *>(&par)) return;
   RooJSONFactoryWSTool::error(std::string("HS3 export of pdf '") + pdfName + "', sample '" + sampleName + "': " +
                               modType + " coefficient '" + par.GetName() + "' is a " + par.ClassName() +
                               ", not a plain RooRealVar; an HS3 modifier \"parameter\" must be a scalar free "
                               "parameter (see OPEN_ISSUES.md, scaled normsys/histosys coefficients).");
}

void exportCoeffModifier(RooJSONFactoryWSTool* tool, JSONNode& mods, RooAbsArg const* c,
                          std::string const& sampleName, std::string const& pdfName) {
  if (!c) {
    // do nothing
  } else if (auto* rcv = dynamic_cast<RooConstVar const*>(c)) {
    // the const'ness of the const parameter is handled elsewhere, we can just pretend like it's a normal NF here
    auto& m = mods.append_child().set_map();
    m["name"] << rcv->GetName();
    m["parameter"] << rcv->GetName();
    m["type"] << "normfactor";
  } else if (auto* rrv = dynamic_cast<RooRealVar const*>(c)) {
    auto& m = mods.append_child().set_map();
    m["name"] << rrv->GetName();
    m["parameter"] << rrv->GetName();
    m["type"] << "normfactor";
    nameConstraintOf(m, rrv->GetName());
    tool->queueExport(*rrv);
  } else if (auto* pn = dynamic_cast<const ProcessNormalization*>(c)) {
    for (auto const* t : pn->thetaList()) requirePlainModifierParameter(*t, pdfName, sampleName, "normsys");
    for (auto const* t : pn->asymmThetaList()) requirePlainModifierParameter(*t, pdfName, sampleName, "normsys");
    double nom = pn->nominalValue();
    auto variations = pn->sigmaVariationsAll();
    for (const auto& v : variations) {
      auto& m = mods.append_child().set_map();
      m["name"] << std::get<0>(v);
      m["parameter"] << std::get<0>(v);
      m["type"] << "normsys";
      nameConstraintOf(m, std::get<0>(v));
      auto& data = m["data"].set_map();
      // sigmaVariationsAll() already gives the relative factors at -1/+1 sigma (the nominal is
      // its own normfactor below), which is what normsys lo/hi mean.
      data["lo"] << std::get<1>(v);
      data["hi"] << std::get<2>(v);
    }

    for (auto& of : pn->otherFactorList()) {
      exportOtherFactorModifiers(tool, mods, static_cast<RooAbsReal*>(of));
    }

    if (nom != 1.) {
      auto& m = mods.append_child().set_map();
      m["name"] << sampleName + "_norm_nominal";
      m["parameter"] << nom;
      m["type"] << "normfactor";
    }
  } else if (auto* ap = dynamic_cast<const AsymPow*>(c)) {
    requirePlainModifierParameter(ap->theta(), pdfName, sampleName, "normsys");
    auto& m = mods.append_child().set_map();
    m["name"] << ap->theta().GetName();
    m["parameter"] << ap->theta().GetName();
    m["type"] << "normsys";
    nameConstraintOf(m, ap->theta().GetName());
    auto& data = m["data"].set_map();
    data["lo"] << ap->kappaLow().getVal();
    data["hi"] << ap->kappaHigh().getVal();
    tool->queueExport(ap->theta());
  } else if (auto* prod = dynamic_cast<RooProduct const*>(c)) {
    for (auto* comp : const_cast<RooProduct*>(prod)->components()) {
      exportCoeffModifier(tool, mods, comp, sampleName, pdfName);
    }
  } else {
    auto& m = mods.append_child().set_map();
    m["name"] << c->GetName();
    m["type"] << "custom";
    tool->queueExport(*c);
  }
}

// The reparameterized autoMCStats term (CombineHS3BBLite.h): a sample with contents `scale` and one
// unconstrained shapefactor, n_j += scale_j * pars_j. Bins without a parameter cannot be written yet.
void exportAdditiveSample(RooJSONFactoryWSTool* tool, JSONNode& samples, std::string const& sampleName,
                          std::string const& modifierName, std::vector<RooAbsReal const*> const& pars,
                          std::vector<double> const& scale) {
  for (auto const* p : pars) {
    if (!p) RooJSONFactoryWSTool::error("shapefactor sample '" + sampleName + "': a bin has no parameter");
  }
  auto &s = samples.append_child().set_map();
  s["name"] << sampleName;
  RooJSONFactoryWSTool::exportArray(scale.size(), scale.data(), s["data"].set_map()["contents"]);
  auto &m = s["modifiers"].set_seq().append_child().set_map();
  m["name"] << modifierName;
  auto &pnode = m["parameters"].set_seq();
  for (auto const* p : pars) {
    pnode.append_child() << p->GetName();
    tool->queueExport(*p);
  }
  m["type"] << "shapefactor";
}

// CMSHistSum's per-process autoMCStats bin parameters (bin types 2, 3; DECISIONS.md 2026-11-04
// hig20013), by bin and process (nullptr where there is none).
//   type 3: n_j += c_p * e_pj * x_jp, x_jp ~ Gauss(x_In | x_jp, 1)
//   type 2: n_j += c_p * (g_jp - 1) * h_pj(theta), g_jp = v_jp * (1 / n_jp), v_jp ~ Poisson(n_jp | v_jp)
// Written as additional samples, each with the process's rate modifiers (c_p) and a shapefactor:
//   mcstat_<p>        contents e_pj in p's type-3 bins (0 elsewhere), shapefactor x_jp
//   mcstat_<p>_bin<j> contents h_pj and p's histosys restricted to bin j, normfactor 1/n_jp
//                     (a literal), shapefactor v_jp in every bin
// and p's own sample has contents and histosys 0 in its type-2 bins, so that bin carries
// c_p * h_pj * v_jp / n_jp = c_p * h_pj + c_p * (g_jp - 1) * h_pj. Bins of the shapefactors where the
// sample has contents 0 borrow a neighbouring parameter, as for the `mcstat` sample.
static std::vector<std::vector<RooAbsReal const *>> perProcessBinPars(const CMSHistSum *hs)
{
  auto const &bt = hs->binTypes();
  RooArgList const &pars = hs->binParList();
  std::vector<std::vector<RooAbsReal const *>> out(hs->nBins(), std::vector<RooAbsReal const *>(hs->nProcs(), nullptr));
  for (unsigned j = 0, r = 0; j < bt.size() && j < hs->nBins(); ++j) {
    for (unsigned i = 0; i < bt[j].size(); ++i) {
      if (bt[j][i] < 1 || bt[j][i] >= 4) continue;
      auto const *par = static_cast<RooAbsReal const *>(pars.at(r++));
      if (bt[j][0] >= 2) out[j].at(i) = par;
    }
  }
  return out;
}

static bool isPerProcessBin(const CMSHistSum *hs, unsigned j, unsigned ip, unsigned type)
{
  auto const &bt = hs->binTypes();
  return j < bt.size() && bt[j][0] >= 2 && bt[j].at(ip) == type;
}

/// v_jp and 1/n_jp of a type-2 bin parameter g_jp = RooProduct(v_jp, RooConstVar(1/n_jp)).
static std::pair<RooRealVar const *, RooConstVar const *> poissonBinPar(RooAbsReal const *g)
{
  auto const *prod = dynamic_cast<RooProduct const *>(g);
  RooRealVar const *v = nullptr;
  RooConstVar const *c = nullptr;
  if (prod) {
    RooArgList comps = const_cast<RooProduct *>(prod)->components();
    if (comps.getSize() == 2) {
      v = dynamic_cast<RooRealVar const *>(comps.at(0));
      c = dynamic_cast<RooConstVar const *>(comps.at(1));
    }
  }
  if (!v || !c)
    RooJSONFactoryWSTool::error(std::string("CMSHistSum: type-2 bin parameter '") + g->GetName() +
                                "' is not RooProduct(RooRealVar, RooConstVar)");
  return {v, c};
}

static void exportShapefactor(RooJSONFactoryWSTool *tool, JSONNode &mods, std::string const &name,
                              std::vector<RooAbsReal const *> pars, RooAbsPdf const *constraint = nullptr)
{
  RooAbsReal const *fill = nullptr;
  for (auto const *p : pars) if (p) { fill = p; break; }
  for (auto &p : pars) {
    if (p) fill = p; else p = fill;
  }
  auto &m = mods.append_child().set_map();
  m["name"] << name;
  auto &pnode = m["parameters"].set_seq();
  for (auto const *p : pars) {
    pnode.append_child() << p->GetName();
    tool->queueExport(*p);
  }
  m["type"] << "shapefactor";
  if (constraint) nameConstraint(m, fill->GetName(), constraint);
}

static void exportPerProcessBinSamples(RooJSONFactoryWSTool *tool, const CMSHistSum *hs, JSONNode &samples)
{
  auto const bp = perProcessBinPars(hs);
  const unsigned nBins = hs->nBins();
  RooWorkspace &ws = *tool->workspace();
  RooArgSet obs(hs->getXVar());
  for (unsigned ip = 0; ip < hs->nProcs(); ++ip) {
    const std::string proc = cms_hist_sum_process_name(*hs, ip);
    auto const errs = hs->binErrors(ip, false);
    auto const nomT = hs->nominalTemplate(ip, false);
    for (unsigned j = 0; j < nBins; ++j) {
      if (isPerProcessBin(hs, j, ip, 3)) {
        auto const *x = dynamic_cast<RooRealVar const *>(bp[j][ip]);
        if (!x)
          RooJSONFactoryWSTool::error(std::string("CMSHistSum '") + hs->GetName() + "': type-3 bin parameter '" +
                                      bp[j][ip]->GetName() + "' is not a RooRealVar");
        if (errs[j] == 0.)
          RooJSONFactoryWSTool::error(std::string("CMSHistSum '") + hs->GetName() + "': type-3 bin with error 0");
        const std::string name = "mcstat_" + proc + "_gaus_bin" + std::to_string(j);
        auto &s = samples.append_child().set_map();
        s["name"] << name;
        std::vector<double> contents(nBins, 0.);
        contents[j] = errs[j];
        RooJSONFactoryWSTool::exportArray(nBins, contents.data(), s["data"].set_map()["contents"]);
        auto &mods = s["modifiers"].set_seq();
        exportCoeffModifier(tool, mods, hs->coeffAt(ip), proc, hs->GetName());
        exportShapefactor(tool, mods, name, std::vector<RooAbsReal const *>(nBins, x), findConstraintPdf(ws, *x, obs));
      } else if (isPerProcessBin(hs, j, ip, 2)) {
        auto [v, c] = poissonBinPar(bp[j][ip]);
        const std::string name = "mcstat_" + proc + "_bin" + std::to_string(j);
        auto &s = samples.append_child().set_map();
        s["name"] << name;
        std::vector<double> one(nBins, 0.);
        one[j] = nomT[j];
        RooJSONFactoryWSTool::exportArray(nBins, one.data(), s["data"].set_map()["contents"]);
        auto &mods = s["modifiers"].set_seq();
        exportCoeffModifier(tool, mods, hs->coeffAt(ip), proc, hs->GetName());
        for (unsigned iv = 0; iv < hs->nMorphs(); ++iv) {
          if (!hs->isMorphUsed(ip, iv)) continue;
          RooAbsReal const *par = hs->morphPar(iv);
          if (!par) continue;
          requirePlainModifierParameter(*par, hs->GetName(), proc, "histosys");
          std::vector<double> lo(nBins, 0.), hi(nBins, 0.);
          lo[j] = hs->downTemplate(ip, iv, false)[j];
          hi[j] = hs->upTemplate(ip, iv, false)[j];
          auto &m = mods.append_child().set_map();
          m["name"] << par->GetName();
          m["parameter"] << par->GetName();
          m["type"] << "histosys";
          nameConstraintOf(m, par->GetName());
          auto &d = m["data"].set_map();
          RooJSONFactoryWSTool::exportArray(nBins, lo.data(), d["lo"].set_map()["contents"]);
          RooJSONFactoryWSTool::exportArray(nBins, hi.data(), d["hi"].set_map()["contents"]);
        }
        auto &nf = mods.append_child().set_map();
        nf["name"] << c->GetName();
        nf["parameter"] << c->getVal();
        nf["type"] << "normfactor";
        exportShapefactor(tool, mods, name, std::vector<RooAbsReal const *>(nBins, v), findConstraintPdf(ws, *v, obs));
      }
    }
  }
}

/// The type-1 term as one sample per bin j, mcstat_bin<j>: contents 1 in bin j (0 elsewhere) and a
/// shapefactor with y_j in every bin, naming x_j's constraint (gaussian_dist with mean y_j / E_j).
static void exportType1BinSamples(RooJSONFactoryWSTool *tool, const CMSHistSum *hs, JSONNode &samples,
                                  std::vector<RooAbsReal const *> const &ys, std::vector<RooAbsPdf *> const *cons)
{
  const unsigned nBins = hs->nBins();
  RooArgSet obs(hs->getXVar());
  for (unsigned j = 0; j < nBins && j < ys.size(); ++j) {
    if (!ys[j]) continue;
    RooAbsPdf const *constraint = cons ? (*cons)[j] : findConstraintPdf(*tool->workspace(), *ys[j], obs);
    const std::string name = "mcstat_bin" + std::to_string(j);
    auto &s = samples.append_child().set_map();
    s["name"] << name;
    std::vector<double> one(nBins, 0.);
    one[j] = 1.0;
    RooJSONFactoryWSTool::exportArray(nBins, one.data(), s["data"].set_map()["contents"]);
    auto &mods = s["modifiers"].set_seq();
    exportShapefactor(tool, mods, name, std::vector<RooAbsReal const *>(nBins, ys[j]), constraint);
  }
}

bool tryExportHistFactory(RooJSONFactoryWSTool* tool, const CMSHistSum* hs, JSONNode& elem, RooAbsArg const* channel) {
  CMSHistSumConstraintCtx ctx{tool, RooArgSet(hs->getXVar()), channel ? channel : hs, {}};
  struct CtxGuard {
    CMSHistSumConstraintCtx *prev = g_hsCtx;
    ~CtxGuard() { g_hsCtx = prev; }
  } guard;
  g_hsCtx = &ctx;

  // if we have external morphs injected, you probably want to error out for now.
  if (hs->hasExternalMorphs()) {  // optional accessor; otherwise remove and always allow
    RooJSONFactoryWSTool::error(
        "CMSHistSum '" + std::string(hs->GetName()) +
        "': exporting with external morphs is not implemented (cannot express as HistFactory modifiers).");
  }

    // Axes: HistFactory exporter uses RooRealVar axis metadata.
    // Use the independent variable name and its range.
    RooAbsReal const &xAbs = hs->getXVar(); 
    auto const *x = dynamic_cast<RooRealVar const *>(&xAbs);
    if (!x) {
      RooJSONFactoryWSTool::error("CMSHistSum '" + std::string(hs->GetName()) +
                                 "': x is not a RooRealVar; cannot export axis metadata.");
    }

    auto &axes = elem["axes"].set_seq();
    auto &ax = axes.append_child().set_map();
    std::string axname = x->GetName();
    RooJSONFactoryWSTool::testValidName(axname, false);
    ax["name"] << axname;
    ax["nbins"] << static_cast<int>(hs->nBins());
    // CMS_TH1obs_<channel> is always combine's own dummy per-bin index, decoupled from the real
    // physical binning (baked into the FastHisto template arrays instead) -- unit-width [0, nBins]
    // is its universal semantic. x->getMin()/getMax() is not reliable: it happens to already equal
    // [0, nBins] for most channels but not all (seen concretely: a channel whose CMSHistSum has 31
    // bins but whose "x" RooRealVar reports range [0, 1]), and disagreeing with nBins here produces
    // an internally-inconsistent axis that silently scales every sample's integral by nBins on
    // reimport.
    ax["min"] << 0.0;
    ax["max"] << static_cast<double>(hs->nBins());

    // Samples array
    auto &samples = elem["samples"].set_seq();

    const unsigned nProcs = hs->nProcs();
    const unsigned nMorphs = hs->nMorphs();
    const unsigned nBins  = hs->nBins();
    auto const *bb = CombineHS3::BBLiteExport::active();
    std::vector<RooRealVar *> const *bbYs = bb ? bb->binParameters(*hs) : nullptr;
    bool anyType1 = false;
    for (auto const &t : hs->binTypes()) anyType1 |= t[0] == 1;

    for (unsigned ip = 0; ip < nProcs; ++ip) {
      const std::string sampleName = cms_hist_sum_process_name(*hs, ip);

      auto &s = samples.append_child().set_map();
      s["name"] << sampleName;

      auto &dataNode = s["data"].set_map();

      auto const &nomT = hs->nominalTemplate(ip, false);
      auto const &errs = hs->binErrors(ip, false);

      // type-2 bins: in the mcstat_<p>_bin<j> sample (exportPerProcessBinSamples)
      std::vector<double> nomVals(nomT.Values().begin(), nomT.Values().begin() + nBins);
      for (unsigned j = 0; j < nBins; ++j) if (isPerProcessBin(hs, j, ip, 2)) nomVals[j] = 0.;
      RooJSONFactoryWSTool::exportArray(nBins, nomVals.data(), dataNode["contents"]);
      // If there are nonzero bin errors, export them and add staterror modifier (HistFactory style).
      bool anyErr = false;
      for (double e : errs.Values()) { if (e != 0.0) { anyErr = true; break; } }
      if (anyErr) {
        RooJSONFactoryWSTool::exportArray(nBins, errs.Values().data(), dataNode["errors"]);
      }

      // --- Modifiers
      auto &mods = s["modifiers"].set_seq();

    auto* coeff = hs->coeffAt(ip);
    exportCoeffModifier(tool, mods, coeff, sampleName, hs->GetName());

      // Vertical morphs -> histosys modifiers
      for (unsigned iv = 0; iv < nMorphs; ++iv) {
        if(!hs->isMorphUsed(ip, iv)) continue;

        RooAbsReal const *par = hs->morphPar(iv);
        if (!par) continue;
        requirePlainModifierParameter(*par, hs->GetName(), sampleName, "histosys");
        const std::string sysName = par->GetName();

        auto const &upT = hs->upTemplate(ip, iv, false); 
        auto const &dnT = hs->downTemplate(ip, iv, false);
	  
        auto &m = mods.append_child().set_map();
        m["name"] << sysName;
        m["parameter"] << sysName;
        m["type"] << "histosys";
        nameConstraintOf(m, sysName);

        std::vector<double> lo(dnT.Values().begin(), dnT.Values().begin() + nBins);
        std::vector<double> hi(upT.Values().begin(), upT.Values().begin() + nBins);
        for (unsigned j = 0; j < nBins; ++j) if (isPerProcessBin(hs, j, ip, 2)) lo[j] = hi[j] = 0.;
        auto &d = m["data"].set_map();
        RooJSONFactoryWSTool::exportArray(nBins, lo.data(), d["lo"].set_map()["contents"]);
        RooJSONFactoryWSTool::exportArray(nBins, hi.data(), d["hi"].set_map()["contents"]);

        tool->queueExport(*par);
      }

      // Bare staterror marker (not exact) only for type-1 bin parameters outside a BBLiteExport,
      // which writes them as the `mcstat` sample: with autoMCStats -1 (no bin parameters) the model
      // has no statistical terms, and per-process bins are written by exportPerProcessBinSamples.
      if (anyErr && anyType1 && !bbYs) {
        auto &m = mods.append_child().set_map();
        m["name"] << "mcstat";
        m["type"] << "staterror";
      }
    }

    exportPerProcessBinSamples(tool, hs, samples);

    auto addPars = hs->additiveBinPars();
    if (bbYs && !addPars.empty())
      RooJSONFactoryWSTool::error("CMSHistSum '" + std::string(hs->GetName()) + "': type-1 and additive bin parameters");
    if (bbYs) {
      exportType1BinSamples(tool, hs, samples, std::vector<RooAbsReal const *>(bbYs->begin(), bbYs->end()),
                            bb->binConstraints(*hs));
    } else if (!addPars.empty() && hs->additiveSampleName().empty()) {
      // imported from mcstat_bin<j> samples: the same, with contents 1
      for (unsigned j = 0; j < nBins; ++j)
        if (addPars[j] && hs->additiveBinScale()[j] != 1.0)
          RooJSONFactoryWSTool::error("CMSHistSum '" + std::string(hs->GetName()) + "': additive bin scale != 1");
      exportType1BinSamples(tool, hs, samples, addPars, nullptr);
    } else if (!addPars.empty()) {
      exportAdditiveSample(tool, samples, hs->additiveSampleName(), hs->additiveModifierName(), addPars,
                           hs->additiveBinScale());
    }
    
    return true;
}

namespace {
// RAII: restores a RooRealVar to its original value when the guard goes out of scope -- used
// while numerically probing CMSHistFunc's shape at fixed vertical-morph parameter points.
struct ParamValueGuard {
  RooRealVar* var;
  double orig;
  explicit ParamValueGuard(RooRealVar* v) : var(v), orig(v ? v->getVal() : 0.0) {}
  ~ParamValueGuard() { if (var) var->setVal(orig); }
  ParamValueGuard(ParamValueGuard const&) = delete;
};

// Numerically read a CMSHistFunc's evaluated shape with its vertical-morph parameters forced to
// `thetas` (one entry per func->vmorphs(), in order; 0 = nominal). f(+-1) is the
// interpolation-scheme-independent anchor point HS3's histosys hi/lo already represents, so this
// works regardless of CMSHistFunc's own interpolation settings (vtype_/vsmooth_par_) -- the same
// convention tryExportHistFactory's own histosys export already relies on for CMSHistSum. Declines
// (false) if a vmorph parameter is not a RooRealVar (should not happen; defensive only).
bool probeCMSHistFuncShape(CMSHistFunc const* func, std::vector<double> const& thetas,
                            std::vector<double>& out) {
  RooArgList const& vm = func->vmorphs();
  if (static_cast<int>(thetas.size()) != vm.getSize()) return false;

  std::vector<std::unique_ptr<ParamValueGuard>> guards;
  guards.reserve(vm.getSize());
  for (int i = 0; i < vm.getSize(); ++i) {
    auto* rrv = dynamic_cast<RooRealVar*>(vm.at(i));
    if (!rrv) return false;
    guards.push_back(std::make_unique<ParamValueGuard>(rrv));
    rrv->setVal(thetas[i]);
  }
  func->updateCache();
  auto const& c = func->cache();
  auto const& vals = c.Values();
  out.assign(vals.data(), vals.data() + c.fullsize());
  return true;
  // guards restore original parameter values here, on scope exit
}

std::string cms_hist_func_process_name(const CMSHistFunc &f)
{
  if (f.stringAttributes().count("combine.process")) {
    return f.getStringAttribute("combine.process");
  }
  return f.GetName();
}
}  // namespace

// CMSHistErrorPropagator holds one CMSHistFunc per process (funcList()) plus a matching rate
// coefficient (coefList()) -- structurally the same "per-process shape + rate" content
// tryExportHistFactory reads from CMSHistSum's fused internal storage, just organized as separate
// objects instead of one combined one. Built whenever text2workspace.py runs its default
// autoMCStats path (i.e. without --use-histsum). The per-process CMSHistFuncWrapper objects
// ShapeTools.py optionally builds alongside it are auxiliary (plotting/shape-breakdown only, not
// part of the likelihood -- see --no-wrappers) and are not needed for this export.
bool tryExportHistFactoryFromPropagator(RooJSONFactoryWSTool* tool, const CMSHistErrorPropagator* prop, JSONNode& elem) {
  RooArgList const& funcs = prop->funcList();
  RooArgList const& coeffs = prop->coefList();
  const unsigned nProcs = funcs.getSize();
  if (nProcs == 0) {
    RooJSONFactoryWSTool::warning("CMSHistErrorPropagator '" + std::string(prop->GetName()) + "': no processes found!");
    return false;
  }

  // Decline pass: every process must be a plain vertical-only CMSHistFunc before a single byte is
  // written into elem, mirroring the contract documented in VerticalInterpHistPdf.cc.
  std::vector<CMSHistFunc const*> hfuncs(nProcs);
  for (unsigned ip = 0; ip < nProcs; ++ip) {
    hfuncs[ip] = dynamic_cast<CMSHistFunc const*>(funcs.at(ip));
    if (!hfuncs[ip]) {
      RooJSONFactoryWSTool::warning("CMSHistErrorPropagator '" + std::string(prop->GetName()) +
                                    "': func '" + funcs.at(ip)->GetName() +
                                    "' is not a CMSHistFunc; declining this export path.");
      return false;
    }
    if (hfuncs[ip]->hmorphs().getSize() != 0) {
      RooJSONFactoryWSTool::warning("CMSHistErrorPropagator '" + std::string(prop->GetName()) +
                                    "': process '" + hfuncs[ip]->GetName() +
                                    "' uses horizontal morphing, which this exporter cannot yet represent; declining.");
      return false;
    }
    for (auto* m : hfuncs[ip]->vmorphs()) {
      if (!dynamic_cast<RooRealVar const*>(m)) {
        RooJSONFactoryWSTool::warning("CMSHistErrorPropagator '" + std::string(prop->GetName()) +
                                      "': process '" + hfuncs[ip]->GetName() +
                                      "' has a vertical-morph parameter that is not a RooRealVar; declining.");
        return false;
      }
    }
  }

  RooAbsReal const &xAbs = hfuncs[0]->getXVar();
  auto const *x = dynamic_cast<RooRealVar const *>(&xAbs);
  if (!x) {
    RooJSONFactoryWSTool::error("CMSHistErrorPropagator '" + std::string(prop->GetName()) +
                                "': x is not a RooRealVar; cannot export axis metadata.");
  }

  auto &axes = elem["axes"].set_seq();
  auto &ax = axes.append_child().set_map();
  std::string axname = x->GetName();
  RooJSONFactoryWSTool::testValidName(axname, false);
  ax["name"] << axname;
  // min/max are written below, once nBins is known from probing the first sample's shape -- see
  // the comment at the equivalent spot in tryExportHistFactory for why x->getMin()/getMax() itself
  // is not a reliable source for this (CMS_TH1obs_<channel> is a dummy per-bin index; its own
  // declared range does not always already equal [0, nBins]).

  auto &samples = elem["samples"].set_seq();
  int nBins = 0;
  auto const *bb = CombineHS3::BBLiteExport::active();
  std::vector<RooRealVar *> const *bbYs = bb ? bb->binParameters(*prop) : nullptr;

  for (unsigned ip = 0; ip < nProcs; ++ip) {
    CMSHistFunc const *func = hfuncs[ip];
    const std::string sampleName = cms_hist_func_process_name(*func);
    RooArgList const &vm = func->vmorphs();
    const unsigned nMorphs = vm.getSize();

    std::vector<double> nomVals;
    if (!probeCMSHistFuncShape(func, std::vector<double>(nMorphs, 0.0), nomVals)) {
      RooJSONFactoryWSTool::error("CMSHistErrorPropagator '" + std::string(prop->GetName()) +
                                  "': failed to probe nominal shape for '" + sampleName + "'.");
    }

    const int proc_nBins = static_cast<int>(nomVals.size());
    if (nBins != 0 && proc_nBins != nBins) {
      RooJSONFactoryWSTool::error("CMSHistErrorPropagator '" + std::string(prop->GetName()) +
                                  "': samples have not the same number of bins!");
    }
    nBins = proc_nBins;

    auto &s = samples.append_child().set_map();
    s["name"] << sampleName;
    auto &dataNode = s["data"].set_map();
    RooJSONFactoryWSTool::exportArray(nBins, nomVals.data(), dataNode["contents"]);

    auto const &errs = func->errors();
    bool anyErr = false;
    for (double e : errs.Values()) { if (e != 0.0) { anyErr = true; break; } }
    if (anyErr) {
      RooJSONFactoryWSTool::exportArray(nBins, errs.Values().data(), dataNode["errors"]);
    }

    auto &mods = s["modifiers"].set_seq();

    // Rate part -- same coefficient dispatch CMSHistSum's own export uses.
    auto* coeff = (ip < static_cast<unsigned>(coeffs.getSize())) ? coeffs.at(ip) : nullptr;
    exportCoeffModifier(tool, mods, coeff, sampleName, prop->GetName());

    // Vertical morphs -> histosys modifiers, probed numerically at theta = +-1. Written in
    // CMSHistSum's canonical order (its std::map of morph names, i.e. sorted by name), so the
    // CMSHistSum the importer builds exports them in the same order.
    std::vector<unsigned> morphOrder(nMorphs);
    for (unsigned iv = 0; iv < nMorphs; ++iv) morphOrder[iv] = iv;
    std::sort(morphOrder.begin(), morphOrder.end(), [&vm](unsigned a, unsigned b) {
      return std::string(vm.at(a)->GetName()) < std::string(vm.at(b)->GetName());
    });
    for (unsigned iv : morphOrder) {
      RooAbsArg const *par = vm.at(iv);
      if (!par) continue;

      std::vector<double> hiThetas(nMorphs, 0.0), loThetas(nMorphs, 0.0);
      hiThetas[iv] = 1.0;
      loThetas[iv] = -1.0;
      std::vector<double> hiVals, loVals;
      if (!probeCMSHistFuncShape(func, hiThetas, hiVals) ||
          !probeCMSHistFuncShape(func, loThetas, loVals)) {
        RooJSONFactoryWSTool::error("CMSHistErrorPropagator '" + std::string(prop->GetName()) +
                                    "': failed to probe shape for '" + sampleName + "' / '" +
                                    par->GetName() + "'.");
      }

      const std::string sysName = par->GetName();
      auto &m = mods.append_child().set_map();
      m["name"] << sysName;
      m["parameter"] << sysName;
      m["type"] << "histosys";
      auto &d = m["data"].set_map();
      RooJSONFactoryWSTool::exportArray(nBins, loVals.data(), d["lo"].set_map()["contents"]);
      RooJSONFactoryWSTool::exportArray(nBins, hiVals.data(), d["hi"].set_map()["contents"]);

      tool->queueExport(*par);
    }

    // Staterror modifier (if errors were exported and there are bin parameters) -- same bare
    // HistFactory-style marker tryExportHistFactory uses. Not when the bin parameters are written
    // as the `mcstat` sample below (CombineHS3BBLite.h).
    if (anyErr && prop->hasBinPars() && !bbYs) {
      auto &m = mods.append_child().set_map();
      m["name"] << "mcstat";
      m["type"] << "staterror";
    }
  }

  if (bbYs) {
    // One parameter per bin of the axis: bins without a term (bin type 0, or beyond the
    // propagator's active bins) get contents 0 and borrow the nearest preceding y_j (the first one
    // for leading bins), so they add 0 * y_j = 0.
    std::vector<RooAbsReal const *> pars(nBins, nullptr);
    std::vector<double> scale(nBins, 0.0);
    for (unsigned j = 0; j < bbYs->size() && j < static_cast<unsigned>(nBins); ++j) {
      if ((*bbYs)[j]) { pars[j] = (*bbYs)[j]; scale[j] = 1.0; }
    }
    RooAbsReal const *fill = nullptr;
    for (auto const *p : pars) if (p) { fill = p; break; }
    for (auto &p : pars) {
      if (p) fill = p; else p = fill;
    }
    exportAdditiveSample(tool, samples, "mcstat", "mcstat", pars, scale);
  }

  ax["nbins"] << nBins;
  ax["min"] << 0.0;
  ax["max"] << static_cast<double>(nBins);
  return true;
}

#include "RooRealSumPdf.h"
#include "RooProdPdf.h"
#include "../interface/CombineHS3Channel.h"

const CMSHistSum* getCMSHistSum(const RooRealSumPdf* pdf){
  CMSHistSum* sum = 0;
  if(pdf->funcList().size() == 1 && pdf->coefList().size() == 1 && pdf->coefList().at(0)->GetName() == std::string("ONE") && (sum = dynamic_cast<CMSHistSum*>(pdf->funcList().at(0))) && sum)
    return sum;
  return nullptr;
}

const CMSHistSum* getCMSHistSum(const RooProdPdf* pdf){
  RooRealSumPdf* sum = 0;
  if(pdf->getComponents()->getSize() == 1 && (sum = dynamic_cast<RooRealSumPdf*>(pdf->getComponents()->first())) && sum){
    return getCMSHistSum(sum);
  }
  return nullptr;
}

// Same "RooRealSumPdf(func=[X], coef=[ONE])" wrapping pattern getCMSHistSum recognizes, but for
// the default autoMCStats path (text2workspace.py without --use-histsum), where X is a
// CMSHistErrorPropagator instead of a CMSHistSum.
const CMSHistErrorPropagator* getCMSHistErrorPropagator(const RooRealSumPdf* pdf){
  CMSHistErrorPropagator* prop = 0;
  if(pdf->funcList().size() == 1 && pdf->coefList().size() == 1 && pdf->coefList().at(0)->GetName() == std::string("ONE") && (prop = dynamic_cast<CMSHistErrorPropagator*>(pdf->funcList().at(0))) && prop)
    return prop;
  return nullptr;
}

const CMSHistErrorPropagator* getCMSHistErrorPropagator(const RooProdPdf* pdf){
  RooRealSumPdf* sum = 0;
  if(pdf->getComponents()->getSize() == 1 && (sum = dynamic_cast<RooRealSumPdf*>(pdf->getComponents()->first())) && sum){
    return getCMSHistErrorPropagator(sum);
  }
  return nullptr;
}

class CMSHistSumStreamer_CMSHistSum : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "histfactory";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
     auto sum = dynamic_cast<const CMSHistSum*>(p);
     if(sum){
       auto retval = tryExportHistFactory(tool,sum,elem,p);
       if(!retval) return false;
       elem["type"] << key();
        std::cout << "exporting CMSHistSum to histfactory" << std::endl;
       return true;
     }
     return false;
   }
};


class CMSHistSumStreamer_SumPdf : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "histfactory_dist";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
     auto* rsp = dynamic_cast<const RooRealSumPdf*>(p);
     if(auto sum = getCMSHistSum(rsp)){
       auto retval = tryExportHistFactory(tool,sum,elem,p);
       if(retval){
	 elem["type"] << key();
  std::cout << "exporting RooRealSumPdf to histfactory" << std::endl;

	 return true;
       }
       std::cout << "failure" << std::endl;
       return false;
     }
     if(auto prop = getCMSHistErrorPropagator(rsp)){
       auto retval = tryExportHistFactoryFromPropagator(tool,prop,elem);
       if(retval){
         elem["type"] << key();
         std::cout << "exporting RooRealSumPdf (CMSHistErrorPropagator) to histfactory" << std::endl;
         return true;
       }
       std::cout << "failure" << std::endl;
     }
     return false;
   }
};


class CMSHistSumStreamer_ProdPdf : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "histfactory_dist";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
     auto* pp = dynamic_cast<const RooProdPdf*>(p);
     if(auto sum = getCMSHistSum(pp)){
       auto retval = tryExportHistFactory(tool,sum,elem,p);
       if(!retval) return false;
       elem["type"] << key();
      std::cout << "exporting RooProdPdf to histfactory" << std::endl;
       return true;
     }
     if(auto prop = getCMSHistErrorPropagator(pp)){
       auto retval = tryExportHistFactoryFromPropagator(tool,prop,elem);
       if(!retval) return false;
       elem["type"] << key();
       std::cout << "exporting RooProdPdf (CMSHistErrorPropagator) to histfactory" << std::endl;
       return true;
     }
     return false;
   }
};

/// RooProdPdf(<name>, [<name>_model, constraints]) as built by Combine's histfactory_dist
/// importers: the histfactory_dist of the channel, whose modifiers name these constraints.
class CombineConstrainedChannelStreamer : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "histfactory_dist";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
      auto const *prod = dynamic_cast<const RooProdPdf *>(p);
      return prod && CombineHS3::exportConstrainedChannel(tool, *prod, elem);
   }
};

/// The constraint of a reparameterized autoMCStats bin parameter x_j (CombineHS3BBLite.h):
/// Gauss(x_In | x_j, 1) is written as gaussian_dist(x = x_In, mean = y_j / E_j, sigma = 1).
class BBLiteConstraintStreamer : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "gaussian_dist";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
      auto const *bb = CombineHS3::BBLiteExport::active();
      auto const *g = dynamic_cast<const SimpleGaussianConstraint *>(p);
      if (!bb || !g) return false;
      RooAbsReal *mean = bb->meanFor(*g);
      if (!mean) return false;
      auto const *glob = dynamic_cast<RooRealVar const *>(&g->getMean());
      if (!glob || !glob->isConstant() || !g->getSigma().isConstant() || g->getSigma().getVal() != 1.0) {
         RooJSONFactoryWSTool::error(std::string("autoMCStats constraint '") + g->GetName() +
                                     "' is not Gauss(x_In | x, 1) with a constant x_In");
      }
      elem["type"] << key();
      elem["x"] << glob->GetName();
      elem["mean"] << mean->GetName();
      elem["sigma"] << 1.0;
      tool->queueExport(*glob);
      tool->queueExport(*mean);
      return true;
   }
};

/// gaussian_dist(x = <variable>, mean = <function>, sigma): Combine's CachingSimNLL groups
/// SimpleGaussianConstraints (SIMNLL_GROUPCONSTRAINTS) and a group only tracks each member's x
/// (SimpleConstraintGroup::add), so a constraint whose parameters all sit in the mean (the
/// reparameterized autoMCStats constraint Gauss(x_In | y_j / E_j, 1), CombineHS3BBLite.h) would keep
/// a stale value during the fit. The Gaussian is symmetric in x and mean, so this importer builds it
/// as SimpleGaussianConstraint(x = <function>, mean = <variable>, sigma), the same value, and
/// SimpleGaussianConstraintSwappedStreamer writes it back unchanged. Every other gaussian_dist is left
/// to SimpleGaussianConstraint's factory expression.
class SimpleGaussianConstraintSwappedFactory : public RooFit::JSONIO::Importer {
public:
   bool importArg(RooJSONFactoryWSTool *tool, const JSONNode &p) const override
   {
      if (!p.has_child("x") || !p.has_child("mean") || !p.has_child("sigma")) return false;
      if (p["x"].is_number() || p["mean"].is_number()) return false;
      std::string const name(RooJSONFactoryWSTool::name(p));
      RooAbsReal *x = tool->requestArg<RooAbsReal>(p, "x");
      RooAbsReal *mean = tool->requestArg<RooAbsReal>(p, "mean");
      if (!dynamic_cast<RooAbsRealLValue *>(x) || dynamic_cast<RooAbsRealLValue *>(mean)) return false;
      RooAbsReal *sigma = p["sigma"].is_number() ? &RooFit::RooConst(p["sigma"].val_double())
                                                 : tool->requestArg<RooAbsReal>(p, "sigma");
      tool->wsEmplace<SimpleGaussianConstraint>(name, *mean, *x, *sigma);
      return true;
   }
};

class SimpleGaussianConstraintSwappedStreamer : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "gaussian_dist";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *arg, JSONNode &elem) const override
   {
      auto const *g = dynamic_cast<const SimpleGaussianConstraint *>(arg);
      if (!g || dynamic_cast<RooAbsRealLValue const *>(&g->getX()) ||
          !dynamic_cast<RooAbsRealLValue const *>(&g->getMean()))
         return false;
      elem["type"] << key();
      elem["x"] << g->getMean().GetName();
      elem["mean"] << g->getX().GetName();
      tool->queueExport(g->getMean());
      tool->queueExport(g->getX());
      if (auto const *c = dynamic_cast<RooConstVar const *>(&g->getSigma())) {
         elem["sigma"] << c->getVal();
      } else {
         elem["sigma"] << g->getSigma().GetName();
         tool->queueExport(g->getSigma());
      }
      return true;
   }
};

#include "static_execute.h"

STATIC_EXECUTE([]() {
  using namespace RooFit::JSONIO;
  
  registerExporter<BBLiteConstraintStreamer>(SimpleGaussianConstraint::Class(), true);
  registerExporter<SimpleGaussianConstraintSwappedStreamer>(SimpleGaussianConstraint::Class(), true);
  registerImporter<SimpleGaussianConstraintSwappedFactory>("gaussian_dist", true);
  registerExporter<CMSHistSumStreamer_SumPdf>(RooRealSumPdf::Class(), true);
  registerExporter<CMSHistSumStreamer_ProdPdf>(RooProdPdf::Class(), true);
  registerExporter<CombineConstrainedChannelStreamer>(RooProdPdf::Class(), true); // ahead of the one above
  registerExporter<CMSHistSumStreamer_CMSHistSum>("CMSHistSum", false);    
 });

namespace {

/// Helper: read an array of doubles
static std::vector<double> readDoubleArray(RooFit::Detail::JSONNode const &arr, std::string const &what)
{
   std::vector<double> out;
   out.reserve(arr.num_children());
   for(const auto& c:arr.children()){
     out.push_back(c.val_double());
   }
   return out;
}

FastHisto readFastHisto(const JSONNode &n,
					 const std::string &name,
					 RooArgSet const &vars)
{
  if (!n.has_child("contents"))
    RooJSONFactoryWSTool::error("no contents given");
  
  JSONNode const &contents = n["contents"];
  if (!contents.is_seq())
    RooJSONFactoryWSTool::error("contents are not in list form");
  
  JSONNode const *errors = nullptr;
  if (n.has_child("errors")) {
    errors = &n["errors"];
    if (!errors->is_seq())
      RooJSONFactoryWSTool::error("errors are not in list form");
  }

  if (vars.size() != 1) {
    std::stringstream errMsg;
    errMsg << "readFastHisto only supports 1D, got vars.size()=" << vars.size();
    RooJSONFactoryWSTool::error(errMsg.str());
  }

  auto *x = dynamic_cast<RooRealVar *>(vars.first());
  if (!x)
    RooJSONFactoryWSTool::error("readFastHisto expects a RooRealVar observable");
  
  const RooAbsBinning &binning = x->getBinning();
  const int nbins = binning.numBins();
  
  if (static_cast<int>(contents.num_children()) != nbins) {
    std::stringstream errMsg;
    errMsg << "inconsistent bin numbers: contents=" << contents.num_children() << ", nbins=" << nbins;
    RooJSONFactoryWSTool::error(errMsg.str());
  }
  if (errors && static_cast<int>(errors->num_children()) != nbins) {
    std::stringstream errMsg;
    errMsg << "inconsistent bin numbers: errors=" << errors->num_children() << ", nbins=" << nbins;
    RooJSONFactoryWSTool::error(errMsg.str());
  }
  
  // Read contents
  std::vector<double> vals;
  vals.reserve(contents.num_children());
  for (auto const &cont : contents.children()) {
    vals.push_back(cont.val_double());
  }
  
  // Build bin edges (nbins+1)
  std::vector<double> edges;
  edges.reserve(static_cast<std::size_t>(nbins) + 1u);
  for (int i = 0; i < nbins; ++i) edges.push_back(binning.binLow(i));
  edges.push_back(binning.binHigh(nbins - 1));
  
  // Construct FastHisto directly
  return FastHisto(edges, vals);
}
  

/// Helper: ensure/lookup observable RooRealVar by name
static RooRealVar &getObservable(RooJSONFactoryWSTool *tool, std::string const &obsName,
                                 double xmin, double xmax)
{
   // If it exists already, use it.
   if (auto *arg = tool->workspace()->arg(obsName.c_str())) {
      auto *rrv = dynamic_cast<RooRealVar *>(arg);
      if (!rrv)
         throw std::runtime_error("histfactory_dist importer: observable '" + obsName + "' exists but is not RooRealVar");
      return *rrv;
   }

   // Otherwise create it (range from axis definition)
   auto *rrv = new RooRealVar(obsName.c_str(), obsName.c_str(), 0.0, xmin, xmax);
   tool->workspace()->import(*rrv, RooFit::RecycleConflictNodes());
   return *dynamic_cast<RooRealVar *>(tool->workspace()->arg(obsName.c_str()));
}

/// Helper: ensure/lookup a parameter (RooRealVar) used by normfactor modifiers
static RooRealVar &getOrMakeParameter(RooJSONFactoryWSTool *tool, std::string const &parName,
                                     double init = 1.0, double lo = 0.0, double hi = 10.0)
{
   if (auto *arg = tool->workspace()->arg(parName.c_str())) {
      auto *rrv = dynamic_cast<RooRealVar *>(arg);
      if (!rrv)
         throw std::runtime_error("histfactory_dist importer: parameter '" + parName + "' exists but is not RooRealVar");
      return *rrv;
   }
   auto *rrv = new RooRealVar(parName.c_str(), parName.c_str(), init, lo, hi);
   tool->workspace()->import(*rrv, RooFit::RecycleConflictNodes());
   return *dynamic_cast<RooRealVar *>(tool->workspace()->arg(parName.c_str()));
}

bool parse_double(const std::string& input, double& out)
{
    const char* begin = input.data();
    const char* end   = input.data() + input.size();

    double tmp;
    auto [ptr, ec] = std::from_chars(begin, end, tmp);

    // Successful conversion requires:
    // 1) No error
    // 2) Entire string was consumed
    if((ec == std::errc() && ptr == end)){
      out = tmp;
      return true;
    }
    return false;
}

inline bool approxSymmetricKappas(double kLo, double kHi,
                                 double absLogTol = 1e-10,
                                 double relLogTol = 1e-10)
{
    // kappas must be positive and finite
    if (!(kLo > 0.0) || !(kHi > 0.0)) return false;
    if (!std::isfinite(kLo) || !std::isfinite(kHi)) return false;

    const double a = std::log(kLo);
    const double b = std::log(kHi);

    // symmetric if a ≈ -b  <=>  (a + b) ≈ 0
    const double diff  = std::abs(a + b);
    const double scale = std::max(std::abs(a), std::abs(b)); // 0 if both ~1
    return diff <= std::max(absLogTol, relLogTol * scale);
}

/// Dispatch predicate, duplicated from JSONFactories_Combine.cxx (and kept in sync with the copy
/// in VerticalInterpHistPdf.cc): true only for nodes that are unambiguously counting-shaped, i.e. a
/// single 1-bin axis, one content value per sample and no shape-like modifiers. Those belong to
/// the single-bin RooPoisson/RooAddition importer, not to this one.
inline bool isSingleBinCountingDist(RooFit::Detail::JSONNode const &elem)
{
    auto const *axes = elem.find("axes");
    if (!axes || axes->num_children() != 1) return false;
    for (auto const &ax : axes->children()) {
        if (ax.has_child("edges")) {
            if (ax["edges"].num_children() != 2) return false; // two edges == one bin
        } else {
            auto const *nb = ax.find("nbins");
            if (!nb || nb->val_int() != 1) return false;
            // the count axis n_obs_bin<b> is unbounded above; a 1-bin shape channel is not
            auto const *mx = ax.find("max");
            if (mx && !mx->is_null()) return false;
        }
    }

    auto const *samples = elem.find("samples");
    if (!samples || samples->num_children() == 0) return false;
    for (auto const &s : samples->children()) {
        auto const *data = s.find("data");
        if (!data || !data->has_child("contents") || (*data)["contents"].num_children() != 1) return false;
        if (auto const *mods = s.find("modifiers")) {
            for (auto const &m : mods->children()) {
                std::string const mtype = m["type"].val();
                if (mtype != "normsys" && mtype != "normfactor" && mtype != "custom") return false;
            }
        }
    }
    return true;
}

/// A FastHisto on x's binning with these bin contents.
static FastHisto fastHistoFrom(std::vector<double> const &vals, RooRealVar const &x, std::string const &what)
{
  const RooAbsBinning &binning = x.getBinning();
  const int nbins = binning.numBins();
  if (static_cast<int>(vals.size()) != nbins)
    RooJSONFactoryWSTool::error("inconsistent bin numbers in " + what);
  std::vector<double> edges;
  edges.reserve(static_cast<std::size_t>(nbins) + 1u);
  for (int i = 0; i < nbins; ++i) edges.push_back(binning.binLow(i));
  edges.push_back(binning.binHigh(nbins - 1));
  return FastHisto(edges, vals);
}

static std::string jsonDump(RooFit::Detail::JSONNode const &n)
{
  std::ostringstream os;
  n.writeJSON(os);
  return os.str();
}

static bool isRateModifier(RooFit::Detail::JSONNode const &m)
{
  std::string const t = m["type"].val();
  return t == "normfactor" || t == "normsys" || t == "custom";
}

/// The per-process autoMCStats samples of one process (exportPerProcessBinSamples).
struct PerProcessBinSamples {
  std::vector<std::pair<unsigned, RooFit::Detail::JSONNode const *>> type3; // (j, mcstat_<p>_gaus_bin<j>)
  std::vector<std::pair<unsigned, RooFit::Detail::JSONNode const *>> type2; // (j, mcstat_<p>_bin<j>)
};

/// "<stem>_bin<j>" -> (stem, j); false if `name` has no such suffix.
static bool splitBinSuffix(std::string const &name, std::string &stem, unsigned &j)
{
  auto k = name.rfind("_bin");
  if (k == std::string::npos || k + 4 >= name.size()) return false;
  std::string const digits = name.substr(k + 4);
  if (digits.find_first_not_of("0123456789") != std::string::npos) return false;
  stem = name.substr(0, k);
  j = std::stoul(digits);
  return true;
}

/// The actual importer
class CMSHistSumHistFactoryImporter final : public RooFit::JSONIO::Importer {
public:
   bool importArg(RooJSONFactoryWSTool *tool, RooFit::Detail::JSONNode const &elem) const override
   {
      using RooFit::Detail::JSONNode;

      // Basic validation
      if (!elem.has_child("name"))
         return false;

      // Single-bin counting nodes belong to the RooPoisson/RooAddition importer.
      if (isSingleBinCountingDist(elem))
         return false;

      const std::string distName = elem["name"].val();

      // Axes: we only support 1D here
      if (!elem.has_child("axes") || elem["axes"].num_children() != 1)
         throw std::runtime_error("CMSHistSum importer: '" + distName + "' must have exactly one axis");

      RooArgSet observables = RooJSONFactoryWSTool::readAxes(elem);
      RooRealVar* x = dynamic_cast<RooRealVar*>(observables.first());
      
      // Samples
      if (!elem.has_child("samples") || elem["samples"].num_children() == 0)
         throw std::runtime_error("CMSHistSum importer: '" + distName + "' has no samples");

      RooArgList funcs;
      RooArgList coeffs;

      // Ensure ONE exists (fixed at 1) for the wrapper RooRealSumPdf
      RooAbsReal *ONE = dynamic_cast<RooAbsReal *>(tool->workspace()->arg("ONE"));
      if (!ONE) {
         auto *one = new RooConstVar("ONE", "ONE", 1.0);
         tool->workspace()->import(*one, RooFit::RecycleConflictNodes());
         ONE = dynamic_cast<RooAbsReal *>(tool->workspace()->arg("ONE"));
      }      

      // An unconstrained shapefactor sample with no other modifiers: the reparameterized autoMCStats
      // term (CombineHS3BBLite.h), added per bin to the process sum (CMSHistSum::setAdditiveBinPars).
      struct Additive {
        std::string sample, modifier;
        std::vector<RooAbsReal *> pars;
        std::vector<double> scale;
      };
      std::unique_ptr<Additive> additive, type1;

      // Per-process autoMCStats samples (exportPerProcessBinSamples): a shapefactor together with other
      // modifiers, named mcstat_<p> (type 3) or mcstat_<p>_bin<j> (type 2) after a process p.
      std::map<std::string, PerProcessBinSamples> perProc;
      {
        std::set<std::string> procs;
        std::vector<RooFit::Detail::JSONNode const *> bbs;
        for (auto const &s : elem["samples"].children()) {
          bool sf = false;
          if (s.has_child("modifiers"))
            for (auto const &m : s["modifiers"].children()) sf |= m["type"].val() == "shapefactor";
          if (!sf) procs.insert(s["name"].val());
          else if (s["modifiers"].num_children() > 1) bbs.push_back(&s);
        }
        for (auto const *s : bbs) {
          std::string const name = (*s)["name"].val();
          std::string const prefix = "mcstat_";
          std::string stem;
          unsigned j = 0;
          if (name.compare(0, prefix.size(), prefix) != 0 || !splitBinSuffix(name.substr(prefix.size()), stem, j))
            throw std::runtime_error("CMSHistSum importer: shapefactor on sample '" + name +
                                     "' together with other modifiers is not supported");
          std::string const gaus = "_gaus";
          bool const t2 = procs.count(stem) > 0;
          bool const t3 = stem.size() > gaus.size() && stem.compare(stem.size() - gaus.size(), gaus.size(), gaus) == 0 &&
                          procs.count(stem.substr(0, stem.size() - gaus.size())) > 0;
          if (t3 == t2)
            throw std::runtime_error("CMSHistSum importer: cannot tell which process autoMCStats sample '" + name +
                                     "' belongs to");
          if (t3) perProc[stem.substr(0, stem.size() - gaus.size())].type3.emplace_back(j, s);
          else perProc[stem].type2.emplace_back(j, s);
        }
      }
      std::map<std::string, unsigned> procIndex;

      for(const auto&s:elem["samples"].children()){
	const std::string sampleName = s["name"].val();

         if (s.has_child("modifiers")) {
            bool hasShapefactor = false;
            for (auto const &m : s["modifiers"].children()) hasShapefactor |= m["type"].val() == "shapefactor";
            if (hasShapefactor) {
               auto const &mods = s["modifiers"];
               if (mods.num_children() != 1) continue; // per-process autoMCStats sample, below
               auto const &m = *mods.children().begin();
               // the type-1 term, one sample per bin (exportType1BinSamples); its `constraint` is
               // applied by wrapChannelInConstraints
               std::string stem;
               unsigned jb = 0;
               if (splitBinSuffix(sampleName, stem, jb) && stem == "mcstat") {
                  if (m["name"].val() != sampleName)
                     throw std::runtime_error("CMSHistSum importer: shapefactor of '" + sampleName + "' has another name");
                  std::set<std::string> pn;
                  for (auto const &par : m["parameters"].children()) pn.insert(par.val());
                  auto const vals = readDoubleArray(s["data"]["contents"], sampleName + ".data.contents");
                  if (pn.size() != 1 || static_cast<int>(m["parameters"].num_children()) != x->numBins() ||
                      static_cast<int>(vals.size()) != x->numBins() || static_cast<int>(jb) >= x->numBins())
                     throw std::runtime_error("CMSHistSum importer: '" + sampleName + "' is not one parameter in every bin");
                  for (unsigned k = 0; k < vals.size(); ++k)
                     if (k != jb && vals[k] != 0.)
                        throw std::runtime_error("CMSHistSum importer: '" + sampleName + "' is not 0 outside its bin");
                  if (!type1) {
                     type1 = std::make_unique<Additive>();
                     type1->pars.assign(x->numBins(), nullptr);
                     type1->scale.assign(x->numBins(), 0.);
                  }
                  if (type1->pars[jb]) throw std::runtime_error("CMSHistSum importer: two samples " + sampleName);
                  type1->pars[jb] = &getOrMakeParameter(tool, *pn.begin(), 0.0, -7.0, 7.0);
                  type1->scale[jb] = vals[jb];
                  continue;
               }
               if (m.has_child("constraint") || m.has_child("constraints"))
                  throw std::runtime_error("CMSHistSum importer: constrained shapefactor on sample '" + sampleName +
                                           "' is not supported");
               if (additive)
                  throw std::runtime_error("CMSHistSum importer: more than one shapefactor sample in '" + distName + "'");
               additive = std::make_unique<Additive>();
               additive->sample = sampleName;
               additive->modifier = m["name"].val();
               additive->scale = readDoubleArray(s["data"]["contents"], sampleName + ".data.contents");
               if (static_cast<int>(additive->scale.size()) != x->numBins() ||
                   static_cast<int>(m["parameters"].num_children()) != x->numBins())
                  throw std::runtime_error("CMSHistSum importer: shapefactor sample '" + sampleName + "' size != nbins");
               for (auto const &par : m["parameters"].children()) {
                  additive->pars.push_back(&getOrMakeParameter(tool, par.val(), 0.0, -7.0, 7.0));
               }
               continue;
            }
         }

         // data.contents (required)
         if (!s.has_child("data") || !s["data"].has_child("contents"))
            throw std::runtime_error("CMSHistSum importer: sample '" + sampleName + "' missing data.contents");
         std::vector<double> contents = readDoubleArray(s["data"]["contents"], sampleName + ".data.contents");
         if (static_cast<int>(contents.size()) != x->numBins())
            throw std::runtime_error("CMSHistSum importer: sample '" + sampleName + "' contents size != nbins");
         // Type-2 bins: the contents (and histosys, below) of bin j are in mcstat_<p>_bin<j>.
         std::vector<std::pair<unsigned, RooFit::Detail::JSONNode const *>> type2;
         if (auto it = perProc.find(sampleName); it != perProc.end()) type2 = it->second.type2;
         auto restoreBin = [&](std::vector<double> &vals, std::vector<double> const &bbVals, unsigned j,
                               std::string const &what) {
           if (bbVals.size() != vals.size() || j >= vals.size())
             throw std::runtime_error("CMSHistSum importer: " + what + ": size != nbins");
           for (unsigned k = 0; k < bbVals.size(); ++k)
             if (k != j && bbVals[k] != 0.)
               throw std::runtime_error("CMSHistSum importer: " + what + ": not 0 outside bin " + std::to_string(j));
           if (vals[j] != 0.)
             throw std::runtime_error("CMSHistSum importer: sample '" + sampleName + "' is not 0 in its type-2 bin " +
                                      std::to_string(j) + " (" + what + ")");
           vals[j] = bbVals[j];
         };
         std::vector<std::string> rateMods, histosysNames;
         for (auto const &[j, bb] : type2) {
           std::string const bbName = (*bb)["name"].val();
           restoreBin(contents, readDoubleArray((*bb)["data"]["contents"], bbName), j, bbName + ".data.contents");
         }

         // data.errors (optional, but we need it for staterror)
         std::vector<double> errors;
         if (s["data"].has_child("errors")) {
            errors = readDoubleArray(s["data"]["errors"], sampleName + ".data.errors");
            if (static_cast<int>(errors.size()) != x->numBins())
               throw std::runtime_error("CMSHistSum importer: sample '" + sampleName + "' errors size != nbins");
         } else {
	   errors.assign(x->numBins(), 0.0);
         }

         // Build a TH1D with that binning (uniform, defined by axis min/max/nbins)
         const std::string hname = distName + "_" + sampleName + "_nominal";
         auto h = std::make_unique<TH1D>(hname.c_str(), hname.c_str(), x->numBins(), x->getMin(), x->getMax());
         for (int ib = 1; ib <= x->numBins(); ++ib) {
            h->SetBinContent(ib, contents[ib - 1]);
            h->SetBinError(ib, errors[ib - 1]);
         }

	 const std::string funcName = distName + "_histosys_" + sampleName;
         auto *func = new CMSHistFunc(funcName.c_str(), funcName.c_str(), *x, *h);

	 struct VMorph {
	   FastHisto low;
	   FastHisto high;
	   VMorph(const FastHisto& low_, const FastHisto& high_):low(std::move(low_)),high(std::move(high_)){}
	 };
	 RooArgList vmorph_pars;
	 std::vector<VMorph> vmorphs;

	 // Rate modifiers, in order, become a product of ProcessNormalizations. exportCoeffModifier
	 // writes one as [symmetric normsys][asymmetric normsys][other factors][<sample>_norm_nominal],
	 // so a new segment starts whenever the next modifier cannot follow in that order; the
	 // re-export then has the same modifiers in the same order.
	 struct RateSegment {
	   enum Phase { Sym, Asym, Other, Nominal } phase = Sym;
	   std::vector<std::pair<RooAbsReal*, std::pair<double,double>>> normsys;
	   RooArgList others;
	   double nominal = 1.;
	 };
	 std::vector<RateSegment> segments;
	 auto segmentFor = [&segments](RateSegment::Phase phase) -> RateSegment& {
	   if (segments.empty() || segments.back().phase > phase) segments.emplace_back();
	   segments.back().phase = phase;
	   return segments.back();
	 };
	 const std::string normNominalName = sampleName + "_norm_nominal";

         if (s.has_child("modifiers")) {
	   for(const auto& m:s["modifiers"].children()){
	     const std::string mtype = m["type"].val();
	     const std::string mname = m["name"].val();
	     if (isRateModifier(m)) rateMods.push_back(jsonDump(m));

               if (mtype == "normfactor") {
                  // Typical: a floating POI or nuisance scaling a process yield
		  double value;
		  if(parse_double(m["parameter"].val(),value)){
		    // numeric: val() keeps only 6 significant digits, so read it with val_double()
		    value = m["parameter"].val_double();
		    if (mname == normNominalName) {
		      segmentFor(RateSegment::Nominal).nominal = value;
		    } else {
		      auto *c = new RooConstVar(mname.c_str(), mname.c_str(), value);
		      tool->workspace()->import(*c, RooFit::RecycleConflictNodes());
		      segmentFor(RateSegment::Other).others.add(*tool->workspace()->arg(mname.c_str()));
		    }
		  } else {
		    RooRealVar &p = getOrMakeParameter(tool, m["parameter"].val(), 1.0, 0.0, 10.0);
		    segmentFor(RateSegment::Other).others.add(p);
		  }
	       } else if (mtype == "custom"){
                  // This might be normalization factor that is given by a ProcessNormalization or RooFormulaVar
                  RooAbsReal *obj = tool->workspace()->function(mname);
		  if(!obj){
		    throw std::runtime_error(
					     "CMSHistSum importer: custom modifier '" + mname + "' on sample '" + sampleName +
					     "' not yet imported, there seems to be a priority/ordering issue");
		  }
		  segmentFor(RateSegment::Other).others.add(*obj);
               } else if (mtype == "staterror") {
                  // Already encoded via per-bin errors; CMSHistSum will handle BBB if configured.
                  continue;
	       } else if (mtype == "normsys"){
		 auto *parameter = m.find("parameter");
		 std::string parname(parameter ? parameter->val() : "alpha_" + mname);
		 RooRealVar &p = getOrMakeParameter(tool, parname, 0.0, -10., 10.0);
		 auto& data = m["data"];
		 double lo = data["lo"].val_double();
		 double hi = data["hi"].val_double();
		 const bool sym = approxSymmetricKappas(lo, hi);
		 segmentFor(sym ? RateSegment::Sym : RateSegment::Asym).normsys.push_back({&p, {lo, hi}});
	       } else if (mtype == "histosys"){
		 // collect all the histosys first, will assign them later when we know how many there are
		 auto *parameter = m.find("parameter");
		 std::string parname(parameter ? parameter->val() : "alpha_" + mname);
		 RooRealVar &p = getOrMakeParameter(tool, parname, 0.0, -10., 10.0);		 
		 vmorph_pars.add(p);
		 auto &data = m["data"];
		 auto lo = readDoubleArray(data["lo"]["contents"], mname);
		 auto hi = readDoubleArray(data["hi"]["contents"], mname);
		 for (auto const &[j, bb] : type2) {
		   RooFit::Detail::JSONNode const *bbm = nullptr;
		   for (auto const &mm : (*bb)["modifiers"].children())
		     if (mm["type"].val() == "histosys" && mm["name"].val() == mname) bbm = &mm;
		   std::string const what = std::string((*bb)["name"].val()) + " histosys " + mname;
		   if (!bbm) throw std::runtime_error("CMSHistSum importer: no " + what);
		   restoreBin(lo, readDoubleArray((*bbm)["data"]["lo"]["contents"], what), j, what + " lo");
		   restoreBin(hi, readDoubleArray((*bbm)["data"]["hi"]["contents"], what), j, what + " hi");
		 }
		 histosysNames.push_back(jsonDump(m["name"]) + jsonDump(m["parameter"]));
		 vmorphs.emplace_back(fastHistoFrom(lo, *x, mname + "Low_" + sampleName),
		                      fastHistoFrom(hi, *x, mname + "High_" + sampleName));
               } else {
                  // Shape and other modifiers: not safely mappable without CMSHistFunc’s public morphing API.
                  throw std::runtime_error(
                      "CMSHistSum importer: modifier type '" + mtype + "' on sample '" + sampleName +
                      "' is not supported by this CMSHistSum importer yet");
               }
	   }
         }

	 assert(vmorphs.size() == vmorph_pars.size());

	 // actually assign the morphs
	 func->setVerticalMorphs(vmorph_pars);
	 func->prepareStorage();
	 for(size_t i=0; i<vmorphs.size(); ++i){
	   func->setShape(0, 0, i+1, 0, vmorphs[i].low);
	   func->setShape(0, 0, i+1, 1, vmorphs[i].high);
	 }
         funcs.add(*func);

	 RooArgList coeffFactors;
	 for (std::size_t iseg = 0; iseg < segments.size(); ++iseg) {
	   RateSegment const &seg = segments[iseg];
	   const std::string pnName = distName + "_overallsys_" + sampleName + (iseg ? "_" + std::to_string(iseg) : "");
	   auto* pn = new ProcessNormalization(pnName.c_str(), pnName.c_str(), seg.nominal);
	   for (auto const &[par, kappas] : seg.normsys) {
	     // lo/hi are the relative factors at theta = -1/+1: symmetric kappa = hi (kappa^theta)
	     if (approxSymmetricKappas(kappas.first, kappas.second)) {
	       pn->addLogNormal(kappas.second, *par);
	     } else {
	       pn->addAsymmLogNormal(kappas.first, kappas.second, *par);
	     }
	   }
	   for (auto *other : seg.others) pn->addOtherFactor(static_cast<RooAbsReal &>(*other));
	   coeffFactors.add(*pn);
	 }

         // One coefficient object per sample, even for a single (possibly shared) factor: it
         // carries the sample name for the exporter (cms_hist_sum_process_name).
         const std::string coeffName = distName + "_" + sampleName + "_rate";
         auto *prod = new RooProduct(coeffName.c_str(), coeffName.c_str(), coeffFactors);
         prod->setStringAttribute("combine.process", sampleName.c_str());
         tool->workspace()->import(*prod, RooFit::RecycleConflictNodes());
         RooAbsReal *coeff = tool->workspace()->function(coeffName.c_str());
         if (!coeff)
            throw std::runtime_error("CMSHistSum importer: failed to import coefficient product '" + coeffName + "'");
         coeffs.add(*coeff);
         procIndex[sampleName] = funcs.getSize() - 1;

         // the per-process autoMCStats samples must carry exactly this process's rate modifiers
         // (and, for type 2, its histosys in the same order), then the normfactor 1/n and the shapefactor
         if (auto it = perProc.find(sampleName); it != perProc.end()) {
           auto check = [&](RooFit::Detail::JSONNode const &bb, bool withHistosys) {
             std::string const bbName = bb["name"].val();
             std::vector<RooFit::Detail::JSONNode const *> mods;
             for (auto const &mm : bb["modifiers"].children()) mods.push_back(&mm);
             std::size_t const tail = withHistosys ? 2 : 1;
             std::size_t const nh = withHistosys ? histosysNames.size() : 0;
             bool ok = mods.size() == rateMods.size() + nh + tail;
             for (std::size_t k = 0; ok && k < rateMods.size(); ++k) ok = jsonDump(*mods[k]) == rateMods[k];
             for (std::size_t k = 0; ok && k < nh; ++k) {
               auto const &mm = *mods[rateMods.size() + k];
               ok = mm["type"].val() == "histosys" && jsonDump(mm["name"]) + jsonDump(mm["parameter"]) == histosysNames[k];
             }
             if (ok && withHistosys) {
               auto const &nf = *mods[mods.size() - 2];
               double v;
               ok = nf["type"].val() == "normfactor" && parse_double(nf["parameter"].val(), v);
             }
             if (ok) ok = mods.back()->operator[]("type").val() == "shapefactor" && mods.back()->operator[]("name").val() == bbName;
             if (!ok)
               throw std::runtime_error("CMSHistSum importer: autoMCStats sample '" + bbName +
                                        "' does not carry the modifiers of process '" + sampleName + "'");
           };
           for (auto const &jb : it->second.type3) check(*jb.second, false);
           for (auto const &jb : it->second.type2) check(*jb.second, true);
         }
      }
      for (auto const &[proc, bb] : perProc) {
        if (!procIndex.count(proc))
          throw std::runtime_error("CMSHistSum importer: autoMCStats samples for unknown process '" + proc + "'");
      }

      // Per-process bin types and parameters (CMSHistSum::setImportedBinPars), bin by bin
      const unsigned nb = x->numBins();
      const unsigned np = funcs.getSize();
      std::vector<std::vector<unsigned>> binTypes(nb, std::vector<unsigned>(1, 0));
      std::vector<std::vector<RooAbsReal *>> binPars(nb, std::vector<RooAbsReal *>(np, nullptr));
      auto mark = [&](unsigned j, unsigned ip, unsigned type, RooAbsReal *par, std::string const &what) {
        if (j >= nb) throw std::runtime_error("CMSHistSum importer: " + what + ": bin out of range");
        if (binTypes[j].size() == 1) binTypes[j].assign(np, 4);
        if (binTypes[j][ip] != 4) throw std::runtime_error("CMSHistSum importer: " + what + ": bin " + std::to_string(j) + " twice");
        binTypes[j][ip] = type;
        binPars[j][ip] = par;
      };
      for (auto const &[proc, bb] : perProc) {
        const unsigned ip = procIndex[proc];
        auto const *hf = static_cast<CMSHistFunc const *>(funcs.at(ip));
        for (auto const &[j, node] : bb.type3) {
          std::string const bbName = (*node)["name"].val();
          auto const vals = readDoubleArray((*node)["data"]["contents"], bbName);
          auto const &sf = (*node)["modifiers"].child((*node)["modifiers"].num_children() - 1);
          std::set<std::string> names;
          for (auto const &pn : sf["parameters"].children()) names.insert(pn.val());
          if (vals.size() != nb || names.size() != 1 || sf["parameters"].num_children() != nb || j >= nb)
            throw std::runtime_error("CMSHistSum importer: " + bbName + ": not one parameter in every bin");
          for (unsigned k = 0; k < nb; ++k)
            if (k != j && vals[k] != 0.)
              throw std::runtime_error("CMSHistSum importer: " + bbName + ": not 0 outside bin " + std::to_string(j));
          if (vals[j] == 0. || vals[j] != hf->errors()[j] * hf->cache().GetWidth(j))
            throw std::runtime_error("CMSHistSum importer: " + bbName + ": contents of bin " + std::to_string(j) +
                                     " are not the error of '" + proc + "'");
          mark(j, ip, 3, &getOrMakeParameter(tool, *names.begin(), 0.0, -7.0, 7.0), bbName);
        }
        for (auto const &[j, node] : bb.type2) {
          std::string const bbName = (*node)["name"].val();
          auto const &sf = (*node)["modifiers"].child((*node)["modifiers"].num_children() - 1);
          std::set<std::string> names;
          for (auto const &pn : sf["parameters"].children()) names.insert(pn.val());
          if (names.size() != 1 || sf["parameters"].num_children() != nb)
            throw std::runtime_error("CMSHistSum importer: " + bbName + ": shapefactor must name one parameter in every bin");
          RooRealVar &v = getOrMakeParameter(tool, *names.begin(), 1.0, 0.0, 10.0);
          std::vector<RooFit::Detail::JSONNode const *> mods;
          for (auto const &mm : (*node)["modifiers"].children()) mods.push_back(&mm);
          auto const &nf = *mods[mods.size() - 2];
          const std::string cname = nf["name"].val();
          const double cval = nf["parameter"].val_double();
          if (!tool->workspace()->arg(cname.c_str())) {
            RooConstVar c(cname.c_str(), "", cval);
            tool->workspace()->import(c, RooFit::RecycleConflictNodes(), RooFit::Silence());
          }
          auto *c = dynamic_cast<RooConstVar *>(tool->workspace()->arg(cname.c_str()));
          if (!c || c->getVal() != cval)
            throw std::runtime_error("CMSHistSum importer: " + bbName + ": constant '" + cname + "' exists with another value");
          const std::string gname = std::string(v.GetName()) + "_prod";
          if (!tool->workspace()->function(gname.c_str())) {
            RooProduct g(gname.c_str(), "", RooArgList(v, *c));
            tool->workspace()->import(g, RooFit::RecycleConflictNodes(), RooFit::Silence());
          }
          auto *g = dynamic_cast<RooProduct *>(tool->workspace()->function(gname.c_str()));
          if (!g) throw std::runtime_error("CMSHistSum importer: failed to import '" + gname + "'");
          mark(j, ip, 2, g, bbName);
        }
      }
      RooArgList importedBinPars;
      for (unsigned j = 0; j < nb; ++j)
        for (auto *par : binPars[j]) if (par) importedBinPars.add(*par);

      // Construct CMSHistSum
      const std::string sumName = distName + "_cmshistsum";
      auto *sum = new CMSHistSum(sumName.c_str(), sumName.c_str(), *x, funcs, coeffs);
      if (importedBinPars.getSize()) sum->setImportedBinPars(binTypes, importedBinPars);
      if (additive && type1)
         throw std::runtime_error("CMSHistSum importer: '" + distName + "' has both an mcstat sample and mcstat_bin<j> samples");
      if (additive) sum->setAdditiveBinPars(additive->pars, additive->scale, additive->sample, additive->modifier);
      if (type1) sum->setAdditiveBinPars(type1->pars, type1->scale, "", "");
      tool->workspace()->import(*sum, RooFit::RecycleConflictNodes());
      sum = dynamic_cast<CMSHistSum *>(tool->workspace()->function(sumName.c_str()));
      if (!sum)
         throw std::runtime_error("CMSHistSum importer: failed to import CMSHistSum '" + sumName + "'");

      // Wrap into RooRealSumPdf with ONE coefficient
      // RooRealSumPdf(name,title, funcs, coeffs) uses:
      //   sum over i coeff_i * func_i ; here: one term: 1 * CMSHistSum
      RooArgList pdfFuncs;
      RooArgList pdfCoeffs;
      pdfFuncs.add(*sum);
      pdfCoeffs.add(*ONE);
      
      // The channel itself; wrapped below into RooProdPdf(distName) with its modifiers' constraints, if any.
      const std::string coreName = CombineHS3::channelCoreName(elem);
      auto *pdf = new RooRealSumPdf(coreName.c_str(), distName.c_str(), pdfFuncs, pdfCoeffs, true /*extended*/);
      tool->workspace()->import(*pdf, RooFit::RecycleConflictNodes());
      CombineHS3::wrapChannelInConstraints(tool, elem);

      return true;
   }
};

STATIC_EXECUTE([]() {
   using namespace RooFit::JSONIO;

   // priority=true so this wins over the standard HistFactoryImporter
   registerImporter<CMSHistSumHistFactoryImporter>("histfactory_dist", /*priority*/ true);
});

} // namespace
