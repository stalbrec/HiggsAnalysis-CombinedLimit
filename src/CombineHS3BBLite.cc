#include "../interface/CombineHS3BBLite.h"

#include "../interface/AsymPow.h"
#include "../interface/CMSHistErrorPropagator.h"
#include "../interface/CMSHistSum.h"
#include "../interface/CMSHistFunc.h"
#include "../interface/ProcessNormalization.h"
#include "../interface/SimpleGaussianConstraint.h"

#include <RooArgList.h>
#include <RooConstVar.h>
#include <RooFormulaVar.h>
#include <RooNumber.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooRealVar.h>
#include <RooStats/ModelConfig.h>
#include <RooWorkspace.h>

#include <cmath>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {

CombineHS3::BBLiteExport const *gActive = nullptr;

std::string literal(double v)
{
   std::ostringstream os;
   os << std::setprecision(std::numeric_limits<double>::max_digits10) << v;
   std::string s = os.str();
   // TFormula needs a floating-point literal where min/max meet a double; HS3 accepts both
   if (s.find_first_of(".eEni") == std::string::npos) s += ".0";
   return v < 0 ? "(" + s + ")" : s;
}

struct Expr {
   std::string text;
   RooArgList deps;
   RooWorkspace *ws = nullptr; // where the shared <theta>_asym_alpha functions go
   void addDep(RooAbsArg const &a)
   {
      if (!deps.find(a.GetName())) deps.add(a);
   }
};

template <class T>
T *importAndGet(RooWorkspace &ws, T const &obj);

/// Combine's logKappaForX (CombineMathFuncs.h) times theta, without branches: for |theta| >= 0.5
/// the clamped polynomial is exactly +-1, so avg +- halfdiff gives logKappaHigh / -logKappaLow
/// (up to rounding of avg +- halfdiff). The polynomial depends on theta only: one shared function
/// <theta>_asym_alpha per nuisance keeps the c_p expressions short (each one is JIT-compiled).
std::string asymTerm(RooAbsArg const &th, double logKappaLow, double logKappaHigh, Expr &e)
{
   double logKhi = logKappaHigh;
   double logKlo = -logKappaLow;
   double avg = 0.5 * (logKhi + logKlo);
   double halfdiff = 0.5 * (logKhi - logKlo);
   std::string const theta = th.GetName();
   std::string const aname = theta + "_asym_alpha";
   RooAbsReal *alpha = e.ws->function(aname.c_str());
   if (!alpha) {
      std::string t = "max(-1.0,min(1.0," + theta + "+" + theta + "))";
      std::string t2 = "(" + t + "*" + t + ")";
      RooFormulaVar f(aname.c_str(), aname.c_str(), ("0.125*" + t + "*(" + t2 + "*(3.0*" + t2 + "-10.0)+15.0)").c_str(),
                      RooArgList(th));
      alpha = importAndGet(*e.ws, f);
   }
   e.addDep(th);
   e.addDep(*alpha);
   return theta + "*(" + literal(avg) + "+" + aname + "*" + literal(halfdiff) + ")";
}

/// c_p as a generic expression, mirroring the evaluation order of the Combine classes.
void coeffExpr(RooAbsArg const *c, Expr &e)
{
   if (auto *rcv = dynamic_cast<RooConstVar const *>(c)) {
      e.text += literal(rcv->getVal());
   } else if (auto *pn = dynamic_cast<ProcessNormalization const *>(c)) {
      // CombineMathFuncs.h processNormalization: nominal * exp(sum) * others
      std::string logVal;
      for (std::size_t i = 0; i < pn->logKappa().size(); ++i) {
         RooAbsArg const *th = pn->thetaList().at(i);
         e.addDep(*th);
         logVal += (logVal.empty() ? "" : "+") + std::string(th->GetName()) + "*" + literal(pn->logKappa()[i]);
      }
      for (std::size_t i = 0; i < pn->logAsymmKappa().size(); ++i) {
         RooAbsArg const *th = pn->asymmThetaList().at(i);
         e.addDep(*th);
         logVal += (logVal.empty() ? "" : "+") +
                   asymTerm(*th, pn->logAsymmKappa()[i].first, pn->logAsymmKappa()[i].second, e);
      }
      e.text += literal(pn->nominalValue());
      if (!logVal.empty()) e.text += "*exp(" + logVal + ")";
      for (RooAbsArg const *o : pn->otherFactorList()) {
         e.text += "*(";
         coeffExpr(o, e);
         e.text += ")";
      }
   } else if (auto *ap = dynamic_cast<AsymPow const *>(c)) {
      // CombineMathFuncs.h asymPow: exp(logKappaForX(theta, log kLo, log kHi) * theta)
      if (!ap->kappaLow().isConstant() || !ap->kappaHigh().isConstant())
         throw std::runtime_error(std::string("HS3 autoMCStats export: AsymPow ") + ap->GetName() +
                                  " has non-constant kappas");
      e.addDep(ap->theta());
      e.text += "exp(" + asymTerm(ap->theta(), std::log(ap->kappaLow().getVal()),
                                  std::log(ap->kappaHigh().getVal()), e) + ")";
   } else if (auto *prod = dynamic_cast<RooProduct const *>(c)) {
      RooArgList comps = const_cast<RooProduct *>(prod)->components();
      if (comps.empty()) {
         e.text += "1.0";
         return;
      }
      for (int i = 0; i < comps.getSize(); ++i) {
         e.text += i ? "*(" : "(";
         coeffExpr(comps.at(i), e);
         e.text += ")";
      }
   } else {
      // a parameter, or a function exported under its own name (rateParams, ...)
      e.addDep(*c);
      e.text += c->GetName();
   }
}

std::string processName(CMSHistFunc const &f)
{
   return f.stringAttributes().count("combine.process") ? f.getStringAttribute("combine.process") : f.GetName();
}

template <class T>
T *importAndGet(RooWorkspace &ws, T const &obj)
{
   ws.import(obj, RooFit::RecycleConflictNodes(), RooFit::Silence());
   auto *out = dynamic_cast<T *>(ws.arg(obj.GetName()));
   if (!out) throw std::runtime_error(std::string("HS3 autoMCStats export: failed to import ") + obj.GetName());
   return out;
}

/// The part of `name` after `primary` + "_" (prop_binch_bonly vs prop_binch -> "bonly"), else `name`.
std::string suffixOf(std::string const &name, std::string const &primary)
{
   if (!primary.empty() && name.size() > primary.size() + 1 && name.compare(0, primary.size() + 1, primary + "_") == 0)
      return name.substr(primary.size() + 1);
   return name;
}

/// `f` with the constraints in `subst` replaced, recursively through RooProdPdfs (new products are
/// named <name><tag>); `f` itself if nothing changes.
RooAbsPdf *substituted(RooWorkspace &ws, RooAbsPdf *f, std::map<std::string, RooAbsPdf *> const &subst,
                       std::string const &tag)
{
   if (auto it = subst.find(f->GetName()); it != subst.end()) return it->second;
   auto *prod = dynamic_cast<RooProdPdf *>(f);
   if (!prod) return f;
   RooArgList factors;
   bool changed = false;
   for (RooAbsArg *a : prod->pdfList()) {
      auto *fa = static_cast<RooAbsPdf *>(a);
      RooAbsPdf *g = substituted(ws, fa, subst, tag);
      changed |= g != fa;
      factors.add(*g);
   }
   if (!changed) return f;
   std::string const name = f->GetName() + tag;
   if (auto *existing = ws.pdf(name.c_str())) return existing;
   RooProdPdf copy(name.c_str(), f->GetTitle(), factors);
   return importAndGet(ws, copy);
}

} // namespace

namespace CombineHS3 {

BBLiteExport::BBLiteExport(RooWorkspace &ws) : previous_(gActive)
{
   // The objects with autoMCStats bin parameters: CMSHistErrorPropagator or CMSHistSum, read through
   // the same accessors (errors in events; both require unit bin widths below).
   struct Source {
      RooAbsReal *obj;
      std::vector<std::vector<unsigned>> const *bt;
      RooArgList const *pars;
      unsigned nProcs;
      std::function<double(unsigned, unsigned)> err; // (process, bin)
      std::function<RooAbsArg const *(unsigned)> coeff;
      std::function<std::string(unsigned)> process;
      std::function<double(unsigned)> width;
   };
   std::vector<Source> props;
   for (RooAbsArg *arg : ws.components()) {
      if (auto *p = dynamic_cast<CMSHistErrorPropagator *>(arg); p && p->hasBinPars()) {
         RooArgList const &funcs = p->funcList();
         std::vector<CMSHistFunc const *> hf(funcs.getSize());
         for (int ip = 0; ip < funcs.getSize(); ++ip) {
            hf[ip] = dynamic_cast<CMSHistFunc const *>(funcs.at(ip));
            if (!hf[ip])
               throw std::runtime_error(std::string("HS3 autoMCStats export: ") + p->GetName() + ": not a CMSHistFunc");
         }
         for (unsigned j = 0; j < p->binTypes().size(); ++j) {
            for (unsigned t : p->binTypes()[j])
               if (t == 2 || t == 3)
                  throw std::runtime_error(std::string("HS3 autoMCStats export: ") + p->GetName() + " bin " +
                                           std::to_string(j) +
                                           " uses per-process bin parameters (types 2/3), written for CMSHistSum only");
         }
         p->getVal(); // fills the caches
         props.push_back({p, &p->binTypes(), &p->binParList(), static_cast<unsigned>(funcs.getSize()),
                          [hf](unsigned ip, unsigned j) { return hf[ip]->errors()[j]; },
                          [p](unsigned ip) -> RooAbsArg const * { return p->coefList().at(ip); },
                          [hf](unsigned ip) { return processName(*hf[ip]); },
                          [p](unsigned j) { return p->cache().GetWidth(j); }});
      } else if (auto *h = dynamic_cast<CMSHistSum *>(arg); h && h->hasBinPars()) {
         h->getVal();
         props.push_back({h, &h->binTypes(), &h->binParList(), h->nProcs(),
                          [h](unsigned ip, unsigned j) { return h->binErrors(ip, true)[j]; },
                          [h](unsigned ip) -> RooAbsArg const * { return h->coeffAt(ip); },
                          [h](unsigned ip) {
                             RooAbsReal const *c = h->coeffAt(ip);
                             return c->stringAttributes().count("combine.process")
                                        ? c->getStringAttribute("combine.process")
                                        : std::string(c->GetName());
                          },
                          [h](unsigned j) { return h->binWidth(j); }});
      }
   }
   // The s+b model's propagators first: they define y_j and keep the plain names.
   if (auto *mc = dynamic_cast<RooStats::ModelConfig *>(ws.genobj("ModelConfig")); mc && mc->GetPdf()) {
      std::unique_ptr<RooArgSet> comps{mc->GetPdf()->getComponents()};
      std::stable_partition(props.begin(), props.end(), [&](auto const &p) { return comps->find(*p.obj) != nullptr; });
   }

   std::map<std::string, RooRealVar *> xToY;
   std::map<std::string, std::string> eExprOf;   // x_j -> sum in E_j of the propagator that defined y_j
   std::map<std::string, std::string> definedBy; // x_j -> that propagator
   std::map<std::string, std::string> normByExpr; // c_p expression -> function name
   std::map<std::string, std::map<std::string, RooAbsPdf *>> substitutes; // propagator -> constraint -> its own copy
   std::map<std::string, std::string> tagOf; // propagator -> suffix of its own copies

   for (Source const &src : props) {
      RooAbsReal *prop = src.obj;
      std::string const pname = prop->GetName();
      auto const &bt = *src.bt;
      RooArgList const &pars = *src.pars;
      // Bins with one parameter for the total (type 1) are written here; per-process bins (types
      // 2, 3, 4) keep their parameters and are written by the histfactory exporter.
      bool anyType1 = false;
      for (unsigned j = 0; j < bt.size(); ++j) anyType1 |= bt[j][0] == 1;
      if (!anyType1) continue;
      for (unsigned j = 0; j < bt.size(); ++j) {
         if (src.width(j) != 1.0)
            throw std::runtime_error("HS3 autoMCStats export: " + pname + " has bins of width != 1");
      }

      // c_p for every process that has a non-zero error somewhere
      std::vector<std::string> cName(src.nProcs);
      for (unsigned ip = 0; ip < src.nProcs; ++ip) {
         bool anyErr = false; // in a type-1 bin
         for (unsigned j = 0; j < bt.size(); ++j) anyErr |= bt[j][0] == 1 && src.err(ip, j) != 0.;
         if (!anyErr) continue;
         Expr e;
         e.ws = &ws;
         coeffExpr(src.coeff(ip), e);
         if (e.deps.empty()) {
            cName[ip] = "(" + e.text + ")";
         } else if (e.deps.getSize() == 1 && e.text == e.deps.at(0)->GetName()) {
            cName[ip] = e.text;
         } else if (auto it = normByExpr.find(e.text); it != normByExpr.end()) {
            cName[ip] = it->second;
         } else {
            std::string const name = pname + "_norm_" + src.process(ip);
            RooFormulaVar f(name.c_str(), name.c_str(), e.text.c_str(), e.deps);
            cName[ip] = normByExpr[e.text] = importAndGet(ws, f)->GetName();
         }
      }

      std::vector<RooRealVar *> &ys = ys_[pname];
      ys.assign(bt.size(), nullptr);
      std::vector<RooAbsPdf *> &cons = cons_[pname];
      cons.assign(bt.size(), nullptr);
      for (unsigned j = 0, r = 0; j < bt.size(); ++j) {
         if (bt[j][0] != 1) { // no term, or per-process parameters
            for (unsigned t : bt[j]) r += t >= 1 && t < 4;
            continue;
         }
         auto *x = dynamic_cast<RooRealVar *>(pars.at(r++));
         if (!x) throw std::runtime_error("HS3 autoMCStats export: " + pname + ": bin parameter is not a RooRealVar");
         std::string const xname = x->GetName();
         // E_j, summed as CMSHistErrorPropagator::updateCache does (vectorized::mul_add_sqr)
         std::string sum;
         RooArgList eDeps;
         for (unsigned ip = 0; ip < src.nProcs; ++ip) {
            double err = src.err(ip, j);
            if (err == 0.) continue;
            std::string const c = cName[ip];
            sum += (sum.empty() ? "" : "+") + c + "*" + c + "*" + literal(err) + "*" + literal(err);
            if (c[0] != '(' && !eDeps.find(c.c_str())) eDeps.add(*ws.arg(c.c_str()));
         }
         // the constraints Gauss(x_In | x_j, 1) of x_j
         std::vector<SimpleGaussianConstraint *> constraints;
         for (RooAbsArg *client : x->clients()) {
            auto *g = dynamic_cast<SimpleGaussianConstraint *>(client);
            if (g && &g->getX() == x && ws.pdf(g->GetName()) == g) constraints.push_back(g);
         }

         auto yIt = xToY.find(xname);
         if (!constraints.empty()) cons[j] = constraints.front();
         if (yIt != xToY.end() && eExprOf[xname] == sum) { // same E_j as where y_j was defined
            ys[j] = yIt->second;
            continue;
         }
         // y_j defined here, or (another model sharing x_j, e.g. the b-only one) its own E_j, mean
         // and constraint copies, which its channel uses instead of the shared ones.
         std::string const tag = yIt == xToY.end() ? "" : "_" + suffixOf(pname, definedBy[xname]);
         std::string const ename = xname + "_E" + tag;
         RooFormulaVar E(ename.c_str(), ename.c_str(), ("sqrt(" + sum + ")").c_str(), eDeps, false);
         RooAbsReal *wsE = importAndGet(ws, E);
         double const E0 = wsE->getVal();
         if (!(E0 > 0.))
            throw std::runtime_error("HS3 autoMCStats export: " + ename + " is not positive at the current values");
         RooRealVar *wsY = nullptr;
         if (yIt == xToY.end()) {
            std::string const yname = xname + "_y";
            // Unbounded: x_j's fixed range is a theta-dependent range in y_j, which no fixed range matches.
            RooRealVar y(yname.c_str(), yname.c_str(), x->getVal() * E0, -RooNumber::infinity(), RooNumber::infinity());
            wsY = importAndGet(ws, y);
            // ROOT writes a step width for every workspace variable with an error: x_j's moves to y_j.
            if (x->hasError()) {
               wsY->setError(x->getError() * E0);
               savedErrors_.emplace_back(x, x->getError());
               x->removeError();
            }
            xToY[xname] = wsY;
            eExprOf[xname] = sum;
            definedBy[xname] = pname;
         } else {
            wsY = yIt->second;
         }
         std::string const mname = xname + "_mean" + tag;
         RooFormulaVar mean(mname.c_str(), mname.c_str(), (std::string(wsY->GetName()) + "/" + ename).c_str(),
                            RooArgList(*wsY, *wsE));
         RooAbsReal *wsMean = importAndGet(ws, mean);
         for (SimpleGaussianConstraint *g : constraints) {
            if (tag.empty()) {
               means_[g->GetName()] = wsMean;
            } else {
               std::string const cname = g->GetName() + tag;
               SimpleGaussianConstraint copy(cname.c_str(), g->GetTitle(), *x, const_cast<RooAbsReal &>(g->getMean()),
                                             const_cast<RooAbsReal &>(g->getSigma()));
               RooAbsPdf *wsCopy = importAndGet(ws, copy);
               means_[cname] = wsMean;
               substitutes[pname][g->GetName()] = wsCopy;
               if (g == constraints.front()) cons[j] = wsCopy;
               if (!tagOf.emplace(pname, tag).second && tagOf[pname] != tag)
                  throw std::runtime_error("HS3 autoMCStats export: " + pname +
                                           " shares bin parameters with two models");
            }
         }
         ys[j] = wsY;
      }
   }

   // Channels of a propagator with its own constraint copies: replace the shared constraints (or the
   // products holding them) among the factors of every RooProdPdf above the propagator.
   for (auto const &[pname, subst] : substitutes) {
      RooAbsArg *prop = ws.arg(pname.c_str());
      std::string const &tag = tagOf[pname];
      for (RooAbsArg *arg : ws.components()) {
         auto *q = dynamic_cast<RooProdPdf *>(arg);
         if (!q || !q->dependsOn(*prop)) continue;
         for (RooAbsArg *f : RooArgList(q->pdfList())) {
            auto *fp = static_cast<RooAbsPdf *>(f);
            if (fp->dependsOn(*prop)) continue;
            RooAbsPdf *g = substituted(ws, fp, subst, tag);
            if (g == fp) continue;
            std::string const attr = std::string("ORIGNAME:") + fp->GetName();
            g->setAttribute(attr.c_str());
            q->redirectServers(RooArgSet(*g), false, true);
            g->setAttribute(attr.c_str(), false);
            redirects_.push_back({q, fp, g});
         }
      }
   }

   // nuisance parameters: x_j -> y_j
   for (TObject *obj : ws.allGenericObjects()) {
      auto *mc = dynamic_cast<RooStats::ModelConfig *>(obj);
      if (!mc || !mc->GetNuisanceParameters()) continue;
      RooArgSet const &nuis = *mc->GetNuisanceParameters();
      bool any = false;
      RooArgSet swapped;
      for (RooAbsArg *a : nuis) {
         auto it = xToY.find(a->GetName());
         any |= it != xToY.end();
         swapped.add(it != xToY.end() ? static_cast<RooAbsArg &>(*it->second) : *a);
      }
      if (!any) continue;
      savedNuisances_.emplace_back(mc, new RooArgSet(nuis));
      mc->SetNuisanceParameters(swapped);
   }
   gActive = this;
}

BBLiteExport::~BBLiteExport()
{
   for (auto it = redirects_.rbegin(); it != redirects_.rend(); ++it) {
      std::string const attr = std::string("ORIGNAME:") + it->replacement->GetName();
      it->original->setAttribute(attr.c_str());
      it->product->redirectServers(RooArgSet(*it->original), false, true);
      it->original->setAttribute(attr.c_str(), false);
   }
   for (auto &[x, err] : savedErrors_) x->setError(err);
   for (auto &[mc, saved] : savedNuisances_) {
      mc->SetNuisanceParameters(*saved);
      delete saved;
   }
   gActive = previous_;
}

BBLiteExport const *BBLiteExport::active()
{
   return gActive;
}

std::vector<RooRealVar *> const *BBLiteExport::binParameters(RooAbsArg const &prop) const
{
   auto it = ys_.find(prop.GetName());
   return it == ys_.end() ? nullptr : &it->second;
}

std::vector<RooAbsPdf *> const *BBLiteExport::binConstraints(RooAbsArg const &prop) const
{
   auto it = cons_.find(prop.GetName());
   return it == cons_.end() ? nullptr : &it->second;
}

RooAbsReal *BBLiteExport::meanFor(RooAbsArg const &constraint) const
{
   auto it = means_.find(constraint.GetName());
   return it == means_.end() ? nullptr : it->second;
}

} // namespace CombineHS3
