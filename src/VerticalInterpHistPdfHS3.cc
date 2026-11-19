// HS3 JSON export and import for FastVerticalInterpHistPdf2 (RooFit::JSONIO).
#include "../interface/VerticalInterpHistPdf.h"

#include <cassert>
#include <memory>

#include "RooFit.h"
#include "Riostream.h"

#include "RooRealVar.h"
#include "RooMsgService.h"
#include "RooAbsData.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "../interface/utils.h"

#include <set>

#include "RooFitHS3/RooJSONFactoryWSTool.h"
#include "RooFitHS3/JSONIO.h"
#include "CombineHS3Helpers.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooProduct.h"
#include "../interface/ProcessNormalization.h"
#include "../interface/AsymPow.h"
#include "../interface/CombineHS3Channel.h"
#include <charconv>
#include <map>
#include <string>
#include <vector>

using RooFit::Detail::JSONNode;


namespace {

using CombineHS3::approxSymmetricKappas;
using CombineHS3::getOrMakeParameter;
using CombineHS3::isSingleBinCountingDist;
using CombineHS3::parse_double;
using CombineHS3::readDoubleArray;
using CombineHS3::requirePlainModifierParameter;
using CombineHS3::writeNormSys;

/// Kept byte-compatible with the copy in JSONFactories_Combine.cxx: the fallback branch must also
/// write "parameter", because the importer below reads it unconditionally via parse_double().
// void exportOtherFactorModifiers(RooJSONFactoryWSTool *tool, JSONNode &mods, RooAbsReal const *factor) {
//   auto &m = mods.append_child().set_map();
//   m["name"] << factor->GetName();
//   m["type"] << "normfactor";
//   if (auto const *rcv = dynamic_cast<RooConstVar const *>(factor)) {
//     m["parameter"] << rcv->getVal();
//   } else {
//     m["parameter"] << factor->GetName();
//     tool->queueExport(*factor);
//   }
// }
void exportOtherFactorModifiers(RooJSONFactoryWSTool* tool, JSONNode& mods, RooAbsReal const* factor) {
  if (auto* rcv = dynamic_cast<RooConstVar const*>(factor)) {
    auto& m = mods.append_child().set_map();
    m["name"] << rcv->GetName();
    m["parameter"] << rcv->getVal();
    m["type"] << "normfactor";
  } else if (auto* rrv = dynamic_cast<RooRealVar const*>(factor)) {
    auto& m = mods.append_child().set_map();
    m["name"] << rrv->GetName();
    m["parameter"] << rrv->GetName();
    m["type"] << "normfactor";
    tool->queueExport(*rrv);
  } else if (auto* prod = dynamic_cast<RooProduct const*>(factor)) {
    for (auto* comp : prod->realComponents()) {
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


// Helper: try to pick a reasonable process name.
static std::string cms_fast_vertinterp_hist_pdf_process_name(const FastVerticalInterpHistPdf2 * pdf, unsigned ip)
{
  
  if (pdf && pdf->stringAttributes().count("combine.process")) {
    return pdf->getStringAttribute("combine.process");
  }
  
  if(pdf){
    return pdf->GetName();
  }
  
  return "proc_"+std::to_string(ip);
}


static std::unique_ptr<TH1D> readTH1D(RooFit::Detail::JSONNode const &arr, RooRealVar x, std::string hname)
{
    std::vector<double> contents = readDoubleArray(arr);
    
    // if (static_cast<int>(contents.size()) != x.numBins()) throw std::runtime_error("Sum of FastVerticalInterpHistPdf2 importer: sample '" + hname + "' contents size != nbins");
    // data.errors can't be used with FastVerticalInterpHistPdf2
    // Avoid registering with gDirectory: these histograms are owned by the
    // returned unique_ptr, and repeated calls with the same hname would
    // otherwise trigger ROOT's "Replacing existing TH1" warnings.
    bool oldAddDirStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
    auto h = std::make_unique<TH1D>(hname.c_str(), hname.c_str(), x.numBins(), x.getMin(), x.getMax());
    TH1::AddDirectory(oldAddDirStatus);
    for(int ib=1; ib<=x.numBins();++ib){
        h->SetBinContent(ib, contents[ib -1]);
    }
    return h;
}


/// Per-export state threaded through the modifier emitters below.
struct HistFactoryExportCtx {
   RooArgSet observables;                                //!< channel observables, for constraint lookup
   RooAbsArg const *channel = nullptr;                   //!< the exported channel pdf
   std::map<std::string, std::string> constraintCache;   //!< theta name -> constraint pdf name ("" = none)
   std::map<std::string, std::string> sampleConstraints; //!< parameter -> constraint name, reset per sample
   std::string sampleName;
};


std::string resolveConstraintPdfName(RooJSONFactoryWSTool *tool, HistFactoryExportCtx &ctx, RooAbsArg const &theta)
{
   return CombineHS3::resolveConstraintPdfName(tool, theta, ctx.observables, ctx.channel, ctx.constraintCache);
}


void registerConstraint(HistFactoryExportCtx &ctx, std::string const &param, std::string const &constraint)
{
   auto res = ctx.sampleConstraints.emplace(param, constraint);
   if (!res.second && res.first->second != constraint) {
      RooJSONFactoryWSTool::error("SumOfFastVerticalInterpHistPdf2: sample '" + ctx.sampleName + "': parameter '" +
                                  param + "' is referenced by modifiers with inconsistent constraints ('" +
                                  res.first->second + "' vs '" + constraint +
                                  "'); refusing to emit inconsistent output.");
   }
}

void exportNormSys(RooJSONFactoryWSTool *tool, HistFactoryExportCtx &ctx, JSONNode &mods, RooAbsArg const &theta,
                   double kappaLo, double kappaHi)
{
   requirePlainModifierParameter(theta, ctx.channel ? ctx.channel->GetName() : "", ctx.sampleName, "normsys");
   const std::string constraint = resolveConstraintPdfName(tool, ctx, theta);
   registerConstraint(ctx, theta.GetName(), constraint);
   writeNormSys(tool, mods, theta, kappaLo, kappaHi, constraint);
}

/// The pieces of a RooAddPdf coefficient this exporter knows how to write.
struct CoefficientParts {
   ProcessNormalization const *pn = nullptr;     //!< at most one, and optional -- see below
   std::vector<AsymPow const *> asymPows;        //!< shape-systematic rate terms (ShapeTools' systeff_*)
   std::vector<RooAbsReal const *> plainFactors; //!< bare RooRealVar / RooConstVar factors
   RooAbsArg const *customCoeff = nullptr;       //!< bare coefficient of an unrecognised type
};


bool decomposeCoefficient(RooAbsArg *coeff, CoefficientParts &out, std::string &err)
{
   auto classifyProductFactor = [&out, &err](RooAbsArg const *factor) {
      if (auto const *pn = dynamic_cast<ProcessNormalization const *>(factor)) {
         if (out.pn) {
            err = "more than one ProcessNormalization ('" + std::string(out.pn->GetName()) + "' and '" +
                  factor->GetName() + "')";
            return false;
         }
         out.pn = pn;
         return true;
      }
      if (auto const *ap = dynamic_cast<AsymPow const *>(factor)) {
         out.asymPows.push_back(ap);
         return true;
      }
      if (dynamic_cast<RooRealVar const *>(factor) || dynamic_cast<RooConstVar const *>(factor)) {
         out.plainFactors.push_back(static_cast<RooAbsReal const *>(factor));
         return true;
      }
      err = "component '" + std::string(factor->GetName()) + "' of type " + factor->ClassName() +
            ", which is neither a ProcessNormalization nor an AsymPow nor a plain variable";
      return false;
   };

   if (auto *prod = dynamic_cast<RooProduct *>(coeff)) {
      for (auto *comp : prod->components()) {
         if (!classifyProductFactor(comp)) return false;
      }
      return true;
   }

   if (auto const *pn = dynamic_cast<ProcessNormalization const *>(coeff)) {
      out.pn = pn;
   } else if (dynamic_cast<RooRealVar const *>(coeff) || dynamic_cast<RooConstVar const *>(coeff)) {
      out.plainFactors.push_back(static_cast<RooAbsReal const *>(coeff));
   } else {
      out.customCoeff = coeff;
   }
   return true;
}

// Number of template bins to export: the smallest active size. text2workspace (ShapeTools.py)
// calls setActiveBins() on morphing samples whose observable has fewer bins than their templates
// (part4: 31-bin templates, observables with 1 and 18 bins); the active size is persisted, and
// FastVerticalInterpHistPdf2 crops (to 1e-9) and normalizes over its active bins only. Samples
// without morphs keep all bins active, so their extra bins do take part in the normalization;
// the exporter writes what they evaluate to instead (see croppedHead). A morphing sample with
// more active bins than that has a parameter-dependent normalization over bins outside the axis,
// which histfactory_dist cannot express, so it is declined.
static int activeBins(RooArgList const &summands, std::string const &name)
{
    int n = -1;
    for (auto const *s : summands) {
        int size = static_cast<FastVerticalInterpHistPdf2 const *>(s)->cacheNominal().size();
        if (n < 0 || size < n) n = size;
    }
    for (auto const *s : summands) {
        auto const &pdf = *static_cast<FastVerticalInterpHistPdf2 const *>(s);
        if (!pdf.morphs().empty() && int(pdf.cacheNominal().size()) > n) {
            RooJSONFactoryWSTool::error("SumOfFastVerticalInterpHistPdf2 '" + name + "': morphing sample '" +
                                        pdf.GetName() + "' has more active bins than the " + std::to_string(n) +
                                        " shared by all samples.");
        }
    }
    return n;
}

// What an unmorphed FastVerticalInterpHistPdf2 evaluates to in its first nBins bins (syncTotal:
// crop to 1e-9 and normalize over all active bins), and the sum of those values. The sum is below
// 1 when active bins lie beyond nBins; the reimported sample is normalized over nBins only, so the
// exporter moves the sum into the sample's yield.
static std::vector<double> croppedHead(FastHisto const &nominal, int nBins, double &headSum)
{
    double total = 0.;
    for (unsigned ib = 0; ib < nominal.size(); ++ib) total += std::max(nominal[ib], 1e-9) * nominal.GetWidth(ib);
    std::vector<double> out(nBins);
    headSum = 0.;
    for (int ib = 0; ib < nBins; ++ib) {
        out[ib] = std::max(nominal[ib], 1e-9) / total;
        headSum += out[ib] * nominal.GetWidth(ib);
    }
    return out;
}

// Sum of a template's first n bins, in the same order wherever it is computed, so that the
// exporter's 1/S and the importer's check of it see the same S.
static double sampleSum(double const *values, int n)
{
    double sum = 0.;
    for (int i = 0; i < n; ++i) sum += values[i];
    return sum;
}

bool tryExportHistFactory(RooJSONFactoryWSTool *tool, const RooAddPdf *sum, JSONNode &elem){
    auto summands = sum->pdfList();
    const unsigned nProcs = summands.getSize();
    if(nProcs==0){
        RooJSONFactoryWSTool::warning("SumOfFastVerticalInterpHistPdf2 '" + std::string(sum->GetName()) + "': no Processes found!");
        return false;
    } 

    std::set<RooRealVar const *> observables;
    for(auto summand : summands){
        auto const *hpdf = dynamic_cast<FastVerticalInterpHistPdf2 const *>(summand);
        // This exports RooAddPdfs of FastVerticalInterpHistPdf2, so if anything else is present, exit here.
        if(!hpdf) {
            RooJSONFactoryWSTool::warning("SumOfFastVerticalInterpHistPdf2 '" + std::string(sum->GetName()) + "': summand was something other than FastVerticalInterpHistPdf2!");
            return false;
        }
        observables.insert(&hpdf->x());
    }
    if(observables.size() != 1){
        RooJSONFactoryWSTool::error("SumOfFastVerticalInterpHistPdf2 '" + std::string(sum->GetName()) +
                                     "': summands do not share a common observable.");
    }

    auto const *x  =  dynamic_cast<RooRealVar const *>(*observables.begin());
    if (!x) {
        RooJSONFactoryWSTool::error("SumOfFastVerticalInterpHistPdf2 '" + std::string(sum->GetName()) +
                                    "': shared observable is not a RooRealVar; cannot export axis metadata.");
    }

    // Decline pass: every coefficient is decomposed before a single byte is written into elem, so
    // that a decline leaves neither a half-filled node nor queued servers behind. The same contract
    // is documented for tryExportCountingHistFactory in JSONFactories_Combine.cxx.
    std::vector<CoefficientParts> coefParts(nProcs);
    for (unsigned ip = 0; ip < nProcs; ++ip) {
        auto *coeff = sum->coefList().at(ip);
        if (!coeff) continue;
        std::string err;
        if (!decomposeCoefficient(coeff, coefParts[ip], err)) {
            RooJSONFactoryWSTool::warning("SumOfFastVerticalInterpHistPdf2 '" + std::string(sum->GetName()) +
                                          "': coefficient '" + coeff->GetName() + "' has " + err +
                                          "; declining this export path.");
            return false;
        }
    }

    HistFactoryExportCtx ctx;
    ctx.observables.add(*x);
    ctx.channel = sum;

    auto &axes = elem["axes"].set_seq();
    auto &ax = axes.append_child().set_map();

    std::string axname(x->GetName());
    RooJSONFactoryWSTool::testValidName(axname,false);
    ax["name"] << axname;
    // min/max are written below once nBins is known -- x->getMin()/getMax() is not a reliable
    // source: CMS_TH1obs_<channel> is combine's own dummy per-bin index, and its own declared
    // range does not always already equal [0, nBins] (concretely: a real 31-bin channel whose "x"
    // reports range [0, 1]). Disagreeing with nBins here silently scales every sample's integral
    // by nBins on reimport -- see CMSHistSum.cc's tryExportHistFactory for the same fix.

    auto &samples = elem["samples"].set_seq();

    const int nBins = activeBins(summands, sum->GetName());

    for (unsigned ip = 0; ip < nProcs; ++ip) {
        auto * fvih = dynamic_cast<const FastVerticalInterpHistPdf2 *>(summands.at(ip));
        const std::string sampleName(cms_fast_vertinterp_hist_pdf_process_name(fvih, ip));
        auto &s = samples.append_child().set_map();
        s["name"] << sampleName;

        const auto nominal = fvih->cacheNominal();

        auto &dataNode = s["data"].set_map();
        double yieldScale = 1.;
        double shapeSum = 1.;
        if (int(nominal.size()) > nBins) {
            // unmorphed, with active bins beyond the axis (activeBins() declines the other case)
            auto const head = croppedHead(nominal, nBins, yieldScale);
            RooJSONFactoryWSTool::exportArray(nBins, head.data(), dataNode["contents"]);
        } else {
            // Built from histograms in this process: write the input nominal back unchanged.
            auto const &input = fvih->inputNominal();
            auto const &contents = input.fullsize() == unsigned(nBins) ? input.Values() : nominal.Values();
            RooJSONFactoryWSTool::exportArray(nBins, contents.data(), dataNode["contents"]);
            shapeSum = sampleSum(contents.data(), nBins);
        }

        auto &mods = s["modifiers"].set_seq();

        // setActiveBins crops a nominal with negative bins to 1e-9 without renormalizing it, so it
        // sums to S > 1, its morphs stay relative to the uncropped shape (hi/lo below sum to S as
        // well), and syncTotal() applies the 1e-9 floor on that scale before it normalizes. The
        // templates are written on that scale, with a literal normfactor 1/S: the histfactory_dist
        // is then exactly coefficient * (nominal + morph) / S, and the importer rebuilds the same
        // object (keepInputScale()), floor included.
        if (std::abs(shapeSum - 1.) > 1e-6) {
            auto &m = mods.append_child().set_map();
            m["name"] << sampleName + "_shape_norm";
            m["parameter"] << 1. / shapeSum;
            m["type"] << "normfactor";
        }

        ctx.sampleName = sampleName;
        ctx.sampleConstraints.clear();

        // Shape part of every nuisance in the coefList.

        for(unsigned icoeff=0; icoeff < fvih->coefList().size(); ++icoeff){
            RooAbsReal * coeff = dynamic_cast<RooAbsReal * >(fvih->coefList().at(icoeff));
            requirePlainModifierParameter(*fvih->coefList().at(icoeff), fvih->GetName(), sampleName, "histosys");
            const std::string constraint = resolveConstraintPdfName(tool, ctx, *coeff);
            registerConstraint(ctx, coeff->GetName(), constraint);
            auto &m = mods.append_child().set_map();
            m["name"] << coeff->GetName();
            m["parameter"] << coeff->GetName();
            m["type"] << "histosys";
            if (!constraint.empty()) m["constraint"] << constraint;
            auto & data = m["data"].set_map();
            FastTemplate morph_sum = fvih->morphs().at(icoeff).sum;
            FastTemplate morph_diff = fvih->morphs().at(icoeff).diff;
            FastTemplate hi(nBins), lo(nBins);
            if (fvih->smoothAlgo() < 0) {
                RooJSONFactoryWSTool::error("smoothalog <0 not implemented yet.");
            }
            // Built from histograms in this process (e.g. by the importer below): write the input
            // templates back unchanged. Reimporting them normalizes them the same way again.
            auto const &inputMorphs = fvih->inputMorphs();
            const bool haveInput = icoeff < inputMorphs.size() && inputMorphs[icoeff].first.fullsize() == unsigned(nBins);
            for (int ib = 0; haveInput && ib<nBins; ++ib){
                hi[ib] = inputMorphs[icoeff].first[ib];
                lo[ib] = inputMorphs[icoeff].second[ib];
            }
            for (int ib = 0; !haveInput && ib<nBins; ++ib){
                hi[ib] = (morph_sum[ib]+morph_diff[ib])/2.;
                lo[ib] = (morph_sum[ib]-morph_diff[ib])/2.;

                if(fvih->smoothAlgo() < 0){
                    hi[ib] = exp(hi[ib])*nominal[ib];
                    lo[ib] = exp(lo[ib])*nominal[ib];
                }else{
                    hi[ib] = hi[ib]+nominal[ib];
                    lo[ib] = lo[ib]+nominal[ib];
                }
            }
            RooJSONFactoryWSTool::exportArray(nBins, hi.Values().data(), data["hi"].set_map()["contents"]);
            RooJSONFactoryWSTool::exportArray(nBins, lo.Values().data(), data["lo"].set_map()["contents"]);
            tool->queueExport(*coeff);
        }
        // Rate part. decomposeCoefficient() already validated this in the decline pass above.
        auto const &parts = coefParts[ip];
        if (auto const *pn = parts.pn) {

            
            auto const &logKappa = pn->logKappa();
            for (int i = 0; i < static_cast<int>(logKappa.size()); ++i) {
                exportNormSys(tool, ctx, mods, pn->thetaList()[i], std::exp(-logKappa[i]), std::exp(logKappa[i]));
            }
            auto const &logAsymmKappa = pn->logAsymmKappa();
            for (int i = 0; i < static_cast<int>(logAsymmKappa.size()); ++i) {
                exportNormSys(tool, ctx, mods, pn->asymmThetaList()[i], std::exp(logAsymmKappa[i].first),
                              std::exp(logAsymmKappa[i].second));
            }

            const double nom = pn->nominalValue() * yieldScale;
            if (nom != 1.) {
                auto &m = mods.append_child().set_map();
                m["name"] << sampleName + "_norm_nominal";
                m["parameter"] << nom;
                m["type"] << "normfactor";
            }
            for (auto const *of : pn->otherFactorList()) {
                exportOtherFactorModifiers(tool, mods, static_cast<RooAbsReal const *>(of));
            }
        } else if (yieldScale != 1.) {
            auto &m = mods.append_child().set_map();
            m["name"] << sampleName + "_norm_nominal";
            m["parameter"] << yieldScale;
            m["type"] << "normfactor";
        }


        for (auto const *ap : parts.asymPows) {
            exportNormSys(tool, ctx, mods, ap->theta(), ap->kappaLow().getVal(), ap->kappaHigh().getVal());
        }

        for (auto const *pf : parts.plainFactors) {
            exportOtherFactorModifiers(tool, mods, pf);
        }

        if (parts.customCoeff) {
            auto &m = mods.append_child().set_map();
            m["name"] << parts.customCoeff->GetName();
            m["type"] << "custom";
            tool->queueExport(*parts.customCoeff);
        }
    }
    

    ax["nbins"] << nBins;
    ax["min"] << 0.0;
    ax["max"] << static_cast<double>(nBins);

    return true;
}

class FastVerticalInterpHistPdf2Streamer : public RooFit::JSONIO::Exporter {
public:

    std::string const & key() const override {
        static const std::string keystring("");
        return keystring;
    }
    bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override {
        // FastVerticalInterpHistPdf2 is exported as part of histfactory_dist below.
        return true;
    }
};


class SumOfFastVerticalInterpHistPdf2Streamer_RooAddition : public RooFit::JSONIO::Exporter {
public:
    bool autoExportDependants() const override {return false;}
    std::string const & key() const override { 
        static const std::string keystring("histfactory_dist");
        return keystring;
    }

    bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override {
        auto sum = static_cast<const RooAddPdf *>(p);
        // std::map<int, std::vector<RooAbsPdf>> sub_sums;
        // unsigned n_pdf = sum->pdfList().size();
        // for(unsigned ipdf=0; ipdf<n_pdf;++ipdf ){
        //     auto *sub_pdf = dynamic_cast<const FastVerticalInterpHistPdf2 *>(sum->pdfList().at(ipdf))
        //     if (!sub_pdf){}
        // }
        // TODO: split AddPdf into sub-sums with same binning
        // TODO: loop ober sub-subsums and export.
        if(sum) {
            bool retval = tryExportHistFactory(tool, sum, elem);
            if(!retval) return false;
            elem["type"] << key();
            std::cout << "exporting Sum of FastVerticalInterpHistPdf2s" << std::endl;

            return true;
        }
        return false;
    }

};

// The literal normfactor 1/S the exporter writes for a sample whose nominal sums to S != 1 (a nominal
// cropped by setActiveBins, see tryExportHistFactory), or nullptr.
RooFit::Detail::JSONNode const *shapeNormModifier(RooFit::Detail::JSONNode const &sample)
{
    auto const *mods = sample.find("modifiers");
    if (!mods) return nullptr;
    const std::string name = sample["name"].val() + "_shape_norm";
    for (auto const &m : mods->children()) {
        if (m["type"].val() == "normfactor" && m["name"].val() == name) return &m;
    }
    return nullptr;
}

// True if, for every sample, the nominal contents and every histosys hi/lo template sum to 1, or
// all to the same S with the sample's "<sample>_shape_norm" normfactor equal to 1/S. The tolerance
// only has to separate normalized shapes from templates in events; this file's exporter can be off
// by the 1e-9 crop floor in empty bins (3e-8 for part4's 31-bin templates).
bool allTemplatesNormalized(RooFit::Detail::JSONNode const &elem)
{
    auto sumOf = [](RooFit::Detail::JSONNode const &contents) {
        std::vector<double> v;
        for (auto const &c : contents.children()) v.push_back(c.val_double());
        return sampleSum(v.data(), v.size());
    };
    auto const *samples = elem.find("samples");
    if (!samples) return false;
    for (auto const &s : samples->children()) {
        auto const *contents = s.find("data", "contents");
        if (!contents) return false;
        double norm = 1.;
        if (auto const *m = shapeNormModifier(s)) {
            double v = 0.;
            if (!parse_double((*m)["parameter"].val(), v)) return false;
            norm = sumOf(*contents);
            if (std::abs((*m)["parameter"].val_double() * norm - 1.) > 1e-12) return false;
        }
        auto normalized = [&](RooFit::Detail::JSONNode const &c) { return std::abs(sumOf(c) / norm - 1.) <= 1e-6; };
        if (!normalized(*contents)) return false;
        if (auto const *mods = s.find("modifiers")) {
            for (auto const &m : mods->children()) {
                if (m["type"].val() != "histosys") continue;
                auto const *hi = m.find("data", "hi", "contents");
                auto const *lo = m.find("data", "lo", "contents");
                if (!hi || !lo || !normalized(*hi) || !normalized(*lo)) return false;
            }
        }
    }
    return true;
}

class SumOfFastVerticalInterpHistPdf2Factory_RooAddition : public RooFit::JSONIO::Importer{
public:
    bool importArg(RooJSONFactoryWSTool *tool, RooFit::Detail::JSONNode const &elem) const override
    {
        using RooFit::Detail::JSONNode;

        if (!elem.has_child("name")) return false;
        const std::string distName = elem["name"].val();

        // Single-bin counting nodes belong to the RooPoisson/RooAddition importer.
        if (isSingleBinCountingDist(elem)) return false;

        // FastVerticalInterpHistPdf2 is a normalized shape with a separate yield, so it only
        // equals the histfactory_dist if every template already sums to 1 (as this file's exporter
        // writes them). Templates in events (e.g. from CMSHistErrorPropagator) go to the CMSHistSum
        // importer, which keeps them as they are.
        if (!allTemplatesNormalized(elem)) return false;

        // Axes: we only support 1D here
        if (!elem.has_child("axes") || elem["axes"].num_children() != 1)
            throw std::runtime_error("Sum of FastVerticalInterpHistPdf2 importer: '" + distName + "' must have exactly one axis");

        RooArgSet observables = RooJSONFactoryWSTool::readAxes(elem);
        RooRealVar* x = dynamic_cast<RooRealVar*>(observables.first());

        // Ensure ONE exists (fixed at 1) in case ProcessNormalization has not a norm of its own.
        RooAbsReal* ONE = dynamic_cast<RooAbsReal*>(tool->workspace()->arg("ONE"));
        if (!ONE) {
            auto* one = new RooConstVar("ONE", "ONE", 1.0);
            tool->workspace()->import(*one, RooFit::RecycleConflictNodes());
            ONE = dynamic_cast<RooAbsReal*>(tool->workspace()->arg("ONE"));
        }

        // Samples
        if (!elem.has_child("samples") || elem["samples"].num_children() == 0)
            throw std::runtime_error("Sum of FastVerticalInterpHistPdf2  importer: '" + distName + "' has no samples");


        // structure to construct: RooAddPdf(summands=FastVerticalPdf(hist;histosys-params), coeff=ProcessNorm(normsys, normfactor))
        RooArgList summand_list;
        RooArgList coeff_list;
        for(const auto &s : elem["samples"].children()){
            RooArgList histosys_coeff;
            TList histosys_funcs;
            const std::string sampleName = s["name"].val();
            // data.contents (required)
            if (!s.has_child("data") || !s["data"].has_child("contents")) {
                throw std::runtime_error("Sum of FastVerticalInterpHistPdf2 importer: sample '" + sampleName + "' missing data.contents");
            }

            const std::string nom_hname = distName + "_" + sampleName + "_nominal";
            auto h = readTH1D(s["data"]["contents"], *x, nom_hname);
            // adding nominal template
            histosys_funcs.Add(h.release());

            RooArgList normsys_pars;
            RooArgList normfactor_pars;
            std::vector<std::pair<double, double>> normsys_vals;
            // Parameters carrying a histosys on this sample. A normsys sharing one of these names is
            // the yield half of a `shape` systematic and gets rebuilt as a standalone AsymPow below.
            std::set<std::string> histosys_par_names;
            double norm_nom = 1.;
            const bool keepScale = shapeNormModifier(s) != nullptr;

            // The exporter writes a sample's rate part as: ProcessNormalization normsys, then
            // "<sample>_norm_nominal" and the other normfactor/custom factors, then the normsys of
            // standalone AsymPows (the yield half of a `shape` systematic). A normsys that comes
            // after a factor can therefore only be such an AsymPow, and one before every factor is a
            // ProcessNormalization entry, even if its parameter also carries a histosys (text2workspace
            // may put it there). Without any factor modifier the two cannot be told apart, and the
            // parameter name alone decides, as before.
            bool sampleHasFactor = false;
            if (s.has_child("modifiers")) {
                for (const auto &m : s["modifiers"].children()) {
                    const std::string mtype = m["type"].val();
                    if ((mtype == "normfactor" && m["name"].val() != sampleName + "_shape_norm") || mtype == "custom")
                        sampleHasFactor = true;
                }
            }
            bool factorSeen = false;
            std::vector<bool> normsys_after_factor;

            if (s.has_child("modifiers")) {
                for (const auto& m: s["modifiers"].children()) {
                    const std::string mtype = m["type"].val();
                    const std::string mname = m["name"].val();
                    if (mtype == "custom" || (mtype == "normfactor" && !(keepScale && mname == sampleName + "_shape_norm")))
                        factorSeen = true;

                    if (keepScale && mtype == "normfactor" && mname == sampleName + "_shape_norm") {
                        // 1/S of a nominal kept on its own scale: the FastVerticalInterpHistPdf2
                        // normalizes itself (keepInputScale() below), so it is not part of the yield.
                    } else if (mtype == "normfactor"){
                        double numeric = 1.;
                        if(parse_double(m["parameter"].val(), numeric)){
                            // A normfactor whose "parameter" is a literal number is a plain
                            // multiplicative constant. Several can occur on one sample (the process
                            // nominal plus any RooConstVar factor), so they accumulate rather than
                            // overwrite each other. val() only has 6 significant digits, so the
                            // value itself is read with val_double().
                            norm_nom *= m["parameter"].val_double();
                        } else {
                            RooRealVar &p = getOrMakeParameter(tool, mname, 1.0, 0.0, 10.0);
                            normfactor_pars.add(p);
                        }
                    }else if (mtype == "custom") {

                        // This might be normalization factor that is given by a ProcessNormalization or RooFormulaVar
                        RooAbsReal* obj = tool->workspace()->function(mname);
                        if (!obj) {
                        throw std::runtime_error("Sum of FastVerticalInterpHistPdf2 importer: custom modifier '" + mname + "' on sample '" +
                                                sampleName +
                                                "' not yet imported, there seems to be a priority/ordering issue");
                        }
                        normfactor_pars.add(*obj);

                    } else if (mtype == "staterror") {
                        return false;
                        RooJSONFactoryWSTool::warning("SumOfFastVerticalInterpHistPdf2 importer: sample '" + sampleName + "' has a staterror modifier " + mname);
                        // throw std::runtime_error("Sum of FastVerticalInterpHistPdf2 importer: sample '" + sampleName + "' has a staterror modifier " + mname);
                    } else if ( mtype == "normsys"){

                        auto * par = m.find("parameter");
                        std::string parname(par ? par->val() : "alpha_"+mname);
                        RooRealVar &p = getOrMakeParameter(tool, parname, 0.0, -10.0, 10.0);
                        auto & data = m["data"];
                        double hi = data["hi"].val_double();
                        double lo = data["lo"].val_double();
                        normsys_pars.add(p);
                        normsys_vals.push_back(std::make_pair(lo,hi));
                        normsys_after_factor.push_back(factorSeen);

                    } else if ( mtype == "histosys" ){

                        auto * par = m.find("parameter");
                        std::string parname(par ? par->val() : "alpha_"+mname);

                        RooRealVar& p = getOrMakeParameter(tool, parname, 0.0, -10., 10.0);

                        auto &data = m["data"];
                        auto low = readTH1D(data["lo"]["contents"], *x, mname+"Low_"+sampleName);
                        auto high = readTH1D(data["hi"]["contents"], *x, mname+"High_"+sampleName);
                        // adding first hi, then lo, since FastVerticalInterpHistPdf2::FastVerticalInterpHistPdf2 dictates this order.
                        histosys_funcs.Add(high.release());
                        histosys_funcs.Add(low.release());
                        histosys_coeff.add(p);
                        histosys_par_names.insert(parname);

                    } else {
                    throw std::runtime_error("Sum of FastVerticalInterpHistPdf2 importer: modifier type '" + mtype + "' on sample '" + sampleName +
                                            "' is not supported by this Sum of FastVerticalInterpHistPdf2 importer yet");
                    }

                } // end modifiers loop
            } // end modifiers safe
            std::string summand_name = distName + "_" +sampleName; 
            FastVerticalInterpHistPdf2 * fvih = new FastVerticalInterpHistPdf2(summand_name.c_str(),summand_name.c_str(), *x, histosys_funcs, histosys_coeff);
            if (keepScale) fvih->keepInputScale();
            fvih->setStringAttribute("combine.process", sampleName.c_str());
            summand_list.add(*fvih);

            // The yield is the ProcessNormalization's own nominal value, as text2workspace.py builds
            // it; the exporter writes that back as the "<sample>_norm_nominal" normfactor.
            const std::string pnName = distName + "_processnorm_" +sampleName;
            auto * pn = new ProcessNormalization(pnName.c_str(), pnName.c_str(), norm_nom);

            // A normsys whose parameter also carries a histosys on this sample is the yield half of
            // a `shape` systematic: rebuild ShapeTools.py's standalone AsymPow for it and multiply it
            // alongside the ProcessNormalization, so the reimported graph matches what
            // text2workspace.py produces. Everything else stays an ordinary (asymm)LogNormal entry.
            //
            // The two encodings are numerically identical -- MathFuncs::asymPow(theta, kLo, kHi) is
            // exactly the asymmThetaList term inside MathFuncs::processNormalization -- so this is
            // about structural fidelity, not about the value.
            RooArgList asymPows;
            std::set<std::string> promoted;
            for(unsigned i = 0; i < normsys_pars.size(); ++i){
                auto par = dynamic_cast<RooAbsReal*>(&normsys_pars[i]);
                if (!par) RooJSONFactoryWSTool::error("Sum of FastVerticalInterpHistPdf2 importer: sample '" +
                                                      sampleName + "': normsys parameter is not a RooAbsReal");
                const std::string parname(par->GetName());

                // Promote at most one normsys per parameter; a parameter carrying both a `shape` and
                // an `lnN` line legitimately produces two, and the remainder belong in the
                // ProcessNormalization. The factors commute, so the split is faithful either way.
                const bool asymPowPosition = !sampleHasFactor || normsys_after_factor[i];
                if (asymPowPosition && histosys_par_names.count(parname) && promoted.insert(parname).second) {
                    const std::string apName = distName + "_asympow_" + sampleName + "_" + parname;
                    auto *kLo = new RooConstVar((apName + "_kappaLo").c_str(), "", normsys_vals[i].first);
                    auto *kHi = new RooConstVar((apName + "_kappaHi").c_str(), "", normsys_vals[i].second);
                    asymPows.add(*new AsymPow(apName.c_str(), apName.c_str(), *kLo, *kHi, *par));
                    continue;
                }

                if(approxSymmetricKappas(normsys_vals[i].first, normsys_vals[i].second)) {
                    pn->addLogNormal(normsys_vals[i].second, *par);
                } else {
                    pn->addAsymmLogNormal(normsys_vals[i].first, normsys_vals[i].second, *par);
                }
            }

            // Export flattens RooProduct siblings and ProcessNormalization::otherFactorList entries
            // into the same "normfactor" modifier, so they cannot be told apart here. Folding them
            // all into the ProcessNormalization is numerically identical to either arrangement.
            for(unsigned i = 0; i < normfactor_pars.size(); ++i){
                auto par = dynamic_cast<RooAbsReal*>(&normfactor_pars[i]);
                if (!par) RooJSONFactoryWSTool::error("Sum of FastVerticalInterpHistPdf2 importer: sample '" +
                                                      sampleName + "': normfactor is not a RooAbsReal");
                pn->addOtherFactor(*par);
            }

            if (asymPows.empty()) {
                coeff_list.add(*pn);
            } else {
                const std::string prodName = distName + "_norm_final_" + sampleName;
                RooArgList factors(*pn);
                factors.add(asymPows);
                coeff_list.add(*new RooProduct(prodName.c_str(), prodName.c_str(), factors));
            }

        } // end samples

        // The channel itself; wrapped below into RooProdPdf(distName) with its modifiers' constraints, if any.
        const std::string coreName = CombineHS3::channelCoreName(elem);
        RooAddPdf * sum = new RooAddPdf(coreName.c_str(), distName.c_str(), summand_list, coeff_list);
        tool->workspace()->import(*sum, RooFit::RecycleConflictNodes());
        sum = dynamic_cast<RooAddPdf*>(tool->workspace()->pdf(coreName.c_str()));
        if (!sum) throw std::runtime_error("Sum of FastVerticalInterpHistPdf2 importer: failed to import Sum of FastVerticalInterpHistPdf2 '" + coreName + "'");
        CombineHS3::wrapChannelInConstraints(tool, elem);
        return true;
        
    } //end import
};

#include "static_execute.h"

STATIC_EXECUTE( []() {
    using namespace RooFit::JSONIO;

    registerExporter<FastVerticalInterpHistPdf2Streamer>("FastVerticalInterpHistPdf2",false);
    registerExporter<SumOfFastVerticalInterpHistPdf2Streamer_RooAddition>("RooAddPdf", true);
    registerImporter<SumOfFastVerticalInterpHistPdf2Factory_RooAddition>("histfactory_dist",true);
});

}
