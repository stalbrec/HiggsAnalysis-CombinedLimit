#ifndef HiggsAnalysis_CombinedLimit_CombineHS3Helpers_h
#define HiggsAnalysis_CombinedLimit_CombineHS3Helpers_h

// Internal helpers shared by the HS3 JSON exporters and importers (JSONFactories_Combine.cxx,
// CMSHistSumHS3.cc, VerticalInterpHistPdfHS3.cc). Not installed, not part of the dictionary.

#include <RooFitHS3/JSONIO.h>
#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include <RooFit/Detail/JSONInterface.h>

#include "../interface/CombineHS3Channel.h"

#include "RooAbsArg.h"
#include "RooAbsCollection.h"
#include "RooAbsPdf.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace CombineHS3 {

/// Parse a whole string as a double; false (and `out` untouched) unless every character is consumed.
inline bool parse_double(const std::string &input, double &out)
{
   const char *begin = input.data();
   const char *end = input.data() + input.size();

   double tmp;
   auto [ptr, ec] = std::from_chars(begin, end, tmp);

   if (ec == std::errc() && ptr == end) {
      out = tmp;
      return true;
   }
   return false;
}

/// True if kLo ~ 1/kHi, i.e. log(kLo) ~ -log(kHi), within the given absolute/relative log tolerance.
inline bool approxSymmetricKappas(double kLo, double kHi, double absLogTol = 1e-10, double relLogTol = 1e-10)
{
   // kappas must be positive and finite
   if (!(kLo > 0.0) || !(kHi > 0.0))
      return false;
   if (!std::isfinite(kLo) || !std::isfinite(kHi))
      return false;

   const double a = std::log(kLo);
   const double b = std::log(kHi);

   // symmetric if a ~ -b  <=>  (a + b) ~ 0
   const double diff = std::abs(a + b);
   const double scale = std::max(std::abs(a), std::abs(b)); // 0 if both ~1
   return diff <= std::max(absLogTol, relLogTol * scale);
}

/// Look up the RooRealVar `parName` in the tool's workspace, or create and import it.
inline RooRealVar &getOrMakeParameter(RooJSONFactoryWSTool *tool, std::string const &parName, double init = 1.0,
                                      double lo = 0.0, double hi = 10.0)
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

/// An HS3 histosys/normsys "parameter" is the name of a scalar free parameter
/// (histfactory_dist in the HS3 standard). A coefficient that is not a plain RooRealVar
/// (e.g. RooProduct(alpha, c)) has no exact representation there, so the export stops instead of
/// writing the compound node's name.
inline void requirePlainModifierParameter(RooAbsArg const &par, std::string const &pdfName,
                                          std::string const &sampleName, const char *modType)
{
   if (dynamic_cast<RooRealVar const *>(&par))
      return;
   RooJSONFactoryWSTool::error(std::string("HS3 export of pdf '") + pdfName + "', sample '" + sampleName + "': " +
                               modType + " coefficient '" + par.GetName() + "' is a " + par.ClassName() +
                               ", not a plain RooRealVar; an HS3 modifier \"parameter\" must be a scalar free "
                               "parameter (see OPEN_ISSUES.md, scaled normsys/histosys coefficients).");
}

/// Name of the constraint pdf of `theta` to write as a modifier's "constraint" ("" = none), cached per theta.
/// The constraint is the workspace's "<theta>_Pdf" (the name ShapeTools.py gives every constraint pdf),
/// else the first pdf that depends on theta but on none of the observables. It is named only if no
/// product enclosing `channel` applies it already (CombineHS3Channel.h), and queued for export either way.
inline std::string resolveConstraintPdfName(RooJSONFactoryWSTool *tool, RooAbsArg const &theta,
                                            RooAbsCollection const &observables, RooAbsArg const *channel,
                                            std::map<std::string, std::string> &cache)
{
   const std::string thetaName(theta.GetName());
   auto cached = cache.find(thetaName);
   if (cached != cache.end())
      return cached->second;

   RooWorkspace *ws = tool->workspace();
   RooAbsPdf *constraint = nullptr;

   if (auto *byName = ws->pdf((thetaName + "_Pdf").c_str())) {
      if (byName->dependsOn(theta))
         constraint = byName;
   }

   if (!constraint) {
      for (auto *arg : ws->allPdfs()) {
         auto *pdf = dynamic_cast<RooAbsPdf *>(arg);
         if (!pdf || !pdf->dependsOn(theta))
            continue;
         if (pdf->dependsOn(observables))
            continue;
         constraint = pdf;
         break;
      }
   }

   std::string name;
   if (constraint) {
      if (!channel || !constraintAppliedOutside(*channel, *constraint))
         name = constraint->GetName();
      // The referenced name has to exist in the output, and these exporters do not
      // autoExportDependants.
      tool->queueExport(*constraint);
   }
   cache[thetaName] = name;
   return name;
}

/// Append a normsys modifier on `theta` with the given kappas and constraint name ("" = none) to `mods`,
/// and queue theta for export.
inline void writeNormSys(RooJSONFactoryWSTool *tool, RooFit::Detail::JSONNode &mods, RooAbsArg const &theta,
                         double kappaLo, double kappaHi, std::string const &constraint)
{
   auto &m = mods.append_child().set_map();
   m["name"] << theta.GetName();
   m["parameter"] << theta.GetName();
   m["type"] << "normsys";
   if (!constraint.empty())
      m["constraint"] << constraint;
   auto &data = m["data"].set_map();
   data["lo"] << kappaLo;
   data["hi"] << kappaHi;
   tool->queueExport(theta);
}

/// The values of a JSON sequence of numbers.
inline std::vector<double> readDoubleArray(RooFit::Detail::JSONNode const &arr)
{
   std::vector<double> out;
   out.reserve(arr.num_children());
   for (const auto &c : arr.children()) {
      out.push_back(c.val_double());
   }
   return out;
}

/// Dispatch predicate between the single-bin counting importer (JSONFactories_Combine.cxx) and the
/// shape importers: true only for nodes that are unambiguously counting-shaped, i.e. a single 1-bin
/// axis, one content value per sample and only rate modifiers.
inline bool isSingleBinCountingDist(RooFit::Detail::JSONNode const &elem)
{
   auto const *axes = elem.find("axes");
   if (!axes || axes->num_children() != 1)
      return false;
   for (auto const &ax : axes->children()) {
      if (ax.has_child("edges")) {
         if (ax["edges"].num_children() != 2) // two edges == one bin
            return false;
      } else {
         auto const *nb = ax.find("nbins");
         if (!nb || nb->val_int() != 1)
            return false;
         // the count axis n_obs_bin<b> is unbounded above (written as null); a 1-bin shape
         // channel (CMS_TH1obs_<channel>, [0, 1]) is not
         auto const *mx = ax.find("max");
         if (mx && !mx->is_null())
            return false;
      }
   }

   auto const *samples = elem.find("samples");
   if (!samples || samples->num_children() == 0)
      return false;
   for (auto const &s : samples->children()) {
      auto const *data = s.find("data");
      if (!data || !data->has_child("contents") || (*data)["contents"].num_children() != 1)
         return false;
      if (auto const *mods = s.find("modifiers")) {
         for (auto const &m : mods->children()) {
            std::string const mtype = m["type"].val();
            if (mtype != "normsys" && mtype != "normfactor" && mtype != "custom")
               return false;
         }
      }
   }
   return true;
}

} // namespace CombineHS3

#endif
