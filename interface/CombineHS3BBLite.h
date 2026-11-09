#ifndef HiggsAnalysis_CombinedLimit_CombineHS3BBLite_h
#define HiggsAnalysis_CombinedLimit_CombineHS3BBLite_h

// HS3 export of autoMCStats (Barlow-Beeston-lite) bin parameters of CMSHistErrorPropagator and
// CMSHistSum (Decision 1 with point 6 (y_j) and point 7 option (a), DECISIONS.md 2026-11-03).
// This covers the bins with one parameter for the total (type 1). CMSHistSum's per-process bins
// (types 2, 3; DECISIONS.md 2026-11-04 hig20013) keep their parameters and are written by its
// histfactory exporter (CMSHistSum.cc, exportPerProcessBinSamples).
//
// Combine, bin type 1:  n_j = S_j(theta) + x_j * E_j(theta),  x_j ~ Gauss(x_j_In | x_j, 1),
//                       E_j(theta) = sqrt(sum_p c_p(theta)^2 e_pj^2).
// HS3:                  n_j = S_j(theta) + y_j  (a `mcstat` sample with contents 1 and an
//                       unconstrained `shapefactor` with parameters y_j), and the constraint
//                       gaussian_dist(x = x_j_In, mean = y_j / E_j(theta), sigma = 1).
// Bins without a term (type 0, or beyond the propagator's active bins) have contents 0 in `mcstat`
// and borrow a neighbouring y_j (shapefactor needs a parameter per bin).
// With y_j = x_j * E_j(theta) both are the same likelihood. The 1e-9 floor Combine applies to the
// bin totals is not written (point 7 option (a)).
//
// Objects written for bin parameter <x> (e.g. prop_binch1_bin3): the unbounded parameter <x>_y, the
// functions <x>_E (E_j) and <x>_mean (y_j / E_j), and per process p of the channel
// <prop>_norm_<p> (c_p, as a generic expression), unless c_p is a plain parameter or a literal.
// The constraint keeps its name (<x>_Pdf).
//
// A propagator that shares x_j with the s+b model's (the b-only model's) but has another E_j gets its own
// <x>_E_<tag>, <x>_mean_<tag> and constraint copies <x>_Pdf_<tag> (tag: "bonly" for prop_bin<ch>_bonly);
// its channel products point at those (and at copies <product>_<tag> of the products holding them)
// for the duration of the export. y_j is shared: in that model it stands for x_j * E_j^<tag>.

#include <map>
#include <string>
#include <vector>

class RooAbsArg;
class RooAbsPdf;
class RooAbsReal;
class RooArgSet;
class RooRealVar;
class RooWorkspace;
namespace RooStats {
class ModelConfig;
}

namespace CombineHS3 {

/// For the duration of an export: builds the objects above in `ws` for every
/// CMSHistErrorPropagator / CMSHistSum with type-1 bin parameters, and swaps x_j for y_j in every
/// ModelConfig's nuisance parameters (restored by the destructor; the new objects stay in the workspace).
/// Throws std::runtime_error for a CMSHistErrorPropagator with bin types 2/3 (not written).
class BBLiteExport {
public:
   explicit BBLiteExport(RooWorkspace &ws);
   ~BBLiteExport();
   BBLiteExport(BBLiteExport const &) = delete;
   BBLiteExport &operator=(BBLiteExport const &) = delete;

   /// The instance of the export in progress, or nullptr.
   static BBLiteExport const *active();

   /// y_j per bin (nullptr for bins without a type-1 term) of this CMSHistErrorPropagator / CMSHistSum, or nullptr if it has no reparameterized bins.
   std::vector<RooRealVar *> const *binParameters(RooAbsArg const &prop) const;
   /// The constraint of x_j per bin (as binParameters(); for a model with its own copies, the copy), or nullptr.
   std::vector<RooAbsPdf *> const *binConstraints(RooAbsArg const &prop) const;
   /// For a constraint Gauss(x_In | x_j, 1) of a bin parameter: its mean function y_j / E_j, or nullptr.
   RooAbsReal *meanFor(RooAbsArg const &constraint) const;

private:
   std::map<std::string, std::vector<RooRealVar *>> ys_; // propagator name -> y_j per bin
   std::map<std::string, std::vector<RooAbsPdf *>> cons_; // propagator name -> constraint of x_j per bin
   std::map<std::string, RooAbsReal *> means_;           // constraint name -> mean function
   struct Redirect {
      RooAbsArg *product, *original, *replacement;
   };
   std::vector<Redirect> redirects_; // channel factors replaced for the export
   std::vector<std::pair<RooStats::ModelConfig *, RooArgSet *>> savedNuisances_;
   std::vector<std::pair<RooRealVar *, double>> savedErrors_; // x_j errors, removed during the export
   BBLiteExport const *previous_;
};

} // namespace CombineHS3

#endif
