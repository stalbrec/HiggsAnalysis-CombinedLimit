#ifndef HiggsAnalysis_CombinedLimit_CombineHS3Input_h
#define HiggsAnalysis_CombinedLimit_CombineHS3Input_h

// Loading of HS3 JSON input files: `combine <file>.json` instead of a .root workspace.

#include <string>

class RooSimultaneous;
class RooSimultaneousOpt;
class RooWorkspace;
namespace RooStats {
class ModelConfig;
}

namespace CombineHS3 {

/// A new workspace "w" with the contents of the HS3 JSON file `fileName` (RooJSONFactoryWSTool::importJSON).
RooWorkspace *importWorkspace(std::string const &fileName);

/// The RooSimultaneousOpt for `pdf`, the model pdf of `mc` imported from `fileName`; replaces
/// `new RooSimultaneousOpt(pdf, ...)` for JSON input. Adds the categories to mc's observables.
///
/// The histfactory_dist importers build a channel whose modifiers name constraints as
/// RooProdPdf(<name>, [<name>_model, constraints]) (CombineHS3Channel.h). Combine's own models keep
/// the channels free of constraints and carry all of them as extra constraints of the
/// RooSimultaneousOpt; CachingSimNLL hands the parameters to the minimizer in the order it meets
/// the constraints, and Minuit's path depends on that order. So the channels go into the
/// RooSimultaneousOpt without their constraints, and all constraints of the model follow as extra
/// constraints in text2workspace's order: the datacard's in ModelConfig nuisance order, then the
/// autoMCStats ones channel by channel. Channel masks are read from the file's misc.CMS.channel_masks.
RooSimultaneousOpt *makeSimultaneousOpt(RooWorkspace &w, RooStats::ModelConfig &mc, RooSimultaneous &pdf,
                                        std::string const &fileName, int verbose);

} // namespace CombineHS3

#endif
