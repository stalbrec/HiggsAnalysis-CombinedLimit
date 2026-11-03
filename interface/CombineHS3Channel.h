#ifndef HiggsAnalysis_CombinedLimit_CombineHS3Channel_h
#define HiggsAnalysis_CombinedLimit_CombineHS3Channel_h

// HS3: where a channel's constraint terms live (hig20013_nobb, "RooSimultaneousOpt items", option A).
//
// A histfactory_dist modifier that names a `constraint` brings that constraint term into the
// distribution (/spec docs/parts/histfactory.md). So Combine's histfactory_dist importers build
// such a channel as RooProdPdf(<name>, [<name>_model, constraints in order of first reference]),
// with the channel itself under <name>_model, as ROOT's own HistFactory importer does.
// exportConstrainedChannel() writes that product back as the same histfactory_dist.
//
// The exporters name a modifier's constraint only if no enclosing RooProdPdf of the channel
// already applies it (e.g. text2workspace's RooProdPdf(channel, constraints) channels, or
// model_b = RooProdPdf(..., nuisancePdf)). Otherwise the JSON would apply it twice.

#include <string>
#include <vector>

class RooAbsArg;
class RooProdPdf;
class RooJSONFactoryWSTool;
namespace RooFit {
namespace Detail {
class JSONNode;
}
} // namespace RooFit

namespace CombineHS3 {

/// The `constraint` names of all modifiers of a histfactory_dist node, in order of first reference.
std::vector<std::string> modifierConstraintNames(RooFit::Detail::JSONNode const &elem);

/// Name under which a histfactory_dist importer builds the channel itself: the distribution name,
/// or "<name>_model" if wrapChannelInConstraints() is going to wrap it.
std::string channelCoreName(RooFit::Detail::JSONNode const &elem);

/// After the channel was imported as channelCoreName(elem): build the RooProdPdf <name> of it and
/// its modifiers' constraints. No-op if no modifier names a constraint.
void wrapChannelInConstraints(RooJSONFactoryWSTool *tool, RooFit::Detail::JSONNode const &elem);

/// Exporter side: true if `constraint` is a factor of a RooProdPdf enclosing `channel` (through
/// RooProdPdf and RooSimultaneous clients), so the channel's modifiers must not name it.
bool constraintAppliedOutside(RooAbsArg const &channel, RooAbsArg const &constraint);

/// Export RooProdPdf(<name>, [<name>_model, constraints]) as the histfactory_dist of <name>_model,
/// if the core's exporter names exactly these constraints, in this order. False otherwise.
bool exportConstrainedChannel(RooJSONFactoryWSTool *tool, RooProdPdf const &prod, RooFit::Detail::JSONNode &elem);

} // namespace CombineHS3

#endif
