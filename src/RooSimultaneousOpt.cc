#include "../interface/RooSimultaneousOpt.h"
#include "../interface/CachingNLL.h"
#include "../interface/Combine.h"
#include "../interface/CombineHS3BBLite.h"

#include <RooCmdConfig.h>
#include <RooWorkspace.h>
#include <RooFit/Detail/JSONInterface.h>
#include <RooFitHS3/RooJSONFactoryWSTool.h>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

#if ROOT_VERSION_CODE < ROOT_VERSION(6,30,0)
RooAbsReal* 
RooSimultaneousOpt::createNLL(RooAbsData& data, const RooLinkedList& cmdList) 
#else
std::unique_ptr<RooAbsReal>
RooSimultaneousOpt::createNLLImpl(RooAbsData& data, const RooLinkedList& cmdList) 
#endif
{
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 36, 0)
  // Match the printout of https://github.com/root-project/root/blob/a6ec2bd88e82f46f131cece8b3cc6b357c756f02/roofit/roofitcore/src/FitHelpers.cxx#L594
  auto timingScope = std::make_unique<ROOT::Math::Util::TimingScope>(
      [this](std::string const& msg) { oocoutI(this, Fitting) << msg << std::endl; }, "Creation of NLL object took");
#endif
  RooCmdConfig pc(Form("RooSimultaneousOpt::createNLL(%s)", GetName()));
  pc.defineSet("cPars", "Constrain", 0, 0);
  RooArgSet* cPars = pc.getSet("cPars");
  auto nll = std::make_unique<cacheutils::CachingSimNLL>(this, &data, cPars);
  nll->setChannelMasks(this->channelMasks());
#if ROOT_VERSION_CODE < ROOT_VERSION(6,30,0)
    return nll.release();
#else
    return nll;
#endif
}

void
RooSimultaneousOpt::addExtraConstraint(const RooAbsPdf &pdf)
{
    if (!_extraConstraints.contains(pdf)) _extraConstraints.add(pdf);
}

void
RooSimultaneousOpt::addExtraConstraints(const RooAbsCollection &pdfs)
{
    if (_extraConstraints.getSize() == 0) {
        // fast (hopefully)
        _extraConstraints.add(pdfs);
    } else {
        // slow
        for (RooAbsArg *a : pdfs) {
            if (!_extraConstraints.contains(*a)) _extraConstraints.add(*a);
        }
    }
}

void
RooSimultaneousOpt::addChannelMasks(const RooArgList &args)
{
    if (args.getSize() != indexCat().numTypes()) {
        std::cerr << "RooSimultaneousOpt: number of channel masks must equal number of channels\n";
        return;
    }
    _channelMasks.removeAll();
    _channelMasks.add(args);
}

namespace {
// CachingSimNLL pairs mask i with category bin i (setBin(i)), whatever the state indices are.
std::map<std::string, std::string> channelMaskNames(const RooSimultaneousOpt &pdf)
{
    std::map<std::string, std::string> out;
    std::unique_ptr<RooAbsCategoryLValue> cat{static_cast<RooAbsCategoryLValue *>(pdf.indexCat().Clone())};
    for (int i = 0, n = pdf.channelMasks().getSize(); i < n; ++i) {
        cat->setBin(i);
        out[cat->getLabel()] = pdf.channelMasks().at(i)->GetName();
    }
    return out;
}
} // namespace

bool
RooSimultaneousOpt::exportJSON(RooWorkspace &ws, const std::string &fileName)
{
    std::map<std::string, std::string> masks;
    for (RooAbsArg *arg : ws.components()) {
        auto *opt = dynamic_cast<RooSimultaneousOpt *>(arg);
        if (!opt || opt->channelMasks().getSize() == 0) continue;
        for (auto const &[channel, mask] : channelMaskNames(*opt)) {
            auto [it, inserted] = masks.emplace(channel, mask);
            if (!inserted && it->second != mask) {
                std::cerr << "RooSimultaneousOpt::exportJSON: channel " << channel << " has masks " << it->second
                          << " and " << mask << "\n";
                return false;
            }
        }
    }
    // The extra constraints are exported where the channels' modifiers name them (CombineHS3Channel.h),
    // which is also where a re-export of the imported JSON finds them: Combine's histfactory_dist
    // importers put them into the channels, and a plain import builds a bare RooSimultaneous. Left
    // attached, ROOT would export them first, as servers of the RooSimultaneousOpt, and the order of
    // distributions, parameters and domains would differ from a re-export. Detached only for the
    // export; unreferenced ones still go out (other clients, or as top-level pdfs).
    struct Detached {
        RooSimultaneousOpt *opt;
        RooArgList constraints;
    };
    std::vector<Detached> detached;
    for (RooAbsArg *arg : ws.components()) {
        auto *opt = dynamic_cast<RooSimultaneousOpt *>(arg);
        if (!opt || opt->_extraConstraints.getSize() == 0) continue;
        detached.push_back({opt, RooArgList(opt->_extraConstraints)});
        opt->_extraConstraints.removeAll();
    }
    struct Restore {
        std::vector<Detached> &d;
        ~Restore() { for (auto &x : d) x.opt->_extraConstraints.add(x.constraints); }
    } restore{detached};

    // autoMCStats bin parameters as y_j (CombineHS3BBLite.h), for the duration of the export
    std::unique_ptr<CombineHS3::BBLiteExport> bblite;
    try {
        bblite = std::make_unique<CombineHS3::BBLiteExport>(ws);
    } catch (std::exception const &e) {
        std::cerr << "RooSimultaneousOpt::exportJSON: " << e.what() << "\n";
        return false;
    }

    RooJSONFactoryWSTool tool(ws);
    if (masks.empty()) return tool.exportJSON(fileName);

    std::unique_ptr<RooFit::Detail::JSONTree> tree = RooFit::Detail::JSONTree::create(tool.exportJSONtoString());
    auto &node = tree->rootnode().get("misc", "CMS", "channel_masks").set_map();
    for (auto const &[channel, mask] : masks) node[channel] << mask;
    std::ofstream out(fileName);
    if (!out) return false;
    tree->rootnode().writeJSON(out);
    return bool(out);
}

bool
RooSimultaneousOpt::importChannelMasksJSON(const std::string &fileName, RooWorkspace &ws)
{
    std::ifstream in(fileName);
    if (!in) return false;
    std::unique_ptr<RooFit::Detail::JSONTree> tree = RooFit::Detail::JSONTree::create(in);
    const RooFit::Detail::JSONNode *node = tree->rootnode().find("misc", "CMS", "channel_masks");
    if (!node) return true;
    if (!node->is_map()) return false;
    std::map<std::string, std::string> masks;
    for (auto const &child : node->children()) masks[child.key()] = child.val();

    RooArgList list;
    std::unique_ptr<RooAbsCategoryLValue> cat{static_cast<RooAbsCategoryLValue *>(indexCat().Clone())};
    for (int i = 0, n = cat->numBins(nullptr); i < n; ++i) {
        cat->setBin(i);
        auto it = masks.find(cat->getLabel());
        RooAbsReal *mask = it == masks.end() ? nullptr : dynamic_cast<RooAbsReal *>(ws.arg(it->second.c_str()));
        if (!mask) {
            std::cerr << "RooSimultaneousOpt::importChannelMasksJSON: no mask for channel " << cat->getLabel() << "\n";
            return false;
        }
        list.add(*mask);
        masks.erase(it);
    }
    if (!masks.empty()) {
        std::cerr << "RooSimultaneousOpt::importChannelMasksJSON: mask for unknown channel "
                  << masks.begin()->first << "\n";
        return false;
    }
    addChannelMasks(list);
    return true;
}
