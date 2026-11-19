// Loading of HS3 JSON input files (combine <file>.json): see interface/CombineHS3Input.h.
#include "../interface/CombineHS3Input.h"

#include "../interface/RooSimultaneousOpt.h"

#include <RooAbsCategoryLValue.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooMsgService.h>
#include <RooProdPdf.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>
#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include <RooStats/ModelConfig.h>

#include <TString.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <typeinfo>
#include <utility>
#include <vector>

RooWorkspace *CombineHS3::importWorkspace(std::string const &fileName)
{
  auto *w = new RooWorkspace("w");
  RooJSONFactoryWSTool tool(*w);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  tool.importJSON(fileName);
  RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
  return w;
}

RooSimultaneousOpt *CombineHS3::makeSimultaneousOpt(RooWorkspace &w, RooStats::ModelConfig &mc, RooSimultaneous &pdf,
                                                    std::string const &fileName, int verbose)
{
  RooArgList extraConstrs;
  // autoMCStats constraints (prop_bin<channel>_bin...) per channel index, in channel order
  // (bin number, constraint) per channel index
  std::map<int, std::vector<std::pair<long, RooAbsArg *>>> bbConstrs;
  // adding cats to observables
  RooArgSet obsWithCats(*mc.GetObservables());
  obsWithCats.add(w.allCats(), /*silent=*/true);
  mc.SetObservables(obsWithCats);
  // The histfactory_dist importers build a channel whose modifiers name constraints as
  // RooProdPdf(<name>, [<name>_model, constraints]) (CombineHS3Channel.h). Combine's own
  // models keep the channels free of constraints and carry all of them as extra
  // constraints of the RooSimultaneousOpt, in ModelConfig nuisance order. CachingSimNLL
  // factorizes the constraints out either way (same likelihood), but it hands the
  // parameters to the minimizer in the order it meets the constraints, and Minuit's path
  // depends on that order. So the channels go into the RooSimultaneousOpt without their
  // constraints, and the constraints follow as extra constraints (sorted below).
  auto *optpdf = new RooSimultaneousOpt(TString(mc.GetPdf()->GetName())+"_opt", pdf.GetTitle(),
                                        const_cast<RooAbsCategoryLValue &>(pdf.indexCat()));
  for (auto const &a : pdf.attributes()) optpdf->setAttribute(a.c_str());
  for (auto const &a : pdf.stringAttributes()) optpdf->setStringAttribute(a.first.c_str(), a.second.c_str());
  std::unique_ptr<RooAbsCategoryLValue> cat{static_cast<RooAbsCategoryLValue *>(pdf.indexCat().Clone())};
  for (int ic = 0, nc = cat->numBins(nullptr); ic < nc; ++ic) {
    cat->setBin(ic);
    RooAbsPdf *chan = pdf.getPdf(cat->getLabel());
    if (!chan) continue;
    RooAbsPdf *core = chan;
    if (auto *prod = dynamic_cast<RooProdPdf *>(chan); prod && typeid(*prod) == typeid(RooProdPdf)) {
      RooAbsPdf *obsPart = nullptr; RooArgList cons; bool flat = true;
      for (RooAbsArg *f : prod->pdfList()) {
        auto *fp = static_cast<RooAbsPdf *>(f);
        if (fp->dependsOn(obsWithCats)) { if (obsPart) flat = false; obsPart = fp; }
        else if (dynamic_cast<RooProdPdf *>(fp)) flat = false;
        else cons.add(*fp);
      }
      if (flat && obsPart) {
        core = obsPart;
        const std::string bbPrefix = std::string("prop_bin") + cat->getCurrentLabel() + "_bin";
        for (RooAbsArg *c : cons) {
          const std::string cname = c->GetName();
          if (cname.rfind(bbPrefix, 0) == 0)
            bbConstrs[cat->getCurrentIndex()].emplace_back(
                std::strtol(cname.c_str() + bbPrefix.size(), nullptr, 10), c);
          else extraConstrs.add(*c, /*silent=*/true);
        }
      }
    }
    optpdf->addPdf(*core, cat->getLabel());
  }
  // Constraints the histfactory_dist importers already put into the channels (modifier
  // `constraint`, CombineHS3Channel.h) are part of the model and must not be added twice.
  std::unique_ptr<RooArgSet> modelComponents{pdf.getComponents()};
  // Nor the constraints that only another model's channels use (e.g. the b-only model's own
  // autoMCStats constraints, CombineHS3BBLite.h): they are not part of this model.
  RooArgSet otherComponents;
  for (RooAbsArg *arg : w.allPdfs()) {
    auto *other = dynamic_cast<RooSimultaneous *>(arg);
    if (!other || other == &pdf) continue;
    std::unique_ptr<RooArgSet> comps{other->getComponents()};
    otherComponents.add(*comps, /*silent=*/true);
  }
  for (RooAbsArg *arg : w.allPdfs()) {
    auto *cpdf = dynamic_cast<RooAbsPdf *>(arg);

    // skip if cpdf is not a pdf
    if (!cpdf || cpdf == optpdf || dynamic_cast<RooProdPdf *>(cpdf)) continue;

    // skip if cpdf depends on a observable
    if (cpdf->dependsOn(*mc.GetObservables())) continue;

    // skip if a channel already contains it (taken out of the channels above)
    if (modelComponents->find(*cpdf)) continue;

    // skip if only another model contains it
    if (otherComponents.find(*cpdf)) continue;

    extraConstrs.add(*cpdf, /*silent=*/true);
  }
  // Order of the imported ModelConfig's nuisance parameters, i.e. of the JSON's
  // nuisance-parameter domain, which ROOT writes sorted by name (as combineCards sorts
  // a datacard's systematics): each constraint by the first nuisance parameter among
  // its direct servers; constraints without one keep their relative order after the others.
  {
    std::map<std::string, std::size_t> rank;
    std::size_t k = 0;
    if (mc.GetNuisanceParameters())
      for (RooAbsArg *p : *mc.GetNuisanceParameters()) rank.emplace(p->GetName(), k++);
    std::vector<std::pair<std::size_t, RooAbsArg *>> keyed;
    for (RooAbsArg *c : extraConstrs) {
      std::size_t key = rank.size();
      for (RooAbsArg *s : c->servers()) {
        auto it = rank.find(s->GetName());
        if (it != rank.end()) key = std::min(key, it->second);
      }
      keyed.emplace_back(key, c);
    }
    std::stable_sort(keyed.begin(), keyed.end(),
                     [](auto const &a, auto const &b) { return a.first < b.first; });
    RooArgList sorted;
    for (auto const &kc : keyed) sorted.add(*kc.second);
    extraConstrs.removeAll();
    extraConstrs.add(sorted);
    // text2workspace adds the autoMCStats constraints after the datacard's, channel by
    // channel in category index order, bin by bin (a bin's per-process terms in the
    // channel's order).
    for (auto &bb : bbConstrs) {
      std::stable_sort(bb.second.begin(), bb.second.end(),
                       [](auto const &a, auto const &b) { return a.first < b.first; });
      for (auto const &bc : bb.second) extraConstrs.add(*bc.second, /*silent=*/true);
    }
  }
  if (extraConstrs.getSize() > 0) {
    if (verbose > 1)
      std::cout << "Attaching " << extraConstrs.getSize() << " constraint PDFs from JSON import." << std::endl;
    optpdf->addExtraConstraints(extraConstrs);
  }
  if (!optpdf->importChannelMasksJSON(fileName.c_str(), w)) {
    throw std::invalid_argument("Invalid misc.CMS.channel_masks in " + fileName);
  }
  return optpdf;
}
