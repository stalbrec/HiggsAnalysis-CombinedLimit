#include <RooFitHS3/RooJSONFactoryWSTool.h>
#include <RooFitHS3/JSONIO.h>
#include "CombineHS3Helpers.h"
#include "static_execute.h"

#include "RooFit.h"

#include "../interface/ProcessNormalization.h"
#include "../interface/CombineHS3Channel.h"

#include "RooAddition.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooNumber.h"
#include "RooPoisson.h"
#include "RooProdPdf.h"
#include "RooSimultaneous.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

#include <charconv>
#include <map>
#include <set>
#include <memory>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

using RooFit::Detail::JSONNode;

namespace {

using CombineHS3::approxSymmetricKappas;
using CombineHS3::getOrMakeParameter;
using CombineHS3::isSingleBinCountingDist;
using CombineHS3::parseDouble;
using CombineHS3::requirePlainModifierParameter;
using CombineHS3::resolveConstraintPdfName;
using CombineHS3::writeNormSys;

// Combine's counting-only models (python/ModelTools.py, CountingModelBuilder) express a channel as
//
//   RooPoisson::pdf_bin<b>[ x=n_obs_bin<b> mean=n_exp_bin<b> ]
//   RooAddition::n_exp_bin<b>[ n_exp_bin<b>_proc_<p> + ... ]
//   ProcessNormalization::n_exp_bin<b>_proc_<p>[ thetaList=(...) asymmThetaList=(...) otherFactorList=(...) ]


void exportOtherFactorModifiers(RooJSONFactoryWSTool *tool, JSONNode &mods, RooAbsReal *factor)
{
   auto &m = mods.append_child().set_map();
   m["name"] << factor->GetName();
   m["type"] << "normfactor";
   if (auto *rcv = dynamic_cast<RooConstVar *>(factor)) {
      m["parameter"] << rcv->getVal();
   } else {
      m["parameter"] << factor->GetName();
      tool->queueExport(*factor);
   }
}


void writeAxisBound(JSONNode &node, double value)
{
   if (RooNumber::isInfinite(value)) {
      node.set_null();
   } else {
      node << value;
   }
}

double readAxisBound(JSONNode const &node, const char *key, double defaultValue)
{
   auto const *bound = node.find(key);
   if (!bound || bound->is_null())
      return defaultValue;
   return bound->val_double();
}

/// Write the single-element data.contents array of a sample.
void writeSingleBinContents(JSONNode &sample, double value)
{
   auto &dataNode = sample["data"].set_map();
   RooJSONFactoryWSTool::exportArray(1, &value, dataNode["contents"]);
}



void exportNormSys(RooJSONFactoryWSTool *tool, JSONNode &mods, RooAbsArg const &theta, double kappaLo, double kappaHi,
                   RooAbsCollection const &observables, RooAbsArg const &channel,
                   std::map<std::string, std::string> &constraintCache, std::string const &sampleName)
{
   requirePlainModifierParameter(theta, channel.GetName(), sampleName, "normsys");
   const std::string constraint = resolveConstraintPdfName(tool, theta, observables, &channel, constraintCache);
   writeNormSys(tool, mods, theta, kappaLo, kappaHi, constraint);
}

// single bin counting experiment exporter
bool tryExportCountingHistFactory(RooJSONFactoryWSTool *tool, const RooPoisson *pois, JSONNode &elem)
{
   auto const *x = dynamic_cast<RooRealVar const *>(&pois->getX());
   if (!x) {
      return false;
   }
   auto const *add = dynamic_cast<RooAddition const *>(&pois->getMean());
   if (!add) {
      return false;
   }

   // Decline pass: everything is checked before a single byte is written into elem, so that a
   // decline leaves neither a half-filled node nor queued servers behind.
   RooArgList const &summands = add->list();
   if (summands.getSize() == 0) {
      RooJSONFactoryWSTool::warning("SumOfRooPoisson '" + std::string(pois->GetName()) + "': no processes found!");
      return false;
   }
   for (auto *summand : summands) {
      if (!dynamic_cast<ProcessNormalization const *>(summand) && !dynamic_cast<RooRealVar const *>(summand) &&
          !dynamic_cast<RooConstVar const *>(summand)) {
         RooJSONFactoryWSTool::warning("SumOfRooPoisson '" + std::string(pois->GetName()) + "': summand '" +
                                       summand->GetName() + "' is a " + summand->ClassName() +
                                       ", not a ProcessNormalization or a constant!");
         return false;
      }
   }

   auto &axes = elem["axes"].set_seq();
   auto &ax = axes.append_child().set_map();
   std::string axname(x->GetName());
   RooJSONFactoryWSTool::testValidName(axname, false);
   ax["name"] << axname;
   ax["nbins"] << 1;

   
   writeAxisBound(ax["min"], x->getMin());
   writeAxisBound(ax["max"], x->getMax());

   auto &samples = elem["samples"].set_seq();

   RooArgSet observables;
   observables.add(*x);
   std::map<std::string, std::string> constraintCache;

   for (auto *summand : summands) {
      std::string sampleName(summand->GetName());
      RooJSONFactoryWSTool::testValidName(sampleName, false);

      auto &s = samples.append_child().set_map();
      s["name"] << sampleName;

      if (auto const *pn = dynamic_cast<ProcessNormalization const *>(summand)) {
         writeSingleBinContents(s, pn->nominalValue());
         auto &mods = s["modifiers"].set_seq();

         // symmetric log-normals: kappaLo == 1/kappaHi
         auto const &logKappa = pn->logKappa();
         for (int i = 0; i < static_cast<int>(logKappa.size()); ++i) {
            exportNormSys(tool, mods, pn->thetaList()[i], std::exp(-logKappa[i]), std::exp(logKappa[i]), observables,
                          *pois, constraintCache, sampleName);
         }
         // asymmetric log-normals
         auto const &logAsymmKappa = pn->logAsymmKappa();
         for (int i = 0; i < static_cast<int>(logAsymmKappa.size()); ++i) {
            exportNormSys(tool, mods, pn->asymmThetaList()[i], std::exp(logAsymmKappa[i].first),
                          std::exp(logAsymmKappa[i].second), observables, *pois, constraintCache, sampleName);
         }
         for (auto *of : pn->otherFactorList()) {
            exportOtherFactorModifiers(tool, mods, static_cast<RooAbsReal *>(of));
         }
      } else if (auto const *rrv = dynamic_cast<RooRealVar const *>(summand); rrv && !rrv->isConstant()) {
         // A floating yield: keep it floating instead of freezing its current value.
         writeSingleBinContents(s, 1.0);
         auto &mods = s["modifiers"].set_seq();
         auto &m = mods.append_child().set_map();
         m["name"] << rrv->GetName();
         m["parameter"] << rrv->GetName();
         m["type"] << "normfactor";
         tool->queueExport(*rrv);
      } else {
         writeSingleBinContents(s, static_cast<RooAbsReal const *>(summand)->getVal());
         s["modifiers"].set_seq();
      }
   }

   return true;
}

class SumOfRooPoissonStreamer_RooAddition : public RooFit::JSONIO::Exporter {
public:
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring("histfactory_dist");
      return keystring;
   }

   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {

      auto const *pois = dynamic_cast<const RooPoisson *>(p);
      if (!pois)
         return false;
      if (!tryExportCountingHistFactory(tool, pois, elem))
         return false;
      elem["type"] << key();

      return true;
   }
};

class SumOfRooPoissonFactory_RooAddition : public RooFit::JSONIO::Importer {
public:
   bool importArg(RooJSONFactoryWSTool *tool, JSONNode const &elem) const override
   {
      if (!elem.has_child("name"))
         return false;
      // Multi-bin (and any shape-carrying) node belongs to the hist-stack importers.
      if (!isSingleBinCountingDist(elem))
         return false;

      const std::string distName = elem["name"].val();

      // Observable: reuse the workspace variable if it is already there (e.g. imported with the
      // observed dataset), otherwise build it from the axis. Read the axis by hand rather than
      // through readAxes(), which cannot deal with the null bounds standing for +/-infinity.
      std::string axName;
      double axMin = -RooNumber::infinity();
      double axMax = RooNumber::infinity();
      for (auto const &ax : elem["axes"].children()) {
         axName = ax["name"].val();
         axMin = readAxisBound(ax, "min", -RooNumber::infinity());
         axMax = readAxisBound(ax, "max", RooNumber::infinity());
      }
      if (axName.empty())
         throw std::runtime_error("SumOfRooPoisson importer: '" + distName + "' has no usable axis");

      auto *x = tool->workspace()->var(axName);
      std::unique_ptr<RooRealVar> ownedX;
      if (!x) {
         // The value-only constructor leaves the range at (-inf, +inf), matching how combine
         // creates n_obs_bin<b>; only finite bounds have to be applied on top.
         ownedX = std::make_unique<RooRealVar>(axName.c_str(), axName.c_str(), 0.0);
         if (!RooNumber::isInfinite(axMin))
            ownedX->setMin(axMin);
         if (!RooNumber::isInfinite(axMax))
            ownedX->setMax(axMax);
         x = ownedX.get();
      }

      RooArgList summands;

      for (const auto &s : elem["samples"].children()) {
         const std::string sampleName = s["name"].val();

         // Samples are named after the ProcessNormalization they came from, so a process shared
         // between the s and the b-only distribution is reconstructed exactly once.
         if (auto *existing = tool->workspace()->function(sampleName)) {
            summands.add(*existing);
            continue;
         }

         double nominal = 0.;
         for (auto const &c : s["data"]["contents"].children()) {
            nominal = c.val_double();
            break;
         }

         RooArgList normsysPars;
         std::vector<std::pair<double, double>> normsysVals;
         RooArgList normfactorPars;

         if (s.has_child("modifiers")) {
            for (const auto &m : s["modifiers"].children()) {
               const std::string mtype = m["type"].val();
               const std::string mname = m["name"].val();

               if (mtype == "normsys") {
                  auto const *par = m.find("parameter");
                  std::string parname(par ? par->val() : "alpha_" + mname);
                  RooRealVar &p = getOrMakeParameter(tool, parname, 0.0, -10.0, 10.0);
                  auto const &data = m["data"];
                  normsysPars.add(p);
                  normsysVals.emplace_back(data["lo"].val_double(), data["hi"].val_double());
               } else if (mtype == "normfactor") {
                  auto const *par = m.find("parameter");
                  const std::string parval(par ? par->val() : mname);
                  double literal = 0.;
                  if (parseDouble(parval, literal)) {
                     // Written by the exporter for a RooConstVar other-factor: recreate it under
                     // its original name so the round trip is name-stable.
                     if (auto *existingCv = tool->workspace()->function(mname)) {
                        normfactorPars.add(*existingCv);
                     } else {
                        RooConstVar cv(mname.c_str(), mname.c_str(), literal);
                        tool->workspace()->import(cv, RooFit::RecycleConflictNodes());
                        normfactorPars.add(*tool->workspace()->function(mname));
                     }
                  } else {
                     normfactorPars.add(*tool->request<RooAbsReal>(parval, sampleName));
                  }
               } else if (mtype == "custom") {
                  normfactorPars.add(*tool->request<RooAbsReal>(mname, sampleName));
               } else {
                  throw std::runtime_error("SumOfRooPoisson importer: modifier type '" + mtype + "' on sample '" +
                                           sampleName + "' is not supported by this importer");
               }
            }
         }

         ProcessNormalization pn(sampleName.c_str(), sampleName.c_str(), nominal);
         for (int i = 0; i < normsysPars.getSize(); ++i) {
            auto &par = static_cast<RooAbsReal &>(normsysPars[i]);
            if (approxSymmetricKappas(normsysVals[i].first, normsysVals[i].second)) {
               pn.addLogNormal(normsysVals[i].second, par);
            } else {
               pn.addAsymmLogNormal(normsysVals[i].first, normsysVals[i].second, par);
            }
         }
         for (int i = 0; i < normfactorPars.getSize(); ++i) {
            pn.addOtherFactor(static_cast<RooAbsReal &>(normfactorPars[i]));
         }
         tool->workspace()->import(pn, RooFit::RecycleConflictNodes());

         auto *imported = tool->workspace()->function(sampleName);
         if (!imported)
            throw std::runtime_error("SumOfRooPoisson importer: failed to import ProcessNormalization '" + sampleName +
                                     "'");
         summands.add(*imported);
      }

      // HS3 has no field to carry the name of the inner sum, so it is derived from the
      // distribution name. Sample names and the distribution name itself round-trip exactly.
      const std::string sumName = distName + "_n_exp";
      RooAddition addition(sumName.c_str(), sumName.c_str(), summands);
      tool->workspace()->import(addition, RooFit::RecycleConflictNodes());
      auto *n_exp = tool->workspace()->function(sumName);
      if (!n_exp)
         throw std::runtime_error("SumOfRooPoisson importer: failed to import RooAddition '" + sumName + "'");

      // combine always builds these counting Poissons with noRounding=1.
      const std::string coreName = CombineHS3::channelCoreName(elem);
      RooPoisson pdf(coreName.c_str(), distName.c_str(), *x, *n_exp, /*noRounding=*/true);
      tool->workspace()->import(pdf, RooFit::RecycleConflictNodes());
      if (!tool->workspace()->pdf(coreName))
         throw std::runtime_error("SumOfRooPoisson importer: failed to import RooPoisson '" + coreName + "'");
      CombineHS3::wrapChannelInConstraints(tool, elem);

      return true;
   }
};

} // namespace

// ---------------------------------------------------------------------------------------------
// Constraint terms inside the channels (CombineHS3Channel.h)
// ---------------------------------------------------------------------------------------------
namespace {
// While exportConstrainedChannel() exports the core of a RooProdPdf(<name>, [<name>_model,
// constraints]), the product stands in for the core in constraintAppliedOutside(): its own
// constraint factors are the channel's, not an enclosing product's.
thread_local RooAbsArg const *g_wrappedCore = nullptr;
thread_local RooAbsArg const *g_wrapper = nullptr;
} // namespace

std::vector<std::string> CombineHS3::modifierConstraintNames(JSONNode const &elem)
{
   std::vector<std::string> names;
   std::set<std::string> seen;
   auto const *samples = elem.find("samples");
   if (!samples)
      return names;
   for (auto const &s : samples->children()) {
      auto const *mods = s.find("modifiers");
      if (!mods)
         continue;
      for (auto const &m : mods->children()) {
         auto const *c = m.find("constraint");
         if (!c || c->is_null())
            continue;
         const std::string name = c->val();
         if (!name.empty() && seen.insert(name).second)
            names.push_back(name);
      }
   }
   return names;
}

std::string CombineHS3::channelCoreName(JSONNode const &elem)
{
   const std::string name = elem["name"].val();
   return modifierConstraintNames(elem).empty() ? name : name + "_model";
}

void CombineHS3::wrapChannelInConstraints(RooJSONFactoryWSTool *tool, JSONNode const &elem)
{
   const std::vector<std::string> names = modifierConstraintNames(elem);
   if (names.empty())
      return;
   const std::string distName = elem["name"].val();
   RooAbsPdf *core = tool->workspace()->pdf(distName + "_model");
   if (!core)
      throw std::runtime_error("histfactory_dist importer: channel '" + distName + "_model' was not imported");
   RooArgList factors(*core);
   for (auto const &n : names)
      factors.add(*tool->request<RooAbsPdf>(n, distName));
   RooProdPdf prod(distName.c_str(), distName.c_str(), factors);
   tool->workspace()->import(prod, RooFit::RecycleConflictNodes(), RooFit::Silence());
   if (!tool->workspace()->pdf(distName))
      throw std::runtime_error("histfactory_dist importer: failed to import RooProdPdf '" + distName + "'");
}

bool CombineHS3::constraintAppliedOutside(RooAbsArg const &channel, RooAbsArg const &constraint)
{
   RooAbsArg const *start = (&channel == g_wrappedCore && g_wrapper) ? g_wrapper : &channel;
   std::vector<RooAbsArg const *> todo{start};
   std::set<RooAbsArg const *> visited{start};
   while (!todo.empty()) {
      RooAbsArg const *node = todo.back();
      todo.pop_back();
      for (RooAbsArg const *client : node->clients()) {
         if (auto const *prod = dynamic_cast<RooProdPdf const *>(client)) {
            if (!prod->pdfList().containsInstance(*node))
               continue;
            for (RooAbsArg const *f : prod->pdfList()) {
               if (f == node)
                  continue;
               if (f == &constraint || f->dependsOn(constraint))
                  return true;
            }
         } else if (!dynamic_cast<RooSimultaneous const *>(client)) {
            // the other channels of a RooSimultaneous are not multiplied with this one
            continue;
         }
         if (visited.insert(client).second)
            todo.push_back(client);
      }
   }
   return false;
}

bool CombineHS3::exportConstrainedChannel(RooJSONFactoryWSTool *tool, RooProdPdf const &prod, JSONNode &elem)
{
   RooArgList const &factors = prod.pdfList();
   if (factors.size() < 2)
      return false;
   RooAbsArg const *core = factors.at(0);
   const std::string name = prod.GetName();
   if (core->GetName() != name + "_model")
      return false;

   auto it = RooFit::JSONIO::exporters().find(core->IsA());
   if (it == RooFit::JSONIO::exporters().end())
      return false;

   auto reset = [&]() {
      elem.clear();
      elem.set_map();
      elem["name"] << name;
   };
   bool done = false;
   g_wrappedCore = core;
   g_wrapper = &prod;
   try {
      for (auto const &exp : it->second) {
         if (exp->key() != "histfactory_dist")
            continue;
         if (exp->exportObject(tool, core, elem)) {
            if (exp->autoExportDependants()) {
               for (RooAbsArg *server : core->servers())
                  tool->queueExport(*server);
            }
            done = true;
            break;
         }
         reset();
      }
   } catch (...) {
      g_wrappedCore = g_wrapper = nullptr;
      throw;
   }
   g_wrappedCore = g_wrapper = nullptr;
   if (!done)
      return false;

   // Only if importing the result gives this product back: same constraints, same order.
   const std::vector<std::string> named = modifierConstraintNames(elem);
   bool same = named.size() + 1 == factors.size();
   for (std::size_t i = 0; same && i < named.size(); ++i)
      same = named[i] == factors.at(i + 1)->GetName();
   if (!same) {
      reset();
      return false;
   }
   for (std::size_t i = 1; i < factors.size(); ++i)
      tool->queueExport(*factors.at(i));
   elem["name"] << name;
   return true;
}

STATIC_EXECUTE([]() {
   using namespace RooFit::JSONIO;

   registerExporter<SumOfRooPoissonStreamer_RooAddition>("RooPoisson", true);
   registerImporter<SumOfRooPoissonFactory_RooAddition>("histfactory_dist", true);
});
