#include "../interface/CMSHistFuncWrapper.h"
#include <vector>
#include <ostream>
#include "RooRealProxy.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "TH1F.h"

#define HFVERBOSE 0

CMSHistFuncWrapper::CMSHistFuncWrapper()
    : pfunc_(nullptr), perr_(nullptr), initialized_(false) {
  idx_ = 0;
}

CMSHistFuncWrapper::CMSHistFuncWrapper(const char* name, const char* title, RooRealVar& x,
              CMSHistFunc & func, CMSHistErrorPropagator & err, unsigned idx)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      func_("func", "", this, func),
      err_("err", "", this, err, true, false),
      idx_(idx),
      sentry_(TString(name) + "_sentry", ""),
      pfunc_(nullptr),
      perr_(nullptr),
      initialized_(false) {
  cache_ = func.cache();
}

CMSHistFuncWrapper::CMSHistFuncWrapper(CMSHistFuncWrapper const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      func_("func", this, other.func_),
      err_("err", this, other.err_),
      cache_(other.cache_),
      idx_(other.idx_),
      sentry_(name ? TString(name) + "_sentry" : TString(other.GetName())+"_sentry", ""),
      pfunc_(nullptr),
      perr_(nullptr),
      initialized_(false) {
}

void CMSHistFuncWrapper::initialize() const {
  if (initialized_) return;
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  pfunc_ = dynamic_cast<CMSHistFunc const*>(&(func_.arg()));
  perr_ = dynamic_cast<CMSHistErrorPropagator *>(err_.absArg());
  auto sentry_args = perr_->getSentryArgs();
  for (RooAbsArg *arg : *sentry_args) {
    sentry_.addArg(*arg);
  }
  sentry_.setValueDirty();
  initialized_ = true;
}

void CMSHistFuncWrapper::updateCache() const {
  initialize();

  // The ErrorPropagator will send a dirty flag whenever we need to update the cache
  if (!sentry_.good()) {
    perr_->applyErrorShifts(idx_, pfunc_->cache(), cache_);
    // cache_.CropUnderflows();
#if HFVERBOSE > 2
    std::cout << "Updated cache from CMSHistFunc:\n";
    pfunc_->cache().Dump();
    std::cout << "After shifts and cropping:\n";
    cache_.Dump();
#endif
  }
  sentry_.reset();
}

Double_t CMSHistFuncWrapper::evaluate() const {
  updateCache();
  return cache_.GetAt(x_);

}

void CMSHistFuncWrapper::printMultiline(std::ostream& os, Int_t contents,
                                 Bool_t verbose, TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  std::cout << ">> Current cache:\n";
  cache_.Dump();
  std::cout << ">> Sentry: " << sentry_.good() << "\n";
  sentry_.Print("v");
}

Int_t CMSHistFuncWrapper::getAnalyticalIntegral(RooArgSet& allVars,
                                         RooArgSet& analVars,
                                         const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_)) return 1;
  return 0;
}

Double_t CMSHistFuncWrapper::analyticalIntegral(Int_t code,
                                         const char* rangeName) const {
  switch (code) {
    case 1: {
      updateCache();
      return cache_.IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

#include "RooFitHS3/RooJSONFactoryWSTool.h"
#include "RooFitHS3/JSONIO.h"
#include "static_execute.h"

namespace {
using RooFit::Detail::JSONNode;

// CMSHistFuncWrapper is a plotting/shape-breakdown convenience (see --no-wrappers): it just
// re-evaluates a process's already-exported CMSHistFunc shape with the shared
// CMSHistErrorPropagator's per-bin Barlow-Beeston correction applied -- no information beyond
// what the parent histfactory_dist already carries, and text2workspace.py imports one into the
// output workspace per process unconditionally (safe_import(), regardless of --no-wrappers'
// existence), so ROOT's own top-level export walks and tries to export every one of them whether
// or not the actual likelihood references it.
//
// A *declining* exporter does not help here (confirmed by reading RooJSONFactoryWSTool::
// exportObject, ~L1010-1024 of root@fbad13f: the {"name": ...} stub is created before the
// exporter is even looked up, and a `false` return explicitly clears any partial content but
// leaves that stub in place) -- so this must actually succeed. There is no way to represent "this
// process's shape with the propagator's live per-bin correction applied" as a reusable HS3 node
// without re-deriving the parent histfactory_dist's own staterror gammas per-bin, which is exactly
// the information histfactory_dist already carries once reimported. Instead, export a plain static
// histogram snapshot at the current (export-time) parameter values, via the same
// RooJSONFactoryWSTool::exportHisto() helper and "step" type key ROOT's own RooHistFunc exporter
// uses (JSONFactories_RooFitCore.cxx) -- reimports as an ordinary RooHistFunc. This is honestly
// lossy (frozen, not tied to the live gamma nuisances) but the object is decorative: nothing in
// the core likelihood depends on it, so a frozen snapshot is strictly better than the current
// "aborts the whole import" status quo, and no CMS-specific HS3 type is needed for it.
class CMSHistFuncWrapperStreamer : public RooFit::JSONIO::Exporter {
public:
   // Don't recurse into func_/err_/x_: they're CMSHistFuncWrapper's *servers* in the RooFit
   // sense, but this exporter already inlined everything it needs from them (a static snapshot).
   // Neither CMSHistFunc nor CMSHistErrorPropagator has a direct-class exporter of their own
   // (CMSHistErrorPropagator is only reachable via the RooRealSumPdf-wrapping pattern
   // tryExportHistFactoryFromPropagator recognizes), so leaving this at the default (true) makes
   // the framework try and fail to export them too.
   bool autoExportDependants() const override { return false; }
   std::string const &key() const override
   {
      static const std::string keystring = "step";
      return keystring;
   }
   bool exportObject(RooJSONFactoryWSTool *tool, const RooAbsArg *p, JSONNode &elem) const override
   {
     auto const *wrapper = dynamic_cast<const CMSHistFuncWrapper *>(p);
     if (!wrapper) return false;

     wrapper->updateCache();
     auto const &c = wrapper->cache();

     // exportHisto() derives the "axes" bin count from the RooRealVar's own numBins(), not from
     // the content array size -- and the shared model observable's binning isn't guaranteed to
     // match this particular wrapper's cache size (multi-channel workspaces have one observable
     // per channel, with per-channel binning). Use a local, throwaway clone with setBins() forced
     // to match, rather than relying on -- or mutating -- the live model variable.
     auto const *xOrig = dynamic_cast<RooRealVar const *>(&wrapper->getXVar());
     if (!xOrig) return false;
     RooRealVar xTmp(xOrig->GetName(), xOrig->GetTitle(), xOrig->getMin(), xOrig->getMax());
     xTmp.setBins(static_cast<int>(c.fullsize()));
     RooArgSet obs(xTmp);
     tool->exportHisto(obs, c.fullsize(), c.Values().data(), elem["data"].set_map());
     elem["type"] << key();
     return true;
   }
};

STATIC_EXECUTE([]() {
  using namespace RooFit::JSONIO;
  registerExporter<CMSHistFuncWrapperStreamer>(CMSHistFuncWrapper::Class(), false);
});
}  // namespace

#undef HFVERBOSE
