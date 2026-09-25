#!/usr/bin/env python3
"""Faithful NLL cross-check: replicates combine's OWN internal NLL
construction (RooSimultaneousOpt + CachingSimNLL via CombineUtils.combineCreateNLL),
not plain RooAbsPdf::createNLL(), since combine does not use the latter for
RooSimultaneousOpt models. This mirrors what Combine.cc's MultiDimFit path
actually evaluates, including the isJSON extraConstraints-reattachment
heuristic, so it's a faithful comparison of the two model's actual likelihoods
at an identical, fully-frozen parameter point.
"""
import argparse
import sys

import ROOT

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


def load_root_workspace(path, ws_name="w"):
    f = ROOT.TFile.Open(path)
    w = f.Get(ws_name)
    if not w:
        raise RuntimeError(f"no workspace '{ws_name}' in {path}")
    return f, w


def build_json_optpdf(json_path):
    """Reimport the JSON and replicate Combine.cc:437-466's fixup:
    wrap the plain RooSimultaneous in a RooSimultaneousOpt and reattach
    extraConstraints by the same structural heuristic."""
    w = ROOT.RooWorkspace("w_json")
    tool = ROOT.RooJSONFactoryWSTool(w)
    if not tool.importJSON(json_path):
        raise RuntimeError(f"failed to import {json_path}")
    mc = w.obj("ModelConfig")
    if not mc:
        raise RuntimeError("no ModelConfig in json-imported workspace")
    pdf = mc.GetPdf()

    optpdf = ROOT.RooSimultaneousOpt(pdf, pdf.GetName() + "_opt")

    # Combine.cc:446-449 folds the category variables into the observables *before* the
    # dependsOn() scan below; without this the heuristic selects a different constraint set than
    # combine itself and the two NLLs are no longer comparable.
    obs_with_cats = ROOT.RooArgSet(mc.GetObservables())
    obs_with_cats.add(w.allCats(), True)
    mc.SetObservables(obs_with_cats)

    extra_constrs = ROOT.RooArgList()
    for cpdf in w.allPdfs():
        if not isinstance(cpdf, ROOT.RooAbsPdf):
            continue
        if cpdf.GetName() == optpdf.GetName():
            continue
        if isinstance(cpdf, ROOT.RooProdPdf):
            continue
        if cpdf.dependsOn(mc.GetObservables()):
            continue
        extra_constrs.add(cpdf)

    print(f"json: attaching {extra_constrs.getSize()} extraConstraints (structural heuristic)")
    if extra_constrs.getSize() > 0:
        optpdf.addExtraConstraints(extra_constrs)

    # Replicate this session's ModelConfig fix: recompute NuisanceParameters /
    # GlobalObservables from the now-complete constraint pdf set.
    new_nuis = ROOT.RooArgSet(mc.GetNuisanceParameters())
    new_globs = ROOT.RooArgSet(mc.GetGlobalObservables())
    for carg in extra_constrs:
        cpars = ROOT.RooArgSet()
        carg.getParameters(mc.GetObservables(), cpars)
        for p in cpars:
            if mc.GetParametersOfInterest().find(p.GetName()):
                continue
            if p.isConstant():
                new_globs.add(p, True)
            else:
                new_nuis.add(p, True)
    mc.SetNuisanceParameters(new_nuis)
    mc.SetGlobalObservables(new_globs)
    mc.SetPdf(optpdf)

    return w, mc, optpdf


def apply_snapshot(w, values):
    n_applied = 0
    n_missing = []
    for name, (val, lo, hi) in values.items():
        var = w.var(name)
        if var:
            if val < var.getMin() or val > var.getMax():
                var.setRange(min(lo, var.getMin()), max(hi, var.getMax()))
            var.setVal(val)
            n_applied += 1
        else:
            n_missing.append(name)
    return n_applied, n_missing


def get_snapshot_values(w):
    values = {}
    for v in w.allVars():
        values[v.GetName()] = (v.getVal(), v.getMin(), v.getMax())
    return values


def get_data(w):
    data = w.data("data_obs")
    if not data:
        datasets = w.allData()
        if datasets.size() == 0:
            raise RuntimeError("no dataset found")
        data = datasets.front()
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--snapshot", required=True)
    args = ap.parse_args()

    _, w_snap = load_root_workspace(args.snapshot)
    snapshot_values = get_snapshot_values(w_snap)
    print(f"Loaded {len(snapshot_values)} parameter values from snapshot workspace.")

    f_root, w_root = load_root_workspace(args.root)
    mc_root = w_root.obj("ModelConfig")
    pdf_root = mc_root.GetPdf()
    data_root = get_data(w_root)

    w_json, mc_json, pdf_json = build_json_optpdf(args.json)
    data_json = get_data(w_json)

    applied_root, missing_root = apply_snapshot(w_root, snapshot_values)
    applied_json, missing_json = apply_snapshot(w_json, snapshot_values)
    print(f"Applied snapshot to .root workspace: {applied_root} vars set, {len(missing_root)} missing")
    print(f"Applied snapshot to .json workspace: {applied_json} vars set, {len(missing_json)} missing")

    print(f"root: nuisance params = {mc_root.GetNuisanceParameters().getSize()}, "
          f"global obs = {mc_root.GetGlobalObservables().getSize()}")
    print(f"json: nuisance params = {mc_json.GetNuisanceParameters().getSize()}, "
          f"global obs = {mc_json.GetGlobalObservables().getSize()}")

    nll_root = ROOT.CombineUtils.combineCreateNLL(
        pdf_root, data_root, mc_root.GetNuisanceParameters(), False, False)
    nll_json = ROOT.CombineUtils.combineCreateNLL(
        pdf_json, data_json, mc_json.GetNuisanceParameters(), False, False)

    val_root = nll_root.getVal()
    val_json = nll_json.getVal()

    print(f"\nNLL(root, combine backend) = {val_root!r}")
    print(f"NLL(json, combine backend) = {val_json!r}")
    print(f"delta                      = {val_json - val_root!r}")

    tol = 1e-3 * max(1.0, abs(val_root))
    if abs(val_json - val_root) < tol:
        print(f"\nPASS: NLL values agree within tolerance ({tol!r}).")
        return 0
    else:
        print(f"\nFAIL: NLL values differ by more than tolerance ({tol!r}).")
        return 1


if __name__ == "__main__":
    sys.exit(main())
