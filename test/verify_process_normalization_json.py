#!/usr/bin/env python3
"""Standalone round-trip check for ProcessNormalization's HS3 JSON export/import.

Builds a ProcessNormalization with both a symmetric and an asymmetric
log-normal term plus an "other factor", exports it through
RooJSONFactoryWSTool, reimports it into a fresh workspace, and checks
that the value agrees across several parameter points. Exits non-zero
on mismatch.

Run inside the dev container after sourcing scripts/env.sh:
    python3 test/verify_process_normalization_json.py
"""
import sys
import tempfile

import ROOT

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")


def build_source_workspace():
    w = ROOT.RooWorkspace("w")
    theta_sym = ROOT.RooRealVar("theta_sym", "theta_sym", 0.0, -5, 5)
    theta_asym = ROOT.RooRealVar("theta_asym", "theta_asym", 0.0, -5, 5)
    other = ROOT.RooRealVar("other_factor", "other_factor", 1.5, 0, 10)

    pn = ROOT.ProcessNormalization("pn_test", "pn_test", 2.0)
    pn.addLogNormal(1.2, theta_sym)
    pn.addAsymmLogNormal(0.8, 1.35, theta_asym)
    pn.addOtherFactor(other)

    getattr(w, "import")(pn, ROOT.RooFit.RecycleConflictNodes())
    return w


def main():
    w = build_source_workspace()
    pn = w.function("pn_test")

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        json_path = f.name

    tool = ROOT.RooJSONFactoryWSTool(w)
    if not tool.exportJSON(json_path):
        print("FAIL: export returned false")
        return 1

    w2 = ROOT.RooWorkspace("w2")
    tool2 = ROOT.RooJSONFactoryWSTool(w2)
    if not tool2.importJSON(json_path):
        print("FAIL: import returned false")
        return 1

    pn2 = w2.function("pn_test")
    if not pn2:
        print("FAIL: pn_test not found in reimported workspace")
        return 1

    test_points = [
        ("theta_sym", "theta_asym", "other_factor"),
    ]
    thetas_sym = [-2.0, -1.0, 0.0, 1.0, 2.0]
    thetas_asym = [-2.0, -1.0, 0.0, 1.0, 2.0]

    ok = True
    for ts in thetas_sym:
        for ta in thetas_asym:
            w.var("theta_sym").setVal(ts)
            w.var("theta_asym").setVal(ta)
            w.var("other_factor").setVal(1.5)
            w2.var("theta_sym").setVal(ts)
            w2.var("theta_asym").setVal(ta)
            w2.var("other_factor").setVal(1.5)

            v1 = pn.getVal()
            v2 = pn2.getVal()
            if abs(v1 - v2) > 1e-9 * max(1.0, abs(v1)):
                print(f"MISMATCH at theta_sym={ts}, theta_asym={ta}: "
                      f"original={v1!r} reimported={v2!r}")
                ok = False

    if ok:
        print("PASS: ProcessNormalization JSON round-trip matches at all test points")
        return 0
    else:
        print("FAIL: ProcessNormalization JSON round-trip diverged (see above)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
