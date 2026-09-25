#!/usr/bin/env python3
"""Round-trip check for `shape`-systematic RooAddPdf terms through HS3 JSON.

Covers the structure text2workspace.py builds for a process carrying `shape` nuisances:

    coefficient = RooProduct(ProcessNormalization, AsymPow_1, AsymPow_2, ...)
    shape       = FastVerticalInterpHistPdf2(coefList=(theta_1, ...), x=obs)

Each such theta maps to *two* modifiers on the same sample -- a histosys (shape) and a normsys
(yield, from the sibling AsymPow) -- sharing one parameter and one constraint. The three cases
below exercise the pairing, the two "half" cases that must keep working, and the decline path.

Run inside the dev container after sourcing scripts/env.sh:
    python3 test/verify_shape_systematics_json.py
"""
import json
import sys
import tempfile

import ROOT

ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
ROOT.TH1.AddDirectory(False)

NBINS = 4
XLO, XHI = 0.0, 4.0
TOL = 1e-9

# Kept alive for the lifetime of the process: PyROOT would otherwise garbage-collect histograms
# still referenced by a TList, or RooFit objects still owned by a workspace under construction.
_keepalive = []


def th1(name, contents):
    h = ROOT.TH1D(name, name, NBINS, XLO, XHI)
    for i, c in enumerate(contents):
        h.SetBinContent(i + 1, c)
    _keepalive.append(h)
    return h


def make_shape_pdf(name, x, nominal, variations):
    """FastVerticalInterpHistPdf2 with one (hi, lo) template pair per theta.

    The constructor wants a TList ordered [nominal, hi0, lo0, hi1, lo1, ...].
    """
    funcs = ROOT.TList()
    funcs.Add(nominal)
    coefs = ROOT.RooArgList()
    for theta, hi, lo in variations:
        funcs.Add(hi)
        funcs.Add(lo)
        coefs.add(theta)
    _keepalive.append(funcs)
    pdf = ROOT.FastVerticalInterpHistPdf2(name, name, x, funcs, coefs)
    pdf.setStringAttribute("combine.process", name)
    _keepalive.append(pdf)
    return pdf


def make_theta(w, name):
    """A nuisance plus the SimpleGaussianConstraint named <theta>_Pdf that combine gives it.

    The exporter resolves the HS3 "constraint" key to this pdf's name.
    """
    theta = ROOT.RooRealVar(name, name, 0.0, -5.0, 5.0)
    glob = ROOT.RooRealVar(name + "_In", name + "_In", 0.0, -5.0, 5.0)
    glob.setConstant(True)
    one = ROOT.RooConstVar("sigma_one", "sigma_one", 1.0)
    constraint = ROOT.SimpleGaussianConstraint(name + "_Pdf", name + "_Pdf", theta, glob, one)
    getattr(w, "import")(constraint, ROOT.RooFit.RecycleConflictNodes())
    _keepalive.extend([theta, glob, one, constraint])
    return w.var(name)


def asym_pow(name, theta, kappa_lo, kappa_hi):
    klo = ROOT.RooConstVar(name + "_kappaLo", "", kappa_lo)
    khi = ROOT.RooConstVar(name + "_kappaHi", "", kappa_hi)
    ap = ROOT.AsymPow(name, name, klo, khi, theta)
    _keepalive.extend([klo, khi, ap])
    return ap


def roundtrip(w):
    """Export w to HS3 JSON and reimport into a fresh workspace. Returns (w2, parsed json)."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        path = f.name
    tool = ROOT.RooJSONFactoryWSTool(w)
    if not tool.exportJSON(path):
        raise RuntimeError("exportJSON returned false")
    w2 = ROOT.RooWorkspace("w2")
    if not ROOT.RooJSONFactoryWSTool(w2).importJSON(path):
        raise RuntimeError("importJSON returned false")
    with open(path) as fh:
        return w2, json.load(fh)


def compare(w, w2, pdf_name, thetas, points):
    """Compare expectedEvents and every bin of the pdf at each parameter point."""
    p1, p2 = w.pdf(pdf_name), w2.pdf(pdf_name)
    if not p2:
        return ["pdf '%s' missing from the reimported workspace" % pdf_name]
    x1, x2 = w.var("x"), w2.var("x")
    obs1, obs2 = ROOT.RooArgSet(x1), ROOT.RooArgSet(x2)

    failures = []
    for point in points:
        for name, val in zip(thetas, point):
            w.var(name).setVal(val)
            w2.var(name).setVal(val)

        n1, n2 = p1.expectedEvents(obs1), p2.expectedEvents(obs2)
        if abs(n1 - n2) > TOL * max(1.0, abs(n1)):
            failures.append("expectedEvents at %s: %r vs %r" % (dict(zip(thetas, point)), n1, n2))

        for ib in range(NBINS):
            xv = XLO + (ib + 0.5) * (XHI - XLO) / NBINS
            x1.setVal(xv)
            x2.setVal(xv)
            v1, v2 = p1.getVal(obs1), p2.getVal(obs2)
            if abs(v1 - v2) > TOL * max(1.0, abs(v1)):
                failures.append(
                    "bin %d at %s: %r vs %r" % (ib, dict(zip(thetas, point)), v1, v2))
    return failures


def modifiers_for(doc, dist_name, sample_name):
    for dist in doc.get("distributions", []):
        if dist.get("name") != dist_name:
            continue
        for sample in dist.get("samples", []):
            if sample.get("name") == sample_name:
                return sample.get("modifiers", [])
    return []


# --------------------------------------------------------------------------------------------
# Case 1: one process, one `shape` nuisance, no lnN at all.
#
# ModelTools.py:984 leaves n_exp_bin<b>_proc_<p> a plain RooRealVar when a process has no rate
# systematics, so the coefficient here is RooProduct(RooRealVar, AsymPow) with no
# ProcessNormalization anywhere -- the decomposition must not require one.
# --------------------------------------------------------------------------------------------
def case_minimal_pairing():
    w = ROOT.RooWorkspace("w")
    x = ROOT.RooRealVar("x", "x", XLO, XHI)
    x.setBins(NBINS)
    theta = make_theta(w, "shape_syst")

    nom = th1("nom", [10.0, 20.0, 30.0, 40.0])
    # Deliberately different integrals from the nominal, so the AsymPow is non-trivial and the
    # yield effect really does have to travel in the normsys rather than the histosys.
    hi = th1("hi", [12.0, 23.0, 31.0, 44.0])
    lo = th1("lo", [9.0, 18.0, 29.0, 37.0])

    pdf = make_shape_pdf("procA", x, nom, [(theta, hi, lo)])
    rate = ROOT.RooRealVar("n_exp_procA", "n_exp_procA", nom.Integral())
    ap = asym_pow("systeff_procA_shape_syst", theta,
                  lo.Integral() / nom.Integral(), hi.Integral() / nom.Integral())
    coef = ROOT.RooProduct("n_exp_final_procA", "n_exp_final_procA", ROOT.RooArgList(rate, ap))
    add = ROOT.RooAddPdf("chan", "chan", ROOT.RooArgList(pdf), ROOT.RooArgList(coef))
    _keepalive.extend([x, rate, coef, add])
    getattr(w, "import")(add, ROOT.RooFit.RecycleConflictNodes())

    w2, doc = roundtrip(w)
    failures = compare(w, w2, "chan", ["shape_syst"], [(0.0,), (1.0,), (-1.0,), (2.0,), (-2.0,)])

    mods = modifiers_for(doc, "chan", "procA")
    histosys = [m for m in mods if m["type"] == "histosys" and m["parameter"] == "shape_syst"]
    normsys = [m for m in mods if m["type"] == "normsys" and m["parameter"] == "shape_syst"]
    if len(histosys) != 1 or len(normsys) != 1:
        failures.append("expected exactly one histosys and one normsys for the shape theta, got "
                        "%d and %d" % (len(histosys), len(normsys)))
    else:
        if histosys[0].get("constraint") != "shape_syst_Pdf":
            failures.append("histosys constraint is %r, expected 'shape_syst_Pdf'"
                            % histosys[0].get("constraint"))
        if histosys[0].get("constraint") != normsys[0].get("constraint"):
            failures.append("paired histosys/normsys disagree on constraint: %r vs %r"
                            % (histosys[0].get("constraint"), normsys[0].get("constraint")))

    coef2 = w2.pdf("chan").coefList().at(0)
    if not isinstance(coef2, ROOT.RooProduct):
        failures.append("reimported coefficient is %s, expected RooProduct" % coef2.ClassName())
    elif not any(isinstance(c, ROOT.AsymPow) for c in coef2.components()):
        failures.append("reimported RooProduct has no AsymPow factor")
    return failures


# --------------------------------------------------------------------------------------------
# Case 2: the same theta is a `shape` systematic on procA and only an lnN on procB.
#
# procB must come back as a bare ProcessNormalization (asymmThetaList entry, no external AsymPow),
# and procA must be unaffected by it. procA also carries a second theta with a histosys but no
# AsymPow -- the pure-shape-morph case.
# --------------------------------------------------------------------------------------------
def case_shared_theta_no_interference():
    w = ROOT.RooWorkspace("w")
    x = ROOT.RooRealVar("x", "x", XLO, XHI)
    x.setBins(NBINS)
    shared = make_theta(w, "shared_syst")
    shape_only = make_theta(w, "shape_only_syst")

    nom_a = th1("nomA", [10.0, 20.0, 30.0, 40.0])
    hi_a = th1("hiA", [12.0, 23.0, 31.0, 44.0])
    lo_a = th1("loA", [9.0, 18.0, 29.0, 37.0])
    # Equal integrals to the nominal: getExtraNorm would drop the AsymPow here, so this theta gets
    # a histosys and no normsys at all.
    hi_b = th1("hiB", [11.0, 21.0, 29.0, 39.0])
    lo_b = th1("loB", [9.0, 19.0, 31.0, 41.0])

    pdf_a = make_shape_pdf("procA", x, nom_a, [(shared, hi_a, lo_a), (shape_only, hi_b, lo_b)])
    rate_a = ROOT.RooRealVar("n_exp_procA", "n_exp_procA", nom_a.Integral())
    ap = asym_pow("systeff_procA_shared", shared,
                  lo_a.Integral() / nom_a.Integral(), hi_a.Integral() / nom_a.Integral())
    coef_a = ROOT.RooProduct("n_exp_final_procA", "n_exp_final_procA", ROOT.RooArgList(rate_a, ap))

    nom_b = th1("nomB", [40.0, 30.0, 20.0, 10.0])
    pdf_b = make_shape_pdf("procB", x, nom_b, [])
    coef_b = ROOT.ProcessNormalization("n_exp_procB", "n_exp_procB", nom_b.Integral())
    coef_b.addAsymmLogNormal(0.92, 1.11, shared)

    add = ROOT.RooAddPdf("chan", "chan", ROOT.RooArgList(pdf_a, pdf_b),
                         ROOT.RooArgList(coef_a, coef_b))
    _keepalive.extend([x, rate_a, coef_a, coef_b, add])
    getattr(w, "import")(add, ROOT.RooFit.RecycleConflictNodes())

    w2, doc = roundtrip(w)
    thetas = ["shared_syst", "shape_only_syst"]
    points = [(0.0, 0.0), (1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0), (1.0, -1.0),
              (-1.5, 2.0)]
    failures = compare(w, w2, "chan", thetas, points)

    mods_b = modifiers_for(doc, "chan", "procB")
    if any(m["type"] == "histosys" for m in mods_b):
        failures.append("procB should have no histosys")
    if not any(m["type"] == "normsys" and m["parameter"] == "shared_syst" for m in mods_b):
        failures.append("procB lost its lnN normsys for the shared theta")

    mods_a = modifiers_for(doc, "chan", "procA")
    shape_only_mods = [m["type"] for m in mods_a if m["parameter"] == "shape_only_syst"]
    if shape_only_mods != ["histosys"]:
        failures.append("shape-only theta on procA produced %r, expected ['histosys']"
                        % shape_only_mods)

    coefs2 = w2.pdf("chan").coefList()
    by_name = {coefs2.at(i).GetName(): coefs2.at(i) for i in range(coefs2.getSize())}
    b = [v for k, v in by_name.items() if k.endswith("procB")]
    if not b:
        failures.append("procB coefficient missing after reimport")
    elif not isinstance(b[0], ROOT.ProcessNormalization):
        failures.append("procB coefficient is %s, expected a bare ProcessNormalization"
                        % b[0].ClassName())
    return failures


# --------------------------------------------------------------------------------------------
# Case 3: a RooProduct coefficient holding something the decomposition does not understand.
# The export must decline cleanly rather than emit a half-written or garbage node.
# --------------------------------------------------------------------------------------------
def case_unknown_product_factor_declines():
    w = ROOT.RooWorkspace("w")
    x = ROOT.RooRealVar("x", "x", XLO, XHI)
    x.setBins(NBINS)
    theta = make_theta(w, "shape_syst")

    nom = th1("nomC", [10.0, 20.0, 30.0, 40.0])
    hi = th1("hiC", [12.0, 23.0, 31.0, 44.0])
    lo = th1("loC", [9.0, 18.0, 29.0, 37.0])

    pdf = make_shape_pdf("procA", x, nom, [(theta, hi, lo)])
    rate = ROOT.RooRealVar("n_exp_procA", "n_exp_procA", nom.Integral())
    # Neither a ProcessNormalization, nor an AsymPow, nor a plain variable.
    weird = ROOT.RooFormulaVar("weird", "weird", "1.0 + 0.1*@0", ROOT.RooArgList(theta))
    coef = ROOT.RooProduct("n_exp_final_procA", "n_exp_final_procA", ROOT.RooArgList(rate, weird))
    add = ROOT.RooAddPdf("chan", "chan", ROOT.RooArgList(pdf), ROOT.RooArgList(coef))
    _keepalive.extend([x, rate, weird, coef, add])
    getattr(w, "import")(add, ROOT.RooFit.RecycleConflictNodes())

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        path = f.name
    try:
        ok = ROOT.RooJSONFactoryWSTool(w).exportJSON(path)
    except Exception:
        return []  # refusing outright is a clean decline

    if not ok:
        return []

    with open(path) as fh:
        doc = json.load(fh)
    for dist in doc.get("distributions", []):
        if dist.get("name") == "chan" and dist.get("type") == "histfactory_dist":
            return ["export produced a histfactory_dist for a coefficient it cannot decompose"]
    return []


# --------------------------------------------------------------------------------------------
# Case 4: one theta carrying *both* a shape effect and an lnN on the same process, i.e. two
# normsys modifiers sharing one parameter alongside a histosys. Exactly one may be promoted to an
# AsymPow; the rest must stay in the ProcessNormalization. The factors commute, so the split is
# faithful either way -- this checks that the importer does not double-count or drop one.
# --------------------------------------------------------------------------------------------
def case_theta_with_shape_and_lnn():
    w = ROOT.RooWorkspace("w")
    x = ROOT.RooRealVar("x", "x", XLO, XHI)
    x.setBins(NBINS)
    theta = make_theta(w, "both_syst")

    nom = th1("nomD", [10.0, 20.0, 30.0, 40.0])
    hi = th1("hiD", [12.0, 23.0, 31.0, 44.0])
    lo = th1("loD", [9.0, 18.0, 29.0, 37.0])

    pdf = make_shape_pdf("procA", x, nom, [(theta, hi, lo)])
    pn = ROOT.ProcessNormalization("n_exp_procA", "n_exp_procA", nom.Integral())
    pn.addAsymmLogNormal(0.95, 1.07, theta)  # the lnN half
    ap = asym_pow("systeff_procA_both", theta,
                  lo.Integral() / nom.Integral(), hi.Integral() / nom.Integral())
    coef = ROOT.RooProduct("n_exp_final_procA", "n_exp_final_procA", ROOT.RooArgList(pn, ap))
    add = ROOT.RooAddPdf("chan", "chan", ROOT.RooArgList(pdf), ROOT.RooArgList(coef))
    _keepalive.extend([x, pn, coef, add])
    getattr(w, "import")(add, ROOT.RooFit.RecycleConflictNodes())

    w2, doc = roundtrip(w)
    failures = compare(w, w2, "chan", ["both_syst"],
                       [(0.0,), (1.0,), (-1.0,), (2.0,), (-2.0,), (0.4,)])

    mods = modifiers_for(doc, "chan", "procA")
    n_normsys = len([m for m in mods if m["type"] == "normsys" and m["parameter"] == "both_syst"])
    if n_normsys != 2:
        failures.append("expected two normsys for the doubled theta, got %d" % n_normsys)
    return failures


def main():
    cases = [
        ("minimal histosys+normsys pairing", case_minimal_pairing),
        ("shared theta, shape vs lnN", case_shared_theta_no_interference),
        ("theta with both a shape and an lnN", case_theta_with_shape_and_lnn),
        ("unknown RooProduct factor declines", case_unknown_product_factor_declines),
    ]
    rc = 0
    for name, fn in cases:
        try:
            failures = fn()
        except Exception as exc:  # noqa: BLE001 - report and keep going
            print("ERROR  %s: %s" % (name, exc))
            rc = 1
            continue
        if failures:
            rc = 1
            print("FAIL   %s" % name)
            for f in failures:
                print("         %s" % f)
        else:
            print("PASS   %s" % name)
    return rc


if __name__ == "__main__":
    sys.exit(main())
