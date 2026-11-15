"""Every number, table and figure of the paper (census, revision 10) and of its supplement, from released outputs.

  nice -n 10 python3 -I make_numbers_v10.py [--draft]

Inputs (checked before anything is written; --draft fills missing macros with '??' for a trial compilation only):
  from the repository ($QV_CERTIFY, default ~/quantavane-spherical-codes/certify):
    census/index.tsv (snapshot.py), paper_v8/v8_numbers.json and paper_v8/v8_verdicts.tsv (make_numbers_v8.py,
    released; the version 8 decisions of record), paper_v8/verify/out_exact*/ (exactcheck.py), paper_v8/verify/out_strict7
    and out_sym7, paper_v8/above/out and out_more (run_above.py), paper_v8/recount/out_flex8 and out_flex8r
    (flexexact8.py, flexrefine8.py), paper_v8/recount/tolscan_close.jsonl, paper_v8/live/out (run_live.py),
    v5/live/live_20261107.tsv (the live table of 7 November 2026), ../codes/MANIFEST.tsv (the author's uploads), and
    paper_v8's figure fig_cost_v8.pdf and table v8_tab_quant.tex (copied unchanged into the supplement);
  from this folder (certify/paper_v9/ of the repository): recog/out_all (recog9.py), recog/sample40.txt,
    recog/out_strict9 and recog/out_strict9_ctl (strict9.py), prestress/out (prestress9.py), cut9/out (cut9.py),
    ranks/out (flintrank9.py), ranks/pari_check9.json (pari_check9.py), certdata/verify.json and export_log.json
    (export9.py, verify_cert9.py), exactdata/ (strict9.py --data, exactdata_above9.py), example/ (example.py).
Writes v10_numbers.tex, v10_numbers.json, the tables v10_tab_*.tex (main text) and v10s_tab_*.tex (supplement), the
figures fig_verdicts_v10.pdf and fig_cost_v10s.pdf, the per-entry file v10_verdicts.tsv and VERIFY_numbers_v10.md.

The version 9 classes are those of version 8 (v8_verdicts.tsv) with these changes, made here and nowhere else:
  - an open hyperstatic entry ('local') that strict9.py proves STRICT becomes strict (certificate strict_recog9);
  - an above-limit entry that prestress9.py proves STRICT becomes strict (strict_second), and one with a second-order
    descent certificate becomes not a local minimum (notlocal_descent);
  - the negative certificate takes the larger of the radii of version 8 and of the intrinsic cut (cut9.py, decided by
    flexexact8.check); an entry is resolved when that radius is set by the separation and is at least 1e-12;
  - for the main tables, 'unresolved' and the two equal-contact non-minima ('equalcontact_notlocal') count as no verdict.
"""
import collections, glob, gzip, json, math, os, re, shutil, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
P8 = os.path.join(CERT, "paper_v8")
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, os.path.join(CERT, "census"))
sys.path.insert(0, os.environ.get("QV_MPLLIB", os.path.expanduser("~/Library/Caches/quantavane-census/mpllib")))
import numpy as np
import snapshot

DRAFT = "--draft" in sys.argv
RANGES = [(3, 5), (6, 10), (11, 17), (18, 24), (25, 32)]
M, SRC = {}, {}


# ============================================================================================ formatting (as version 8)
def put(name, value, src):
    assert name.isalpha(), "a LaTeX macro name must be letters only: " + name
    M[name] = value
    SRC[name] = src


def num(k):
    return f"{k:,}".replace(",", "{,}")


def sci(x, digits=1):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "--"
    if x == 0:
        return "$0$"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, digits) >= 10:
        m, e = m / 10, e + 1
    if 0 <= e <= 2:
        return f"${x:.{max(digits - e, 0)}f}$"
    return f"${m:.{digits}f}\\cdot10^{{{e}}}$"


def sci_up(x, digits=1):
    """An upper bound, rounded up in the printed digit."""
    if x is None or x == 0:
        return sci(x, digits)
    e = int(math.floor(math.log10(abs(x))))
    m = math.ceil(x / 10 ** e * 10 ** digits - 1e-9) / 10 ** digits
    if m >= 10:
        m, e = m / 10, e + 1
    if 0 <= e <= 2:
        return sci(m * 10 ** e, digits)
    return f"${m:.{digits}f}\\cdot10^{{{e}}}$"


def sci_dn(x, digits=1):
    """A lower bound, rounded down in the printed digit."""
    if x is None or x == 0:
        return sci(x, digits)
    e = int(math.floor(math.log10(abs(x))))
    m = math.floor(x / 10 ** e * 10 ** digits + 1e-9) / 10 ** digits
    if 0 <= e <= 2:
        return sci(m * 10 ** e, digits)
    return f"${m:.{digits}f}\\cdot10^{{{e}}}$"


def quant(vals, ps=(0, 0.1, 0.5, 0.9, 1)):
    """Order statistics: the value of rank floor(p (n - 1)); the median of an even number of values is the lower one."""
    v = sorted(vals)
    if not v:
        return [None] * len(ps)
    return [v[int(math.floor(p * (len(v) - 1) + 1e-12))] for p in ps]


def med(vals):
    return quant(vals, (0.5,))[0]


def pct(a, b, digits=0):
    if not b:
        return "--"
    v = 100.0 * a / b
    return f"{v:.{digits}f}"


def elist(ks):
    return ", ".join(f"$({d},{n})$" for d, n in sorted(ks))


def load_dir(path, pattern="*.json"):
    out = {}
    for f in glob.glob(os.path.join(path, pattern)):
        j = json.load(open(f))
        out[(j["d"], j["n"])] = j
    return out


def jsonl(path):
    return [json.loads(l) for l in open(path) if l.strip()]


def band(d):
    for i, (a, b) in enumerate(RANGES):
        if a <= d <= b:
            return i


# ============================================================================================ inputs
REQ_CERT = ["census/index.tsv", "paper_v8/v8_numbers.json", "paper_v8/v8_verdicts.tsv", "paper_v8/verify/out_exact",
            "paper_v8/verify/out_exact2", "paper_v8/verify/out_strict7", "paper_v8/verify/out_sym7", "paper_v8/above/out",
            "paper_v8/above/out_more", "paper_v8/recount/out_flex8", "paper_v8/recount/out_flex8r",
            "paper_v8/recount/tolscan_close.jsonl", "paper_v8/live/out", "v5/live/live_20261107.tsv",
            "paper_v8/fig_cost_v8.pdf", "paper_v8/v8_tab_quant.tex"]
REQ_HERE = ["recog/out_all", "recog/out_pass2", "recog/sample40.txt", "recog/ctl_plan.json", "recog/out_strict9", "recog/out_strict9_ctl", "prestress/out",
            "cut9/out", "ranks/out", "ranks/pari_check9.json", "certdata/verify.json", "certdata/export_log.json",
            "exactdata", "example/example_d5_n16.json"]


def check_inputs():
    missing = [r for r in REQ_CERT if not os.path.exists(os.path.join(CERT, r))]
    missing += [r for r in REQ_HERE if not os.path.exists(os.path.join(HERE, r)) or
                (os.path.isdir(os.path.join(HERE, r)) and not os.listdir(os.path.join(HERE, r)))]
    if not os.path.exists(os.path.join(os.path.dirname(CERT), "codes", "MANIFEST.tsv")):
        missing.append("codes/MANIFEST.tsv")
    if missing and not DRAFT:
        sys.exit("missing inputs, nothing written: " + ", ".join(missing))
    return missing


def maybe_dir(rel, pattern="*.json"):
    p = os.path.join(HERE, rel)
    return load_dir(p, pattern) if os.path.isdir(p) else {}


def main():
    missing = check_inputs()
    M8 = json.load(open(os.path.join(P8, "v8_numbers.json")))
    S8 = "paper_v8/v8_numbers.json (make_numbers_v8.py, released; unchanged in version 9)"
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(P8, "v8_verdicts.tsv")) if not l.startswith("#")]
    head8 = rows[0]
    R8 = {(int(r[0]), int(r[1])): dict(zip(head8, r)) for r in rows[1:]}
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index() if e["listed"]}
    keys = sorted(k for k in idx if 3 <= k[0] <= 32 and k[1] <= 1024)
    assert sorted(R8) == keys, "v8_verdicts.tsv does not cover the census"
    S_TSV = "paper_v8/v8_verdicts.tsv (released per-entry verdicts of version 8)"

    # ------------------------------------------------------------------- pass-through macros (meaning unchanged)
    PASS8 = """AboveGanContacts AboveGanDegs AboveGanDist AboveGanMB AboveGanNext AboveGanRank AboveGanSec AboveMaxCols
    AboveMaxColsDone AboveNonStrictClause BigHyperAll BigHyperGood BigIsoAll BigIsoGood BigIsoMaxU BigStrictAll
    BigStrictGood CJKTList CPUHoursCensus CPUMedLoc CPUMedPass CloseRun CloseSentence DiagRowsumMed
    ExBetaK ExCDigits ExCPU ExDegrees ExDist
    ExEpsGamma ExEpsMax ExGap ExKappa ExKappaHalf ExLambda ExMuP ExNContacts ExNextGap ExNormC ExOnlyDistCoordMax
    ExOnlyDistLSMax ExOnlyN ExR ExRho ExRu ExSpan ExStressMin ExStressRel ExTableCos ExTheta ExU ExZeta ExactCPUMax
    FXEDenseOps FXEHours FXELimbOps FXEMaxCols FXEMaxRows FXEMaxSec FXEQR FXERatioMax FXERatioMed FXERatioMin
    FXERefineList FamilyList FlatColoursMax FlatColoursMin FlatList FlatLiveImpRange FlatOtherList FlatUnionList GapCertExistFail
    GapCertMaxU GapCertRun GapFound GapFoundSmall GapScanN ImpLowerA ImpLowerB LiveCertBeaten LiveCertListed LiveDate
    LiveDistLowMedBeaten LiveDistLowMinBeaten LiveFetchFlex LiveFetchPos LiveListed LiveRerunBeatenN LiveRerunBeatenPos
    LiveRerunFlexFlex LiveRerunFlexNone LiveRerunFlexPos LiveRerunFlexPosPct LiveRerunPosFlex LiveRerunPosOther
    LiveRerunPosPos LiveRerunPosStrictIso LiveStrictImpMed LocExtraMax LocExtraMed LocExtraMin LowCondFamily
    MaxUnknownsLoc MaxUnknownsPass NAboveAll NAboveCore NAboveDeg NAboveField NAboveFlat NAboveFlatUnion NAboveLPNone
    NAboveMoreRun NAboveMoreStopped NAboveNoStress NAboveNoStressZero NAboveNotRigid NAboveNotRigidList
    NAboveNotRigidStress NAboveSecondPass NAboveStrict NAboveTooLarge NAboveTried NAttFailExist NAttFailNoStress
    NAttFailOther NAttemptFail NCensus NCtlRun NCtlStrict NDecidedExact NDiagClose NDiagSep NDiagSepSing NDispBelow
    NExactRecognised NExactRootInside NExactStrict NExactTested NExactTestedHyper NExactTestedHyperRecog NExactTestedIso
    NExactTried NFlagAll NFlatLive NFlatLiveImproved NFlexCensus NListed NLocal NLocalAll NNotLocalFlat NNotLocalNonopt
    NOpenDegree NOpenExact NOpenNoPoly NOpenNotInField NOpenTooMany NOptAll NPassRattlers NPolyCensus NPolyListed
    NPropNever NRattlerMoved NStrictAbove NStrictHyper NStrictHyperFieldHigh NStrictHyperFieldOne NStrictHyperFieldTwo
    NStrictIso NStrictSym NSurvey NSymNoSym NSymNontrivial NSymSquare NSymStrict NUnflagAll NewColsA NewColsB NewGapA
    NewRankA NewRankB NonoptDistA NonoptDistB NonoptStressA NonoptStressB OpenContacts OpenLPMaxMin OpenRankNeeded
    OpenRankStressed OpenStressed OptAboveText PctPropNever RattlerMovedMax SecondFlexGood SecondFlexN SizeCap
    SplitCoreA SplitCoreB SplitDistMax SplitRankA SplitRankB StrictDistMax StrictDistMed StrictNaNList StrictNewList
    StructShareLarge StructShareMid StructShareSmall SymExEntry SymExExtra SymExGroup SymExRed SymExU
    SymGroupMax SymGroupMin SymRedMax TiltMax TiltMaxEntry TiltMin TiltMinEntry TolConflictSentence TolOtherEntries
    TolOtherProposals TolOtherSentence ValChiMax XCErrMax XCPUSlope XCSame XChiMax XDistRatioMax XDistRatioMed XDistTinyN
    XGapHypEleven XGapHypNine XGapIsoEleven XGapIsoNine XGapMin XGapMinEntry XGapMinExact XGapMinTol XHours XInsideHyp
    XInsideIso XKappaHalfMed XKappaRatioMax XKappaRatioMin XLocDistTubeMax XMaxU XNChecked XNLocOK XNNonoptOK XNPassOK
    XNRhoDefined XNoRhoIso XOutsideIso XRUniqueMed XRepeatList XRhoMed XSecMax XStageTwoList XStageTwoN
    XStressRatioMax XStressRatioMin XStrictOK XStrictRun XSymOK XSymRun XTubeMax XUniqueCovers ZetaN ZetaRatioMax
    ZetaRatioMed ZetaRatioMin""".split()
    for k in PASS8:
        put(k, M8[k], S8)
    # minor 10: the smallest ratio tau'/tau, printed to three digits rounded down (version 8 printed 1.00)
    fl8 = {**load_dir(os.path.join(P8, "recount", "out_flex8")), **load_dir(os.path.join(P8, "recount", "out_flex8r"))}
    fl8 = {k: j for k, j in fl8.items() if j.get("status") == "FLEX"}
    ratios = [j["tau_prime"] / j["tol"] for j in fl8.values()]
    put("FlexTauRatioMin", f"{math.floor(min(ratios) * 1000) / 1000:.3f}",
        "paper_v8/recount/out_flex8*/: min tau_prime/tol over the negative certificates, rounded down to 3 digits")
    put("NTauBelowTol", num(sum(1 for r in ratios if r < 1)), "paper_v8/recount/out_flex8*/: tau_prime < tol")

    # ------------------------------------------------------------------- the worked example (5,16), regenerated (minor 9)
    def dec(x, sig=2, how="up"):
        from decimal import Decimal as D_, ROUND_CEILING, ROUND_FLOOR
        if not (1e-3 <= abs(x) < 1e3):
            return (sci_up if how == "up" else sci_dn)(x, sig - 1)
        q = D_(repr(x))
        step = D_(1).scaleb(q.adjusted() - sig + 1)
        return "$" + format(q.quantize(step, rounding=ROUND_CEILING if how == "up" else ROUND_FLOOR), "f") + "$"
    exf = os.path.join(HERE, "example", "example_d5_n16.json")
    if os.path.exists(exf):
        S_E2 = "example/example_d5_n16.json (paper_v8/example/example.py 5 16, regenerated in version 9)"
        E2 = json.load(open(exf))
        x2 = E2["exactcheck"]
        rec516 = [j for j in (json.load(open(f)) for f in glob.glob(os.path.join(P8, "verify", "out_exact", "d5_n16_LOCALMIN_*.json")))][0]
        for key in ("chi", "tube", "beta_K", "kappa", "rho", "dist_ls", "stress_min_rel"):
            assert x2[key] == rec516[key], ("the regenerated example differs from the decision of record", key)
        put("ExBNContacts", num(E2["contacts"]), S_E2)
        put("ExBExtra", num(E2["extra"]), S_E2)
        put("ExBIso", num(E2["iso_count"]), S_E2)
        put("ExBU", num(E2["unknowns"]), S_E2)
        put("ExBBetaK", sci_up(x2["beta_K"]), S_E2 + ": beta_K")
        put("ExBStressRel", dec(x2["stress_min_rel"], 3, "dn"), S_E2 + ": smallest relative stress")
        put("ExBChi", sci_up(x2["chi"]), S_E2 + ": chi")
        put("ExBKappa", dec(x2["kappa"], 3, "dn"), S_E2 + ": kappa")
        put("ExBRho", dec(x2["rho"], 3, "dn"), S_E2 + ": rho")
        put("ExBTube", sci_up(x2["tube"]), S_E2 + ": tube")
        put("ExBDist", sci_up(x2["dist_ls"]), S_E2 + ": dist_ls")

    # ------------------------------------------------------------------- per-entry records of version 8
    ex = {}
    for sub in ("out_exact", "out_exact2"):
        for f in glob.glob(os.path.join(P8, "verify", sub, "*.json")):
            j = json.load(open(f))
            ex[(j["d"], j["n"], j["status"])] = j            # out_exact2 overrides for its three records
    s7 = load_dir(os.path.join(P8, "verify", "out_strict7"))
    sym7 = load_dir(os.path.join(P8, "verify", "out_sym7"))
    above = load_dir(os.path.join(P8, "above", "out"))
    for k, j in load_dir(os.path.join(P8, "above", "out_more")).items():
        if str(j.get("status", "")).startswith(("STRICT", "LOCAL", "NOT LOCAL")):
            above[k] = dict(above[k], **j, second_pass=True)

    # ------------------------------------------------------------------- version 9 outputs
    rec_all = maybe_dir("recog/out_all")
    st9 = maybe_dir("recog/out_strict9")
    st9c = maybe_dir("recog/out_strict9_ctl")
    ps9 = maybe_dir("prestress/out")
    cut = maybe_dir("cut9/out")
    fr9 = maybe_dir("ranks/out")

    # ------------------------------------------------------------------- the version 9 classes
    cls = {k: [R8[k]["verdict"], R8[k]["certificate"]] for k in keys}
    for k, j in st9.items():
        if j.get("status") == "STRICT":
            assert cls[k][0] == "local", k
            cls[k] = ["strict", "strict_recog9"]
    for k, j in ps9.items():
        st = str(j.get("status", ""))
        if st.startswith("STRICT"):
            cls[k] = ["strict", "strict_second"]
        elif st.startswith("NOT A LOCAL MINIMUM"):
            cls[k] = ["notlocal", "notlocal_descent"]
    # the negative certificate: the larger radius of version 8 and of the cut
    neg = {}
    for k in keys:
        if R8[k]["verdict"] in ("flex", "unresolved") and k in fl8:
            j8 = fl8[k]
            d8, by8 = j8["delta"], j8.get("delta_by")
            dcut, bycut, tcut = None, None, None
            jc = cut.get(k)
            if jc and jc.get("status") == "FLEX":
                e = jc["exact"]
                dcut, bycut, tcut = e["delta"], e["delta_by"], e["tau_prime"]
            # a certificate is resolved when its radius is set by the separation and is at least 1e-12; the version 9
            # radius is the largest resolved one (the cut never demotes a resolved verdict of version 8)
            res8 = R8[k]["verdict"] == "flex"
            rescut = dcut is not None and dcut >= 1e-12 and bycut == "tau'"
            cand = ([d8] if res8 else []) + ([dcut] if rescut else [])
            res9 = bool(cand)
            best = max(cand) if cand else max(d8, dcut or 0.0)
            from_cut = rescut and dcut >= (d8 if res8 else 0.0)
            neg[k] = dict(d8=d8, by8=by8, res8=res8, dcut=dcut, bycut=bycut, tcut=tcut, d9=best, from_cut=from_cut,
                          res9=res9, cut_done=k in cut, flagged=idx[k]["flagged"])
            cls[k] = ["flex", "flex_cut9" if from_cut else "flex"] if res9 else ["unresolved", "flex_unres"]

    def main_class(k):
        v = cls[k][0]
        if v in ("unresolved", "equalcontact_notlocal", "none"):
            return "none"
        return v

    C = collections.Counter(main_class(k) for k in keys)
    put("NStrictV", num(C["strict"]), "v9 classes: strict")
    put("NLocalV", num(C["local"]), "v9 classes: local (open hyperstatic)")
    put("NNonStrictV", num(C["nonstrict"]), "v9 classes: nonstrict")
    put("NNotLocalV", num(C["notlocal"]), "v9 classes: notlocal (flat cores; second-order descent)")
    put("NNegResV", num(C["flex"]), "v9 classes: resolved negative certificate (version 8 or the cut)")
    put("NNoVerdictV", num(C["none"]), "v9 classes: no verdict (incl. unresolved and equal-contact non-minima)")
    put("NUnresolvedV", num(sum(1 for k in keys if cls[k][0] == "unresolved")), "v9 classes: unresolved")
    put("NStrictHyperV", num(sum(1 for k in keys if cls[k][1] in ("strict_hyper", "strict_recog9"))),
        "v9 classes: strict by recognition (version 8 and strict9.py)")
    put("NStrictExactAllV", num(sum(1 for k in keys if cls[k][1] in ("strict_exact", "strict_new", "strict_above",
                                                                       "strict_second"))),
        "v9 classes: strict by exact arithmetic (strict_exact, strict_new, strict_above, strict_second)")
    put("NRecogNewStrict", num(sum(1 for k in keys if cls[k][1] == "strict_recog9")), "recog/out_strict9: STRICT")

    # ------------------------------------------------------------------- distances, tiers, neighbourhoods
    def exrec(k, status):
        return ex.get((k[0], k[1], status))

    # Tier A for a hyperstatic strict certificate (M2 of the review of v9): the exact equal-contact configuration P* of the
    # strictness certificate (Theorem 3.13 or 3.15), rotated to the reference vectors xi by Eckart's condition, is a zero
    # of the square system F of the hyperstatic certificate (unit norms, the contacts of E' at c*, psi = 0). If it lies in
    # the box |z - z~| <= r_unique where the Krawczyk test (K) holds, it is the certificate's zero P', so chi = 0 at P' and
    # Theorem 3.8 with the certificate's omega_- and zeta gives rho* >= rho: P* is the only minimiser in {D < rho*}, which
    # contains the stored code when the certificate's D < rho. Bounds used: |R P* - xi|_inf <= D(xi, P*) <= D(stored, P*)
    # + sqrt(nV d) 2^-50 (xi is the stored code rounded and normalised in double precision); |xi - p~|_inf <= upsilon, the
    # certificate's bound max|xi - p'| + r, recovered as theta / (2 nV mu_p) with a relative margin of 2^-40; and
    # |c* - c~| <= |c* - c_mid| + 10^-39 + (width of c_interval_30) + r.
    s7b = load_dir(os.path.join(P8, "verify", "out_strict7b"))
    same_zero_log = {}

    def same_zero(k, c, e):
        src = None
        if c == "strict_hyper":
            src = next((j for j in (s7.get(k), s7b.get(k)) if j and str(j.get("status")).startswith("STRICT")), None)
        elif c == "strict_sym":
            src = sym7.get(k)
        elif c == "strict_recog9":
            src = st9.get(k)
        if not src or not e or not e.get("r_unique"):
            return False
        nV, d = int(e["nV"]), int(e["d"])
        dP = src.get("dist_ls_upper")
        if dP is None or not math.isfinite(float(dP)):
            dP = math.sqrt(nV * d) * float(src["dist_coord_upper"]) * (1 + 2 ** -40)
        ups = float(e["theta"]) / (2 * nV * float(e["mu_p"])) * (1 + 2 ** -40)
        s_p = float(dP) + math.sqrt(nV * d) * 2 ** -50 + ups
        lo, hi = (Fraction(int(x), 10 ** 30) for x in e["c_interval_30"])
        cs = Fraction(src["gamma_mid"]) if c == "strict_sym" else Fraction(str(src["c_star"]))
        s_c = abs(cs - (lo + hi) / 2) + Fraction(1, 10 ** 39) + (hi - lo) + Fraction(float(e["r"]))
        ok = s_p <= float(e["r_unique"]) and s_c <= Fraction(float(e["r_unique"]))
        same_zero_log[k] = dict(ok=ok, ratio=s_p / float(e["r_unique"]))
        return ok

    info = {}
    for k in keys:
        v, c = cls[k]
        r = R8[k]
        dist = float(r["dist_ls"]) if r["dist_ls"] not in ("", "--") else None
        rho = kap = inside = None
        if c == "strict_iso":
            e = exrec(k, "PASS")
            rho, kap, inside = e.get("rho"), e.get("kappa_half"), e.get("code_inside_rho")
        elif c in ("strict_hyper", "strict_sym", "strict_recog9", "local"):
            e = exrec(k, "LOCALMIN")
            if e:
                rho, kap, inside = e.get("rho"), e.get("kappa_half"), e.get("code_inside_rho")
        if c == "strict_recog9":
            dist = st9[k].get("dist_ls_upper")
        if c in ("strict_second", "notlocal_descent"):
            dist = ps9[k].get("dist_ls_upper")
        tier = None
        if v in ("strict", "local", "nonstrict"):
            if dist is not None and dist <= 1e-12:
                # tier A only for strict verdicts: isostatic (chi = 0 at P' itself), or hyperstatic with P' = P*
                okA = v == "strict" and bool(inside) and (c == "strict_iso" or
                                                          (c in ("strict_hyper", "strict_sym", "strict_recog9")
                                                           and same_zero(k, c, exrec(k, "LOCALMIN"))))
                tier = "A" if okA else "B"
            else:
                tier = "C"
        info[k] = dict(dist=dist, rho=rho, kappa_half=kap, inside=inside, tier=tier)
    strict_k = [k for k in keys if cls[k][0] == "strict"]
    dists = [info[k]["dist"] for k in strict_k if info[k]["dist"] is not None]
    put("DistStrictMax", sci_up(max(dists)), "max certified distance over the strict verdicts (v8_verdicts.tsv dist_ls; strict9, prestress9)")
    put("DistStrictMin", sci_up(min(dists)), "min certified distance over the strict verdicts (an upper bound, rounded up)")
    put("NStrictDistSmall", num(sum(1 for x in dists if x <= 1e-12)), "strict verdicts with certified distance <= 1e-12")
    tiers = collections.Counter(info[k]["tier"] for k in strict_k)
    hypA = [k for k in strict_k if info[k]["tier"] == "A" and cls[k][1] != "strict_iso"]
    hyp_tried = [k for k in same_zero_log]
    put("NTierAHyper", num(len(hypA)), "tier A hyperstatic strict certificates: P* in the uniqueness box of the hyperstatic certificate")
    put("NTierAHyperFail", num(sum(1 for k in hyp_tried if not same_zero_log[k]["ok"])),
        "hyperstatic strict certificates inside rho whose P* was not shown to be the certificate's zero (tier B)")
    put("TierAHyperRatioMax", sci_up(max(same_zero_log[k]["ratio"] for k in hyp_tried)) if hyp_tried else "--",
        "largest (distance bound)/r_unique among the hyperstatic tier-A checks")
    put("NTierAIso", num(sum(1 for k in strict_k if info[k]["tier"] == "A" and cls[k][1] == "strict_iso")),
        "tier A isostatic certificates")
    put("NStrictTierA", num(tiers["A"]), "strict, distance <= 1e-12 and code inside rho (out_exact code_inside_rho)")
    put("NStrictTierB", num(tiers["B"]), "strict, distance <= 1e-12, code not inside an explicit neighbourhood")
    put("NStrictTierC", num(tiers["C"]), "strict, distance > 1e-12")
    kh = [ex[(k[0], k[1], "PASS")]["kappa_half"] for k in keys if cls[k][1] == "strict_iso" and ex[(k[0], k[1], "PASS")].get("kappa_half")]
    rh = [ex[(k[0], k[1], "PASS")]["rho"] for k in keys if cls[k][1] == "strict_iso" and ex[(k[0], k[1], "PASS")].get("rho")]
    put("XKappaHalfMin", sci_dn(min(kh)), "paper_v8/verify/out_exact: min kappa/(2 sqrt d) over the isostatic certificates")
    put("XRhoMin", sci_dn(min(rh)), "paper_v8/verify/out_exact: min rho over the isostatic certificates")
    # close calls and double certificates (minor 12)
    close = [k for k in keys if cls[k][0] in ("strict", "local") and R8[k]["gap_first_pass"] not in ("", "--")
             and float(R8[k]["gap_first_pass"]) < 1e-9]
    put("NCloseCalls", num(len(close)), "v8_verdicts.tsv: positive verdicts with gap_first_pass < 1e-9")
    dbl = {}
    for j in jsonl(os.path.join(P8, "recount", "tolscan_close.jsonl")):
        k = (j["d"], j["n"])
        for a in j.get("attempts", []):
            e = a.get("exact") or {}
            if e.get("status") in ("PASS", "LOCALMIN") and info[k]["dist"] is not None:
                dbl[k] = max(dbl.get(k, 0.0), info[k]["dist"] + float(e.get("dist_ls_plus_tube") or e["dist_ls"]))
    put("NDoubleCert", num(len(dbl)), "paper_v8/recount/tolscan_close.jsonl: entries with the larger graph certified too")
    put("DoubleCertList", elist(dbl), "the same entries")
    put("DoubleCertDistMax", sci_up(max(dbl.values())) if dbl else "--",
        "max over them of dist(stored, cert 1) + dist(stored, cert 2) (triangle inequality modulo rotations)")

    # ------------------------------------------------------------------- what was known (Table known)
    # optimal and unique up to isometry by published theorems (statements read: MT12 Thm 1, MT15 Thm 1, CdLL24 Thm 1.1
    # with its Table 1.2, Cohn-Kumar JAMS 2007 Appendix A; the others as in the earlier versions). An optimal code that
    # is unique up to isometry is a strict local minimum modulo rotations, so these strict verdicts were known.
    K1 = {(3, 12): "\\cite{FT43}", (4, 10): "\\cite{BV09}", (4, 24): "\\cite{dLLM24}", (4, 120): "\\cite{And99,BD01}",
          (7, 56): "\\cite{Lev79,BS81}", (8, 240): "\\cite{Lev79,OS79,BS81}", (22, 891): "\\cite{Lev79,CK07}",
          (3, 13): "\\cite[Thm.~1]{MT12}", (3, 14): "\\cite[Thm.~1]{MT15}",
          (16, 288): "\\cite[Thm.~1.1]{CdLL24}", (20, 56): "\\cite[Thm.~1.1]{CdLL24}",
          (21, 50): "\\cite[Thm.~1.1]{CdLL24}", (21, 77): "\\cite[Thm.~1.1]{CdLL24}",
          (5, 16): "\\cite[App.~A]{CK07J}", (6, 27): "\\cite[App.~A]{CK07J}", (21, 162): "\\cite[App.~A]{CK07J}",
          (22, 100): "\\cite[App.~A]{CK07J}", (22, 275): "\\cite[App.~A]{CK07J}", (23, 552): "\\cite[App.~A]{CK07J}"}
    # optimal, uniqueness not established by us: Robinson's paper could not be read; for (21,112) uniqueness would follow
    # from the argument of CK07J App. A and the uniqueness of the strongly regular graph (112,30,2,10), not read by us
    UNSETTLED = {(3, 24): "\\cite{Rob61}", (21, 112): "\\cite{DGS77,Lev79}"}
    cjkt = [(4, 24), (5, 40), (6, 72), (7, 126), (8, 240), (9, 306)]
    mt14 = [k for k in keys if k[0] == 3 and k[1] <= 11]
    opt = [k for k in keys if idx[k]["optimal"]]
    known = set(K1) | set(cjkt) | set(mt14)
    for k in known | set(UNSETTLED):
        assert cls[k][0] == "strict", ("a control is not strict", k)
    k4 = sorted(k for k in opt if k not in known and k not in UNSETTLED)
    assert DRAFT or not k4, ("a provably optimal entry neither known nor unsettled: the text has no proposition for it", k4)
    put("NKnownStrict", num(len(known)), "entries whose strictness or jamming was known: optimal and unique, CJKT, MT14 (d=3, N<=11)")
    put("NOptUnique", num(len(K1)), "optimal and unique by published theorems")
    put("NUnsettledOpt", num(len(UNSETTLED)), "optimal; uniqueness not established by us (not claimed as new)")
    put("UnsettledOptList", elist(UNSETTLED), "the same entries")
    put("NNotClaimedNew", num(len(known | set(UNSETTLED))), "strict verdicts not claimed as new")
    put("NNewStrictAll", num(len(strict_k) - len(known | set(UNSETTLED))), "strict verdicts minus the known and unsettled ones")
    put("NOptNewStrict", num(len(k4)), "provably optimal census entries outside the known and unsettled lists")
    put("OptNewStrictList", elist(k4) if k4 else "none", "the same entries")
    certname = {"strict_iso": "iso", "strict_hyper": "recog", "strict_recog9": "recog-hd", "strict_sym": "sym",
                "strict_exact": "exact", "strict_new": "exact", "strict_above": "exact", "strict_second": "second order"}
    certlong = {"iso": "the isostatic certificate (Theorem~\\ref{thm:check})",
                "recog": "recognition (Theorem~\\ref{thm:strict} with Proposition~\\ref{prop:gram})",
                "recog-hd": "recognition at higher degree", "sym": "symmetry (Theorem~\\ref{thm:sym})",
                "exact": "exact arithmetic (Proposition~\\ref{prop:exact})", "second order": "the second-order criterion"}
    how = collections.defaultdict(list)
    for k in k4:
        how[certname[cls[k][1]]].append(k)
    put("OptNewStrictCert", "; ".join(f"{elist(v)} by {certlong[nm]}" for nm, v in sorted(how.items())) + ".",
        "certificates of the provably optimal entries")
    TK = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{p{0.24\\linewidth}rp{0.56\\linewidth}}", "\\toprule",
          "status before this census & entries & entries (certificate) \\\\", "\\midrule"]
    def row(title, ks, refs=None):
        cell = ", ".join(f"$({d},{n})$" + (f" {refs[(d, n)]}" if refs and (d, n) in refs else "") + f" ({certname[cls[(d, n)][1]]})"
                         for d, n in sorted(ks))
        TK.append(f"{title} & {len(ks)} & {cell} \\\\")
    row("optimal and unique (published)", sorted(K1), K1)
    row("jammed~\\cite{CJKT11}", [k for k in cjkt if k not in K1])
    row("$d=3$, $N\\le11$ (contact graphs in~\\cite{MT14}); not claimed as new", [k for k in mt14 if k not in K1])
    row("optimal; uniqueness not established by us; not claimed as new", sorted(UNSETTLED), UNSETTLED)
    if k4:
        row("optimal, uniqueness not known to us", k4)
    rest = [k for k in strict_k if k not in known and k not in UNSETTLED and k not in k4]
    TK.append(f"not known to us to be strict & {num(len(rest))} & " +
              ", ".join(f"{v} {nm}" for nm, v in sorted(collections.Counter(certname[cls[k][1]] for k in rest).items())) + " \\\\")
    TK += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_known.tex"), "w").write("\n".join(TK) + "\n")

    # ------------------------------------------------------------------- recognition at higher degree (A1)
    sample = [tuple(map(int, s.split("_"))) for s in open(os.path.join(HERE, "recog", "sample40.txt")).read().strip().split(",")] \
        if os.path.exists(os.path.join(HERE, "recog", "sample40.txt")) else []
    if rec_all:
        st = collections.Counter(j["status"] for j in rec_all.values())
        full = [k for k, j in rec_all.items() if str(j["status"]).startswith("recognised")]
        cosf = [k for k, j in rec_all.items() if j.get("c_degree")]
        put("RecogFull", num(len(full)), "recog/out_all: recognised (every inner product)")
        put("RecogFullTheta", num(sum(1 for k in full if rec_all[k]["status"] == "recognised in Q(theta)")),
            "recog/out_all: recognised in Q(theta)")
        put("RecogCosFound", num(len(cosf)), "recog/out_all: the cosine has a polynomial of degree <= 30")
        put("RecogCosOnly", num(len(cosf) - len(full)), "recog/out_all: cosine found, some inner product not")
        put("RecogCosDegMin", str(min(rec_all[k]["c_degree"] for k in cosf)), "recog/out_all: min degree of the cosine")
        put("RecogCosDegMax", str(max(rec_all[k]["c_degree"] for k in cosf)), "recog/out_all: max degree of the cosine")
        devs = [j["extra_dev_log2"] for j in rec_all.values() if j.get("extra_dev_log2") is not None]
        assert len(devs) == len(rec_all) == int(M8["NLocal"].replace("{,}", ""))
        put("ExtraDevDigits", str(int(math.floor(-max(devs) * math.log10(2)))),
            "recog/out_all: min over the open entries of -log10 |extra contact - c'| at the refined zero (digits)")
        put("RecogResidMinBits", str(min(-int(j["refine_log2_residuals"][-1]) for j in rec_all.values())),
            "recog/out_all: the largest final residual (as a power of two) over the open entries")
        put("RecogPBits", str(int(med([j["accuracy_bits"] for j in rec_all.values()]) - 16)),
            "recog/out_all: median accuracy_bits - 16 (the precision P of the lattice reduction)")
        sf = [k for k in sample if k in full]
        sc = [k for k in sample if k in cosf and k not in full]
        put("RecogSampleFull", num(len(sf)), "recog/out_all restricted to recog/sample40.txt: recognised")
        put("RecogSampleCos", num(len(sc)), "the same: cosine only")
        put("RecogSampleSentence", f"{len(sf)} of the 40 were recognised completely and {len(sc)} more had a cosine of "
            f"degree at most 30",
            "recog/out_all on the 40 sampled entries")
        put("RecogOutcomeSentence", ", ".join(f"{v} {s}" for s, v in sorted(st.items(), key=lambda t: -t[1])),
            "recog/out_all: statuses")
    if st9 or rec_all:
        ok = [k for k, j in st9.items() if j.get("status") == "STRICT"]
        degs = [st9[k]["field_degree"] for k in ok]
        put("RecogNewDegMin", str(min(degs)) if degs else "--", "recog/out_strict9: min field degree of the STRICT ones")
        put("RecogNewDegMax", str(max(degs)) if degs else "--", "recog/out_strict9: max field degree")
        put("RecogNewDistMax", sci_up(max(st9[k]["dist_ls_upper"] for k in ok)) if ok else "--",
            "recog/out_strict9: max certified distance")
        fails = collections.Counter(j.get("status") for k, j in st9.items() if j.get("status") != "STRICT")
        put("RecogFailClause", ("" if not fails else "; for the others " + ", ".join(f"{v}: {s}" for s, v in fails.items())),
            "recog/out_strict9: statuses other than STRICT")
        T = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{lrrrr@{\\qquad}lrrrr}", "\\toprule",
             "entry & field & cosine & values & distance & entry & field & cosine & values & distance \\\\", "\\midrule"]
        ok = sorted(ok)
        cells = [f"$({k[0]},{k[1]})$ & {st9[k]['field_degree']} & {st9[k].get('c_degree') or rec_all.get(k, {}).get('c_degree', '')} & "
                 f"{st9[k].get('distinct_inner_products')} & {sci_up(st9[k]['dist_ls_upper'])}" for k in ok]
        half = (len(cells) + 1) // 2
        for i in range(half):
            T.append(cells[i] + " & " + (cells[half + i] if half + i < len(cells) else "& & & &") + " \\\\")
        T += ["\\bottomrule", "\\end{tabular}"]
        open(os.path.join(HERE, "v10_tab_recog9.tex"), "w").write("\n".join(T) + "\n")
    rec2 = maybe_dir("recog/out_pass2")
    if rec_all:
        full1 = [k for k, j in rec_all.items() if str(j["status"]).startswith("recognised")]
        put("RecogFullStrict", num(sum(1 for k in full1 if st9.get(k, {}).get("status") == "STRICT")),
            "recog/out_strict9: STRICT among the first-pass recognitions")
    if rec2:
        f2 = [k for k, j in rec2.items() if str(j["status"]).startswith("recognised")]
        new2 = [k for k in f2 if st9.get(k, {}).get("status") == "STRICT"]
        d2 = [rec2[k]["field_degree"] for k in f2]
        nall = sum(1 for j in st9.values() if j.get("status") == "STRICT")
        put("RecogPassTwoSentence", f"A second pass at 4{{,}}000 bits (about 1{{,}}200 digits) on the {len(rec2)} entries with "
            f"at most 20 distinct inner products that the first left unrecognised recognised {len(f2)} more" +
            (f" ({elist(f2)}; fields of degree {min(d2)} to {max(d2)}), and " +
             (f"all {len(new2)}" if len(new2) == len(f2) else f"{len(new2)}") + " of them are proved strict, "
             f"so that {nall} codes in all are newly proved strict by recognition" if f2 else "") + ".",
            "recog/out_pass2 (recog9.py --pass2 --bits 4000) and recog/out_strict9")
        put("RecogPassTwoN", num(len(f2)), "recog/out_pass2: recognised")
    else:
        put("RecogPassTwoSentence", "", "recog/out_pass2 (not present)")
        put("RecogPassTwoN", "0", "recog/out_pass2 (not present)")
    if st9c:
        okc = sum(1 for j in st9c.values() if j.get("status") == "STRICT")
        plan = json.load(open(os.path.join(HERE, "recog", "ctl_plan.json")))
        prun = {(e["d"], e["n"]) for e in plan["entries"] if e["run"]}
        assert set(st9c) <= prun, "a control outside the plan"
        assert DRAFT or set(st9c) == prun, "controls missing"
        nleft = len(plan["entries"]) - len(prun)
        put("StrictNineCtlSentence", f"run on the {len(st9c)} with at most {num(plan['max_unknowns'])} unknowns" +
            (f" (the other {nleft} need more memory than our runs allowed)" if nleft else "") + ", it confirms " +
            ("all of them" if okc == len(st9c) else f"{okc}; for the others the recognition of this version takes another "
             "route (the version 8 decision stands)"), "recog/out_strict9_ctl: STRICT; recog/ctl_plan.json")

    # ------------------------------------------------------------------- second order (A2)
    second_strict = sorted(k for k in keys if cls[k][1] == "strict_second")
    second_desc = sorted(k for k in keys if cls[k][1] == "notlocal_descent")
    notinf = sorted(k for k, j in ps9.items() if j.get("flex_space_exact") and j.get("flexes_beyond_rotations", 0) > 0
                    and all(j.get("lifted_flexes_exact", [])) and j.get("stress_equilibrium_exact") and j.get("stress_positive_exact"))
    fk = {k: ps9[k]["flexes_beyond_rotations"] for k in notinf}
    names_ = {(11, 604): "the kissing configuration of 604 points in $\\R^{11}$ of Bianchi, Kwon, Pappu and Zou~\\cite{BKPZ26}",
              (22, 166): "the code of 166 points in $\\R^{22}$ with cosine $1/7$, Cohn, 2026"}
    if notinf:
        dl = [sci_up(ps9[k]["dist_ls_upper"]) for k in notinf]
        resp = "" if len(notinf) == 1 else ", respectively"
        put("NotInfStatement", "The exact code" + ("s " if len(notinf) > 1 else " ") +
            " and ".join(f"$({k[0]},{k[1]})$ ({names_[k]})" if k in names_ else f"$({k[0]},{k[1]})$" for k in notinf) +
            " lie within " + (dl[0] if len(set(dl)) == 1 else " and ".join(dl) + resp.replace(", r", " r") + ",") +
            " of the stored codes and are not infinitesimally jammed: their contact frameworks carry exact strictly positive "
            "equilibrium stresses and have exactly " + " and ".join(str(fk[k]) for k in notinf) + " independent infinitesimal "
            "flexes beyond the rotations" + resp + ".",
            "prestress/out (prestress9.py): flex_space_exact, flexes_beyond_rotations, dist_ls_upper")
        put("NotInfCert", "For each code the lifted vectors are exact flexes ($\\mathcal M f=0$ over the field), they and the "
            "rotations have rank $g+k$ modulo a prime, and $\\mathcal M^{\\mathsf T}\\mathcal M$ has rank $|V|d-g-k$ over "
            "$\\F_p$ (" + ", ".join(f"$p={num(ps9[k]['rank_trials'][-1]['prime'])}$" for k in notinf) + "); the stresses are "
            "checked exactly as in Section~\\ref{sec:exactres}. The lifts needed " +
            " and ".join(f"{ps9[k]['primes_used']} primes" for k in notinf) + ".", "prestress/out")
    else:
        put("NotInfStatement", "\\textbf{(no exact flexes were obtained)}", "prestress/out")
        put("NotInfCert", "--", "prestress/out")
    # the untrusted explorations (logs of explore9_*.py)
    import re as _re
    det = []
    lp = os.path.join(HERE, "prestress", "explore", "descent_d22_n166.log")
    if os.path.exists(lp):
        etas = [float(x) for x in _re.findall(r"eta (-?[0-9.e-]+) ", open(lp).read())]
        det.append(f"For $(22,166)$, the linear programs of Proposition~\\ref{{prop:descent}} (\\nolinkurl{{explore9_descent.py}}, "
                   f"untrusted) find no descent in the {len(etas)} flex directions tried (largest margin $\\eta={max(etas):.2f}$, "
                   "where a descent needs $\\eta>0$), so by duality each of these directions has a nonnegative stress with "
                   "$H_\\omega\\ge0$ on it")
    def hsign(j):
        H = j.get("H_float") or []
        if H and all(abs(x) == 0 for row in H for x in row):
            return "identically zero"
        piv = [Fraction(x.split(" + ")[0]) if "sqrt" not in x else None for x in j.get("H_pivots", [])]
        if piv and all(q is not None and q < 0 for q in piv):
            return "negative definite"
        return "positive definite" if j.get("H_pos_def") else "not positive definite"
    ex_txt = [f"$({k[0]},{k[1]})$ it is {hsign(ps9[k])}" for k in notinf if ps9[k].get("H_float") is not None]
    if ex_txt:
        det.append("for the exact strictly positive stress constant on the classes of contacts, $H_\\omega$ is computed "
                   "exactly on the exact flexes: for " + " and for ".join(ex_txt))
    cl22 = os.path.join(HERE, "prestress", "explore", "classes_d22_n166.log")
    if os.path.exists(cl22):
        eig = _re.findall(r"H eig \[([^\]]*)\]", open(cl22).read())
        ev = [[float(x) for x in e.split()] for e in eig]
        signs = "negative definite" if all(max(e) < 0 for e in ev) or all(min(e) > 0 for e in ev) else "indefinite"
        det.append("numerically (\\nolinkurl{explore9_classes.py}), every strictly positive stress constant on the classes of "
                   "contacts makes $H_\\omega$ negative definite on the flexes, so Theorem~\\ref{thm:second} cannot be applied "
                   "with a symmetric stress")
    cl11 = os.path.join(HERE, "prestress", "explore", "classes_d11_n604.log")
    sh11 = os.path.join(HERE, "prestress", "explore", "shape_d11_n604.log")
    if os.path.exists(cl11) and os.path.exists(sh11):
        eig = _re.findall(r"H eig \[([^\]]*)\]", open(cl11).read())
        hmax = max(abs(float(x)) for e in eig for x in e.split())
        txt = open(sh11).read()
        sup = _re.search(r"union support (\d+)", txt).group(1)
        fixed = _re.search(r"fixed points (\d+)", txt).group(1)
        rk = _re.search(r"rank of the fixed points (\d+)", txt).group(1)
        det.append(f"for $(11,604)$ the flexes move {sup} of the 604 points (the other {fixed} span a subspace of dimension "
                   f"{rk}), and numerically $H_\\omega$ vanishes on them, to {sci_up(hmax)}, for every stress constant on the "
                   "classes of contacts; whether the code is a local minimum, and whether it lies in a continuous family of "
                   "kissing configurations, are open")
    fam = os.path.join(HERE, "prestress", "explore", "family_d11_n604.log")
    if os.path.exists(fam):
        rows_f = _re.findall(r"h ([0-9.e-]+) residual ([0-9.e+-]+) iterations \d+ max non-contact ([0-9.e+-]+) max cos "
                             r"([0-9.e+-]+) distance to the orbit ([0-9.e+-]+)", open(fam).read())
        good = [(float(h_), float(r_), float(nc_), float(dd)) for h_, r_, nc_, mc, dd in rows_f
                if float(r_) < 1e-13 and float(nc_) < 0.5]
        if good:
            hb = max(good, key=lambda t: t[3])
            # "at least" rounds down, "at most" rounds up
            fdist = f"{math.floor(hb[3] * 100) / 100:.2f}"
            fnc = f"{math.ceil(hb[2] * 10 ** 4) / 10 ** 4:.4f}"
            put("FamilyDist", fdist, "prestress/explore/family_d11_n604.log: largest distance with every contact at 1/2 (down)")
            put("FamilyNonContact", fnc, "the same: largest non-contact inner product there (up)")
            put("FamilyResid", sci_up(max(t[1] for t in good)), "the same: largest residual of the contact equations")
            det.append(f"numerically (\\nolinkurl{{explore9_family.py}}), Gauss--Newton projections along a flex of $(11,604)$ "
                       f"reach configurations at least {fdist} away from it (least squares, modulo $O(11)$) with every "
                       f"contact at cosine $1/2$ to {sci_up(max(t[1] for t in good))} and every other inner product at most "
                       f"{fnc}, so the code appears to lie in a continuous family of kissing configurations "
                       "(Conjecture~\\ref{conj:kiss})")
    put("SecondOrderDetail", (" ".join(x[0].upper() + x[1:] + "." for x in det)) if det else "",
        "prestress/out and prestress/explore/*.log (explore9_*.py, untrusted)")
    r30 = above.get((30, 656), {})
    put("FlexThirtyCount", str(int(r30.get("kernel_dim_upper", 0)) - int(r30.get("kernel_vectors", 0))),
        "paper_v8/above/out/d30_n656.json: kernel_dim_upper - kernel_vectors")
    if notinf:
        kiss_txt = (", the first of them the 604-point kissing configuration in $\\R^{11}$ found in 2026~\\cite{BKPZ26},"
                    if notinf and notinf[0] == (11, 604) else "")
        two = " and ".join(f"$({k[0]},{k[1]})$" for k in notinf)
        put("SecondOrderFinding", "the exact codes " + two + kiss_txt + " are not infinitesimally jammed, and the "
            "second-order tests, with the stresses we tried, do not decide whether they are jammed"
            + ("; numerically the kissing configuration lies in a continuous family of kissing configurations "
               "(Conjecture~\\ref{conj:kiss})." if "FamilyDist" in M else "."), "prestress/out, prestress/explore")
        put("AbstractSecondOrder", "that the exact codes " + " and ".join(f"$({k[0]},{k[1]})$" for k in notinf) +
            " (the first a 2026 kissing configuration) are not infinitesimally jammed, ", "prestress/out")
    else:
        put("SecondOrderFinding", "the second-order tests leave the three flexing exact codes open.", "prestress/out")
        put("AbstractSecondOrder", "", "prestress/out")
    put("SecondOrderStrictClause", (f", and {len(second_strict)} by the second-order criterion" if second_strict else ""),
        "prestress/out")
    para = [f"$({k[0]},{k[1]})$: " + str(ps9[k].get("status")) + (f", {ps9[k].get('flexes_beyond_rotations')} flexes beyond "
            f"the rotations, lifts with {ps9[k].get('primes_used')} primes ({ps9[k].get('modulus_bits')} bits), "
            f"$H_\\omega$ = {ps9[k].get('H_float')}" if ps9[k].get('flex_space_exact') else "") for k in sorted(ps9)]
    put("SecondOrderSupp", ("\\nolinkurl{prestress/prestress9.py} gives: " + "; ".join(para) + ". ") if para else "--",
        "prestress/out")
    if os.path.exists(fam):
        TF2 = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{rrrrr}", "\\toprule",
               "step $h$ & residual & largest non-contact & largest cosine & distance modulo $O(11)$ \\\\", "\\midrule"]
        for h_, r_, nc_, mc, dd in rows_f:
            TF2.append(f"{h_} & {sci_up(float(r_))} & {float(nc_):.6f} & {float(mc):.6f} & {float(dd):.4f} \\\\")
        TF2 += ["\\bottomrule", "\\end{tabular}"]
        open(os.path.join(HERE, "v10s_tab_family.tex"), "w").write("\n".join(TF2) + "\n")
    if os.path.exists(lp):
        put("DescentEtaList", ", ".join(f"{x:.3f}" for x in etas), "prestress/explore/descent_d22_n166.log: eta of each direction")

    # ------------------------------------------------------------------- exact descriptions (A5)
    exact_routes = {"strict_hyper": "recognition (v8)", "strict_recog9": "recognition (v9)", "strict_above": "exact core",
                    "strict_exact": "exact arithmetic", "strict_new": "exact core", "strict_second": "exact core"}
    exk = [k for k in strict_k if cls[k][1] in exact_routes]
    put("NExactDesc", num(len(exk)), "strict entries with an exact description (recognition, exact cores)")
    deg = {}
    for k in exk:
        c = cls[k][1]
        if c == "strict_hyper":
            deg[k] = s7[k].get("field_degree")
        elif c == "strict_recog9":
            deg[k] = st9[k]["field_degree"]
        elif c in ("strict_above", "strict_second"):
            deg[k] = 1 if int(above[k].get("disc", 1)) == 1 else 2
        else:
            deg[k] = 1
    put("FieldDegMax", str(max(v for v in deg.values() if v)), "max degree of the field over the exact descriptions")
    nf = {sub: len(glob.glob(os.path.join(HERE, "exactdata", sub, "*.json.gz"))) for sub in ("open", "v8recog", "core")}
    have = set()
    for sub in ("open", "v8recog", "core"):
        for f_ in glob.glob(os.path.join(HERE, "exactdata", sub, "*.json.gz")):
            m_ = _re.match(r"d(\d+)_n(\d+)\.json\.gz$", os.path.basename(f_))
            have.add((int(m_.group(1)), int(m_.group(2))))
    nofile = sorted(k for k in exk if k not in have)
    plan_ = json.load(open(os.path.join(HERE, "recog", "ctl_plan.json")))
    notrun_ = {(e["d"], e["n"]) for e in plan_["entries"] if not e["run"]}
    assert DRAFT or all(cls[k][1] == "strict_hyper" and k in notrun_ for k in nofile), "an exact description without a file"
    put("NExactFiles", num(len(exk) - len(nofile)), "exact descriptions with a file in exactdata/")
    put("NExactNoFile", num(len(nofile)), "exact descriptions without a file (version 8 recognitions not rerun, memory)")
    put("ExactDataSentence", f"{nf['open']} files for the codes recognised at higher degree, {nf['v8recog']} for codes "
        f"recognised at lower degree (written when \\nolinkurl{{strict9.py}} decided them again) and {nf['core']} for the "
        f"exact cores", "exactdata/{open,v8recog,core}/*.json.gz")
    dc = collections.Counter(deg.values())
    put("FieldSentence", f"{dc[1]} have rational inner products, {dc[2]} lie in quadratic fields and "
        f"{sum(v for d_, v in dc.items() if d_ and d_ > 2)} in fields of degree 3 to {max(v for v in deg.values() if v)}",
        "the same")
    TF = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{l" + "r" * 6 + "}", "\\toprule",
          "route & $\\Q$ & degree 2 & 3--4 & 5--8 & 9--" + str(max(v for v in deg.values() if v)) + " & total \\\\", "\\midrule"]
    bins = [(1, 1), (2, 2), (3, 4), (5, 8), (9, 99)]
    for rt in sorted(set(exact_routes.values())):
        ks = [k for k in exk if exact_routes[cls[k][1]] == rt]
        if not ks:
            continue
        TF.append(rt + " & " + " & ".join(str(sum(1 for k in ks if lo <= (deg[k] or 0) <= hi)) for lo, hi in bins) +
                  f" & {len(ks)} \\\\")
    TF.append("\\midrule all & " + " & ".join(str(sum(1 for k in exk if lo <= (deg[k] or 0) <= hi)) for lo, hi in bins) +
              f" & {len(exk)} \\\\")
    TF += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_fields.tex"), "w").write("\n".join(TF) + "\n")

    # ------------------------------------------------------------------- FLINT ranks (B)
    agree = [k for k, j in fr9.items() if j.get("status") == "AGREE"]
    dis = [k for k, j in fr9.items() if j.get("status") not in ("AGREE",) and not str(j.get("status", "")).startswith("skipped")]
    # every strict entry above the size limit without a FLINT rank: not recomputed because of memory (dense matrix size)
    sa = [k for k in keys if cls[k][1] == "strict_above"]
    skip = sorted(k for k in sa if k not in set(agree) and k not in set(dis))
    gb = lambda k: int(above[k]["cols"]) ** 2 * 8 / 1e9     # the dense FLINT matrix (8 bytes per entry), in GB
    if skip and agree:
        # the run went through the entries in the order of their size; the ones left are exactly the largest
        assert DRAFT or min(gb(k) for k in skip) > max(gb(k) for k in agree), "a skipped FLINT rank is not among the largest"
    put("NFlintRankAgree", num(len(agree)), "ranks/out (flintrank9.py): AGREE")
    put("NFlintRankSkip", num(len(skip)), "strict entries above the size limit without a FLINT rank")
    if skip:
        put("FlintRankSkipGB", f"{min(gb(k) for k in skip):.2f} to {max(gb(k) for k in skip):.2f}",
            "dense FLINT matrix of the entries left (cols^2 * 8 bytes)")
        put("FlintRankSkipNeedGB", f"{2 * min(gb(k) for k in skip):.1f}", "twice the smallest dense matrix left")
    if agree:
        put("FlintRankDoneGB", f"{max(gb(k) for k in agree):.2f}", "largest dense FLINT matrix computed")
    put("FlintRankSkipList", (", ".join(f"$({k[0]},{k[1]})$ ({gb(k):.2f}~GB)" for k in sorted(skip, key=gb)))
        if skip else "none", "strict entries above the size limit without a FLINT rank (dense size)")
    primes = sorted(set(j.get("prime") for j in fr9.values() if j.get("prime")))
    put("FlintRankPrime", num(primes[0]) if len(primes) == 1 else ", ".join(num(p) for p in primes[:-1]) + " and " +
        num(primes[-1]), "ranks/out: primes")
    put("FlintRankPrimeClause", (f"$p={num(primes[0])}$" if len(primes) == 1 else
         f"a prime $p<2^{{26}}$ in which $D$ is a square (one of {', '.join(num(q) for q in primes[:-1])} and "
         f"{num(primes[-1])})"), "ranks/out: primes (kfield.prime_for below 2^26)")
    put("FlintRankSentence", f"the rank agrees with run\\_above.py's and equals $|V|d-g$ for all {len(agree)} entries "
        "computed, whose dense matrices take up to " + (f"{max(gb(k) for k in agree):.2f}" if agree else "0") + "~GB" +
        (f"; the {len(skip)} largest were not recomputed, because FLINT's rank routine works on a copy of the dense matrix, "
         f"so that each of them needs more than {2 * min(gb(k) for k in skip):.1f}~GB of memory, which our runs did not "
         "allow: " + ", ".join(f"$({k[0]},{k[1]})$ ({gb(k):.2f}~GB)" for k in sorted(skip, key=gb)) if skip else "") +
        (f"; {len(dis)} disagree" if dis else ""), "ranks/out")
    # PARI
    pp = os.path.join(HERE, "ranks", "pari_check9.json")
    if os.path.exists(pp):
        pj = json.load(open(pp))
        put("PariSentence", ", and the exact Gram matrices of " + pj["summary_short"] + " were checked again with PARI/GP",
            "ranks/pari_check9.json")
        put("PariSupp", pj["summary_long"], "ranks/pari_check9.json")

    # ------------------------------------------------------------------- the intrinsic cut (A3)
    nk = list(neg)
    done = [k for k in nk if neg[k]["cut_done"]]
    put("CutDone", num(len(done)), "cut9/out: entries with the negative certificate searched")
    samp9 = [k for k in done if cut[k].get("group") == "sample" and neg[k]["res8"]]
    g = [neg[k]["d9"] / neg[k]["d8"] for k in samp9]
    put("CutSampleN", num(len(samp9)), "cut9/out: the random sample (group 'sample') with a resolved verdict of version 8")
    put("CutGainMed", f"{med(g):.1f}" if g else "--", "cut9/out, random sample: median of delta9/delta8")
    put("CutDeltaMed", sci(med([neg[k]["d9"] for k in samp9])) if samp9 else "--",
        "cut9/out, random sample: median version 9 radius")
    put("CutDeltaMedGrid", sci(med([neg[k]["d8"] for k in samp9])) if samp9 else "--",
        "cut9/out, random sample: median version 8 radius")
    unfl = [k for k in nk if not neg[k]["flagged"]]
    unfl_res = [k for k in unfl if neg[k]["res9"]]
    put("NNegUnflagAll", num(len(unfl)), "negative-certificate entries the table does not flag")
    put("NNegUnflag", num(len(unfl_res)), "the same, resolved in version 9")
    put("CutDeltaMedUnflag", sci(med([neg[k]["d9"] for k in unfl_res])) if unfl_res else "--",
        "median version 9 radius over them")
    put("CutNewlyResolved", num(sum(1 for k in nk if neg[k]["res9"] and not neg[k]["res8"])),
        "entries unresolved in version 8 and resolved at the cut")
    rest_not = len([k for k in nk if not neg[k]["cut_done"]])
    put("CutRestClause", "" if rest_not == 0 else f"; the other {num(rest_not)} keep their radius on the tolerance grid",
        "entries not searched")
    pop = [k for k in nk if neg[k]["flagged"] and neg[k]["res8"] and cut.get(k, {}).get("group") in ("sample", None)
           and not (k in cut and cut[k].get("group") in ("unflagged", "unresolved"))]
    put("CutSamplePop", num(len(pop)), "the population of the random sample: flagged entries resolved on the grid")
    # the size of the window E of the certificate that sets the radius (review of v9, M4)
    win = {}
    for k in nk:
        if not neg[k]["res9"]:
            continue
        if neg[k]["from_cut"]:
            win[k] = int(cut[k]["exact"]["contacts"])
        else:
            win[k] = int(fl8[k]["contacts"])
    wv = sorted(win.values())
    put("NegWinN", num(len(wv)), "resolved negative verdicts with the window of the certificate that sets the radius")
    put("NegWinMed", num(int(med(wv))), "median number of pairs in that window")
    put("NegWinMax", num(wv[-1]), "largest window")
    put("NegWinOne", num(sum(1 for w in wv if w == 1)), "radius set by a certificate on a single pair")
    put("NegWinLeThree", num(sum(1 for w in wv if w <= 3)), "radius set by a certificate on at most 3 pairs")
    put("NegWinLeTen", num(sum(1 for w in wv if w <= 10)), "radius set by a certificate on at most 10 pairs")
    gw = [int(fl8[k]["contacts"]) for k in nk if neg[k]["res8"] and k in fl8]
    put("NegGridWinOne", num(sum(1 for w in gw if w == 1)), "grid certificates (resolved) on a single pair")
    put("NegGridWinLeTen", num(sum(1 for w in gw if w <= 10)), "grid certificates (resolved) on at most 10 pairs")
    put("NegGridWinN", num(len(gw)), "resolved grid certificates")
    ex1 = sorted(k for k, w in win.items() if w == 1)
    put("NegWinOneEx", elist(ex1[:3]) if ex1 else "none", "examples of single-pair certificates")

    # ------------------------------------------------------------------- the taxonomy (Table 1)
    def rng(ks, key="dist"):
        v = [info[k][key] for k in ks if info[k].get(key) is not None]
        return (sci_up(min(v)) + " to " + sci_up(max(v))) if v else "--"
    loc_k = [k for k in keys if cls[k][0] == "local"]
    ns_k = [k for k in keys if cls[k][0] == "nonstrict"]
    nl_k = [k for k in keys if cls[k][0] == "notlocal"]
    fx_k = [k for k in keys if main_class(k) == "flex"]
    nv_k = [k for k in keys if main_class(k) == "none"]
    tl = collections.Counter(info[k]["tier"] for k in loc_k)
    tn = collections.Counter(info[k]["tier"] for k in ns_k)
    d9 = [neg[k]["d9"] for k in fx_k]
    TT = ["% written by make_numbers_v10.py; do not edit",
          "\\begin{tabular}{p{0.16\\linewidth}p{0.35\\linewidth}rrrrp{0.15\\linewidth}}", "\\toprule",
          "verdict & what is certified near the stored code & entries & A & B & C & distance \\\\", "\\midrule",
          f"strict local minimum & an exact strict local minimum modulo rotations and rattlers & {num(len(strict_k))} & "
          f"{num(tiers['A'])} & {num(tiers['B'])} & {num(tiers['C'])} & {rng(strict_k)} \\\\",
          f"local minimum, strictness open & an exact local minimum in an explicit tube (hyperstatic) & {num(len(loc_k))} & "
          f"{num(tl['A'])} & {num(tl['B'])} & {num(tl['C'])} & {rng(loc_k)} \\\\",
          f"local minimum, not strict & an exact local minimum, neither strict nor jammed & {num(len(ns_k))} & "
          f"{num(tn['A'])} & {num(tn['B'])} & {num(tn['C'])} & {rng(ns_k)} \\\\",
          f"not a local minimum & an exact configuration that is not a local minimum (flat core"
          + (", second-order descent" if second_desc else "") + f") & {num(len(nl_k))} & & & & {rng(nl_k)} \\\\",
          f"no local minimum within $\\delta$ & no local minimum within $\\delta$ (largest move of a point, modulo rotations), "
          f"$\\delta$ set by the separation, $\\delta\\ge10^{{-12}}$ & {num(len(fx_k))} & & & & $\\delta$ {sci_dn(min(d9))} to "
          f"{sci_dn(max(d9))} \\\\",
          f"no verdict & every attempt failed or was above the limits; or the negative certificate is unresolved & "
          f"{num(len(nv_k))} & & & & -- \\\\", "\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_taxonomy.tex"), "w").write("\n".join(TT) + "\n")
    put("NotLocalDistMax", rng(nl_k), "distances of the not-local-minimum verdicts")

    # ------------------------------------------------------------------- census table by dimension
    cols = ["strict", "local", "nonstrict", "notlocal", "flex", "none"]
    names = ["strict", "open", "not strict", "not loc.\\ min.", "negative", "none"]
    TC = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{l" + "r" * (len(cols) + 2) + "}", "\\toprule",
          "$d$ & " + " & ".join(names) + " & all & unmarked \\\\", "\\midrule"]
    for i, (a, b) in enumerate(RANGES):
        ks = [k for k in keys if a <= k[0] <= b]
        cc = collections.Counter(main_class(k) for k in ks)
        TC.append(f"{a}--{b} & " + " & ".join(num(cc[c]) for c in cols) + f" & {num(len(ks))} & "
                  f"{num(sum(1 for k in ks if not idx[k]['flagged']))} \\\\")
    TC.append("\\midrule all & " + " & ".join(num(C[c]) for c in cols) + f" & {num(len(keys))} & "
              f"{num(sum(1 for k in keys if not idx[k]['flagged']))} \\\\")
    TC += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_census.tex"), "w").write("\n".join(TC) + "\n")

    # ------------------------------------------------------------------- structure (A4)
    pos = [k for k in keys if cls[k][0] in ("strict", "local", "nonstrict")]
    def iso(k):
        return cls[k][1] == "strict_iso"
    def ratt(k):
        try:
            return int(R8[k]["nonrattlers"]) < k[1]
        except ValueError:
            return False
    def cpp(k):
        try:
            return 2 * int(R8[k]["contacts"]) / int(R8[k]["nonrattlers"])
        except (ValueError, ZeroDivisionError):
            return None
    TS = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{lrrrrrrrr}", "\\toprule",
          "$d$ & entries & cert.\\ loc.\\ min. & isostatic & rattlers & contacts/point & $N/d$ & $N/d$, neg. & exact (strict) \\\\",
          "\\midrule"]
    for a, b in RANGES:
        ks = [k for k in keys if a <= k[0] <= b]
        ps_ = [k for k in pos if a <= k[0] <= b]
        ng = [k for k in fx_k if a <= k[0] <= b]
        stk = [k for k in strict_k if a <= k[0] <= b]
        cp = [cpp(k) for k in ps_ if cpp(k) is not None]
        TS.append(f"{a}--{b} & {num(len(ks))} & {num(len(ps_))} & {pct(sum(1 for k in ps_ if iso(k)), len(ps_))}\\% & "
                  f"{pct(sum(1 for k in ps_ if ratt(k)), len(ps_))}\\% & {med(cp):.1f} & {med([k[1] / k[0] for k in ps_]):.1f} & "
                  f"{med([k[1] / k[0] for k in ng]):.1f} & {pct(sum(1 for k in stk if k in set(exk)), len(stk))}\\% \\\\")
    TS += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_structure.tex"), "w").write("\n".join(TS) + "\n")
    lowp = [k for k in pos if k[0] <= 5]
    highp = [k for k in pos if k[0] >= 25]
    put("StructIsoLowPct", pct(sum(1 for k in lowp if iso(k)), len(lowp)), "certified local minima with d<=5: isostatic share")
    put("StructIsoHighPct", pct(sum(1 for k in highp if iso(k)), len(highp)), "the same, d>=25")
    put("StructRatPct", pct(sum(1 for k in pos if ratt(k)), len(pos)), "certified local minima with a rattler")
    for name, lo, hi in (("Small", 0, 4), ("Mid", 4, 16), ("Large", 16, 10 ** 9)):
        ks = [k for k in keys if lo * k[0] < k[1] <= hi * k[0]] if lo else [k for k in keys if k[1] <= hi * k[0]]
        put("StructShare" + name, pct(sum(1 for k in ks if k in set(pos)), len(ks)), f"share of certified local minima, {lo}d<N<={hi}d")
    isohigh = [k for k in keys if k[0] >= 25 and iso(k)]
    fam = [(d, 2 * d + 1) for d in range(3, 33) if (d, 2 * d + 1) in idx]
    fam_ok = [k for k in fam if iso(k) and not ratt(k)]
    put("StructFamilySentence", f"All {len(isohigh)} isostatic certificates with $d\\ge25$ have $2d<N\\le"
        f"{max(k[1] / k[0] for k in isohigh):.2f}\\,d$. The codes with $N=2d+1$ form a uniform family: for every $d$ from "
        f"{min(k[0] for k in fam)} to {max(k[0] for k in fam)} the table's code of $2d+1$ points is an isostatic strict local "
        f"minimum without rattlers ({len(fam_ok)} of {len(fam)}), with exactly $\\iso=(d-1)(3d+2)/2+1$ contacts.",
        "v9 classes and v8_verdicts.tsv contacts/nonrattlers")
    put("ConjFamilyText", "For every $d\\ge3$, the best known code of $2d+1$ points in $\\R^d$ is an isostatic strict local "
        "minimum of the largest cosine without rattlers: its contact graph has exactly $(d-1)(3d+2)/2+1$ edges, carries a "
        f"unique strictly positive stress and is bar-rigid. (Verified for the table's codes for $3\\le d\\le32$, "
        f"{len(fam_ok)} of {len(fam)}.)", "the same")
    for k in fam:
        assert int(R8[k]["contacts"]) == (k[0] - 1) * (3 * k[0] + 2) // 2 + 1 or not iso(k), k

    # ------------------------------------------------------------------- the live table, with and without the author
    live = {}
    for line in open(os.path.join(CERT, "v5", "live", "live_20261107.tsv")):
        if line.startswith("#"):
            continue
        p_ = line.rstrip("\n").split("\t")
        live[(int(p_[0]), int(p_[1]))] = (p_[2], p_[4] if len(p_) > 4 else "")
    from decimal import Decimal
    def lower(k):
        return k in live and Decimal(live[k][0]) < Decimal(idx[k]["cos"])
    def author(k):
        return k in live and "Arhan Canli" in live[k][1]
    listed = [k for k in keys if k in live]
    put("LiveLowerAll", num(sum(1 for k in listed if lower(k))), "v5/live/live_20261107.tsv: listed census entries with a lower displayed cosine")
    auth = [k for k in listed if author(k)]
    put("LiveAuthorAll", num(len(auth)), "the same: credited to the author")
    byc = collections.Counter({"strict": "strict", "local": "open", "nonstrict": "not strict", "notlocal": "not local",
                               "flex": "with the negative verdict", "none": "without a verdict"}[main_class(k)] for k in auth)
    put("LiveAuthorByClass", ", ".join(f"{v} {nm}" for nm, v in sorted(byc.items(), key=lambda t: -t[1])) +
        f"; {sum(1 for k in auth if lower(k))} of them display a lower cosine now", "credited census entries by v9 class")
    TL = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{lrrrrrr}", "\\toprule",
          " & \\multicolumn{3}{c}{all unmarked, listed} & \\multicolumn{3}{c}{without the author's} \\\\",
          "verdict & listed & lower now & share & listed & lower now & share \\\\", "\\midrule"]
    stats = {}
    for c, nm in zip(cols, names):
        ks = [k for k in listed if main_class(k) == c and not idx[k]["flagged"]]
        lw = [k for k in ks if lower(k)]
        ks2 = [k for k in ks if not author(k)]
        lw2 = [k for k in ks2 if lower(k)]
        stats[c] = (len(ks), len(lw), len(ks2), len(lw2))
        TL.append(f"{nm} & {num(len(ks))} & {num(len(lw))} & {pct(len(lw), len(ks), 1)}\\% & {num(len(ks2))} & "
                  f"{num(len(lw2))} & {pct(len(lw2), len(ks2), 1)}\\% \\\\")
    TL += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_live.tex"), "w").write("\n".join(TL) + "\n")
    a, b_, a2, b2 = stats["flex"]
    put("LiveUnmarkedFlexN", num(a), "unmarked listed entries with the negative verdict (v9)")
    put("LiveUnmarkedFlexImp", num(b_), "the same, lower displayed cosine now")
    put("LiveUnmarkedFlexPct", pct(b_, a, 1), "the share")
    put("LiveUnmarkedFlexNNA", num(a2), "the same without the author-credited entries")
    put("LiveUnmarkedFlexImpNA", num(b2), "the same, lower now")
    put("LiveUnmarkedFlexPctNA", pct(b2, a2, 1), "the share")
    a, b_, a2, b2 = stats["strict"]
    put("LiveUnmarkedStrictPct", pct(b_, a, 1), "unmarked listed strict entries with a lower displayed cosine")
    put("LiveUnmarkedStrictPctNA", pct(b2, a2, 1), "the same without the author-credited entries")
    gains = [float(Decimal(idx[k]["cos"]) - Decimal(live[k][0])) for k in listed if main_class(k) == "flex" and lower(k)]
    put("LiveFlexGainMed", sci(med(gains)), "median (archived displayed - live displayed) over improved negative-verdict entries")
    lv = load_dir(os.path.join(P8, "live", "out"))
    fsamp = [j for j in lv.values() if j.get("census_class") == "flex" and j.get("dist_lower_point") is not None]
    put("LiveFlexNewConf", num(sum(1 for j in fsamp if j["dist_lower_point"] >= 1e-4)),
        "paper_v8/live/out: sampled negative entries with dist_lower_point >= 1e-4")
    put("LiveFlexPolish", num(sum(1 for j in fsamp if j["dist_lower_point"] < 1e-6)),
        "the same with dist_lower_point < 1e-6")
    put("LiveFlexSampleN", num(len(fsamp)), "paper_v8/live/out: sampled negative entries")
    # minor 16: the beaten certified minima of version 8 that were not fetched (the definition of numbers_v8_extra.py)
    def cstar_lo8(k):
        v8c = R8[k]["verdict"]
        st = "PASS" if R8[k]["certificate"] == "strict_iso" else "LOCALMIN"
        e_ = ex.get((k[0], k[1], st))
        if e_ and R8[k]["certificate"] in ("strict_iso", "strict_hyper", "strict_sym", "local"):
            lo_ = Fraction(int(e_["c_interval_30"][0]), 10 ** 30)
            if v8c == "local":
                lo_ -= Fraction(e_["chi"])
            return lo_
        return Fraction(idx[k]["cos"]) - Fraction(1, 10 ** 12)
    cert8 = [k for k in keys if R8[k]["verdict"] in ("strict", "local", "nonstrict") and k in live]
    beat8 = [k for k in cert8 if Decimal(live[k][0]) < Decimal(idx[k]["cos"]) and Fraction(live[k][0]) < cstar_lo8(k)]
    notf = sorted(k for k in beat8 if k not in lv)
    put("LiveBeatenNotFetched", (elist(notf) + (", whose certificate was obtained after the codes to fetch had been selected"
        if all(R8[k]["certificate"] in ("strict_above", "strict_exact", "strict_new") for k in notf) else "")) if notf else "none",
        "beaten version 8 certified minima (numbers_v8_extra.py definition) not in paper_v8/live/out")
    # the author's uploads
    man = os.path.join(os.path.dirname(CERT), "codes", "MANIFEST.tsv")
    if os.path.exists(man):
        ups = [l.rstrip("\n").split("\t") for l in open(man)][1:]
        dates = sorted(u[4][:10] for u in ups if u and len(u) > 4)
        mon = {"10": "October", "11": "November"}
        fmt = lambda x: f"{int(x[8:10])} {mon[x[5:7]]} {x[:4]}"
        put("AuthorUploads", num(len(ups)), "codes/MANIFEST.tsv (the author's uploads to the live table)")
        put("AuthorUploadFirst", fmt(dates[0]), "codes/MANIFEST.tsv: first upload date")
        put("AuthorUploadLast", fmt(dates[-1]), "codes/MANIFEST.tsv: last upload date")
    # the kissing entries
    kiss = [(11, 592), (11, 593), (12, 841)]
    put("KissNegSentence", ", ".join(f"{sci_dn(neg[k]['d9'])}" for k in kiss if k in neg) + " respectively",
        "version 9 radii of (11,592), (11,593), (12,841)")

    # ------------------------------------------------------------------- the negative certificate table
    TN = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{lrrrrrr}", "\\toprule",
          "entries & number & searched & median $\\delta$, grid & median $\\delta$, cut & median gain & resolved at the cut \\\\",
          "\\midrule"]
    def qtxt(v):
        return sci(med(v)) if v else "--"
    groups = [("all", nk), ("unmarked", [k for k in nk if not neg[k]["flagged"]]),
              ("marked", [k for k in nk if neg[k]["flagged"]]), ("unresolved on the grid", [k for k in nk if not neg[k]["res8"]])]
    for nm, ks in groups:
        ds = [k for k in ks if neg[k]["cut_done"]]
        TN.append(f"{nm} & {num(len(ks))} & {num(len(ds))} & {qtxt([neg[k]['d8'] for k in ds])} & "
                  f"{qtxt([neg[k]['d9'] for k in ds])} & {med([neg[k]['d9'] / neg[k]['d8'] for k in ds]) if ds else 0:.1f} & "
                  f"{num(sum(1 for k in ds if neg[k]['res9']))} \\\\")
    TN += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_neg.tex"), "w").write("\n".join(TN) + "\n")
    kv8 = [cut[k]["v8"]["contacts"] for k in done if cut[k].get("v8", {}).get("contacts")]
    kc = [cut[k]["k"] for k in done if cut[k].get("k")]
    put("CutKVEightMed", num(int(med(kv8))) if kv8 else "--", "cut9/out: median |E| of version 8")
    put("CutKMed", num(int(med(kc))) if kc else "--", "cut9/out: median |E| at the cut (where decided)")
    rv8 = [cut[k]["v8"]["rows"] for k in done if cut[k].get("v8", {}).get("rows")]
    rc = [cut[k]["exact"]["rows"] for k in done if cut[k].get("exact")]
    put("CutRowsVEightMed", num(int(med(rv8))) if rv8 else "--", "cut9/out: median rows of version 8")
    put("CutRowsMed", num(int(med(rc))) if rc else "--", "cut9/out: median rows at the cut")
    hb = sum(1 for k in done if neg[k]["bycut"] == "h" and neg[k]["dcut"] is not None)
    put("CutHBindSentence", f"At the cut, $\\delta$ is set by $h$ for {num(hb)} entries and by the separation for the others.",
        "cut9/out: delta_by of the exact decisions")
    st_c = collections.Counter(cut[k].get("status") for k in done)
    plain_st = {"FLEX": "decided at the cut", "no larger radius than version 8": "no larger radius than on the grid",
                "no admissible cut": "no admissible cut"}
    put("CutFailSentence", ", ".join(f"{num(v)} {plain_st.get(s, s)}" for s, v in sorted(st_c.items(), key=lambda t: -t[1])),
        "cut9/out: statuses (plain names)")
    gr = collections.Counter(cut[k].get("group") for k in done)
    put("CutGroupsSentence", ", ".join(f"{num(v)} {s}" for s, v in sorted(gr.items(), key=lambda t: -t[1])), "cut9/out: groups")

    # ------------------------------------------------------------------- certificate data (B)
    cdv = os.path.join(HERE, "certdata", "verify.json")
    cde = os.path.join(HERE, "certdata", "export_log.json")
    if os.path.exists(cdv) and os.path.exists(cde):
        V9 = json.load(open(cdv))["files"]
        E9 = json.load(open(cde))["files"]
        posf = {k: v for k, v in V9.items() if v["kind"] == "positive"}
        negf = {k: v for k, v in V9.items() if v["kind"] == "negative"}
        put("CertDataNPos", num(len(posf)), "certdata/verify.json: positive certificate files")
        put("CertDataNNeg", num(len(negf)), "certdata/verify.json: negative certificate files")
        dims = [v["d"] for v in V9.values()]
        put("CertDataDimRange", f"{min(dims)} to {max(dims)}", "certdata/verify.json: dimensions")
        put("CertDataMaxU", num(max(e.get("U", 0) for e in E9.values())), "certdata/export_log.json: largest U")
        put("CertDataMB", f"{sum(e['bytes'] for e in E9.values()) / 2 ** 20:.0f}", "certdata/export_log.json: total size")
        same = sum(1 for v in V9.values() if v["all_same"])
        put("CertDataSentence", f"the verifier reproduces the decision of record in every field compared for "
            f"{same} of the {len(V9)} files", "certdata/verify.json: all_same")
        negh = [v for v in negf.values() if v["result"].get("delta_by") == "h"]
        put("CertDataNegSentence", f"including {len(negh)} limited by $h$", "certdata/verify.json")
    us = [int(ex[(k[0], k[1], st)]["unknowns"]) for k in keys for st in ("PASS", "LOCALMIN", "NONOPT")
          if (k[0], k[1], st) in ex]
    put("CertDataFullGB", f"{sum(u * u * 8 for u in us) / 2 ** 30:.0f}", "sum of 8 U^2 bytes over the exact decisions")

    # ------------------------------------------------------------------- CPU of the version 9 runs
    cpu = 0.0
    for sub in ("recog/out_all", "recog/out_strict9", "recog/out_strict9_ctl", "prestress/out", "cut9/out"):
        for j in maybe_dir(sub).values():
            cpu += float(j.get("cpu_seconds") or 0)
    for j in fr9.values():
        cpu += float(j.get("wall") or 0)
    put("VNineCPUHours", f"{cpu / 3600:.1f}", "sum of cpu_seconds (wall for ranks/out) over the version 9 outputs")

    # ------------------------------------------------------------------- above the size limit (supplement table)
    def above_outcome(k, j):
        c = cls[k][1]
        if c == "strict_above":
            return "strict, " + ("$\\Q$" if int(j.get("disc", 1)) == 1 else f"$\\Q(\\sqrt{{{j.get('disc')}}})$")
        if c == "strict_second":
            return "strict (second order)"
        if c == "notlocal_descent":
            return "not a local minimum (descent)"
        if c == "notlocal_flat_above":
            notes = idx[k]["notes"].lower()
            return "flat core" + ("$^u$" if ("orthogonal" in notes or "union" in notes) else "")
        return "undecided: " + str(j.get("status", "")).replace("undecided: ", "").split("(")[0].strip()
    groups_a = collections.defaultdict(list)
    for k, j in sorted(above.items()):
        groups_a[above_outcome(k, j)].append(k)
    TA = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{p{0.24\\linewidth}rp{0.66\\linewidth}}",
          "\\toprule", "outcome & entries & entries \\\\", "\\midrule"]
    for o in sorted(groups_a, key=lambda o: (not o.startswith("strict"), o)):
        TA.append(f"{o} & {len(groups_a[o])} & {elist(groups_a[o])} \\\\")
    TA += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_above.tex"), "w").write("\n".join(TA) + "\n")

    # ------------------------------------------------------------------- notable entries (supplement; minor 20)
    optref = {(3, 7): "\\cite{SvdW51}", (3, 8): "\\cite{SvdW51}", (3, 9): "\\cite{SvdW51}", (3, 10): "\\cite{Dan86}",
              (3, 12): "\\cite{FT43}", (3, 13): "\\cite{MT12}", (3, 14): "\\cite{MT15}", (3, 24): "\\cite{Rob61}",
              (4, 10): "\\cite{BV09}", (4, 24): "\\cite{dLLM24}", (4, 120): "\\cite{And99,BD01}", (5, 16): "\\cite{Lev79}",
              (6, 27): "\\cite{Lev79}", (7, 56): "\\cite{Lev79,BS81}", (8, 240): "\\cite{Lev79,OS79,BS81}",
              (16, 288): "\\cite{CdLL24}", (20, 56): "\\cite{CdLL24}", (21, 50): "\\cite{CdLL24}", (21, 77): "\\cite{CdLL24}",
              (21, 112): "\\cite{DGS77,Lev79}", (21, 162): "\\cite{DGS77,Lev79}", (22, 100): "\\cite{DGS77,Lev79}",
              (22, 275): "\\cite{DGS77,Lev79}", (22, 891): "\\cite{DGS77,Lev79,CK07}", (23, 552): "\\cite{DGS77,Lev79}"}
    for k in opt:
        assert k in optref, ("a provably optimal entry without a reference", k)
    notable = {}
    for k in opt:
        notable[k] = "optimal " + optref[k]
    for k in cjkt:
        notable[k] = (notable.get(k, "") + "; " if k in notable else "") + "CJKT"
    for k in [(10, 510), (8, 72), (9, 330), (11, 592), (11, 593), (11, 604), (12, 841), (22, 116), (22, 166), (30, 656),
              (15, 63), (3, 576), (3, 791), (17, 317)] + second_strict + second_desc:
        notable.setdefault(k, "")
    vname = {"strict": "strict", "local": "local, open", "nonstrict": "not strict", "notlocal": "not a local min.",
             "flex": "no loc.\\ min.\\ within $\\delta$", "none": "no verdict"}
    TNb = ["% written by make_numbers_v10.py; do not edit", "\\begin{tabular}{lp{0.2\\linewidth}lp{0.38\\linewidth}}",
           "\\toprule", "entry & verdict & certificate & known \\\\", "\\midrule"]
    for k in sorted(notable):
        cert_tex = cls[k][1].replace("_", "\\_")
        TNb.append(f"$({k[0]},{k[1]})$ & {vname[main_class(k)]} & {cert_tex} & {notable[k]} \\\\")
    TNb += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v10_tab_notable.tex"), "w").write("\n".join(TNb) + "\n")

    # ------------------------------------------------------------------- the figure of the verdicts (minor 24)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.4), sharex=True, sharey=True, constrained_layout=True)
    nb = np.logspace(np.log10(4), np.log10(1024), 25)
    db = np.arange(2.5, 33.5, 1.0)
    for ax, c, nm in zip(axes.ravel(), cols, ["strict", "local minimum, open", "not strict", "not a local minimum",
                                                "no local minimum within $\\delta$", "no verdict"]):
        ks = [k for k in keys if main_class(k) == c]
        H_, _, _ = np.histogram2d([k[0] for k in ks], [k[1] for k in ks], bins=[db, nb])
        Hm = np.ma.masked_where(H_.T == 0, H_.T)
        im = ax.pcolormesh(db, nb, Hm, cmap="viridis", shading="flat", vmin=1)
        ax.set_yscale("log")
        ax.set_title(f"{nm} ({len(ks):,})", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.8)
    for ax in axes[-1]:
        ax.set_xlabel("$d$")
    for ax in axes[:, 0]:
        ax.set_ylabel("$N$")
    fig.savefig(os.path.join(HERE, "fig_verdicts_v10.pdf"), metadata={"CreationDate": None})
    plt.close(fig)
    shutil.copy(os.path.join(P8, "fig_cost_v8.pdf"), os.path.join(HERE, "fig_cost_v10s.pdf"))
    shutil.copy(os.path.join(P8, "v8_tab_quant.tex"), os.path.join(HERE, "v10s_tab_quant.tex"))

    # ------------------------------------------------------------------- the per-entry file (minor 23, major 4)
    head9 = ["d", "n", "flagged", "optimal", "verdict", "certificate", "tier", "tol", "contacts", "nonrattlers", "c_lo_30",
             "c_hi_30", "dist_ls", "rho", "kappa_half", "code_inside_rho", "gap_lower", "gap_first_pass", "close_call",
             "double_cert", "double_cert_dist", "tau_prime", "h", "delta_v8", "delta_cut", "delta", "resolved",
             "field_degree", "live_cos", "live_lower", "live_author", "note"]
    with open(os.path.join(HERE, "v10_verdicts.tsv"), "w") as f:
        f.write("# Census version 9: one line per entry of the archived table with 3<=d<=32, N<=1024 (make_numbers_v10.py).\n")
        f.write("# verdict: strict | local (strictness open) | nonstrict | notlocal | flex (no local minimum within delta, "
                "resolved) | unresolved (delta < 1e-12 or set by h) | equalcontact_notlocal | none\n")
        f.write("# certificate: strict_iso, strict_hyper (recognition, version 8), strict_recog9 (recognition at higher degree, "
                "version 9), strict_sym, strict_exact, strict_new, strict_above, strict_second (second-order criterion), "
                "local, nonstrict_split, notlocal_flat, notlocal_flat_above, notlocal_descent, notlocal_nonopt, flex "
                "(version 8 radius), flex_cut9 (radius of the intrinsic cut), flex_unres, none_size, failed\n")
        f.write("# tier (positive verdicts): A = strict, dist_ls <= 1e-12 and code_inside_rho, and for a hyperstatic "
                "certificate the exact configuration shown to be the certificate's zero (Remark rem:same); B = dist_ls <= "
                "1e-12 otherwise (including every entry of open strictness); C = dist_ls > 1e-12\n")
        f.write("# c_lo_30, c_hi_30: integers; the certified cosine (of the zero of the square system) lies in "
                "[c_lo_30, c_hi_30] * 10^-30. dist_ls: certified least-squares distance to the stored code, modulo "
                "rotations, on the non-rattlers (local minima: plus the tube). rho, kappa_half: the radius and growth rate "
                "kappa/(2 sqrt d) of Theorem 3.9 (isostatic) or of the hyperstatic certificate; code_inside_rho: the stored "
                "code lies within rho (exact decision).\n")
        f.write("# close_call: the first pass's gap to the nearest non-contact pair is below 1e-9; double_cert: the larger "
                "contact graph is certified as well (paper_v8/recount/tolscan_close.jsonl); double_cert_dist: an upper bound "
                "for the distance between the two certified configurations (sum of their distances to the stored code).\n")
        f.write("# tau_prime, h: version 8 negative certificate; delta_v8, delta_cut: radius on the grid and at the intrinsic "
                "cut (exact; empty if not searched or not decided); delta: the version 9 radius; resolved: 1 if delta >= 1e-12 "
                "and set by the separation.\n")
        f.write("# field_degree: degree of the field of the exact description (strict entries described exactly). "
                "live_cos: displayed cosine on 7 November 2026; live_lower: 1 if below the archived displayed cosine; "
                "live_author: 1 if the live table credits the entry to the author. note: further exact results (prestress9.py).\n")
        f.write("\t".join(head9) + "\n")
        for k in keys:
            r = R8[k]
            ii = info[k]
            ng = neg.get(k, {})
            row9 = [str(k[0]), str(k[1]), r["flagged"], r["optimal"], cls[k][0], cls[k][1], ii["tier"] or "", r["tol"],
                    r["contacts"], r["nonrattlers"], r["c_lo_30"], r["c_hi_30"],
                    "" if ii["dist"] is None else f"{ii['dist']:.4g}", "" if ii["rho"] is None else f"{ii['rho']:.4g}",
                    "" if ii["kappa_half"] is None else f"{ii['kappa_half']:.4g}",
                    "" if ii["inside"] is None else str(int(bool(ii["inside"]))), r["gap_lower"], r["gap_first_pass"],
                    str(int(k in set(close))), str(int(k in dbl)), f"{dbl[k]:.3g}" if k in dbl else "",
                    r["tau_prime"], r["h"], f"{ng['d8']:.6g}" if ng else "",
                    f"{ng['dcut']:.6g}" if ng and ng.get("dcut") is not None else "", f"{ng['d9']:.6g}" if ng else "",
                    str(int(ng["res9"])) if ng else "", str(deg.get(k) or "") if k in set(exk) else "",
                    live[k][0] if k in live else "", str(int(lower(k))) if k in live else "",
                    str(int(author(k))) if k in live else "",
                    (f"not infinitesimally jammed: {ps9[k]['flexes_beyond_rotations']} exact flexes beyond the rotations"
                     if k in set(notinf) else "")]
            f.write("\t".join(row9) + "\n")

    # ------------------------------------------------------------------- write
    import re
    # the paper's sources (not part of the code release): every macro they use must be defined
    texfiles = [f for f in ("tex/body_v10.tex", "abstract_v10.tex", "supp_v10.tex") if os.path.exists(os.path.join(HERE, f))]
    texall = "".join(open(os.path.join(HERE, f)).read() for f in texfiles)
    known_latex = {"R", "Q", "F", "Sph", "QVrepo", "QVCommit", "VenueEndMatter", "Gamma", "Phi", "Delta", "Lambda", "Theta",
                   "Pi", "Psi", "Omega", "Bigl", "Bigr", "Big", "TeX", "LaTeX", "S"}
    undefined = sorted(n for n in set(re.findall(r"\\([A-Z][A-Za-z]+)", texall)) if n not in M and n not in known_latex)
    if undefined and not DRAFT:
        sys.exit("macros used in the text but not defined, nothing written: " + ", ".join(undefined))
    if undefined:
        print("DRAFT: undefined macros:", ", ".join(undefined))
    if DRAFT:
        tex = texall
        import re
        for name in set(re.findall(r"\\([A-Za-z]+)", tex)):
            if name not in M and name[0].isupper() and name not in ("R", "Q", "F", "Sph", "QVrepo", "QVCommit",
                                                                    "VenueEndMatter", "Gamma", "Phi", "Delta", "Lambda",
                                                                    "Theta", "Pi", "Psi", "Omega"):
                put(name, "\\textbf{??}", "DRAFT placeholder")
        for t in ("taxonomy", "census", "structure", "neg", "live", "above", "notable", "recog9", "known", "fields"):
            p = os.path.join(HERE, f"v10_tab_{t}.tex")
            if not os.path.exists(p):
                open(p, "w").write("% draft\n\\begin{tabular}{l}draft table " + t + "\\end{tabular}\n")
    with open(os.path.join(HERE, "v10_numbers.tex"), "w") as f:
        f.write("% written by make_numbers_v10.py; do not edit\n")
        for k_ in sorted(M):
            cmd = "providecommand" if SRC[k_] == "DRAFT placeholder" else "newcommand"
            f.write(f"\\{cmd}{{\\{k_}}}{{{M[k_]}}}\n")
    json.dump(M, open(os.path.join(HERE, "v10_numbers.json"), "w"), indent=1, sort_keys=True)
    with open(os.path.join(HERE, "VERIFY_numbers_v10.md"), "w") as f:
        f.write("# Every macro of v10_numbers.tex: value and source\n\nWritten by make_numbers_v10.py; do not edit.\n\n")
        f.write("| macro | value | source |\n|---|---|---|\n")
        for k_ in sorted(M):
            v = str(M[k_]).replace("|", "\\|")
            f.write(f"| `\\{k_}` | {v} | {SRC[k_]} |\n")
    print(len(M), "macros; classes", dict(C), "missing" if missing else "", missing)


if __name__ == "__main__":
    main()
