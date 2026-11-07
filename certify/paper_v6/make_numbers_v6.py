"""Every number and table of paper_v6.tex (census version 6), from released outputs and the census_v6 runs.

  nice -n 10 python3 -I make_numbers_v6.py

Reads, from the repository (~/quantavane-spherical-codes/certify, or $QV_CERTIFY):
  census/index.tsv, census/results_v4/*.json (one file per entry), census/strict_v4/*.json,
  census/exact_v4/{strict,big}/*.json, census/validate_v4.tsv, census/summary_v4.json (only for the macros listed in
  FROM_V4, written by census/summarize.py from the same per-entry files), v5/out/*.json, v5/out/summary_v5.json (the
  class of every entry, written by v5/summarize_v5.py), v5/live/live_20261107.tsv;
and, from this folder (census_v6, proposed for certify/paper_v6/ of the repository):
  checks/same_gram.json, checks/split_halves.json, recount/diag_fail*.jsonl (diag_fail.py),
  phase2/survey_*.jsonl (sym_survey.py), phase2/out_open/*.json, phase2/out_controls/*.json, phase2/out_red_*/*.json
  (symstrict.py).
Writes v6_numbers.tex, v6_numbers.json, v6_tab_census.tex, v6_tab_quant.tex, v6_tab_flow.tex, v6_tab_verdicts.tex,
and VERIFY_numbers.md (every macro, its value and its source).

Class changes against version 5 (summary_v5.json):
  - (22,116): the second-order verdict of the separate note is not counted here (class 'no verdict, attempts failed');
  - every open hyperstatic entry for which symstrict.py reports STRICT moves to 'strict' (kind 'strict_sym');
  - the entries with a FLEX certificate get the negative verdict of the flexibility theorem (no local minimum within
    delta of the stored code); this needs no new run, only the census outputs.
"""
import collections, glob, json, math, os, sys
from fractions import Fraction
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
CEN = os.path.join(CERT, "census")
OUT5 = os.path.join(CERT, "v5", "out")
sys.path.insert(0, os.environ.get("QV_PYLIB", os.path.expanduser("~/Library/Caches/quantavane-census/pylib")))
sys.path.insert(0, CEN)
import numpy as np
import snapshot

RANGES = [(3, 5), (6, 10), (11, 17), (18, 24), (25, 32)]
M, SRC = {}, {}
U53 = 2.0 ** -53


def put(name, value, src):
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
    """x rounded UP in the printed digit (for upper bounds)."""
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
    """x rounded DOWN in the printed digit (for lower bounds)."""
    if x is None or x == 0:
        return sci(x, digits)
    e = int(math.floor(math.log10(abs(x))))
    m = math.floor(x / 10 ** e * 10 ** digits + 1e-9) / 10 ** digits
    if 0 <= e <= 2:
        return sci(m * 10 ** e, digits)
    return f"${m:.{digits}f}\\cdot10^{{{e}}}$"


def quant(vals, ps=(0, 0.1, 0.5, 0.9, 1)):
    """Order statistics: the value of rank floor(p (n - 1)) (0-based) of the sorted values; the median of an even
    number of values is the lower middle one."""
    v = sorted(vals)
    if not v:
        return [None] * len(ps)
    return [v[int(math.floor(p * (len(v) - 1) + 1e-12))] for p in ps]


def elist(ks):
    return ", ".join(f"$({d},{n})$" for d, n in ks)


def load_dir(path):
    out = {}
    for f in glob.glob(os.path.join(path, "*.json")):
        j = json.load(open(f))
        out[(j["d"], j["n"])] = j
    return out


def main():
    m4 = json.load(open(os.path.join(CEN, "summary_v4.json")))["macros"]
    s5 = json.load(open(os.path.join(OUT5, "summary_v5.json")))
    m5 = s5["macros"]
    rows = list(snapshot.load_index())
    idx = {(e["d"], e["n"]): e for e in rows if e["listed"]}
    keys = sorted(k for k in idx if 3 <= k[0] <= 32 and k[1] <= 1024)
    S_IDX = "census/index.tsv (snapshot.py from the archived data set)"
    # ------------------------------------------------------------------------------------------- the data set
    put("NCodesData", num(len(rows)), S_IDX + ": codes in data.txt")
    put("NListed", num(len(idx)), S_IDX + ": codes listed in table.txt")
    put("NNotListed", num(len(rows) - len(idx)), S_IDX + ": codes not listed")
    put("NPolyListed", num(sum(1 for e in idx.values() if e["poly"])), S_IDX + ": listed codes with a minimal polynomial")
    put("NCensus", num(len(keys)), S_IDX + ": listed, 3<=d<=32, N<=1024")
    put("NFlagAll", num(sum(1 for k in keys if idx[k]["flagged"])), S_IDX + ": census entries marked 'needs more optimization'")
    put("NUnflagAll", num(sum(1 for k in keys if not idx[k]["flagged"])), S_IDX + ": census entries not so marked")
    put("NOptAll", num(sum(1 for k in keys if idx[k]["optimal"])), S_IDX + ": census entries marked 'Provably optimal'")
    put("NPolyCensus", num(sum(1 for k in keys if idx[k]["poly"])), S_IDX + ": census entries with a minimal polynomial")
    assert M["NCensus"] == m4["NCensus"] and M["NPolyCensus"] == m4["NExact"]
    # ------------------------------------------------------------------------------------------- per-entry outputs
    res = load_dir(os.path.join(CEN, "results_v4"))
    S_RES_ROOT = "census/results_v4/*.json (top-level record)"
    strict4 = load_dir(os.path.join(CEN, "strict_v4"))
    sym = load_dir(os.path.join(HERE, "phase2", "out_open"))
    assert sorted(res) == sorted(set(res)) and set(keys) <= set(res)
    # ------------------------------------------------------------------------------------------- classes (version 6)
    cls5 = {tuple(map(int, k.split(","))): tuple(v) for k, v in s5["class"].items()}
    assert sorted(cls5) == keys
    cls = dict(cls5)
    cls[(22, 116)] = ("none_tried", "failed")            # the note's verdict is not counted (item M5)
    symstrict = sorted(k for k, j in sym.items() if j.get("status") == "STRICT")
    for k in symstrict:
        assert cls5[k][0] == "local", k
        cls[k] = ("strict", "strict_sym")
    C = collections.Counter(v[0] for v in cls.values())
    KC = collections.Counter(v[1] for v in cls.values())
    S_CLS = ("class of every entry: v5/out/summary_v5.json (v5/summarize_v5.py), with (22,116) moved to 'no verdict' "
             "and the STRICT outputs of phase2/out_open (symstrict.py) moved to 'strict'")
    for name, val in (("NStrict", C["strict"]), ("NLocal", C["local"]), ("NNonStrict", C["nonstrict"]),
                      ("NNotLocal", C["notlocal"]), ("NNoneTried", C["none_tried"]), ("NTooLarge", C["none_size"]),
                      ("NStrictIso", KC["strict_iso"]), ("NStrictHyper", KC["strict_hyper"]),
                      ("NStrictSym", KC["strict_sym"]), ("NStrictExact", KC["strict_exact"]),
                      ("NStrictNew", KC["strict_new"]), ("NNotLocalFlat", KC["notlocal_flat"]),
                      ("NNotLocalNonopt", KC["notlocal_nonopt"]), ("NFlex", KC["flex"]), ("NFailedAll", KC["failed"])):
        put(name, num(val), S_CLS)
    assert KC["strict_second"] == 0
    put("NStrictVFive", m5["VNStrictAll"], "v5/out/summary_v5.json VNStrictAll (version 5 count, for the comparison)")
    put("NLocalVFive", m5["VNLocal"], "v5/out/summary_v5.json VNLocal")
    dec_cls = ("strict", "local", "nonstrict", "notlocal")
    ndec = sum(C[c] for c in dec_cls)
    put("NAnyVerdict", num(ndec), S_CLS)
    put("PctAnyVerdict", f"{100 * ndec / len(keys):.1f}", S_CLS)
    put("NLocalAtLeast", num(C["strict"] + C["local"] + C["nonstrict"]), S_CLS)
    put("NStrictOpenSec", num(KC["strict_new"]), S_CLS + " (exact core, bar-rigid)")
    unm = [k for k in keys if not idx[k]["flagged"]]
    mar = [k for k in keys if idx[k]["flagged"]]
    CU = collections.Counter(cls[k][0] for k in unm)
    CM = collections.Counter(cls[k][0] for k in mar)
    nu_dec = sum(CU[c] for c in dec_cls)
    put("UAnyVerdict", num(nu_dec), S_CLS + " restricted to unmarked entries")
    put("PctUAnyVerdict", f"{100 * nu_dec / len(unm):.0f}", S_CLS + " restricted to unmarked entries")
    put("UStrict", num(CU["strict"]), S_CLS + " (unmarked)")
    put("UFlex", num(sum(1 for k in unm if cls[k][1] == "flex")), S_CLS + " (unmarked)")
    put("MAnyVerdict", num(sum(CM[c] for c in dec_cls)), S_CLS + " (marked)")
    opt = [k for k in keys if idx[k]["optimal"]]
    put("NOptStrict", num(sum(1 for k in opt if cls[k][0] == "strict")), S_CLS + " and index 'optimal'")
    put("OptMissingList", elist([k for k in opt if cls[k][0] != "strict"]), S_CLS + " and index 'optimal'")
    put("StrictExactList", elist(sorted(k for k in keys if cls[k][1] == "strict_exact")), S_CLS)
    put("StrictNewList", elist(sorted(k for k in keys if cls[k][1] == "strict_new")), S_CLS)
    # ------------------------------------------------------------------------------------------- the census table
    ranges = RANGES + [(3, 32)]
    sets = [[k for k in keys if lo <= k[0] <= hi] for lo, hi in ranges]

    def cnt(ks, c, kd):
        return sum(1 for k in ks if cls[k][0] == c and (kd is None or cls[k][1] == kd))

    def line(lab, c, kd):
        return lab + " & " + " & ".join(num(cnt(s, c, kd)) for s in sets) + f" & {num(cnt(unm, c, kd))} \\\\"
    L = ["% written by make_numbers_v6.py; do not edit", "\\begin{tabular}{lrrrrrrr}", "\\toprule",
         " & \\multicolumn{5}{c}{dimension $d$} & & \\\\ \\cmidrule(lr){2-6}",
         " & 3--5 & 6--10 & 11--17 & 18--24 & 25--32 & all & unmarked \\\\", "\\midrule",
         "entries & " + " & ".join(num(len(s)) for s in sets) + f" & {num(len(unm))} \\\\", "\\midrule",
         line("strict local minimum", "strict", None),
         line("\\quad isostatic (Thm.~\\ref{thm:check})", "strict", "strict_iso"),
         line("\\quad recognised Gram matrix (Sec.~\\ref{sec:recog})", "strict", "strict_hyper"),
         line("\\quad symmetry (Thm.~\\ref{thm:sym})", "strict", "strict_sym"),
         line("\\quad exact arithmetic (Prop.~\\ref{prop:exact})", "strict", "strict_exact"),
         line("\\quad exact core, bar-rigid (Sec.~\\ref{sec:open})", "strict", "strict_new"),
         line("local min., strictness open (Thm.~\\ref{thm:hyper})", "local", None),
         line("local min., not strict (Prop.~\\ref{prop:split})", "nonstrict", None),
         line("not a local minimum", "notlocal", None),
         line("\\quad negative stress (Prop.~\\ref{prop:nonopt})", "notlocal", "notlocal_nonopt"),
         line("\\quad flat core (Prop.~\\ref{prop:flat})", "notlocal", "notlocal_flat"),
         "\\midrule",
         line("no local min.\\ within $\\delta$ (Thm.~\\ref{thm:flex})", "none_tried", "flex"),
         "\\multicolumn{8}{l}{no verdict:} \\\\",
         line("\\quad attempts failed or did not apply", "none_tried", "failed"),
         line("\\quad above the size limit", "none_size", None),
         "\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v6_tab_census.tex"), "w").write("\n".join(L) + "\n")
    # program versions: the certifiers have one SHA-256 each; the driver two
    shas = collections.Counter()
    drv = collections.Counter()
    for k in keys:
        cs = res[k]["code_sha256"]
        shas[(cs["cq.py"], cs["sphcert.py"])] += 1
        drv[cs["run_census.py"]] += 1
    assert len(shas) == 1, shas
    (da, na), (db, nb) = sorted(drv.items(), key=lambda x: -x[1])
    put("DriverShaA", da[:6] + "\\ldots", S_RES_ROOT + ": code_sha256['run_census.py'], the more frequent")
    put("NDriverA", num(na), S_RES_ROOT + ": entries with DriverShaA")
    put("DriverShaB", db[:6] + "\\ldots", S_RES_ROOT + ": code_sha256['run_census.py'], the other")
    put("NDriverB", num(nb), S_RES_ROOT + ": entries with DriverShaB")
    assert M["NDriverB"] == m4["NWatchdogEntries"]
    put("CertShaCq", list(shas)[0][0][:6] + "\\ldots", S_RES_ROOT + ": code_sha256['cq.py'] (one value in all entries)")
    # ------------------------------------------------------------------------------------------- isostatic certificates
    S_RES = "census/results_v4/*.json (cq.py via run_census.py), the attempt with status "
    first = {}
    for k in keys:
        st = res[k]["status"]
        att = {a.get("status"): a for a in res[k]["attempts"]}
        if st in ("PASS", "LOCALMIN"):
            first[k] = att[st]
    passes = [k for k in keys if res[k]["status"] == "PASS"]
    locs = [k for k in keys if res[k]["status"] == "LOCALMIN"]
    assert num(len(passes)) == M["NStrictIso"] and num(len(locs)) == m4["NLocalAll"]
    A = {k: first[k] for k in passes}
    qrows = [("kappa_half", "growth rate $\\kappa/(2\\sqrt d)$"), ("rho", "radius $\\rho$"),
             ("dist_ls", "distance $D(X_V)$ to the stored code"),
             ("r_unique", "uniqueness radius $r_u$"),
             ("krawczyk_beta", "largest row sum $\\beta_K$ in (K)"), ("stress_beta", "$\\beta_S$ in (S)"),
             ("cond_inf_upper", "condition bound"), ("stress_min_rel", "smallest relative stress $|E|\\,\\omega_-$"),
             ("unknowns", "unknowns $U$"), ("cpu_seconds", "CPU time (s)")]
    for k in passes:
        kp = A[k].get("kappa_point")
        A[k]["kappa_half"] = None if kp is None else kp / 2
    T = ["% written by make_numbers_v6.py; do not edit", "\\begin{tabular}{lrrrrr}", "\\toprule",
         " & min & 10\\% & median & 90\\% & max \\\\", "\\midrule"]
    for key, lab in qrows:
        vals = [float(A[k][key]) for k in passes if A[k].get(key) is not None]
        qq = quant(vals)
        if key == "unknowns":
            cells = [num(int(x)) for x in qq]
        elif key == "cpu_seconds":
            cells = [f"${x:.2f}$" if x < 10 else f"${x:.0f}$" for x in qq]
        else:
            cells = [sci(x) for x in qq]
        T.append(lab + " & " + " & ".join(cells) + " \\\\")
    T += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v6_tab_quant.tex"), "w").write("\n".join(T) + "\n")
    put("NRhoDefined", num(sum(1 for k in passes if A[k].get("rho") is not None)), S_RES + "PASS: rho present")
    put("NNoRho", num(sum(1 for k in passes if A[k].get("rho") is None)), S_RES + "PASS: rho absent")
    put("NInside", num(sum(1 for k in passes if A[k].get("code_inside_rho") is True)), S_RES + "PASS: code_inside_rho")
    put("NOutside", num(sum(1 for k in passes if A[k].get("code_inside_rho") is False)), S_RES + "PASS: code_inside_rho false")
    put("NUniqueCoversCode", num(sum(1 for k in passes if A[k]["r_unique"] > A[k]["dist_nonrattlers"])),
        S_RES + "PASS: r_unique > dist_nonrattlers")
    for a in ("NRhoDefined", "NInside", "NOutside", "NNoRho", "NUniqueCoversCode"):
        assert M[a] == m4[a], a
    q = quant([A[k]["kappa_half"] for k in passes if A[k].get("kappa_half") is not None])
    put("QKappaHalfMed", sci(q[2]), S_RES + "PASS: kappa_point/2, median")
    q = quant([A[k]["rho"] for k in passes if A[k].get("rho") is not None])
    put("QRhoMed", sci(q[2]), S_RES + "PASS: rho, median")
    put("QCErrMax", sci_up(max(A[k]["c_err"] for k in passes)), S_RES + "PASS: max c_err")
    put("QRUniqueMed", sci(quant([A[k]["r_unique"] for k in passes])[2]), S_RES + "PASS: r_unique, median")
    dls = [A[k]["dist_ls"] for k in passes]
    put("PassDistMax", sci_up(max(dls)), S_RES + "PASS: max dist_ls")
    put("PassDistMed", sci(quant(dls)[2]), S_RES + "PASS: dist_ls, median")
    put("NPassDistAbove", num(sum(1 for x in dls if x > 1e-8)), S_RES + "PASS: dist_ls > 1e-8")
    # rattlers (item m1): moves recorded by sphcert.place_rattlers in 'rattler_placement'
    withr = [k for k in passes if (A[k].get("rattler_placement") or {}).get("rattlers", 0) > 0]
    moved = [k for k in passes if (A[k].get("rattler_placement") or {}).get("moved_first_order", 0)
             + (A[k].get("rattler_placement") or {}).get("moved_second_order", 0) > 0]
    put("NPassRattlers", num(len(withr)), S_RES + "PASS: rattler_placement.rattlers > 0")
    put("NRattlerMoved", num(len(moved)), S_RES + "PASS: rattler_placement moved_first_order + moved_second_order > 0")
    put("RattlerMovedMax", sci_up(max(A[k]["dist_rattlers"] for k in moved), 2), S_RES + "PASS: max dist_rattlers over moved (2 digits, rounded up)")
    put("NRattlerStuck", num(sum(1 for k in passes if (A[k].get("rattler_placement") or {}).get("stuck", 0) > 0)),
        S_RES + "PASS: rattler_placement.stuck > 0")
    # timings
    cpuP = [A[k]["cpu_seconds"] for k in passes]
    cpuL = [first[k]["cpu_seconds"] for k in locs]
    put("CPUMedPass", f"{quant(cpuP)[2]:.2f}", S_RES + "PASS: cpu_seconds median")
    put("CPUMaxPass", f"{max(cpuP):.0f}", S_RES + "PASS: cpu_seconds max")
    put("CPUMedLoc", f"{quant(cpuL)[2]:.2f}", S_RES + "LOCALMIN: cpu_seconds median")
    put("CPUMaxLoc", f"{max(cpuL):.0f}", S_RES + "LOCALMIN: cpu_seconds max")
    big = [k for k in passes + locs if first[k]["unknowns"] >= 1000]
    Uu = np.array([first[k]["unknowns"] for k in big], float)
    tt = np.array([first[k]["cpu_seconds"] for k in big], float)
    bfit, afit = np.polyfit(np.log(Uu), np.log(tt), 1)
    put("CPUSlope", f"{bfit:.1f}", S_RES + "PASS/LOCALMIN with >=1000 unknowns: least-squares slope of log cpu on log U")
    put("NCPUSlope", num(len(big)), S_RES + "PASS/LOCALMIN with >=1000 unknowns")
    sh = [first[k]["seconds_inverse"] / first[k]["wall_seconds"] for k in big
          if first[k].get("seconds_inverse") and first[k].get("wall_seconds")]
    put("InvShareMed", f"{100 * quant(sh)[2]:.0f}", S_RES + "PASS/LOCALMIN >=1000 unknowns: seconds_inverse/wall_seconds, median (wall-clock share)")
    kbig = max(passes + locs, key=lambda k: first[k]["unknowns"])
    put("BigU", num(first[kbig]["unknowns"]), S_RES + "largest certified system")
    put("BigCPU", f"{first[kbig]['cpu_seconds']:.0f}", S_RES + "largest certified system: cpu_seconds")
    put("BigEntry", f"$({kbig[0]},{kbig[1]})$", S_RES + "largest certified system")
    put("MaxUnknownsPass", num(max(A[k]["unknowns"] for k in passes)), S_RES + "PASS: max unknowns")
    put("MaxUnknownsLoc", num(max(first[k]["unknowns"] for k in locs)), S_RES + "LOCALMIN: max unknowns")
    tolP = collections.Counter(A[k]["tol"] for k in passes)
    tolL = collections.Counter(first[k]["tol"] for k in locs)
    put("NPassTolOther", num(len(passes) - tolP[1e-10]), S_RES + "PASS: tol != 1e-10")
    put("NLocTolOther", num(len(locs) - tolL[1e-10]), S_RES + "LOCALMIN: tol != 1e-10")
    assert M["NPassTolOther"] == m4["NPassTolOther"] and M["NLocTolOther"] == m4["NLocTolOther"]
    # verdicts at several tolerances (question 4): a NONOPT attempt followed by PASS or LOCALMIN at another tolerance
    multi = [k for k in keys if any(a.get("status") == "NONOPT" for a in res[k]["attempts"])
             and res[k]["status"] in ("PASS", "LOCALMIN")]
    put("NNonoptThenCert", num(len(multi)), S_RES + "an attempt NONOPT and a later PASS/LOCALMIN")
    assert not multi          # the text says that no entry has both
    # ------------------------------------------------------------------------------------------- proposals and attempts
    prop = {}
    for k in keys:
        att = [a for a in res[k]["attempts"] if a.get("kind") in ("isostatic", "hyperstatic")]
        sk = [a for a in res[k].get("skipped_above_cap", []) if a.get("kind") in ("iso", "hyper")]
        prop[k] = "attempted" if att else ("abovecap" if sk else "never")
    PC = collections.Counter(prop.values())
    S_PROP = S_RES + "(attempts of kind isostatic/hyperstatic; skipped_above_cap of kind iso/hyper)"
    put("NPropNever", num(PC["never"]), S_PROP + ": none, none")
    put("NPropAboveCap", num(PC["abovecap"]), S_PROP + ": none attempted, some skipped")
    put("NPropAttempted", num(PC["attempted"]), S_PROP + ": some attempted")
    put("PctPropNever", f"{100 * PC['never'] / len(keys):.0f}", S_PROP)
    put("MPropNever", num(sum(1 for k in keys if prop[k] == "never" and idx[k]["flagged"])), S_PROP + " (marked)")
    put("UPropNever", num(sum(1 for k in keys if prop[k] == "never" and not idx[k]["flagged"])), S_PROP + " (unmarked)")
    att_k = [k for k in keys if prop[k] == "attempted"]
    verdict = lambda k: cls[k][0] in dec_cls
    nv = sum(1 for k in att_k if verdict(k))
    put("NPropAttVerdict", num(nv), S_PROP + " and class")
    put("PctPropAttVerdict", f"{100 * nv / len(att_k):.0f}", S_PROP + " and class")
    ua = [k for k in att_k if not idx[k]["flagged"]]
    put("UPropAttempted", num(len(ua)), S_PROP + " (unmarked)")
    put("UPropAttVerdict", num(sum(1 for k in ua if verdict(k))), S_PROP + " and class (unmarked)")
    put("PctUPropAttVerdict", f"{100 * sum(1 for k in ua if verdict(k)) / len(ua):.0f}", S_PROP + " and class (unmarked)")
    assert all(not verdict(k) for k in keys if prop[k] == "never")
    put("NTooLargeNever", num(sum(1 for k in keys if prop[k] == "never" and cls[k][0] == "none_size")), S_PROP)
    fc = collections.Counter()
    failed_exist = []
    for k in att_k:
        if verdict(k):
            continue
        a0 = [a for a in res[k]["attempts"] if a.get("kind") in ("isostatic", "hyperstatic")][0]["status"]
        tag = ("exist" if a0 in ("Jacobian numerically singular", "FAIL existence (Krawczyk)") else
               ("nostress" if "no positive stress" in a0 else "other"))
        fc[tag] += 1
        if tag == "exist":
            failed_exist.append(k)
    put("NAttemptFail", num(len(att_k) - nv), S_PROP + " without a verdict")
    put("NAttFailExist", num(fc["exist"]), S_PROP + ": first attempt failed existence")
    put("NAttFailNoStress", num(fc["nostress"]), S_PROP + ": first attempt found no positive stress")
    put("NAttFailOther", num(fc["other"]), S_PROP + ": first attempt failed otherwise")
    ab = [k for k in keys if prop[k] == "abovecap" and cls[k][0] in ("none_size", "none_tried")]
    put("NAboveCapOpen", num(len(ab)), S_PROP)
    put("NAboveCapExact", num(sum(1 for k in keys if prop[k] == "abovecap" and cls[k][0] == "strict")), S_PROP)
    Uab = [min(a["size"] for a in res[k]["skipped_above_cap"] if a.get("kind") in ("iso", "hyper")) for k in ab]
    put("AboveCapMaxU", num(max(Uab)), S_PROP + ": smallest skipped system, max over entries")
    put("AboveCapGB", f"{8 * max(Uab) ** 2 / 1e9:.1f}", "8 AboveCapMaxU^2 bytes")
    bigk = [k for k in keys if k in first and first[k]["unknowns"] >= 1000]
    UU = np.array([first[k]["unknowns"] for k in bigk], float)
    t2 = np.array([first[k]["cpu_seconds"] for k in bigk], float)
    bf, af = np.polyfit(np.log(UU), np.log(t2), 1)
    put("AboveCapCPUh", f"{sum(math.exp(af) * u ** bf for u in Uab) / 3600:.0f}", "extrapolation of the fit CPUSlope to the skipped systems")
    # flow table (item m13): proposals against verdicts
    def vgroup(k):
        c, kd = cls[k]
        if c in dec_cls:
            return "verdict"
        if kd == "flex":
            return "flex"
        return "none"
    FL = collections.Counter((prop[k], vgroup(k)) for k in keys)
    labels = [("attempted", "a certificate attempted"), ("abovecap", "certifiable proposal only above the size limit"),
              ("never", "no isostatic or hyperstatic proposal")]
    TF = ["% written by make_numbers_v6.py; do not edit", "\\begin{tabular}{lrrrr}", "\\toprule",
          "proposal at some $\\tau\\in\\{10^{-13},\\dots,10^{-7}\\}$ & verdict & no local min.\\ within $\\delta$ & none & all \\\\",
          "\\midrule"]
    for p_, lab in labels:
        TF.append(f"{lab} & {num(FL[(p_, 'verdict')])} & {num(FL[(p_, 'flex')])} & {num(FL[(p_, 'none')])} & "
                  f"{num(sum(FL[(p_, g)] for g in ('verdict', 'flex', 'none')))} \\\\")
    TF.append("\\midrule")
    TF.append(f"all & {num(sum(FL[(p_, 'verdict')] for p_, _ in labels))} & {num(sum(FL[(p_, 'flex')] for p_, _ in labels))} & "
              f"{num(sum(FL[(p_, 'none')] for p_, _ in labels))} & {num(len(keys))} \\\\")
    TF += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v6_tab_flow.tex"), "w").write("\n".join(TF) + "\n")
    for (p_, g), v in FL.items():
        put(f"Flow{p_.capitalize()}{g.capitalize()}", num(v), S_PROP + " x class")
    # ------------------------------------------------------------------------------------------- diagnostics of failures (M3)
    diag = {}
    for f in glob.glob(os.path.join(HERE, "recount", "diag_fail*.jsonl")):
        for line_ in open(f):
            j = json.loads(line_)
            diag[(j["d"], j["n"])] = j
    S_D = "recount/diag_fail*.jsonl (recount/diag_fail.py, double precision, untrusted)"
    fe = [k for k in failed_exist]
    assert set(fe) <= set(diag), sorted(set(fe) - set(diag))[:5]
    D_ = [diag[k] for k in fe]
    put("NDiag", num(len(D_)), S_D)
    sing = lambda j: j["rel_sv_min"] < 1e-12
    sep = lambda j: j.get("gap_next") is not None and j["gap_next"] > 10 * j["tol"]
    conv = lambda j: j.get("newton_last") is not None and j["newton_last"] < 1e-25
    put("NDiagSingular", num(sum(1 for j in D_ if sing(j))), S_D + ": smallest/largest singular value < 1e-12")
    put("NDiagNonsing", num(sum(1 for j in D_ if not sing(j))), S_D + ": smallest/largest singular value >= 1e-12")
    put("NDiagSep", num(sum(1 for j in D_ if sep(j))), S_D + ": next pair more than 10 tau below the largest cosine")
    put("NDiagSepSing", num(sum(1 for j in D_ if sep(j) and sing(j))), S_D + ": separated and numerically singular")
    put("NDiagSepNonsing", num(sum(1 for j in D_ if sep(j) and not sing(j))), S_D + ": separated and nonsingular")
    put("NDiagClose", num(sum(1 for j in D_ if not sep(j))), S_D + ": next pair within 10 tau")
    put("NDiagCloseSing", num(sum(1 for j in D_ if not sep(j) and sing(j))), S_D + ": not separated and numerically singular")
    close_ = [j for j in D_ if not sep(j)]
    put("DiagGapRatioMed", f"{quant([j['gap_next'] / j['tol'] for j in close_ if j.get('gap_next') is not None])[2]:.2f}",
        S_D + ": over the unseparated proposals, median of (largest cosine - next pair) / tau")
    put("DiagSpreadRatioMed", f"{quant([j['spread_E'] / j['tol'] for j in close_])[2]:.2f}",
        S_D + ": over the unseparated proposals, median of (largest - smallest contact cosine) / tau")
    put("NDiagConv", num(sum(1 for j in D_ if conv(j))), S_D + ": census Newton residual below 1e-25")
    put("NDiagNoConv", num(sum(1 for j in D_ if not conv(j))), S_D + ": census Newton residual not below 1e-25")
    put("NDiagCorankOne", num(sum(1 for j in D_ if sing(j) and j["n_below_1e12"] == 1)), S_D + ": exactly one singular value below 1e-12 relative")
    rs = [j["krawczyk_rowsum"] for j in D_ if j.get("krawczyk_rowsum") is not None]
    put("DiagRowsumMin", sci_dn(min(rs)), S_D + " (census record krawczyk_attempts[0].max_rowsum): min")
    put("DiagRowsumMed", sci(quant(rs)[2]), S_D + ": median")
    put("NDiagRowsumBelowOne", num(sum(1 for x in rs if x < 1)), S_D + ": contraction factor below 1")
    # ------------------------------------------------------------------------------------------- hyperstatic certificates
    put("NLocalAll", num(len(locs)), S_RES + "LOCALMIN")
    put("ChiMax", sci_up(max(first[k]["tau_extra"] for k in locs)), S_RES + "LOCALMIN: max tau_extra")
    put("TubeMax", sci_up(max(first[k]["localmin_displacement"] for k in locs)), S_RES + "LOCALMIN: max localmin_displacement")
    ldist = [first[k]["dist_ls"] + first[k]["localmin_displacement"] for k in locs]
    put("LocDistMax", sci_up(max(first[k]["dist_ls"] for k in locs)), S_RES + "LOCALMIN: max dist_ls (to P')")
    put("LocDistTubeMax", sci_up(max(ldist)), S_RES + "LOCALMIN: max dist_ls + localmin_displacement")
    ex_ = [first[k]["contacts"] - first[k]["isostatic_count"] for k in locs]
    put("LocExtraMin", num(min(ex_)), S_RES + "LOCALMIN: min contacts - isostatic_count")
    put("LocExtraMed", num(quant(ex_)[2]), S_RES + "LOCALMIN: median")
    put("LocExtraMax", num(max(ex_)), S_RES + "LOCALMIN: max")
    S_ST = "census/strict_v4/*.json (strict.py)"
    st_ok = {k: j for k, j in strict4.items() if j.get("status") == "STRICT"}
    assert num(len(st_ok)) == M["NStrictHyper"]
    deg = collections.Counter(j.get("field_degree") for j in st_ok.values())
    put("NStrictHyperFieldOne", num(deg[1]), S_ST + ": field_degree 1")
    put("NStrictHyperFieldTwo", num(deg[2]), S_ST + ": field_degree 2")
    put("NStrictHyperFieldHigh", num(sum(v for d_, v in deg.items() if d_ and d_ >= 3)), S_ST + ": field_degree >= 3")
    sd = []
    nan_fix = []
    for k, j in st_ok.items():
        dl = j.get("dist_ls_upper")
        if dl is None or not math.isfinite(dl):
            dl = math.sqrt(j["nV"] * k[0]) * j["dist_coord_upper"] * (1 + 1e-12)
            nan_fix.append(k)
        sd.append(dl)
    put("StrictDistMax", sci_up(max(sd)), S_ST + ": max dist_ls_upper (sqrt(n d) times dist_coord_upper where NaN)")
    put("StrictDistMed", sci(quant(sd)[2]), S_ST + ": dist_ls_upper, median (order statistic)")
    put("StrictNaNList", elist(sorted(nan_fix)), S_ST + ": dist_ls_upper NaN")
    put("NStrictRotExact", num(sum(1 for j in st_ok.values() if j.get("rotation_exactly_orthogonal"))), S_ST + ": rotation_exactly_orthogonal")
    reasons = collections.Counter()
    for k in keys:
        if cls[k][0] != "local":
            continue
        s_ = strict4[k]["status"]
        reasons["nopoly" if s_.startswith("no polynomial") else "notinfield" if s_.startswith("an inner product")
                else "toomany" if s_.startswith("too many") else "degree" if s_.startswith("cosine of degree")
                else "other"] += 1
    put("NOpenNoPoly", num(reasons["nopoly"]), S_ST + " status of the entries still open")
    put("NOpenNotInField", num(reasons["notinfield"]), S_ST + " status of the entries still open")
    put("NOpenTooMany", num(reasons["toomany"]), S_ST + " status of the entries still open")
    put("NOpenDegree", num(reasons["degree"]), S_ST + " status of the entries still open")
    assert reasons["other"] == 0
    # ------------------------------------------------------------------------------------------- validation (m6 of M6)
    vt = [l.rstrip("\n").split("\t") for l in open(os.path.join(CEN, "validate_v4.tsv"))][1:]
    tested = [r for r in vt if r[-1]]
    S_V = "census/validate_v4.tsv (validate_poly.py)"
    put("NExactTested", num(len(tested)), S_V + ": rows with a verdict")
    put("NExactRootInside", num(sum(1 for r in tested if r[-1] == "root inside")), S_V + ": 'root inside'")
    put("NExactTestedIso", num(sum(1 for r in tested if r[5] == "isostatic")), S_V + ": kind isostatic")
    put("NExactTestedHyper", num(sum(1 for r in tested if r[5] == "hyperstatic")), S_V + ": kind hyperstatic")
    hyper_t = [(int(r[0]), int(r[1])) for r in tested if r[5] == "hyperstatic"]
    put("NExactTestedHyperStrict", num(sum(1 for k in hyper_t if cls[k][0] == "strict")), S_V + " and class")
    put("NExactTestedHyperRecog", num(sum(1 for k in hyper_t if cls[k][1] == "strict_hyper")), S_V + " and class strict_hyper (recognition with the table's polynomial)")
    put("ValChiMax", sci_up(max(first[k]["tau_extra"] for k in hyper_t)), S_RES + "LOCALMIN of the validated hyperstatic entries: max tau_extra")
    # ------------------------------------------------------------------------------------------- exact arithmetic, second program
    S_EX = "census/exact_v4/strict/*.json (exact_strict.py)"
    ex = [json.load(open(f)) for f in glob.glob(os.path.join(CEN, "exact_v4", "strict", "*.json"))]
    put("NExactTried", num(len(ex)), S_EX)
    rec_ = [j for j in ex if j["status"] != "coordinates not recognised as exact"]
    put("NExactRecognised", num(len(rec_)), S_EX + ": recognised")
    okx = [j for j in ex if str(j["status"]).startswith("STRICT")]
    put("NExactStrict", num(len(okx)), S_EX + ": STRICT")
    put("ExactCPUMax", f"{max(float(j['cpu_seconds']) for j in okx):.0f}", S_EX + ": max cpu_seconds of STRICT")
    S_BIG = "census/exact_v4/big/*.json (exact_big.py)"
    bigj = [json.load(open(f)) for f in glob.glob(os.path.join(CEN, "exact_v4", "big", "*.json"))]
    for mode, pre in (("iso", "Iso"), ("hyper", "Hyper"), ("strict", "Strict")):
        js = [j for j in bigj if j.get("mode") == mode]
        good = [j for j in js if j["status"].startswith(("PASS", "LOCALMIN confirmed", "STRICT confirmed")) and "FAIL" not in j["status"]]
        put(f"Big{pre}Good", num(len(good)), S_BIG + f": mode {mode}, confirmed")
        put(f"Big{pre}All", num(len(js)), S_BIG + f": mode {mode}")
        uu = [int(j.get("unknowns") or 0) for j in js]
        put(f"Big{pre}MaxU", num(max(uu)) if any(uu) else "--", S_BIG + f": mode {mode}, max unknowns")
    put("BigSecondsMax", f"{max(float(j.get('seconds', 0)) for j in bigj if j.get('mode') != 'flex'):.0f}", S_BIG + ": max seconds")
    # ------------------------------------------------------------------------------------------- version 5 analyses
    S5 = "v5/out/summary_v5.json (v5/summarize_v5.py) macro "
    for a in ("VExOnlyN", "VExOnlyDistLSMax", "VExOnlyDistCoordMax", "VFlatList", "VFlatDistMax", "VFlatColoursMin",
              "VFlatColoursMax", "VNUnion", "VSplitDistMax", "VImpLowerA", "VImpLowerB", "VFamilyList", "VNFlatLive"):
        put(a[1:], m5[a], S5 + a)
    put("NOpenExact", m5["VNOpenExact"], S5 + "VNOpenExact")
    put("NDecidedExact", num(int(m5["VNDecidedExact"]) - 1), S5 + "VNDecidedExact minus (22,116)")
    sub = {}
    for f in glob.glob(os.path.join(OUT5, "subspace_d*_n*.json")):
        j = json.load(open(f))
        sub[(j["d"], j["n"])] = j
    flat = sorted(k for k in sub if sub[k]["status"].startswith("NOT A LOCAL MINIMUM"))
    assert num(len(flat)) == M["NNotLocalFlat"]
    gains = {k: -sub[k]["core_max_minus_c"] for k in flat}
    kmin, kmax = min(gains, key=gains.get), max(gains, key=gains.get)
    put("TiltMin", sci(gains[kmin]), "v5/out/subspace_*.json: -core_max_minus_c, min")
    put("TiltMinEntry", f"$({kmin[0]},{kmin[1]})$", "v5/out/subspace_*.json")
    put("TiltMax", sci(gains[kmax]), "v5/out/subspace_*.json: -core_max_minus_c, max")
    put("TiltMaxEntry", f"$({kmax[0]},{kmax[1]})$", "v5/out/subspace_*.json")
    rig = {}
    for f in glob.glob(os.path.join(OUT5, "rigidity_d*_n*.json")):
        j = json.load(open(f))
        rig[(j["d"], j["n"])] = j
    S_RIG = "v5/out/rigidity_*.json (v5/rigidity.py)"
    for (d, n), nm in (((8, 72), "A"), ((9, 330), "B")):
        j = rig[(d, n)]
        assert j["status"].startswith("LOCAL, NOT STRICT")
        put(f"SplitRank{nm}", num(j["rank_mod_p"]), S_RIG + f" ({d},{n}): rank_mod_p")
        put(f"SplitCore{nm}", num(j["core_size"]), S_RIG + f" ({d},{n}): core_size")
    assert int(rig[(8, 72)]["rank_mod_p"]) == 2 * (36 * 8 - 28) and int(rig[(8, 72)]["core_size"]) == 72
    assert int(rig[(9, 330)]["rank_mod_p"]) == 2 * (120 * 9 - 36) and int(rig[(9, 330)]["core_size"]) == 240
    assert int(rig[(10, 238)]["rank_mod_p"]) == int(rig[(10, 238)]["cols"]) - 45
    assert int(rig[(12, 210)]["rank_mod_p"]) == int(rig[(12, 210)]["cols"]) - 66
    for (d, n), nm in (((10, 238), "A"), ((12, 210), "B")):
        j = rig[(d, n)]
        assert j["status"].startswith("STRICT")
        put(f"NewRank{nm}", num(j["rank_mod_p"]), S_RIG + f" ({d},{n}): rank_mod_p")
        put(f"NewCols{nm}", num(j["cols"]), S_RIG + f" ({d},{n}): cols")
        put(f"NewGap{nm}", sci_dn(-float(j["rattler_max_minus_c_upper"])), S_RIG + f" ({d},{n}): -rattler_max_minus_c_upper")
        put(f"NewDist{nm}", sci_up(float(j["dist_ls_upper"])), S_RIG + f" ({d},{n}): dist_ls_upper")
    put("NewDistMax", sci_up(max(float(rig[k]["dist_ls_upper"]) for k in ((10, 238), (12, 210)))), S_RIG + ": max dist_ls_upper of (10,238), (12,210)")
    j = rig[(19, 166)]
    put("OpenContacts", num(j["contacts"]), S_RIG + " (19,166)")
    put("OpenStressed", num(j["stressed_contacts"]), S_RIG + " (19,166)")
    put("OpenRankStressed", num(j["rank_mod_p_stressed"]), S_RIG + " (19,166)")
    put("OpenRankNeeded", num(j["cols"] - 19 * 18 // 2), S_RIG + " (19,166): cols - g")
    put("OpenLPMaxMin", sci(abs(j["lp_maxmin_stress"])), S_RIG + " (19,166)")
    d116 = diag[(22, 116)]
    assert d116["rel_sv_min"] < 1e-12 and d116["n_below_1e12"] == 1
    nonopt = sorted(k for k in keys if any(a.get("status") == "NONOPT" for a in res[k]["attempts"]))
    put("NonoptList", elist(nonopt), S_RES + "NONOPT")
    nod = []
    for k, nm in zip(nonopt, ("A", "B")):
        a = [x for x in res[k]["attempts"] if x.get("status") == "NONOPT"][0]
        put(f"NonoptDist{nm}", sci_up(a["dist_ls"]), S_RES + f"NONOPT {k}: dist_ls")
        put(f"NonoptStress{nm}", sci((a.get("stress_quantiles_over_mean") or [None])[0]), S_RES + f"NONOPT {k}: smallest stress over mean")
        nod.append(a["dist_ls"])
    sg = json.load(open(os.path.join(HERE, "checks", "same_gram.json")))["entries"]
    assert num(sum(1 for e in sg if e["status"] == "SAME GRAM MATRIX")) == M["ExOnlyN"]
    shv = json.load(open(os.path.join(HERE, "checks", "split_halves.json")))["entries"]
    assert sorted((e["d"], e["n"]) for e in shv if e["status"] == "HALVES ANTIPODAL, CORE RANK d") == [(8, 72), (9, 330)]
    # ------------------------------------------------------------------------------------------- Arb re-check (m15)
    ap = os.path.join(HERE, "recount", "arbcheck.jsonl")
    if os.path.exists(ap):
        AR = [json.loads(l_) for l_ in open(ap)]
        S_A = "recount/arbcheck.jsonl (recount/arbcheck.py, Arb ball arithmetic)"
        put("NArbRun", num(len(AR)), S_A)
        put("NArbIso", num(sum(1 for r in AR if r["census_status"] == "PASS")), S_A + ": census PASS")
        put("NArbHyper", num(sum(1 for r in AR if r["census_status"] == "LOCALMIN")), S_A + ": census LOCALMIN")
        put("NArbPass", num(sum(1 for r in AR if r.get("arb_pass"))), S_A + ": Krawczyk test passes in Arb")
        put("NArbOverlap", num(sum(1 for r in AR if r.get("overlaps_census_interval"))), S_A + ": cosine ball overlaps the census 30-digit interval")
        put("ArbMaxU", num(max(r["unknowns"] for r in AR)), S_A + ": max unknowns")
    # ------------------------------------------------------------------------------------------- the flexibility theorem (M4)
    S_F = S_RES + "FLEX (cq.certify_flex) and recount/flexgram*.jsonl (recount/flexgram.py): "
    fg = {}
    for f in glob.glob(os.path.join(HERE, "recount", "flexgram*.jsonl")):
        for line_ in open(f):
            r_ = json.loads(line_)
            fg[(r_["d"], r_["n"])] = r_
    fl = []
    for k in keys:
        if cls[k][1] != "flex":
            continue
        a = [x for x in res[k]["attempts"] if x.get("status") == "FLEX"][0]
        if os.environ.get("QV_PARTIAL") and k not in fg:
            continue
        r_ = fg[k]
        assert r_["proposal_matches_census"], k
        assert r_["delta"] > 0 and r_["tau_prime"] > 0, k
        assert r_["h"] == a["flex_radius_l1_stored"] and r_["tol"] == a["tol"]
        fl.append((k, a["tol"], r_["h"], r_["delta"], r_["eps"], r_["tau_prime"]))
    if not fl and os.environ.get("QV_PARTIAL"):            # layout runs only: placeholders
        fl = [(keys[0], 1e-7, 1.0, 1e-9, 1e-15, 1e-7)]
    put("NFlexTheorem", num(len(fl)), S_F + "count")
    put("FlexEpsMax", sci_up(max(x[4] for x in fl)), S_F + "max of the Gram error bound eps")
    tc = collections.Counter(x[1] for x in fl)
    for t_, nm in ((1e-7, "Seven"), (1e-8, "Eight"), (1e-9, "Nine"), (1e-10, "Ten")):
        put(f"NFlexTol{nm}", num(tc[t_]), S_F + f"tol {t_}")
    put("NFlexTolSmaller", num(sum(v for t_, v in tc.items() if t_ < 1e-10)), S_F + "tol < 1e-10")
    dls_ = [x[3] for x in fl]
    q = quant(dls_)
    put("FlexDeltaMin", sci_dn(q[0]), S_F + "delta = min(h/sqrt d, tau'/4), min")
    put("FlexDeltaMed", sci_dn(q[2]), S_F + "delta, median")
    put("FlexDeltaMax", sci_dn(q[4]), S_F + "delta, max")
    put("NFlexDeltaAtTauQuarter", num(sum(1 for k, t_, h, dl, e, tp in fl if h / math.sqrt(k[0]) >= tp / 4)),
        S_F + "delta attained by tau'/4")
    put("FlexTauRatioMin", f"{min(x[5] / x[1] for x in fl):.2f}", S_F + "min of tau'/tau")
    put("NFlexTauBelowTau", num(sum(1 for x in fl if x[5] < x[1])), S_F + "entries with tau' < tau")
    hq = quant([x[2] / math.sqrt(x[0][0]) for x in fl])
    put("FlexHMed", sci(hq[2]), S_F + "h/sqrt d, median")
    put("NFlexUnmarked", num(sum(1 for x in fl if not idx[x[0]]["flagged"])), S_F + "unmarked entries")
    put("NFlexDeltaAboveNine", num(sum(1 for x in fl if x[3] >= 1e-9)), S_F + "delta >= 1e-9")
    # ------------------------------------------------------------------------------------------- displayed values (M7e, Q6)
    below = []
    for k in passes + locs:
        a = first[k]
        hi = Fraction(int(a["c_interval_30"][1]), 10 ** 30) + Fraction(a.get("tau_extra") or 0)
        c12 = Fraction(-(-hi.numerator * 10 ** 12 // hi.denominator), 10 ** 12)
        if c12 < Fraction(Decimal(idx[k]["cos"])):
            below.append(k)
    put("NDispBelow", num(len(below)), S_RES + "PASS/LOCALMIN: certified value (+ tau_extra) rounded up to 12 decimals below the displayed cosine")
    assert M["NDispBelow"] == m4["NDispBelow"]
    live = {}
    for line_ in open(os.path.join(CERT, "v5", "live", "live_20261107.tsv")):
        if line_.startswith("#"):
            continue
        p = line_.rstrip("\n").split("\t")
        live[(int(p[0]), int(p[1]))] = p
    S_L = "v5/live/live_20261107.tsv (the live table read on 7 November 2026)"
    put("NLiveAuthor", num(sum(1 for p in live.values() if "Arhan Canli" in p[4])), S_L + ": credited to the author")
    bl = [k for k in below if k in live]
    put("NDispBelowLiveAuthor", num(sum(1 for k in bl if "Arhan Canli" in live[k][4])), S_L + ": entries of NDispBelow credited to the author")
    assert M["NDispBelowLiveAuthor"] == "0"   # the text says 'credits none to the author'
    put("NDispBelowLiveImproved", num(sum(1 for k in bl if Decimal(live[k][2]) < Decimal(idx[k]["cos"]))),
        S_L + ": entries of NDispBelow whose live cosine is below the archived one")
    # ------------------------------------------------------------------------------------------- Phase 2: symmetry
    S_SV = "phase2/survey_*.jsonl (phase2/sym_survey.py, untrusted)"
    sv = {}
    for f in glob.glob(os.path.join(HERE, "phase2", "survey_*.jsonl")):
        for line_ in open(f):
            j = json.loads(line_)
            sv[(j["d"], j["n"])] = j
    openv5 = [k for k in keys if cls5[k][0] == "local"]
    assert set(openv5) <= set(sv)
    put("NSurvey", num(len(openv5)), S_SV)
    put("NSymNontrivial", num(sum(1 for k in openv5 if (sv[k].get("n_auts_found") or 0) > 1)), S_SV + ": an automorphism other than the identity found")
    put("NSymSquare", num(sum(1 for k in openv5 if sv[k].get("square"))), S_SV + ": reduced system square")
    S_SY = "phase2/out_open/*.json (phase2/symstrict.py)"
    if not os.environ.get("QV_PARTIAL"):
        assert set(openv5) <= set(sym), sorted(set(openv5) - set(sym))[:5]
    else:
        openv5 = [k for k in openv5 if k in sym]
    stc = collections.Counter(sym[k]["status"].split(" (")[0] for k in openv5)
    put("NSymRun", num(len(openv5)), S_SY)
    put("NSymNoSym", num(stc["no symmetry"]), S_SY + ": status 'no symmetry' (search complete, only the identity)")
    put("NSymTimeout", num(stc["no symmetry found within the time limit"]), S_SY + ": no automorphism found within the time limit")
    put("NSymCap", num(stc["group above the enumeration cap"]), S_SY + ": group above the enumeration cap")
    put("NSymNotSquare", num(sum(1 for k in openv5 if sym[k]["status"].startswith("reduced system not square"))), S_SY)
    put("NSymNoModel", num(stc["no exact model"]), S_SY + ": status 'no exact model'")
    put("NSymStrict", num(len(symstrict)), S_SY + ": status STRICT")
    other = [k for k in openv5 if sym[k]["status"] not in ("STRICT", "no symmetry", "no exact model",
             "no symmetry found within the time limit", "group above the enumeration cap")
             and not sym[k]["status"].startswith("reduced system not square")]
    put("NSymOther", num(len(other)), S_SY + ": other failures")
    assert all(sym[k]["status"] == "not LOCALMIN on re-run" for k in other), [sym[k]["status"] for k in other]
    put("SymOtherStatuses", "; ".join(f"{v} {k}" for k, v in collections.Counter(sym[k]["status"] for k in other).most_common()),
        S_SY + ": statuses of the other failures (for the text, not printed if empty)")
    sqs = [k for k in openv5 if sv[k].get("square")]
    put("NSymSquareStrict", num(sum(1 for k in sqs if sym[k]["status"] == "STRICT")), S_SY + " and " + S_SV + ": square in the survey and STRICT")
    put("NSymSquareNoModel", num(sum(1 for k in sqs if sym[k]["status"] == "no exact model")), S_SY + ": square in the survey, no exact model")
    put("SymSquareNoModelList", elist(sorted(k for k in sqs if sym[k]["status"] == "no exact model")), S_SY)
    if symstrict:
        SJ = [sym[k] for k in symstrict]
        put("SymStrictList", elist(symstrict), S_SY)
        put("SymDistMax", sci_up(max(j["dist_ls_upper"] for j in SJ)), S_SY + ": max dist_ls_upper")
        put("SymGammaRadMax", sci_up(max(j["gamma_rad"] for j in SJ)), S_SY + ": max radius of gamma*")
        put("SymGroupMin", num(min(j["group_order_exact"] for j in SJ)), S_SY + ": group_order_exact")
        put("SymGroupMax", num(max(j["group_order_exact"] for j in SJ)), S_SY + ": group_order_exact")
        put("SymRedMax", num(max(j["reduced_unknowns"] for j in SJ)), S_SY + ": reduced_unknowns max")
        put("SymUMax", num(max(j["unknowns"] for j in SJ)), S_SY + ": unknowns of the full system, max")
        fd = collections.Counter(j["field_D"] for j in SJ)
        put("NSymFieldQ", num(fd[1]), S_SY + ": field_D = 1")
        put("NSymFieldQuad", num(sum(v for D_, v in fd.items() if D_ != 1)), S_SY + ": field_D > 1")
        put("SymFieldList", ", ".join(f"$\\mathbb Q(\\sqrt{{{D_}}})$" for D_ in sorted(fd) if D_ != 1), S_SY)
        md = collections.Counter(j["model"] for j in SJ)
        put("NSymModelInv", num(md["involutions"]), S_SY + ": model involutions")
        put("NSymModelCyc", num(md["cyclic"]), S_SY + ": model cyclic")
        put("NSymModelLines", num(md["lines"]), S_SY + ": model lines")
        put("NSymModelIso", num(md["isotypic"]), S_SY + ": model isotypic")
        put("SymCPUMax", f"{max(j['cpu_seconds'] for j in SJ):.0f}", S_SY + ": max cpu_seconds")
        put("SymStressMinRel", sci_dn(min(j["stress_min_rel"] for j in SJ)), S_SY + ": min stress_min_rel")
        put("NSymUnmarked", num(sum(1 for k in symstrict if not idx[k]["flagged"])), S_SY + " (unmarked)")
        put("NSymDimThree", num(sum(1 for k in symstrict if k[0] == 3)), S_SY + " (d = 3)")
        kx = max(symstrict, key=lambda k: (sym[k]["unknowns"] / sym[k]["reduced_unknowns"], -k[0]))
        jx = sym[kx]
        put("SymExEntry", f"$({kx[0]},{kx[1]})$", S_SY + ": the entry with the largest ratio of unknowns to reduced unknowns")
        put("SymExGroup", num(jx["group_order_exact"]), S_SY + ": its group order")
        put("SymExU", num(jx["unknowns"]), S_SY + ": its full system")
        put("SymExRed", num(jx["reduced_unknowns"]), S_SY + ": its reduced system")
        put("SymExExtra", num(jx["extra"]), S_SY + ": its extra contacts")
    ctl = load_dir(os.path.join(HERE, "phase2", "out_controls"))
    if ctl:
        S_C = "phase2/out_controls/*.json (symstrict.py on entries certified strict by recognition)"
        cst = [k for k, j in ctl.items() if j["status"] == "STRICT"]
        put("NCtlRun", num(len(ctl)), S_C)
        put("NCtlStrict", num(len(cst)), S_C + ": STRICT")
        agree = 0
        for k in cst:
            cs = strict4[k].get("c_star")
            g = ctl[k]["gamma_mid"]
            if cs is not None and abs(Decimal(g) - Decimal(cs)) < Decimal(10) ** -35:
                agree += 1
        put("NCtlAgree", num(agree), S_C + ": gamma* agrees with strict.py's c* to 35 digits")
    for kind in ("model", "perm", "swap"):
        rd = load_dir(os.path.join(HERE, "phase2", f"out_red_{kind}"))
        if rd:
            S_R = f"phase2/out_red_{kind}/*.json (symstrict.py --red {kind})"
            put(f"NRed{kind.capitalize()}Run", num(len(rd)), S_R)
            put(f"NRed{kind.capitalize()}Strict", num(sum(1 for j in rd.values() if j["status"] == "STRICT")), S_R + ": STRICT (must be 0)")
    # ------------------------------------------------------------------------------------------- verdict table (M6)
    TV = ["% written by make_numbers_v6.py; do not edit", "\\begin{tabular}{p{0.16\\linewidth}p{0.35\\linewidth}p{0.17\\linewidth}p{0.17\\linewidth}}", "\\toprule",
          "verdict & configuration certified & distance to the stored code & value \\\\", "\\midrule",
          f"strict, isostatic & the zero $P^*$ of the square system, with exactly the proposed contacts; rattlers placed "
          f"(moved, by at most {M['RattlerMovedMax']}, in {M['NRattlerMoved']} cases) & $\\le$ {M['PassDistMax']} "
          f"({M['NPassDistAbove']} above $10^{{-8}}$) & $c^*$ in a ball of radius $\\le$ {M['QCErrMax']} \\\\",
          f"strict, exact Gram matrix & the configuration with the recognised Gram matrix & $\\le$ {M['StrictDistMax']} & "
          f"$c^*$ exact (algebraic) \\\\"]
    if symstrict:
        TV.append(f"strict, symmetry & the zero of the reduced system, all contacts exactly equal & $\\le$ {M['SymDistMax']} & "
                  f"$\\gamma^*$ in a ball of radius $\\le$ {M['SymGammaRadMax']} \\\\")
    TV += [f"strict, exact arithmetic & the stored or recognised exact configuration & $\\le$ {M['ExOnlyDistLSMax']} & exact \\\\",
           f"strict, exact core & the exact core, stored rattlers & $\\le$ {M['NewDistMax']} & exact \\\\",
           f"local minimum & a local minimum $Q^0$ within the tube around $P'$ & $\\le$ {M['LocDistTubeMax']} & "
           f"$M(Q^0)$ within {M['ChiMax']} of $c'$ \\\\",
           f"not strict & the exact core (two rigid halves) & $\\le$ {M['SplitDistMax']} & exact \\\\",
           f"not a local min., negative stress & the equal-contact zero $Q^*$ of the \\emph{{proposed}} graph & $\\le$ "
           f"{sci_up(max(nod))} & $c^*$ enclosed \\\\",
           f"not a local min., flat core & the exact core, stored rattlers & $\\le$ {M['FlatDistMax']} & exact \\\\",
           f"no local min.\\ within $\\delta$ & every configuration within $\\delta$ of the stored code (max.\\ per point) & "
           f"$\\delta$ from {M['FlexDeltaMin']} to {M['FlexDeltaMax']} & --- \\\\",
           "\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v6_tab_verdicts.tex"), "w").write("\n".join(TV) + "\n")
    # ------------------------------------------------------------------------------------------- copied from summary_v4.json
    FROM_V4 = ("NClassVaries", "SizeCap", "NWatchdogEntries", "NStoppedMemory", "NEntriesStoppedMemory",
               "CPUHoursCensus", "CPUHoursStrict")
    for a in FROM_V4:
        put(a, m4[a], "census/summary_v4.json macro " + a + " (census/summarize.py from the per-entry outputs)")
    put("NHypoEmptyTen", num(int(m4["NClassEmptyTen"].replace("{,}", "")) + int(m4["NClassHypoTen"].replace("{,}", ""))),
        "census/summary_v4.json NClassEmptyTen + NClassHypoTen")
    put("NIsoTen", m4["NClassIsoTen"], "census/summary_v4.json NClassIsoTen")
    put("NHyperTen", m4["NClassHyperTen"], "census/summary_v4.json NClassHyperTen")
    # ------------------------------------------------------------------------------------------- write
    with open(os.path.join(HERE, "v6_numbers.tex"), "w") as f:
        f.write("% written by make_numbers_v6.py; do not edit\n")
        for k_ in sorted(M):
            f.write(f"\\newcommand{{\\{k_}}}{{{M[k_]}}}\n")
    json.dump(M, open(os.path.join(HERE, "v6_numbers.json"), "w"), indent=1, sort_keys=True)
    with open(os.path.join(HERE, "VERIFY_numbers.md"), "w") as f:
        f.write("# Every macro of v6_numbers.tex: value and source\n\nWritten by make_numbers_v6.py; do not edit.\n\n")
        f.write("| macro | value | source |\n|---|---|---|\n")
        for k_ in sorted(M):
            v = str(M[k_]).replace("|", "\\|")
            f.write(f"| `\\{k_}` | {v} | {SRC[k_]} |\n")
    print(len(M), "macros; classes", dict(C))


if __name__ == "__main__":
    main()
