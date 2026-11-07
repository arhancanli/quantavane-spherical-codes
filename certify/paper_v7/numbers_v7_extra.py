"""The numbers, tables and figures that version 7 adds (called by make_numbers_v7.py; see its docstring).

Reads, from this folder: verify/out_exact/*.json (exactcheck.py), verify/out_strict7/*.json and verify/out_sym7/*.json
(strict7.py, symstrict7.py), recount/flexgram7_*.jsonl, recount/flexexact.jsonl, recount/gapcut_*.jsonl,
recount/gapcut_cert.jsonl, example/example_d3_n7.json, example/example_d5_n16.json, live/out/*.json (run_live.py),
live/live_select.tsv; and from the repository census/exact_v4/big/flex_*.json and v5/live/live_20261107.tsv.
Writes v7_tab_quant.tex, v7_tab_flexcat.tex, v7_tab_live.tex, v7_tab_notable.tex, fig_*_v7.pdf and v7_verdicts.tsv.
"""
import collections, glob, json, math, os, sys
from fractions import Fraction


def run(g):
    put, num, sci, sci_up, sci_dn, quant, elist = (g[k] for k in ("put", "num", "sci", "sci_up", "sci_dn", "quant", "elist"))
    HERE, CERT, CEN = g["HERE"], g["CERT"], g["CEN"]
    keys, idx, cls, res, first, passes, locs = (g[k] for k in ("keys", "idx", "cls", "res", "first", "passes", "locs"))
    strict4, sym = g["strict4"], g["sym"]
    np = g["np"]
    M = g["M"]

    def jl(pattern):
        out = []
        for p in sorted(glob.glob(os.path.join(HERE, pattern))):
            out += [json.loads(x) for x in open(p) if x.strip()]
        return out

    def jd(pattern):
        return [json.load(open(p)) for p in sorted(glob.glob(os.path.join(HERE, pattern)))]

    def f2(x, nd=2):
        return f"${x:.{nd}f}$"

    def dec(x, sig=2, how="up"):
        """A plain decimal with sig significant digits, rounded up or down (for values from 0.001 to 1000); else sci."""
        from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR
        if not (1e-3 <= abs(x) < 1e3):
            return (sci_up if how == "up" else sci_dn)(x, sig - 1)
        q = Decimal(repr(x))
        e = q.adjusted()
        step = Decimal(1).scaleb(e - sig + 1)
        r = q.quantize(step, rounding=ROUND_CEILING if how == "up" else ROUND_FLOOR)
        return "$" + format(r, "f") + "$"

    # ======================================================================== exact decisions (exactcheck.py)
    S_X = "verify/out_exact/*.json (exactcheck.py)"
    ex = {}
    for j in jd("verify/out_exact/*.json"):
        ex[(j["d"], j["n"], j["census_status"])] = j
    v1_fail = sorted(k for k, j in ex.items() if j["status"] != k[2])
    ex2 = {}
    for j in jd("verify/out_exact2/*.json"):            # version 2 (exactcheck2.py), run where version 1 failed
        ex2[(j["d"], j["n"], j["census_status"])] = j
    assert set(ex2) == set(v1_fail), ("version 2 must cover exactly the version 1 failures", set(ex2) ^ set(v1_fail))
    ex.update(ex2)
    S_X = S_X + " and verify/out_exact2/*.json (exactcheck2.py, where exactcheck.py left a stress sign undecided)"
    agree = collections.Counter((k[2], j["status"]) for k, j in ex.items())
    bad = {k: j["status"] for k, j in ex.items() if j["status"] != k[2]}
    put("XNChecked", num(len(ex)), S_X + ": files")
    put("XNPassOK", num(agree[("PASS", "PASS")]), S_X + ": census PASS, exact PASS")
    put("XNLocOK", num(agree[("LOCALMIN", "LOCALMIN")]), S_X + ": census LOCALMIN, exact LOCALMIN")
    put("XNNonoptOK", num(agree[("NONOPT", "NONOPT")]), S_X + ": census NONOPT, exact NONOPT")
    put("XNDisagree", num(len(bad)), S_X + ": exact status differs from the census status (" + repr(sorted(bad.items()))[:300] + ")")
    exp = {k: ex.get((k[0], k[1], "PASS")) for k in passes}
    exl = {k: ex.get((k[0], k[1], "LOCALMIN")) for k in locs}
    assert all(v is not None for v in exp.values()) and all(v is not None for v in exl.values()), "exact outputs missing"
    allx = list(ex.values())
    put("XMaxU", num(max(j["unknowns"] for j in allx)), S_X + ": max unknowns")
    cpu = [j["cpu_seconds"] for j in allx]
    put("XHours", f"{sum(cpu) / 3600:.1f}", S_X + ": sum of cpu_seconds / 3600")
    put("XSecMed", f"{quant(cpu)[2]:.1f}", S_X + ": cpu_seconds, median")
    jm = max(allx, key=lambda j: j["cpu_seconds"])
    put("XSecMax", f"{jm['cpu_seconds']:.0f}", S_X + ": cpu_seconds, max")
    put("XSecMaxEntry", f"$({jm['d']},{jm['n']})$", S_X + ": entry of the max cpu_seconds")
    put("XSecMaxU", num(jm["unknowns"]), S_X + ": its unknowns")
    big = [j for j in allx if j["unknowns"] >= 1000]
    b_, a_ = np.polyfit(np.log([j["unknowns"] for j in big]), np.log([j["cpu_seconds"] for j in big]), 1)
    put("XCPUSlope", f"{b_:.1f}", S_X + ": least-squares slope of log cpu_seconds on log unknowns, unknowns >= 1000")
    # comparison with the first pass
    same, ratios_s, ratios_k, ratios_d = 0, [], [], []
    for k in passes + locs:
        st = res[k]["status"]
        e_ = ex[(k[0], k[1], st)]
        c_ = first[k]
        if e_["c_interval_30"] == c_.get("c_interval_30"):
            same += 1
        else:
            lo1, hi1 = map(int, e_["c_interval_30"])
            lo2, hi2 = map(int, c_["c_interval_30"])
            assert lo1 <= hi2 and lo2 <= hi1, ("c intervals disjoint", k)
        if c_.get("stress_min_lower") and e_.get("stress_min_lower"):
            ratios_s.append(e_["stress_min_lower"] / c_["stress_min_lower"])
        if c_.get("kappa") and e_.get("kappa"):
            ratios_k.append(e_["kappa"] / c_["kappa"])
        if c_.get("dist_ls") and e_.get("dist_ls"):
            ratios_d.append(e_["dist_ls"] / c_["dist_ls"])
    put("XCSame", num(same), S_X + " vs census/results_v4: identical 30-digit intervals (all others overlap; asserted)")
    put("XStressRatioMin", f2(min(ratios_s), 3), S_X + ": min of exact/census stress_min_lower")
    put("XStressRatioMax", f2(max(ratios_s), 3), S_X + ": max of exact/census stress_min_lower")
    put("XKappaRatioMin", f2(min(ratios_k), 3), S_X + ": min of exact/census kappa")
    put("XKappaRatioMax", f2(max(ratios_k), 3), S_X + ": max of exact/census kappa")
    put("XDistRatioMax", f2(max(ratios_d), 3), S_X + ": max of exact/census dist_ls")
    put("XDistRatioMed", f2(quant(ratios_d)[2], 3), S_X + ": median of exact/census dist_ls")
    # constants (isostatic)
    P_ = [exp[k] for k in passes]
    put("XNRhoDefined", num(sum(1 for j in P_ if j.get("rho") is not None)), S_X + ": PASS with rho")
    put("XNoRhoIso", num(sum(1 for j in P_ if j.get("rho") is None)), S_X + ": PASS without rho")
    put("XInsideIso", num(sum(1 for j in P_ if j.get("code_inside_rho") is True)), S_X + ": PASS, code_inside_rho")
    put("XOutsideIso", num(sum(1 for j in P_ if j.get("code_inside_rho") is False)), S_X + ": PASS, not code_inside_rho")
    put("XKappaHalfMed", sci(quant([j["kappa_half"] for j in P_ if j.get("kappa_half")])[2]), S_X + ": PASS kappa_half median")
    put("XRhoMed", sci(quant([j["rho"] for j in P_ if j.get("rho")])[2]), S_X + ": PASS rho median")
    put("XRUniqueMed", sci(quant([j["r_unique"] for j in P_])[2]), S_X + ": PASS r_unique median")
    put("XUniqueCovers", num(sum(1 for j in P_ if j["r_unique"] > j["dist_coord"])), S_X + ": PASS, r_unique > dist_coord")
    wid = max((int(j["c_interval_30"][1]) - int(j["c_interval_30"][0])) for j in allx)
    put("XCErrMax", f"$10^{{-30}}$" if wid == 1 else f"${wid}\\cdot10^{{-30}}$", S_X + ": widest 30-digit interval")
    L_ = [exl[k] for k in locs]
    put("XInsideHyp", num(sum(1 for j in L_ if j.get("code_inside_rho") is True)), S_X + ": LOCALMIN, code_inside_rho")
    put("XChiMax", sci_up(max(j["chi"] for j in L_)), S_X + ": LOCALMIN max chi")
    put("XTubeMax", sci_up(max(j["tube"] for j in L_)), S_X + ": LOCALMIN max tube")
    # the gap to the nearest pair outside the contacts, rattlers included: the first pass's value (min_gap_float);
    # the exact decision proves every such gap positive, and records the exact lower bound over the non-rattlers
    S_G1 = "census/results_v4: gaps.min_gap_float of the deciding attempt (first pass, rattlers included)"
    gi = [(first[k]["gaps"]["min_gap_float"], k) for k in passes]
    gh = [(first[k]["gaps"]["min_gap_float"], k) for k in locs]
    put("XGapIsoNine", num(sum(1 for x, _ in gi if x < 1e-9)), S_G1 + ": PASS < 1e-9")
    put("XGapIsoEleven", num(sum(1 for x, _ in gi if x < 1e-11)), S_G1 + ": PASS < 1e-11")
    put("XGapHypNine", num(sum(1 for x, _ in gh if x < 1e-9)), S_G1 + ": LOCALMIN < 1e-9")
    put("XGapHypEleven", num(sum(1 for x, _ in gh if x < 1e-11)), S_G1 + ": LOCALMIN < 1e-11")
    gmin, kg = min(gi + gh)
    jg = ex[(kg[0], kg[1], res[kg]["status"])]
    put("XGapMin", sci_dn(gmin, 2), S_G1 + ": smallest")
    put("XGapMinEntry", f"$({kg[0]},{kg[1]})$", S_G1 + ": its entry")
    put("XGapMinTol", f"$10^{{{int(round(math.log10(jg['tol'])))}}}$", S_X + ": its tolerance")
    put("XGapMinExact", sci_dn(jg["gap_lower"], 2), S_X + ": its exact gap_lower (non-rattlers)")
    T = ["% written by make_numbers_v7.py; do not edit", "\\begin{tabular}{lrrrrr}", "\\toprule",
         " & min & 10\\% & median & 90\\% & max \\\\", "\\midrule"]
    for key, lab in (("kappa_half", "growth rate $\\kappa/(2\\sqrt d)$"), ("rho", "radius $\\rho$"),
                     ("dist_ls", "distance $D(X_V)$ to the stored code"), ("r_unique", "uniqueness radius $r_u$"),
                     ("beta_K", "largest row sum $\\beta_K$ in (K)"),
                     ("stress_min_rel", "smallest relative stress $|E|\\,\\omega_-$"),
                     ("unknowns", "unknowns $U$"), ("cpu_seconds", "CPU time of the exact decision (s)")):
        qq = quant([float(j[key]) for j in P_ if j.get(key) is not None])
        if key == "unknowns":
            cells = [num(int(x)) for x in qq]
        elif key == "cpu_seconds":
            cells = [f"${x:.2f}$" if x < 10 else f"${x:.0f}$" for x in qq]
        else:
            cells = [sci(x) for x in qq]
        T.append(lab + " & " + " & ".join(cells) + " \\\\")
    T += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v7_tab_quant.tex"), "w").write("\n".join(T) + "\n")
    # strictness re-checks
    S_S = "verify/out_strict7/*.json (strict7.py: strict.py with exact (S') and span)"
    st7 = {(j["d"], j["n"]): j for j in jd("verify/out_strict7/*.json")}
    st7b = {(j["d"], j["n"]): j for j in jd("verify/out_strict7b/*.json")}
    s7_fail = sorted(k for k, j in st7.items() if j.get("status") != "STRICT")
    assert set(st7b) == set(s7_fail), ("strict7b must cover exactly the strict7 failures", set(st7b) ^ set(s7_fail))
    st7.update(st7b)
    st7 = list(st7.values())
    sy7 = {(j["d"], j["n"]): j for j in jd("verify/out_sym7/*.json")}
    sy7b = {(j["d"], j["n"]): j for j in jd("verify/out_sym7b/*.json")}
    y7_fail = sorted(k for k, j in sy7.items() if j.get("status") != "STRICT")
    assert set(sy7b) == set(y7_fail), ("sym7b must cover exactly the sym7 failures", set(sy7b) ^ set(y7_fail))
    sy7.update(sy7b)
    sy7 = list(sy7.values())
    used2 = sorted(set((j["d"], j["n"]) for j in list(ex2.values()) if j.get("stress_enclosure") == "both")
                   | set((j["d"], j["n"]) for j in st7 + sy7 if j.get("exact_stress_enclosure") == "both"))
    put("NNoVerdict", num(int(M["NFailedAll"].replace("{,}", "")) + int(M["NTooLarge"].replace("{,}", ""))), "NFailedAll + NTooLarge")
    put("XStageTwoN", num(len(used2)), "certificates whose exact stress test needed the second enclosure (stress_enclosure 'both')")
    put("XStageTwoList", elist(used2) if used2 else "none", "their entries")
    rep_lp = sorted(k for k in s7_fail + y7_fail if (st7b.get(k) or sy7b.get(k) or {}).get("status") == "STRICT")
    put("XRepeatList", elist(rep_lp) if rep_lp else "none", "strictness re-checks repeated with version 2 after the untrusted LP hit its time limit (stress_lp_status 1)")
    put("XRepeatN", num(len(rep_lp)), "their number")
    put("XStrictRun", num(len(st7)), S_S + ": files")
    put("XStrictOK", num(sum(1 for j in st7 if j.get("status") == "STRICT")), S_S + ": STRICT")
    S_Y = "verify/out_sym7/*.json (symstrict7.py: symstrict.py with exact (S') and span)"
    put("XSymRun", num(len(sy7)), S_Y + ": files")
    put("XSymOK", num(sum(1 for j in sy7 if j.get("status") == "STRICT")), S_Y + ": STRICT")
    put("XStrictHours", f"{sum(j.get('cpu_seconds', 0) for j in st7 + sy7) / 3600:.1f}", S_S + " and " + S_Y + ": cpu_seconds / 3600")
    # verdict-table cells
    put("XPassDistMax", sci_up(max(j["dist_ls"] for j in P_)), S_X + ": PASS max dist_ls")
    put("XNPassDistAbove", num(sum(1 for j in P_ if j["dist_ls"] > 1e-8)), S_X + ": PASS dist_ls > 1e-8")
    lopen = [exl[k] for k in locs if cls[k][0] == "local"]
    put("XLocDistTubeMax", sci_up(max(j["dist_ls_plus_tube"] for j in lopen)), S_X + ": open LOCALMIN max dist_ls_plus_tube")
    put("XNonoptDistMax", sci_up(max(j["dist_ls"] for k, j in ex.items() if k[2] == "NONOPT")), S_X + ": NONOPT max dist_ls")

    # ======================================================================== negative certificate
    S_F = "recount/flexgram7_*.jsonl (flexgram7.py)"
    fg = {(r["d"], r["n"]): r for r in jl("recount/flexgram7_*.jsonl")}
    flex = [k for k in keys if cls[k][1] == "flex"]
    assert set(flex) <= set(fg), ("flexgram7 incomplete", len(set(flex) - set(fg)))
    F_ = [fg[k] for k in flex]
    assert all(r["proposal_matches_census"] for r in F_)
    dl = [r["delta"] for r in F_]
    put("FlexDeltaMin", sci_dn(min(dl)), S_F + ": min delta")
    put("FlexDeltaMax", sci_dn(max(dl)), S_F + ": max delta")
    put("FlexDeltaMed", sci(quant(dl)[2]), S_F + ": median delta")
    put("FlexEpsMax", sci_up(max(r["eps"] for r in F_)), S_F + ": max eps")
    put("FlexTauRatioMin", f"{min(r['tau_prime'] / r['tol'] for r in F_):.2f}", S_F + ": min tau'/tol")
    put("NFlexTauBelowTau", num(sum(1 for r in F_ if r["tau_prime"] < r["tol"])), S_F + ": tau' < tol")
    put("FlexHMed", sci(quant([r["h"] / math.sqrt(r["d"]) for r in F_])[2]), S_F + ": h/sqrt(d), median")
    put("NFlexDeltaAtTauQuarter", num(sum(1 for r in F_ if r["tau_prime"] / 4 <= r["h"] / math.sqrt(r["d"]))), S_F + ": tau'/4 <= h/sqrt d")
    put("NFlexDeltaAboveNine", num(sum(1 for x in dl if x >= 1e-9)), S_F + ": delta >= 1e-9")
    put("NFlexDeltaBelowNine", num(sum(1 for x in dl if x < 1e-9)), S_F + ": delta < 1e-9")
    put("NFlexDeltaBelowTwelve", num(sum(1 for x in dl if x < 1e-12)), S_F + ": delta < 1e-12")
    rmin = min(F_, key=lambda r: r["delta"])
    put("FlexDeltaMinEntry", f"$({rmin['d']},{rmin['n']})$", S_F + ": entry of the min delta")
    put("FlexDeltaMinH", sci_up(rmin["h"]), S_F + ": its h")
    put("NFlexFlagged", num(sum(1 for k in flex if idx[k]["flagged"])), "census/index.tsv: FLEX entries marked")

    def category(k):
        j = res[k]
        if any(a.get("planned_kind") in ("iso", "hyper") for a in j["attempts"]):
            return "after-failed"
        if any(s["kind"] in ("iso", "hyper") for s in j.get("skipped_above_cap", [])):
            return "above-cap"
        return "no-proposal-flagged" if idx[k]["flagged"] else "no-proposal-unflagged"
    cat = {k: category(k) for k in flex}
    labels = [("no-proposal-unflagged", "no proposal, unmarked"), ("no-proposal-flagged", "no proposal, marked"),
              ("after-failed", "attempt failed"), ("above-cap", "above the size limit"), (None, "all")]
    TF = ["% written by make_numbers_v7.py; do not edit", "\\begin{tabular}{lrrrrrr}", "\\toprule",
          "category & entries & min $\\delta$ & median & max & $\\delta=\\tau'/4$ & $\\delta<10^{-9}$ \\\\", "\\midrule"]
    curves = {}
    for c_, lab in labels:
        ks = [k for k in flex if c_ is None or cat[k] == c_]
        ds = [fg[k]["delta"] for k in ks]
        q = quant(ds)
        tq = sum(1 for k in ks if fg[k]["tau_prime"] / 4 <= fg[k]["h"] / math.sqrt(k[0]))
        if c_ is None:
            TF.append("\\midrule")
        else:
            curves[lab] = sorted(ds)
        TF.append(f"{lab} & {num(len(ks))} & {sci_dn(q[0])} & {sci(q[2])} & {sci_dn(q[4])} & "
                  f"{100 * tq / len(ks):.1f}\\% & {num(sum(1 for x in ds if x < 1e-9))} \\\\")
        if c_ is not None:
            put("FlexCat" + "".join(w.capitalize() for w in c_.split("-")) + "N", num(len(ks)), S_F + ": category " + c_)
            put("FlexCat" + "".join(w.capitalize() for w in c_.split("-")) + "Med", sci(q[2]), S_F + ": category " + c_ + ", median delta")
    TF += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v7_tab_flexcat.tex"), "w").write("\n".join(TF) + "\n")
    # exact re-check on the sample
    S_FX = "recount/flexexact.jsonl (flexexact.py)"
    fx = jl("recount/flexexact.jsonl")
    att = [r for r in fx if not str(r.get("status", "")).startswith("not attempted")]
    put("FXRun", num(len(fx)), S_FX + ": entries")
    put("FXOK", num(sum(1 for r in fx if r.get("status") == "FLEX")), S_FX + ": status FLEX")
    put("FXHashMatch", num(sum(1 for r in att if r.get("E_sha256") and fg.get((r["d"], r["n"]), {}).get("E_sha256") == r["E_sha256"])),
        S_FX + " vs " + S_F + ": same E_sha256")
    put("FXMaxRows", num(max(r["rows"] for r in att)), S_FX + ": max rows")
    put("FXMaxCols", num(max(r["columns"] for r in att)), S_FX + ": max columns")
    rr = [r["delta"] / fg[(r["d"], r["n"])]["delta"] for r in att if r.get("status") == "FLEX"]
    put("FXDeltaRatioMin", f2(min(rr), 3), S_FX + ": min of exact delta / flexgram7 delta")
    put("FXDeltaRatioMax", f2(max(rr), 3), S_FX + ": max of exact delta / flexgram7 delta")
    put("FXHours", f"{sum(r.get('cpu_seconds', 0) for r in fx) / 3600:.1f}", S_FX + ": cpu hours")
    # the second program's negative re-checks
    S_2 = "census/exact_v4/big/flex_*.json (exact_big.py --flex)"
    sp = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(CEN, "exact_v4", "big", "flex_*.json")))]
    put("SecondFlexN", num(len(sp)), S_2)
    put("SecondFlexGood", num(sum(1 for j in sp if j["status"] == "FLEX confirmed")), S_2 + ": FLEX confirmed")
    put("SecondFlexRowsMin", num(min(j["rows"] for j in sp)), S_2 + ": min rows")
    put("SecondFlexRowsMax", num(max(j["rows"] for j in sp)), S_2 + ": max rows")

    # ======================================================================== worked example
    S_E = "example/example_d3_n7.json (example.py 3 7)"
    E1 = json.load(open(os.path.join(HERE, "example", "example_d3_n7.json")))
    x1 = E1["exactcheck"]
    put("ExTableCos", "$" + E1["table_cos"] + "$", S_E + ": table cosine")
    put("ExNContacts", num(E1["contacts"]), S_E)
    dc = collections.Counter(E1["degrees"])
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven"}
    put("ExDegrees", " and ".join(f"{words[c]} points of degree {dg}" for dg, c in sorted(dc.items(), reverse=True)), S_E + ": degrees")
    put("ExU", num(E1["unknowns"]), S_E)
    put("ExNextGap", dec(E1["next_gap_double"], 3, "dn"), S_E + ": next pair below the largest cosine (double)")
    put("ExNormC", dec(E1["normC"], 3, "up"), S_E + ": ||C||_inf")
    put("ExR", sci(x1["r"]), S_E + ": r")
    put("ExBetaK", sci_up(x1["beta_K"]), S_E + ": beta_K")
    put("ExEpsMax", sci_up(E1["eps_max"]), S_E + ": max eps")
    put("ExEpsGamma", sci_up(E1["eps_gamma"]), S_E + ": eps_gamma")
    put("ExRu", dec(x1["r_unique"], 3, "dn"), S_E + ": r_unique")
    lo, hi = x1["c_interval_30"]
    pre = os.path.commonprefix([lo, hi])
    put("ExCDigits", "\\,".join(pre[i:i + 5] for i in range(0, len(pre), 5)), S_E + ": common digits of the 30-digit interval, in groups of five")
    put("ExStressMin", dec(x1["stress_min_lower"], 3, "dn"), S_E + ": smallest stress, lower bound")
    put("ExStressRel", dec(x1["stress_min_rel"], 3, "dn"), S_E + ": smallest stress times |E|")
    put("ExSpan", sci_up(x1["span_value"]), S_E + ": d r ||B^-1||")
    put("ExGap", dec(x1["gap_lower"], 3, "dn"), S_E + ": gap_lower")
    put("ExMuP", dec(x1["mu_p"], 3, "up"), S_E + ": mu_p")
    put("ExTheta", sci_up(x1["theta"]), S_E + ": vartheta")
    put("ExZeta", dec(x1["zeta"], 3, "dn"), S_E + ": zeta")
    put("ExKappa", dec(x1["kappa"], 3, "dn"), S_E + ": kappa")
    put("ExLambda", dec(x1["Lambda"], 3, "up"), S_E + ": Lambda")
    put("ExRho", dec(x1["rho"], 3, "dn"), S_E + ": rho")
    put("ExKappaHalf", dec(x1["kappa_half"], 3, "dn"), S_E + ": kappa/(2 sqrt d)")
    put("ExDist", sci_up(x1["dist_ls"]), S_E + ": dist_ls")
    put("ExCPU", f"{ex[(3, 7, 'PASS')]['cpu_seconds']:.1f}", S_X + ": (3,7) cpu_seconds")
    S_E2 = "example/example_d5_n16.json (example.py 5 16)"
    E2 = json.load(open(os.path.join(HERE, "example", "example_d5_n16.json")))
    x2 = E2["exactcheck"]
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

    # ======================================================================== the live table
    S_L = "v5/live/live_20261107.tsv (the live table's index read on 7 November 2026)"
    live = {}
    for line in open(os.path.join(CERT, "v5", "live", "live_20261107.tsv")):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        live[(int(p[0]), int(p[1]))] = p
    put("LiveDate", "7 November 2026", S_L)
    put("LiveListed", num(sum(1 for k in keys if k in live)), S_L + ": census entries still listed")
    groups = [("strict", "strict", lambda k: cls[k][0] == "strict"),
              ("local", "local, strictness open", lambda k: cls[k][0] == "local"),
              ("nonstrict", "local, not strict", lambda k: cls[k][0] == "nonstrict"),
              ("notlocal", "not local", lambda k: cls[k][0] == "notlocal"),
              ("flex", "none within $\\delta$", lambda k: cls[k][1] == "flex"),
              ("none", "no verdict", lambda k: cls[k][0] in ("none_tried", "none_size") and cls[k][1] != "flex")]
    # a lower bound of the certified value: the exact interval of the zero of the square system (c* = c' for the strict
    # hyperstatic entries; minus chi for the open ones); for the few exact-only verdicts the displayed value minus 1e-12
    # (the table rounds up)
    def cstar_lo(k):
        st = res[k]["status"]
        e_ = ex.get((k[0], k[1], st)) if st in ("PASS", "LOCALMIN") else None
        if e_:
            lo_ = Fraction(int(e_["c_interval_30"][0]), 10 ** 30)
            if cls[k][0] == "local":
                lo_ -= Fraction(e_["chi"])
            return lo_
        return Fraction(idx[k]["cos"]) - Fraction(1, 10 ** 12)
    TL = ["% written by make_numbers_v7.py; do not edit", "\\begin{tabular}{lrrrrrr}", "\\toprule",
          "verdict & entries & listed & lower now & below $c^*$ & median gain & not beaten \\\\", "\\midrule"]
    imp_of = {}
    for k in keys:
        if k in live:
            a, l = Fraction(idx[k]["cos"]), Fraction(live[k][2])
            imp_of[k] = float(a - l)
    stat = {}
    beaten_of = {}
    for tag, lab, f in groups:
        ks = [k for k in keys if f(k)]
        lk = [k for k in ks if k in live]
        im = [k for k in lk if imp_of[k] > 0]
        if tag in ("strict", "local", "nonstrict"):
            bt = [k for k in im if Fraction(live[k][2]) < cstar_lo(k)]
        else:
            bt = im
        for k in bt:
            beaten_of[k] = True
        med = quant([imp_of[k] for k in bt])[2] if bt else None
        stat[tag] = (ks, lk, im, med, bt)
        btc = num(len(bt)) if tag in ("strict", "local", "nonstrict") else "--"
        TL.append(f"{lab} & {num(len(ks))} & {num(len(lk))} & {num(len(im))} & {btc} & {sci(med) if med else '--'} & "
                  f"{num(len(lk) - len(bt))} \\\\")
    TL += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v7_tab_live.tex"), "w").write("\n".join(TL) + "\n")
    ks, lk, im, med, bt = stat["flex"]
    put("LiveFlexListed", num(len(lk)), S_L + ": flex entries listed")
    put("LiveFlexImproved", num(len(im)), S_L + ": flex entries listed with a lower cosine")
    put("LiveFlexPct", f"{100 * len(im) / len(lk):.1f}", S_L + ": flex entries listed with a lower cosine, percent")
    put("LiveFlexImpMed", sci(med), S_L + ": flex entries, median improvement")
    ks, lk, im, med, bt = stat["strict"]
    put("LiveStrictListed", num(len(lk)), S_L + ": strict entries listed")
    put("LiveStrictLower", num(len(im)), S_L + ": strict entries listed with a lower displayed cosine")
    put("LiveStrictBeaten", num(len(bt)), S_L + ": strict entries listed with a cosine below the certified value")
    put("LiveStrictPct", f"{100 * len(bt) / len(lk):.0f}", S_L + ": strict entries beaten, percent of listed")
    put("LiveStrictImpMed", sci(med), S_L + ": beaten strict entries, median improvement of the displayed cosine")
    put("LiveStrictAuthor", num(sum(1 for k in bt if "Arhan Canli" in live[k][4])), S_L + ": beaten strict entries credited to the author")
    put("LiveStrictRefined", num(len(im) - len(bt)), S_L + ": strict entries with a lower displayed cosine not below the certified value")
    cert = [k for k in keys if cls[k][0] in ("strict", "local", "nonstrict")]
    cl_ = [k for k in cert if k in live]
    ci = [k for k in cl_ if beaten_of.get(k)]
    put("LiveCertListed", num(len(cl_)), S_L + ": certified local minima still listed")
    put("LiveCertStillBest", num(len(cl_) - len(ci)), S_L + ": of these, not beaten (no listed cosine below the certified value)")
    put("LiveCertNotBestPct", f"{100 * len(ci) / len(cl_):.0f}", S_L + ": of these, beaten, percent")
    put("LiveCertBeaten", num(len(ci)), S_L + ": of these, beaten")
    flat = [k for k in keys if cls[k][1] == "notlocal_flat"]
    fi = [k for k in flat if k in live and imp_of[k] > 0]
    put("NFlatLiveImproved", num(len(fi)), S_L + ": flat-core entries with a lower cosine now")
    put("FlatLiveImpRange", (f"{sci(min(imp_of[k] for k in fi))} to {sci(max(imp_of[k] for k in fi))}" if fi else "--"),
        S_L + ": their improvements")
    # rerun of the rule on fetched codes
    S_LR = "live/out/*.json (run_live.py on codes fetched by fetch_live.py)"
    lr = {(j["d"], j["n"]): j for j in jd("live/out/*.json")}
    pos = [j for j in lr.values() if j["census_class"] != "flex"]
    fl = [j for j in lr.values() if j["census_class"] == "flex"]

    def verdict(j):
        if j.get("status_exact") in ("PASS", "LOCALMIN"):
            return j["status_exact"]
        if j.get("status_first_pass") in ("PASS", "LOCALMIN"):
            return "unconfirmed"
        return j.get("status_first_pass")
    put("LiveFetchPos", num(len(pos)), S_LR + ": positive-class entries")
    put("LiveFetchFlex", num(len(fl)), S_LR + ": sampled flex entries")
    put("LiveRerunPosPos", num(sum(1 for j in pos if verdict(j) in ("PASS", "LOCALMIN"))), S_LR + ": positive-class, exact PASS/LOCALMIN")
    put("LiveRerunPosStrictIso", num(sum(1 for j in pos if verdict(j) == "PASS")), S_LR + ": positive-class, exact PASS")
    put("LiveRerunPosFlex", num(sum(1 for j in pos if verdict(j) == "FLEX")), S_LR + ": positive-class, FLEX")
    put("LiveRerunPosOther", num(sum(1 for j in pos if verdict(j) not in ("PASS", "LOCALMIN", "FLEX"))), S_LR + ": positive-class, other")
    put("LiveRerunFlexPos", num(sum(1 for j in fl if verdict(j) in ("PASS", "LOCALMIN"))), S_LR + ": flex sample, exact PASS/LOCALMIN")
    put("LiveRerunFlexFlex", num(sum(1 for j in fl if verdict(j) == "FLEX")), S_LR + ": flex sample, FLEX")
    put("LiveRerunFlexNone", num(sum(1 for j in fl if verdict(j) not in ("PASS", "LOCALMIN", "FLEX"))), S_LR + ": flex sample, no verdict")
    put("LiveRerunFlexPosPct", f"{100 * sum(1 for j in fl if verdict(j) in ('PASS', 'LOCALMIN')) / max(1, len(fl)):.0f}", S_LR + ": flex sample, positive, percent")
    put("LiveRerunUnconfirmed", num(sum(1 for j in lr.values() if verdict(j) == "unconfirmed")), S_LR + ": first pass positive, exact not")
    dlow = [j["dist_lower_point"] for j in pos if j.get("dist_lower_point") is not None]
    put("LiveDistLowMed", sci_dn(quant(dlow)[2]) if dlow else "--", S_LR + ": positive-class, median dist_lower_point")
    put("LiveDistLowMin", sci_dn(min(dlow)) if dlow else "--", S_LR + ": positive-class, min dist_lower_point")
    pb = [j for j in pos if j["census_class"] in ("strict", "local", "nonstrict") and beaten_of.get((j["d"], j["n"]))]
    db = [j["dist_lower_point"] for j in pb if j.get("dist_lower_point") is not None]
    put("LiveRerunBeatenN", num(len(pb)), S_LR + ": fetched positive-class entries beaten (live cosine below the certified value)")
    put("LiveDistLowMedBeaten", sci_dn(quant(db)[2]) if db else "--", S_LR + ": beaten ones, median dist_lower_point")
    put("LiveDistLowMinBeaten", sci_dn(min(db)) if db else "--", S_LR + ": beaten ones, min dist_lower_point")
    put("LiveRerunBeatenPos", num(sum(1 for j in pb if verdict(j) in ("PASS", "LOCALMIN"))), S_LR + ": beaten ones whose new code gets an exact PASS/LOCALMIN")
    put("LiveRerunBeatenFlex", num(sum(1 for j in pb if verdict(j) == "FLEX")), S_LR + ": beaten ones whose new code gets FLEX")

    # ======================================================================== gap-based proposals
    S_G = "recount/gapcut_*.jsonl (gapcut.py scan)"
    gs = jl("recount/gapcut_[0-9]*.jsonl")
    put("GapScanN", num(len(gs)), S_G + ": entries scanned")
    put("GapFound", num(sum(1 for r in gs if r.get("found"))), S_G + ": a non-hypostatic cut found")
    put("GapFoundSmall", num(sum(1 for r in gs if r.get("found") and r["spread"] <= 1e-7)), S_G + ": with spread <= 1e-7")
    S_GC = "recount/gapcut_cert.jsonl (gapcut.py certify)"
    gc = jl("recount/gapcut_cert.jsonl")
    gpos = [r for r in gc if (r.get("exact") or {}).get("status") in ("PASS", "LOCALMIN")]
    put("GapCertRun", num(len(gc)), S_GC + ": entries")
    put("GapCertPos", num(len(gpos)), S_GC + ": exact PASS or LOCALMIN")
    put("GapCertList", elist(sorted((r["d"], r["n"]) for r in gpos)) if gpos else "none", S_GC)
    put("GapCertDistMax", sci_up(max(r["exact"]["dist_ls"] for r in gpos)) if gpos else "--", S_GC + ": max dist_ls")
    put("GapCertExistFail", num(sum(1 for r in gc if str(r.get("first_pass", "")).startswith("FAIL existence"))), S_GC + ": existence test failed")
    put("GapCertMaxU", num(1500), "gapsample.py: at most 1,500 unknowns")

    # ======================================================================== notable entries
    def vname(k):
        c, kd = cls[k]
        return {"strict_iso": "strict, isostatic", "strict_hyper": "strict, recognition", "strict_sym": "strict, symmetry",
                "strict_exact": "strict, exact", "strict_new": "strict, exact core", "local": "local, strictness open",
                "nonstrict_split": "local, not strict", "notlocal_flat": "not local, flat core",
                "notlocal_nonopt": "not local, negative stress", "flex": "no local min.\\ within $\\delta$",
                "failed": "no verdict", "none_size": "no verdict (size)"}.get(kd, kd)
    opt = [k for k in keys if idx[k]["optimal"]]
    cjkt = [k for k in [(4, 24), (5, 40), (6, 72), (7, 126), (8, 240), (9, 306)] if k in cls]
    assert all(cls[k][0] == "strict" for k in cjkt)
    put("CJKTList", elist(cjkt), "kissing configurations of [CJKT11] that are census entries (all strict)")
    rows_ = []

    def grp(title, ks, verdict_txt=None):
        byv = collections.defaultdict(list)
        for k in ks:
            byv[verdict_txt or vname(k)].append(k)
        for v, kk in byv.items():
            rows_.append(f"{title} & {elist(sorted(kk))} & {v} \\\\")
    grp("provably optimal", opt)
    grp("kissing configurations~\\cite{CJKT11}", cjkt)
    grp("worked examples (Sec.~\\ref{sec:example})", [(3, 7), (5, 16)])
    grp("strict by symmetry", [k for k in keys if cls[k][1] == "strict_sym"])
    grp("strict by exact arithmetic", [k for k in keys if cls[k][1] in ("strict_exact", "strict_new")])
    grp("local minima, not strict", [k for k in keys if cls[k][1] == "nonstrict_split"])
    grp("not local minima", [k for k in keys if cls[k][0] == "notlocal"])
    grp("not decided (Sec.~\\ref{sec:open})", [(22, 116), (19, 166), (22, 47), (22, 48), (23, 50)])
    ext = [(rmin["d"], rmin["n"])]
    rmax = max(F_, key=lambda r: r["delta"])
    rows_.append(f"largest and smallest $\\delta$ & {elist([(rmax['d'], rmax['n'])])}, {elist(ext)} & no local min.\\ within $\\delta$ \\\\")
    rows_.append(f"smallest proven gap & {M['XGapMinEntry']} & {vname((jg['d'], jg['n']))} \\\\")
    TN = ["% written by make_numbers_v7.py; do not edit", "\\begin{tabular}{p{0.2\\linewidth}p{0.55\\linewidth}p{0.19\\linewidth}}",
          "\\toprule", "group & entries $(d,N)$ & verdict \\\\", "\\midrule"] + rows_ + ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(HERE, "v7_tab_notable.tex"), "w").write("\n".join(TN) + "\n")

    # ======================================================================== per-entry file
    with open(os.path.join(HERE, "v7_verdicts.tsv"), "w") as f:
        f.write("# Census version 7: one line per entry of the archived table with 3<=d<=32, N<=1024 (make_numbers_v7.py).\n"
                "# verdict: strict | local | nonstrict | notlocal | flex (no local minimum within delta) | none\n"
                "# c_lo_30, c_hi_30: the 30-digit decimal interval of the certified cosine (zero of the square system); "
                "dist_ls: certified least-squares distance to the stored code (non-rattlers; local minima: plus the tube); "
                "gap_lower: exact lower bound of the gap to the nearest non-contact pair of non-rattlers; gap_first_pass: the "
                "first pass's value of the gap to the nearest non-contact pair, rattlers included (double precision); "
                "tau_prime, h, delta: the negative certificate.\n")
        f.write("d\tn\tflagged\toptimal\tverdict\tcertificate\ttol\tc_lo_30\tc_hi_30\tdist_ls\tgap_lower\tgap_first_pass\ttau_prime\th\t"
                "delta\tlive_cos\tlive_lower\n")
        for k in keys:
            c, kd = cls[k]
            v = {"strict": "strict", "local": "local", "nonstrict": "nonstrict", "notlocal": "notlocal"}.get(c, "flex" if kd == "flex" else "none")
            st = res[k]["status"]
            e_ = ex.get((k[0], k[1], st)) if st in ("PASS", "LOCALMIN") else ex.get((k[0], k[1], "NONOPT"))
            tol = e_["tol"] if e_ else (fg[k]["tol"] if k in fg and kd == "flex" else "")
            clo, chi_ = (e_["c_interval_30"] if e_ else ("", ""))
            dist = e_.get("dist_ls_plus_tube", e_.get("dist_ls")) if e_ else ""
            gap = e_.get("gap_lower", "") if e_ else ""
            g1 = first[k]["gaps"]["min_gap_float"] if k in first else ""
            tp, hh, dd = ((fg[k]["tau_prime"], fg[k]["h"], fg[k]["delta"]) if kd == "flex" else ("", "", ""))
            lc = live[k][2] if k in live else ""
            ll = (1 if imp_of.get(k, 0) > 0 else 0) if k in live else ""
            f.write(f"{k[0]}\t{k[1]}\t{idx[k]['flagged']}\t{idx[k]['optimal']}\t{v}\t{kd}\t{tol}\t{clo}\t{chi_}\t{dist}\t{gap}\t{g1}\t"
                    f"{tp}\t{hh}\t{dd}\t{lc}\t{ll}\n")

    # ======================================================================== figures
    sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/mpllib"))
    import matplotlib
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
                         "xtick.color": "#52514e", "ytick.color": "#52514e", "axes.linewidth": 0.6,
                         "pdf.fonttype": 42, "font.family": "serif"})
    SER = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
    GREY = "#b5b4ae"
    fig, ax = plt.subplots(figsize=(6.3, 3.4))
    order = [("none", "no verdict", GREY, ".", 2.0), ("flex", "no local minimum within $\\delta$", SER[4], ".", 2.0),
             ("notlocal", "not a local minimum", SER[3], "D", 9), ("nonstrict", "local minimum, not strict", SER[2], "s", 12),
             ("local", "local minimum, strictness open", SER[1], "x", 5), ("strict", "strict local minimum", SER[0], "o", 4)]
    for tag, lab, col, mk, sz in order:
        ks = stat[tag][0]
        ax.scatter([k[1] for k in ks], [k[0] for k in ks], s=sz, c=col, marker=mk, linewidths=0.6 if mk == "x" else 0,
                   label=f"{lab} ({len(ks):,})", rasterized=(len(ks) > 2000))
    ax.set_xscale("log")
    ax.set_xlabel("number of points $N$")
    ax.set_ylabel("dimension $d$")
    ax.set_ylim(2, 33)
    ax.grid(True, color="#e6e5e0", linewidth=0.4)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, markerscale=2.5, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_verdicts_v7.pdf"), dpi=300)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.0, 2.8))
    styles = ["-", "--", "-.", ":"]
    for i, (lab, ds) in enumerate(curves.items()):
        y = np.arange(1, len(ds) + 1) / len(ds)
        ax.step(ds, y, where="post", color=SER[i], linestyle=styles[i], linewidth=1.4, label=f"{lab} ({len(ds):,})")
    ax.set_xscale("log")
    ax.set_xlabel("radius $\\delta$ of the negative certificate")
    ax.set_ylabel("fraction of entries")
    ax.grid(True, color="#e6e5e0", linewidth=0.4)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_delta_v7.pdf"))
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.6))
    for ks, col, mk, lab in ((passes, SER[0], "o", "isostatic"), (locs, SER[1], "x", "hyperstatic")):
        U = [first[k]["unknowns"] for k in ks]
        e_ = exp if ks is passes else exl
        axs[0].scatter(U, [max(e_[k]["dist_ls"], 1e-19) for k in ks], s=4 if mk == "o" else 6, c=col, marker=mk,
                       linewidths=0 if mk == "o" else 0.6, label=lab)
        axs[1].scatter(U, [first[k]["cond_inf_upper"] for k in ks], s=4 if mk == "o" else 6, c=col, marker=mk,
                       linewidths=0 if mk == "o" else 0.6, label=lab)
    for a_, yl in ((axs[0], "certified distance $D(X_V)$"), (axs[1], "condition bound (first pass)")):
        a_.set_xscale("log")
        a_.set_yscale("log")
        a_.set_xlabel("unknowns $U$")
        a_.set_ylabel(yl)
        a_.grid(True, color="#e6e5e0", linewidth=0.4)
        a_.set_axisbelow(True)
    axs[0].set_ylim(5e-20, 1e-5)
    axs[0].legend(frameon=False, fontsize=7, markerscale=2, loc="lower right")
    put("XDistTinyN", num(sum(1 for k in passes if exp[k]["dist_ls"] < 1e-19) + sum(1 for k in locs if exl[k]["dist_ls"] < 1e-19)),
        S_X + ": certificates with dist_ls < 1e-19 (drawn on the lower edge of Figure 3)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_cost_v7.pdf"))
    plt.close(fig)
