"""Every number and table of the paper (version 5) that is new in version 5, from the census outputs (version 4) and the
version 5 outputs (certify/v5/out/).

  python3 -I certify/v5/summarize_v5.py PAPERDIR

Writes PAPERDIR/v5_numbers.tex, v5_tab_census.tex, v5_tab_census_dims.tex, v5_tab_open.tex and v5_tab_quant.tex, and
certify/v5/out/summary_v5.json (the class of every entry and every macro). The census numbers of version 4 that are
unchanged come from v4_numbers.tex (census/summarize.py).
"""
import collections, glob, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CEN = os.path.join(os.path.dirname(HERE), "census")
sys.path.insert(0, os.path.expanduser("~/Library/Caches/quantavane-census/pylib"))
sys.path.insert(0, CEN)
import numpy as np
import snapshot
from summarize import num, sci, quant, word

OUT = os.path.join(HERE, "out")
RANGES = [(3, 5), (6, 10), (11, 17), (18, 24), (25, 32)]
OPEN22 = [(8, 72), (9, 330), (10, 238), (12, 210), (16, 117), (17, 162), (19, 166), (20, 92), (22, 47), (22, 48),
          (22, 79), (22, 114), (22, 116), (22, 127), (23, 50), (23, 102), (23, 129), (23, 147), (24, 121), (24, 149),
          (26, 57), (27, 96)]
COLS = ["strict", "local", "nonstrict", "notlocal", "none_tried", "none_size"]


def entry_list(lst):
    return ", ".join(f"$({d},{n})$" for d, n in lst)


def main():
    pdir = sys.argv[1]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index() if e["listed"]}
    keys = sorted(k for k in idx if 3 <= k[0] <= 32 and k[1] <= 1024)
    res, strict, exact4 = {}, {}, {}
    for f in glob.glob(os.path.join(CEN, "results_v4", "*.json")):
        j = json.load(open(f))
        res[(j["d"], j["n"])] = j
    for f in glob.glob(os.path.join(CEN, "strict_v4", "*.json")):
        j = json.load(open(f))
        strict[(j["d"], j["n"])] = j
    for f in glob.glob(os.path.join(CEN, "exact_v4", "strict", "*.json")):
        j = json.load(open(f))
        exact4[(j["d"], j["n"])] = j
    # ------------------------------------------------------------------ version 4 statuses (as census/summarize.py)
    v4, first = {}, {}
    for k in keys:
        j = res[k]
        st = j["status"]
        att = {a.get("status"): a for a in j["attempts"]}
        if st == "PASS":
            v4[k], first[k] = "strict_iso", att["PASS"]
        elif st == "LOCALMIN":
            first[k] = att["LOCALMIN"]
            v4[k] = "strict_hyper" if strict.get(k, {}).get("status") == "STRICT" else "local"
        elif st == "FLEX":
            v4[k], first[k] = "flex", att["FLEX"]
        elif st == "not attempted":
            v4[k] = "too_large"
        else:
            v4[k] = "failed"
    nonopt = sorted(k for k in keys if any(a.get("status") == "NONOPT" for a in res[k]["attempts"]))
    # ------------------------------------------------------------------ version 5 outputs
    sub = {}
    for f in glob.glob(os.path.join(OUT, "subspace_d*_n*.json")):
        j = json.load(open(f))
        sub[(j["d"], j["n"])] = j
    rig = {}
    for f in glob.glob(os.path.join(OUT, "rigidity_d*_n*.json")):
        j = json.load(open(f))
        rig[(j["d"], j["n"])] = j
    jam = json.load(open(os.path.join(OUT, "jam_d22_n116.json")))
    dist = {(e["d"], e["n"]): e for e in json.load(open(os.path.join(OUT, "distances_exact_only.json")))["entries"]}
    live = {}
    for line in open(os.path.join(HERE, "live", "live_20261107.tsv")):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        live[(int(p[0]), int(p[1]))] = (p[2], p[4])
    # ------------------------------------------------------------------ version 5 classes
    cls, kind = {}, {}
    for k in keys:
        s = v4[k]
        if s in ("strict_iso", "strict_hyper"):
            cls[k], kind[k] = "strict", s
        elif s == "local":
            cls[k], kind[k] = "local", "local"
        elif s == "too_large":
            if exact4.get(k, {}).get("status", "").startswith("STRICT"):
                cls[k], kind[k] = "strict", "strict_exact"
            else:
                cls[k], kind[k] = "none_size", "none_size"
        else:
            cls[k], kind[k] = "none_tried", s
    for k in nonopt:
        cls[k], kind[k] = "notlocal", "notlocal_nonopt"
    for k, j in sub.items():
        if j.get("status", "").startswith("NOT A LOCAL MINIMUM"):
            cls[k], kind[k] = "notlocal", "notlocal_flat"
    for k, j in rig.items():
        st = j.get("status", "")
        if st.startswith("STRICT"):
            cls[k], kind[k] = "strict", "strict_new"
        elif st.startswith("LOCAL, NOT STRICT"):
            cls[k], kind[k] = "nonstrict", "nonstrict_split"
    if jam.get("status", "").startswith("JAMMED"):
        cls[(22, 116)], kind[(22, 116)] = "strict", "strict_second"
    cnt = collections.Counter(cls.values())
    kc = collections.Counter(kind.values())
    M = {}
    M["VNStrictAll"] = num(cnt["strict"])
    M["VNStrictIso"] = num(kc["strict_iso"])
    M["VNStrictHyper"] = num(kc["strict_hyper"])
    M["VNStrictExact"] = num(kc["strict_exact"])
    M["VNStrictNew"] = num(kc["strict_new"])
    M["VNStrictSecond"] = num(kc["strict_second"])
    M["VNStrictOther"] = num(kc["strict_exact"] + kc["strict_new"] + kc["strict_second"])
    M["VNLocal"] = num(cnt["local"])
    M["VNNonStrict"] = num(cnt["nonstrict"])
    M["VNNotLocal"] = num(cnt["notlocal"])
    M["VNNotLocalFlat"] = num(kc["notlocal_flat"])
    M["VNNotLocalNonopt"] = num(kc["notlocal_nonopt"])
    M["VNNoneTried"] = num(cnt["none_tried"])
    M["VNTooLarge"] = num(cnt["none_size"])
    M["VNNoneAll"] = num(cnt["none_tried"] + cnt["none_size"])
    M["VNCertAll"] = num(cnt["strict"] + cnt["local"] + cnt["nonstrict"] + cnt["notlocal"])
    M["VPctDecided"] = f"{100 * (cnt['strict'] + cnt['local'] + cnt['nonstrict'] + cnt['notlocal']) / len(keys):.1f}"
    M["VNUnflagCert"] = num(sum(1 for k in keys if not idx[k]["flagged"] and cls[k] in ("strict", "local")))
    M["VNFlagCert"] = num(sum(1 for k in keys if idx[k]["flagged"] and cls[k] in ("strict", "local")))
    M["VNStrictExactList"] = entry_list(sorted(k for k in keys if kind[k] == "strict_exact"))
    M["VNStrictNewList"] = entry_list(sorted(k for k in keys if kind[k] == "strict_new"))
    M["VNonoptList"] = entry_list(nonopt)
    # the exact-arithmetic sentence of version 4, with the counting of version 5
    import re as _re
    v4t = open(os.path.join(pdir, "v4_text.tex")).read()
    m_ = _re.search(r"\\newcommand\{\\ExactStrictText\}\{(.*)\}\s*$", v4t, flags=_re.M)
    est = m_.group(1)
    cut = est.find(". Of the confirmed entries")
    if cut > 0:
        above_ = sorted(k for k in keys if kind[k] == "strict_exact")
        est = est[:cut] + (". Of the confirmed entries, " + entry_list(above_) + " are above the size limit of the "
                           "census rule; for them this is the only certificate, and Table~\\ref{tab:census} counts them "
                           "as strict.")
    M["VExactStrictText"] = est
    EXPT = {1e-13: "-13", 1e-12: "-12", 1e-11: "-11", 1e-10: "-10", 1e-9: "-9", 1e-8: "-8", 1e-7: "-7"}
    rows_ = []
    for k in nonopt:
        atts = [a for a in res[k]["attempts"] if a.get("status") == "NONOPT"]
        taus = ", ".join(f"10^{{{EXPT[a['tol']]}}}" for a in atts)
        a = atts[0]
        q = a.get("stress_quantiles_over_mean") or [None]
        far = "; a configuration rather far from the code" if a["dist_ls"] > 1e-6 else ""
        rows_.append(f"$({k[0]},{k[1]})$ at $\\tau={taus}$ ($Q^*$ within {sci(a['dist_ls'])} of the stored code{far}; "
                     f"smallest relative stress {sci(q[0])})")
    M["VNonoptText"] = ("Proposition~\\ref{prop:nonopt} applies to " + " and ".join(rows_) + ".") if rows_ else \
        "Proposition~\\ref{prop:nonopt} applies to no entry."
    # the diagnostic, reported apart: the version 4 flex statements of entries without a certificate
    diag = [k for k in keys if kind[k] == "flex"]
    F = {k: first[k] for k in diag}
    M["VNDiag"] = num(len(diag))
    M["VNDiagForest"] = num(sum(1 for k in diag if F[k].get("forest")))
    rr = [F[k]["flex_radius_l1_stored"] for k in diag]
    M["VNDiagSmallRad"] = num(sum(1 for x in rr if x < 1e-6))
    M["VNDiagFlagged"] = num(sum(1 for k in diag if idx[k]["flagged"]))
    M["VNDiagUseful"] = num(sum(1 for k in diag if not F[k].get("forest") and F[k]["flex_radius_l1_stored"] >= 1e-6))
    q = quant(rr)
    M["VDiagRadMin"], M["VDiagRadMed"], M["VDiagRadMax"] = sci(q[0]), sci(q[2]), sci(q[4])
    M["VNFailed"] = num(sum(1 for k in keys if kind[k] == "failed"))
    # ------------------------------------------------------------------ the 22 open entries with exact cosines
    rows, verdict = [], {}
    for k in OPEN22:
        d, n = k
        kd = kind[k]
        if kd == "notlocal_flat":
            j = sub[k]
            v = "not a local minimum: flat core"
            info = (j["core_size"], j["rattlers"], j["core_rank"], j["contacts"], j["dist_ls_upper"])
        elif kd == "nonstrict_split":
            j = rig[k]
            v = "local minimum, not strict: split contact graph"
            info = (j["core_size"], j["rattlers"], j["core_rank"], j["contacts"], j["dist_ls_upper"])
        elif kd == "strict_new":
            j = rig[k]
            v = "strict local minimum (bar-rigid core)"
            info = (j["core_size"], j["rattlers"], j["core_rank"], j["contacts"], j["dist_ls_upper"])
        elif kd == "strict_second":
            v = "strict local minimum, jammed, not infinitesimally jammed"
            info = (116, 0, 22, jam["contacts"], jam["dist_ls_upper"])
        else:
            j = rig.get(k)
            v = "undecided"
            info = (None, None, None, None, None)
            if "continuous family" in idx[k]["notes"]:
                v = "undecided (continuous family)"
        verdict[k] = v
        poly = idx[k]["poly"].replace("*", "")
        lv = live.get(k)
        lstat = ("same" if lv and lv[0] == idx[k]["cos"] else ("lower" if lv else "not listed"))
        rows.append((k, poly, info, v, lstat))
    M["VNOpenExact"] = num(len(OPEN22))
    dec = [k for k in OPEN22 if not verdict[k].startswith("undecided")]
    M["VNDecidedExact"] = num(len(dec))
    flat = [k for k in OPEN22 if kind[k] == "notlocal_flat"]
    M["VNFlat"] = num(len(flat))
    flat_live = [k for k in flat if live.get(k) and live[k][0] == idx[k]["cos"]]
    M["VNFlatLive"] = num(len(flat_live))
    M["VFlatLiveList"] = entry_list(flat_live)
    M["VFlatList"] = entry_list(flat)
    M["VFlatDistMax"] = sci(max(sub[k]["dist_ls_upper"] for k in flat))
    M["VFlatColoursMin"] = num(min(sub[k]["colours"] for k in flat))
    M["VFlatColoursMax"] = num(max(sub[k]["colours"] for k in flat))
    M["VFlatRattlerGapMin"] = sci(min(-sub[k]["rattler_max_minus_c_upper"] for k in flat))
    split = [k for k in OPEN22 if kind[k] == "nonstrict_split"]
    M["VNUnion"] = num(sum(1 for k in flat if "orthogonal" in idx[k]["notes"].lower()))
    M["VNSplit"] = num(len(split))
    M["VSplitList"] = entry_list(split)
    M["VSplitDistMax"] = sci(max(rig[k]["dist_ls_upper"] for k in split)) if split else "--"
    newstrict = [k for k in OPEN22 if kind[k] == "strict_new"]
    fam = [k for k in OPEN22 if verdict[k] == "undecided (continuous family)"]
    und = [k for k in OPEN22 if verdict[k] == "undecided"]
    M["VNFamily"] = num(len(fam))
    M["VFamilyList"] = entry_list(fam)
    M["VUndecidedList"] = entry_list(und)
    M["VNUndecidedOther"] = num(len(und))
    # sentences for the abstract and the introduction
    ns = len(newstrict) + 1
    if newstrict:
        M["VNStrictNewSentence"] = ("\\emph{" + word_n(len(newstrict)).capitalize() + " more " +
                                    ("is a strict local minimum" if len(newstrict) == 1 else "are strict local minima") +
                                    "}, " + entry_list(newstrict) + ", by Theorem~\\ref{thm:main} checked in exact "
                                    "arithmetic, where the census subsystem had failed;")
    else:
        M["VNStrictNewSentence"] = ""
    M["VNFamilyItem"] = ("\\emph{" + word_n(len(fam) + len(und)).capitalize() + " remain undecided}: " + entry_list(fam) +
                         ", which the table notes as members of continuous families of codes without rattlers, and our "
                         "numerical analysis agrees (one infinitesimal flex, on which $H_\\omega$ vanishes)" +
                         (("; and " + entry_list(und)) if und else "") + ".")
    M["LiveDate"] = "7 November 2026"
    sp = [k for k in newstrict if rig[k].get("stressed_contacts")]
    M["VStressedSentence"] = " ".join(
        f"For $({k[0]},{k[1]})$ no stress is positive on every contact, but the {num(rig[k]['stressed_contacts'])} "
        f"contacts (of {num(rig[k]['contacts'])}) that carry a positive entry in some stress form a bar-rigid framework "
        "with a strictly positive stress, which suffices by the remark after Theorem~\\ref{thm:main}; the other "
        "contacts carry no stress." for k in sp)
    M["VStressedSentenceSpace"] = " " if M.get("VStressedSentence") else ""
    try:
        famj = json.load(open(os.path.join(OUT, "families.json")))["entries"]
    except FileNotFoundError:
        famj = [{"flexes": 0, "lp_maxmin_stress": float("nan"), "H_max_over_stresses": float("nan"),
                 "H_min_over_stresses": float("nan")}]
    M["VFamFlexes"] = ", ".join(str(e["flexes"]) for e in famj)
    M["VFamFlexesText"] = ("one infinitesimal flex each" if all(e["flexes"] == 1 for e in famj) else
                           "respectively " + M["VFamFlexes"] + " infinitesimal flexes")
    try:
        r48 = json.load(open(os.path.join(OUT, "rigidity_d22_n48.json")))
        M["VFamRecog"] = r48.get("status", "")
    except FileNotFoundError:
        M["VFamRecog"] = "not run"
    M["VFamHAbsMax"] = sci(max(max(abs(e.get("H_max_over_stresses") or 0), abs(e.get("H_min_over_stresses") or 0))
                               for e in famj))
    M["VFamStressMin"] = sci(min(e["lp_maxmin_stress"] for e in famj))
    M["VNStrictNewShort"] = num(len(newstrict))
    M["VNOpenLeft"] = num(len(fam) + len(und))
    us = []
    for k in und:
        j = rig.get(k)
        if j is None:
            us.append(f"$({k[0]},{k[1]})$ was not analysed in exact arithmetic within our time limits")
        else:
            kd_ = j.get("kernel_dim_upper")
            g_ = k[0] * (k[0] - 1) // 2
            if kd_ is not None and kd_ == g_ and j.get("stressed_contacts"):
                us.append(f"the exact core of $({k[0]},{k[1]})$ is bar-rigid, with {num(j['contacts'])} contacts, but "
                          "numerically no stress is positive on every contact (the largest smallest entry of a stress of "
                          f"total one is {sci(abs(j['lp_maxmin_stress']))} in absolute value), and the "
                          f"{num(j['stressed_contacts'])} contacts that carry stress fail our test for bar-rigidity (rank "
                          f"{num(j['rank_mod_p_stressed'])} modulo $q$ against {num(j['cols'] - g_)} needed), so the "
                          "first-order criterion does not apply; it remains open")
            elif kd_ is not None and kd_ == g_:
                us.append(f"the exact core of $({k[0]},{k[1]})$ is bar-rigid, but no exact strictly positive stress "
                          "was found")
            else:
                us.append(f"for $({k[0]},{k[1]})$ the exact analysis was inconclusive ({j.get('status')})")
    M["VUndecidedSentence"] = (" Finally, " + "; ".join(us) + ".") if us else ""
    # the jamming certificate
    M["JamDistLS"] = sci(jam["dist_ls_upper"])
    M["JamDistCoord"] = sci(jam["dist_coord_upper"])
    M["JamContacts"] = num(jam["contacts"])
    M["JamMinContacts"] = num(jam["min_contacts_per_point"])
    M["JamRank"] = num(jam["rank_mod_p"])
    M["JamCols"] = num(jam["cols"])
    M["JamPrime"] = num(jam["prime"])
    M["JamFlexPairs"] = num(jam["flex_S_nonzero_pairs"])
    M["JamStressTypes"] = num(jam["stress_types"])
    M["JamStressDim"] = num(jam["stress_space_dim_on_types"])
    M["JamStressMinRel"] = f"{jam['stress_min_relative']:.2f}"
    M["JamH"] = sci(jam["H_over_norm2_total_stress"])
    M["JamSeconds"] = f"{jam['seconds']:.0f}"
    # exact-only distances
    M["VExOnlyN"] = num(len(dist))
    M["VExOnlyDistLSMax"] = sci(max(e["dist_ls_upper"] for e in dist.values()))
    M["VExOnlyDistCoordMax"] = sci(max(e["dist_coord_upper"] for e in dist.values()))
    # moved rattlers of the isostatic certificates
    passes = [k for k in keys if v4[k] == "strict_iso"]
    moved = [first[k]["dist_rattlers"] for k in passes if (first[k].get("dist_rattlers") or 0) > 0]
    M["VNRattlerMoved"] = num(len(moved))
    M["VRattlerMovedMax"] = sci(max(moved)) if moved else "--"
    M["VRattlerMovedMed"] = sci(float(np.median(moved))) if moved else "--"
    # optimal entries not certified strict
    optk = [k for k in keys if idx[k]["optimal"]]
    miss = [k for k in optk if cls[k] != "strict"]
    M["VNOptMissing"] = num(len(miss))
    M["VOptMissingList"] = entry_list(miss)
    M["VNStrictSecFive"] = num(kc["strict_new"] + kc["strict_second"])
    # the local searches from the NONOPT certificates, evaluated exactly (maxcos_exact.py)
    imp = {}
    for n_ in (576, 791):
        a = json.load(open(os.path.join(OUT, "improved", f"maxcos_d3_n{n_}.json")))
        b = json.load(open(os.path.join(OUT, "improved", f"maxcos_stored_d3_n{n_}.json")))
        from decimal import Decimal, ROUND_CEILING
        va, vb = Decimal(a["max_cos_30"]), Decimal(b["max_cos_30"])
        imp[n_] = (vb - va, va.quantize(Decimal("1e-12"), rounding=ROUND_CEILING))
    M["VImpLowerA"] = sci(float(imp[576][0]))
    M["VImpDispA"] = str(imp[576][1])
    M["VImpLowerB"] = sci(float(imp[791][0]))
    M["VImpLiveB"] = live[(3, 791)][0]
    # the table of the thirteen inner products of the jammed code
    rec = json.load(open(os.path.join(OUT, "recognition_d22_n116.json")))
    pcnt = collections.Counter(rec["pair_class"])
    from fractions import Fraction as Fr
    lines = ["% written by certify/v5/summarize_v5.py; do not edit", "\\begin{tabular}{rrl}", "\\toprule",
             "value & pairs & exact value \\\\", "\\midrule"]

    def fmt(rel):
        r0 = rel[0]
        q = [Fr(-x, r0) for x in rel[1:]]
        den = 1
        for x in q:
            den = den * x.denominator // math.gcd(den, x.denominator)
        z = [int(x * den) for x in q]

        def lin(a, b):
            if a == 0 and b == 0:
                return None
            s = ""
            if a:
                s += str(a)
            if b:
                s += ("+" if (b > 0 and s) else ("-" if b < 0 else "")) + (str(abs(b)) if abs(b) != 1 else "") + "\\sqrt{145}"
            return s
        parts = []
        t0 = lin(z[0], z[1])
        if t0:
            parts.append(t0)
        for (a, b), sym in (((z[2], z[3]), "s_1"), ((z[4], z[5]), "s_2"), ((z[6], z[7]), "s_1s_2")):
            t = lin(a, b)
            if t:
                if a and b:
                    if a < 0 and b < 0:
                        parts.append("-(" + lin(-a, -b) + ")" + sym)
                    else:
                        parts.append("+(" + t + ")" + sym)
                else:
                    tt = t
                    sign = "-" if tt.startswith("-") else "+"
                    tt = tt.lstrip("-")
                    parts.append(sign + (tt if tt not in ("1",) else "") + sym)
        num_ = "".join(parts).replace("+-", "-")
        return "\\frac{" + num_ + "}{" + str(den) + "}" if den != 1 else num_
    for k in range(len(rec["classes"])):
        lines.append(f"${rec['classes'][k]:.6f}$ & {num(pcnt[k])} & ${fmt(rec['class_rel'][k])}$ \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v5_tab_jam.tex"), "w").write("\n".join(lines) + "\n")
    D1, D2 = rec["Delta1_rel"], rec["Delta2_rel"]
    M["JamDeltaOne"] = f"\\frac{{{-D1[2]}\\sqrt{{145}}{'+' if -D1[1] > 0 else '-'}{abs(D1[1])}}}{{{D1[0]}}}"
    M["JamDeltaTwo"] = f"\\frac{{{-D2[2]}\\sqrt{{145}}{'+' if -D2[1] > 0 else '-'}{abs(D2[1])}}}{{{D2[0]}}}"
    M["JamClasses"] = num(len(rec["classes"]))
    # ------------------------------------------------------------------ tables
    lines = ["% written by certify/v5/summarize_v5.py; do not edit", "\\begin{tabular}{lrrrrrrr}", "\\toprule",
             " & & \\multicolumn{3}{c}{exact local minimum nearby} & not a local & \\multicolumn{2}{c}{no certificate} \\\\",
             "\\cmidrule(lr){3-5}\\cmidrule(lr){7-8}",
             "$d$ & entries & strict & open & not strict & minimum & attempted & size limit \\\\", "\\midrule"]
    for lo, hi in RANGES + [(3, 32)]:
        ks = [k for k in keys if lo <= k[0] <= hi]
        c = collections.Counter(cls[k] for k in ks)
        if (lo, hi) == (3, 32):
            lines.append("\\midrule")
        lab = "all" if (lo, hi) == (3, 32) else f"{lo}--{hi}"
        lines.append(f"{lab} & {num(len(ks))} & " + " & ".join(num(c.get(x, 0)) for x in COLS) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v5_tab_census.tex"), "w").write("\n".join(lines) + "\n")
    lines = ["% written by certify/v5/summarize_v5.py; do not edit", "\\begin{tabular}{rrrrrrrrrr}", "\\toprule",
             "$d$ & entries & strict iso. & strict hyper. & strict other & open & not strict & not local & attempted"
             " & too large \\\\", "\\midrule"]
    for d in range(3, 33):
        ks = [k for k in keys if k[0] == d]
        kc_ = collections.Counter(kind[k] for k in ks)
        c = collections.Counter(cls[k] for k in ks)
        other = kc_["strict_exact"] + kc_["strict_new"] + kc_["strict_second"]
        lines.append(f"{d} & {num(len(ks))} & {num(kc_['strict_iso'])} & {num(kc_['strict_hyper'])} & {num(other)} & "
                     f"{num(c['local'])} & {num(c['nonstrict'])} & {num(c['notlocal'])} & {num(c['none_tried'])} & "
                     f"{num(c['none_size'])} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(OUT, "census_by_dimension.tex"), "w").write("\n".join(lines) + "\n")
    lines = ["% written by certify/v5/summarize_v5.py; do not edit", "\\begin{tabular}{llrrrl}", "\\toprule",
             "entry & polynomial of $c$ & core & rattlers & core rank & verdict \\\\", "\\midrule"]
    short = {"not a local minimum: flat core": "not a local min.\\ (flat core)",
             "local minimum, not strict: split contact graph": "local min., not strict (split)",
             "strict local minimum (bar-rigid core)": "strict (bar-rigid)",
             "strict local minimum, jammed, not infinitesimally jammed": "strict, not inf.\\ jammed",
             "undecided (continuous family)": "open (continuous family)", "undecided": "open"}
    for k, poly, info, v, lstat in rows:
        core, rat, rk, ct, dl = info
        lines.append(f"$({k[0]},{k[1]})$ & ${poly}$ & {num(core) if core is not None else '--'} & "
                     f"{num(rat) if rat is not None else '--'} & {rk if rk is not None else '--'} & {short[v]} \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v5_tab_open.tex"), "w").write("\n".join(lines) + "\n")
    # quantitative table of the isostatic certificates, version 4 rows plus the moved rattlers
    A = {k: first[k] for k in passes}
    qrows = [("kappa_point", "sharpness $\\kappa/\\sqrt d$"), ("rho", "radius $\\rho$"),
             ("dist_ls", "distance $D(X_V)$ (least squares)"), ("dist_nonrattlers", "largest coordinate distance"),
             ("r_unique", "uniqueness radius $r_u$ (max norm)"), ("c_err", "radius of the enclosure of $c^*$"),
             ("krawczyk_beta", "largest row sum of (K)"), ("stress_beta", "$\\beta_S$ of (S)"),
             ("cond_inf_upper", "condition bound"), ("stress_min_rel", "smallest relative stress $|E|\\,\\omega_-$"),
             ("unknowns", "unknowns $U$")]
    lines = ["% written by certify/v5/summarize_v5.py; do not edit", "\\begin{tabular}{lrrrrr}",
             "\\toprule", " & min & 10\\% & median & 90\\% & max \\\\", "\\midrule"]
    for key, lab in qrows:
        vals = [float(A[k][key]) for k in passes if A[k].get(key) is not None]
        qq = quant(vals)
        lines.append(lab + " & " + " & ".join(num(int(x)) if key == "unknowns" and x is not None else sci(x) for x in qq)
                     + " \\\\")
    qq = quant(moved)
    lines.append(f"moved rattlers ({num(len(moved))} codes): distance & " + " & ".join(sci(x) for x in qq) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    open(os.path.join(pdir, "v5_tab_quant.tex"), "w").write("\n".join(lines) + "\n")
    with open(os.path.join(pdir, "v5_numbers.tex"), "w") as f:
        f.write("% written by certify/v5/summarize_v5.py; do not edit\n")
        for k_ in sorted(M):
            f.write(f"\\newcommand{{\\{k_}}}{{{M[k_]}}}\n")
    json.dump({"macros": M, "class": {f"{k[0]},{k[1]}": [cls[k], kind[k]] for k in keys},
               "open22": {f"{k[0]},{k[1]}": verdict[k] for k in OPEN22}},
              open(os.path.join(OUT, "summary_v5.json"), "w"), indent=1)
    print(len(M), "macros;", dict(cnt))


def word_n(k):
    w = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve",
         "thirteen", "fourteen", "fifteen"]
    return w[k] if k < len(w) else str(k)


if __name__ == "__main__":
    main()
