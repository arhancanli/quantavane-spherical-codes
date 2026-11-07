# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Numbers and a table for Section 7 of the note, from the outputs of the searches (family.py, family_table.py) and of
the exact programs run on the second code X' (xprime.gp, structure.gp and check_xprime.gp on data_xprime/).

  python3 -I certify/note/family_summary.py NOTEDIR

Reads out/family/family_q9.json, family_q25.json, family_qother.json, family_q13.json, family_qmore.json,
xstar_restarts.json, table_near.json, table_near_b.json and
out/structure_xprime.json, out/check_xprime.json, out/xprime.json, out/jam_xprime.json and the red controls
out/jam_xprime_red_gram.json, out/jam_xprime_red_ratio.json. Writes NOTEDIR/family_numbers.tex,
NOTEDIR/tab_family.tex and out/family/summary_family.json. Floating-point numbers are only displayed; every statement
about X' that the note calls proved comes from check_xprime.gp (exact).
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
O = os.path.join(HERE, "out")
FD = os.path.join(O, "family")


def sci(x, k=1):
    e = math.floor(math.log10(abs(x)))
    return f"{x / 10 ** e:.{k}f}" + r"\cdot10^{%d}" % e


def outcome(an):
    """a short verdict for the best configuration of a search case"""
    if an is None:
        return "---"
    if an.get("analysis_error"):
        return "analysis failed"
    if an["rattler_orbits"]:
        return "rattlers"
    if not an.get("strictly_positive_stress"):
        return "no positive stress"
    fl = an.get("flexes_mod_rotations")
    if fl is None:
        lin = an.get("linear_sectors", {})
        if all(v["flex_dim_mod_rotations"] == 0 for v in lin.values()):
            return "no flex in a linear character"
        return "flex in a linear character"
    if fl == 0:
        return "infinitesimally jammed"
    h = an.get("H_on_flexes_maxmin_stress_eigs") or []
    return "%d flex%s, $H%s0$" % (fl, "" if fl == 1 else "es", ">" if h and min(h) > 0 else "\\not>")


def main():
    notedir = sys.argv[1]
    mac, rows, summ = {}, [], {}
    cases = []
    for f in ("family_q9.json", "family_q25.json", "family_qother.json", "family_q13.json", "family_qmore.json",
              "family_q16.json"):
        p = os.path.join(FD, f)
        if os.path.exists(p):
            cases += json.load(open(p))["results"]
    cases.sort(key=lambda c: c["case"])
    names = {"PGammaL": r"$\mathrm{P\Gamma L}(2,%d)$", "PGL": r"$\mathrm{PGL}(2,%d)$", "PSigmaL": r"$\mathrm{P\Sigma L}(2,%d)$"}
    eps_name = {"M": r"$\varepsilon$", "det": r"$\chi\circ\det$", "frob": r"$(-1)^{\phi}$"}
    for c in cases:
        if c.get("skipped"):
            continue
        an = c["analysed"][0] if c.get("analysed") else None
        best = c["optimised"][0] if c.get("optimised") else None
        kinds = "$" + "\\cup ".join(k if k in ("P", "Q") else ("A" if k == "A" else "B_{%s}" % k[1:]) for k in c["kinds"]) + "$"
        sizes = "+".join(str(s) for s in c["orbit_sizes"])
        isprime = all(c["q"] % k for k in range(2, c["q"]))
        grp = (r"$\mathrm{PGL}(2,%d)$" % c["q"]) if (isprime and c["group"] == "PGammaL") else names[c["group"]] % c["q"]
        dims = "+".join(str(x) for x in best["profile_dims"]) if best else "---"
        mu = f"{best['largest_cosine']:.6f}" if best else "---"
        verdict = outcome(an)
        if c["q"] == 9 and an and abs(an.get("largest_cosine", 0) - (math.sqrt(145) - 5) / 60) < 1e-9:
            mm = an.get("maxmin_invariant_stress_times_E", 0)       # 0.2658 for X*, 0.2817 for X'
            verdict = "$X^*$" if abs(mm - 0.26581) < 1e-4 else "$X'$" if abs(mm - 0.28167) < 1e-4 else verdict
        rows.append(f"{grp} & {kinds} & {sizes} & {eps_name[c['eps']]} & {c['profiles']} & {dims} & ${mu}$ & {verdict} \\\\")
        summ[str(c["case"])] = {"q": c["q"], "group": c["group"], "kinds": c["kinds"], "N": c["N"], "profiles": c["profiles"],
                                "best": best, "verdict": verdict}
    with open(os.path.join(notedir, "tab_family.tex"), "w") as f:
        f.write("% written by certify/note/family_summary.py; do not edit\n\\begin{tabular}{llllrlll}\n\\toprule\n"
                "group & $\\Omega$ & sizes & sign & profiles & best $\\R^d$ & largest cosine & best configuration \\\\\n"
                "\\midrule\n" + "\n".join(rows) + "\n\\bottomrule\n\\end{tabular}\n")
    # ------------------------------------------------------------ restarts and the analyses of X* and X'
    rs = json.load(open(os.path.join(FD, "xstar_restarts.json")))
    mac["FamStarts"] = rs["starts"]
    mac["FamReachStar"] = rs["starts_reaching_Xstar"]
    mac["FamReachPrime"] = rs["starts_reaching_Xprime"]
    mac["FamOther"] = rs["starts"] - rs["starts_reaching_Xstar"] - rs["starts_reaching_Xprime"]
    others = rs.get("other_local_minima", [])
    mac["FamOtherMin"] = f"{min(others):.4f}" if others else "---"
    for kind, tag in (("X*", "Star"), ("X'", "Prime")):
        a = rs["analysis"][kind]
        g = a["singular_values_around_the_gap"]
        mac[f"FamGapLow{tag}"] = sci(g[0])
        mac[f"FamGapHigh{tag}"] = f"{g[1]:.3f}"
        mac[f"FamH{tag}"] = sci(a["H_on_flexes_maxmin_stress_eigs"][0], 3)
        mac[f"FamMaxmin{tag}"] = f"{a['maxmin_invariant_stress_times_E']:.4f}"
        mac[f"FamFlex{tag}"] = a["flexes_mod_rotations"]
    mac["FamGapLow"], mac["FamGapHigh"], mac["FamH"] = mac["FamGapLowStar"], mac["FamGapHighStar"], mac["FamHStar"]
    # ------------------------------------------------------------ the exact results on X'
    xp = json.load(open(os.path.join(O, "xprime.json")))["results"]
    sx = json.load(open(os.path.join(O, "structure_xprime.json")))["results"]
    cx = json.load(open(os.path.join(O, "check_xprime.json")))
    assert cx["results"].get("STATUS") == "ALL CHECKS PASSED" and not cx["error"]
    mac["XpChanged"] = xp["changed_pairs"]
    mac["XpAut"] = cx["results"]["automorphisms_by_backtracking"]
    mac["XpRankModp"] = cx["results"]["rank_mod_p"]
    mac["XpPrime"] = f'{int(cx["results"]["prime"]):,}'.replace(",", "{,}")
    mac["XpMinDeg"] = cx["results"]["min_contacts_per_point"]
    mac["XpH"] = f'{float(cx["results"]["H_numeric"].replace(" E", "E")):.4f}'
    mac["XpRatio"] = f'{float(sx["maxmin_ratio_QQ_over_PP"].replace(" E", "E")):.3f}'
    mac["XpThreshold"] = f'{float(sx["threshold_ratio_QQ_over_PP_numeric"].replace(" E", "E")):.4f}'
    mac["XpMaxmin"] = f'{float(sx["maxmin_min_entry_times_E"].replace(" E", "E")):.4f}'
    mac["XpAlphaSign"] = "<" if sx["alpha_sign"].strip() == "-1" else ">"
    mac["XpBetaSign"] = "<" if sx["beta_sign"].strip() == "-1" else ">"
    mac["XpSeconds"] = round(float(cx["seconds"]))
    # the second, independent check of X' (jam_xprime.py, Python/Arb) and its red controls
    jx = json.load(open(os.path.join(O, "jam_xprime.json")))
    assert jx["status"].startswith("JAMMED") and jx["agrees_with_xprime_gp"]
    mac["XpPyRank"] = f'{jx["rank_mod_p"]:,}'.replace(",", "{,}")
    mac["XpPyPrime"] = f'{jx["prime"]:,}'.replace(",", "{,}")
    for m_ in ("gram", "ratio"):
        jr = json.load(open(os.path.join(O, "jam_xprime_red_%s.json" % m_)))
        assert jr["status"].startswith("FAILED")
    mac["XpRedGram"] = json.load(open(os.path.join(O, "jam_xprime_red_gram.json")))["status"].replace("FAILED: ", "")
    mac["XpRedRatio"] = json.load(open(os.path.join(O, "jam_xprime_red_ratio.json")))["status"].replace("FAILED: ", "")
    # ------------------------------------------------------------ the entries of the table near (22,116)
    tp = os.path.join(FD, "table_near.json")
    if os.path.exists(tp):
        t = json.load(open(tp))
        nmax = t["range"]["nmax"]
        tb = os.path.join(FD, "table_near_b.json")          # the entries with 101 <= N <= 130 (time limit per group)
        if os.path.exists(tb) and "seconds" in json.load(open(tb)):          # only a completed run is used
            t2 = json.load(open(tb))
            nmax = t2["range"]["nmax"]
            t = {"entries": t["entries"] + t2["entries"], "structured": t["structured"] + t2["structured"],
                 "grouped": t["grouped"] + t2["grouped"]}
        mac["TabNmax"] = nmax
        g = t["grouped"]
        mac["TabAutTimeout"] = sum(1 for r in t["structured"] if r.get("aut_timeout_seconds"))
        mac["TabAutTimeoutList"] = ", ".join("$(%d,%d)$" % (r["d"], r["n"]) for r in t["structured"]
                                             if r.get("aut_timeout_seconds")) or "none"
        mac["TabAutText"] = ("" if not mac["TabAutTimeout"] else
                             " For %s the symmetry group was not found within the time limit (300 seconds), so %s not "
                             "analysed." % (mac["TabAutTimeoutList"], "it was" if mac["TabAutTimeout"] == 1 else "they were"))
        mac["TabEntries"] = t["entries"]
        mac["TabStructured"] = len(t["structured"])
        mac["TabGrouped"] = len(g)
        mac["TabFlex"] = sum(1 for r in g if r.get("flexes_mod_rotations", 0) > 0)
        mac["TabRattlers"] = sum(1 for r in g if r.get("rattlers", 0) > 0)
        mac["TabNoStress"] = sum(1 for r in g if not r.get("strictly_positive_stress", False))
        flex = [r for r in g if r.get("flexes_mod_rotations", 0) > 0]
        mac["TabFlexList"] = ", ".join("$(%d,%d)$" % (r["d"], r["n"]) for r in flex) or "none"
        def one(r):
            bits = []
            if "continuous" in r["notes"]:
                bits.append("a member of a continuous family, as the table notes")
            if r["flexes_mod_rotations"] > 1:
                bits.append("%d flexes" % r["flexes_mod_rotations"])
            if r.get("rattlers"):
                bits.append("%d rattler%s" % (r["rattlers"], "" if r["rattlers"] == 1 else "s"))
            if not r.get("strictly_positive_stress"):
                bits.append("no strictly positive stress")
            h = r.get("H_eigs_min_max")
            if h and r.get("strictly_positive_stress"):
                bits.append("$H>0$" if h[0] > 1e-9 else "$H=0$ on a flex" if abs(h[0]) <= 1e-9 else "$H<0$ on a flex")
            return "$(%d,%d)$%s" % (r["d"], r["n"], " (" + ", ".join(bits) + ")" if bits else "")
        if not flex:
            txt = "no infinitesimal flex modulo rotations of the contact framework of the core (the points that are not rattlers)"
        else:
            txt = ("infinitesimal flexes modulo rotations of the contact framework of the core (the points that are not "
                   "rattlers) only for " + ", ".join(one(r) for r in flex) + "; $H$ is computed for a max-min stress")
        mac["TabFlexText"] = txt
        summ["table"] = {"entries": t["entries"], "structured": len(t["structured"]), "grouped": len(g),
                         "with_flexes": [(r["d"], r["n"], r["flexes_mod_rotations"], r.get("rattlers"),
                                          r.get("strictly_positive_stress"), r.get("H_eigs_min_max")) for r in flex]}
    with open(os.path.join(notedir, "family_numbers.tex"), "w") as f:
        f.write("% written by certify/note/family_summary.py from the outputs of the searches; do not edit\n")
        for k, v in mac.items():
            f.write("\\newcommand{\\N%s}{%s}\n" % (k, v))
    summ["macros"] = mac
    json.dump(summ, open(os.path.join(FD, "summary_family.json"), "w"), indent=1, default=str)
    print(json.dumps(mac))


if __name__ == "__main__":
    main()
