# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""The numbers of the note's count across the table (Section 4; floating point), from the outputs of symcount.py.

  python3 -I certify/note/symcount_summary.py NOTEDIR

Reads out/symcount/symcount_*.jsonl (run_symcount.sh) and writes NOTEDIR/symcount_numbers.tex (macros \\NSc...),
NOTEDIR/tab_symcount.tex and out/symcount/summary_symcount.json. Nothing is computed here beyond counting.
"""
import collections, glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
O = os.path.join(HERE, "out", "symcount")
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "census"))
from snapshot import load_index


def f(n):
    return "{:,}".format(n).replace(",", "{,}")


def entries(lst):
    return ", ".join("$(%d,%d)$" % k for k in sorted(lst))


def main():
    notedir = sys.argv[1]
    recs = {}
    for p in sorted(glob.glob(os.path.join(O, "symcount_*.jsonl"))):
        for line in open(p):
            r = json.loads(line)
            recs[(r["d"], r["n"])] = r
    index = [(e["d"], e["n"]) for e in load_index()]
    missing = [k for k in index if k not in recs]
    cat = collections.Counter(r["category"] for r in recs.values())
    by = lambda c: sorted(k for k, r in recs.items() if r["category"] == c)
    trivial, forced = by("forced, trivial group"), by("forced by symmetry")
    notforced, norank, undecided = by("flexes not forced"), by("not forced (rank not computed)"), by("undecided")

    def examined(r):                        # the stress test of symcount.py: unflagged codes with N d <= 2600
        return not r["flagged"] and r["d"] * r["n"] <= 2600

    def positive_stress(r):                 # smallest entry at least 1e-3 of the average entry: positive beyond the numerics
        return examined(r) and r.get("maxmin_times_E", 0) > 1e-3

    def weak_stress(r):                     # positive, but only at the level of the numerics
        return examined(r) and 1e-6 < r.get("maxmin_times_E", 0) <= 1e-3

    def prestress(r):                       # H > 0 on every flex direction for the max-min stress, or (f = 1) for some
        e = r.get("H_maxmin_eigs")        # (beyond 1e-6: with sum omega = 1 and unit flexes; zero energy is ~1e-15)
        return positive_stress(r) and e is not None and (e[0] > 1e-6 or r.get("H_best_positive_stress", -1) > 1e-6)
    flex = [k for k, r in recs.items() if r.get("f", 0) > 0]
    stressed = sorted(k for k in flex if positive_stress(recs[k]))
    cands = sorted(k for k in stressed if prestress(recs[k]))
    weak = sorted(k for k in flex if weak_stress(recs[k]))
    weak_pos = sorted(k for k in weak if (recs[k].get("H_maxmin_eigs") or [-1])[0] > 1e-6)
    tested = sorted(k for k in flex if "H_maxmin_eigs" in recs[k] or recs[k].get("s", 1) == 0)
    soft = collections.defaultdict(list)
    for k in stressed:
        if "soft_modes" in recs[k]:
            soft[k[0]].append(recs[k]["soft_modes"])
    star = recs.get((22, 116), {})
    other = sorted((g, m) for g, m in star.get("counts_by_degree", []) if g > 1)
    bydeg = collections.defaultdict(list)
    for g, m in other:
        bydeg[g].append(m)
    star_text = "; ".join("$%s$ (degree %d)" % (",".join(str(m) for m in sorted(v)), g) for g, v in sorted(bydeg.items()))
    soft_entries = sorted((k, recs[k]["soft_modes"]) for k in stressed if "soft_modes" in recs[k])
    only_const = [k for k, s in soft_entries if s == 1]
    more = [k for k, s in soft_entries if s > 1]
    lowd_more = [k for k in more if k[0] <= 4]
    def sci(x):
        mant, ex = ("%.1e" % x).split("e")
        return "%s\\cdot10^{%d}" % (mant, int(ex))
    wt = []
    for k in weak_pos:
        r = recs[k]
        forced_by = (("its group of order %d" % r["group_order"]) if r.get("group_order", 1) > 1
                     else "symmetry (it has none)")
        txt = ("$(%d,%d)$ has %s flex%s modulo rotations, not forced by %s; for the stress found by linear programming $H$ "
               "is positive on %s ($\\approx%.4f$ with $\\sum_E\\omega=1$ and unit flexes), but the smallest entry of that "
               "stress is only $%s$ of the average"
               % (k[0], k[1], r["f"] if r["f"] > 1 else "one", "es" if r["f"] > 1 else "", forced_by,
                  "the flexes" if r["f"] > 1 else "the flex", r["H_maxmin_eigs"][0], sci(r["maxmin_times_E"])))
        wt.append(txt)
    if wt:
        words = {2: "Two", 3: "Three", 4: "Four", 5: "Five"}
        weak_text = ("One more entry comes close: " if len(wt) == 1 else
                     "%s more entries come close: " % words.get(len(wt), str(len(wt)))) + \
            "; ".join(wt) + ". We leave open whether %s jammed." % ("it is" if len(wt) == 1 else "they are")
    else:
        weak_text = ""
    if stressed:
        dmin = min(k[0] for k in stressed)
        # flexes in lower dimensions whose stresses were not examined (no linear program solved, or too large)
        unknown_low = sorted(k for k in flex if k[0] < dmin and examined(recs[k]) and recs[k].get("s", 1) > 0
                             and "maxmin_times_E" not in recs[k])
        unknown_low_text = (" (apart from %s, where no stress was examined)" % entries(unknown_low)) if unknown_low else ""
        soft_text = ("In the archived table (Section~\\ref{sec:flex}), of the entries with a flex whose stresses we "
                     "examined, the %d with a stress whose smallest entry is at least $10^{-3}$ of the average all lie in "
                     "dimension %d or more: %s (with the number of soft modes, for the stress with the largest smallest "
                     "entry). In lower dimensions the flexes we examined come only with stresses that vanish, or nearly "
                     "vanish, on some contacts, so that Theorem~\\ref{thm:second} cannot apply%s."
                     % (len(stressed), dmin,
                        ", ".join("$(%d,%d)$%s" % (k[0], k[1], (" (%d)" % recs[k]["soft_modes"]) if "soft_modes" in recs[k]
                                                    else "") for k in stressed),
                        unknown_low_text))
    else:
        soft_text = ""
    mac = {
        "ScEntries": f(len(recs)), "ScMissing": f(len(missing)),
        "ScRattlers": f(cat["rattlers"]), "ScNotSpanning": f(cat["not spanning"]),
        "ScTooLarge": f(cat["too large"]),
        "ScNoRattlers": f(len(recs) - cat["rattlers"] - cat["not spanning"] - cat["too large"]),
        "ScBarRigid": f(cat["bar-rigid"]), "ScTrivial": f(len(trivial)),
        "ScTrivialFlagged": f(sum(1 for k in trivial if recs[k]["flagged"])),
        "ScTrivialNoStress": f(sum(1 for k in trivial if recs[k].get("s", 1) == 0)),
        "ScForcedSym": f(len(forced)), "ScForcedSymList": entries(forced),
        "ScNotForced": f(len(notforced)), "ScNotForcedList": entries(notforced),
        "ScNoRank": f(len(norank)), "ScUndecided": f(len(undecided)), "ScErrors": f(cat["error"]),
        "ScFlex": f(len(flex)), "ScFlexStressed": f(len(stressed)), "ScTested": f(len(tested)),
        "ScExamined": f(sum(1 for k in flex if examined(recs[k]))),
        "ScCandidates": f(len(cands)), "ScCandidateList": entries(cands) if cands else "none",
        "ScCandidateOtherList": entries([k for k in cands if k != (22, 116)]) or "none",
        "ScCandidateOtherCount": f(len([k for k in cands if k != (22, 116)])),
        "ScCandidateDetail": "; ".join(
            "$(%d,%d)$, with %s flex%s modulo rotations%s" % (
                k[0], k[1], recs[k]["f"] if recs[k]["f"] > 1 else "one", "es" if recs[k]["f"] > 1 else "",
                (", forced by its symmetry group of order %s" % f(recs[k]["group_order"])) if recs[k].get("forced", 0) > 0
                else (" and no symmetry" if recs[k].get("group_order", 1) == 1 else
                      ", not forced by its symmetry group of order %s" % f(recs[k]["group_order"])))
            for k in cands if k != (22, 116)) or "none",
        "ScStarExamined": "yes" if (22, 116) in cands else "no",
        "ScWeak": f(len(weak)), "ScWeakList": entries(weak) if weak else "none",
        "ScWeakPos": f(len(weak_pos)), "ScWeakPosList": entries(weak_pos) if weak_pos else "none",
        "ScSoftEntries": f(len(soft_entries)), "ScSoftOnlyConst": f(len(only_const)), "ScSoftMore": f(len(more)),
        "ScSoftMoreList": entries(more) if more else "none",
        "StarOtherCounts": star_text,
        "ScWeakText": weak_text,
        "ScSoftText": soft_text,
    }
    rows = []
    for k in sorted(set(forced) | set(cands)):
        r = recs[k]
        types = ", ".join("%d$\\times$%d" % (m, g) if m > 1 else "%d" % g for g, m in r.get("forced_types", []))
        e = r.get("H_maxmin_eigs")
        rows.append("$(%d,%d)$ & %s & %s & %s & %d & %s & %s & %s\\\\" % (
            k[0], k[1], f(r["group_order"]) if "group_order" in r else "--", types or "--", r.get("m0", "--"),
            r.get("f", -1), "yes" if positive_stress(r) else "no",
            ("%.2g" % e[0]) if e else "--", "yes" if r["flagged"] else ""))
    tab = ["\\begin{tabular}{@{}lrrrrccc@{}}", "\\toprule",
           "entry & $|H|$ & forced (degree) & Maxwell & flexes & stress $>0$ & min $H$ & flagged\\\\", "\\midrule"]
    tab += rows + ["\\bottomrule", "\\end{tabular}"]
    with open(os.path.join(notedir, "tab_symcount.tex"), "w") as fo:
        fo.write("% Written by certify/note/symcount_summary.py from out/symcount/ (floating point).\n" + "\n".join(tab) + "\n")
    with open(os.path.join(notedir, "symcount_numbers.tex"), "w") as fo:
        fo.write("% Written by certify/note/symcount_summary.py from out/symcount/ (floating point).\n")
        for kk, v in sorted(mac.items()):
            fo.write("\\newcommand{\\N%s}{%s}\n" % (kk, v))
    summ = {"program": "symcount_summary.py", "records": len(recs), "missing": missing, "categories": dict(cat),
            "forced_by_symmetry": forced, "flexes_not_forced": notforced, "undecided": undecided,
            "rank_not_computed": norank, "flex_with_positive_stress": stressed, "prestress_candidates": cands,
            "weak_stress": weak, "weak_stress_positive_energy": weak_pos,
            "soft_modes": {"%d_%d" % k: s for k, s in soft_entries}, "low_dimension_more_soft_modes": lowd_more,
            "star_counts_by_degree": star.get("counts_by_degree")}
    json.dump(summ, open(os.path.join(O, "summary_symcount.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k not in ("soft_modes",)}))


if __name__ == "__main__":
    main()
