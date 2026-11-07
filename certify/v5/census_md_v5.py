"""CENSUS.md for the paper (version 5): a plain summary of the census and of the version 5 analysis, from the outputs.

  python3 -I certify/v5/census_md_v5.py > CENSUS.md
"""
import collections, glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CEN = os.path.join(os.path.dirname(HERE), "census")
sys.path.insert(0, CEN)
import snapshot

LAB = [("strict", "exact strict local minimum nearby"), ("local", "exact local minimum nearby, strictness open"),
       ("nonstrict", "exact local minimum nearby, not strict (split contact graph)"),
       ("notlocal", "the exact configuration nearby is not a local minimum"),
       ("none_tried", "no certificate: every attempted certificate failed or did not apply"),
       ("none_size", "no certificate: every system above the size limit (3,100 unknowns)")]
KINDS = [("strict_iso", "isostatic certificate"), ("strict_hyper", "hyperstatic + exact Gram matrix"),
         ("strict_exact", "exact arithmetic alone (above the size limit)"),
         ("strict_new", "exact arithmetic, version 5 (bar-rigid core)"),
         ("strict_second", "second-order criterion, version 5 (jammed, not infinitesimally jammed)")]


def f(k):
    return f"{k:,}"


def main():
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index() if e["listed"]}
    S = json.load(open(os.path.join(HERE, "out", "summary_v5.json")))
    cls = {tuple(map(int, k.split(","))): v for k, v in S["class"].items()}
    keys = sorted(cls)
    snap = json.load(open(os.path.join(CEN, "snapshot_files.json")))
    out = ["# Census of Cohn's archived table of spherical codes (paper version 5)", "",
           "Data set: DSpace@MIT handle 1721.1/153543, data as of 31 July 2026 (downloaded 6 November 2026; DSpace lists "
           "the item as last modified 4 August 2026). Every file's MD5 agrees with DSpace's.", "",
           "| file | bytes | SHA-256 |", "|---|---:|---|"]
    for name in ["table.txt"] + [f"data-{k}.txt" for k in range(1, 13)] + ["data.txt"]:
        out.append(f"| {name} | {f(snap[name]['bytes'])} | `{snap[name]['sha256']}` |")
    nfl = sum(idx[k]["flagged"] for k in keys)
    nopt = sum(idx[k]["optimal"] for k in keys)
    out += ["", f"Census: the {f(len(keys))} entries that `table.txt` lists with 3 <= d <= 32 and N <= 1024 ({f(nfl)} marked "
            f"\"needs more optimization\", {nopt} marked \"provably optimal\"). The census rule and certificates are those of "
            "version 4 (paper, Sections 3 and 4); version 5 adds the exact analysis of the 22 open entries with an exact "
            "cosine (paper, Section 5), certified distances for the exact-only confirmations, and reports the diagnostic "
            "of Proposition 4.2 apart from the classes.", "",
            "The objects of the theorems are exact configurations near the stored codes, never the stored codes; every "
            "distance between them is in the outputs.", "", "## Counts by class", "",
            "| class | all | marked | unmarked |", "|---|---:|---:|---:|"]
    for c, lab in LAB:
        ks = [k for k in keys if cls[k][0] == c]
        fl = sum(idx[k]["flagged"] for k in ks)
        out.append(f"| {lab} | {f(len(ks))} | {f(fl)} | {f(len(ks) - fl)} |")
    out.append(f"| total | {f(len(keys))} | {f(nfl)} | {f(len(keys) - nfl)} |")
    out += ["", "Strict entries by certificate:", ""]
    for kd, lab in KINDS:
        out.append(f"- {lab}: {f(sum(1 for k in keys if cls[k][1] == kd))}")
    M = S["macros"]
    out += ["", f"Diagnostic (Proposition 4.2, a neighbourhood statement, not a class): holds for {M['VNDiag']} entries "
            f"without a certificate; elementary (forest) for {M['VNDiagForest']}; radius below 1e-6 for "
            f"{M['VNDiagSmallRad']}; marked \"needs more optimization\" {M['VNDiagFlagged']}.", "",
            "## Counts by dimension", "",
            "| d | entries | strict | open | not strict | not local | attempted | too large |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for d in range(3, 33):
        ks = [k for k in keys if k[0] == d]
        c = collections.Counter(cls[k][0] for k in ks)
        out.append(f"| {d} | {f(len(ks))} | " + " | ".join(f(c.get(x, 0)) for x, _ in LAB) + " |")
    out += ["", "## The 22 open entries with an exact cosine (version 5)", "", "| entry | verdict |", "|---|---|"]
    for k, v in S["open22"].items():
        out.append(f"| ({k}) | {v} |")
    jam = json.load(open(os.path.join(HERE, "out", "jam_d22_n116.json")))
    out += ["", "## The jammed code that is not infinitesimally jammed", "",
            f"(22,116): {jam['status']}. Field Q(sqrt145, s1, s2) of degree 8; {jam['contacts']} contacts; rank of the "
            f"rigidity matrix at least {jam['rank_mod_p']} = 116*22 - 231 - 1 (mod {jam['prime']}); one flex (nonzero on "
            f"{jam['flex_S_nonzero_pairs']} pairs of the Gram derivative); a positive stress on {jam['stress_types']} contact "
            f"types; H_omega(F)/(|F|^2 sum omega) = {jam['H_over_norm2_total_stress']:.3g}; distance to the stored code "
            f"{jam['dist_ls_upper']:.2g} (least squares), {jam['dist_coord_upper']:.2g} (coordinates).", "",
            "## Entries marked provably optimal", "", "| entry | class |", "|---|---|"]
    for k in keys:
        if idx[k]["optimal"]:
            out.append(f"| ({k[0]},{k[1]}) | {dict(LAB)[cls[k][0]]} ({cls[k][1]}) |")
    dist = json.load(open(os.path.join(HERE, "out", "distances_exact_only.json")))["entries"]
    out += ["", "## Certified distances of the exact-only confirmations (version 5)", "",
            "| entry | least-squares distance | coordinate distance |", "|---|---:|---:|"]
    for e in dist:
        out.append(f"| ({e['d']},{e['n']}) | {e['dist_ls_upper']:.2g} | {e['dist_coord_upper']:.2g} |")
    out += ["", "## Programs", "", "Every output records the SHA-256 of the programs that produced it. Version 4: "
            "`certify/census/MANIFEST_v4.tsv`; version 5: `certify/v5/MANIFEST_v5.tsv`.", "",
            f"Macros of the paper: certify/census/summary_v4.json and certify/v5/out/summary_v5.json ({len(M)} values)."]
    print("\n".join(out).replace("{,}", ","))


if __name__ == "__main__":
    main()
