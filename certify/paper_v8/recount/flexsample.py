"""The stratified sample of negative certificates for flexexact.py (census version 7). Untrusted bookkeeping only.

  python3 -I flexsample.py FLEXGRAM.jsonl[,FLEXGRAM2.jsonl...] OUT.tsv

Every FLEX entry gets one of four categories (from the census record and the snapshot flags):
  after-failed  an isostatic or hyperstatic attempt was made and failed;
  above-cap     no such attempt was made, but one was planned and skipped above the size cap;
  no-proposal-unflagged / no-proposal-flagged   the proposal was hypostatic or empty at every tolerance.
In each category and each dimension band (3-5, 6-10, 11-17, 18-24, 25-32) the entries with the smallest, the median
and the largest delta (of flexgram) among those with rows * columns <= MAX_CELLS are taken, plus five entries checked
independently. The categories of every FLEX entry are also written to OUT.tsv.categories (for the paper).
"""
import json, os, sys
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
sys.path.insert(0, os.path.join(CERT, "census"))
import snapshot

MAX_CELLS = 12_000_000
BANDS = [(3, 5), (6, 10), (11, 17), (18, 24), (25, 32)]
NAMED = [(12, 841), (14, 73), (26, 55), (24, 75), (11, 105)]


def category(j, flagged):
    made = any(a.get("planned_kind") in ("iso", "hyper") for a in j["attempts"])
    skipped = any(s["kind"] in ("iso", "hyper") for s in j.get("skipped_above_cap", []))
    if made:
        return "after-failed"
    if skipped:
        return "above-cap"
    return "no-proposal-flagged" if flagged else "no-proposal-unflagged"


def main():
    fg = {}
    for p in sys.argv[1].split(","):
        for line in open(p):
            r = json.loads(line)
            fg[(r["d"], r["n"])] = r
    out = sys.argv[2]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    cats = {}
    for (d, n), r in sorted(fg.items()):
        j = json.load(open(os.path.join(CERT, "census", "results_v4", f"d{d}_n{n}.json")))
        cats[(d, n)] = category(j, idx[(d, n)]["flagged"])
    with open(out + ".categories", "w") as f:
        f.write("# d\tn\tcategory\tdelta\ttau_prime\th\trows\tcolumns\n")
        for (d, n), c in sorted(cats.items()):
            r = fg[(d, n)]
            f.write(f"{d}\t{n}\t{c}\t{r['delta']!r}\t{r['tau_prime']!r}\t{r['h']!r}\t"
                    f"{r['contacts'] + r['points_in_contacts']}\t{r['points_in_contacts'] * d}\n")
    pick = {}
    for c in sorted(set(cats.values())):
        for lo, hi in BANDS:
            cand = [(fg[k]["delta"], k) for k, cc in cats.items() if cc == c and lo <= k[0] <= hi
                    and (fg[k]["contacts"] + fg[k]["points_in_contacts"]) * fg[k]["points_in_contacts"] * k[0]
                    <= MAX_CELLS]
            if not cand:
                continue
            cand.sort()
            for why, (dl, k) in (("smallest delta", cand[0]), ("median delta", cand[len(cand) // 2]),
                                 ("largest delta", cand[-1])):
                pick.setdefault(k, (c, f"{why} in d {lo}-{hi}"))
    for k in NAMED:
        if k in cats:
            pick.setdefault(k, (cats[k], "checked independently"))
    with open(out, "w") as f:
        f.write("# d\tn\tcategory\treason\n")
        for k in sorted(pick, key=lambda k: (fg[k]["contacts"] + fg[k]["points_in_contacts"]) * fg[k]["points_in_contacts"] * k[0]):
            f.write(f"{k[0]}\t{k[1]}\t{pick[k][0]}\t{pick[k][1]}\n")
    from collections import Counter
    print(Counter(cats.values()), len(pick))


if __name__ == "__main__":
    main()
