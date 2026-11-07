"""Which entries of the live table to fetch for the currency check of census version 7 (untrusted bookkeeping).

  python3 -I select_live.py LIVE.tsv CATEGORIES.tsv OUT.tsv

Every entry with a positive verdict in the census (classes strict, local, nonstrict, notlocal of summary_v5.json; the
symmetry upgrades of version 6 are inside 'local') whose live cosine differs from the archived one, and a stratified
random sample (seed 7, ten per category and dimension band) of the negative-certificate entries whose live cosine is
lower than the archived one. The live cosine and the archived cosine are compared as the 12-decimal strings shown.
"""
import json, os, random, sys
from fractions import Fraction
CERT = os.environ.get("QV_CERTIFY", os.path.expanduser("~/quantavane-spherical-codes/certify"))
sys.path.insert(0, os.path.join(CERT, "census"))
import snapshot

BANDS = [(3, 5), (6, 10), (11, 17), (18, 24), (25, 32)]


def main():
    live = {}
    for line in open(sys.argv[1]):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        live[(int(p[0]), int(p[1]))] = p[2]
    cats = {}
    for line in open(sys.argv[2]):
        if line.startswith("#"):
            continue
        p = line.split("\t")
        cats[(int(p[0]), int(p[1]))] = p[2]
    cls = json.load(open(os.path.join(CERT, "v5", "out", "summary_v5.json")))["class"]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    pos, neg = [], {}
    for k, v in cls.items():
        d, n = map(int, k.split(","))
        if (d, n) not in live:
            continue
        a, l = Fraction(idx[(d, n)]["cos"]), Fraction(live[(d, n)])
        if v[0] in ("strict", "local", "nonstrict", "notlocal") and a != l:
            pos.append((d, n, v[0], "positive verdict, live cosine differs"))
        elif v[1] == "flex" and l < a:
            band = next(i for i, (lo, hi) in enumerate(BANDS) if lo <= d <= hi)
            neg.setdefault((cats[(d, n)], band), []).append((d, n))
    rng = random.Random(7)
    sample = []
    for key in sorted(neg):
        lst = sorted(neg[key])
        for d, n in sorted(rng.sample(lst, min(10, len(lst)))):
            sample.append((d, n, "flex", f"sample: {key[0]}, d {BANDS[key[1]][0]}-{BANDS[key[1]][1]}"))
    with open(sys.argv[3], "w") as f:
        f.write("# d\tn\tcensus class\treason\n")
        for r in sorted(pos) + sample:
            f.write("\t".join(map(str, r)) + "\n")
    print(len(pos), len(sample))


if __name__ == "__main__":
    main()
