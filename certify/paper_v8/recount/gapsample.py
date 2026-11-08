"""The sample of gap-based proposals for gapcut.py certify (census version 7). Untrusted bookkeeping only.

  python3 -I gapsample.py SCAN.jsonl[,SCAN2.jsonl] OUT.tsv

Among the scanned entries with a non-hypostatic cut and at most 1,500 unknowns: for each census status (FLEX, failed)
and each spread band (<= 1e-7, 1e-7 to 1e-6, 1e-6 to 1e-5, 1e-5 to 1e-4), ten entries at random (seed 7), or all if
fewer. Columns: d, n, tol (inside the gap), stratum.
"""
import json, random, sys

BANDS = [(0, 1e-7), (1e-7, 1e-6), (1e-6, 1e-5), (1e-5, 1e-4)]


def main():
    rs = []
    for p in sys.argv[1].split(","):
        rs += [json.loads(x) for x in open(p) if x.strip()]
    rng = random.Random(7)
    out = []
    for st in ("FLEX", "failed"):
        for lo, hi in BANDS:
            c = [r for r in rs if r.get("found") and r["census_status"] == st and r["unknowns"] <= 1500
                 and (lo < r["spread"] <= hi or (lo == 0 and r["spread"] == 0))]
            c.sort(key=lambda r: (r["d"], r["n"]))
            for r in sorted(rng.sample(c, min(10, len(c))), key=lambda r: (r["d"], r["n"])):
                out.append((r["d"], r["n"], r["tol"], f"{st}, spread {lo:g}-{hi:g}"))
    with open(sys.argv[2], "w") as f:
        f.write("# d\tn\ttol\tstratum\n")
        for d, n, t, s in out:
            f.write(f"{d}\t{n}\t{t!r}\t{s}\n")
    print(len(out))


if __name__ == "__main__":
    main()
