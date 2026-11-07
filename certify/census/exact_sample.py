"""The sample of certificates re-checked by exact_big.py (paper version 4), drawn reproducibly from the census outputs.

  python3 -I certify/census/exact_sample.py SNAPDIR CENSUSDIR > exact_v4/sample.tsv

For each certificate type: the largest systems certified (isostatic: the four largest; hyperstatic: the largest) and a
random sample (Python's random with seed 2026, over the entries below a size bound that keeps the exact arithmetic to
minutes): isostatic PASS (8, at most 1,500 unknowns), hyperstatic LOCALMIN (8, at most 1,200), flexibility FLEX (8, at
most 800 rows), strictness with rational inner products (the E8 code (8,240) and 8 more, at most 1,200 unknowns). The
codes are written to CENSUSDIR/exact_v4/codes/. Output lines: mode, d, n, tolerance, code file.
"""
import glob, json, os, random, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import snapshot


def main():
    snapdir, cdir = sys.argv[1], sys.argv[2]
    idx = {(e["d"], e["n"]): e for e in snapshot.load_index()}
    res, strict = {}, {}
    for f in glob.glob(os.path.join(cdir, "results_v4", "*.json")):
        j = json.load(open(f))
        res[(j["d"], j["n"])] = j
    for f in glob.glob(os.path.join(cdir, "strict_v4", "*.json")):
        j = json.load(open(f))
        strict[(j["d"], j["n"])] = j

    def att(k, st):
        return [a for a in res[k]["attempts"] if a.get("status") == st][0]
    rng = random.Random(2026)
    pick = []
    passes = sorted((att(k, "PASS")["unknowns"], k) for k in res if res[k]["status"] == "PASS")
    big = [k for _, k in passes[-4:]]
    small = sorted(k for u, k in passes if u <= 1500 and k not in big)
    pick += [("iso", k, att(k, "PASS")["tol"]) for k in big + rng.sample(small, min(8, len(small)))]
    locs = sorted((att(k, "LOCALMIN")["unknowns"], k) for k in res if res[k]["status"] == "LOCALMIN")
    big = [k for _, k in locs[-1:]]
    small = sorted(k for u, k in locs if u <= 1200 and k not in big)
    pick += [("hyper", k, att(k, "LOCALMIN")["tol"]) for k in big + rng.sample(small, min(8, len(small)))]
    flex = sorted(k for k in res if res[k]["status"] == "FLEX" and att(k, "FLEX")["unknowns"] <= 800)
    pick += [("flex", k, att(k, "FLEX")["tol"]) for k in rng.sample(flex, min(8, len(flex)))]
    rat = sorted(k for k, s in strict.items() if s.get("status") == "STRICT" and s.get("field_degree") == 1
                 and s.get("unknowns", 10 ** 9) <= 1200 and k != (8, 240))
    sel = ([(8, 240)] if strict.get((8, 240), {}).get("status") == "STRICT" else []) + rng.sample(rat, min(8, len(rat)))
    pick += [("strict", k, strict[k]["tol"]) for k in sel]
    out = os.path.join(cdir, "exact_v4", "codes")
    os.makedirs(out, exist_ok=True)
    for mode, (d, n), tol in pick:
        path = os.path.join(out, f"d{d}_n{n}.txt")
        if not os.path.exists(path):
            with open(path, "w") as f:
                for p in snapshot.read_points(snapdir, idx[(d, n)]):
                    f.write(",".join(t.strip() for t in p) + "\n")
        print(f"{mode}\t{d}\t{n}\t{tol:g}\t{path}")


if __name__ == "__main__":
    main()
