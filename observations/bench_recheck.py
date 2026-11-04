"""Recompute, from the saved codes, the 31 Oct annealing benchmark (8 flagged entries, data/bench/picks.txt; seeds 31-33,
Python annealer 30/60 s and C annealer 12/30/60 s; outputs bench/s3x/{py30,py60,c12,c30,c60}_dD_nN.txt): is each
annealing output the start code (the table's code), and what is its largest cosine?

  OMP_NUM_THREADS=1 python3 -I bench_recheck.py > bench_recheck.out          (needs numpy)
Every one of the 120 outputs is the table's code of 25-26 Oct (to 2.2e-16 per coordinate), and neither those outputs
nor the table's codes are in this repository (see README; data/unpublished_outputs.tsv and data/table_codes.tsv list
them with their SHA-256). The values are read from data/measured/bench_recheck.tsv, unless table_codes.py has put
matching copies in cache/, when they are recomputed and checked against it. --measure writes that file.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
import numpy as np
from common import DATA, code, output, table_code, measured, write_measured, same, note, MEASURE

stored = measured("bench_recheck", {"maxcos": np.float64, "dev": float})
src, diff, out_rows = Counter(), [], []
picks = [l.split() for l in open(os.path.join(DATA, "bench", "picks.txt")) if l.strip()]
tot = same_n = 0
for d, n, table, _, _ in picks:
    d, n = int(d), int(n)
    start = table_code("needs_opt", f"d{d}_n{n}.txt")
    if MEASURE and start is None:
        sys.exit(f"--measure: the table's code needs_opt/d{d}_n{n}.txt is not in cache/ (see README)")
    key0 = f"start_{d}_{n}"
    if start is not None:
        Y = np.array(code(start, d))
        Y /= np.linalg.norm(Y, axis=1, keepdims=True)
        G = Y @ Y.T; np.fill_diagonal(G, -2.0); m0 = G.max()
        src["start recomputed"] += 1
        if key0 in stored and not same(float(m0), stored[key0]["maxcos"]):
            diff.append(key0)
    else:
        Y, m0 = None, stored[key0]["maxcos"]
        src["start from data/measured"] += 1
    out_rows.append({"key": key0, "maxcos": float(m0), "dev": ""})
    worst = 0.0
    for s in (31, 32, 33):
        for k in ("py30", "py60", "c12", "c30", "c60"):
            key = f"s{s}_{k}_{d}_{n}"
            p = output(f"bench/s{s}/{k}_d{d}_n{n}.txt")
            if p is None or Y is None:
                if MEASURE:
                    sys.exit(f"--measure: bench/s{s}/{k}_d{d}_n{n}.txt is not in cache/ (see README)")
                if key not in stored:
                    continue
                dev, gx = stored[key]["dev"], stored[key]["maxcos"]
                src["output from data/measured"] += 1
            else:
                X = np.array(code(p, d))
                if X.shape != (n, d):
                    continue
                X /= np.linalg.norm(X, axis=1, keepdims=True)
                G = X @ X.T; np.fill_diagonal(G, -2.0)
                dev, gx = float(np.abs(X - Y).max()), G.max()
                src["output recomputed"] += 1
                if key in stored and not (same(dev, stored[key]["dev"]) and same(float(gx), stored[key]["maxcos"])):
                    diff.append(key)
            out_rows.append({"key": key, "maxcos": float(gx), "dev": dev})
            tot += 1
            same_n += int(dev < 1e-12 and abs(gx - m0) < 1e-15)
            worst = max(worst, dev)
    print(f"d{d} N{n}: start (table code) max cos {m0:.15f}; largest coordinate difference of any annealing output from "
          f"the start {worst:.1e}")
print(f"annealing outputs found {tot}; identical to the start code (coordinates to 1e-12, max cos to 1e-15): {same_n}")
if MEASURE:
    write_measured("bench_recheck", ["maxcos", "dev"], out_rows)
note(f"{dict(src)}; recomputed values that differ from data/measured/bench_recheck.tsv: {len(diff)} {diff[:5]}")
