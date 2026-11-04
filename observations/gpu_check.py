"""Largest cosine of the 26 Oct GPU-annealer test outputs (data/runs/gpu) vs the table and the CPU sweep.

  OMP_NUM_THREADS=1 python3 -I gpu_check.py > gpu_check.out          (needs numpy)
Three of the four outputs did not beat the table of 25 Oct by 1e-10 and are not in this repository (see README); their
largest cosine is read from data/measured/gpu_check.tsv unless table_codes.py has put matching copies in cache/.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
import numpy as np
from common import R, code, snapshot, lane_runs, output, listed, measured, write_measured, same, note, MEASURE
snap = snapshot()
stored = measured("gpu_check", {"maxcos": np.float64})
src, diff, out_rows = Counter(), [], []
sw = {(r["d"], r["n"]): r for r in lane_runs("sweep")}
names = set(os.listdir(os.path.join(R, "gpu"))) | {f.split("/")[-1] for f in listed("unpublished_outputs")
                                                    if f.startswith("runs/gpu/")}
for f in sorted(names):
    d, n = int(f[1:f.index("_")]), int(f.split("_n")[1].split("_")[0].split(".")[0])
    p = output(f"runs/gpu/{f}")
    if p is None:
        if MEASURE:
            sys.exit(f"--measure: runs/gpu/{f} is not in cache/ (see README)")
        gx = stored[f]["maxcos"]
        src["from data/measured"] += 1
    else:
        X = np.array(code(p, d)); X /= np.linalg.norm(X, axis=1, keepdims=True)
        G = X @ X.T; np.fill_diagonal(G, -2.0)
        gx = G.max()
        src["recomputed"] += 1
        if f in stored and not same(float(gx), stored[f]["maxcos"]):
            diff.append(f)
    out_rows.append({"key": f, "maxcos": float(gx)})
    s = sw.get((d, n))
    print(f"{f}: max cos {gx:.12f}; table {snap[(d, n)][0]}; CPU sweep {s['checked'] if s else '-'}")
if MEASURE:
    write_measured("gpu_check", ["maxcos"], out_rows)
note(f"{dict(src)}; recomputed values that differ from data/measured/gpu_check.tsv: {len(diff)} {diff[:5]}")
