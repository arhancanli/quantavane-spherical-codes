"""Look at the gap below the largest cosine in a few converged codes, to choose the contact tolerance.

  OMP_NUM_THREADS=1 python3 -I gaps.py > gaps.out          (needs numpy)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from common import lane_runs, code, output

runs = [r for r in lane_runs("finish") if r["polish"] == "converged"]
runs.sort(key=lambda r: r["d"] * r["n"])
pick = runs[::max(1, len(runs) // 12)][:12]
for r in pick:
    X = np.array(code(output(f"runs/finish/{r['file']}"), r["d"]))
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    G = X @ X.T
    iu = np.triu_indices(len(X), 1)
    g = np.sort(G[iu])[::-1]
    M = g[0]
    counts = {t: int(np.sum(g >= M - t)) for t in (1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5)}
    print(f"d{r['d']} N{r['n']}: M {M:.15f}; pairs within t of M: " + " ".join(f"{t:.0e}:{c}" for t, c in counts.items()))
