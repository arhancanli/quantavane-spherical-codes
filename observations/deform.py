"""Is each improved code a small deformation of the table's code it started from, or a different configuration?

  OMP_NUM_THREADS=1 python3 -I deform.py > deform.out          (needs numpy)
For every entry with a filed gain >= 1e-10 by the paper's page time: our best code filed by then (data/best_codes.tsv:
the uploaded code in ../codes, or data/best/dD_nN.txt) against the table's code we started from (downloaded on 25-26 Oct for flagged
entries, 31 Oct for the others), same point order, after the best rotation. The per-point move is divided by the
minimal distance of our code, sqrt(2 - 2 max cos). The median over points is used, so that a few free points
(rattlers), which wander during the polish, do not decide it.
The table's codes are not in this repository (see README): the moves are read from data/measured/deform.tsv, unless
table_codes.py has put matching copies in cache/, when they are recomputed and checked against it. --measure writes
that file (every file needed).
"""
import math, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datetime as dt
from collections import Counter
from decimal import Decimal
import numpy as np
from common import C, snapshot, code, med, q, best_code, table_code, measured, write_measured, same, note, MEASURE

TZ = dt.timezone(dt.timedelta(hours=4))
PAGE = dt.datetime(2026, 11, 2, 18, 46, tzinfo=TZ).timestamp()
snap = snapshot()


def ts(s):
    m = re.match(r"(\d{4}-\d\d-\d\d) (\d\d):(\d\d)(?::(\d\d))?", s.strip())
    y, mo, d = (int(x) for x in m.group(1).split("-"))
    return dt.datetime(y, mo, d, int(m.group(2)), int(m.group(3)), int(m.group(4) or 0), tzinfo=TZ).timestamp()


best = {}
for line in open(os.path.join(C, "IMPROVED.tsv")).read().splitlines()[1:]:
    p = line.split("\t")
    if len(p) >= 5 and p[0].isdigit() and ts(p[4]) <= PAGE:
        k = (int(p[0]), int(p[1]))
        best[k] = min(best.get(k, Decimal(p[3])), Decimal(p[3]))
keys = sorted(k for k, v in best.items() if snap[k][0] - v >= Decimal("1e-10"))
COLS = ["M", "median", "q90", "max"]
stored = measured("deform", {c: float for c in COLS})
src, diff, out_rows = Counter(), [], []
rows, miss = [], Counter()
for k in keys:
    d, n = k
    f = f"d{d}_n{n}.txt"
    key = f"{d}_{n}"
    start = table_code("needs_opt" if snap[k][2] else "unflag", f)
    ours = best_code(f)
    s = stored.get(key)
    if MEASURE and (start is None or ours is None):
        sys.exit(f"--measure: the table's code or our best code {f} is missing (see README)")
    M = stats = None
    if ours is not None:
        X = np.array(code(ours, d))
        if X.shape != (n, d):
            miss["shape"] += 1
            continue
        X /= np.linalg.norm(X, axis=1, keepdims=True)
        G = X @ X.T; np.fill_diagonal(G, -2.0); M = G.max()
        src["code recomputed"] += 1
        if start is not None:
            Y = np.array(code(start, d))
            if Y.shape != (n, d):
                miss["shape"] += 1
                continue
            Y /= np.linalg.norm(Y, axis=1, keepdims=True)
            U, _, Vt = np.linalg.svd(Y.T @ X)
            mv = np.linalg.norm(Y @ (U @ Vt) - X, axis=1) / math.sqrt(2 - 2 * M)
            stats = (float(np.median(mv)), float(np.quantile(mv, 0.9)), float(mv.max()))
            src["move recomputed"] += 1
    if s is not None:
        if M is None:
            M = s["M"]
            src["code from data/measured"] += 1
        elif not same(float(M), s["M"]):
            diff.append((key, "M"))
        if stats is None:
            stats = (s["median"], s["q90"], s["max"])
            src["move from data/measured"] += 1
        elif not all(same(a, b) for a, b in zip(stats, (s["median"], s["q90"], s["max"]))):
            diff.append((key, "move"))
    if M is None or stats is None:
        miss["missing"] += 1
        continue
    if abs(float(M) - float(best[k])) > 1e-13:
        miss["code is not the best filed code"] += 1
    out_rows.append({"key": key, "M": float(M), "median": stats[0], "q90": stats[1], "max": stats[2]})
    rows.append((k, float(snap[k][0]) - float(M), stats[0], stats[1], stats[2]))
if MEASURE:
    write_measured("deform", COLS, out_rows)
note(f"{dict(src)}; recomputed values that differ from data/measured/deform.tsv: {len(diff)} {diff[:5]}")
print(f"improved entries compared: {len(rows)} of {len(keys)} {dict(miss)}")
medmv = [r[2] for r in rows]
print(f"median per-point move / minimal distance: median {med(medmv):.1e}, q90 {q(medmv, .9):.1e}, q99 {q(medmv, .99):.1e},"
      f" max {max(medmv):.2f}")
for t in (0.01, 0.05, 0.1, 0.25, 0.5):
    print(f"  entries whose median point moved by more than {t} of the minimal distance: {sum(1 for x in medmv if x > t)}")
for t in (1e-4, 1e-5, 1e-6, 1e-7):
    s = [r for r in rows if r[1] >= t]
    print(f"  gain >= {t:.0e}: {len(s)} entries, median point move > 0.1 in {sum(1 for r in s if r[2] > 0.1)}")
s = [r for r in rows if r[1] < 1e-5]
print(f"  gain < 1e-5: {len(s)} entries, median point move > 0.1 in {sum(1 for r in s if r[2] > 0.1)}")
# rank correlation between gain and move
def rank(x):
    o = sorted(range(len(x)), key=lambda i: x[i]); r = [0] * len(x)
    for i, j in enumerate(o):
        r[j] = i
    return r
a, b = rank([r[1] for r in rows]), rank(medmv)
nn = len(rows)
print(f"Spearman rank correlation of gain and median move: {1 - 6 * sum((x - y) ** 2 for x, y in zip(a, b)) / (nn * (nn * nn - 1)):.2f}")
far = sorted(rows, key=lambda r: -r[2])[:15]
print("largest median moves (entry, gain, median move, q90, max):", [(r[0], f"{r[1]:.1e}", f"{r[2]:.2f}", f"{r[3]:.2f}", f"{r[4]:.2f}") for r in far])
print("by d of the entries with median move > 0.1:", dict(Counter(r[0][0] for r in rows if r[2] > 0.1)))
print("d8 entries with median move > 0.1 of the minimal distance, N:", sorted(r[0][1] for r in rows if r[2] > 0.1 and r[0][0] == 8))
print("entries with median move > 0.1, gains:", sorted(f"{r[1]:.1e}" for r in rows if r[2] > 0.1)[:5], "...")
