"""Entries run both by an anneal-then-LP-polish lane (sweep, sweep3, hitrate, lowtest) and later by finish.py straight
from the table's code with no annealing (fast lane after 1 Nov 00:30): same values? same local optimum?

  OMP_NUM_THREADS=1 python3 -I same_basin.py > same_basin.out          (needs numpy)
Same point order in both outputs (both descend from the same downloaded table code), so the two codes are compared
point by point after the best rotation, over the points that are not rattlers in the finish.py code.
An output that did not beat the table of 25 Oct by 1e-10 is the table's code unchanged or polished in the last digits
and is not in this repository (see README): for a pair with such a code the move is read from
data/measured/same_basin.tsv, unless table_codes.py has put a matching copy in cache/. --measure writes that file.
"""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter
import numpy as np
from common import lane_runs, code, med, q, FAST_SPLIT, output, measured, write_measured, same, note, MEASURE

TOL = 1e-10
stored = measured("same_basin", {"mv": float})
src, diff, out_rows = Counter(), [], []
ff = {(r["d"], r["n"]): r for r in lane_runs("fast") if r["mtime"] >= FAST_SPLIT}
rows = []
for lane in ("sweep", "sweep3", "lowtest"):
    for r in lane_runs(lane):
        k = (r["d"], r["n"])
        if k not in ff:
            continue
        f = ff[k]
        key = f"{lane}_{k[0]}_{k[1]}"
        pa, pb = output(f"runs/{lane}/{r['file']}"), output(f"runs/fast/{f['file']}")
        if pa is None or pb is None:
            if MEASURE:
                sys.exit(f"--measure: runs/{lane}/{r['file']} or runs/fast/{f['file']} is missing (see README)")
            if key not in stored:
                continue
            mv = stored[key]["mv"]
            src["move from data/measured"] += 1
        else:
            A = np.array(code(pa, k[0])); B = np.array(code(pb, k[0]))
            if A.shape != (k[1], k[0]) or B.shape != A.shape:
                continue
            A /= np.linalg.norm(A, axis=1, keepdims=True); B /= np.linalg.norm(B, axis=1, keepdims=True)
            G = B @ B.T; np.fill_diagonal(G, -2.0); M = G.max()
            C = G >= M - TOL
            alive = np.ones(k[1], bool)
            while True:
                deg = (C & alive[None, :] & alive[:, None]).sum(1)
                drop = alive & (deg < k[0])
                if not drop.any():
                    break
                alive &= ~drop
            if alive.sum() < k[0]:
                alive[:] = True  # no rigid core at this tolerance (an unconverged run): compare all points
            U, _, Vt = np.linalg.svd(A[alive].T @ B[alive])
            mv = float(np.linalg.norm(A[alive] @ (U @ Vt) - B[alive], axis=1).max()) / math.sqrt(2 - 2 * M)
            src["move recomputed"] += 1
            if key in stored and not same(mv, stored[key]["mv"]):
                diff.append(key)
        out_rows.append({"key": key, "mv": mv})
        rows.append((lane, k, float(r["checked"]), float(f["checked"]), float(r["table"]), mv, f["polish"]))
if MEASURE:
    write_measured("same_basin", ["mv"], out_rows)
note(f"{dict(src)}; recomputed moves that differ from data/measured/same_basin.tsv: {len(diff)} {diff[:5]}")
print(f"entries run by an annealing lane and by finish.py from the table's code: {len(rows)}")
for lane in ("sweep", "sweep3", "lowtest"):
    s = [x for x in rows if x[0] == lane]
    if not s:
        continue
    diff = [x[2] - x[3] for x in s]  # > 0: finish.py from the table code is better
    print(f"  {lane}: {len(s)}; finish.py-from-table better by >= 1e-10 in {sum(1 for x in diff if x >= 1e-10)}, worse by"
          f" >= 1e-10 in {sum(1 for x in diff if x <= -1e-10)}, within 1e-10 in {sum(1 for x in diff if abs(x) < 1e-10)}")
    anneal_gain = [x for x in s if x[4] - x[2] >= 1e-10]
    print(f"    of the {len(anneal_gain)} where the annealing lane gained >= 1e-10: finish.py from the table code reached at "
          f"least as low (within 1e-10) in {sum(1 for x in anneal_gain if x[3] <= x[2] + 1e-10)}")
    mv = [x[5] for x in s]
    print(f"    non-rattler max move between the two codes / minimal distance: median {med(mv):.1e}, q90 {q(mv, .9):.1e},"
          f" max {max(mv):.2e}; above 0.1: {sum(1 for x in mv if x > 0.1)}")
    conv = [x for x in s if x[6] == "converged"]
    print(f"    where finish.py converged ({len(conv)}): values within 1e-10 of the annealing lane's "
          f"{sum(1 for x in conv if abs(x[2] - x[3]) < 1e-10)}, finish.py lower by >= 1e-10 {sum(1 for x in conv if x[2] - x[3] >= 1e-10)}"
          f", higher by >= 1e-10 {sum(1 for x in conv if x[3] - x[2] >= 1e-10)}")
    worse = sorted([x for x in s if x[3] - x[2] >= 1e-10], key=lambda x: x[2] - x[3])[:8]
    print("    finish.py worse (lane, entry, anneal-lane value, finish.py value, move):",
          [(x[1], f"{x[2]:.12f}", f"{x[3]:.12f}", f"{x[5]:.2f}", x[6]) for x in worse])
print("\nwhere finish.py from the table's code converged but the annealing lane's code is lower by >= 1e-10:")
for x in sorted(rows, key=lambda x: x[1]):
    if x[6] == "converged" and x[3] - x[2] >= 1e-10:
        print(f"  d{x[1][0]} N{x[1][1]}: table {x[4]:.12f}, anneal+LP {x[2]:.12f}, finish.py from table {x[3]:.12f}, "
              f"difference {x[3] - x[2]:.2e}, non-rattler move between them / minimal distance {x[5]:.2f}")
big = [x for x in rows if x[4] - x[2] >= 2.5e-4]
print(f"entries where the annealing lane gained >= 2.5e-4 over the table: {len(big)}; of them finish.py from the table "
      f"code reached within 1e-10 or lower: {sum(1 for x in big if x[3] <= x[2] + 1e-10)}")
