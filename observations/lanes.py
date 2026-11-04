"""Per-lane yield: completed runs, how many beat the 25 Oct table at 12 decimals and by >= 1e-10, and the gains.

  python3 -I lanes.py > lanes.out
A run counts as completed when its .done file has the independent checker's lines. The gain is the table value of
25 Oct minus the checker's largest cosine of the run's output code.
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import defaultdict
from decimal import Decimal
from common import R, snapshot, lane_runs, med, q, FAST_SPLIT

snap = snapshot()
M10 = Decimal("1e-10")


def summary(name, runs):
    if not runs:
        print(f"{name:16s} no completed runs")
        return
    fl = sum(1 for r in runs if snap[(r["d"], r["n"])][2])
    gains = [float(r["table"] - r["checked"]) for r in runs]
    hit12 = [r for r in runs if r["yes"]]
    hit10 = [g for g in gains if g >= 1e-10]
    conv = sum(1 for r in runs if r["polish"] == "converged")
    ds = sorted({r["d"] for r in runs})
    print(f"{name:16s} runs {len(runs):5d} (flagged {fl:5d}) d {ds[0]}-{ds[-1]}  beats@12dec {len(hit12):5d} "
          f"({100 * len(hit12) / len(runs):5.1f}%)  gain>=1e-10 {len(hit10):5d} ({100 * len(hit10) / len(runs):5.1f}%)  "
          f"median gain of those {med(hit10):.2e} [q25 {q(hit10, .25):.2e}, q75 {q(hit10, .75):.2e}]  "
          f"polish converged {conv}  median secs {med([r['secs'] for r in runs]):.0f}")


# the hit-rate test of 25-26 Oct (results.txt: d n best table improved secs; no checker lines there)
hr = []
for line in open(os.path.join(R, "hitrate", "results.txt")):
    p = line.split()
    if len(p) >= 6 and p[0].isdigit():
        hr.append((int(p[0]), int(p[1]), Decimal(p[2]), Decimal(p[3]), int(p[4]), float(p[5])))
hits = [x for x in hr if x[3] - x[2] >= M10]
print(f"hitrate test: targets 16, completed {len(hr)}, gain >= 1e-10: {len(hits)}; "
      f"by d: " + ", ".join(f"d{x[0]} N{x[1]} {float(x[3] - x[2]):.2e}" for x in hr))
print(f"   completed with d >= 18: {sum(1 for x in hr if x[0] >= 18)}, of which improved {sum(1 for x in hits if x[0] >= 18)};"
      f" d <= 17: {sum(1 for x in hr if x[0] <= 17)}, improved {sum(1 for x in hits if x[0] <= 17)}")

lanes = {}
for lane in ["abtest", "lowtest", "sweep", "sweep3", "sweep4", "deep", "repolish", "converge", "unflag", "uall"]:
    lanes[lane] = lane_runs(lane)
fast = lane_runs("fast")
lanes["fast-LP(31Oct)"] = [r for r in fast if r["mtime"] < FAST_SPLIT]
lanes["fast-finish(1-2Nov)"] = [r for r in fast if r["mtime"] >= FAST_SPLIT]
lanes["finish(converged)"] = lane_runs("finish")
lanes["finish(cut)"] = [r for r in lane_runs("finish", include_cut=True) if ".cut" in r["path"]]
for k, v in lanes.items():
    summary(k, v)

# the start-from-the-table lanes split by dimension range and by flag
print("\nby population (runs that start from the table's own code):")
for name in ["sweep", "sweep3", "fast-LP(31Oct)", "fast-finish(1-2Nov)", "unflag", "uall", "lowtest"]:
    runs = lanes[name]
    for lo, hi in [(3, 10), (11, 17), (18, 24), (25, 32)]:
        sub = [r for r in runs if lo <= r["d"] <= hi]
        if sub:
            g = [float(r["table"] - r["checked"]) for r in sub]
            print(f"  {name:20s} d {lo:2d}-{hi:2d}: runs {len(sub):4d}, gain>=1e-10 {sum(1 for x in g if x >= 1e-10):4d}"
                  f" ({100 * sum(1 for x in g if x >= 1e-10) / len(sub):5.1f}%)")

# fast-finish vs table: how many of its runs converged, and its hits by flag
ff = lanes["fast-finish(1-2Nov)"]
print(f"\nfast-finish: flagged runs {sum(1 for r in ff if snap[(r['d'], r['n'])][2])}, unflagged {sum(1 for r in ff if not snap[(r['d'], r['n'])][2])}")
fl_tier2 = sum(1 for r in fast if r["tier2"])
print(f"fast (all): TIER2 marks (no gain >= 1e-10) {fl_tier2} of {len(fast)}")
