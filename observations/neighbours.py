"""Entries whose 25 Oct value nearly equals the value for one more point (a code that may be the larger code with a
point removed): are they more open to improvement?

  python3 -I neighbours.py > neighbours.out
"""
import csv, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datetime as dt
from decimal import Decimal
from common import C, R, snapshot, lane_runs, med

TZ = dt.timezone(dt.timedelta(hours=4))
PAGE = dt.datetime(2026, 11, 2, 18, 46, tzinfo=TZ).timestamp()
snap = snapshot()


def ts(s):
    m = re.match(r"(\d{4}-\d\d-\d\d) (\d\d):(\d\d)(?::(\d\d))?", s.strip())
    y, mo, d = (int(x) for x in m.group(1).split("-"))
    return dt.datetime(y, mo, d, int(m.group(2)), int(m.group(3)), int(m.group(4) or 0), tzinfo=TZ).timestamp()


tried = set()
for lane in ["lowtest", "sweep", "sweep3", "fast"]:
    for r in lane_runs(lane):
        tried.add((r["d"], r["n"]))
best = {}
for line in open(os.path.join(C, "IMPROVED.tsv")).read().splitlines()[1:]:
    p = line.split("\t")
    if len(p) >= 5 and p[0].isdigit() and ts(p[4]) <= PAGE:
        k = (int(p[0]), int(p[1]))
        best[k] = min(best.get(k, Decimal(p[3])), Decimal(p[3]))
ft = [k for k in tried if snap[k][2]]
for tol in (Decimal("1e-8"), Decimal("1e-6")):
    near = [k for k in ft if (k[0], k[1] + 1) in snap and snap[(k[0], k[1] + 1)][0] - snap[k][0] <= tol]
    rest = [k for k in ft if k not in near]
    for name, s in (("value within %s of N+1's" % tol, near), ("others", rest)):
        imp = [k for k in s if k in best and snap[k][0] - best[k] >= Decimal("1e-10")]
        g = [float(snap[k][0] - best[k]) for k in imp]
        print(f"flagged tried, {name}: {len(s)}, improved {len(imp)} ({100 * len(imp) / max(1, len(s)):.1f}%), "
              f"median gain {med(g):.2e}, gain >= 1e-4: {sum(1 for x in g if x >= 1e-4)}")
