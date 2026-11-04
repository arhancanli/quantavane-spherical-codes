"""Which entries were open to improvement: by dimension, by size, by minimal angle; the gains.

  python3 -I coverage.py > coverage.out
"Tried": a completed run (checker lines present) that starts from the table's own code: lanes hitrate, lowtest,
sweep, sweep3, fast (both programs), unflag, uall. "Improved": some filed code of ours (IMPROVED.tsv, to the paper's
page time 2 Nov 18:46) beats the 25 Oct value by >= 1e-10. "Ours": the 447 entries the paper counts (live_table.tsv
credit, same rule as make_table.py: credited to Canli and displaying one of our uploaded cosines rounded up).
"""
import datetime as dt
import math, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter, defaultdict
from decimal import Decimal, ROUND_CEILING
from common import C, LIVE_TABLE, R, snapshot, lane_runs, med, q

TZ = dt.timezone(dt.timedelta(hours=4))
PAGE = dt.datetime(2026, 11, 2, 18, 46, tzinfo=TZ).timestamp()
snap = snapshot()
M10 = Decimal("1e-10")


def ts(s):
    m = re.match(r"(\d{4}-\d\d-\d\d) (\d\d):(\d\d)(?::(\d\d))?", s.strip())
    y, mo, d = (int(x) for x in m.group(1).split("-"))
    return dt.datetime(y, mo, d, int(m.group(2)), int(m.group(3)), int(m.group(4) or 0), tzinfo=TZ).timestamp()


tried = set()
for x in open(os.path.join(R, "hitrate", "results.txt")):
    p = x.split()
    if len(p) >= 6 and p[0].isdigit():
        tried.add((int(p[0]), int(p[1])))
for lane in ["lowtest", "sweep", "sweep3", "fast", "unflag", "uall"]:
    for r in lane_runs(lane):
        tried.add((r["d"], r["n"]))
best = {}
for line in open(os.path.join(C, "IMPROVED.tsv")).read().splitlines()[1:]:
    p = line.split("\t")
    if len(p) >= 5 and p[0].isdigit() and ts(p[4]) <= PAGE:
        k = (int(p[0]), int(p[1]))
        best[k] = min(best.get(k, Decimal(p[3])), Decimal(p[3]))
improved = {k for k, v in best.items() if snap[k][0] - v >= M10}
print(f"tried {len(tried)} (flagged {sum(1 for k in tried if snap[k][2])}); improved {len(improved)} "
      f"(flagged {sum(1 for k in improved if snap[k][2])}); improved but never tried by a table-start lane: "
      f"{len(improved - tried)}")

# the paper's 447 (same rule as make_table.py)
uploaded = defaultdict(list)
for line in open(os.path.join(C, "UPLOADED.tsv")).read().splitlines()[1:]:
    a, b = line.split("\t")[:2]
    d, n = (int(x) for x in re.findall(r"\d+", a))
    uploaded[(d, n)].append(Decimal(b))
live = {}
for line in open(LIVE_TABLE, encoding="utf-8"):
    if not line.startswith("#"):
        p = line.rstrip("\n").split("\t")
        live[(int(p[0]), int(p[1]))] = (Decimal(p[2]), p[3] if len(p) > 3 else "")
ours = {k for k, v in uploaded.items() if k in live and "Canli" in live[k][1]
        and any(-Decimal("1e-15") <= live[k][0] - x < Decimal("1e-12") + Decimal("1e-15") for x in v)}
print(f"ours on the 2 Nov page: {len(ours)}")

# by dimension
print("\nd | table entries (proved, flagged) | flagged tried, improved (rate) | unflagged tried, improved | ours "
      "| improved / table entries")
rows = []
for d in range(3, 33):
    ent = [k for k in snap if k[0] == d]
    fl = [k for k in ent if snap[k][2]]
    pr = [k for k in ent if snap[k][1]]
    ft = [k for k in fl if k in tried]
    fi = [k for k in ft if k in improved]
    ut = [k for k in ent if k in tried and not snap[k][2]]
    ui = [k for k in ut if k in improved]
    oc = sum(1 for k in ours if k[0] == d)
    imp_all = sum(1 for k in improved if k[0] == d)
    rows.append((d, len(ent), len(fl), len(ft), len(fi), oc))
    print(f"{d:2d} | {len(ent):4d} ({len(pr):2d}, {len(fl):4d} = {100 * len(fl) / max(1, len(ent)):4.1f}%) | "
          f"{len(ft):4d}, {len(fi):4d} ({100 * len(fi) / max(1, len(ft)):5.1f}%) | {len(ut):3d}, {len(ui):3d} | "
          f"{oc:3d} | {imp_all:4d} / {len(ent):4d} = {100 * imp_all / max(1, len(ent)):4.1f}%")
nmax = Counter()
for k in snap:
    nmax[k[0]] = max(nmax[k[0]], k[1])
print("largest N listed per d:", dict(nmax))

# rank correlation across d between the flagged share and our improved share of the table's entries
def spearman(a, b):
    def rank(x):
        o = sorted(range(len(x)), key=lambda i: x[i])
        r = [0] * len(x)
        for i, j in enumerate(o):
            r[j] = i
        return r
    ra, rb = rank(a), rank(b)
    n = len(a)
    return 1 - 6 * sum((x - y) ** 2 for x, y in zip(ra, rb)) / (n * (n * n - 1))

# flagged rate among tried entries, overall, by d range, by size
ft = [k for k in tried if snap[k][2]]
print(f"\nflagged tried {len(ft)}, improved {sum(1 for k in ft if k in improved)} "
      f"({100 * sum(1 for k in ft if k in improved) / len(ft):.1f}%)")
for lo, hi in [(3, 5), (6, 10), (11, 17), (18, 24), (25, 32)]:
    s = [k for k in ft if lo <= k[0] <= hi]
    a = [k for k in snap if lo <= k[0] <= hi and snap[k][2]]
    if s:
        print(f"  d {lo:2d}-{hi:2d}: flagged in table {len(a):4d}, tried {len(s):4d}, improved "
              f"{sum(1 for k in s if k in improved):4d} ({100 * sum(1 for k in s if k in improved) / len(s):5.1f}%)")
print("by N/d (flagged tried entries):")
for lo, hi in [(2, 4), (4, 8), (8, 16), (16, 32), (32, 64), (64, 400)]:
    s = [k for k in ft if lo < k[1] / k[0] <= hi]
    a = [k for k in snap if lo < k[1] / k[0] <= hi and snap[k][2]]
    if s:
        print(f"  N/d in ({lo},{hi}]: flagged in table {len(a):4d}, tried {len(s):4d}, improved "
              f"{sum(1 for k in s if k in improved):4d} ({100 * sum(1 for k in s if k in improved) / len(s):5.1f}%)")
print("by minimal angle of the table's code (flagged tried entries):")
for name, f in [("theta > 60 (cos < 1/2)", lambda c: c < Decimal("0.5")), ("theta <= 60", lambda c: c >= Decimal("0.5"))]:
    s = [k for k in ft if f(snap[k][0])]
    a = [k for k in snap if f(snap[k][0]) and snap[k][2]]
    print(f"  {name}: flagged in table {len(a)}, tried {len(s)}, improved {sum(1 for k in s if k in improved)} "
          f"({100 * sum(1 for k in s if k in improved) / max(1, len(s)):.1f}%)")
print("unflagged tried:", len([k for k in tried if not snap[k][2]]), "improved:",
      len([k for k in tried if not snap[k][2] and k in improved]),
      sorted((k, float(snap[k][0] - best[k])) for k in tried if not snap[k][2] and k in improved))

# where the table's flags sit, by d (share of a dimension's entries flagged) vs where we improved
fl_share = [r[2] / r[1] for r in rows]
imp_share = [sum(1 for k in improved if k[0] == r[0]) / r[1] for r in rows]
print(f"\nSpearman rank correlation over d = 3..32 of (flagged share, improved share): {spearman(fl_share, imp_share):.2f}")

# gains
def angle_gain(before, after):
    one = Decimal(1)
    s = before * (one - after * after).sqrt() - after * (one - before * before).sqrt()
    return math.degrees(math.asin(float(s)))

for name, keys in (("improved (gain >= 1e-10, to 2 Nov 18:46)", improved),):
    g = [float(snap[k][0] - best[k]) for k in keys]
    a = [angle_gain(snap[k][0], best[k]) for k in keys]
    print(f"\ngains, {name}: n={len(g)}")
    print(f"  cosine: median {med(g):.2e}, q10 {q(g, .1):.2e}, q25 {q(g, .25):.2e}, q75 {q(g, .75):.2e}, "
          f"q90 {q(g, .9):.2e}, max {max(g):.2e}")
    for t in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
        print(f"    below {t:.0e}: {sum(1 for x in g if x < t)} ({100 * sum(1 for x in g if x < t) / len(g):.1f}%)")
    print(f"  angle (degrees): median {med(a):.2e}, q25 {q(a, .25):.2e}, q75 {q(a, .75):.2e}, max {max(a):.3f}")
    for t in (1e-8, 1e-6, 1e-4, 1e-3, 1e-2):
        print(f"    below {t:.0e} deg: {sum(1 for x in a if x < t)} ({100 * sum(1 for x in a if x < t) / len(a):.1f}%)")
    for lo, hi in [(3, 10), (11, 17), (18, 24), (25, 32)]:
        s = [float(snap[k][0] - best[k]) for k in keys if lo <= k[0] <= hi]
        print(f"    d {lo}-{hi}: n={len(s)} median cosine gain {med(s):.2e}, >= 1e-6: {sum(1 for x in s if x >= 1e-6)}")
    for lo, hi in [(2, 4), (4, 8), (8, 16), (16, 32), (32, 64), (64, 400)]:
        s = [float(snap[k][0] - best[k]) for k in keys if lo < k[1] / k[0] <= hi]
        if s:
            print(f"    N/d ({lo},{hi}]: n={len(s)} median cosine gain {med(s):.2e}")
    big = sorted(keys, key=lambda k: -float(snap[k][0] - best[k]))[:12]
    print("  largest cosine gains:", [(k, f"{float(snap[k][0] - best[k]):.2e}") for k in big])

# how the improvements relate to the size of the table in each dimension
ent_d = [r[1] for r in rows]
print(f"\nSpearman over d = 3..32: (entries in the table, improved share) {spearman(ent_d, imp_share):.2f}; "
      f"(entries, flagged share) {spearman(ent_d, fl_share):.2f}")
lo_d = [r for r in rows if r[0] <= 6]
print(f"d 3-6: entries {sum(r[1] for r in lo_d)}, flagged {sum(r[2] for r in lo_d)}, improved "
      f"{sum(1 for k in improved if k[0] <= 6)}")
mid = [r for r in rows if 7 <= r[0] <= 21]
print(f"d 7-21: entries {sum(r[1] for r in mid)}, flagged {sum(r[2] for r in mid)} ({100 * sum(r[2] for r in mid) / sum(r[1] for r in mid):.0f}%),"
      f" improved {sum(1 for k in improved if 7 <= k[0] <= 21)}")
hi = [r for r in rows if r[0] >= 22]
print(f"d 22-32: entries {sum(r[1] for r in hi)}, flagged {sum(r[2] for r in hi)} ({100 * sum(r[2] for r in hi) / sum(r[1] for r in hi):.0f}%),"
      f" improved {sum(1 for k in improved if k[0] >= 22)}")
print(f"unflagged tried by d range: d<=17 {sum(1 for k in tried if not snap[k][2] and k[0] <= 17)} (improved "
      f"{sum(1 for k in tried if not snap[k][2] and k[0] <= 17 and k in improved)}), d>=18 "
      f"{sum(1 for k in tried if not snap[k][2] and k[0] >= 18)} (improved {sum(1 for k in tried if not snap[k][2] and k[0] >= 18 and k in improved)})")
