"""Which program produced each filed improvement and each upload, and how much the strict finisher moved values.

  python3 -I history.py > history.out
Sources: ../records/IMPROVED.tsv (every filing), UPLOADED.tsv + UPLOADED.log (every confirmed upload, in order),
the lanes' .done files in data/run_records.jsonl (each run's summary line and the independent checker's 20-digit
largest cosine; a filed or uploaded value is attributed to the run whose checked value is the same string) and
../records/live_table.tsv, the paper's saved page. Program of a run: finish.py for the finish lane after 1 Nov 01:35
(+04) with secs < 1700 and for the fast lane after 1 Nov 00:30; sphopt.py (LP polish, optionally after annealing)
otherwise.
"""
import datetime as dt
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collections import Counter, defaultdict
from decimal import Decimal
from common import C, LIVE_TABLE, R, snapshot, lane_runs, med, q, FAST_SPLIT

TZ = dt.timezone(dt.timedelta(hours=4))
FIN_SPLIT = dt.datetime(2026, 11, 1, 1, 35, tzinfo=TZ).timestamp()
GATE = dt.datetime(2026, 10, 31, 21, 48, tzinfo=TZ).timestamp()      # the uploader's gate (engine commit fcf130e): converged codes only
PAGE = dt.datetime(2026, 11, 2, 18, 46, tzinfo=TZ).timestamp()       # the paper's saved page
snap = snapshot()


def ts(s):
    m = re.match(r"(\d{4}-\d\d-\d\d) (\d\d):(\d\d)(?::(\d\d))?", s.strip())
    if not m:
        return None
    y, mo, d = (int(x) for x in m.group(1).split("-"))
    return dt.datetime(y, mo, d, int(m.group(2)), int(m.group(3)), int(m.group(4) or 0), tzinfo=TZ).timestamp()


# every lane run (completed ones, cut ones included) indexed by (d, n, checked value)
runs = []
for lane in ["hitrate", "lowtest", "sweep", "sweep3", "sweep4", "deep", "repolish", "converge", "unflag", "uall",
             "fast", "finish"]:
    for r in lane_runs(lane, include_cut=True):
        prog = "sphopt"
        if lane == "fast" and r["mtime"] >= FAST_SPLIT:
            prog = "finish.py"
        if lane == "finish" and r["mtime"] >= FIN_SPLIT and r["secs"] < 1700:
            prog = "finish.py"
        r["lane"], r["prog"] = lane, prog
        runs.append(r)
by_val = defaultdict(list)
for r in runs:
    by_val[(r["d"], r["n"], str(r["checked"]))].append(r)


def source(d, n, val, when=None):
    """The run that produced this checked value (the earliest one, at or before `when` + 5 min if given)."""
    cands = by_val.get((d, n, str(val)), [])
    if when is not None:
        cands = [r for r in cands if r["mtime"] <= when + 300] or cands
    return min(cands, key=lambda r: r["mtime"]) if cands else None


# 1. filings
filings = []
for line in open(os.path.join(C, "IMPROVED.tsv")).read().splitlines()[1:]:
    p = line.split("\t")
    if len(p) < 5 or not p[0].isdigit():
        continue
    d, n, table, ours = int(p[0]), int(p[1]), Decimal(p[2]), Decimal(p[3])
    filings.append((ts(p[4]), d, n, table, ours, p[4]))
print(f"filings in IMPROVED.tsv: {len(filings)}, distinct entries {len({(f[1], f[2]) for f in filings})}")
pre = [f for f in filings if f[0] is not None and f[0] <= PAGE]
ent_pre = {}
for f in pre:
    k = (f[1], f[2])
    ent_pre[k] = min(ent_pre.get(k, f[4]), f[4])
print(f"  filed by the paper's page time (2 Nov 18:46): {len(ent_pre)} entries; with gain >= 1e-10 over 25 Oct: "
      f"{sum(1 for k, v in ent_pre.items() if snap[k][0] - v >= Decimal('1e-10'))}")
ent_all = {}
for f in filings:
    k = (f[1], f[2])
    ent_all[k] = min(ent_all.get(k, f[4]), f[4])
print(f"  all filings (to the end of the runs): {len(ent_all)} entries; gain >= 1e-10: "
      f"{sum(1 for k, v in ent_all.items() if snap[k][0] - v >= Decimal('1e-10'))}")

# first filing per entry with gain >= 1e-10: which lane/program found it
first = {}
for f in sorted(pre, key=lambda f: f[0]):
    k = (f[1], f[2])
    if k in first or snap[k][0] - f[4] < Decimal("1e-10"):
        continue
    s = source(f[1], f[2], f[4], f[0])
    first[k] = (s["lane"] + ("/" + s["prog"] if s["lane"] in ("fast", "finish") else "")) if s else "unmatched:" + f[5][17:]
print("\nfirst filing with gain >= 1e-10 (to 2 Nov 18:46), by lane:", dict(Counter(first.values()).most_common()))

# 2. uploads, in order, with their time and the run that produced the uploaded code
up_vals = []
for line in open(os.path.join(C, "UPLOADED.tsv")).read().splitlines()[1:]:
    a, b = line.split("\t")[:2]
    up_vals.append((a, Decimal(b)))
# times from UPLOADED.log, matched to each UPLOADED.tsv row by entry and value (the log line carries the live value,
# equal to ours within 1e-14); rows with no such line are the hand uploads of 26 Oct (logged without values)
times = defaultdict(list)
for line in open(os.path.join(C, "UPLOADED.log")):
    m = re.match(r"# (?:attended curl POST|autoupload) (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d): (d\d+_n\d+) live (?:= ours )?([0-9.]+)", line)
    if m:
        times[m.group(2)].append((ts(m.group(1)), Decimal(m.group(3))))
HAND = dt.datetime(2026, 10, 26, 12, 0, tzinfo=TZ).timestamp()
ups, unmatched_time = [], 0
for name, v in up_vals:
    cand = [t for t, x in times[name] if abs(x - v) < Decimal("1e-14")]
    t = min(cand) if cand else None
    if t is None:
        unmatched_time += 1
        t = HAND
    d, n = (int(x) for x in re.findall(r"\d+", name))
    s = source(d, n, v, t)
    ups.append((name, d, n, v, t, s))
print(f"\nupload rows with no timed log line (hand uploads of 26 Oct, given 26 Oct 12:00): {unmatched_time}")
after = [u for u in ups if u[4] is not None and u[4] >= GATE]
print(f"\nuploads (lines of UPLOADED.tsv): {len(ups)}; after the converged-only gate (31 Oct 21:48): {len(after)}")
kinds = Counter()
for u in after:
    s = u[5]
    if s is None:
        kinds["no matching run"] += 1
    else:
        kinds[f"{s['prog']} {s['polish']}"] += 1
print("  program and polish status of the run that produced each uploaded code after the gate:", dict(kinds))
# the gate accepts an entry when ANY converged run of it has the same cosine to 1e-13; count by that rule too
conv = defaultdict(list)
for r in runs:
    if r["polish"] == "converged":
        conv[(r["d"], r["n"])].append(r)
k2 = Counter()
for u in after:
    c = [r for r in conv[(u[1], u[2])] if abs(r["best"] - u[3]) <= Decimal("1e-13")]
    k2["finish.py" if any(r["prog"] == "finish.py" for r in c) else "sphopt only" if c else "none"] += 1
print("  converged run within 1e-13 of the uploaded cosine (the uploader's gate):", dict(k2))
before = [u for u in ups if u[4] is not None and u[4] < GATE]
kb = Counter((u[5]["prog"] + " " + u[5]["polish"]) if u[5] else "no matching run" for u in before)
print(f"  before the gate: {len(before)} uploads; producing run:", dict(kb))

# 3. how much finish.py moved values: for each entry, best filed value before its first finish.py run vs the best
# value after its finish.py runs (filings tagged 'finish' or from the fast lane's finish.py runs)
fin_runs = [r for r in runs if r["prog"] == "finish.py"]
gain_conv, gain_any, n_conv = [], [], 0
for k in {(r["d"], r["n"]) for r in fin_runs if r["lane"] == "finish"}:
    rs = sorted([r for r in fin_runs if (r["d"], r["n"]) == k and r["lane"] == "finish"], key=lambda r: r["mtime"])
    t0 = rs[0]["mtime"] - rs[0]["secs"] * 4 - 600  # the run started at most 4x its CPU budget earlier
    prior = [f[4] for f in filings if (f[1], f[2]) == k and f[0] is not None and f[0] < t0]
    if not prior:
        continue
    b0 = min(prior)
    b1 = min(r["checked"] for r in rs)
    g = float(b0 - b1)
    gain_any.append(g)
    if any(r["polish"] == "converged" for r in rs):
        n_conv += 1
        gain_conv.append(g)
print(f"\nfinish lane (finish.py) on entries we had already filed: {len(gain_any)} entries; converged {n_conv}")
for name, g in (("all", gain_any), ("converged", gain_conv)):
    print(f"  gain over our best earlier code ({name}): median {med(g):.2e}, q25 {q(g, .25):.2e}, q75 {q(g, .75):.2e}, "
          f"max {max(g):.2e}; >= 1e-8: {sum(1 for x in g if x >= 1e-8)}, >= 1e-10: {sum(1 for x in g if x >= 1e-10)},"
          f" < 1e-12: {sum(1 for x in g if x < 1e-12)}")

# 4. the paper's breakdown: superseded uploads vs the program that produced the last code we uploaded
live = {}
for line in open(LIVE_TABLE, encoding="utf-8"):
    if line.startswith("#"):
        continue
    p = line.rstrip("\n").split("\t")
    live[(int(p[0]), int(p[1]))] = (Decimal(p[2]), p[3] if len(p) > 3 else "")
last = {}
for u in ups:
    last[(u[1], u[2])] = u
stat = Counter()
for k, u in last.items():
    if k not in live:
        st = "gone"
    else:
        shown, who = live[k]
        if "Canli" in who:
            st = "ours"
        elif "needs more optimization" in who:
            st = "maintainer"
        else:
            st = "other"
    s = u[5]
    conv_any = any(abs(r["best"] - u[3]) <= Decimal("1e-13") for r in conv[k])
    conv_fin = any(abs(r["best"] - u[3]) <= Decimal("1e-13") and r["prog"] == "finish.py" for r in conv[k])
    tag = "finish.py-converged" if conv_fin else "sphopt-converged" if conv_any else "not converged"
    stat[(st, tag)] += 1
print("\nlast uploaded code of each entry, by its status on the 2 Nov 18:46 page and by convergence of that code:")
for st in ("ours", "maintainer", "other", "gone"):
    row = {tag: stat[(st, tag)] for tag in ("finish.py-converged", "sphopt-converged", "not converged")}
    print(f"  {st:10s} {row}")
for tag in ("finish.py-converged", "sphopt-converged", "not converged"):
    tot = sum(stat[(st, tag)] for st in ("ours", "maintainer", "other", "gone"))
    if tot:
        print(f"  {tag:20s}: {tot} entries, still ours {stat[('ours', tag)]} ({100 * stat[('ours', tag)] / tot:.0f}%)")

# 5. exposure control: survival to the 2 Nov page by the date of the entry's last upload and its convergence class
print("\nsurvival by date of the last upload (ours on the 2 Nov 18:46 page / entries):")
tab = defaultdict(Counter)
for k, u in last.items():
    if u[4] is None:
        continue
    day = dt.datetime.fromtimestamp(u[4], TZ).strftime("%m-%d")
    conv_any = any(abs(r["best"] - u[3]) <= Decimal("1e-13") for r in conv[k])
    conv_fin = any(abs(r["best"] - u[3]) <= Decimal("1e-13") and r["prog"] == "finish.py" for r in conv[k])
    tag = "finish.py" if conv_fin else "sphopt-conv" if conv_any else "not-conv"
    ours = k in live and "Canli" in live[k][1]
    tab[(day, tag)]["n"] += 1
    tab[(day, tag)]["ours"] += int(ours)
for day in sorted({x[0] for x in tab}):
    print("  " + day + ": " + "  ".join(f"{tag} {tab[(day, tag)]['ours']}/{tab[(day, tag)]['n']}"
                                    for tag in ("not-conv", "sphopt-conv", "finish.py") if tab[(day, tag)]["n"]))

# 6. held back by the gate: entries with a filed gain >= 1e-10 (to the page time) that were never uploaded
upk = {(u[1], u[2]) for u in ups}
found10 = {k for k, v in ent_pre.items() if snap[k][0] - v >= Decimal("1e-10")}
never = found10 - upk
hb = Counter()
for k in never:
    v = ent_pre[k]
    c = [r for r in conv[k] if abs(r["best"] - v) <= Decimal("1e-13")]
    anyconv = bool(conv[k])
    hb["best code converged" if c else "an older code converged" if anyconv else "never converged"] += 1
print(f"\nentries improved by >= 1e-10 by the page time: {len(found10)}; uploaded {len(found10 & upk)}; never uploaded "
      f"{len(never)}: {dict(hb)}")
seen_live = {}
for line in open(os.path.join(C, "LIVE_SEEN.tsv")):
    p = line.split("\t")
    if len(p) >= 2 and p[0].startswith("d"):
        seen_live[p[0]] = p[1]
nb = sum(1 for k in never if f"d{k[0]}_n{k[1]}" in seen_live)
print(f"  of those, with a live value seen better than ours (LIVE_SEEN.tsv): {nb}")
# the uploads after the gate whose code had no converged run within 1e-13 (the uploader's own rule)
odd = []
for u in after:
    c = [r for r in conv[(u[1], u[2])] if abs(r["best"] - u[3]) <= Decimal("1e-13")]
    if not c:
        cl = sorted(conv[(u[1], u[2])], key=lambda r: r["best"])
        odd.append((u[0], str(u[3])[:16], dt.datetime.fromtimestamp(u[4], TZ).strftime("%m-%d %H:%M"),
                    [(r["lane"], str(r["best"])[:16]) for r in cl[:2]], u[5]["lane"] + "/" + u[5]["polish"] if u[5] else None))
print(f"  uploads after the gate with no converged run within 1e-13 now: {len(odd)}")
for o in odd:
    print("   ", o)

# 7. the never-uploaded entries by the time their first >= 1e-10 code was filed
STOP = dt.datetime(2026, 11, 2, 7, 45, tzinfo=TZ).timestamp()  # the uploader's stop (AUTOUPLOAD_STOP)
firstt = {}
for f in sorted(pre, key=lambda f: f[0]):
    k = (f[1], f[2])
    if k not in firstt and snap[k][0] - f[4] >= Decimal("1e-10"):
        firstt[k] = f[0]
bins = Counter()
for k in never:
    t = firstt[k]
    conv_best = any(abs(r["best"] - ent_pre[k]) <= Decimal("1e-13") for r in conv[k])
    w = "before the gate" if t < GATE else "gate to uploader stop" if t < STOP else "after the stop"
    bins[(w, "converged" if conv_best else "not converged")] += 1
print("never uploaded, by when first found and convergence of our best code:", dict(bins))
found_window = [k for k in found10 if GATE <= firstt[k] < STOP]
print(f"entries first found (>= 1e-10) between the gate and the stop: {len(found_window)}, uploaded "
      f"{sum(1 for k in found_window if k in upk)}")
held = [k for k in found_window if k not in upk and not any(abs(r["best"] - ent_pre[k]) <= Decimal("1e-13") for r in conv[k])]
print(f"held back (first found in the window, never uploaded, best code not converged): {len(held)}; of them with a "
      f"live value seen better than ours: {sum(1 for k in held if f'd{k[0]}_n{k[1]}' in seen_live)}")
