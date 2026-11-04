"""Every number in paper.tex, from our claim records and the last saved copy of the table page (all in this folder).

  python3 -I make_table.py                         writes numbers.tex, table.tex, examples.tex and headline.json
  python3 -I make_table.py --live-page PAGE.html   first rebuilds live_table.tsv from a saved copy of the front page
                                                   of https://spherical-codes.org (one plain GET), then as above

Inputs (read only, this folder):
  targets_20261025.csv  the table on 25 Oct 2026 (d, n, cos, proved_optimal, needs_opt, ref, coords_link)
  IMPROVED.tsv          every checked code that beat that table (d, n, table_cos, our_cos, time filed); rows filed
                        after the page was saved are left out, so the numbers are those of the saved page's time
  UPLOADED.tsv          every upload the live table served just after the post (code, cos)
  BEATEN.tsv            entries others improved past ours (optional)
  LIVE_SEEN.tsv         live values the uploader saw when it held back (optional, no header)
  live_table.tsv        the table as displayed when the page was saved (d, n, displayed cosine, credit); its first
                        line gives the time it was saved, which dates the numbers
  before_upload.tsv (this folder)  for the entries where the uploader kept its copy of the live code it replaced, that
                                   code's largest cosine (written once, on 31 Oct 2026, from the copies autoupload.py
                                   and attended.sh had saved; optional)

Counting.
  records rule: distinct codes in UPLOADED.tsv, minus BEATEN.tsv, minus those whose latest LIVE_SEEN value is below our
    latest uploaded cosine by more than 1e-14 (the site re-serialises coordinates, which moves a cosine by about 1e-17).
  live rule (the count in the paper): of those, the entries whose displayed cosine on the saved page equals one of our
    uploaded cosines rounded up to 12 decimals. The site displays cosines rounded up: on 31 Oct 2026 this held for all
    205 entries it credited to us. A displayed value below our code means a better code was posted after ours; an entry
    missing from the page is dominated by a code with more points; a displayed value above every code of ours means our
    upload came after the page was saved (not counted until the page is saved again).
Gains run from the value the table held before our first upload to our latest uploaded code: in cosine, before - ours;
in angle, arccos(ours) - arccos(before), in degrees. "Before" is the replaced code's value from before_upload.tsv when it
is below the 25 Oct value (someone improved the entry in between), and otherwise the 25 Oct value.
"""
import csv
import html
import json
import math
import os
import re
import statistics
import sys
import time
from collections import Counter, defaultdict
from decimal import ROUND_CEILING, Decimal, InvalidOperation, getcontext

getcontext().prec = 50
HERE = os.path.dirname(os.path.abspath(__file__))
CLAIMS = HERE
SNAPSHOT = os.path.join(HERE, "targets_20261025.csv")
LIVE = os.path.join(HERE, "live_table.tsv")
BEFORE = os.path.join(HERE, "before_upload.tsv")
SEEN_TOL = Decimal("1e-14")
SLACK = Decimal("1e-15")
ULP12 = Decimal("1e-12")
CHANGED = Decimal("1e-11")  # a replaced code this far below the 25 Oct value means the entry had changed meanwhile
MONTHS = "January February March April May June July August September October November December".split()


def key_of(name):
    m = re.fullmatch(r"d(\d+)_n(\d+)", name.strip())
    return (int(m.group(1)), int(m.group(2))) if m else None


def number(s):
    try:
        return Decimal(s.strip())
    except (InvalidOperation, AttributeError):
        return None


def rows(path):
    """Tab-separated rows of a file, [] if it is missing; header and malformed lines are left to the caller."""
    if not os.path.exists(path):
        return []
    return [line.rstrip("\n").split("\t") for line in open(path, encoding="utf-8") if line.strip()]


def keyed(path):
    """[(key, value)] from lines 'dD_nN<TAB>value...', skipping headers and malformed lines, in file order."""
    out = []
    for r in rows(path):
        k, v = key_of(r[0]), number(r[1]) if len(r) > 1 else None
        if k and v is not None:
            out.append((k, v))
    return out


def extract_live(page):
    """live_table.tsv from a saved copy of the table's front page: one line per entry (d, n, displayed cos, credit)."""
    text = open(page, encoding="utf-8", errors="replace").read()
    out = []
    for tr in re.findall(r"<tr>(.*?)</tr>", text, re.S):
        cells = [re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", c))).strip()
                 for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)]
        if len(cells) < 4 or not cells[0].isdigit():
            continue
        m = re.fullmatch(r"(\d+)\*?", cells[1])
        if not m or not re.fullmatch(r"-?\d\.\d+", cells[2]):
            continue
        out.append((cells[0], m.group(1), cells[2], cells[3].replace("\t", " ")[:200]))
    saved = time.strftime("%Y-%m-%d %H:%M %z", time.localtime(os.path.getmtime(page)))
    with open(LIVE, "w", encoding="utf-8") as f:
        f.write(f"# saved\t{saved}\tthe table at https://spherical-codes.org as displayed; d, n, cosine, credit\n")
        f.writelines("\t".join(r) + "\n" for r in out)
    print(f"{len(out)} entries -> live_table.tsv (page saved {saved})")


def angle_gain(before, ours):
    """arccos(ours) - arccos(before) in degrees, from the sine of the difference computed in 50 digits."""
    one = Decimal(1)
    sin_diff = before * (one - ours * ours).sqrt() - ours * (one - before * before).sqrt()
    return math.degrees(math.asin(float(sin_diff)))


def sci(x, digits=2):
    """x as $m\\cdot10^{e}$ (LaTeX), m with `digits` significant digits."""
    if x == 0:
        return "$0$"
    e = math.floor(math.log10(abs(x)))
    m = round(x / 10 ** e, digits - 1)
    if abs(m) >= 10:
        m, e = m / 10, e + 1
    return f"${m:.{digits - 1}f}\\cdot10^{{{e}}}$"


def plain(x, digits=2):
    """x with `digits` significant digits, as a plain decimal between 0.001 and 1000, else as sci()."""
    if not 1e-3 <= abs(x) < 1e3:
        return sci(x, digits)
    places = max(0, digits - 1 - math.floor(math.log10(abs(x))))
    return f"{x:.{places}f}"


def grouped(n):
    return f"{n:,}".replace(",", "{,}")


def long_date(stamp):
    y, mo, d = (int(v) for v in stamp[:10].split("-"))
    return f"{d} {MONTHS[mo - 1]} {y}"


def period(start, end):
    """'25 to 31 October 2026', '25 October to 3 November 2026' or with both years."""
    a, b = long_date(start).split(), long_date(end).split()
    if a[2] != b[2]:
        return f"{' '.join(a)} to {' '.join(b)}"
    if a[1] != b[1]:
        return f"{a[0]} {a[1]} to {' '.join(b)}"
    return f"{a[0]} to {' '.join(b)}"


def ceil12(x):
    """A cosine as the table displays it: rounded up to 12 decimals."""
    return x.quantize(Decimal("1e-12"), rounding=ROUND_CEILING)


def page_saved():
    """'YYYY-MM-DD HH:MM +zzzz' from the first line of live_table.tsv: the numbers are to this time."""
    assert os.path.exists(LIVE), "no live_table.tsv: save the table's front page and run with --live-page PAGE.html"
    first = open(LIVE, encoding="utf-8").readline().rstrip("\n").split("\t")
    assert first[0] == "# saved", "live_table.tsv does not start with its '# saved' line"
    return first[1]


def main():
    if "--live-page" in sys.argv:
        extract_live(sys.argv[sys.argv.index("--live-page") + 1])
    cut = page_saved()[:16]  # filings up to the minute the page was saved

    # the table of 25 Oct 2026: (cos, proved optimal, needs more optimization)
    snap = {}
    for r in csv.reader(open(SNAPSHOT, encoding="utf-8")):
        if r[0] != "d":
            snap[(int(r[0]), int(r[1]))] = (Decimal(r[2]), r[3] == "1", r[4] == "1")

    # every checked code that beat it; the table value it was compared with must be the snapshot's
    found = {}
    for r in rows(os.path.join(CLAIMS, "IMPROVED.tsv")):
        if len(r) < 4 or not r[0].isdigit() or (len(r) > 4 and r[4][:16] > cut):
            continue
        k, table, ours = (int(r[0]), int(r[1])), number(r[2]), number(r[3])
        assert k in snap and table == snap[k][0], f"IMPROVED.tsv {k}: table value {table} is not the snapshot's"
        found[k] = min(found.get(k, ours), ours)

    # uploads confirmed live, in file order (later lines are re-uploads of better codes), and the records rule
    uploaded = defaultdict(list)
    for k, v in keyed(os.path.join(CLAIMS, "UPLOADED.tsv")):
        uploaded[k].append(v)
    assert all(k in snap for k in uploaded), "an uploaded entry is missing from the snapshot"
    margin = Decimal("1e-10")  # the paper says every upload gained at least this much; check it
    for k, vals in uploaded.items():
        assert snap[k][0] - vals[0] >= margin, f"{k}: first upload gains less than 1e-10 over the 25 Oct table"
        # re-uploads (31 Oct 22:xx rule): an entry the table credits to us is uploaded again when the displayed value
        # (rounded up to 12 decimals) improves, so each re-upload must lower the displayed value
        assert all(ceil12(b) < ceil12(a) for a, b in zip(vals, vals[1:])), f"{k}: a re-upload does not lower the displayed value"
    beaten = {k for k, _ in keyed(os.path.join(CLAIMS, "BEATEN.tsv"))}
    seen = dict(keyed(os.path.join(CLAIMS, "LIVE_SEEN.tsv")))  # the latest line for each entry wins
    records = {k: v for k, v in uploaded.items()
               if k not in beaten and not (k in seen and seen[k] < v[-1] - SEEN_TOL)}

    # the live rule, against the saved page
    assert os.path.exists(LIVE), "no live_table.tsv: save the table's front page and run with --live-page PAGE.html"
    live, saved = {}, None
    for r in rows(LIVE):
        if r[0].startswith("#"):
            saved = r[1]
        else:
            live[(int(r[0]), int(r[1]))] = (Decimal(r[2]), r[3] if len(r) > 3 else "")
    status, credit, maint_gap = {}, Counter(), []
    for k, vals in uploaded.items():  # every upload is classified; only records can count as ours
        if k not in live:
            status[k] = "dominated"
            continue
        shown, who = live[k]
        if any(-SLACK <= shown - v < ULP12 + SLACK for v in vals):
            status[k] = "ours" if k in records else "ours on the page but not by the records"
            if k in records:
                credit["us" if "Canli" in who else "no name" if "anonymous" in who else "other"] += 1
        elif shown < vals[-1] - SLACK:
            status[k] = "maintainer" if "needs more optimization" in who else "others"
            if status[k] == "maintainer":
                maint_gap.append(float(vals[-1] - shown))
        else:
            status[k] = "pending"
    counts = Counter(status.values())
    ours = sorted(k for k in records if status[k] == "ours")
    assert ours, "no entry is ours on the saved page"
    credited = [k for k in uploaded if k in live and "Canli" in live[k][1]]
    if any(not any(-SLACK <= live[k][0] - v < ULP12 + SLACK for v in uploaded[k]) for k in credited):
        print("WARNING: an entry credited to us does not display our cosine rounded up; the paper says it always does",
              file=sys.stderr)
    unexplained = len(uploaded) - sum(counts[s] for s in ("ours", "maintainer", "others", "dominated", "pending"))
    if unexplained:  # the paper's breakdown of the uploads would not add up: read the status counts before publishing
        print(f"WARNING: {unexplained} uploads fall outside the paper's breakdown: {dict(counts)}", file=sys.stderr)

    # gains, from the value before our first upload
    copy = {key_of(r[0]): number(r[2]) for r in rows(BEFORE)  # code, copy_saved, live_cos, saved_by
            if key_of(r[0]) and len(r) > 2 and number(r[2]) is not None}
    gains, changed = {}, []
    for k in ours:
        before, after = snap[k][0], records[k][-1]
        if k in copy and copy[k] < before - CHANGED:
            before = copy[k]
            changed.append(k)
        assert before > after, f"{k}: our code {after} is not below the value before it, {before}"
        gains[k] = (float(before - after), angle_gain(before, after), before, after)
    with_copy = [k for k in ours if k in copy]

    # table.tex: one row per dimension with an upload
    by_d = defaultdict(list)
    for k in ours:
        by_d[k[0]].append(gains[k])
    up_d = Counter(k[0] for k in uploaded)

    def stats(g):
        return (f"{sci(max(x[0] for x in g))} & {sci(max(x[1] for x in g))} & "
                f"{sci(statistics.median(x[0] for x in g))} & {sci(statistics.median(x[1] for x in g))}")
    lines = [r"\begin{tabular}{rrrllll}", r"\toprule",
             r" & \multicolumn{2}{c}{entries} & \multicolumn{2}{c}{largest gain} & \multicolumn{2}{c}{median gain} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
             r"$d$ & uploaded & ours now & cosine & angle ($^\circ$) & cosine & angle ($^\circ$) \\", r"\midrule"]
    for d in sorted(up_d):
        g = by_d.get(d, [])
        rest = stats(g) if g else "-- & -- & -- & --"
        lines.append(f"{d} & {grouped(up_d[d])} & {grouped(len(g))} & {rest} \\\\")
    allg = list(gains.values())
    lines += [r"\midrule", f"all & {grouped(len(uploaded))} & {grouped(len(ours))} & {stats(allg)} \\\\",
              r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(HERE, "table.tex"), "w").write("\n".join(lines) + "\n")

    # examples.tex: among entries whose replaced code we kept (else all of ours), the largest angular gain, the largest
    # in the lowest dimension and the largest in the highest dimension
    pool = with_copy if len(with_copy) >= 3 else ours
    best = sorted(pool, key=lambda k: -gains[k][1])
    dims = sorted({k[0] for k in pool})
    lowest = max((k for k in pool if k[0] == dims[0]), key=lambda k: gains[k][1])
    highest = max((k for k in pool if k[0] == dims[-1]), key=lambda k: gains[k][1])
    picks = []
    for k in [best[0], lowest, highest] + best[1:]:
        if k not in picks:
            picks.append(k)
        if len(picks) == 3:
            break
    picks.sort()
    ex = [r"\begin{tabular}{rrllll}", r"\toprule",
          r" & & \multicolumn{2}{c}{largest cosine} & \multicolumn{2}{c}{minimal angle ($^\circ$)} \\",
          r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
          r"$d$ & $N$ & before & ours & before & ours \\", r"\midrule"]
    for k in picks:
        before, after = gains[k][2], gains[k][3]
        places = max(6, 2 - math.floor(math.log10(gains[k][1])))
        ex.append(f"{k[0]} & {k[1]} & {before.quantize(ULP12, rounding=ROUND_CEILING)} & "
                  f"{after.quantize(ULP12, rounding=ROUND_CEILING)} & "
                  f"{math.degrees(math.acos(float(before))):.{places}f} & "
                  f"{math.degrees(math.acos(float(after))):.{places}f} \\\\")
    ex += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(HERE, "examples.tex"), "w").write("\n".join(ex) + "\n")

    # numbers.tex (macros used by paper.tex) and headline.json
    kc = max(ours, key=lambda k: gains[k][0])
    ka = max(ours, key=lambda k: gains[k][1])
    flagged = sum(1 for k in ours if snap[k][2])
    nums = {
        "DataDate": long_date(saved),
        "Period": period("2026-10-25", saved[:10]),
        "LiveDate": long_date(saved),
        "LiveTime": f"{saved[11:16]} (UTC{saved[17:20]}:{saved[20:22]})",
        "NSnap": grouped(len(snap)),
        "NFlagged": grouped(sum(1 for v in snap.values() if v[2])),
        "NProved": grouped(sum(1 for v in snap.values() if v[1])),
        "NFound": grouped(len(found)),
        "NFoundMargin": grouped(sum(1 for k, v in found.items() if snap[k][0] - v >= Decimal("1e-10"))),
        "NUploaded": grouped(len(uploaded)),
        "NRecords": grouped(len(records)),
        "NLive": grouped(len(ours)),
        "NLiveFlagged": grouped(flagged),
        "FlaggedPhrase": (f"All {grouped(len(ours))} entries were" if flagged == len(ours)
                          else f"{grouped(flagged)} of the {grouped(len(ours))} entries were"),
        "NLiveCredited": grouped(credit["us"]),
        "NLiveCopy": grouped(len(with_copy)),
        "NLiveNoCopy": grouped(len(ours) - len(with_copy)),
        "NChanged": grouped(len(changed)),
        "ChangedDims": ", ".join(str(d) for d in sorted({k[0] for k in changed})) or "none",
        "NBig": grouped(sum(1 for k in ours if gains[k][0] >= 1e-6)),
        "NSuperseded": grouped(counts["maintainer"] + counts["others"] + counts["dominated"]),
        "NMaintainer": grouped(counts["maintainer"]),
        "MedMaintGap": sci(statistics.median(maint_gap)) if maint_gap else "--",
        "NOthers": grouped(counts["others"]),
        "NDominated": grouped(counts["dominated"]),
        "NPending": grouped(counts["pending"]),
        "NDims": str(len(by_d)),
        "DMin": str(min(by_d)),
        "DMax": str(max(by_d)),
        "MaxDcos": sci(gains[kc][0]), "MaxDcosD": str(kc[0]), "MaxDcosN": str(kc[1]),
        "MaxDang": plain(gains[ka][1]), "MaxDangD": str(ka[0]), "MaxDangN": str(ka[1]),
        "MedDcos": sci(statistics.median(x[0] for x in allg)),
        "MedDang": sci(statistics.median(x[1] for x in allg)),
    }
    with open(os.path.join(HERE, "numbers.tex"), "w") as f:
        f.write("% written by make_table.py; do not edit\n")
        f.writelines(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in nums.items())
    head = {
        "computed": saved,  # the numbers are to the time the page was saved
        "live_page_saved": saved,
        "snapshot_entries": len(snap),
        "found_better_than_25oct": len(found),
        "uploaded_confirmed": len(uploaded),
        "count_records_rule": len(records),
        "count_live_rule": len(ours),
        "records_rule_breakdown_on_page": dict(counts),
        "ours_on_page_by_credit": dict(credit),
        "median_gap_to_maintainer_codes": statistics.median(maint_gap) if maint_gap else None,
        "dimensions": {"min": min(by_d), "max": max(by_d), "count": len(by_d)},
        "before_value": {"replaced_code_kept": len(with_copy), "changed_since_25oct": [f"d{d}_n{n}" for d, n in changed]},
        "largest_cos_gain": {"d": kc[0], "n": kc[1], "gain": gains[kc][0]},
        "largest_angle_gain_deg": {"d": ka[0], "n": ka[1], "gain": gains[ka][1]},
        "median_cos_gain": statistics.median(x[0] for x in allg),
        "median_angle_gain_deg": statistics.median(x[1] for x in allg),
        "examples": [{"d": k[0], "n": k[1], "before": str(gains[k][2]), "ours": str(gains[k][3])} for k in picks],
        "per_dimension": {d: {"uploaded": up_d[d], "ours_now": len(by_d.get(d, []))} for d in sorted(up_d)},
    }
    json.dump(head, open(os.path.join(HERE, "headline.json"), "w"), indent=1)
    print(f"uploaded {len(uploaded)}; records rule {len(records)}; ours on the page saved {saved}: {len(ours)} "
          f"({dict(counts)}); dimensions {min(by_d)}-{max(by_d)}; before value changed since 25 Oct: {len(changed)}")


if __name__ == "__main__":
    main()
