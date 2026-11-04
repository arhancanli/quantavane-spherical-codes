"""Shared readers for the observation scripts. Every path is relative to this folder (observations/); the claim records,
the table of 25 October 2026 and the saved page are read from ../records, and the uploaded codes from ../codes, as in
the rest of the repository. Nothing here writes outside data/measured/.
"""
import csv, fnmatch, hashlib, json, os, re, statistics, sys
from decimal import Decimal, getcontext

getcontext().prec = 50
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
R = os.path.join(DATA, "runs")         # run outputs (our codes), the hit-rate results and the random-start logs
C = os.path.join(HERE, os.pardir, "records")  # IMPROVED.tsv, UPLOADED.tsv, UPLOADED.log, LIVE_SEEN.tsv
TARGETS = os.path.join(C, "targets_20261025.csv")  # the table of 25 Oct 2026
LIVE_TABLE = os.path.join(C, "live_table.tsv")     # the page saved 2 Nov 2026 18:46 +04
BEST = os.path.join(DATA, "best_codes.tsv")        # our best code filed by the page time, per improved entry
RECORDS = os.path.join(DATA, "run_records.jsonl")    # each run's .done file: its text and modification time
MEASURED = os.path.join(DATA, "measured")            # values that need files this repository does not hold
CACHE = os.path.join(HERE, "cache")                  # those files, when table_codes.py has found them (not in git)
FAST_SPLIT = 1793478600  # 2026-11-01 00:30 +04: fast.sh switched from the LP polish (sphopt.py) to finish.py
PAGE = 1793630760        # 2026-11-02 18:46 +04: the paper's saved page; runs completed later are left out
MEASURE = "--measure" in sys.argv[1:]  # write data/measured/<script>.tsv (needs every file, see README)


def snapshot():
    """(d, n) -> (cos Decimal, proved, flagged) from the table of 25 Oct 2026."""
    snap = {}
    for r in csv.reader(open(TARGETS, encoding="utf-8")):
        if r[0] != "d":
            snap[(int(r[0]), int(r[1]))] = (Decimal(r[2]), r[3] == "1", r[4] == "1")
    return snap


def parse_done(text, path, mtime):
    """A lane's .done file -> dict or None (run never finished / no check line)."""
    lines = text.splitlines()
    if not lines:
        return None
    first = lines[0].split()
    if len(first) < 7 or not first[0].isdigit():
        return None
    out = {"d": int(first[0]), "n": int(first[1]), "best": Decimal(first[2]), "table": Decimal(first[3]),
           "secs": float(first[5]), "polish": first[6].replace("polish=", ""), "path": path,
           "mtime": mtime, "tier2": any(l.startswith("TIER2") for l in lines)}
    chk = [l for l in lines if l.startswith("points ")]
    yes = [l for l in lines if "improves at 12 decimals" in l]
    if not chk or not yes:
        return None
    m = re.search(r"max_cos ([0-9.]+)", chk[-1])
    out["checked"] = Decimal(m.group(1))
    out["yes"] = yes[-1].strip().endswith("YES")
    return out


_RECS = None


def _records():
    """lane -> [(file name, modification time, text)], from data/run_records.jsonl (git keeps no file times)."""
    global _RECS
    if _RECS is None:
        _RECS = {}
        for line in open(RECORDS, encoding="utf-8"):
            x = json.loads(line)
            lane, name = x["path"].split("/")
            _RECS.setdefault(lane, []).append((name, x["mtime"], x["text"]))
    return _RECS


def lane_runs(lane, include_cut=False, until=PAGE):
    """The lane's completed runs (in the order of a sorted glob of the lane's *.txt.done[.cut*] files)."""
    pats = ["*.txt.done"] + (["*.txt.done.cut*"] if include_cut else [])
    recs = _records().get(lane, [])
    runs = []
    for pat in pats:
        for name, mtime, text in sorted((x for x in recs if fnmatch.fnmatchcase(x[0], pat)), key=lambda x: x[0]):
            r = parse_done(text, os.path.join(R, lane, name), mtime)
            if r and (until is None or r["mtime"] <= until):
                r["file"] = name.split(".txt")[0] + ".txt"
                runs.append(r)
    return runs


def code(path, d=None):
    pts = []
    for line in open(path, errors="replace"):
        if line.startswith("#") or not line.strip():
            continue
        try:
            v = [float(x) for x in line.strip().split(",")]
        except ValueError:
            continue
        if d is None or len(v) == d:
            pts.append(v)
    return pts


def q(xs, p):
    xs = sorted(xs)
    if not xs:
        return float("nan")
    k = (len(xs) - 1) * p
    lo, hi = int(k), min(int(k) + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def med(xs):
    return statistics.median(xs) if xs else float("nan")


# --- files this repository lists but does not hold -------------------------------------------------------------------
# data/table_codes.tsv: the table's codes of 25-26 Oct 2026 that the scripts compare with (Henry Cohn's data, not
# redistributed here). data/unpublished_outputs.tsv: our run outputs that did not beat the table of 25 Oct by 1e-10,
# which are the table's codes unchanged or polished in the last digits. table_codes.py puts copies with the listed
# SHA-256 into cache/; without them the scripts take the values computed from those files from data/measured/.

def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def listed(name):
    """data/<name>.tsv -> {file: row dict}."""
    rows = list(csv.DictReader(open(os.path.join(DATA, name + ".tsv"), encoding="utf-8"), delimiter="\t"))
    return {r["file"]: r for r in rows}


_LISTS = {}


def _cached(kind, rel):
    if kind not in _LISTS:
        _LISTS[kind] = listed(kind)
    row = _LISTS[kind].get(rel)
    p = os.path.join(CACHE, "table" if kind == "table_codes" else "outputs", rel)
    return p if row and os.path.exists(p) and _sha256(p) == row["sha256"] else None


def is_listed(kind, rel):
    """rel is in data/<kind>.tsv (kind: table_codes or unpublished_outputs)."""
    if kind not in _LISTS:
        _LISTS[kind] = listed(kind)
    return rel in _LISTS[kind]


def table_code(folder, name):
    """Path of the table's code of 25-26 Oct (folder needs_opt or unflag), or None when cache/ has no matching copy."""
    return _cached("table_codes", folder + "/" + name)


_BEST = None


def best_code(name):
    """Path of our best code filed by the page time for this entry (data/best_codes.tsv: ../codes/<name> when it is the
    uploaded code, else data/best/<name>), or None when that file is missing or its SHA-256 is not the listed one."""
    global _BEST
    if _BEST is None:
        _BEST = {r["entry"] + ".txt": r for r in csv.DictReader(open(BEST, encoding="utf-8"), delimiter="\t")}
    row = _BEST.get(name)
    p = os.path.normpath(os.path.join(HERE, row["file"])) if row else None
    return p if p and os.path.exists(p) and _sha256(p) == row["sha256"] else None


def output(rel):
    """Path of one of our codes (rel: runs/<lane>/<file>, claims/<file> or bench/s<seed>/<file>): data/<rel> when
    published, else a matching copy in cache/outputs/<rel>, else None."""
    p = os.path.join(DATA, rel)
    return p if os.path.exists(p) else _cached("unpublished_outputs", rel)


def measured(script, types):
    """data/measured/<script>.tsv -> {key: row}; types maps a column to its type (others stay str)."""
    path = os.path.join(MEASURED, script + ".tsv")
    if MEASURE or not os.path.exists(path):
        return {}
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        out[r["key"]] = {k: (types[k](v) if k in types and v != "" else v) for k, v in r.items()}
    return out


def write_measured(script, cols, rows):
    """--measure: write data/measured/<script>.tsv (rows: dicts with a 'key' and the columns), floats as repr."""
    os.makedirs(MEASURED, exist_ok=True)
    with open(os.path.join(MEASURED, script + ".tsv"), "w", encoding="utf-8") as f:
        f.write("\t".join(["key"] + cols) + "\n")
        for r in rows:
            f.write("\t".join([r["key"]] + [repr(float(r[c])) if isinstance(r[c], float) else str(r[c])
                                             for c in cols]) + "\n")


def same(a, b):
    """Recomputed value equals the measured one (floats compared exactly; nan equals nan)."""
    if isinstance(a, float) and isinstance(b, float):
        return a == b or (a != a and b != b)
    return a == b


def note(msg):
    """Where the values came from: printed to stderr, so that stdout stays the paper's output."""
    print(os.path.basename(sys.argv[0]) + ": " + msg, file=sys.stderr)
