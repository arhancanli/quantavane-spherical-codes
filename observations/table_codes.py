"""Find copies of the files this folder lists but does not hold, and put them in cache/ (not in git):
  - data/table_codes.tsv: the codes of Henry Cohn's table (spherical-codes.org) as we downloaded them on 25-26 Oct 2026
    (marked entries) and 31 Oct 2026 (unmarked entries), with each file's URL and SHA-256;
  - data/unpublished_outputs.tsv: our run outputs that did not beat the table of 25 Oct by 1e-10 (the table's codes
    unchanged or polished in the last digits), with each file's SHA-256.
A file is kept only when its SHA-256 is the listed one, so the scripts then recompute exactly what the paper used.

  python3 -I table_codes.py fetch [dD_nN ...]   download listed table codes from spherical-codes.org (one per second)
  python3 -I table_codes.py copy DIR [DIR ...]  copy matching files found anywhere under DIR
  python3 -I table_codes.py status              how many listed files cache/ holds

The table serves only the current code of each entry, so an entry improved since 25 Oct (by us or by the table's own
polishing) no longer matches; the archive at https://hdl.handle.net/1721.1/153543 is a snapshot of 31 July 2026.
"""
import hashlib, os, sys, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import CACHE, listed

UA = "quantavane-spherical-codes observations/table_codes.py (https://github.com/arhancanli/quantavane-spherical-codes)"


def sha(b):
    return hashlib.sha256(b).hexdigest()


def put(kind, rel, data):
    p = os.path.join(CACHE, "table" if kind == "table_codes" else "outputs", rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "wb") as f:
        f.write(data)


def fetch(names):
    rows = listed("table_codes")
    todo = [r for r in rows.values() if not names or r["file"].split("/")[1][:-4] in names]
    got = changed = failed = 0
    for i, r in enumerate(todo):
        if i:
            time.sleep(1)
        try:
            req = urllib.request.Request(r["url"], headers={"User-Agent": UA})
            data = urllib.request.urlopen(req, timeout=60).read()
        except OSError as e:
            failed += 1
            print(f"{r['file']}: {e}")
            continue
        if sha(data) == r["sha256"]:
            put("table_codes", r["file"], data)
            got += 1
        else:
            changed += 1
        print(f"{r['file']}: {'same as on ' + r['fetched'][:10] if sha(data) == r['sha256'] else 'changed since ' + r['fetched'][:10]}")
    print(f"fetched {len(todo)}: kept {got}, changed since our download {changed}, failed {failed}")


def copy(dirs):
    want = {}
    for kind in ("table_codes", "unpublished_outputs"):
        for r in listed(kind).values():
            want.setdefault(os.path.basename(r["file"]), []).append((kind, r))
    kept = set()
    for top in dirs:
        for dp, dns, fns in os.walk(top):
            dns[:] = [d for d in dns if not d.startswith(".")]
            for fn in fns:
                if fn not in want:
                    continue
                with open(os.path.join(dp, fn), "rb") as f:
                    data = f.read()
                h = sha(data)
                for kind, r in want[fn]:
                    if r["sha256"] == h and (kind, r["file"]) not in kept:
                        put(kind, r["file"], data)
                        kept.add((kind, r["file"]))
    for kind in ("table_codes", "unpublished_outputs"):
        print(f"{kind}: {sum(1 for k, _ in kept if k == kind)} of {len(listed(kind))} found with the listed SHA-256")


def status():
    for kind in ("table_codes", "unpublished_outputs"):
        rows = listed(kind)
        n = 0
        for r in rows.values():
            p = os.path.join(CACHE, "table" if kind == "table_codes" else "outputs", r["file"])
            if os.path.exists(p) and sha(open(p, "rb").read()) == r["sha256"]:
                n += 1
        print(f"{kind}: cache/ holds {n} of {len(rows)}")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("fetch", "copy", "status") or (sys.argv[1] == "copy" and len(sys.argv) < 3):
        sys.exit(__doc__)
    {"fetch": lambda: fetch(set(sys.argv[2:])), "copy": lambda: copy(sys.argv[2:]), "status": status}[sys.argv[1]]()
