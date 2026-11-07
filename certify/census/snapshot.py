"""The archived table: Cohn's Table of spherical codes, DSpace@MIT handle 1721.1/153543 (data as of 31 July 2026).

  python3 -I certify/census/snapshot.py SNAPDIR          # writes SNAPDIR/data.txt (if missing) and census/index.tsv

SNAPDIR holds the files of the data set as downloaded (table.txt, data-1.txt ... data-12.txt). The readme of the data
set says that data-1.txt ... data-12.txt are to be concatenated into data.txt. Each code in data.txt is a block of
lines "Dimension: d", "Number of points: N", "Cosine: t" (rounded up), optionally "Minimal polynomial: p" and
"Provably optimal", then "References and notes: ...", "Points:" and N lines of d comma-separated coordinates.
The files are read as data only. The index records, for every code, its byte offset in data.txt.
"""
import hashlib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
INDEX = os.environ.get("CENSUS_INDEX", os.path.join(HERE, "index.tsv"))
COLS = ["d", "n", "cos", "poly", "optimal", "flagged", "listed", "offset", "notes"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def build(snapdir):
    data = os.path.join(snapdir, "data.txt")
    if not os.path.exists(data):
        with open(data + ".part", "wb") as out:
            for k in range(1, 13):
                with open(os.path.join(snapdir, f"data-{k}.txt"), "rb") as f:
                    for chunk in iter(lambda: f.read(1 << 22), b""):
                        out.write(chunk)
        os.rename(data + ".part", data)
    rows, cur = [], None
    with open(data, "rb") as f:
        off = 0
        for raw in f:
            line = raw.decode("utf-8").rstrip("\n")
            if line.startswith("Dimension: "):
                cur = {"d": int(line[11:]), "offset": off, "poly": "", "optimal": 0, "notes": ""}
            elif line.startswith("Number of points: "):
                cur["n"] = int(line[18:])
            elif line.startswith("Cosine: "):
                cur["cos"] = line[8:].strip()
            elif line.startswith("Minimal polynomial: "):
                cur["poly"] = line[20:].strip()
            elif line.strip() == "Provably optimal":
                cur["optimal"] = 1
            elif line.startswith("References and notes: "):
                cur["notes"] = line[22:].strip().replace("\t", " ")
            elif line.startswith("Points:"):
                cur["flagged"] = int("needs more optimization" in cur["notes"])
                rows.append(cur)
            off += len(raw)
    listed = set()
    for line in open(os.path.join(snapdir, "table.txt")):
        p = line.strip().split(",", 3)
        if len(p) >= 3:
            listed.add((int(p[0]), int(p[1])))
    for r in rows:  # table.txt omits the codes of data.txt for which a code with more points is as good
        r["listed"] = int((r["d"], r["n"]) in listed)
    with open(INDEX, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in COLS) + "\n")
    meta = {"data.txt": {"bytes": os.path.getsize(data), "sha256": sha256(data)}, "codes": len(rows)}
    for name in ["table.txt", "readme.pdf"] + [f"data-{k}.txt" for k in range(1, 13)]:
        p = os.path.join(snapdir, name)
        meta[name] = {"bytes": os.path.getsize(p), "sha256": sha256(p)}
    json.dump(meta, open(os.path.join(HERE, "snapshot_files.json"), "w"), indent=1)
    print(len(rows), "codes indexed")


def load_index():
    out = []
    with open(INDEX) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            v = line.rstrip("\n").split("\t")
            r = dict(zip(head, v))
            for c in ("d", "n", "optimal", "flagged", "listed", "offset"):
                r[c] = int(r[c])
            out.append(r)
    return out


def read_points(snapdir, entry):
    """The N points of an entry as decimal strings (lists of d strings), read from data.txt at the indexed offset."""
    d, n = entry["d"], entry["n"]
    pts = []
    with open(os.path.join(snapdir, "data.txt"), "rb") as f:
        f.seek(entry["offset"])
        seen = False
        for raw in f:
            line = raw.decode("utf-8").strip()
            if not seen:
                seen = line.startswith("Points:")
                continue
            parts = line.split(",")
            if len(parts) != d:
                break
            pts.append(parts)
            if len(pts) == n:
                break
    if len(pts) != n:
        raise ValueError(f"d={d} n={n}: read {len(pts)} points")
    return pts


if __name__ == "__main__":
    build(sys.argv[1])
