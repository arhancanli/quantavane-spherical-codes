"""MANIFEST_v4.tsv and COMMIT_LIST_v4.txt: the files of this version of certify/ (paper version 4), with their sizes
and SHA-256.

  python3 -I certify/census/manifest.py
"""
import glob, hashlib, os

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.dirname(HERE)
ROOT = os.path.dirname(CERT)

PROGRAMS = ["certify/README.md", "certify/.gitignore", "certify/sphcert.py",
            "certify/census/cq.py", "certify/census/strict.py", "certify/census/exact_strict.py",
            "certify/census/exact_big.py", "certify/census/exact_sample.py", "certify/census/snapshot.py",
            "certify/census/classify.py", "certify/census/run_census.py", "certify/census/run_all.sh",
            "certify/census/run_exact.sh", "certify/census/fetch_snapshot.sh", "certify/census/summarize.py",
            "certify/census/make_text.py", "certify/census/census_md.py", "certify/census/validate_poly.py",
            "certify/census/manifest.py"]
INPUTS = ["certify/census/index.tsv", "certify/census/classify.jsonl", "certify/census/snapshot_files.json"]
OUTPUTS = ["certify/census/validate_v4.tsv", "certify/census/summary_v4.json"]
DIRS = ["certify/census/results_v4", "certify/census/strict_v4", "certify/census/exact_v4"]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    files = PROGRAMS + INPUTS + OUTPUTS
    for d in DIRS:
        files += sorted(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(ROOT, d, "**", "*"), recursive=True)
                        if os.path.isfile(p) and not p.endswith(".part"))
    files += sorted(os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(HERE, "logs", "*v4*")))
    missing = [f for f in files if not os.path.exists(os.path.join(ROOT, f))]
    files = [f for f in files if f not in missing]
    with open(os.path.join(HERE, "MANIFEST_v4.tsv"), "w") as out:
        out.write("path\tbytes\tsha256\n")
        for f in files:
            p = os.path.join(ROOT, f)
            out.write(f"{f}\t{os.path.getsize(p)}\t{sha(p)}\n")
    with open(os.path.join(CERT, "COMMIT_LIST_v4.txt"), "w") as out:
        out.write("# Files of certify/ for paper version 4 (paths from the repository root). Directories stand for all\n"
                  "# their files; certify/census/MANIFEST_v4.tsv lists every file with its size and SHA-256.\n")
        for f in PROGRAMS + INPUTS + OUTPUTS + ["certify/census/MANIFEST_v4.tsv", "certify/COMMIT_LIST_v4.txt"]:
            if f not in missing:
                out.write(f + "\n")
        for d in DIRS:
            out.write(d + "/\n")
        out.write("certify/census/logs/*v4*\n")
    print(len(files), "files;", "missing:", missing)


if __name__ == "__main__":
    main()
