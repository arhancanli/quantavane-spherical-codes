"""The release plan of certify/ (paper, version 5): which files go to the git repository and which to the Zenodo
deposit, with sizes; and a manifest (SHA-256, bytes) of the version 5 files.

  python3 -I certify/v5/release_list.py          # writes certify/RELEASE_FILES.tsv and certify/v5/MANIFEST_v5.tsv

Rule: programs, documentation and small outputs go to git (aim: under 20 MB in all); the one-file-per-entry census
outputs (census/results_v4/, 16,265 files), the proposal file census/classify.jsonl and the outputs of version 2 of the
paper (results/, lists/, logs/ at this level) go to Zenodo. Superseded version 3 outputs (see .gitignore) go nowhere.
"""
import hashlib, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.dirname(HERE)

ZENODO_PREFIXES = ("census/results_v4/", "census/classify.jsonl", "results/", "lists/", "logs/")
IGNORED_PREFIXES = ("census/results/", "census/results_over_cap/", "census/exact/", "census/candidates/",
                    "census/improve/", "census/uploaded_results/", "census/paper_out/", "v5/out/tilt/")
IGNORED_FILES = {"census/candidates.tsv", "census/summary.json", "census/text.json", "census/records.json",
                 "census/validate.tsv", "census/cap.json", "census/run_codes.py", "census/run_improve.sh",
                 "census/run_exact2.sh", "census/improve.py", "census/check_candidates.py", "census/records_v3.py",
                 "census/post.sh"}


def classify(rel):
    base = os.path.basename(rel)
    if "__pycache__" in rel or base.endswith(".pyc") or base in (".DS_Store", "ruvector.db", "agentdb.rvf",
                                                                    "agentdb.rvf.idmap.json"):
        return None
    if rel.startswith((".claude-flow", ".swarm")) or "/.claude-flow/" in rel or "/.swarm/" in rel:
        return None
    if rel in IGNORED_FILES or rel.startswith(IGNORED_PREFIXES):
        return None
    if rel.startswith("census/classify_") and rel.endswith(".jsonl"):
        return None
    if rel.startswith("census/logs/") and "v4" not in base:
        return None
    if rel.startswith(ZENODO_PREFIXES):
        return "zenodo"
    return "git"


def main():
    rows = []
    for root, dirs, files in os.walk(CERT):
        dirs[:] = [x for x in dirs if x not in (".claude-flow", ".swarm", "__pycache__")]
        for f in files:
            p = os.path.join(root, f)
            rel = os.path.relpath(p, CERT)
            dest = classify(rel)
            if dest is None:
                continue
            rows.append((rel, os.path.getsize(p), dest))
    rows.sort()
    tot = {"git": 0, "zenodo": 0}
    cnt = {"git": 0, "zenodo": 0}
    for rel, size, dest in rows:
        tot[dest] += size
        cnt[dest] += 1
    with open(os.path.join(CERT, "RELEASE_FILES.tsv"), "w") as f:
        f.write("# certify/ release plan (paper version 5): path relative to certify/, bytes, destination\n")
        f.write(f"# git: {cnt['git']} files, {tot['git']} bytes; zenodo: {cnt['zenodo']} files, {tot['zenodo']} bytes\n")
        # directories with many files are listed once, with their totals
        agg = {}
        for rel, size, dest in rows:
            key = rel
            for pre in ("census/results_v4/", "census/strict_v4/", "census/exact_v4/strict/", "census/exact_v4/big/",
                        "census/exact_v4/codes/", "results/"):
                if rel.startswith(pre):
                    key = pre
            a = agg.setdefault((key, dest), [0, 0])
            a[0] += size
            a[1] += 1
        for (key, dest), (size, n) in sorted(agg.items()):
            f.write(f"{key}\t{size}\t{dest}" + (f"\t({n} files)" if n > 1 else "") + "\n")
    with open(os.path.join(HERE, "MANIFEST_v5.tsv"), "w") as f:
        f.write("# certify/v5: path, bytes, sha256\n")
        for rel, size, dest in rows:
            if rel.startswith("v5/") and not rel.endswith("MANIFEST_v5.tsv"):
                h = hashlib.sha256(open(os.path.join(CERT, rel), "rb").read()).hexdigest()
                f.write(f"{rel}\t{size}\t{h}\n")
    print(f"git: {cnt['git']} files, {tot['git'] / 1e6:.1f} MB; zenodo: {cnt['zenodo']} files, {tot['zenodo'] / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
