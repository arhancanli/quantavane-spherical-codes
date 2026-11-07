# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""The release list for the note: every file of certify/note and certify/v5, and the files of certify/ that they use,
with size, SHA-256 and destination (git or Zenodo), in certify/note/RELEASE_note.tsv. Nothing is committed or pushed.

  python3 -I certify/note/release_note.py

Rules: those of certify/v5/release_list.py (tool-state files of editor or agent plugins, caches and superseded outputs
are left out; everything here goes to git, as no file is large). Also checks that certify/v5/MANIFEST_v5.tsv still
matches the version 5 files it lists, and certify/note/MANIFEST_note.tsv the note files.
"""
import hashlib, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(CERT, "v5"))
from release_list import classify          # the version 5 rules (path relative to certify/ -> 'git', 'zenodo' or None)

EXTRA = ["LICENSE", ".gitignore", "README.md", "census/snapshot.py", "census/index.tsv"]   # used by v5/jam116.py


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def manifest_ok(path):
    bad, n = [], 0
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        rel, size, h = line.rstrip("\n").split("\t")
        p = os.path.join(CERT, rel)                  # manifest paths are relative to certify/ (note/..., v5/...)
        n += 1
        if not os.path.exists(p) or os.path.getsize(p) != int(size) or sha(p) != h:
            bad.append(rel)
    return n, bad


def main():
    rows = []
    for sub in ("note", "v5"):
        for root, dirs, files in os.walk(os.path.join(CERT, sub)):
            dirs[:] = sorted(x for x in dirs if x not in (".claude-flow", ".swarm", "__pycache__"))
            for f in sorted(files):
                rel = os.path.relpath(os.path.join(root, f), CERT)
                if rel == "note/RELEASE_note.tsv":
                    continue
                dest = classify(rel)
                if dest:
                    rows.append(rel)
    rows += [x for x in EXTRA if os.path.exists(os.path.join(CERT, x))]
    rows.append("note/RELEASE_note.tsv")
    n5, bad5 = manifest_ok(os.path.join(CERT, "v5", "MANIFEST_v5.tsv"))
    nn, badn = manifest_ok(os.path.join(HERE, "MANIFEST_note.tsv"))
    with open(os.path.join(HERE, "RELEASE_note.tsv"), "w") as out:
        out.write("# Release list for the note (certify/note, certify/v5 and the files of certify/ they use): path from the\n"
                  "# repository root, bytes, SHA-256, destination. Written by certify/note/release_note.py; nothing has been\n"
                  "# committed. To commit and push these files in one step, from the repository root:\n"
                  "#   git add -- $(awk -F'\\t' '!/^#/ {print $1}' certify/note/RELEASE_note.tsv) && "
                  "git commit -m 'certify: note and version 5 certificates' && git push\n"
                  "# The note cites the public commit in the macros \\CertCommit and \\CertCommitFull of note.tex;\n"
                  "# after pushing these files, update them to the new commit if it is to be cited.\n"
                  "# Manifests checked: v5/MANIFEST_v5.tsv %d files, %d mismatches; note/MANIFEST_note.tsv %d files, %d "
                  "mismatches.\n" % (n5, len(bad5), nn, len(badn)))
        total = 0
        for rel in rows:
            p = os.path.join(CERT, rel)
            if rel == "note/RELEASE_note.tsv":
                out.write("certify/%s\t-\t-\tgit\n" % rel)
                continue
            size = os.path.getsize(p)
            total += size
            out.write("certify/%s\t%d\t%s\tgit\n" % (rel, size, sha(p)))
        out.write("# %d files, %d bytes\n" % (len(rows), total))
    print({"files": len(rows), "bytes": total, "v5_manifest_mismatches": bad5, "note_manifest_mismatches": badn})


if __name__ == "__main__":
    main()
