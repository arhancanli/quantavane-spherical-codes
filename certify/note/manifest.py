# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
"""Write MANIFEST_note.tsv: path, bytes and SHA-256 of every program, data file and output of certify/note/.
  python3 -I certify/note/manifest.py"""
import hashlib, os
HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for sub in ("", "data", "out"):
    base = os.path.join(HERE, sub)
    for f in sorted(os.listdir(base)):
        p = os.path.join(base, f)
        if os.path.isdir(p) or f.startswith(".") or f in ("MANIFEST_note.tsv", "manifest.log", "RELEASE_note.tsv") or f.endswith(".db"):
            continue
        b = open(p, "rb").read()
        rows.append((os.path.join("note", sub, f) if sub else os.path.join("note", f), len(b), hashlib.sha256(b).hexdigest()))
with open(os.path.join(HERE, "MANIFEST_note.tsv"), "w") as out:
    out.write("# certify/note: path, bytes, sha256\n")
    for r in rows:
        out.write("%s\t%d\t%s\n" % r)
print(len(rows), "files")
