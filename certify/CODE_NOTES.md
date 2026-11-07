# Notes on the released census programs (proposed as certify/CODE_NOTES.md)

The census programs below are kept byte-identical to the versions that produced the released outputs, because every
output records their SHA-256 (cq.py c12bf1..., sphcert.py 807930..., run_census.py 4a7435... and de9662..., strict.py
1b2fba...). Changing a comment would break that link, so the following are recorded here instead of fixed in place.

- `census/strict.py`, docstring: "Theorem 2.5" refers to the numbering of paper version 4; it is Theorem 2.4 (the
  first-order criterion) of paper versions 6 and 7.
- `census/run_census.py`, docstring: "paper, Section 4.2" refers to version 4; the rule is in the section "The data
  and the run" of the census section (Section 5.1 of version 6, Section 6.1 of version 7).
- `census/cq.py`, `census/strict.py`, `census/classify.py`, `census/exact_big.py`, `census/exact_strict.py`,
  `census/summarize.py` and others prepend `~/Library/Caches/quantavane-census/pylib` (and `flintlib`) to `sys.path`
  before importing NumPy, SciPy and python-flint. On another machine, install the versions of
  `requirements-lock.txt` into that directory (`pip install --target ~/Library/Caches/quantavane-census/pylib -r
  requirements-lock.txt`), or into the interpreter's own site-packages, which Python searches after it. The programs of
  `paper_v7/` read the paths from `QV_PYLIB` and `QV_FLINTLIB` (default: the same directories), the repository from
  `QV_CERTIFY` and the archived table from `QV_SNAP`.
