# VERIFY: where every number of paper_v9 comes from

Paper: `em_v9.tex` (Experimental Mathematics) with the body `tex/body_v9.tex`, the abstract `abstract_v9.tex`, the
macros `macros_v9.tex`, and the supplement `supp_v9.tex`. Written 8 November 2026. Nothing here was committed, pushed,
uploaded or sent.

## 1. Numbers, tables and figures produced by programs

Every count, quantile, distance, radius, ratio and list of entries in the paper and the supplement is a macro of
`v9_numbers.tex`, a cell of a table `v9_tab_*.tex` or `v9s_tab_*.tex`, or a point of `fig_verdicts_v9.pdf` or
`fig_cost_v9s.pdf`. All are written by `make_numbers_v9.py`, which also writes **`VERIFY_numbers_v9.md` (every macro,
its value and its source file and field)**, `v9_numbers.json` and the per-entry file `v9_verdicts.tsv`.
`make_numbers_v9.py` first checks that every input exists and writes nothing if one is missing (`--draft`, for trial
compilations only, fills missing macros with ?? and is not used for the paper).

Inputs:

| output | program | what |
|---|---|---|
| `paper_v8/v8_numbers.json`, `paper_v8/v8_verdicts.tsv` | `make_numbers_v8.py` (released) | the decisions of record of version 8 and their macros (passed through where their meaning is unchanged) |
| `paper_v8/verify/out_exact*`, `out_strict7`, `out_sym7`, `above/out*`, `recount/out_flex8*`, `recount/tolscan_close.jsonl`, `live/out` | version 8 programs (released) | per-entry constants (rho, kappa, inside rho), distances, close calls and double certificates, live reruns |
| `v5/live/live_20261107.tsv` | released | the live table of 7 November 2026, with its credits (author-credited entries: "Arhan Canli") |
| `../codes/MANIFEST.tsv` | released | the author's uploads to the live table (dates) |
| `recog/out_all/` | `recog/recog9.py` | the survey of the 630 open entries (2,000 bits, LLL degree <= 30) |
| `recog/out_pass2/` | `recog/recog9.py --pass2` | the second pass (4,000 bits) on entries with at most 20 distinct inner products |
| `recog/out_strict9/` | `recog/strict9.py` | **the exact strictness decision of the recognised open entries** |
| `recog/out_strict9_ctl/` | `recog/strict9.py --class strict_hyper` (the shards of `recog/ctl_plan9.py`) | the version 8 certificates by recognition with at most 2,300 unknowns, decided again (controls; 234 of 258) |
| `exactdata/open`, `exactdata/v8recog` | `strict9.py --data` | the exact descriptions of the recognised codes (field, elements, index matrix) |
| `exactdata/core` | `exactdata/exactdata_core9.py` | the exact descriptions of the strict exact cores |
| `exactdata/INDEX.tsv` | `exactdata/index9.py` | one line per exact description (field, exact minimal polynomial of the cosine, surveyed symmetry order) |
| `prestress/out/` | `prestress/prestress9.py` | exact flexes and the second-order form for (22,166), (11,604) |
| `prestress/explore/*.log` | `prestress/explore9_*.py` | untrusted numerical explorations, quoted as such |
| `cut9/out/` | `cut9/cut9.py` (decision: `paper_v8/recount/flexexact8.check`) | the negative certificate at the intrinsic cut |
| `ranks/out/` | `ranks/flintrank9.py` | FLINT ranks of the strict entries above the size limit (72 computed; the 14 largest recorded as skipped for memory) |
| `ranks/pari_check9.json` | `ranks/pari_check9.py` | PARI/GP spot check of exact descriptions |
| `certdata/files/`, `export_log.json`, `verify.json` | `certdata/export9.py`, `verify_cert9.py` | certificate data of the subset, and their stand-alone decision |
| `example/example_d5_n16.json`, `example_d3_n7.json` | `paper_v8/example/example.py` (released) | the worked examples, regenerated |

Assertions in `make_numbers_v9.py`: the version 8 file covers the census; every entry moved to strict by `strict9.py`
was open ('local') in version 8; every control (optimal and unique, CJKT, d = 3 with N <= 11) is strict; every provably
optimal entry has an optimality reference; the open entries surveyed are exactly the 630 of version 8; for the
N = 2d + 1 family the number of contacts is the isostatic count.

## 2. Program versions

Every output records the SHA-256 of its programs (`program_sha256`). `recog9.py`, `rec9lib.py`, `strict9.py`,
`prestress9.py`, `cut9.py`, `flintrank9.py`, `export9.py` and `verify_cert9.py` were not changed after their
production runs started; three earlier trial runs of `recog9.py` (a random sample of 40 with a version that accepted LLL
relations too readily, and two runs stopped within minutes when the relation test was fixed) and one of an exploratory
linear program were discarded. Environment: Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1, python-flint 0.6.0, mpmath 1.3.0,
matplotlib 3.9.4, cypari 2.5.7 (PARI 2.15.4) (`public_changes/v9/requirements-lock.txt`).

## 3. Numbers written in the text (not macros)

| text | value | source |
|---|---|---|
| refinement grid and recognition | 2^-2000 (about 600 digits), LLL degree at most 30, relations of at most 3P/4 bits | `recog/rec9lib.py`, `recog9.py` |
| second pass | 4,000 bits (about 1,200 digits), at most 20 distinct inner products | `run_v9.sh pass2` |
| controls | at most 2,300 unknowns | `recog/ctl_plan9.py --max-unknowns 2300` |
| primes of the rank cross-check | the largest primes below 2^26 (with D a square) | `ranks/flintrank9.py`, `kfield.prime_for` |
| dense FLINT matrices computed | up to 1,280 MiB (`--max-mb 1280`); FLINT's rank routine works on a copy | `run_v9.sh ranks`, `ranks/flintrank9.py` |
| intrinsic cut search | windows moved by 10, 10^(1/2), 10^(1/4); cuts at the largest gap within 10^(1/8); at most 3,100 rows; beta < 0.9 | `cut9/cut9.py` |
| sample of the cut search | 400 random entries of the others (seed 20261108) | `cut9/cut9.py` |
| sample of the recognition survey | 40 random open entries (seed 20261108) | `recog/sample40.txt` |
| tiers | 10^-12 | `make_numbers_v9.py` |
| close calls | first-pass gap below 10^-9 | `make_numbers_v9.py` (v8_verdicts.tsv) |
| polishing vs new configuration | point move at least 10^-4 / below 10^-6 (lower bound from sorted inner products) | `make_numbers_v9.py` (paper_v8/live/out) |
| dates | archived snapshot 31 July 2026; census run 6 and 7 November 2026; live table 7 November 2026 | DSpace record; repository history |
