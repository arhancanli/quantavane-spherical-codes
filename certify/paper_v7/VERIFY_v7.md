# VERIFY: where every number of paper_v7 comes from

Paper: `em_v7.tex` (Experimental Mathematics) and `mcom_v7.tex` (Mathematics of Computation), one body `body_v7.tex`,
abstract `abstract_v7.tex`, macros `macros_v7.tex`. Written 7-8 November 2026. Nothing here was committed, pushed,
uploaded or sent.

## 1. Numbers, tables and figures produced by programs

Every count, quantile, distance, radius, ratio and list of entries in the paper is a macro of `v7_numbers.tex`, a cell
of `v7_tab_{census,verdicts,flow,quant,flexcat,live,notable}.tex`, or a point of `fig_{verdicts,delta,cost}_v7.pdf`.
All are written by `make_numbers_v7.py` (the version 6 script with renamed outputs) together with
`numbers_v7_extra.py`; together they also write **`VERIFY_numbers_v7.md`: every macro, its value and its source file
and field**, `v7_numbers.json` and the per-entry file `v7_verdicts.tsv`. They read only:

- released census outputs in `~/quantavane-spherical-codes/certify/` (as in version 6: `census/index.tsv`,
  `census/results_v4/`, `census/strict_v4/`, `census/exact_v4/{strict,big}/`, `census/validate_v4.tsv`,
  `census/summary_v4.json` for seven macros, `v5/out/`, `v5/live/live_20261107.tsv`);
- the version 6 runs of this folder (`phase2/`, `recount/diag_fail*.jsonl`, `checks/`), unchanged;
- the version 7 runs of this folder:

| output | program | what |
|---|---|---|
| `verify/out_exact/*.json` | `verify/exactcheck.py` | the exact decision of every PASS, LOCALMIN and NONOPT certificate (decision of record) |
| `verify/out_exact2/*.json` | `verify/exactcheck2.py` (`run_v2.sh`) | version 2 (second stress enclosure, Lemma 3.4(b)) on the 3 certificates where version 1 left a sign undecided: (3,613), (3,877) PASS, (3,791) NONOPT |
| `verify/out_strict7/*.json` | `verify/strict7.py` (= `census/strict.py` + `exactpatch.py`) | strictness by recognition with exact (S') and span |
| `verify/out_strict7b/d19_n102.json` | `verify/strict7b.py` (= `strict.py` + `exactpatch2.py`) | the one strictness re-check repeated after the untrusted LP hit its time limit (stress_lp_status 1) |
| `verify/out_sym7/*.json` | `verify/symstrict7.py` (= `phase2/symstrict.py` + `exactpatch.py`) | strictness by symmetry with exact (S') and span |
| `recount/flexgram7_*.jsonl` | `recount/flexgram7.py` | tau', delta and the SHA-256 of E for every FLEX entry |
| `recount/flexsample.tsv(.categories)` | `recount/flexsample.py` | the categories and the stratified sample (untrusted) |
| `recount/flexexact.jsonl` | `recount/flexexact.py` | the negative certificate in exact arithmetic on the sample |
| `recount/gapcut_*.jsonl`, `gapcut_sample.tsv`, `gapcut_cert.jsonl` | `recount/gapcut.py`, `gapsample.py` | gap-based proposals (scan untrusted; certificates exact) |
| `example/example_d3_n7.json`, `example_d5_n16.json` | `example/example.py` | the worked examples (same exact code, more fields) |
| `live/live_select.tsv`, `live/out/*.json`, `live/fetch_manifest.tsv` (815 files, all HTTP 200, fetched 7 Nov 2026 21:26-21:47 UTC) | `live/select_live.py`, `fetch_live.py`, `run_live.py` | the currency check on fetched codes |

`make_numbers_v7.py` asserts: every PASS/LOCALMIN of the census has an exact output; every exact 30-digit value interval
overlaps the first pass's; every flexgram7 proposal matches the census record; the CJKT kissing configurations are
certified strict; and the version 6 assertions (class counts against the version 4 and 5 summaries, ranks, SHA-256 of the
certifiers in every output).

## 2. Program versions

Every output records the SHA-256 of the programs that produced it (`program_sha256`, `exact_kernel`, `code_sha256`).
`exactcheck.py` was not changed after the runs started: the SHA-256 recorded by `exactpatch.py` in the strictness
outputs and by `exactcheck.py` in `out_exact` is the same one. `exactcheck2.py` is `exactcheck.py` plus the second
stress enclosure, tried only when the first leaves a sign undecided; it was run only where version 1 failed (3 entries),
and `make_numbers_v7.py` asserts that its outputs cover exactly those entries. The environment is in `public_changes/v7/requirements-lock.txt`
(Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1, python-flint 0.6.0, mpmath 1.3.0, matplotlib 3.9.4).

## 3. Numbers written in the text (not macros)

| text | value | source |
|---|---|---|
| tolerances | 1e-13 ... 1e-7, order 1e-10, 1e-11, 1e-12, 1e-13, 1e-9, 1e-8, 1e-7 | `census/run_census.py` TOLS, FLEX_TOLS |
| size limit | 3,100 unknowns (`\SizeCap`) | `run_census.py` CAP |
| grid of z~ | 2^-100 | `cq.py` K=100 |
| rows of C in the exact decision | 60 bits; blocks of 384 rows; Arb at 256 bits | `exactcheck.py` CBITS, BLOCK, APREC |
| r tried | first power of two from about 2 max abs(CF), six tries | `exactcheck.py` square_test |
| flexexact | 60-bit right inverse, column blocks of 256, at most 1.2e7 entries of L | `flexexact.py` CBITS, BLOCK, MAX_CELLS |
| gap scan | spread at most 1e-4, 40 candidate cuts, top 6 n d pairs | `gapcut.py` S_MAX, NCAND |
| gap sample | ten per (status, spread band), at most 1,500 unknowns, seed 7 | `gapsample.py` |
| live sample | ten per (category, dimension band), seed 7; one request every 1.5 s | `select_live.py`, `fetch_live.py` |
| recognition caps | table polynomial up to degree 12; PSLQ up to degree 10 from 110 digits | `strict.py` MAX_DEGREE, find_poly |
| worked example | (3,7): 7 points, isostatic count 12, U = 22; (5,16): 80 contacts, isostatic count 55, U = 81 | `example/*.json` (also macros) |

## 4. What a reader can check without our machine

- Every exact decision: `zsh run_v7.sh SNAP exact` (about the CPU hours of `\XHours`), or a single entry with
  `python3 -I verify/exactcheck.py OUT --only d_n`.
- Every number: `zsh run_v7.sh SNAP numbers` reproduces `v7_numbers.tex` and the tables and figures from the outputs.
- The live comparison of Table 7 needs only `v5/live/live_20261107.tsv`; the rerun on fetched codes needs the codes,
  which are not redistributed (their SHA-256 are in `live/fetch_manifest.tsv`; the live site may have changed since).
