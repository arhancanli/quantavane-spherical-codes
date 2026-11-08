# VERIFY: where every number of paper_v8 comes from

Paper: `em_v8.tex` (Experimental Mathematics) and `mcom_v8.tex` (Mathematics of Computation), one body `body_v8.tex`,
abstract `abstract_v8.tex`, macros `macros_v8.tex`. Written 8 November 2026. Nothing here was committed, pushed,
uploaded or sent.

## 1. Numbers, tables and figures produced by programs

Every count, quantile, distance, radius, ratio and list of entries in the paper is a macro of `v8_numbers.tex`, a cell
of `v8_tab_{census,verdicts,flow,quant,flexcat,live,notable,above,structure}.tex`, or a point of
`fig_{verdicts,delta,cost}_v8.pdf`. All are written by `make_numbers_v8.py` with `numbers_v8_extra.py`, which also write
**`VERIFY_numbers_v8.md` (every macro, its value and its source file and field)**, `v8_numbers.json` and the per-entry
file `v8_verdicts.tsv`. `make_numbers_v8.py` first checks that every input exists (`REQUIRED`) and writes nothing if
one is missing. It reads only the released census outputs in `certify/` (`census/`, `v5/`) and the files of this
directory:

| output | program | what |
|---|---|---|
| `verify/out_exact/`, `verify/out_exact2/` | `verify/exactcheck.py`, `exactcheck2.py` | exact decision of every PASS, LOCALMIN, NONOPT certificate |
| `verify/out_strict7/`, `out_strict7b/`, `out_sym7/` | `verify/strict7.py`, `strict7b.py`, `symstrict7.py` | strictness by recognition and symmetry, (S') and span exact |
| `recount/out_flex8/` | `recount/flexexact8.py` | **the exact negative certificate of every FLEX entry** (decision of record) |
| `recount/out_flex8r/` | `recount/flexrefine8.py` | the same with a right inverse from an SVD, where flexexact8's fail |
| `recount/flexgram7_*.jsonl` | `recount/flexgram7.py` | the double-precision separation (comparison only) |
| `above/out/`, `above/out_more/`, `above/gan510.json` | `above/run_above.py`, `above_more.py`, `gan510.py` | the structured entries above the size limit (second pass: larger cap on the classes) |
| `recount/tolscan_other_*.jsonl`, `tolscan_close.jsonl` | `recount/tolscan.py` | other tolerances; close calls |
| `recount/zetasample.jsonl` | `recount/zetasample.py` | sharper rigidity constant (untrusted) |
| `recount/gapcut_*.jsonl`, `gapcut_cert.jsonl` | `recount/gapcut.py`, `gapsample.py` | gap-based proposals |
| `example/*.json` | `example/example.py` | the worked examples |
| `live/out/`, `live/fetch_manifest.tsv` | `live/run_live.py`, `fetch_live.py`, `select_live.py` | the currency check |
| `phase2/`, `checks/`, `recount/diag_fail*.jsonl`, `recount/arbcheck.jsonl` | version 6 programs (unchanged) | symmetry certificate, checks, diagnosis, Arb sample |

Assertions in the scripts: every census PASS/LOCALMIN has an exact output and the 30-digit intervals overlap; version 2
outputs cover exactly the version 1 failures; every FLEX entry has an exact flexexact8 output, and every one that is not FLEX a flexrefine8 output
(an entry that fails both loses its verdict); every confirmed one has the census's proposal; the above-limit outputs cover exactly the 155 entries; the CJKT controls are strict; the seven flat
orthogonal unions; and the version 6 assertions (class counts against the version 4 and 5 summaries).

## 2. Program versions

Every output records the SHA-256 of its programs (`program_sha256`, `exact_kernel`, `code_sha256`). `run_above.py`,
`above_more.py`, `flexexact8.py` and `flexrefine8.py` were not changed after their production runs started; two
earlier runs of earlier versions of `run_above.py` (16 and 41 entries) were stopped and their outputs discarded. Every
output was checked to record the SHA-256 of the released file (12,482 + 1 negative certificates, 155 above-limit
records and the 19 of the second pass, whose 8 completed ones also record run_above.py's SHA-256 while the 11 stopped
at the time or memory limit carry only above_more.py's, 1,778 + 3 exact decisions, 815 live records). The whole numbers step was rerun from a fresh local clone of
c4a3dcc8 with `public_changes/v8/make_copy_v8.sh` and reproduced every output byte for byte. Environment: Python 3.9.6, NumPy 2.0.2, SciPy
1.13.1, python-flint 0.6.0, mpmath 1.3.0, matplotlib 3.9.4 (`public_changes/v8/requirements-lock.txt`).

## 3. Numbers written in the text (not macros)

| text | value | source |
|---|---|---|
| tolerances | 1e-13 ... 1e-7, order 1e-10, 1e-11, 1e-12, 1e-13, 1e-9, 1e-8, 1e-7 | `census/run_census.py` |
| size limit | 3,100 unknowns (`\SizeCap`) | `run_census.py` CAP |
| grid of z~ | 2^-100 | `cq.py` |
| exact decision | rows of C at 60 bits; blocks of 384; Arb at 256 bits | `exactcheck.py` |
| negative certificate | grid 2^-52 for X0; C at 52 bits; limbs of 26/27 bits; QR when beta >= 1/2 and L has at most 1.2e7 entries | `flexexact8.py` |
| negative certificate, third right inverse | SVD (double precision), only where flexexact8's two give beta >= 1 or none; at most 1.2e7 entries | `flexrefine8.py` |
| resolved | delta >= 1e-12 and delta = tau'/4 | `make_numbers_v8.py` |
| above the limit | p < 2^16 (65521, or the largest prime below 2^16 in which D is a square); S = M^T, then a random S with the next prime; matrix stored as 16-bit integers; panels of 256 columns, inner blocks of 16; at most 20,000 columns; 1.5 GB per entry; class-constant stresses up to 300 classes (FLINT nullspace of B^T B); exact square solve up to 800 columns; LP classes up to 4,000 | `run_above.py` |
| tolerance scan | one proposal per distinct (contacts, non-rattlers), at most 3,100 unknowns | `tolscan.py` |
| close calls | next pairs within 10 times the gap of the next one | `tolscan.py` |
| zeta sample | isostatic certificates with at most 1,600 unknowns | `zetasample.py` |
