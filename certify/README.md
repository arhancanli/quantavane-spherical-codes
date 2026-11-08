# Certified local optimality of spherical codes: the programs

The programs and outputs of the paper *Certified local optimality of spherical codes: a census of Cohn's table*
(A. Canli, census version 8, November 2026; the version submitted to Mathematics of Computation on 6 November 2026 is an
earlier version). They apply to the archived version of Cohn's table of spherical codes (DSpace@MIT handle
1721.1/153543, data as of 31 July 2026), entries with 3 <= d <= 32 and N <= 1024. Section, theorem and table numbers
below are those of version 8 of the paper.

**Licence.** Everything in this directory (`certify/`) is released under the MIT licence (`certify/LICENSE`,
Copyright (c) 2026 Arhan Canli (Quantavane)). The table's data remain Cohn's; the repository's other directories have
their own licence (`LICENSE` at the top).

**What is trusted.** Every function is marked TRUSTED (its result enters a decision) or UNTRUSTED (it only proposes
candidates: contacts, rattlers, Newton steps, approximate inverses, subsets, coefficients and classes of stresses,
recognised numbers, automorphisms and models of representations, colourings). Whatever a program reports as proved is
proved whichever candidates it was given. In version 8 every verdict of the paper is decided in exact integer and
rational arithmetic and Arb ball arithmetic; the double-precision kernel of the census (`sphcert.py`, `census/cq.py`)
only proposes. Table 1 of the paper lists the arithmetic of every step.

## Environment

Python 3.9 with the versions in `requirements-lock.txt` (NumPy 2.0.2, SciPy 1.13.1, python-flint 0.6.0, mpmath 1.3.0;
matplotlib 3.9.4 for the figures). The census programs of version 4 put `~/Library/Caches/quantavane-census/pylib` and
`flintlib` first on `sys.path` (see `CODE_NOTES.md`); install the libraries there, or into the interpreter's
site-packages. The programs of `paper_v8/` read the paths from `QV_PYLIB` and `QV_FLINTLIB`, the repository from
`QV_CERTIFY` and the archived table from `QV_SNAP`. Use one BLAS thread (`VECLIB_MAXIMUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`).

## Rerunning everything from a clean checkout

    git clone https://github.com/arhancanli/quantavane-spherical-codes && cd quantavane-spherical-codes
    pip install --target ~/Library/Caches/quantavane-census/pylib -r certify/requirements-lock.txt
    zsh certify/census/fetch_snapshot.sh SNAPDIR          # table.txt, data-1..12.txt, every MD5 checked
    python3 -I certify/census/snapshot.py SNAPDIR         # data.txt and census/index.tsv
    export QV_CERTIFY=$PWD/certify QV_SNAP=SNAPDIR
    zsh certify/paper_v8/run_v8.sh SNAPDIR numbers         # every number, table and figure (about a minute)

`run_v8.sh SNAPDIR all` reruns every computation of version 8 from the released census outputs, in roughly fifteen
hours on one core of our machine (shared with other work): the exact decisions of the positive certificates, about 2.7
CPU hours; the strictness re-checks, about 1.5; the exact negative certificates, about 2.3; the structured codes above
the size limit, about 8 hours of wall-clock time (`run_above.py` and the second pass `above_more.py`, whose entries stop
at 15 minutes or 1.5 GB); the other tolerances and close calls, under half an hour; the rigidity-constant sample, the
gap-based scans and the worked examples, about an hour. Each
step can be run alone (`exact`, `strict`, `flex8`, `above`, `tolscan`, `zeta`, `gapscan`, `gapcert`, `example`, `live`,
`numbers`). The census itself (version 4, about five CPU hours) and the version 5 analyses are rerun with
`certify/census/run_all.sh` and `certify/v5/run_v5.sh`; their outputs are released, so this is optional. The fetch of
current codes from the live table (`run_v8.sh SNAPDIR livefetch`, about 800 requests, one every 1.5 s) is not part of
`all`; its outputs and the SHA-256 of the fetched files are released, the files themselves are not.

## What each program proves

The census (`census/`, outputs of version 4, used unchanged): `cq.py` + `sphcert.py` (the first pass: `PASS`,
`LOCALMIN`, `NONOPT`, `FLEX` in double precision with rigorous bounds; every success is decided again exactly in
`paper_v8/`), `strict.py` (strictness by recognition, Section 3.6, Theorem 3.14), `exact_strict.py` (Proposition 3.19),
`exact_big.py` (a second program in exact integers), `run_census.py` (the rule of Section 6.1), `classify.py` (the
proposals at every tolerance). The analyses of the entries with an exact cosine (`v5/`, Section 6.7): `core.py`,
`subspace.py`, `rigidity.py` and their helpers.

Version 8 (`paper_v8/`, self-contained; it imports only `census/`, `v5/` and its own files):

- `verify/exactcheck.py`: the decision of record of every PASS, LOCALMIN and NONOPT certificate (Lemma 3.4, Theorem 3.5,
  Theorem 3.11), exact. `verify/exactcheck2.py`: the same with the second stress enclosure (Lemma 3.4(b)), used where
  the first left a sign undecided (three certificates).
- `verify/strict7.py`, `verify/symstrict7.py`, `verify/strict7b.py` (with `exactpatch.py`, `exactpatch2.py`): the
  strictness certificates of `census/strict.py` and `phase2/symstrict.py` with (S') and the span decided exactly.
- `phase2/symstrict.py`: strictness by symmetry (Section 3.7, Theorem 3.16); `phase2/sym_survey.py`: the untrusted survey.
- `recount/flexexact8.py`: the negative certificate (Theorem 4.1) decided exactly, with sparse exact products;
  `recount/flexrefine8.py`: the same with a third right inverse (from a singular value decomposition) where the census's
  two fail the exact test.
- `above/run_above.py`: the structured codes above the size limit (Section 6.8; Proposition 3.19 with the rank over
  F_p of S M, S = M^T or random, and an exact stress); `above/above_more.py`: a second pass with a larger cap on the
  classes of contacts where no stress was found. `above/gan510.py`: the exact inner products of (10,510).
- `recount/tolscan.py`: the positive certificates at the other tolerances and at the close calls (Section 6.3).
- `recount/zetasample.py`: the comparison of rigidity constants (Section 6.4; untrusted).
- `recount/gapcut.py`, `recount/gapsample.py`: gap-based proposals (Section 6.9).
- `recount/flexgram7.py`: the double-precision separation of the negative certificate (version 7, kept for comparison).
- `recount/diag_fail.py`, `recount/arbcheck.py`, `checks/check_*.py`: the diagnosis of the failures, an Arb re-check of a
  sample, two exact identity checks (version 6, unchanged).
- `live/select_live.py`, `live/fetch_live.py`, `live/run_live.py`: the currency check (Section 6.10).
- `example/example.py`: the worked example (Section 5).
- `make_numbers_v8.py` (+ `numbers_v8_extra.py`): every number, table and figure, `v8_verdicts.tsv` (the verdict of every
  entry) and `VERIFY_numbers_v8.md` (every number with its source). It checks that every input exists before writing
  anything.

Every output records the SHA-256 of the programs that produced it.

## Files: git and Zenodo

Everything the paper uses is in git: the census outputs in `census/`, all of `v5/`, and all of `paper_v8/`. A copy of
the directory at the version of the paper will be deposited at Zenodo. `paper_v6/` and `paper_v7/` hold the earlier
versions and are not needed by version 8.

## Earlier versions

`run_batch.py`, `analysis.py`, `exact_check.py`, `lp_dual.py`, `make_lists.py`, `dedup.py`, `lists/` and `results/` at this
level belong to version 2 of the paper and are kept for that record. Version 3 outputs are superseded by version 4.

## The note: a jammed spherical code that is not infinitesimally jammed (`note/`)


The programs of the note *A jammed spherical code that is not infinitesimally jammed* (A. Canli). The theorem that the
entry (22,116) is jammed but not infinitesimally jammed appeared in version 5 of the census paper (its Theorem 5.4,
with `v5/recognise116.py` and `v5/jam116.py`); from version 6 on the census paper does not contain it, and the note is
its only place. Rerun everything with

    zsh certify/note/run_note.sh SNAPDIR NOTEDIR     # about twenty-five minutes on one core; peak memory about 2.3 GB

which runs, in order: `extract.py` (the stored entry), `classes.py` (classes and the 13 exact values, from the
recognition `v5/out/recognition_d22_n116.json`), `aut116.py` (the group by individualisation and refinement),
`construct.py` (the model in the circle geometry of PG(1,9)), `model_check.py` (a second, separately written program
that builds the model from the text of the note and finds the same bijection with the deposited class matrices of X*
and X'; adapted from an anonymous referee check), `orbit_model.py` (each orbit on pairs described in the model; Tables 1
and 5), `structure.gp` (the flex, the stresses, the characters), `check.gp` (the second exact verification, including the
count of the 1440 automorphisms by backtracking: 147,129 nodes), `v5/jam116.py` (the first exact verification),
`field.gp` (the field in reduced form), then for the second code X' of Section 7 `xprime.gp` (the exact data of X' in
`data_xprime/`), `structure.gp` on that directory (`QV_NOTE_DATA=data_xprime`), `galois_xprime.py` (the three Galois
relabellings of X', `data_xprime/galois.gp`), `check_xprime.gp` (which checks them), `jam_xprime.py` (the second check of
X': an adapted copy of `v5/jam116.py` in Python with ball arithmetic, which recomputes X' from X* itself and must agree
with `data_xprime/`) and its two red controls (`--red gram`, `--red ratio`, which must fail: run_note.sh stops if either
passes), and finally `optima.py`, `summarize.py`, `family_summary.py`, `symcount_summary.py` (numbers and tables) and
`manifest.py`. The last lines printed are the STATUS lines, the two Python statuses, the model check and the red-control
line. On the author's machine (loaded) the GP scripts take two to four minutes each, `jam116.py` about a minute and
`model_check.py` about two minutes.

Two floating-point computations are run separately (hours). `zsh certify/note/run_symcount.sh SNAPDIR` runs
`symcount.py` in three shards over every entry of the table: contacts, rattlers, the rank of the rigidity matrix, the
symmetry group and its character table where needed, and the symmetry-forced count of flexes (Section 4 of the note);
outputs in `note/out/symcount/`, summarised by `symcount_summary.py`. Some steps have time limits (the backtracking
search for a group, each linear program), so on another machine a few entries may move between "undecided" and a
decided class; the outputs here were written by successive versions of `symcount.py` with the same decision rules
(the last ones only computed less where the result cannot change the class), and the entries whose stresses the note
examines were recomputed with the last two versions (which differ only in solving the max-min linear program over the invariant stresses when the group is known). `zsh certify/note/run_family.sh SNAPDIR` runs `family.py` (Gamma-invariant codes from
PGammaL(2,q) and subgroups: the restarts for X* and X', then the cases of its table `CASES`) and `family_table.py`;
outputs in `note/out/family/`. Nothing in them is a certificate; the note states which results are floating point.
