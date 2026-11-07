# Certified local optimality of spherical codes: the programs

The programs and outputs of the paper *Certified local optimality of spherical codes* (A. Canli, census version 6,
November 2026; the version submitted to Mathematics of Computation on 6 November 2026 is an earlier version). They apply
to the archived version of Cohn's table of spherical codes (DSpace@MIT handle 1721.1/153543, data as of 31 July 2026),
entries with 3 <= d <= 32 and N <= 1024. Section, theorem and table numbers below are those of the paper.

**Licence.** Everything in this directory (`certify/`) is released under the MIT licence (`certify/LICENSE`,
Copyright (c) 2026 Arhan Canli (Quantavane)). The table's data remain Cohn's; the repository's other directories have
their own licence (`LICENSE` at the top).

**What is trusted.** Every function is marked TRUSTED (its result enters a decision, with exact integer or rational
arithmetic, arithmetic in a number field, ball arithmetic with rigorous radii, or floating point with rigorous error
bounds formed only from individual products) or UNTRUSTED (it only proposes candidates: contacts, rattlers, Newton steps,
approximate inverses, subsets, coefficients of stresses, recognised numbers, automorphisms and models of
representations, colourings, flex patterns). Whatever a program reports as proved is proved whichever candidates it was
given. Table 1 of the paper lists the arithmetic of every step.

## What each program proves

The census (`census/`, outputs of version 4, used unchanged by the paper):

- `cq.py` + `sphcert.py`: `PASS` (isostatic graph, Theorem 3.4: an exact configuration with exactly the proposed
  contacts exists within the printed distance; a strict local minimum, with the explicit constants of Theorem 3.8),
  `LOCALMIN` (hyperstatic graph, Theorem 3.10: a local minimum within the printed distance of the zero P' of the square
  subsystem, plus the tube `localmin_displacement`), `NONOPT` (Proposition 3.5: the exact equal-contact configuration
  of the proposed graph is not a local minimum), `FLEX` (Theorem 4.1, the negative certificate: no configuration within
  the printed l1 radius h carries a nonzero stress on the pairs within tau; with the tolerance this proves that no local
  minimum lies within delta = min(h/sqrt d, tau'/4) of the stored code; `paper_v6/recount/flexgram.py` computes tau').
- `strict.py`: `STRICT` on a LOCALMIN entry (Section 3.6: the Gram matrix recognised in Q(c*), verified exactly, then
  Theorem 3.13).
- `exact_strict.py`: Theorem 2.4 with no floating point in any decision (Proposition 3.18), for exact coordinates or a
  rational Gram matrix. `exact_big.py`: a second program (exact integers and rationals) that re-checks a sample.
- `run_census.py`: the rule of Section 5.1. It has two SHA-256 values in the outputs: 4a7435... for the first 15,992
  entries and de9662... for the last 273, which added the memory watchdog (`--mem-limit-mb`); the certifiers are the same.

The analyses of the entries with an exact cosine (`v5/`, Section 5.6): `core.py` (the exact core), `subspace.py` (flat
cores, Proposition 2.7), `rigidity.py` (split contact graphs, Proposition 2.8; bar-rigid cores, Theorem 2.4),
`distances.py` (certified distances of the exact-only confirmations), `maxcos_exact.py`, `families.py` (untrusted),
`summarize_v5.py`, and the helpers `kfield.py`, `tower.py`, `dist.py`, `exact2.py`. (`recognise116.py` and `jam116.py`
belong to the separate note on (22,116); the paper does not use them.)

Version 6 (`paper_v6/`):

- `phase2/symstrict.py`: `STRICT` on a LOCALMIN entry by symmetry (Section 3.7, Theorem 3.15): an exact model of the
  code's symmetry over Q(sqrt D), checked exactly, a Krawczyk test on the symmetry-reduced system in Arb ball
  arithmetic, then (S'), the span and the rattler gaps by the kernel of `cq.py` and `strict.py`. `--red model|perm|swap`
  runs the red controls of Section 5.4. `phase2/sym_survey.py`: the untrusted survey of the symmetries.
- `recount/flexgram.py`: the separation tau' of the negative certificate (Theorem 4.1), recomputed with no dense BLAS
  product in any bound; `make_numbers_v6.py` takes delta = min(h/sqrt d, tau'/4) from its output.
- `recount/arbcheck.py`: an independent re-check of a sample of existence proofs in Arb ball arithmetic (Section 5.5).
- `recount/diag_fail.py`: the untrusted diagnosis of the existence failures (Section 5.7).
- `checks/check_same_gram.py`, `checks/check_split_halves.py`: two exact identity checks used in Sections 5.5 and 5.6.
- `make_numbers_v6.py`: every number and table of the paper, from released files only, and `VERIFY_numbers.md`.

Every output records the SHA-256 of the programs that produced it (`code_sha256`).

## Rerunning every number of the paper

You need Python 3.9 with NumPy 2.0.2, SciPy 1.13.1, python-flint 0.6.0 and mpmath 1.3.0, and the data set:

    zsh certify/census/fetch_snapshot.sh SNAPDIR          # downloads table.txt, data-1..12.txt, checks every MD5
    python3 -I certify/census/snapshot.py SNAPDIR         # data.txt and census/index.tsv

The census (about 4 CPU hours, then 1 for the strictness stage; outputs are released, so this is optional):

    zsh certify/census/run_all.sh SNAPDIR classify
    zsh certify/census/run_all.sh SNAPDIR certify
    zsh certify/census/run_all.sh SNAPDIR strict
    zsh certify/census/run_all.sh SNAPDIR exact
    zsh certify/census/run_all.sh SNAPDIR post
    zsh certify/v5/run_v5.sh SNAPDIR OUTDIR

Version 6 (about one CPU hour) and the numbers of the paper:

    export QV_CERTIFY=$PWD/certify QV_SNAP=SNAPDIR
    zsh certify/paper_v6/run_v6.sh SNAPDIR               # survey, diag, sym, controls, red, checks, numbers

`make_numbers_v6.py` writes `v6_numbers.tex` and the tables `v6_tab_*.tex`, which the paper's source inputs; every
number of the paper comes from these files. One certificate alone:

    python3 -I certify/census/cq.py CODE.txt D N --tol 1e-10 [--flex] [--json OUT.json]
    python3 -I certify/paper_v6/phase2/symstrict.py OUTDIR --only D_N

The floating-point steps can differ in the last bits between machines, so a rerun can give slightly different
enclosures; whatever is reported as proved is proved either way.

## The trusted kernel

The pseudocode of the floating-point kernel (steps K1-K5, S1-S3, P, G, Q, F) is Appendix A of the paper; it is part of
the proofs. `up(s) = fl(fl(s(1 + 2^-29)) + 2^-996)` bounds from above any sum of at most 2^20 nonnegative products
computed in any order (Lemma 3.19(a)).

## Files: git and Zenodo

Everything the paper uses is in git: the programs, `census/index.tsv`, `census/snapshot_files.json`, the per-entry
outputs `census/results_v4/` (16,265 files), `census/strict_v4/`, `census/exact_v4/`, `census/summary_v4.json`,
`census/validate_v4.tsv`, `census/classify.jsonl`, `census/MANIFEST_v4.tsv`, `census/logs/*v4*`, all of `v5/`, and all of
`paper_v6/`. The Zenodo deposit **[DOI]** holds a copy of the directory at the version of the paper. Superseded version
3 outputs (see `.gitignore`) go nowhere.

## Earlier versions

`run_batch.py`, `analysis.py`, `exact_check.py`, `lp_dual.py`, `make_lists.py`, `dedup.py`, `lists/` and `results/` at this
level belong to version 2 of the paper (certificates of the codes we uploaded) and are kept for that record. Version 3
outputs are superseded by version 4 and are not part of the release.

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
