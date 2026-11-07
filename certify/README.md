# Certified local optimality in a table of spherical codes: the programs

The programs and outputs of the paper *Certified local optimality in a table of spherical codes, and a jammed code that
is not infinitesimally jammed* (A. Canli, version 5, November 2026). They apply to the archived version of Cohn's table
of spherical codes (DSpace@MIT handle 1721.1/153543, data as of 31 July 2026), entries with 3 <= d <= 32 and N <= 1024.

**Licence.** Everything in this directory (`certify/`) is released under the MIT licence (`certify/LICENSE`,
Copyright (c) 2026 Arhan Canli (Quantavane)). The table's data remain Cohn's; the repository's other directories have
their own licence (`LICENSE` at the top).

**What is trusted.** Every function is marked TRUSTED (its result enters a decision, with exact integer or rational
arithmetic, arithmetic in a number field, ball arithmetic with rigorous radii, or floating point with rigorous error
bounds formed only from individual products) or UNTRUSTED (it only proposes candidates: contacts, rattlers, Newton steps,
approximate inverses, subsets, coefficients of stresses, recognised numbers, colourings, flex patterns). Whatever a
program reports as proved is proved whichever candidates it was given.

## What each program proves

Version 4 (the census, `census/`):

- `cq.py` + `sphcert.py`: `PASS` (isostatic graph: an exact configuration with exactly the proposed contacts exists within
  the printed distance, a strict local minimum, explicit constants), `LOCALMIN` (hyperstatic graph: a local minimum
  within the printed distance), `NONOPT` (the exact equal-contact configuration is not a local minimum), `FLEX` (the
  diagnostic of Proposition 4.2: no configuration within the printed l1 radius carries a nonzero stress on the pairs
  within tau; a statement about a neighbourhood, not about an exact configuration).
- `strict.py`: `STRICT` on a LOCALMIN entry (the Gram matrix recognised in Q(c*) and verified exactly).
- `exact_strict.py`: the classical criterion with no floating point in any decision, for exact coordinates or a rational
  Gram matrix. `exact_big.py`: a second program (exact integers and rationals) that re-checks a sample.

Version 5 (`v5/`):

- `core.py`: the exact contact core of an entry: the core's Gram matrix recognised in Q or Q(sqrt D) (untrusted) and
  verified exactly with its rank, the rattlers kept at their stored positions and proved to touch nothing (ball
  arithmetic), and the certified distance to the stored code.
- `subspace.py`: `NOT A LOCAL MINIMUM (flat core ...)`: the core spans a proper subspace and the cosine is positive
  (Proposition 2.10); one explicit tilted configuration is checked below c.
- `rigidity.py`: `STRICT` (bar-rigid core with an exact strictly positive stress; Theorem 2.5) or `LOCAL, NOT STRICT`
  (split contact graph, each component bar-rigid with a positive stress; Proposition 2.11).
- `recognise116.py` (untrusted) and `jam116.py` (trusted): the entry (22,116) is jammed but not infinitesimally jammed
  (Theorem 5.1): exact arithmetic in the field Q(sqrt145, s1, s2) of degree 8 (`tower.py`).
- `distances.py`: certified distances for the exact-only confirmations of version 4.
- `maxcos_exact.py`: the exact largest cosine of a code given by decimals. `families.py`: numerical (untrusted) evidence
  on the continuous families. `summarize_v5.py`: every version 5 number and table of the paper. `release_list.py`: the
  release plan and `v5/MANIFEST_v5.tsv`.
- `kfield.py` (Q(sqrt D)), `tower.py` (a tower of quadratic fields), `dist.py` (exactly orthogonal rational rotations),
  `exact2.py` (shared helpers: recognition, ranks over F_p, the second-order form).

Every output records the SHA-256 of the programs that produced it (`code_sha256`).

## Rerunning every number of the paper

You need Python 3.9 with NumPy 2.0.2, SciPy 1.13.1, python-flint 0.6.0 and mpmath 1.3.0 (the paths used here are set at
the top of each program; edit them for your installation), and the data set:

    zsh certify/census/fetch_snapshot.sh SNAPDIR          # downloads table.txt, data-1..12.txt, checks every MD5
    python3 -I certify/census/snapshot.py SNAPDIR         # data.txt and census/index.tsv

Then:

    zsh certify/census/run_all.sh SNAPDIR classify        # contacts and rattlers at every tolerance (classify.jsonl)
    zsh certify/census/run_all.sh SNAPDIR certify         # the census certificates (results_v4/), about 4 CPU hours
    zsh certify/census/run_all.sh SNAPDIR strict          # the strictness stage (strict_v4/)
    zsh certify/census/run_all.sh SNAPDIR exact           # exact_strict.py and the exact_big.py sample (exact_v4/)
    PAPER=path/to/paper zsh certify/census/run_all.sh SNAPDIR post    # v4_numbers.tex, v4_tab_*.tex, validate_v4.tsv
    zsh certify/v5/run_v5.sh SNAPDIR path/to/paper        # everything of version 5 (about half an hour) and
                                                          # v5_numbers.tex, v5_tab_*.tex, RELEASE_FILES.tsv

The paper's LaTeX source inputs `v4_numbers.tex`, `v4_text.tex`, `v4_tab_strict.tex`, `v5_numbers.tex` and the
`v5_tab_*.tex` tables; every number in it comes from these files. One certificate alone:

    python3 -I certify/census/cq.py CODE.txt D N --tol 1e-10 [--flex] [--json OUT.json]
    python3 -I certify/v5/subspace.py SNAPDIR D N OUT.json [--disc D]
    python3 -I certify/v5/rigidity.py SNAPDIR D N OUT.json [--disc D]
    python3 -I certify/v5/recognise116.py SNAPDIR REC.json && python3 -I certify/v5/jam116.py SNAPDIR REC.json OUT.json

The floating-point steps can differ in the last bits between machines, so a rerun can give slightly different
enclosures; whatever is reported as proved is proved either way. NumPy built against Apple's Accelerate may print
"divide by zero / overflow / invalid value encountered in matmul" warnings from untrusted matrix-vector products; they
come from the library's handling of subnormal intermediate values in code whose results are only candidates, and they
do not affect any decision.

## The trusted kernel (the paper, Section 3.7)

`up(s) = fl(fl(s(1 + 2^-29)) + 2^-996)` bounds from above any sum of at most 2^20 nonnegative products computed in any
order (the paper's Lemma 3.14(a)); part (b) bounds the error of a computed sum of products by gamma_k sum|a_i b_i| + k eta0.
"sum" is a floating-point sum in any order; arrays are in double precision unless marked exact; C is processed in blocks
of 128 rows or columns.

- **K1** Exact: F(z~) as rationals. F_f = round to nearest; F_err = |F_f| 2^-52 + eta0. z^ = fl(z~),
  delta = |z^| 2^-52 + eta0, J0 = DF(z^) (exact).
- **K2** u_J = up(row sums of |J0|); a_Delta = up(Delta(delta) 1); b_Delta = Delta(1) 1 (2d or 0, exact).
- **K3** For each row i of C: Ybar_i = up(|fl(sum_j C_ij F_f,j)| + up(gamma_{U+1} up(sum_j |C_ij||F_f,j|)) +
  up(sum_j |C_ij| F_err,j)) [(b) for the first term]; P_i = fl(C_i J0) (sparse), E^_i = fl(e_i - P_i);
  a_i = up(up(up(sum_j |E^_ij|)) + up(gamma_{U+1} up(sum_j |C_ij| u_J,j)) + up(sum_j |C_ij| a_Delta,j))
  [(b) for the rounding of P]; b_i = up(2d up(sum over the rows F_i, F_ij of |C_ij|));
  nu_i = up(sum_j |C_ij|), nuE_i = up(sum over the rows F_ij, F_st of |C_ij|).
- **K4** r = up(2 max Ybar); beta_i = up(a_i + up(r b_i)); eps_i = up(Ybar_i + up(beta_i r)). Require eps_i < r for all i
  (otherwise retry with r = up(4 max Ybar / (1 - max beta)), at most three times). This is (K); beta_K = max beta_i.
- **K5** (epsilon-inflation) For r' = q/2, q/4, ... with q ~ min_i (1 - a_i)/b_i: the first r' that passes K4's test is r_u.
- **S1** y~ ~ J0^-T b (untrusted), b = -e_gamma, or for a hyperstatic graph b = -e_gamma - sum_f lambda_f grad F_f(p^)
  with a radius b_rad covering its rounding and the enclosure of p. Delta_r = Delta(up(delta + r)) (sparse).
  Res = up(up(up(|fl(b - J0^T y~)|) + up(gamma_{U+1} up(|J0|^T |y~|))) + up(Delta_r^T |y~|) + b_rad) [(b)].
- **S2** For each column j of C: w_j = up(sum_i |C_ij| Res_i); Q_j = fl(J0 C_j), E^ = fl(e_j - Q_j);
  s_j = up(up(up(sum_i |E^_ij|)) + up(gamma_{U+1} up(sum_i v_J,i |C_ij|)) + up(sum_i v_Delta,i |C_ij|)) with v_J, v_Delta
  the upper bounds of the column sums of |J0| and Delta_r.
- **S3** beta_S = max s, require beta_S < 1; wbar = up(up(max w) / ((1 - beta_S)(1 - 2^-50))); s = up(w + up(s wbar)).
  Require y~_e > s_e for every contact e (taken as the next float below y~_e - s_e). This is (S); omega_- is the smallest
  of these lower bounds (and of the lambda_f).
- **P** Choose d points by pivoted QR (untrusted), B their rows of z^, C_B ~ B^-1; require every row sum of
  up(up(|fl(I - C_B B)|) + up(gamma_d up(|C_B||B|)) + up(|C_B| (delta + eps)_B)) to be below 1.
- **G** Every inner product is bounded in floating point by (b) with the distances eps; a pair not placed below
  gamma~ - eps_gamma is decided in rational arithmetic; a rattler enters only through the exact rational |Y_r|^2.
- **Q** mu = up(max nu / ((1 - beta_K)(1 - 2^-50))), mu_p = max_{i <= nd} up(nuE_i + up(beta_i mu)),
  delta_xi = up(max |p^ - xi| + max(delta + eps)), theta = up(2 n delta_xi mu_p); if theta < 1: sigma, kappa, rho of
  Theorem 3.6 rounded down; D(X_V) <= up(sqrt(sum(|p^ - xi| + e_X + delta + eps)^2)), rounded up, where e_X bounds the
  distance between xi and the exactly normalised stored code.
- **F** (Proposition 4.2) L exact in floating point; for each block of columns of C: P = fl(L C_blk) (sparse),
  E^ = fl(I - P); the row sums of |I - LC| are bounded by up(up(sum |E^|)) + up(gamma_2d up(|L| up(|C| 1))) [(b)];
  ||C||_inf <= max up(|C| 1); h rounded down.

## Files: git and Zenodo

`RELEASE_FILES.tsv` lists every file of `certify/` with its size and destination. Programs, documentation and the small
outputs go to git (under 20 MB in all): `census/strict_v4/`, `census/exact_v4/`, `census/index.tsv`,
`census/snapshot_files.json`, `census/summary_v4.json`, `census/validate_v4.tsv`, `census/MANIFEST_v4.tsv`,
`census/logs/*v4*`, and all of `v5/` (programs, `v5/out/`, `v5/live/`, `v5/MANIFEST_v5.tsv`). The one-file-per-entry
census outputs `census/results_v4/` (16,265 files, about 23 MB), the proposal file `census/classify.jsonl` (16 MB) and the
outputs of version 2 of the paper (`results/`, `lists/`, `logs/` at this level) go to the Zenodo deposit
**[DOI to be inserted by the author]**, which also holds an archive of the whole directory at the version of the paper.
Superseded version 3 outputs (see `.gitignore`) go nowhere.

## Earlier versions

`run_batch.py`, `analysis.py`, `exact_check.py`, `lp_dual.py`, `make_lists.py`, `dedup.py`, `lists/` and `results/` at this
level belong to version 2 of the paper (certificates of the codes we uploaded) and are kept for that record. Version 3
outputs are superseded by version 4 and are not part of the release.

## The note: a jammed spherical code that is not infinitesimally jammed (`note/`)

The programs of the note *A jammed spherical code that is not infinitesimally jammed* (A. Canli). Rerun everything with

    zsh certify/note/run_note.sh SNAPDIR NOTEDIR     # about fifteen minutes on one core

which runs, in order: `extract.py` (the stored entry), `classes.py` (classes and the 13 exact values, from the version 5
recognition), `aut116.py` (the group by individualisation and refinement), `construct.py` (the model in the circle
geometry of PG(1,9)), `orbit_model.py` (each orbit on pairs described in the model; Tables 1 and 5), `structure.gp` (the
flex, the stresses, the characters), `check.gp` (the second exact verification, including the count of the 1440
automorphisms by backtracking: 147,129 nodes), `v5/jam116.py` (the first exact verification), then for the second code
X' of Section 7 `xprime.gp` (the exact data of X' in `data_xprime/`), `structure.gp` and `check_xprime.gp` on that
directory (`QV_NOTE_DATA=data_xprime`), `jam_xprime.py` (the second, independent check of X': an adapted copy of
`v5/jam116.py` in Python with ball arithmetic, which recomputes X' from X* itself and must agree with `data_xprime/`)
and its two red controls (`--red gram`, `--red ratio`, which must fail: run_note.sh stops if either passes), and finally
`summarize.py`, `family_summary.py` (numbers and tables) and `manifest.py`. The last lines printed are the two GP
STATUS lines, the two Python statuses and the red-control line. On the author's machine (loaded) the GP scripts take two to four minutes each and `jam116.py` about a
minute.

The floating-point searches of Section 7 are run separately (hours): `zsh certify/note/run_family.sh SNAPDIR` runs
`family.py` (Gamma-invariant codes from PGammaL(2,q) and subgroups: the restarts for X* and X', then the cases of its
table `CASES`) and `family_table.py` (the entries of the table in dimension 22 with N <= 100); their outputs are in
`note/out/family/`. Nothing in them is a certificate; the note states which results are floating point.
