# Spherical codes uploaded to Henry Cohn's table

1,039 spherical codes (N points in dimension 3 to 32) whose largest cosine was below the value listed on Henry Cohn's
table of spherical codes, spherical-codes.org, when we uploaded them (October and November 2026). On 1 November 2026
the table showed 414 of them as the best known code for their entry; for most of the others the table's own polishing
pass has since listed a slightly better code.

- `codes/dD_nN.txt`: our best checked code for each uploaded entry, one point per line (comma-separated coordinates,
  the table's format); the first line gives D, N and the largest cosine.
- `checker/sph_check.py`: an independent checker (standard library only, 50-digit decimal arithmetic) that normalises
  every point and reports the largest cosine: `python3 checker/sph_check.py codes/d10_n147.txt 10 147`.
- `records/`: the upload records (UPLOADED.tsv: each upload the table served; IMPROVED.tsv: every checked code that
  beat the table of 25 October 2026; BEATEN.tsv, LIVE_SEEN.tsv, REPLACED.tsv), the table as read on 25 October
  2026 (targets_20261025.csv) and on 1 November 2026 (live_table.tsv), and make_table.py, the script that computes
  every number in the paper from these files.

Paper: A. Canli, Improved spherical codes for Cohn's table (in preparation).

By Arhan Canli (Quantavane). CC BY 4.0 (see LICENSE). Corrections are welcome as GitHub issues.
