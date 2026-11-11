# Spherical codes uploaded to Henry Cohn's table

1,042 spherical codes (N points in dimension 3 to 32) whose largest cosine was below the value listed on Henry Cohn's
table of spherical codes, spherical-codes.org, when we uploaded them (26 October to 2 November 2026). On the page saved
on 2 November 2026 at 18:46 (+04) the table showed 447 of them as the best known code for their entry; for most of the
others the table's own polishing pass has since listed a slightly better code.

- `codes/dD_nN.txt`: for each uploaded entry, the code we uploaded last, byte for byte as posted (one point per line,
  comma-separated coordinates, the table's format). `codes/MANIFEST.tsv` gives each file's dimension, N, largest cosine
  (our checker's 20-digit value, which the table showed after the upload), upload time and SHA-256, and where the file
  was recovered from. For one entry (d14_n143) the posted file was not kept: the folder holds the same code from a
  later filing with the same 20-digit largest cosine.
- `checker/sph_check.py`: an independent checker (standard library only, 50-digit decimal arithmetic) that normalises
  every point and reports the largest cosine: `python3 checker/sph_check.py codes/d10_n147.txt 10 147`.
- `records/`: the upload records (UPLOADED.tsv and UPLOADED.log: each upload the table served, with its time;
  IMPROVED.tsv: every checked code that beat the table of 25 October 2026, with its time; BEATEN.tsv, LIVE_SEEN.tsv,
  REPLACED.tsv), the table as read on 25 October 2026 (targets_20261025.csv) and on 2 November 2026 at 18:46
  (live_table.tsv), and make_table.py, the script that computes every number in the paper from these files:
  `python3 -I records/make_table.py` writes numbers.tex, table.tex, examples.tex and headline.json in records/.
- `observations/`: the scripts and run records behind the paper's section "Observations", with the outputs they give
  (see observations/README.md).
- `certify/`: computer-assisted certificates that the codes lie next to exact strict local optima of the largest
  cosine (interval Newton test with rigorous floating-point error bounds), one per code, with the programs that make
  and analyse them (see certify/README.md).

Paper: A. Canli, Improved spherical codes in dimensions 3 to 32 (in preparation).

By Arhan Canli (Quantavane). CC BY 4.0 (see LICENSE). Corrections are welcome as GitHub issues.

## How to cite

Cite this repository for the codes and the checker, and Henry Cohn's table of spherical codes (spherical-codes.org),
which lists the best known code for each entry. For one code, give its entry as well, for example d4_n310.

```bibtex
@misc{quantavane-spherical-codes,
  author       = {Canli, Arhan},
  title        = {Spherical codes uploaded to Henry Cohn's table, with an independent checker},
  year         = {2026},
  howpublished = {\url{https://github.com/arhancanli/quantavane-spherical-codes}},
  note         = {Quantavane; CC BY 4.0}
}
```

[CITATION.cff](CITATION.cff) gives the same reference for GitHub's "Cite this repository" button.
