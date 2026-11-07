"""census/strict.py with its two floating-point kernel calls replaced by the exact tests (exactpatch2.py: version 2, with the second stress enclosure).
  nice -n 10 python3 -I strict7b.py SNAPDIR RESULTS_DIR OUTDIR [--only d_n,...]"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exactpatch2 as exactpatch          # noqa: F401  (applies the patches)
import strict
strict.SHA.update(exactpatch.SHA)
if __name__ == "__main__":
    strict.main()
