"""paper_v6/symstrict.py with its two floating-point kernel calls replaced by the exact tests (exactpatch2.py: version 2, with the second stress enclosure).
  nice -n 10 python3 -I symstrict7b.py OUTDIR --entries open --only d_n,..."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "phase2"))
import exactpatch2 as exactpatch          # noqa: F401
import symstrict
symstrict.SHA.update(exactpatch.SHA)
if __name__ == "__main__":
    symstrict.main()
