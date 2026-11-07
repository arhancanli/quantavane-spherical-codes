#!/bin/zsh
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
# Rerun everything behind the note "A jammed spherical code that is not infinitesimally jammed".
#   zsh certify/note/run_note.sh SNAPDIR NOTEDIR
# SNAPDIR: the archived table (the directory holding data.txt, DSpace@MIT handle 1721.1/153543); NOTEDIR: where the LaTeX
# numbers and tables go. About twenty-five minutes on one core (on the author's loaded machine); peak resident memory about
# 2.3 GB (PARI's stack may grow to 4 GB, run_gp.py).
#
# Environment (see requirements_note.txt): Python 3.9 or later with numpy, scipy, mpmath, python-flint 0.6.0 (for the
# version 5 program certify/v5/jam116.py) and cypari 2.5.7, which bundles PARI/GP 2.15.4 (for the GP scripts). The
# scripts run in isolated mode (python3 -I), which ignores PYTHONPATH; extra library folders can be given in QV_LIBS
# (colon-separated), and PYTHON selects the interpreter. The defaults below are the author's machine; with the packages
# installed in the interpreter itself, run as   QV_LIBS= PYTHON=python3 zsh certify/note/run_note.sh SNAPDIR NOTEDIR
#
# The floating-point computations that take hours are run separately: the count across the table (symcount.py, by
# run_symcount.sh) and the searches of Appendix B (family.py, family_table.py, by run_family.sh); this script only reads
# their outputs in out/symcount/ and out/family/ (symcount_summary.py, family_summary.py).
set -e
HERE=${0:A:h}
SNAP=$1
NOTE=$2
CACHE=$HOME/Library/Caches/quantavane-census
export QV_LIBS=${QV_LIBS-$CACHE/pylib:$CACHE/parilib:$CACHE/flintlib}
PY=(nice -n 10 ${PYTHON:-/usr/bin/python3} -I)
cd $HERE
mkdir -p out data_xprime
# the code X* (Theorem 1.1)
$PY extract.py $SNAP                                   > out/extract.log
$PY classes.py ../v5/out/recognition_d22_n116.json     > out/classes.log
$PY aut116.py                                          > out/aut116.log
$PY construct.py                                       > out/construct.log
$PY orbit_model.py                                     > out/orbit_model.log
$PY run_gp.py structure.gp structure flex.gp stress.gp > out/run_structure.txt 2>&1
$PY run_gp.py check.gp check                           > out/run_check.txt 2>&1
(cd ../v5 && $PY jam116.py $SNAP out/recognition_d22_n116.json ../note/out/jam116_rerun.json > ../note/out/jam116_rerun.log)
$PY run_gp.py field.gp field                           > out/run_field.txt 2>&1
# the second code X' (Theorem 7.1): X* with Q reflected in the hyperplane of A
$PY run_gp.py xprime.gp xprime                         > out/run_xprime.txt 2>&1
QV_NOTE_DATA=data_xprime $PY run_gp.py structure.gp structure_xprime flex.gp stress.gp > out/run_structure_xprime.txt 2>&1
$PY galois_xprime.py                                   > out/galois_xprime.log
QV_NOTE_DATA=data_xprime $PY run_gp.py check_xprime.gp check_xprime                    > out/run_check_xprime.txt 2>&1
# the model of Definition 3.1, built a second time by a separately written program, against data/ and data_xprime/
$PY model_check.py out/model_check.json                > out/model_check.log
# the second, independent check of X' (Python/Arb: an adapted copy of v5/jam116.py), and two red controls that must fail
$PY jam_xprime.py out/jam_xprime.json                  > out/jam_xprime.log
for m in gram ratio; do
  if $PY jam_xprime.py out/jam_xprime_red_$m.json --red $m > out/jam_xprime_red_$m.log 2>&1; then
    echo "RED CONTROL $m DID NOT FAIL"; exit 1
  fi
done
echo "RED CONTROLS: 2 of 2 failed as they must (gram: (G) fails; ratio: (H) fails)" > out/red_controls.txt
# numbers and tables
$PY optima.py $SNAP out/optima.json                    > out/optima.log
$PY summarize.py $NOTE                                 > out/summary.log
$PY family_summary.py $NOTE                            > out/family_summary.log
$PY symcount_summary.py $NOTE                          > out/symcount_summary.log
$PY manifest.py                                        > out/manifest.log
grep -h STATUS out/check.log out/check_xprime.log out/field.log
grep -ho '"status": "[^"]*"' out/model_check.json
cat out/galois_xprime.log
grep -ho '"status": "[^"]*"' out/jam116_rerun.json out/jam_xprime.json
cat out/red_controls.txt
