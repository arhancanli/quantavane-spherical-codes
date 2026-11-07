#!/bin/zsh
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
# Rerun everything behind the note "A jammed spherical code that is not infinitesimally jammed".
#   zsh certify/note/run_note.sh SNAPDIR NOTEDIR
# SNAPDIR: the archived table (the directory holding data.txt, DSpace@MIT handle 1721.1/153543); NOTEDIR: where the LaTeX
# numbers and tables go. About 8 minutes on one core.
#
# Environment (see requirements_note.txt): Python 3.9 or later with numpy, mpmath, python-flint 0.6.0 (for the version 5
# program certify/v5/jam116.py) and cypari 2.5.7, which bundles PARI/GP 2.15.4 (for the GP scripts). The scripts run in
# isolated mode (python3 -I), which ignores PYTHONPATH; extra library folders can be given in QV_LIBS (colon-separated),
# and PYTHON selects the interpreter. The defaults below are the author's machine; with the packages installed in the
# interpreter itself, run as   QV_LIBS= PYTHON=python3 zsh certify/note/run_note.sh SNAPDIR NOTEDIR
set -e
HERE=${0:A:h}
SNAP=$1
NOTE=$2
CACHE=$HOME/Library/Caches/quantavane-census
export QV_LIBS=${QV_LIBS-$CACHE/pylib:$CACHE/parilib}
PY=(nice -n 10 ${PYTHON:-/usr/bin/python3} -I)
cd $HERE
mkdir -p out
$PY extract.py $SNAP                                   > out/extract.log
$PY classes.py ../v5/out/recognition_d22_n116.json     > out/classes.log
$PY aut116.py                                          > out/aut116.log
$PY construct.py                                       > out/construct.log
$PY run_gp.py structure.gp structure flex.gp stress.gp > out/run_structure.txt 2>&1
$PY run_gp.py check.gp check                           > out/run_check.txt 2>&1
(cd ../v5 && $PY jam116.py $SNAP out/recognition_d22_n116.json ../note/out/jam116_rerun.json > ../note/out/jam116_rerun.log)
$PY summarize.py $NOTE                                 > out/summary.log
$PY manifest.py                                        > out/manifest.log
grep -h STATUS out/check.log
