#!/bin/zsh
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
# The floating-point searches of Section 7 of the note (hours on one core; run them in parallel if you can):
#   zsh certify/note/run_family.sh SNAPDIR
# Writes out/family/*.json and *.log, which family_summary.py (run by run_note.sh) reads. Environment as run_note.sh.
set -e
HERE=${0:A:h}
SNAP=$1
CACHE=$HOME/Library/Caches/quantavane-census
export QV_LIBS=${QV_LIBS-$CACHE/pylib:$CACHE/parilib:$CACHE/flintlib}
PY=(nice -n 10 ${PYTHON:-/usr/bin/python3} -I)
cd $HERE
mkdir -p out/family
$PY family.py out/family/xstar_restarts.json --xstar 100               > out/family/xstar_restarts.log 2>&1
$PY family.py out/family/family_q9.json --cases 0,1,2,3,4              > out/family/family_q9.log 2>&1
$PY family.py out/family/family_q25.json --cases 5,6,7,8,9             > out/family/family_q25.log 2>&1
$PY family.py out/family/family_qother.json --cases 10,11               > out/family/family_qother.log 2>&1
$PY family.py out/family/family_q13.json --cases 12 --starts 4          > out/family/family_q13.log 2>&1
$PY family.py out/family/family_qmore.json --cases 13,14 --starts 4    > out/family/family_qmore.log 2>&1
$PY family.py out/family/family_q16.json --cases 16,15 --starts 4       > out/family/family_q16.log 2>&1
$PY family_table.py $SNAP out/family/table_near.json --dmin 22 --dmax 22 --nmax 100 > out/family/table_near.log 2>&1
