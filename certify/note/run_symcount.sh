#!/bin/zsh
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Arhan Canli (Quantavane)
# The screening of the archived table for symmetry-forced flexes (Section 4 of the note; floating point, hours).
#   zsh certify/note/run_symcount.sh SNAPDIR
# Three shards of symcount.py in parallel (nice 10, one BLAS thread each, peak resident memory below 1.5 GB each); a
# shard that is stopped resumes where it stopped. run_note.sh reads the outputs (symcount_summary.py).
set -e
HERE=${0:A:h}
SNAP=$1
CACHE=$HOME/Library/Caches/quantavane-census
export QV_LIBS=${QV_LIBS-$CACHE/pylib:$CACHE/parilib:$CACHE/flintlib}
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
PY=(nice -n 10 ${PYTHON:-/usr/bin/python3} -I)
cd $HERE
mkdir -p out/symcount
for i in 0 1 2; do
  $PY symcount.py $SNAP out/symcount/symcount_$i.jsonl --shard $i/3 > out/symcount/symcount_$i.log 2>&1 &
done
wait
wc -l out/symcount/symcount_*.jsonl
