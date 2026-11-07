#!/bin/zsh
# The census, end to end (paper version 4). Light background jobs: nice 10, at most 4 processes, one BLAS thread each.
# SNAP is the folder with the downloaded data set (table.txt, data-1.txt ... data-12.txt).
#   zsh certify/census/run_all.sh SNAP [classify|certify|strict|exact|post|all]
set -e
SNAP=${1:?snapshot folder}
STEP=${2:-all}
C=$(cd "$(dirname "$0")" && pwd)
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
PY=(nice -n 10 /usr/bin/python3 -I)   # an array: zsh does not split a string
mkdir -p $C/logs
stamp() { echo "$(date -u +%FT%TZ) $1" >> $C/logs/run_all_v4.log }
if [[ $STEP == all || $STEP == classify ]]; then
  stamp "classify start"
  [[ -f $C/index.tsv ]] || $PY $C/snapshot.py $SNAP
  for k in 0 1 2 3; do $PY $C/classify.py $SNAP $C/classify_$k.jsonl --shard $k/4 > $C/logs/classify_$k.log 2>&1 & done
  wait
  cat $C/classify_[0-3].jsonl > $C/classify.jsonl
  stamp "classify done"
fi
if [[ $STEP == all || $STEP == certify ]]; then
  stamp "certify start"
  for k in 0 1 2 3; do $PY $C/run_census.py $SNAP $C/classify.jsonl $C/results_v4 --shard $k/4 --mem-limit-mb 1800 \
      > $C/logs/census_v4_$k.log 2>&1 & done
  wait
  stamp "certify done"
fi
if [[ $STEP == all || $STEP == strict ]]; then
  stamp "strict start"
  for k in 0 1 2 3; do $PY $C/strict.py $SNAP $C/results_v4 $C/strict_v4 --shard $k/4 > $C/logs/strict_v4_$k.log 2>&1 & done
  wait
  stamp "strict done"
fi
if [[ $STEP == all || $STEP == exact ]]; then
  stamp "exact start"
  zsh $C/run_exact.sh $SNAP
  stamp "exact done"
fi
if [[ $STEP == all || $STEP == post ]]; then
  $PY $C/validate_poly.py $C/results_v4 $C/validate_v4.tsv > $C/logs/validate_v4.log 2>&1
  $PY $C/summarize.py $C ${PAPER:-$C/paper_out} > $C/logs/summarize_v4.log 2>&1
  $PY $C/make_text.py $C ${PAPER:-$C/paper_out} >> $C/logs/summarize_v4.log 2>&1
  $PY $C/census_md.py $C > ${PAPER:-$C/paper_out}/CENSUS.md
  stamp "post done"
fi
