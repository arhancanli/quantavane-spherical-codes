#!/bin/zsh
# The independent exact re-checks of the paper (version 4), one after the other (nice 10, one BLAS thread):
# exact_strict.py on every hyperstatic certificate whose coordinates are exact, then exact_big.py on the sample of
# exact_sample.py (every certificate type).
#   zsh certify/census/run_exact.sh SNAP
SNAP=${1:?snapshot folder}
C=$(cd "$(dirname "$0")" && pwd)
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
PY=(nice -n 10 /usr/bin/python3 -I)   # an array: zsh does not split a string
mkdir -p $C/exact_v4/big $C/logs
$PY $C/exact_strict.py $SNAP $C/results_v4 $C/exact_v4/strict --with-optimal >> $C/logs/exact_v4.log 2>&1
[[ -f $C/exact_v4/sample.tsv ]] || $PY $C/exact_sample.py $SNAP $C > $C/exact_v4/sample.tsv
while IFS=$'\t' read -r mode d n tol cfile; do
  out=$C/exact_v4/big/${mode}_d${d}_n${n}.json
  [[ -f $out ]] && continue
  flag=""; [[ $mode != iso ]] && flag="--$mode"
  $PY $C/exact_big.py $cfile $d $n $tol ${=flag} --json $out >> $C/logs/exact_v4.log 2>&1
done < $C/exact_v4/sample.tsv
echo "$(date -u +%FT%TZ) exact done" >> $C/logs/exact_v4.log
