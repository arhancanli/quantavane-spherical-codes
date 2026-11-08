#!/bin/zsh
# Census version 8: every run behind paper_v8, end to end, from a clean checkout of quantavane-spherical-codes.
# Light background jobs: nice 10, one BLAS thread, at most three processes at a time. Needs the archived data set (SNAP:
# the folder with table.txt, data.txt, data-1..12.txt) and the released census outputs in certify/ (QV_CERTIFY, default
# ~/quantavane-spherical-codes/certify). Everything else this script reads is in this directory.
#   zsh run_v8.sh SNAP [exact|strict|flex8|above|tolscan|zeta|gapscan|gapcert|example|live|numbers|all]
#   zsh run_v8.sh SNAP livefetch     (fetches about 800 codes from the live web site; not part of 'all')
set -e
SNAP=${1:?snapshot folder}
STEP=${2:-all}
H=$(cd "$(dirname "$0")" && pwd)
export QV_SNAP=$SNAP OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
CERT=${QV_CERTIFY:-$HOME/quantavane-spherical-codes/certify}
export QV_CERTIFY=$CERT
PY=(nice -n 10 python3 -I)
mkdir -p $H/verify/logs $H/recount/logs $H/live/logs $H/above/logs
if [[ $STEP == all || $STEP == exact ]]; then           # the decision of record of every positive certificate
  $PY $H/verify/exactcheck.py $H/verify/out_exact --shard 0/2 > $H/verify/logs/exact_0.log 2>&1 &
  $PY $H/verify/exactcheck.py $H/verify/out_exact --shard 1/2 > $H/verify/logs/exact_1.log 2>&1 &
  wait
  # version 2 (the second stress enclosure) where version 1 left a stress sign undecided
  $PY $H/verify/exactcheck2.py $H/verify/out_exact2 --only 3_613,3_877 --status PASS > $H/verify/logs/exact2_pass.log 2>&1
  $PY $H/verify/exactcheck2.py $H/verify/out_exact2 --only 3_791 --status NONOPT > $H/verify/logs/exact2_nonopt.log 2>&1
fi
if [[ $STEP == all || $STEP == strict ]]; then          # strictness certificates with exact (S') and span
  ONLY=$(python3 -I -c "
import json
s = json.load(open('$CERT/v5/out/summary_v5.json'))['class']
print(','.join(k.replace(',', '_') for k, v in s.items() if v[1] == 'strict_hyper'))")
  $PY $H/verify/strict7.py $SNAP $CERT/census/results_v4 $H/verify/out_strict7 --only $ONLY > $H/verify/logs/strict7.log 2>&1
  ONLY=$(python3 -I -c "
import json, glob
print(','.join(sorted(f\"{j['d']}_{j['n']}\" for j in (json.load(open(f)) for f in glob.glob('$H/phase2/out_open/*.json')) if j['status'] == 'STRICT')))")
  $PY $H/verify/symstrict7.py $H/verify/out_sym7 --entries open --only $ONLY > $H/verify/logs/sym7.log 2>&1
  # repeated with version 2 after the untrusted linear program of the first run reached its time limit
  $PY $H/verify/strict7b.py $SNAP $CERT/census/results_v4 $H/verify/out_strict7b --only 19_102 > $H/verify/logs/strict7b.log 2>&1
fi
if [[ $STEP == all || $STEP == flex8 ]]; then           # the negative certificate, exactly, for every entry
  $PY $H/recount/flexgram7.py $H/recount/flexgram7_0.jsonl --shard 0/2 > $H/recount/logs/flexgram7_0.log 2>&1 &
  $PY $H/recount/flexgram7.py $H/recount/flexgram7_1.jsonl --shard 1/2 > $H/recount/logs/flexgram7_1.log 2>&1 &
  wait
  $PY $H/recount/flexexact8.py $H/recount/out_flex8 > $H/recount/logs/flexexact8_0.log 2>&1 &
  $PY $H/recount/flexexact8.py $H/recount/out_flex8 > $H/recount/logs/flexexact8_1.log 2>&1 &
  wait
  # a third right inverse (singular value decomposition) where the first two fail the exact test
  $PY $H/recount/flexrefine8.py $H/recount/out_flex8 $H/recount/out_flex8r > $H/recount/logs/flexrefine8.log 2>&1
fi
if [[ $STEP == all || $STEP == above ]]; then           # the structured codes above the size limit, exactly
  $PY $H/above/run_above.py $H/above/out > $H/above/logs/above.log 2>&1
  # a second pass with a larger cap on the classes of contacts, where no stress was found
  $PY $H/above/above_more.py $H/above/out $H/above/out_more > $H/above/logs/above_more.log 2>&1
  $PY $H/above/gan510.py $H/above/gan510.json > $H/above/logs/gan510.log 2>&1
fi
if [[ $STEP == all || $STEP == tolscan ]]; then         # the other tolerances and the close calls
  $PY $H/recount/tolscan.py other $H/recount/tolscan_other_0.jsonl --shard 0/2 > $H/recount/logs/tolscan_other_0.log 2>&1 &
  $PY $H/recount/tolscan.py other $H/recount/tolscan_other_1.jsonl --shard 1/2 > $H/recount/logs/tolscan_other_1.log 2>&1 &
  wait
  $PY $H/recount/tolscan.py close $H/recount/tolscan_close.jsonl > $H/recount/logs/tolscan_close.log 2>&1
fi
if [[ $STEP == all || $STEP == zeta ]]; then            # a sharper rigidity constant (untrusted)
  $PY $H/recount/zetasample.py $H/recount/zetasample.jsonl > $H/recount/logs/zetasample.log 2>&1
fi
if [[ $STEP == all || $STEP == gapscan ]]; then         # gap-based proposals (untrusted scan)
  $PY $H/recount/gapcut.py scan $H/recount/gapcut_0.jsonl --shard 0/1 > $H/recount/logs/gapcut_0.log 2>&1
fi
if [[ $STEP == all || $STEP == gapcert ]]; then         # certificates on a sample of the gap-based proposals
  $PY $H/recount/gapsample.py $H/recount/gapcut_0.jsonl $H/recount/gapcut_sample.tsv
  $PY $H/recount/gapcut.py certify $H/recount/gapcut_sample.tsv $H/recount/gapcut_cert.jsonl > $H/recount/logs/gapcut_cert.log 2>&1
fi
if [[ $STEP == all || $STEP == example ]]; then         # the worked examples of Section 5
  $PY $H/example/example.py 3 7 $H/example/example_d3_n7.json > $H/example/example_d3_n7.log 2>&1
  $PY $H/example/example.py 5 16 $H/example/example_d5_n16.json > $H/example/example_d5_n16.log 2>&1
fi
if [[ $STEP == livefetch ]]; then                       # current codes from the live table (network)
  FETCH=${QV_LIVEFETCH:?folder for the fetched codes}
  $PY $H/live/select_live.py $CERT/v5/live/live_20261107.tsv $H/recount/flexsample.tsv.categories $H/live/live_select.tsv
  $PY $H/live/fetch_live.py $H/live/live_select.tsv $FETCH > $H/live/logs/fetch.log 2>&1
fi
if [[ $STEP == all || $STEP == live ]]; then            # the census rule and the exact decision on the fetched codes
  if [[ -n ${QV_LIVEFETCH:-} ]]; then
    $PY $H/live/run_live.py $QV_LIVEFETCH $H/live/live_select.tsv $H/live/out --shard 0/2 > $H/live/logs/live_0.log 2>&1 &
    $PY $H/live/run_live.py $QV_LIVEFETCH $H/live/live_select.tsv $H/live/out --shard 1/2 > $H/live/logs/live_1.log 2>&1 &
    wait
  else
    echo "live: QV_LIVEFETCH not set; using the released outputs in live/out"
  fi
fi
if [[ $STEP == all || $STEP == numbers ]]; then         # every number, table and figure, and VERIFY_numbers_v8.md
  $PY $H/make_numbers_v8.py
fi
