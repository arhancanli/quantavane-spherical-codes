#!/bin/zsh
# Census version 7: every run behind paper_v7, end to end. Light background jobs: nice 10, one BLAS thread, at most
# three processes at a time. Needs the archived data set (SNAP: the folder with table.txt, data.txt, data-1..12.txt) and
# the released census outputs in certify/ (QV_CERTIFY, default ~/quantavane-spherical-codes/certify). The version 6
# runs (run_v6.sh) are inputs and are not repeated here.
#   zsh run_v7.sh SNAP [exact|strict|flexgram|flexexact|gapscan|gapcert|example|livefetch|live|numbers|pdf|all]
# 'livefetch' fetches codes from the live web site (about 800 requests, one every 1.5 s); it is not part of 'all'.
set -e
SNAP=${1:?snapshot folder}
STEP=${2:-all}
H=$(cd "$(dirname "$0")" && pwd)
export QV_SNAP=$SNAP OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
CERT=${QV_CERTIFY:-$HOME/quantavane-spherical-codes/certify}
PY=(nice -n 10 /usr/bin/python3 -I)
mkdir -p $H/verify/logs $H/recount/logs $H/live/logs
if [[ $STEP == all || $STEP == exact ]]; then           # the decision of record for every positive certificate
  $PY $H/verify/exactcheck.py $H/verify/out_exact --shard 0/2 > $H/verify/logs/exact_0.log 2>&1 &
  $PY $H/verify/exactcheck.py $H/verify/out_exact --shard 1/2 > $H/verify/logs/exact_1.log 2>&1 &
  wait
  # version 2 (the second stress enclosure) where version 1 left a stress sign undecided
  $PY $H/verify/exactcheck2.py $H/verify/out_exact2 --only 3_613,3_877 --status PASS > $H/verify/logs/exact2_pass.log 2>&1
  $PY $H/verify/exactcheck2.py $H/verify/out_exact2 --only 3_791 --status NONOPT > $H/verify/logs/exact2_nonopt.log 2>&1
fi
if [[ $STEP == all || $STEP == strict ]]; then          # strictness certificates with exact (S') and span
  ONLY=$(/usr/bin/python3 -I -c "
import json
s = json.load(open('$CERT/v5/out/summary_v5.json'))['class']
print(','.join(k.replace(',', '_') for k, v in s.items() if v[1] == 'strict_hyper'))")
  $PY $H/verify/strict7.py $SNAP $CERT/census/results_v4 $H/verify/out_strict7 --only $ONLY > $H/verify/logs/strict7.log 2>&1
  ONLY=$(/usr/bin/python3 -I -c "
import json, glob
print(','.join(sorted(f\"{j['d']}_{j['n']}\" for j in (json.load(open(f)) for f in glob.glob('$H/phase2/out_open/*.json')) if j['status'] == 'STRICT')))")
  $PY $H/verify/symstrict7.py $H/verify/out_sym7 --entries open --only $ONLY > $H/verify/logs/sym7.log 2>&1
  # repeated with version 2 after the untrusted LP of the first run reached its time limit (stress_lp_status 1)
  $PY $H/verify/strict7b.py $SNAP $CERT/census/results_v4 $H/verify/out_strict7b --only 19_102 > $H/verify/logs/strict7b.log 2>&1
fi
if [[ $STEP == all || $STEP == flexgram ]]; then        # tau' of the negative certificate, with the hash of E
  $PY $H/recount/flexgram7.py $H/recount/flexgram7_0.jsonl --shard 0/2 > $H/recount/logs/flexgram7_0.log 2>&1 &
  $PY $H/recount/flexgram7.py $H/recount/flexgram7_1.jsonl --shard 1/2 > $H/recount/logs/flexgram7_1.log 2>&1 &
  wait
fi
if [[ $STEP == all || $STEP == flexexact ]]; then       # the negative certificate in exact arithmetic, on a sample
  $PY $H/recount/flexsample.py $H/recount/flexgram7_0.jsonl,$H/recount/flexgram7_1.jsonl $H/recount/flexsample.tsv
  $PY $H/recount/flexexact.py $H/recount/flexexact.jsonl $H/recount/flexsample.tsv > $H/recount/logs/flexexact.log 2>&1
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
  FETCH=${QV_LIVEFETCH:?folder for the fetched codes}
  $PY $H/live/run_live.py $FETCH $H/live/live_select.tsv $H/live/out --shard 0/2 > $H/live/logs/live_0.log 2>&1 &
  $PY $H/live/run_live.py $FETCH $H/live/live_select.tsv $H/live/out --shard 1/2 > $H/live/logs/live_1.log 2>&1 &
  wait
fi
if [[ $STEP == all || $STEP == numbers ]]; then         # every number, table and figure, and VERIFY_numbers_v7.md
  $PY $H/make_numbers_v7.py
  $PY $H/fill_response_v7.py
fi
if [[ $STEP == all || $STEP == pdf ]]; then
  cd $H
  for f in em_v7 mcom_v7; do pdflatex -interaction=nonstopmode $f.tex > /dev/null; pdflatex -interaction=nonstopmode $f.tex > build_$f.log; done
fi
