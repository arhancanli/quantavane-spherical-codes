#!/bin/zsh
# Census version 6: every run behind paper_v6, end to end. Light background jobs: nice 10, one BLAS thread, at most
# three processes. Needs the archived data set (SNAP: the folder with table.txt, data.txt, data-1..12.txt) and the
# released census outputs in certify/ (QV_CERTIFY, default ~/quantavane-spherical-codes/certify).
#   zsh run_v6.sh SNAP [survey|diag|flex|arb|sym|controls|red|checks|numbers|pdf|all]
set -e
SNAP=${1:?snapshot folder}
STEP=${2:-all}
H=$(cd "$(dirname "$0")" && pwd)
export QV_SNAP=$SNAP OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
PY=(nice -n 10 /usr/bin/python3 -I)
mkdir -p $H/phase2/logs $H/recount/logs
if [[ $STEP == all || $STEP == survey ]]; then          # untrusted: symmetry of the open entries
  for k in 0 1 2; do $PY $H/phase2/sym_survey.py $H/phase2/survey_$k.jsonl --shard $k/3 > $H/phase2/logs/survey_$k.log 2>&1 & done
  wait
fi
if [[ $STEP == all || $STEP == diag ]]; then            # untrusted: what the existence failures are
  $PY $H/recount/diag_fail.py $H/recount/diag_fail.jsonl --shard 0/2 > $H/recount/logs/diag_0.log 2>&1 &
  $PY $H/recount/diag_fail.py $H/recount/diag_fail_1.jsonl --shard 1/2 > $H/recount/logs/diag_1.log 2>&1 &
  wait
fi
if [[ $STEP == all || $STEP == flex ]]; then            # the separation tau' of the negative certificate (Theorem 4.1)
  $PY $H/recount/flexgram.py $H/recount/flexgram_0.jsonl --shard 0/2 > $H/recount/logs/flexgram_0.log 2>&1 &
  $PY $H/recount/flexgram.py $H/recount/flexgram_1.jsonl --shard 1/2 > $H/recount/logs/flexgram_1.log 2>&1 &
  wait
fi
if [[ $STEP == all || $STEP == arb ]]; then             # independent re-check of a sample in Arb ball arithmetic
  $PY $H/recount/arbcheck.py $H/recount/arbcheck.jsonl --max-unknowns 400 --every 10 > $H/recount/logs/arbcheck.log 2>&1
fi
if [[ $STEP == all || $STEP == sym ]]; then             # the symmetry certificate on the open entries
  for k in 0 1; do $PY $H/phase2/symstrict.py $H/phase2/out_open --entries open --shard $k/2 > $H/phase2/logs/symstrict_open_$k.log 2>&1 & done
  wait
fi
if [[ $STEP == all || $STEP == controls ]]; then        # positive control: the entries certified by recognition
  $PY $H/phase2/symstrict.py $H/phase2/out_controls --entries controls > $H/phase2/logs/symstrict_controls.log 2>&1
fi
if [[ $STEP == all || $STEP == red ]]; then             # red controls on every entry certified by symmetry
  ONLY=$(/usr/bin/python3 -I -c "
import json, glob
print(','.join(sorted(f\"{j['d']}_{j['n']}\" for j in (json.load(open(f)) for f in glob.glob('$H/phase2/out_open/*.json')) if j['status'] == 'STRICT')))")
  for kind in model perm swap; do
    $PY $H/phase2/symstrict.py $H/phase2/out_red_$kind --entries open --only $ONLY --red $kind > $H/phase2/logs/red_$kind.log 2>&1 &
  done
  wait
fi
if [[ $STEP == all || $STEP == checks ]]; then          # exact identity checks of version 5 (unchanged)
  $PY $H/checks/check_same_gram.py $SNAP $H/checks/same_gram.json > $H/checks/same_gram.log 2>&1
  $PY $H/checks/check_split_halves.py $SNAP $H/checks/split_halves.json > $H/checks/split_halves.log 2>&1
fi
if [[ $STEP == all || $STEP == numbers ]]; then         # every number and table of the paper, and VERIFY_numbers.md
  $PY $H/make_numbers_v6.py
fi
if [[ $STEP == all || $STEP == pdf ]]; then
  cd $H
  for f in em_v6 mcom_v6; do pdflatex -interaction=nonstopmode $f.tex > /dev/null; pdflatex -interaction=nonstopmode $f.tex > build_$f.log; done
fi
