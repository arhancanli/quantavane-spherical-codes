#!/bin/zsh
# Census version 9: every run behind paper_v9, from a clean checkout of quantavane-spherical-codes. It builds on the
# released outputs of version 8 (certify/paper_v8/, unchanged). Light background jobs: nice 10, one BLAS thread.
# Needs the archived data set (SNAP: the folder with table.txt, data.txt) and the repository (QV_CERTIFY, default
# ~/quantavane-spherical-codes/certify). Everything else this script reads is in this directory or in paper_v8/.
#   zsh run_v9.sh SNAP [recog|strict9|pass2|prestress|explore|cut|ranks|pari|certdata|exactdata|example|numbers|all]
set -e
SNAP=${1:?snapshot folder}
STEP=${2:-all}
H=$(cd "$(dirname "$0")" && pwd)
export QV_SNAP=$SNAP OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
CERT=${QV_CERTIFY:-$HOME/quantavane-spherical-codes/certify}
export QV_CERTIFY=$CERT
PY=(nice -n 10 python3 -I -W ignore)
mkdir -p $H/recog/logs $H/prestress/logs $H/prestress/explore $H/cut9/logs $H/ranks/logs $H/certdata/files $H/example
if [[ $STEP == all || $STEP == recog ]]; then           # A1: the survey of the open entries (2,000 bits, degree <= 30)
  $PY $H/recog/recog9.py $H/recog/out_all --shard 0/2 > $H/recog/logs/all_0.log 2>&1 &
  $PY $H/recog/recog9.py $H/recog/out_all --shard 1/2 > $H/recog/logs/all_1.log 2>&1 &
  wait
fi
if [[ $STEP == all || $STEP == strict9 ]]; then         # A1: the exact decision; the controls (version 8 recognitions)
  ONLY=$(python3 -I -c "
import json, glob
print(','.join(f\"{j['d']}_{j['n']}\" for j in (json.load(open(f)) for f in sorted(glob.glob('$H/recog/out_all/*.json'))) if str(j['status']).startswith('recognised')))")
  $PY $H/recog/strict9.py $H/recog/out_strict9 --only $ONLY --data $H/exactdata/open > $H/recog/logs/strict9.log 2>&1
  # the controls: the version 8 certificates by recognition with at most 2,300 unknowns (recog/ctl_plan.json), two shards
  $PY $H/recog/ctl_plan9.py --max-unknowns 2300 > $H/recog/logs/ctl_plan.log 2>&1
  for S in A B; do
    $PY $H/recog/strict9.py $H/recog/out_strict9_ctl --class strict_hyper --only $(cat $H/recog/ctl_shard$S.txt) --data $H/exactdata/v8recog > $H/recog/logs/strict9_ctl_$S.log 2>&1 &
  done
  wait
fi
if [[ $STEP == all || $STEP == pass2 ]]; then           # A1: a second pass at 4,000 bits for entries with <= 20 values
  $PY $H/recog/recog9.py $H/recog/out_pass2 --pass2 $H/recog/out_all --fewvals 20 --bits 4000 --maxdeg 30 > $H/recog/logs/pass2.log 2>&1
  ONLY=$(python3 -I -c "
import json, glob
print(','.join(f\"{j['d']}_{j['n']}\" for j in (json.load(open(f)) for f in sorted(glob.glob('$H/recog/out_pass2/*.json'))) if str(j['status']).startswith('recognised')))")
  if [[ -n $ONLY ]]; then
    $PY $H/recog/strict9.py $H/recog/out_strict9 --only $ONLY --bits 4000 --data $H/exactdata/open > $H/recog/logs/strict9_pass2.log 2>&1
  fi
fi
if [[ $STEP == all || $STEP == prestress ]]; then       # A2: exact flexes and the second-order form
  $PY $H/prestress/prestress9.py $H/prestress/out --only 22_166,11_604 --max-mb 1500 > $H/prestress/logs/prestress.log 2>&1
fi
if [[ $STEP == all || $STEP == explore ]]; then         # A2: untrusted numerical explorations (logs only)
  $PY $H/prestress/explore9_descent.py 22 166 1 > $H/prestress/explore/descent_d22_n166.log 2>&1
  $PY $H/prestress/explore9_classes.py 22 166 1 > $H/prestress/explore/classes_d22_n166.log 2>&1
  $PY $H/prestress/explore9_classes.py 11 604 2 > $H/prestress/explore/classes_d11_n604.log 2>&1
  $PY $H/prestress/explore9_shape.py 11 604 2 > $H/prestress/explore/shape_d11_n604.log 2>&1
fi
if [[ $STEP == all || $STEP == cut ]]; then             # A3: the intrinsic cut (unflagged, unresolved, the random sample)
  for G in unflagged unresolved sample; do              # (--group rest would search the others; not run for the paper)
    $PY $H/cut9/cut9.py $H/cut9/out --group $G > $H/cut9/logs/cut_$G.log 2>&1
  done
fi
if [[ $STEP == all || $STEP == ranks ]]; then           # B: FLINT ranks of the strict entries above the size limit
  # dense matrices up to 1,280 MiB; the job needs several times that (FLINT's rank routine copies the matrix and the
  # build step holds the nonzeros of K, tens of millions, as Python lists: 5.4 GB at 1,277 MiB). The 14 larger entries
  # are recorded as skipped; --max-mb 3000 computes them all on a machine with enough memory.
  $PY $H/ranks/flintrank9.py $H/ranks/out --max-mb 1280 > $H/ranks/logs/ranks.log 2>&1
fi
if [[ $STEP == all || $STEP == pari ]]; then            # B: PARI/GP spot check of the exact descriptions
  $PY $H/ranks/pari_check9.py $H/ranks/pari_check9.json $H/exactdata/open/*.json.gz $(ls $H/exactdata/v8recog/*.json.gz | head -20) > $H/ranks/logs/pari.log 2>&1
fi
if [[ $STEP == all || $STEP == certdata ]]; then        # B: certificate data of the subset and the stand-alone verifier
  $PY $H/certdata/export9.py $H/certdata/files $(cat $H/certdata/subset.txt) > $H/certdata/export.log 2>&1
  mv $H/certdata/files/export_log.json $H/certdata/export_log.json
  $PY $H/certdata/verify_cert9.py $H/certdata/verify.json $H/certdata/files/*.npz > $H/certdata/verify.log 2>&1
fi
if [[ $STEP == all || $STEP == exactdata ]]; then       # A5: the exact cores and the index of the data set
  $PY $H/exactdata/exactdata_core9.py $H/exactdata/core > $H/exactdata/core.log 2>&1
  $PY $H/exactdata/index9.py > $H/exactdata/index.log 2>&1
fi
if [[ $STEP == all || $STEP == example ]]; then         # the worked examples, regenerated with paper_v8/example/example.py
  $PY $CERT/paper_v8/example/example.py 3 7 $H/example/example_d3_n7.json > $H/example/example_d3_n7.log 2>&1
  $PY $CERT/paper_v8/example/example.py 5 16 $H/example/example_d5_n16.json > $H/example/example_d5_n16.log 2>&1
fi
if [[ $STEP == all || $STEP == numbers ]]; then         # every number, table and figure, and VERIFY_numbers_v9.md
  $PY $H/make_numbers_v9.py
fi
