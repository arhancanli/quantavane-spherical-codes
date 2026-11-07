#!/bin/zsh
# Version 2 of the exact decision (exactcheck2.py: the second stress enclosure) on the certificates for which the first
# enclosure of exactcheck.py left a sign undecided. Outputs in out_exact2/ and out_strict7b/.
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
H=$(cd "$(dirname "$0")" && pwd)
SNAP=${QV_SNAP:-$HOME/cohn-table-dspace-153543/files}
CERT=${QV_CERTIFY:-$HOME/quantavane-spherical-codes/certify}
nice -n 10 /usr/bin/python3 -I $H/exactcheck2.py $H/out_exact2 --only 3_613,3_877 --status PASS > $H/logs/exact2_pass.log 2>&1
nice -n 10 /usr/bin/python3 -I $H/exactcheck2.py $H/out_exact2 --only 3_791 --status NONOPT > $H/logs/exact2_nonopt.log 2>&1
nice -n 10 /usr/bin/python3 -I $H/strict7b.py $SNAP $CERT/census/results_v4 $H/out_strict7b --only 19_102 > $H/logs/strict7b.log 2>&1
echo done > $H/logs/v2_done
