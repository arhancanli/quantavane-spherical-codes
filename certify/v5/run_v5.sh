#!/bin/zsh
# Every computation of version 5 of the paper, end to end (light background jobs: nice 10, two BLAS threads, at most two
# processes at a time). SNAP is the folder with the downloaded data set; PAPER the folder of the paper's LaTeX sources.
#   zsh certify/v5/run_v5.sh SNAP PAPER
set -e
SNAP=${1:?snapshot folder}
PAPER=${2:?paper folder}
V=${0:A:h}
cd $V
export OPENBLAS_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 OMP_NUM_THREADS=2
PY=(nice -n 10 /usr/bin/python3 -I)
mkdir -p out out/improved
stamp() { echo "$(date -u +%FT%TZ) $1" >> out/run_v5.log }
stamp "start"
# flat cores (Section 5.2)
for e in "16 117 3" "17 162 3" "20 92 3" "22 79 1" "22 114 1" "22 127 1" "23 102 1" "23 129 1" "23 147 1" "24 121 1" \
         "24 149 1" "26 57 1" "27 96 6"; do
  parts=(${=e})
  $PY subspace.py $SNAP $parts[1] $parts[2] out/subspace_d$parts[1]_n$parts[2].json --disc $parts[3] \
      > out/subspace_$parts[1]_$parts[2].log 2>&1
done
stamp "subspace done"
# split contact graphs and strict cores (Sections 5.3 and 5.4); two at a time
$PY rigidity.py $SNAP 8 72 out/rigidity_d8_n72.json > out/rigidity_8_72.log 2>&1 &
$PY rigidity.py $SNAP 12 210 out/rigidity_d12_n210.json > out/rigidity_12_210.log 2>&1 &
wait
$PY rigidity.py $SNAP 10 238 out/rigidity_d10_n238.json --disc 3 > out/rigidity_10_238.log 2>&1 &
$PY rigidity.py $SNAP 19 166 out/rigidity_d19_n166.json > out/rigidity_19_166.log 2>&1 &
wait
$PY rigidity.py $SNAP 9 330 out/rigidity_d9_n330.json > out/rigidity_9_330.log 2>&1
stamp "rigidity done"
# (22,116), the jammed code that is not infinitesimally jammed: now the separate note (certify/note/), no longer in
# the census paper
$PY recognise116.py $SNAP out/recognition_d22_n116.json > out/recognise116.log 2>&1
$PY jam116.py $SNAP out/recognition_d22_n116.json out/jam_d22_n116.json > out/jam116.log 2>&1
stamp "jam116 done"
# certified distances of the exact-only confirmations of version 4 (Section 4.4)
$PY distances.py $SNAP out/distances_exact_only.json 3,12,5 4,120,5 6,27,3 6,72,3 4,10,1 14,42,1 14,64,1 14,70,1 \
    20,56,1 20,81,1 20,105,1 20,120,1 21,50,1 21,77,1 21,112,1 22,100,1 21,162,1 22,275,1 > out/distances_exact_only.log 2>&1
stamp "distances done"
# numerical evidence on the continuous families (Section 5.5; untrusted)
$PY families.py $SNAP out/families.json 22,47 22,48 23,50 > out/families.log 2>&1
$PY rigidity.py $SNAP 22 48 out/rigidity_d22_n48.json --disc 145 > out/rigidity_22_48.log 2>&1 || true
stamp "families done"
# the local searches from the NONOPT certificates (version 3 search), evaluated exactly
for n in 576 791; do
  $PY extract_code.py $SNAP 3 $n out/improved/stored_d3_n$n.txt
  $PY maxcos_exact.py out/improved/d3_n$n.txt 3 > out/improved/maxcos_d3_n$n.json
  $PY maxcos_exact.py out/improved/stored_d3_n$n.txt 3 > out/improved/maxcos_stored_d3_n$n.json
done
stamp "maxcos done"
$PY summarize_v5.py $PAPER > out/summarize_v5.log 2>&1
$PY census_md_v5.py > $PAPER/CENSUS.md
$PY release_list.py >> out/summarize_v5.log 2>&1
stamp "done"
