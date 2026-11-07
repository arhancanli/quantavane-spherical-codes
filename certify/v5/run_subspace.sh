#!/bin/zsh
# exact flat-core certificates (subspace.py) for the census entries whose contact core is flat
S=${1:?snapdir}
cd ${0:A:h}
run() { perl -e 'alarm shift; exec @ARGV' 3600 nice -n 10 python3 -I subspace.py $S $1 $2 out/subspace_d$1_n$2.json --disc $3 > out/subspace_$1_$2.log 2>&1 }
run 22 79 1; run 22 114 1; run 23 102 1; run 23 129 1; run 24 149 1; run 26 57 1; run 23 147 1; run 24 121 1
run 16 117 3; run 17 162 3; run 20 92 3; run 27 96 6
echo done
