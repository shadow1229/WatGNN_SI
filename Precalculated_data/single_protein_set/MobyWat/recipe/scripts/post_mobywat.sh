#!/bin/bash
# usage: post_mobywat.sh <workdir> <mobywat binary>
# Post-MD processing (MobyWat manual Sec. 9.1.2.1) + MobyWat prediction (IDa/IDe/POS/MER)
set -euo pipefail
W=$1; MW=$2; GMX=${GMX:-gmx}
cd "$W"
t0=$(date +%s)
echo 0   | $GMX trjconv -f md.xtc    -s md.tpr -o pbc_1.xtc -pbc whole > trj1.log 2>&1
printf "1\n0\n" | $GMX trjconv -f pbc_1.xtc -s md.tpr -o pbc_2.xtc -pbc cluster > trj2.log 2>&1
printf "1\n0\n" | $GMX trjconv -f pbc_2.xtc -s md.tpr -o pbc_3.xtc -center -pbc mol -ur compact > trj3.log 2>&1
printf "3\n0\n" | $GMX trjconv -f pbc_3.xtc -s md.tpr -o system.xtc -fit progressive > trj4.log 2>&1
printf "3\n0\n" | $GMX trjconv -f pbc_3.xtc -s md.tpr -o system_tpy.pdb -b 0 -e 0 -fit progressive > trj5.log 2>&1
NFR=$($GMX check -f system.xtc 2>&1 | awk '/^Step/ {print $2}')
NPROT=$(awk 'NR==2{print $1}' conf.gro)
LAST=$((NFR-1))
echo "frames=$NFR protein atoms=$NPROT"
t1=$(date +%s)
"$MW" -f system.xtc -tpy system_tpy.pdb -t 1-$NPROT -w Auto -n 0-$LAST -cls MER -v Verbose > mobywat.out 2>&1
t2=$(date +%s)
echo "post-MD wall time: $((t1-t0)) s; MobyWat wall time: $((t2-t1)) s"
printf "postprocess\t%s\nmobywat\t%s\n" $((t1-t0)) $((t2-t1)) >> timing.tsv
