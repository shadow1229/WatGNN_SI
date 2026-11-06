#!/bin/bash
#SBATCH -J mw_prep_bench
#SBATCH -A swpark
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 01:00:00
#SBATCH -o %x_%j.out
set -eo pipefail
export GMX=$HOME/.conda/envs/gromacs/bin/gmx
export OMP_NUM_THREADS=8
H=$HOME/mobywat_test
cd $H
t0=$(date +%s)
STOP_BEFORE_MD=1 ./scripts/run_md.sh 1byi_A $H/ref/1byi_A.pdb $H/run_1byi 10 1
echo "prep wall: $(( $(date +%s)-t0 )) s"
cd $H/run_1byi
# 50 ps benchmarks, three thread counts at once on disjoint cores
for spec in "8 0" "16 8" "32 24"; do
  set -- $spec
  ( $GMX mdrun -s md.tpr -deffnm bench_$1 -nsteps 25000 -ntmpi 1 -ntomp $1 -pin on -pinoffset $2 -pinstride 1 -noconfout > bench_$1.out 2>&1 ) &
done
wait
grep -H "Performance" bench_*.log
