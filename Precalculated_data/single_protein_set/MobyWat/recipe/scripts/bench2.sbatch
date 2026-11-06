#!/bin/bash
#SBATCH -J mw_bench2
#SBATCH -A swpark
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 112
#SBATCH -t 00:30:00
#SBATCH -o %x_%j.out
export GMX=$HOME/.conda/envs/gromacs/bin/gmx
unset OMP_NUM_THREADS
cd $HOME/mobywat_test/run_1byi
t0=$(date +%s)
for spec in "16 0" "32 16" "64 48"; do
  set -- $spec
  ( s=$(date +%s); $GMX mdrun -s md.tpr -deffnm bench_$1 -nsteps 25000 -ntmpi 1 -ntomp $1 -pin on -pinoffset $2 -pinstride 1 -noconfout > bench_$1.out 2>&1; echo "ntomp=$1 wall=$(( $(date +%s)-s )) s" ) &
done
wait
echo "total wall: $(( $(date +%s)-t0 )) s"
grep -H "Performance" bench_*.log
