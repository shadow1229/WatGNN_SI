#!/bin/bash
#SBATCH -J mw_1byi_10ns
#SBATCH -A swpark
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 08:00:00
#SBATCH -o %x_%j.out
export GMX=$HOME/.conda/envs/gromacs/bin/gmx
unset OMP_NUM_THREADS
H=$HOME/mobywat_test
cd $H/run_1byi
printf "prep_min\t846\n" >> timing.tsv   # job 2959
s=$(date +%s)
$GMX mdrun -s md.tpr -deffnm md -ntmpi 1 -ntomp 64 -pin on -pinoffset 0 -pinstride 1 > mdrun_md.log 2>&1
e=$(date +%s)
echo "MD 10 ns wall: $((e-s)) s"
printf "md_10ns_ntomp64\t%s\n" $((e-s)) >> timing.tsv
grep Performance md.log >> timing.tsv
$H/scripts/post_mobywat.sh $H/run_1byi $H/bin/mobywat
cat timing.tsv
