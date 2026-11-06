#!/bin/bash
# usage: sbatch scripts/vesper_post_batch.sbatch   (reads $H/from_pm/READY.txt: one verified target per line)
#SBATCH -J mw_post_batch
#SBATCH -A swpark
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -t 12:00:00
#SBATCH -o %x_%j.out
unset OMP_NUM_THREADS
H=$HOME/mobywat_test
cd $H
# 2 MobyWat runs per target are sequential inside a task; 48 targets in parallel, each ~single-core
grep -v '^$' ${LIST:-from_pm/READY.txt} | while read t; do [ -f from_pm/$t/STATUS_POST ] && grep -q '^OK' from_pm/$t/STATUS_POST || echo $t; done | xargs -P ${PAR:-48} -I{} $H/scripts/vesper_post_task.sh {}
