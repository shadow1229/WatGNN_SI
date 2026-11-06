#!/bin/bash
#SBATCH -J mw_post_test
#SBATCH -A swpark
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 00:20:00
#SBATCH -o %x_%j.out
export GMX=$HOME/.conda/envs/gromacs/bin/gmx
H=$HOME/mobywat_test
rm -rf $H/post_test && mkdir -p $H/post_test && cd $H/post_test
cp ../run_1byi/{md.tpr,conf.gro,bench_8.xtc} . && mv bench_8.xtc md.xtc
$H/scripts/post_mobywat.sh $H/post_test $H/bin/mobywat
cat timing.tsv; ls; head -12 O_system_prIDa.pdb 2>/dev/null
