#!/bin/bash
# vesper (CPU) stage, one target: gmx trjconv chain (MobyWat manual 9.1.2.1) + MobyWat on all frames (10 ns) and first 1 ns.
# input dir: $B/<target>/ with md.xtc, md.tpr, conf.gro (copied from Perlmutter, sha256-verified)
T=$1
H=${H:-$HOME/mobywat_test}
W=$H/from_pm/$T
GMX=${GMX:-$HOME/.conda/envs/gromacs/bin/gmx}
MW=$H/bin/mobywat
export GMX
cd "$W" || exit 1
fail() { echo "FAIL $1 $(date +%FT%T)" > STATUS_POST; exit 1; }
t0=$(date +%s)
printf "0\n"    | $GMX trjconv -f md.xtc    -s md.tpr -o pbc_1.xtc -pbc whole                    > trj1.log 2>&1 || fail trjconv1
printf "1\n0\n" | $GMX trjconv -f pbc_1.xtc -s md.tpr -o pbc_2.xtc -pbc cluster                  > trj2.log 2>&1 || fail trjconv2
rm -f pbc_1.xtc
printf "1\n0\n" | $GMX trjconv -f pbc_2.xtc -s md.tpr -o pbc_3.xtc -center -pbc mol -ur compact   > trj3.log 2>&1 || fail trjconv3
rm -f pbc_2.xtc
printf "3\n0\n" | $GMX trjconv -f pbc_3.xtc -s md.tpr -o system.xtc -fit progressive              > trj4.log 2>&1 || fail trjconv4
printf "3\n0\n" | $GMX trjconv -f pbc_3.xtc -s md.tpr -o system_tpy.pdb -b 0 -e 0 -fit progressive > trj5.log 2>&1 || fail trjconv5
rm -f pbc_3.xtc
printf "postprocess\t%s\n" $(( $(date +%s)-t0 )) >> timing_post.tsv
NFR=$($GMX check -f system.xtc 2>&1 | awk '/^Step/ {print $2}')
NPROT=$(awk 'NR==2{print $1}' conf.gro)
echo "frames=$NFR protein_atoms=$NPROT" > frames.txt
NTOT=$(sed -n 2p md.gro | tr -d ' ')
if [ "$NTOT" -gt 99999 ]; then
  # PDB serials wrap at 99999 -> MobyWat atom ranges become ambiguous: keep protein + ions + nearest waters (<=99999 atoms)
  s=$(date +%s)
  python3 $H/scripts/crop_waters.py system_tpy.pdb crop.ndx crop_info.txt 99999 > crop.log 2>&1 || fail crop_select
  echo Crop | $GMX trjconv -f system.xtc -s md.tpr -n crop.ndx -o system_crop.xtc > trj_crop1.log 2>&1 || fail crop_xtc
  echo Crop | $GMX trjconv -f system.xtc -s md.tpr -n crop.ndx -o system_tpy_crop.pdb -b 0 -e 0 > trj_crop2.log 2>&1 || fail crop_pdb
  mv system_tpy.pdb system_tpy_full.pdb; mv system_tpy_crop.pdb system_tpy.pdb
  mv system_crop.xtc system.xtc
  printf "crop_waters\t%s\n" $(( $(date +%s)-s )) >> timing_post.tsv
  echo "frames=$NFR protein_atoms=$NPROT cropped_from=$NTOT" > frames.txt
fi
for N in 10000 1000; do
  LAST=$N; [ $((NFR-1)) -lt $N ] && LAST=$((NFR-1))
  mkdir -p mw_n$N && cd mw_n$N
  ln -sf ../system.xtc system.xtc; ln -sf ../system_tpy.pdb system_tpy.pdb
  s=$(date +%s)
  "$MW" -f system.xtc -tpy system_tpy.pdb -t 1-$NPROT -w Auto -n 0-$LAST -cls MER -v Verbose > mobywat.out 2>&1 || { cd ..; fail mobywat_$N; }
  printf "mobywat_n%s\t%s\n" $N $(( $(date +%s)-s )) >> ../timing_post.tsv
  rm -f system.xtc system_tpy.pdb O_system.xtc O_system.plw
  cd ..
done
rm -f system.xtc
printf "total_post\t%s\n" $(( $(date +%s)-t0 )) >> timing_post.tsv
echo "OK $(date +%FT%T)" > STATUS_POST
