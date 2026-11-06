#!/bin/bash
# usage: run_md.sh <PDB id_chain> <ref pdb> <workdir> <ns> <seed>
set -euo pipefail
ID=$1; REF=$2; W=$3; NS=$4; SEED=${5:-1}
GMX=${GMX:-gmx}
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$W" && cd "$W"
t0=$(date +%s)
python3 "$HERE/scripts/strip_water.py" "$REF" dry.pdb
$GMX pdb2gmx -f dry.pdb -water tip3p -ff amber99sb-ildn -ignh -o conf.gro -p topol.top -i posre.itp > pdb2gmx.log 2>&1
$GMX editconf -f conf.gro -o box.gro -d 1.5 -bt cubic > editconf.log 2>&1
$GMX solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top > solvate.log 2>&1
$GMX grompp -f "$HERE/mdp/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 > grompp_ions.log 2>&1
CHG=$(grep -E "non-zero total charge" grompp_ions.log | awk '{print $NF}')
CHG=${CHG:-0}
echo "net charge before ions: $CHG"
NPOS=$(python3 -c "c=float('$CHG'); print(max(0,int(round(-c))))")
NNEG=$(python3 -c "c=float('$CHG'); print(max(0,int(round(c))))")
echo SOL | $GMX genion -s ions.tpr -o ion.gro -p topol.top -pname NA -nname CL -np $NPOS -nn $NNEG > genion.log 2>&1
$GMX grompp -f "$HERE/mdp/steep.mdp" -c ion.gro -r ion.gro -p topol.top -o st.tpr -maxwarn 2 > grompp_st.log 2>&1
$GMX mdrun -deffnm st -v > mdrun_st.log 2>&1
$GMX grompp -f "$HERE/mdp/cg.mdp" -c st.gro -r st.gro -p topol.top -o cg.tpr -maxwarn 2 > grompp_cg.log 2>&1
$GMX mdrun -deffnm cg -v > mdrun_cg.log 2>&1
NSTEPS=$(python3 -c "print(int($NS*1000/0.002))")
sed -e "s/@NSTEPS@/$NSTEPS/" -e "s/@NSTXTC@/500/" -e "s/@SEED@/$SEED/" "$HERE/mdp/md.mdp.in" > md.mdp
$GMX grompp -f md.mdp -c cg.gro -r cg.gro -p topol.top -o md.tpr -maxwarn 2 > grompp_md.log 2>&1
t1=$(date +%s)
echo "prep+minimization wall time: $((t1-t0)) s"
printf "prep_min\t%s\n" $((t1-t0)) >> timing.tsv
if [ "${STOP_BEFORE_MD:-0}" = 1 ]; then exit 0; fi
$GMX mdrun -deffnm md -ntmpi 1 -ntomp ${NTOMP:-8} > mdrun_md.log 2>&1
t2=$(date +%s)
echo "MD wall time: $((t2-t1)) s"
printf "md_ntomp%s\t%s\n" ${NTOMP:-8} $((t2-t1)) >> timing.tsv
grep -E "Performance" md.log >> timing.tsv || true
