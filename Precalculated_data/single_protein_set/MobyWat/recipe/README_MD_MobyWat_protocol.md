# MobyWat benchmark on the WatGNN 91-target single-protein set — MD / preparation / post-processing record

## 1. Targets
91 single-protein targets from `WatGNN_SI/Dataset/single protein set/single_protein_set.txt`; 1byi_A (test) and 1mn8_A (pilot) were run first, then the other 89 as a Slurm array (concurrency 6). Input = reference crystal structure with all HOH removed (`strip_water.py`).

## 2. MD protocol (follows MobyWat Bioinformatics 2015 Suppl. S1.2 / JCIM 2016)
| item | setting |
|---|---|
| input | reference PDB, HOH and non-protein records removed, first altloc (` `/`A`) kept; PRO `OE` atom removed (1i1w: pyroglutamate-like N-terminus not in AMBER template) |
| software | GROMACS 2026.3 (conda-forge), mixed precision; CUDA build on Perlmutter, CPU build on vesper |
| force field / water | AMBER99SB-ILDN, TIP3P, `pdb2gmx -ignh` |
| box | cubic, 1.5 nm solute–box distance (paper: 15 Å; manual example: dodecahedron 10 Å) |
| neutralization | genion NA/CL, net charge only (no extra salt); neutral systems add no ions |
| minimization | steepest descent (emtol 1000, 50000 steps) then conjugate gradient (emtol 10, 50000 steps), both with `-DPOSRES`; CG with `-DFLEXIBLE` (rigid water fails in CG) |
| MD | NPT, 10 ns, dt 2 fs, h-bond constraints (LINCS), 300 K v-rescale τ=0.1 ps (Protein / Non-Protein), Parrinello–Rahman τ=0.5 ps, compressibility 4.5e-5, 1 bar, PME, cut-offs 1.1 nm, `refcoord-scaling=com` |
| restraints | `-DPOSRES`, 1000 kJ/mol/nm² on all protein heavy atoms, whole run; no separate equilibration (as in the paper) |
| velocities | gen-vel yes, 300 K, seed = replicate number (replicate 1: seed 1) |
| output | compressed xtc every 500 steps = 1 ps → 10,001 frames |
| GPU flags (Perlmutter) | A100 `-q shared -C gpu -G 1 -c 32`, `-nb gpu -pme gpu -bonded gpu -pin on` (update on CPU because of Parrinello–Rahman) |

Deviations from the published protocol (to be disclosed): mixed instead of double precision (CG does not reach Fmax<10 in mixed precision; steep/CG are only for preparation), modern GROMACS version, `-DFLEXIBLE` in CG, cubic 1.5 nm box, one 10 ns replicate so far (the paper uses repeats with different velocity seeds).

## 3. Post-processing and MobyWat (vesper CPU; `scripts/vesper_post_task.sh`)
GROMACS trjconv chain from MobyWat manual §9.1.2.1:
1. `-pbc whole` (System) → 2. `-pbc cluster` (Protein) → 3. `-center -pbc mol -ur compact` (Protein, System) → 4. `-fit progressive` on Cα (group 3) → `system.xtc`; topology frame `system_tpy.pdb` via `-b 0 -e 0`.

MobyWat: `mobywat -f system.xtc -tpy system_tpy.pdb -t 1-<N protein atoms> -w Auto -n 0-<last frame> -cls MER`, run twice per target: all frames (`mw_n10000`) and the first 1,000 frames (`mw_n1000`, ≈1 ns). Output: IDa / IDe / POS / MER hydration-site lists ordered by mobility.
Intermediates (pbc_*.xtc) are deleted after each step; `system.xtc` is deleted after MobyWat.

## 4. Systems with >= 100,000 atoms (water cropping before MobyWat)
MobyWat selects atoms by PDB serial number (`-t 1-N`) and the PDB serial field wraps after 99,999, so for systems with &gt;= 100,000 atoms the serial range 1..N also matches wrapped water atoms (warning "Found atom number in pdb: 2N" in `mobywat.out`; protein selection wrong, run extremely slow). The first attempts for 2fvy\_A, 2jhf\_A, 1muw\_A, 2bf6\_A were therefore invalid and were discarded.  
Affected targets (8): 4ga2\_A, 2fvy\_A, 2jhf\_A, 3vla\_A, 2bf6\_A, 3l8w\_A, 1muw\_A, 3qr7\_A (1jfb\_A, 99,065 atoms, is unaffected).  
Fix (`scripts/crop_waters.py`, applied automatically by `vesper_post_task.sh` if md.gro has &gt; 99,999 atoms): after the trjconv chain, keep protein + ions + the waters whose O is closest to the protein (min-image distance in frame 0) so that &lt;= 99,999 atoms remain; same atoms selected from system.xtc and system\_tpy.pdb; protein atoms stay first. Details per target: `crop_info.txt`, `crop.ndx`, full topology PDB kept as `system_tpy_full.pdb`.
Resulting shell (all waters within this distance of the protein kept): 1muw\_A 33.9 A, 2bf6\_A 35.3 A, 2jhf\_A 40.2 A, 3l8w\_A 33.7 A, 3qr7\_A 39.8 A, 3vla\_A 38.4 A, 2fvy\_A 45.1 A, 4ga2\_A 54.2 A, i.e. far beyond the hydration layers; only bulk water is removed. MD itself is unchanged (full boxes).
