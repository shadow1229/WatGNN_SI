#!/usr/bin/env python3
"""Remove water molecules (residue name SOL: OW, HW1, HW2) from MobyWat_aligned/system_tpy/<target>.pdb, replacing each file in place.
Protein atoms, ions (NA, CL) and all non-atom lines (REMARK, CRYST1, MODEL, TER, ENDMDL) are kept unchanged; atom serial numbers are not renumbered.
usage: remove_water_system_tpy.py [system_tpy dir]
The untouched originals remain in vesper_mobywat_test/from_pm/<target>/system_tpy.pdb."""
import os, sys

d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "MobyWat_aligned", "system_tpy")
tot_before = tot_after = 0
files = sorted(f for f in os.listdir(d) if f.endswith(".pdb"))
for f in files:
    p = os.path.join(d, f)
    keep, nb = [], 0
    for line in open(p):
        if line.startswith(("ATOM", "HETATM")):
            nb += 1
            if line[17:20].strip() == "SOL":
                continue
        keep.append(line)
    na = sum(1 for l in keep if l.startswith(("ATOM", "HETATM")))
    tmp = p + ".tmp"
    open(tmp, "w").writelines(keep)
    os.replace(tmp, p)
    tot_before += nb; tot_after += na
print(f"{len(files)} files: atoms {tot_before} -> {tot_after} (removed {tot_before - tot_after} water atoms)")
