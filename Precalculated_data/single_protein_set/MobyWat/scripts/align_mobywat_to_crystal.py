#!/usr/bin/env python3
"""MobyWat predictions are in the coordinate frame of the processed MD trajectory (system_tpy.pdb: protein centred in the simulation box).
Superimpose that frame onto the crystal structure (dry.pdb) with a Calpha Kabsch fit and apply the same rigid transform to the predicted waters.
usage: align_mobywat_to_crystal.py <from_pm dir> <targets.txt> <MobyWat dir (box frame)> <output dir>
output: <output dir>/<target>.pdb (crystal frame) and alignment_transforms.csv (R, t, Calpha RMSD before/after)"""
import csv, os, sys
import numpy as np

src, tfile, indir, outdir = sys.argv[1:5]
os.makedirs(outdir, exist_ok=True)

def ca(path, first_model=False):
    xyz = []
    for l in open(path):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        elif first_model and l.startswith("ENDMDL"):
            break
    return np.array(xyz)

def kabsch(P, Q):   # R, t minimising |R p + t - q|
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    R = Vt.T @ D @ U.T
    return R, qc - R @ pc

rows = []
for t in [l.strip() for l in open(tfile) if l.strip()]:
    d = os.path.join(src, t)
    P, Q = ca(os.path.join(d, "dry.pdb")), ca(os.path.join(d, "system_tpy.pdb"), True)   # crystal, MD frame 0
    assert len(P) == len(Q), t
    R, tv = kabsch(Q, P)
    rms = lambda A, B: float(np.sqrt(((A - B) ** 2).sum(1).mean()))
    with open(os.path.join(indir, f"{t}.pdb")) as fi, open(os.path.join(outdir, f"{t}.pdb"), "w") as fo:
        for l in fi:
            x = R @ np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) + tv
            fo.write(l[:30] + "%8.3f%8.3f%8.3f" % tuple(x) + l[54:])
    rows.append([t, round(rms(Q, P), 3), round(rms(Q @ R.T + tv, P), 3)] + [round(v, 6) for v in R.ravel()] + [round(v, 3) for v in tv])
with open(os.path.join(outdir, "alignment_transforms.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["target", "CA_rmsd_before_A", "CA_rmsd_after_A"] + [f"R{i}{j}" for i in range(3) for j in range(3)] + ["tx", "ty", "tz"])
    w.writerows(rows)
print(len(rows), "targets ->", outdir)
