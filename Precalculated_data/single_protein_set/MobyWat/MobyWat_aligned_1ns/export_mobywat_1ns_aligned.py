#!/usr/bin/env python3
"""First-1-ns MobyWat predictions (mw_n1000/O_system_prIDa.pdb, frames 0-1000) in the same format and crystal frame as MobyWat_aligned/<target>.pdb.
Used for a limited sampling-length check (1 ns vs 10 ns).
usage: export_mobywat_1ns_aligned.py <from_pm dir> <targets.txt> <output dir>"""
import os, sys
import numpy as np

src, tfile, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)

def ca(p, first_model=False):
    x = []
    for l in open(p):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            x.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        elif first_model and l.startswith("ENDMDL"):
            break
    return np.array(x)

def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc

n = sites = 0
for t in [l.strip() for l in open(tfile) if l.strip()]:
    d = os.path.join(src, t)
    R, tv = kabsch(ca(os.path.join(d, "system_tpy.pdb"), True), ca(os.path.join(d, "dry.pdb")))
    with open(os.path.join(d, "mw_n1000", "O_system_prIDa.pdb")) as fi, open(os.path.join(out, f"{t}.pdb"), "w") as fo:
        for l in fi:
            if l.startswith(("HETATM", "ATOM")):
                x = R @ np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) + tv
                fo.write(l[:12] + " O  " + l[16] + "HOH" + l[20:30] + "%8.3f%8.3f%8.3f" % tuple(x) + l[54:])
                sites += 1
    n += 1
print(f"{n} targets, {sites} sites -> {out}")
