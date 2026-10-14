#!/usr/bin/env python3
"""Build water-prediction inputs from original_PDB following the MobyWat protocol.

protein/<ID>.pdb : receptor + ligand (+ reagents within 5 A of ligand), no HOH,
                   altloc 'A' (or blank) only, no fallback; hydrogens present in the input are KEPT
water/<ID>.pdb    : HOH of all chains, no B-factor limit, same altloc rule
water_45/<ID>.pdb : as water/, and within 4.5 A of BOTH a ligand and a receptor heavy atom
water_35/<ID>.pdb : as water/, and within 3.5 A of BOTH a ligand and a receptor heavy atom

Rules: reagents = the crystallization-reagent list below (ignored as receptor/ligand atoms, kept in
protein/ only if within CUTOFF of the ligand); no B-factor upper limit; hydrogens are excluded
only when choosing interfacial waters (distance test), not from protein/<ID>.pdb.
"""
import os
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
REAGENTS = {"ACT", "CL", "EDO", "EOH", "GOL", "IUM", "U1", "PO4", "SO4", "TFA"}
# residues never used as ligand/receptor atoms for the interfacial-water test.
# ACE is excluded from the test only; it stays in protein/ as a normal residue.
DIST_EXCLUDE = REAGENTS | {"ACE"}
CUTOFF = 5.0
BFAC_MAX = 40


def read_targets():
    out = []
    with open(os.path.join(ROOT, "trg_list_rec_lig.txt")) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            pid, rec, lig = line.split(",")
            out.append((pid, rec.split("+"), lig.split("+")))
    return out


def filter_altloc(lines):
    """Keep only blank or 'A' altloc atoms (no fallback to other altlocs)."""
    kept = [l for l in lines if l[16] in (" ", "A")]
    return [l[:16] + " " + l[17:] for l in kept], 0


def is_h(l):
    e = l[76:78].strip().upper()
    if e:
        return e in ("H", "D")
    return l[12:16].strip().lstrip("0123456789")[:1].upper() in ("H", "D")


def residue_id(l):
    return (l[21], l[22:27], l[17:20])


def main():
    os.makedirs(os.path.join(ROOT, "protein"), exist_ok=True)
    for d in ("water", "water_45", "water_35"):
        os.makedirs(os.path.join(ROOT, d), exist_ok=True)
    for pid, rec, lig in read_targets():
        with open(os.path.join(ROOT, "original_PDB", pid + ".pdb")) as f:
            raw = [l.rstrip("\n").ljust(80) for l in f
                   if l.startswith(("ATOM  ", "HETATM"))]
        atoms, nfb = filter_altloc(raw)
        keep_chains = set(rec) | set(lig)
        waters = [l for l in atoms if l[17:20] == "HOH"]
        atoms = [l for l in atoms if l[17:20] != "HOH"]

        # ligand atoms: non-reagent atoms in ligand chains
        lig_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                            for l in atoms
                            if l[21] in lig and l[17:20].strip() not in DIST_EXCLUDE and not is_h(l)])
        
        rec_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                            for l in atoms
                            if l[21] in rec and l[17:20].strip() not in DIST_EXCLUDE and not is_h(l)])
        def xyz_of(l):
            return [float(l[30:38]), float(l[38:46]), float(l[46:54])]

        water_b = [l for l in waters if BFAC_MAX is None or float(l[60:66]) < BFAC_MAX]
        if water_b:
            wd_l = np.linalg.norm(np.array([xyz_of(l) for l in water_b])[:, None]
                                - lig_xyz[None], axis=2).min(axis=1)
            wd_r = np.linalg.norm(np.array([xyz_of(l) for l in water_b])[:, None]
                                - rec_xyz[None], axis=2).min(axis=1)
        else:
            wd_l = np.array([])
            wd_r = np.array([])
        water_45 = [l for l, dl, dr in zip(water_b, wd_l, wd_r) if (dl <= 4.5 and dr <= 4.5)]
        water_35 = [l for l, dl, dr in zip(water_b, wd_l, wd_r) if (dl <= 3.5 and dr <= 3.5)]

        # reagent residues (any chain) within CUTOFF of ligand
        res_atoms = {}
        for l in atoms:
            if l[17:20].strip() in REAGENTS:
                res_atoms.setdefault(residue_id(l), []).append(l)
        near = set()
        for rid, ls in res_atoms.items():
            xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in ls])
            d = np.linalg.norm(xyz[:, None] - lig_xyz[None], axis=2).min()
            if d <= CUTOFF:
                near.add(rid)

        out, dropped = [], {}
        for l in atoms:
            name = l[17:20].strip()
            if name in REAGENTS:
                if residue_id(l) in near:
                    out.append(l)
                else:
                    dropped[name] = dropped.get(name, set()) | {residue_id(l)}
            elif l[21] in keep_chains:
                out.append(l)

        def write(path, ls):
            with open(path, "w") as f:
                prev = None
                for l in ls:
                    if prev is not None and l[21] != prev:
                        f.write("TER\n")
                    f.write(l.rstrip() + "\n")
                    prev = l[21]
                f.write("TER\nEND\n")

        write(os.path.join(ROOT, "protein", pid + ".pdb"), out)
        write(os.path.join(ROOT, "water", pid + ".pdb"), water_b)
        write(os.path.join(ROOT, "water_45", pid + ".pdb"), water_45)
        write(os.path.join(ROOT, "water_35", pid + ".pdb"), water_35)
        nr = sum(1 for rid in near)
        print(f"{pid}: protein atoms={len(out)} reagent residues kept={nr} "
              f"dropped={ {k: len(v) for k, v in dropped.items()} } "
              f"water={len(water_b)}/{len(waters)} 4.5A={len(water_45)} 3.5A={len(water_35)} altloc-fallback={nfb}")


if __name__ == "__main__":
    main()
