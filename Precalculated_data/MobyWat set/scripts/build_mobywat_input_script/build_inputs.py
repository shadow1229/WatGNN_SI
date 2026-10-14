#!/usr/bin/env python3
"""Build water-prediction inputs from original_PDB following the MobyWat protocol.

protein/<ID>.pdb : receptor + ligand (+ reagents within 5 A of ligand), no HOH,
                   altloc 'A' (or blank) only
water/<ID>.pdb    : HOH of all chains with B-factor < 40, same altloc rule
water_45/<ID>.pdb : as water/, and within 4.5 A of a ligand atom
water_35/<ID>.pdb : as water/, and within 3.5 A of a ligand atom
"""
import os
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
REAGENTS = {"ACT", "CL", "EDO", "EOH", "GOL", "IUM", "U1", "PO4", "SO4", "TFA"}
CUTOFF = 5.0
BFAC_MAX = 40.0


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
    """Keep blank/'A' altlocs; if an atom has altlocs but no 'A', keep its first one."""
    def key(l):
        return (l[21], l[22:27], l[17:20], l[12:16])
    has_a, first = set(), {}
    for l in lines:
        k = key(l)
        if l[16] == "A":
            has_a.add(k)
        if l[16] != " " and k not in first:
            first[k] = l[16]
    kept, fallback = [], 0
    for l in lines:
        a = l[16]
        if a == " " or a == "A":
            kept.append(l)
        elif key(l) not in has_a and first[key(l)] == a:
            kept.append(l)
            fallback += 1
    return [l[:16] + " " + l[17:] for l in kept], fallback


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
                            if l[21] in lig and l[17:20].strip() not in REAGENTS])
        def xyz_of(l):
            return [float(l[30:38]), float(l[38:46]), float(l[46:54])]

        water_b = [l for l in waters if float(l[60:66]) < BFAC_MAX]
        if water_b:
            wd = np.linalg.norm(np.array([xyz_of(l) for l in water_b])[:, None]
                                - lig_xyz[None], axis=2).min(axis=1)
        else:
            wd = np.array([])
        water_45 = [l for l, d in zip(water_b, wd) if d <= 4.5]
        water_35 = [l for l, d in zip(water_b, wd) if d <= 3.5]

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
              f"water B<40={len(water_b)}/{len(waters)} 4.5A={len(water_45)} 3.5A={len(water_35)} altloc-fallback={nfb}")


if __name__ == "__main__":
    main()
