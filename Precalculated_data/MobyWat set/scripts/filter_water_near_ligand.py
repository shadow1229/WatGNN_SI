#!/usr/bin/env python3
"""
Keep only predicted water positions within 3.5 / 4.5 A of ligand atoms,
sorted by the B-factor column (= prediction score) in descending order
(MobyWat: ascending, since the column holds mobility values).

Methods / folders (relative to this script):
  CNN_mobywat/<T>.pdb, GNN_mobywat/<T>_pred.pdb, hydraprot_mobywat/<T>_waters.pdb,
  rism_mobywat/<T>_pred.pdb, mobywat_mtol15/<T>/O_<T>_prIDa.pdb
Ligand: chains listed in trg_list_rec_lig.txt (target,protein,ligand; '+' joins chains)
        taken from original_PDB/<T>.pdb (same source as build_inputs.py).
        - residues in REAGENTS and ACE, and HOH, are ignored (same as build_inputs.py)
        - a predicted water must be within the cutoff of BOTH ligand and receptor heavy atoms
        - only altloc ' ' / 'A' are used
        - hydrogens are ignored
Output: water_cutoff_<cutoff>/<method>/<target>_pred.pdb  (HETATM, atom name O, resname HOH, chain X)
"""
import os, sys
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
CUTOFFS = [3.5,4.5]
PADDING = 1.0 #PADDING distance for 1.0A recall
ASCENDING = {"mobywat_mtol15"}
# only these residue names are used as predicted waters (others: all records)
WATER_RESNAME = {"GNN_mobywat": {"PRD"}}
REAGENTS = {"ACT", "CL", "EDO", "EOH", "GOL", "IUM", "U1", "PO4", "SO4", "TFA"}
# same rule as build_inputs.py: reagents + ACE are not used as receptor/ligand atoms
# in the interfacial-water distance test (ACE stays a normal residue elsewhere)
DIST_EXCLUDE = REAGENTS | {"ACE"}
METHODS = {  # name: (folder, filename pattern)
    "CNN_mobywat":       "{t}.pdb",
    "GNN_mobywat":       "{t}_pred.pdb",
    "hydraprot_mobywat": "{t}_waters.pdb",
    "rism_mobywat":      "{t}_pred.pdb",
    "mobywat_mtol15":    "{t}/O_{t}_prIDa.pdb",
}

def parse_atom(line):
    try:
        return dict(
            alt=line[16], res=line[17:20].strip(), chain=line[21],
            xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
            elem=line[76:78].strip() if len(line) >= 78 else "",
            name=line[12:16].strip(),
            b=float(line[60:66]),
        )
    except ValueError:
        return None

def is_h(a):
    e = a["elem"].upper()
    return e in ("H", "D") if e else a["name"].lstrip("0123456789")[:1].upper() in ("H", "D")

def read_targets():
    out = []
    with open(os.path.join(ROOT, "trg_list_rec_lig.txt")) as f:
        for ln in f:
            p = ln.strip().split(",")
            if len(p) >= 3:
                out.append((p[0], p[1].split("+"), p[2].split("+")))
    return out

def reclig_coords(target, chains):
    xyz = []
    with open(os.path.join(ROOT, "original_PDB", target + ".pdb")) as f:
        for ln in f:
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            a = parse_atom(ln)
            if a is None or a["chain"] not in chains:
                continue
            if a["alt"] not in (" ", "A"):
                continue
            if a["res"] in DIST_EXCLUDE or a["res"] == "HOH" or is_h(a):
                continue
            xyz.append(a["xyz"])
    return np.array(xyz)

def main():
    summary = []
    for target, rec, lig in read_targets():
        lxyz = reclig_coords(target, set(lig))
        rxyz = reclig_coords(target, set(rec))
        if len(lxyz) == 0:
            print(f"[WARN] {target}: no ligand atoms found for chains {lig}")
            continue
        if len(rxyz) == 0:
            print(f"[WARN] {target}: no receptor atoms found for chains {lig}")
            continue
        for method, pat in METHODS.items():
            path = os.path.join(ROOT, method, pat.format(t=target))
            if not os.path.exists(path):
                print(f"[WARN] missing {path}")
                continue
            lines, coords = [], []
            with open(path, newline="") as f:
                for ln in f:
                    ln = ln.rstrip("\r\n")
                    if ln.startswith(("ATOM", "HETATM")):
                        a = parse_atom(ln)
                        if a is None:
                            continue
                        if method in WATER_RESNAME and a["res"] not in WATER_RESNAME[method]:
                            continue
                        lines.append((a["b"], ln))
                        coords.append(a["xyz"])
            if not lines:
                print(f"[WARN] no atoms in {path}"); continue
            c = np.array(coords)
            dmin_l = np.sqrt(((c[:, None, :] - lxyz[None, :, :]) ** 2).sum(-1)).min(axis=1)
            dmin_r = np.sqrt(((c[:, None, :] - rxyz[None, :, :]) ** 2).sum(-1)).min(axis=1)
            for cut in CUTOFFS:
                idx = [i for i in range(len(lines)) if (dmin_l[i] <= (cut+PADDING) and dmin_r[i]<=(cut+PADDING))]
                # B-factor column = prediction score: descending, except MobyWat
                # (mobility value) which is ascending. Stable: ties keep original order.
                sign = 1 if method in ASCENDING else -1
                idx.sort(key=lambda i: sign * lines[i][0])
                od = os.path.join(ROOT, f"water_cutoff_{cut}", method)
                os.makedirs(od, exist_ok=True)
                with open(os.path.join(od, f"{target}_pred.pdb"), "w", newline="\n") as o:
                    for i in idx:
                        ln = lines[i][1].ljust(54)
                        # fix residue name -> HOH (cols 18-20), chain -> X (col 22)
                        ln = "HETATM" + ln[6:12] + " O  " + ln[16] + "HOH" + ln[20] + "X" + ln[22:]
                        # element column (77-78), if present, -> O
                        if len(ln.rstrip()) >= 78:
                            ln = ln[:76] + " O" + ln[78:]
                        o.write(ln.rstrip() + "\n")
                    o.write("END\n")
                summary.append((target, method, cut, len(lines), len(idx)))
    with open(os.path.join(ROOT, "water_cutoff_summary.tsv"), "w") as o:
        o.write("target\tmethod\tcutoff\tn_predicted\tn_kept\n")
        for s in summary:
            o.write("\t".join(map(str, s)) + "\n")
    print(f"done: {len(summary)} outputs")

if __name__ == "__main__":
    main()
