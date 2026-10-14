#!/usr/bin/env python3
"""
  1) ref/[ID]_ligand.mol2        -> ref_ligand_pdb/[ID]_ligand.pdb
  2) ref/[ID]_protein.pdb + ligand PDB from 1) -> ref_complex_pdb/[ID]_complex.pdb

Usage:  python3 make_complex_pdb.py
options: --list, --ref, --resname, --chain
"""
import argparse, os, sys

ELEM2 = {"CL", "BR", "FE", "ZN", "MG", "MN", "CU", "CO", "NI", "NA", "CA", "SE", "LI", "AL", "SI", "HG", "CD"}

def parse_mol2(path):
    atoms, bonds, sec = [], [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("@<TRIPOS>"):
                sec = s[9:].upper()
                continue
            if not s or s.startswith("#"):
                continue
            t = s.split()
            if sec == "ATOM" and len(t) >= 6:
                name, x, y, z, atype = t[1], float(t[2]), float(t[3]), float(t[4]), t[5]
                elem = atype.split(".")[0].upper()
                if len(elem) > 2 or elem not in ELEM2 and len(elem) == 2:
                    elem = name.rstrip("0123456789")[:1].upper()
                atoms.append((int(t[0]), name, x, y, z, elem))
            elif sec == "BOND" and len(t) >= 3:
                bonds.append((int(t[1]), int(t[2]), t[3] if len(t) > 3 else "1"))
    return atoms, bonds

def fmt_name(name, elem):
    name = name[:4]
    return name.ljust(4) if (len(elem) == 2 or len(name) == 4) else (" " + name).ljust(4)

def ligand_lines(atoms, bonds, resname, chain, start_serial=1, resseq=1):
    out, idmap = [], {}
    for i, (aid, name, x, y, z, elem) in enumerate(atoms):
        serial = start_serial + i
        idmap[aid] = serial
        out.append("HETATM%5d %s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
                   % (serial, fmt_name(name, elem), resname[:3], chain, resseq, x, y, z, 1.0, 0.0, elem.capitalize().rjust(2)))
    nbr = {}
    for a, b, _ in bonds:
        nbr.setdefault(idmap[a], []).append(idmap[b])
        nbr.setdefault(idmap[b], []).append(idmap[a])
    conect = []
    for s in sorted(nbr):
        l = sorted(nbr[s])
        for k in range(0, len(l), 4):
            conect.append("CONECT%5d" % s + "".join("%5d" % n for n in l[k:k+4]))
    return out, conect

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", default="pdbbind_exclude30.txt")
    ap.add_argument("--ref", default="ref")
    ap.add_argument("--resname", default="LIG")
    ap.add_argument("--chain", default="L", help="ligand chain ID")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    os.chdir(here)
    os.makedirs("ref_ligand_pdb", exist_ok=True)
    os.makedirs("ref_complex_pdb", exist_ok=True)
    ids = [l.split()[0] for l in open(a.list) if l.strip() and not l.startswith("#")]
    ok, fail = 0, []
    for t in ids:
        mol2 = os.path.join(a.ref, f"{t}_ligand.mol2")
        prot = os.path.join(a.ref, f"{t}_protein.pdb")
        try:
            atoms, bonds = parse_mol2(mol2)
            if not atoms:
                raise ValueError("no atom in mol2 file")
            # 1. ligand pdb
            lig, con = ligand_lines(atoms, bonds, a.resname, a.chain)
            with open(f"ref_ligand_pdb/{t}_ligand.pdb", "w") as f:
                f.write("\n".join(lig + con + ["END"]) + "\n")
            # 2. complex pdb
            plines = [l.rstrip("\n") for l in open(prot)]
            body = [l for l in plines if l.startswith(("ATOM", "HETATM", "TER"))]
            serials = [int(l[6:11]) for l in body if l.startswith(("ATOM", "HETATM"))]
            chains = {l[21] for l in body if l.startswith(("ATOM", "HETATM"))}
            chain = a.chain
            if chain in chains:
                chain = next(c for c in "LXYZMNOPQRSTUVWKJIHGFEDCBA" if c not in chains)
            last = max(serials) if serials else 0
            resseqs = [int(l[22:26]) for l in body if l.startswith(("ATOM", "HETATM"))]
            last_res = max(resseqs) if resseqs else 0
            if body and not body[-1].startswith("TER"):
                body.append("TER")
            lig2, con2 = ligand_lines(atoms, bonds, a.resname, chain, start_serial=last + 2, resseq=last_res + 2)
            with open(f"ref_complex_pdb/{t}_complex.pdb", "w") as f:
                f.write("\n".join(body + lig2 + ["TER"] + con2 + ["END"]) + "\n")
            ok += 1
        except Exception as e:
            fail.append((t, str(e)))
    print(f"success: {ok}/{len(ids)}  failure: {len(fail)}")
    for t, e in fail:
        print("  FAIL", t, e)

if __name__ == "__main__":
    main()
