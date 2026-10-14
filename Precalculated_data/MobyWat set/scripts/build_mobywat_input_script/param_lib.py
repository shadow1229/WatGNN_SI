"""GAFF / AM1-BCC parameterisation of non-protein residues with antechamber."""
import os
import subprocess
import numpy as np
from pdbtools import Atom, bonded, write_pdb, COV_R

AMBERHOME = os.environ.get("AMBERHOME", "/opt/anaconda3/envs/AmberTools23")
BIN = os.path.join(AMBERHOME, "bin")
ENV = dict(os.environ, AMBERHOME=AMBERHOME, PATH=BIN + os.pathsep + os.environ.get("PATH", ""))
CH_LEN = {"C": 1.09, "N": 1.01, "O": 0.96}

# ff14SB backbone types used for in-chain modified residues
BB_TYPES = {"N": "N", "H": "H", "CA": "CT", "HA": "H1", "HA2": "H1", "HA3": "H1",
            "C": "C", "O": "O"}


def run(cmd, cwd):
    r = subprocess.run(cmd, cwd=cwd, env=ENV, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("%s failed:\n%s\n%s" % (cmd[0], r.stdout[-1500:], r.stderr[-1500:]))
    return r


def make_atom(name, element, xyz):
    a = Atom("HETATM    1  X   UNK A   1       0.000   0.000   0.000  1.00  0.00          X ")
    a.name, a.element, a.xyz = name, element, np.asarray(xyz, float)
    return a


def cap_open_bonds(frag, full):
    """Add H atoms along every bond from a fragment atom to an atom outside it."""
    ids = {id(a) for a in frag}
    caps = []
    for a in frag:
        if a.is_h():
            continue
        for b in full:
            if id(b) in ids or b is a or not bonded(a, b):
                continue
            v = b.xyz - a.xyz
            v /= np.linalg.norm(v)
            caps.append(make_atom("HC", "H", a.xyz + v * CH_LEN.get(a.element, 1.0)))
    return frag + caps


def backbone_neighbor(res_atoms, nb_atoms, which):
    """Return the neighbour's cap atoms (acetyl-like for prev, N-methylamide-like for next)."""
    if not nb_atoms:
        return []
    names = {a.name: a for a in nb_atoms}
    if which == "prev":
        if "C" not in names or "N" not in {a.name for a in res_atoms}:
            return []
        n = next(a for a in res_atoms if a.name == "N")
        if not bonded(names["C"], n):
            return []
        if nb_atoms[0].resname == "ACE":
            return list(nb_atoms)
        keep = {"CA", "C", "O", "HA", "HA2", "HA3"}
    else:
        if "N" not in names or "C" not in {a.name for a in res_atoms}:
            return []
        c = next(a for a in res_atoms if a.name == "C")
        if not bonded(c, names["N"]):
            return []
        keep = {"N", "H", "CA", "HA", "HA2", "HA3", "CD"} if nb_atoms[0].resname == "PRO" \
            else {"N", "H", "CA", "HA", "HA2", "HA3"}
    return [a for a in nb_atoms if a.name in keep]


def read_mol2(path):
    lines = open(path).read().splitlines()
    sec, head, atoms, bonds = None, [], [], []
    for l in lines:
        if l.startswith("@<TRIPOS>"):
            sec = l.strip()
            continue
        if sec == "@<TRIPOS>MOLECULE":
            head.append(l)
        elif sec == "@<TRIPOS>ATOM" and l.strip():
            t = l.split()
            atoms.append([int(t[0]), t[1], float(t[2]), float(t[3]), float(t[4]), t[5], t[6], t[7], float(t[8])])
        elif sec == "@<TRIPOS>BOND" and l.strip():
            bonds.append(l.split())
    return head, atoms, bonds


def write_mol2(path, resname, atoms, bonds):
    with open(path, "w") as f:
        f.write("@<TRIPOS>MOLECULE\n%s\n %d %d 1 0 0\nSMALL\nUSER_CHARGES\n\n" % (resname, len(atoms), len(bonds)))
        f.write("@<TRIPOS>ATOM\n")
        for i, a in enumerate(atoms, 1):
            f.write("%7d %-8s %9.4f %9.4f %9.4f %-5s %3d %-4s %10.6f\n" % (i, a[1], a[2], a[3], a[4], a[5], 1, resname, a[8]))
        f.write("@<TRIPOS>BOND\n")
        for i, b in enumerate(bonds, 1):
            f.write("%6d %5s %5s %s\n" % (i, b[1], b[2], b[3]))
        f.write("@<TRIPOS>SUBSTRUCTURE\n     1 %s 1 TEMP 0 **** **** 0 ROOT\n" % resname)


def antechamber_bcc(workdir, name, atoms, net_charge, ftype="gaff"):
    """Run AM1-BCC on `atoms` (unique names assigned internally). Returns mol2 path."""
    os.makedirs(workdir, exist_ok=True)
    counts = {}
    tmp = []
    for a in atoms:
        counts[a.element] = counts.get(a.element, 0) + 1
        b = make_atom("%s%d" % (a.element, counts[a.element]), a.element, a.xyz)
        tmp.append(b)
    write_pdb(os.path.join(workdir, name + "_in.pdb"), tmp, resname=name[:3].upper())
    run([os.path.join(BIN, "antechamber"), "-i", name + "_in.pdb", "-fi", "pdb", "-o", name + "_bcc.mol2",
         "-fo", "mol2", "-c", "bcc", "-nc", str(net_charge), "-at", ftype, "-rn", name[:3].upper(),
         "-pf", "y", "-s", "2"], workdir)
    return os.path.join(workdir, name + "_bcc.mol2"), [t.name for t in tmp]


FF_TYPES = {"N", "H", "CT", "H1", "C", "O"}  # protein types already defined by ff99SB
BOND_RE = None


def filter_frcmod(src, dst):
    """Drop MASS/NONBON/bonded lines that only involve existing protein atom types,
    so loading the frcmod never overrides ff99SB parameters. Returns #ATTN lines."""
    import re
    pat = re.compile(r"^([A-Za-z0-9*]{1,2} ?(?:-[A-Za-z0-9*]{1,2} ?){1,3})\s")
    out, sec, attn = [], None, 0
    for l in open(src).read().splitlines():
        s = l.strip()
        if s in ("MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"):
            sec = s
            out.append(l)
            continue
        if sec is None or not s:
            out.append(l)
            continue
        if sec in ("MASS", "NONBON"):
            t = [s.split()[0]]
        else:
            m = pat.match(l)
            if not m:
                out.append(l)
                continue
            t = [x.strip() for x in m.group(1).split("-")]
        if all(x in FF_TYPES for x in t):
            continue
        if "ATTN" in l:
            attn += 1
        out.append(l)
    open(dst, "w").write("\n".join(out) + "\n")
    return attn


def parmchk(workdir, mol2, frcmod_out):
    raw = frcmod_out + ".raw"
    run([os.path.join(BIN, "parmchk2"), "-i", mol2, "-f", "mol2", "-o", raw, "-s", "1"], workdir)
    return filter_frcmod(os.path.join(workdir, raw) if not os.path.isabs(raw) else raw,
                         os.path.join(workdir, frcmod_out) if not os.path.isabs(frcmod_out) else frcmod_out)


def normalise_charges(atoms, target):
    q = sum(a[8] for a in atoms)
    d = (target - q) / len(atoms)
    for a in atoms:
        a[8] = round(a[8] + d, 6)
    # absorb rounding error
    atoms[0][8] = round(atoms[0][8] + (target - sum(a[8] for a in atoms)), 6)
    return q


def residue_mol2(workdir, mol2_bcc, resatoms, resname, formal_q, retype_backbone, full_atoms=None):
    """Build the final residue mol2 (original names, protein backbone types if requested)."""
    _, at, bo = read_mol2(mol2_bcc)
    n = len(resatoms)
    res = [list(a) for a in at[:n]]
    for a, r in zip(res, resatoms):
        a[1] = r.name
    if retype_backbone:
        ca = next(r for r in resatoms if r.name == "CA")
        nn = next(r for r in resatoms if r.name == "N")
        types = {r.name: a for a, r in zip(res, resatoms)}
        for a, r in zip(res, resatoms):
            if r.name == "N":
                a[5] = "N"
            elif r.is_h() and bonded(r, nn):
                a[5] = "H"
            elif r.name == "CA":
                a[5] = "CT"
            elif r.is_h() and bonded(r, ca):
                a[5] = "H1"
            elif r.name == "C" and a[5] == "c":
                a[5] = "C"
            elif r.name == "O" and a[5] == "o":
                a[5] = "O"
    q0 = normalise_charges(res, formal_q)
    ids = {a[0] for a in res}
    b = [x for x in bo if int(x[1]) in ids and int(x[2]) in ids]
    path = os.path.join(workdir, resname + ".mol2")
    write_mol2(path, resname, res, b)
    return path, q0
