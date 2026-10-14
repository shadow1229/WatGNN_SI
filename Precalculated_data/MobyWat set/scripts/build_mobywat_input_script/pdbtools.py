"""Small PDB helpers shared by the placevent input scripts."""
import numpy as np

COV_R = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07}
STD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
          "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


class Atom:
    def __init__(self, line):
        self.record = line[:6]
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21]
        self.resseq = int(line[22:26])
        self.icode = line[26]
        self.xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        el = line[76:78].strip() if len(line) >= 78 else ""
        if not el:
            el = "".join(c for c in self.name if c.isalpha())[:1]
        self.element = el.capitalize()

    @property
    def rid(self):
        return (self.chain, self.resseq, self.icode)

    def is_h(self):
        return self.element in ("H", "D")


def read_atoms(path):
    with open(path) as f:
        return [Atom(l.rstrip("\n").ljust(80)) for l in f if l.startswith(("ATOM  ", "HETATM"))]


def group_residues(atoms):
    """Ordered dict rid -> list of atoms (file order)."""
    res = {}
    for a in atoms:
        res.setdefault(a.rid, []).append(a)
    return res


def write_pdb(path, atoms, resname=None):
    """Write atoms with sequential serials; atom names are used as given."""
    with open(path, "w") as f:
        for i, a in enumerate(atoms, 1):
            nm = a.name if len(a.name) >= 4 else (" " + a.name).ljust(4)
            f.write("%-6s%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (
                "HETATM", i, nm, resname or a.resname, a.chain, a.resseq % 10000,
                a.xyz[0], a.xyz[1], a.xyz[2], a.element))
        f.write("END\n")


def bonded(a, b, tol=0.45):
    d = np.linalg.norm(a.xyz - b.xyz)
    return 0.4 < d < COV_R.get(a.element, 0.8) + COV_R.get(b.element, 0.8) + tol
