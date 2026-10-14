#!/usr/bin/env python3

import glob, os, collections

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)

files = sorted(glob.glob("ref/*_ligand.mol2"))
name_cnt, atype_cnt = collections.Counter(), collections.Counter()
name_files, atype_files = collections.defaultdict(set), collections.defaultdict(set)
n_atoms = 0

for path in files:
    tid = os.path.basename(path)[:-len("_ligand.mol2")]
    sec = None
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
                name_cnt[name] += 1; name_files[name].add(tid)
                atype_cnt[atype] += 1; atype_files[atype].add(tid)
                n_atoms += 1

def dump(fn, cnt, fset, label):
    with open(fn, "w") as f:
        f.write(f"# {label}: {len(cnt)} unique values, {n_atoms} atoms, {len(files)} mol2 files\n")
        f.write(f"# value\tn_atoms\tn_files\n")
        for k in sorted(cnt):
            f.write(f"{k}\t{cnt[k]}\t{len(fset[k])}\n")

dump("name.log", name_cnt, name_files, "atom name")
dump("atype.log", atype_cnt, atype_files, "atom type")
print(f"files={len(files)} atoms={n_atoms} unique_name={len(name_cnt)} unique_atype={len(atype_cnt)}")
