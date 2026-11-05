#!/usr/bin/env python3
"""Copy from_pm/<target>/dry.pdb -> MobyWat_aligned/dry/<target>.pdb (byte-for-byte) for every target of the single protein set.
usage: copy_dry_pdb.py [from_pm dir] [targets.txt] [output dir]"""
import os, shutil, sys, hashlib

src = sys.argv[1] if len(sys.argv) > 1 else "/Users/shadow1229/claude_project/Mobywat_revision/vesper_mobywat_test/from_pm"
tfile = sys.argv[2] if len(sys.argv) > 2 else "/Users/shadow1229/claude_project/WatGNN/WatGNN_SI/Dataset/single protein set/single_protein_set.txt"
out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "MobyWat_aligned", "dry")
os.makedirs(out, exist_ok=True)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
n = 0
for t in [l.strip() for l in open(tfile) if l.strip()]:
    s, d = os.path.join(src, t, "dry.pdb"), os.path.join(out, f"{t}.pdb")
    shutil.copyfile(s, d)
    assert sha(s) == sha(d), t
    n += 1
print(f"{n} files copied -> {out}")
