#!/usr/bin/env python3
"""Prepare Amber (prmtop/inpcrd) inputs for 3D-RISM + placevent from protein/<ID>.pdb.

* standard residues : ff99SB (lab convention, see ~/rism/leaprc.ff99SB); hydrogens are
  rebuilt by tleap, His/Cys protonation states are taken from `reduce -BUILD`
* non-protein / modified residues (DR7, ALY, M3L, PTR, 2A1, SO4, PO4, ACT) : GAFF
  atom types + AM1-BCC charges (antechamber/sqm) ; parmchk2 for missing terms
* Zn2+/Ca2+ : Li-Merz 12-6 (SPC/E) divalent ions (matches the SPC/E solvent of the 3D-RISM run)
* missing atoms inside residues are rebuilt by tleap; missing residues are not built.

usage: prepare_placevent_inputs.py [ID ...]   (default: every file in protein/)
"""
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pdbtools import Atom, read_atoms, group_residues, bonded, STD_AA  # noqa: E402
import param_lib as pl  # noqa: E402

OUT = os.path.join(HERE, "rism_input")
AMBERHOME = pl.AMBERHOME
BIN = pl.BIN
HET_DICT = os.path.join(AMBERHOME, "dat", "reduce_wwPDB_het_dict.txt")

IN_CHAIN = {"ALY": 0, "M3L": 1, "PTR": -2, "2A1": 0}           # formal charge
SEPARATE = {"DR7": 0, "SO4": -2, "PO4": -3, "ACT": -1}           # PO4 modelled fully deprotonated
SHARED = {"SO4", "PO4", "ACT"}                                   # geometry independent -> cache
IONS = {"ZN", "CA"}
KEEP_REDUCE_H = True   # keep reduce-optimised hydrogens of standard residues
BREAK_DIST = 2.2   # C(i)-N(i+1) distance (A) above which a chain break is assumed


def sh(cmd, cwd=None, out=None):
    r = subprocess.run(cmd, cwd=cwd, env=pl.ENV, capture_output=True, text=True)
    if out:
        open(out, "w").write(r.stdout)
    return r


def dist(a, b):
    return float(np.linalg.norm(a.xyz - b.xyz))


def heavy_template(resname):
    """Heavy atom names of a het residue from the reduce dictionary (None if absent)."""
    names, on = [], False
    for l in open(HET_DICT):
        if l.startswith("RESIDUE   %-3s" % resname):
            on = True
            continue
        if on:
            if l.startswith("CONECT"):
                n = l.split()[1]
                if not n.startswith("H"):
                    names.append(n)
            elif l.startswith("END"):
                break
    return names or None


def his_state(R):
    names = {a.name for a in R}
    d, e = "HD1" in names, "HE2" in names
    return "HIP" if d and e else "HID" if d else "HIE"


def process(pid, log):
    wd = os.path.join(OUT, pid)
    os.makedirs(os.path.join(wd, "params"), exist_ok=True)
    src = os.path.join(HERE, "protein", pid + ".pdb")

    # 1. hydrogens / His flips etc. with reduce
    hpdb = os.path.join(wd, pid + "_reduce.pdb")
    env = dict(pl.ENV, REDUCE_HET_DICT=HET_DICT)
    r = subprocess.run([os.path.join(BIN, "reduce"), "-BUILD", "-NUClear", "-Quiet", src], env=env,
                       capture_output=True, text=True)
    open(hpdb, "w").write(r.stdout)
    atoms = read_atoms(hpdb)
    res = group_residues(atoms)
    keys = list(res)
    zn = [a for a in atoms if a.element.upper() in ("ZN", "CA") and a.resname in IONS]

    # 2. residue naming (CYS/HIS states, disulfides, metal sites)
    rename = {}
    sg = [(rid, next(a for a in R if a.name == "SG")) for rid, R in res.items()
          if R[0].resname == "CYS" and any(a.name == "SG" for a in R)]
    for i, (r1, a1) in enumerate(sg):
        for r2, a2 in sg[i + 1:]:
            if dist(a1, a2) < 2.4:
                rename[r1] = rename[r2] = "CYX"
        if r1 not in rename and any(dist(a1, m) < 3.0 for m in zn):
            rename[r1] = "CYM"
        elif r1 not in rename:
            rename[r1] = "CYS" if any(a.name == "HG" for a in res[r1]) else "CYS"
    drop_h = {}
    for rid, R in res.items():
        if R[0].resname == "HIS":
            nm = his_state(R)
            for a in R:
                if a.name in ("ND1", "NE2") and any(dist(a, m) < 2.8 for m in zn):
                    # metal-bound N must be unprotonated
                    hn = "HD1" if a.name == "ND1" else "HE2"
                    if hn in {x.name for x in R}:
                        drop_h.setdefault(rid, set()).add(hn)
            gone = drop_h.get(rid, set())
            has_d = "HD1" in {x.name for x in R} and "HD1" not in gone
            has_e = "HE2" in {x.name for x in R} and "HE2" not in gone
            rename[rid] = "HIP" if has_d and has_e else "HID" if has_d else "HIE"

    # 3. parameterise non-protein parts
    templ, notes = {}, []
    for idx, rid in enumerate(keys):
        R = res[rid]
        rn = R[0].resname
        if rn in IN_CHAIN and rn not in templ:
            pass
    params = {}   # resname -> dict(mol2, frcmod, head, tail)
    for idx, rid in enumerate(keys):
        R = res[rid]
        rn = R[0].resname
        if rn in params or (rn not in IN_CHAIN and rn not in SEPARATE):
            continue
        tpl = heavy_template(rn)
        missing = sorted(set(tpl) - {a.name for a in R} - {'OXT'}) if tpl else []
        if missing:
            notes.append("%s %s%d is missing heavy atoms %s" % (rn, rid[0], rid[1], missing))
        pdir = os.path.join(wd, "params", rn)
        if rn in SHARED:
            sd = os.path.join(OUT, "_shared", rn)
            if os.path.exists(os.path.join(sd, rn + ".mol2")):
                shutil.copytree(sd, pdir, dirs_exist_ok=True)
                params[rn] = dict(mol2=os.path.join(pdir, rn + ".mol2"), frcmod=os.path.join(pdir, rn + ".frcmod"))
                continue
        os.makedirs(pdir, exist_ok=True)
        if rn in IN_CHAIN:
            P = res[keys[idx - 1]] if idx > 0 and keys[idx - 1][0] == rid[0] else []
            N = res[keys[idx + 1]] if idx + 1 < len(keys) and keys[idx + 1][0] == rid[0] else []
            frag = list(R) + pl.backbone_neighbor(R, P, "prev") + pl.backbone_neighbor(R, N, "next")
            frag = pl.cap_open_bonds(frag, atoms)
            q = IN_CHAIN[rn]
            m, _ = pl.antechamber_bcc(pdir, rn, frag, q)
            mol2, q0 = pl.residue_mol2(pdir, m, R, rn, q, retype_backbone=True)
            has_tail = any(a.name == "C" for a in R) and any(a.name == "O" for a in R) and rn != "2A1"
            params[rn] = dict(head="N", tail="C" if has_tail else None)
            notes.append("%s: capped-fragment AM1-BCC (nc=%+d), residue charge before renormalisation %.3f" % (rn, q, q0))
        else:
            q = SEPARATE[rn]
            m, _ = pl.antechamber_bcc(pdir, rn, R, q)
            _, at, bo = pl.read_mol2(m)
            for a, r_ in zip(at, R):
                a[1] = r_.name
            q0 = pl.normalise_charges(at, q)
            mol2 = os.path.join(pdir, rn + ".mol2")
            pl.write_mol2(mol2, rn, at, bo)
            params.setdefault(rn, {})
            notes.append("%s: AM1-BCC (nc=%+d), charge before renormalisation %.3f" % (rn, q, q0))
        attn = pl.parmchk(pdir, rn + ".mol2", rn + ".frcmod")
        params[rn].update(mol2=mol2, frcmod=os.path.join(pdir, rn + ".frcmod"))
        if attn:
            notes.append("%s: %d parmchk2 terms flagged ATTN" % (rn, attn))
        if rn in SHARED:
            shutil.copytree(pdir, os.path.join(OUT, "_shared", rn), dirs_exist_ok=True)

    # 4. write tleap input PDB
    lines, prev = [], None
    serial = 0
    for idx, rid in enumerate(keys):
        R = res[rid]
        rn = R[0].resname
        poly = rn in STD_AA or rn == "ACE"
        if prev is not None:
            prid, pR = prev
            brk = rid[0] != prid[0]
            ppoly = pR[0].resname in STD_AA or pR[0].resname in ("ACE",) or pR[0].resname in IN_CHAIN
            cpoly = poly or rn in IN_CHAIN
            if not brk and ppoly and cpoly:
                c = next((a for a in pR if a.name == "C"), None)
                n = next((a for a in R if a.name == "N"), None)
                if c is not None and n is not None and dist(c, n) > BREAK_DIST:
                    brk = True
                    notes.append("chain break between %s%d and %s%d (%.1f A) -> TER"
                                 % (prid[0], prid[1], rid[0], rid[1], dist(c, n)))
            elif not (ppoly and cpoly):
                brk = True
            if brk:
                lines.append("TER\n")
        std_like = rn in STD_AA or rn == "ACE"
        out_rn = rename.get(rid, rn)
        for a in R:
            if std_like and a.is_h() and (not KEEP_REDUCE_H or rn == "ACE"):
                continue
            if rid in drop_h and a.name in drop_h[rid]:
                continue
            if rn == "CYS" and out_rn in ("CYX", "CYM") and a.name == "HG":
                continue
            serial += 1
            nm = a.name if len(a.name) >= 4 else (" " + a.name).ljust(4)
            lines.append("%-6s%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (
                "ATOM" if std_like or rn in IN_CHAIN else "HETATM", serial % 100000, nm, out_rn, a.chain,
                a.resseq % 10000, a.xyz[0], a.xyz[1], a.xyz[2], a.element))
        prev = (rid, R)
    lines.append("TER\nEND\n")
    open(os.path.join(wd, "input.pdb"), "w").writelines(lines)

    # 5. tleap
    li = ["source oldff/leaprc.ff99SB", "source leaprc.gaff", "loadamberparams frcmod.ions234lm_126_spce",
          "loadOff atomic_ions.lib"]
    for rn, p in params.items():
        li.append("loadamberparams %s" % os.path.relpath(p["frcmod"], wd))
        li.append("%s = loadmol2 %s" % (rn, os.path.relpath(p["mol2"], wd)))
        if "head" in p and p["head"]:
            li.append("set %s head %s.1.%s" % (rn, rn, p["head"]))
            if p.get("tail"):
                li.append("set %s tail %s.1.%s" % (rn, rn, p["tail"]))
    li += ["mol = loadpdb input.pdb", "check mol", "charge mol",
           "savepdb mol %s.pdb" % pid, "saveamberparm mol %s.prmtop %s.rst7" % (pid, pid), "quit"]
    open(os.path.join(wd, "leap.in"), "w").write("\n".join(li) + "\n")
    r = subprocess.run([os.path.join(BIN, "tleap"), "-f", "leap.in"], cwd=wd, env=pl.ENV,
                       capture_output=True, text=True)
    open(os.path.join(wd, "leap.log"), "w").write(r.stdout + r.stderr)
    txt = r.stdout
    added = len(re.findall(r"Added missing heavy atom", txt))
    charge = re.findall(r"Total unperturbed charge:\s*(-?[\d.]+)", txt)
    errs = [l for l in txt.splitlines() if re.search(r"FATAL|Could not find|Parameter file was not saved|Unknown residue|ERROR", l)]
    ok = os.path.exists(os.path.join(wd, pid + ".prmtop")) and not errs
    msg = "%s: %s netcharge=%s added-heavy-atoms=%d gaps=%d" % (
        pid, "OK" if ok else "FAILED", charge[-1] if charge else "?", added,
        sum(1 for n in notes if n.startswith("chain break")))
    log.write(msg + "\n")
    for n in notes:
        log.write("    " + n + "\n")
    for e in errs[:10]:
        log.write("    leap: " + e + "\n")
    log.flush()
    print(msg)
    return ok


def main():
    ids = sys.argv[1:] or sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(HERE, "protein", "*.pdb")))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "summary.txt"), "a") as log:
        for pid in ids:
            process(pid, log)


if __name__ == "__main__":
    main()
