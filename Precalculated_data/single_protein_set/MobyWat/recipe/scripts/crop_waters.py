"""Keep protein + ions + the water molecules closest to the protein so that the system has <= max_atoms atoms.
Needed because MobyWat selects atoms by PDB serial number and PDB serials wrap after 99999.
usage: crop_waters.py system_tpy.pdb out.ndx out_info.txt [max_atoms=99999]
Distance = minimum image distance of water O to any protein atom, frame 0 of the fitted trajectory."""
import sys
import numpy as np

pdb, ndx, info = sys.argv[1:4]
maxat = int(sys.argv[4]) if len(sys.argv) > 4 else 99999
box = None
res, xyz = [], []
for l in open(pdb):
    if l.startswith("CRYST1"):
        box = np.array([float(l[6:15]), float(l[15:24]), float(l[24:33])])
    elif l.startswith(("ATOM", "HETATM")):
        res.append(l[17:20].strip())
        xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    elif l.startswith("ENDMDL"):
        break
res = np.array(res); xyz = np.array(xyz)
prot = ~np.isin(res, ["SOL", "NA", "CL"])
ion = np.isin(res, ["NA", "CL"])
ow = np.where(res == "SOL")[0][::3]            # OW is first atom of each SOL
P = xyz[prot]
d = np.empty(len(ow))
for i in range(0, len(ow), 500):
    diff = xyz[ow[i:i + 500], None, :] - P[None, :, :]
    diff -= box * np.round(diff / box)
    d[i:i + 500] = np.sqrt((diff ** 2).sum(-1)).min(1)
nwat = (maxat - prot.sum() - ion.sum()) // 3
nwat = min(nwat, len(ow))
keep_w = np.sort(np.argsort(d)[:nwat])
keep = np.zeros(len(res), bool)
keep[prot] = True; keep[ion] = True
for w in keep_w:
    keep[ow[w]:ow[w] + 3] = True
idx = np.where(keep)[0] + 1
with open(ndx, "w") as f:
    f.write("[ Crop ]\n")
    for i in range(0, len(idx), 15):
        f.write(" ".join(str(x) for x in idx[i:i + 15]) + "\n")
with open(info, "w") as f:
    f.write(f"atoms_total_before\t{len(res)}\natoms_kept\t{len(idx)}\nprotein_atoms\t{int(prot.sum())}\nions_kept\t{int(ion.sum())}\n"
            f"waters_total\t{len(ow)}\nwaters_kept\t{nwat}\nmax_kept_water_distance_A\t{d[keep_w].max():.2f}\nmin_dropped_water_distance_A\t{np.sort(d)[nwat] if nwat < len(ow) else float('nan'):.2f}\n")
print(open(info).read())
