"""Remove crystallographic HOH (and other non-protein records) from a reference PDB."""
import sys

src, dst = sys.argv[1], sys.argv[2]
keep = []
for line in open(src):
    if line.startswith("ATOM"):
        if line[16] not in (" ", "A"):  # keep only the first alternate location
            continue
        keep.append(line[:16] + " " + line[17:])
    elif line.startswith(("TER", "END")):
        keep.append(line)
open(dst, "w").writelines(keep)
n_hoh = sum(1 for l in open(src) if l[17:20] == "HOH")
print(f"removed {n_hoh} HOH; kept {sum(l.startswith('ATOM') for l in keep)} protein atoms")
