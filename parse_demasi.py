"""Read the CLAS database export of the De Masi beam-spin asymmetry.

Sixty measurement blocks, each one (xB, Q2, -t) bin with the asymmetry as a
function of phi.  This replaces the figure extraction we had been using: the
alpha values in clas6_demasi_alu.json were digitised from figure 5 of the 2008
Rapid Communication, and comparing them with the CLAS12 supplemental suggested
they were high by a factor near three.
"""
import re
import numpy as np

SRC = "2024_CLAS6_pi0_BSA.txt"
NUM = r"-?\d+\.?\d*(?:[eE][-+]?\d+)?"
blocks, cur, name, indata = [], [], None, False
for line in open(SRC, encoding="utf-8", errors="replace"):
    m = re.match(r"Measurement\s+(\S+)", line)
    if m:
        if cur: blocks.append((name, cur))
        name, cur, indata = m.group(1), [], False
        continue
    if line.startswith("Data:"): indata = True; continue
    if not indata: continue
    v = re.findall(NUM, line)
    if len(v) == 6: cur.append([float(z) for z in v])
if cur: blocks.append((name, cur))

print(f"{len(blocks)} measurement blocks")
rows = []
for nm, b in blocks:
    d = np.array(b)
    if len(d) < 4: print(f"   {nm}: only {len(d)} phi points, skipped"); continue
    rows.append((nm, d))
print(f"{len(rows)} usable, {sum(len(d) for _, d in rows)} phi points total")
print(f"\n{'block':>8}{'npts':>6}{'<xB>':>8}{'<Q2>':>7}{'<-t>':>7}{'phi range':>16}")
for nm, d in rows[:8]:
    print(f"{nm:>8}{len(d):6d}{d[:,0].mean():8.3f}{d[:,1].mean():7.2f}{d[:,2].mean():7.3f}"
          f"{f'{d[:,3].min():.0f} - {d[:,3].max():.0f}':>16}")
np.save("/tmp/demasi_blocks.npy", np.array([1]))
import pickle
pickle.dump(rows, open("data/demasi_phi.pkl", "wb"))
print("\nsaved -> data/demasi_phi.pkl")
