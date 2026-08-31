"""Read tables IV, V, VI and VII of PRC 83 025201 (2011).

Four kinematics, seven t' bins each, and four structure functions per bin --
sigma_U, sigma_LT, sigma_TT and sigma_LT', the last of which we did not have.
The two tables are typeset side by side, table VI in the left column of the page
and table VII in the right, so each line is split at column 72 and the two
halves are walked independently.
"""
import re
import numpy as np

SRC = "/tmp/2011.txt"
SPLIT = 72
KIN = {                     # label: (xB, Q2, eps, tmin, column index in its table)
    "Kin3":  (0.369, 2.350, 0.649, 0.173, 0),
    "Kin2":  (0.368, 1.941, 0.769, 0.170, 1),
    "KinX3": (0.335, 2.155, 0.648, 0.137, 0),
    "KinX2": (0.394, 2.073, 0.768, 0.199, 1),
}
# table V: <tmin - t> weighted by the cross section
TP = {"Kin3":  [0.0095, 0.0298, 0.0546, 0.0844, 0.1188, 0.1583, 0.2063],
      "Kin2":  [0.0094, 0.0296, 0.0541, 0.0839, 0.1179, 0.1576, 0.2050],
      "KinX3": [0.0095, 0.0297, 0.0545, 0.0843, 0.1188, 0.1579, 0.2057],
      "KinX2": [0.0094, 0.0296, 0.0542, 0.0840, 0.1181, 0.1579, 0.2051]}

def block_of(seg):
    if "dσT T" in seg: return "TT"
    if "dσT L\x02" in seg or "dσTL’" in seg: return "LTp"
    if "dσT L" in seg: return "LT"
    if "dσT /dt" in seg: return "U"
    return None

ROW = re.compile(r"^\s*(0\.\d{3})\s+(.*)$")
NUM = r"[-−]?\d+"
TRIPLE = re.compile(NUM + r"\s*±\s*(\d+)\s*±\s*(\d+)")

def parse(side):
    cur, out = None, {}
    for line in open(SRC, encoding="utf-8"):
        seg = line[:SPLIT] if side == 0 else line[SPLIT:]
        b = block_of(seg)
        if b: cur = b
        m = ROW.match(seg)
        if not m or cur is None: continue
        vals = [float(z.replace("−", "-")) for z in re.findall(NUM, m.group(2))]
        if len(vals) != 6: continue          # two kinematics x (value, stat, syst)
        out.setdefault(cur, []).append([float(m.group(1))] + vals)
    return out

rows = []
for side, names in ((0, ("Kin3", "Kin2")), (1, ("KinX3", "KinX2"))):
    got = parse(side)
    for b in ("U", "LT", "TT", "LTp"):
        if b not in got or len(got[b]) != 7:
            print(f"  side {side} block {b}: {len(got.get(b, []))} rows (expected 7)")
    for j, name in enumerate(names):
        xB, Q2, eps, tm, _ = KIN[name]
        for i in range(7):
            r = [Q2, xB, eps, tm, TP[name][i]]
            for b in ("U", "LT", "TT", "LTp"):
                r += got[b][i][1 + 3*j : 4 + 3*j]
            rows.append(r)
d = np.array(rows)
np.savetxt("data/halla_y11.data", d, fmt="%9.4f", header=(
    "Hall-A pi0 off the proton, PRC 83 025201 (2011), tables IV-VII\n"
    "four kinematics Kin3 Kin2 KinX3 KinX2, seven t' bins each; sigma in nb/GeV^2\n"
    "t' = tmin - t taken from table V (cross-section weighted); tmin from table IV\n"
    "Q2 xB eps tmin tp  sigU stat syst  sigLT stat syst  sigTT stat syst  sigLTp stat syst"))
print(f"parsed {len(d)} rows -> data/halla_y11.data")
for name in KIN:
    xB, Q2 = KIN[name][:2]
    m = (np.abs(d[:,0]-Q2) < 0.001) & (np.abs(d[:,1]-xB) < 0.001)
    print(f"   {name:6s} Q2={Q2:.3f} xB={xB:.3f}: {int(m.sum())} bins, "
          f"sigma_U {d[m,5].min():.0f}..{d[m,5].max():.0f}")
