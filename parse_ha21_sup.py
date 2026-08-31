"""Read table II of the PRL 127 152301 supplemental: the four structure functions
of the nine Hall-A settings, including dsigma_LT', which the main paper only
shows as a figure.

Layout note: the setting is printed down the left margin, xB on the first row of
the block and Q2 on the second, so the block's Q2 is not yet known when its first
data row is read.  Rows are therefore buffered per block and stamped afterwards.
"""
import re
import numpy as np

SRC = "/tmp/ha21sup.txt"          # pdftotext -layout of 2021_HallA_pi0_sup.pdf
NUM = r"-?\d+\.?\d*"
out, buf, xb, q2 = [], [], None, None

def flush():
    if buf and xb is not None and q2 is not None:
        out.extend([[q2, xb] + r for r in buf])
    buf.clear()

for line in open(SRC, encoding="utf-8"):
    m = re.search(r"xB\s*=\s*(" + NUM + ")", line)
    if m:
        flush(); xb, q2 = float(m.group(1)), None
    m = re.search(r"Q2\s*=\s*(" + NUM + ")", line)
    if m: q2 = float(m.group(1))
    m = re.search(r"\[\s*(" + NUM + r")\s*,\s*(" + NUM + r")\s*\]\s+(" + NUM + r")\s+(.*)$", line)
    if not m or xb is None: continue
    v = [float(z) for z in re.findall(NUM, m.group(4))]
    if len(v) != 12: continue
    buf.append([float(m.group(1)), float(m.group(2)), float(m.group(3))] + v)
flush()

d = np.array(out)
np.savetxt("data/halla_y21.data", d, fmt="%9.4f", header=(
    "Hall-A pi0 off the proton, PRL 127 152301 (2021), supplemental table II\n"
    "sigma in nb/GeV^2; t' = tmin - t; <t'> is the event-weighted average of the bin\n"
    "Q2 xB tp_low tp_up tp  sigU stat syst  sigLT stat syst  sigTT stat syst  sigLTp stat syst"))
print(f"parsed {len(d)} rows -> data/halla_y21.data")
print("settings: " + ", ".join(f"({xb:.2f},{q:.2f})x{n}" for (xb, q), n in
      sorted({(k[1], k[0]): int(((d[:,1]==k[1]) & (d[:,0]==k[0])).sum())
              for k in set(zip(d[:,0], d[:,1]))}.items())))

old = np.array([[float(x) for x in l.split()[3:16]] for l in open("data/halla_more.data")
                if l.split() and l.split()[0] == "HallA_y21"])
nm, mx = 0, 0.0
for r in d:
    m = ((np.abs(old[:,0]-r[0]) < 0.02) & (np.abs(old[:,1]-r[1]) < 0.01)
         & (np.abs(old[:,3]-r[5]) < 0.02))
    if m.sum() == 1:
        nm += 1; o = old[m][0]
        mx = max(mx, *(abs(o[i]-r[j]) for i, j in
                       ((3,5),(4,6),(5,7),(6,8),(7,9),(8,10),(9,11),(10,12),(11,13))))
print(f"cross-check against the numbers we already had: {nm}/36 rows matched, "
      f"max |difference| = {mx:.4f}")
s, e = d[:,14], np.hypot(d[:,15], d[:,16])
print(f"sigma_LT' is new: {s.min():.2f} .. {s.max():.2f} nb/GeV^2, median error {np.median(e):.2f}")
print(f"   nonzero at more than 2 sigma: {int((np.abs(s) > 2*e).sum())}/36, "
      f"{int((s > 0).sum())} of 36 positive")
print(f"   chi2 against zero: {float(np.sum((s/e)**2)):.1f}/36")
