"""sf_<ch>_iN.txt  ->  strfun format used by the amplitude fit:
Q2 xB |t| sigU stat sys sigLT stat sys sigTT stat sys
Systematics: carry the PUBLISHED relative systematic across (the RC ratio is a
model quantity with no statistical error; the published relative systematics
are preserved by construction of the correction)."""
import sys
import numpy as np
ch, src, out = sys.argv[1], sys.argv[2], sys.argv[3]
pub = np.loadtxt(f"/Users/vpk/pi0_eta_amplitude_model/data/strfun_{ch}.data")
n = np.loadtxt(src)
rows = []
for row in n:
    Q2, xB, mt = row[0], row[1], row[2]
    sU, eU, sTT, eTT, sLT, eLT = row[6], row[7], row[8], row[9], row[10], row[11]
    m = (abs(pub[:,0]-Q2) < 0.08) & (abs(pub[:,1]-xB) < 0.01) & (abs(pub[:,2]-mt) < 0.05)
    if m.sum():
        j = np.argmax(m)
        fU  = pub[j,5]/pub[j,3]  if pub[j,3]  else 0.0
        fLT = pub[j,8]/abs(pub[j,6])  if pub[j,6]  else 0.0
        fTT = pub[j,11]/abs(pub[j,9]) if pub[j,9] else 0.0
    else:
        fU = fLT = fTT = 0.0
    rows.append([Q2, xB, mt, sU, eU, abs(sU*fU),
                 sLT, eLT, abs(sLT*fLT), sTT, eTT, abs(sTT*fTT)])
np.savetxt(out, np.array(rows), fmt="%8.3f")
print(f"{ch}: {len(rows)} rows -> {out}")
