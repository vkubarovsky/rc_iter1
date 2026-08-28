"""Iteration 1: RC ratio r = eta(amp2021) / eta(published model), both computed
with exclurad_py so that the CODE cancels in the ratio and only the MODEL
difference survives.

    pi0: published model = pi0.vpk2013   (PRC 90 025205), v_cut = 0.18
    eta: published model = eta.tableB1   (PRC 95 035202), v_cut = 0.094

Output: q2 xb t phi eta_pub eta_new r      (r = eta_new/eta_pub)
Correction applied downstream:  sigma_new = sigma_pub / r.
"""
import math, os, sys
from concurrent.futures import ProcessPoolExecutor
import numpy as np
sys.path.insert(0, "/Users/vpk/exclurad_py")
from exclurad_py import models
from exclurad_py.core.rc import rc_factor
from exclurad_py.core.constants import M_P

CH = sys.argv[1]
JOBS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
CFG = dict(
    pi0=dict(pub="pi0.vpk2013", new="pi0.amp2021", vcut=0.18, m=0.1349768,
             table="/Users/vpk/exclurad-2026/refit/xs_table.txt", cols=(0,1,2,3)),
    eta=dict(pub="eta.tableB1", new="eta.amp2021", vcut=0.094, m=0.547862,
             table="/Users/vpk/eta-rc-2026/data/xs_corrected_eta_v020.txt", cols=(0,1,2,3)),
)[CH]

d = np.loadtxt(CFG["table"])
q2, xb, t, phi = (d[:, c] for c in CFG["cols"])
mx2 = CFG["m"]**2 + CFG["vcut"]
print(f"{CH}: {len(d)} points, v_cut={CFG['vcut']}, MX2cut={mx2:.4f}", flush=True)

def one(i):
    Q2, xB, mt, ph = float(q2[i]), float(xb[i]), float(t[i]), float(phi[i])
    W2 = M_P*M_P + Q2*(1-xB)/xB
    out = []
    for key in (CFG["pub"], CFG["new"]):
        try:
            m = models.get(key, t_nucl=-mt, Ebeam=5.75)
            r = rc_factor(5.75, Q2, W2, -mt, math.radians(ph), m,
                          h=+1, mx2_cut=mx2)
            out.append(float(r["eta"]))
        except Exception:
            out.append(float('nan'))
    return (Q2, xB, mt, ph, out[0], out[1])

if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=JOBS) as ex:
        res = list(ex.map(one, range(len(d)), chunksize=8))
    a = np.array(res)
    ratio = a[:, 5] / a[:, 4]
    a = np.column_stack([a, ratio])
    good = np.isfinite(ratio)
    np.savetxt(f"/Users/vpk/rc_iter1/rc_{CH}.txt", a,
               header="q2 xb t phi eta_pub eta_new r=eta_new/eta_pub",
               fmt="%9.4f %8.4f %7.3f %8.2f %10.5f %10.5f %9.5f")
    print(f"{CH}: {good.sum()}/{len(a)} finite; "
          f"r median {np.nanmedian(ratio):.4f}, "
          f"16-84% [{np.nanpercentile(ratio,16):.4f}, {np.nanpercentile(ratio,84):.4f}]",
          flush=True)
