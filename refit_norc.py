"""Same extraction as refit_phi.py but with r = 1: the ORIGINAL structure
functions, obtained from the published phi distributions with the published RC.

Two purposes: it is the 'as published' column of the database, and it validates
our machinery -- these numbers should reproduce the published structure-function
tables, since only the RC differs from the production set.
"""
import sys
import numpy as np
MP, EBEAM = 0.938272, 5.75
def epsilon(q2, xb, e=EBEAM):
    nu = q2/(2*MP*xb); y = nu/e; g2 = q2/nu**2
    return (1 - y - y*y*g2/4) / (1 - y + y*y/2 + y*y*g2/4)
def wfit(ph, y, ey):
    X = np.column_stack([np.ones_like(ph), np.cos(np.radians(ph)), np.cos(2*np.radians(ph))])
    W = 1/ey**2
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    p = cov @ (X.T @ (y*W)); chi2 = np.sum(W*(y - X@p)**2)
    return p, np.sqrt(np.diag(cov)), chi2, len(y)-3
CH, OUT = sys.argv[1], sys.argv[2]
if CH == "pi0":
    d = np.loadtxt("/Users/vpk/exclurad-2026/refit/xs_table.txt")
    T_CENTERS = [0.12, 0.175, 0.25, 0.35, 0.49, 0.77, 1.21, 1.71]
else:
    d = np.loadtxt("/Users/vpk/eta-rc-2026/data/xs_corrected_eta_v020.txt")
    T_CENTERS = [0.12, 0.17, 0.25, 0.35, 0.50, 0.80, 1.25, 1.75]
q2, xb, t, phi, sig, stat = (d[:, i] for i in range(6))
syst = d[:, 6]
keys = {}
for i in range(len(d)):
    for (kq, kx) in keys:
        if abs(q2[i]-kq) < 0.06 and abs(xb[i]-kx) < 0.008:
            keys[(kq, kx)].append(i); break
    else:
        keys[(q2[i], xb[i])] = [i]
rows = []
for (kq, kx), idx in sorted(keys.items()):
    idx = np.array(idx); eps = epsilon(np.mean(q2[idx]), np.mean(xb[idx]))
    for tc in T_CENTERS:
        s = idx[np.abs(t[idx]-tc) < 0.045]
        if len(s) < 8: continue
        p, e, c, n = wfit(phi[s], sig[s], stat[s])
        # systematic on the structure functions: refit with sigma shifted by +syst
        p2, _, _, _ = wfit(phi[s], sig[s]+syst[s], stat[s])
        rows.append([kq, kx, tc, eps, n+3, c/max(n,1),
                     2*np.pi*p[0], 2*np.pi*e[0], 2*np.pi*abs(p2[0]-p[0]),
                     2*np.pi*p[2]/eps, 2*np.pi*e[2]/eps, 2*np.pi*abs(p2[2]-p[2])/eps,
                     2*np.pi*p[1]/np.sqrt(2*eps*(1+eps)), 2*np.pi*e[1]/np.sqrt(2*eps*(1+eps)),
                     2*np.pi*abs(p2[1]-p[1])/np.sqrt(2*eps*(1+eps))])
rows = np.array(rows)
np.savetxt(OUT, rows, fmt="%10.4f",
           header="Q2 xB t eps npts chi2ndf sigU stat syst sigTT stat syst sigLT stat syst")
print(f"{CH}: {len(rows)} bins -> {OUT}   median chi2/ndf {np.median(rows[:,5]):.2f}")
