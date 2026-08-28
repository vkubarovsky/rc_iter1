"""Refit sigma(phi) = A + B cos(phi) + C cos(2phi) on RC-recorrected 4-fold
cross sections, for pi0 and eta.  Same fit as ~/exclurad-2026/refit/refit.py;
generalised over channel and over the RC-ratio file of the current iteration.

    sigma_new = sigma_pub / r,   r = eta(new model) / eta(published model)

Output: <out>  with columns
    Q2 xB t eps npts chi2ndf  sigU err sigTT err sigLT err
"""
import sys
import numpy as np

MP, EBEAM = 0.938272, 5.75

def epsilon(q2, xb, e=EBEAM):
    nu = q2/(2*MP*xb); y = nu/e; g2 = q2/nu**2
    return (1 - y - y*y*g2/4) / (1 - y + y*y/2 + y*y*g2/4)

def wfit(ph, y, ey):
    X = np.column_stack([np.ones_like(ph), np.cos(np.radians(ph)),
                         np.cos(2*np.radians(ph))])
    W = 1/ey**2
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    p = cov @ (X.T @ (y*W))
    chi2 = np.sum(W*(y - X@p)**2)
    return p, np.sqrt(np.diag(cov)), chi2, len(y)-3

CH   = sys.argv[1]
RCF  = sys.argv[2]
OUT  = sys.argv[3]

if CH == "pi0":
    d = np.loadtxt("/Users/vpk/exclurad-2026/refit/xs_table.txt")
    q2, xb, t, phi, sig, stat = (d[:, i] for i in range(6))
    T_CENTERS = [0.12, 0.175, 0.25, 0.35, 0.49, 0.77, 1.21, 1.71]
else:
    d = np.loadtxt("/Users/vpk/eta-rc-2026/data/xs_corrected_eta_v020.txt")
    q2, xb, t, phi, sig, stat = (d[:, i] for i in range(6))
    T_CENTERS = [0.12, 0.17, 0.25, 0.35, 0.50, 0.80, 1.25, 1.75]

g = np.loadtxt(RCF)                      # q2 xb t phi eta_pub eta_new r
r = np.ones(len(d)); nm = 0
for i in range(len(d)):
    m = ((np.abs(g[:,3]-phi[i]) < 1.0) & (np.abs(g[:,2]-t[i]) < 0.06)
         & (np.abs(g[:,0]-q2[i]) < 0.10) & (np.abs(g[:,1]-xb[i]) < 0.012))
    if m.sum():
        j = np.argmax(m)
        if np.isfinite(g[j,6]) and g[j,6] > 0:
            r[i] = g[j,6]; nm += 1
print(f"{CH}: matched {nm}/{len(d)} points to the RC ratio")

keys = {}
for i in range(len(d)):
    for (kq, kx) in keys:
        if abs(q2[i]-kq) < 0.06 and abs(xb[i]-kx) < 0.008:
            keys[(kq, kx)].append(i); break
    else:
        keys[(q2[i], xb[i])] = [i]

rows = []
for (kq, kx), idx in sorted(keys.items()):
    idx = np.array(idx)
    eps = epsilon(np.mean(q2[idx]), np.mean(xb[idx]))
    for tc in T_CENTERS:
        s = idx[np.abs(t[idx]-tc) < 0.045]
        if len(s) < 8: continue
        ph, y, ey, rr = phi[s], sig[s], stat[s], r[s]
        p, e, c, n = wfit(ph, y/rr, ey/rr)
        sigU  = 2*np.pi*p[0];                        esigU  = 2*np.pi*e[0]
        sigTT = 2*np.pi*p[2]/eps;                    esigTT = 2*np.pi*e[2]/eps
        sigLT = 2*np.pi*p[1]/np.sqrt(2*eps*(1+eps)); esigLT = 2*np.pi*e[1]/np.sqrt(2*eps*(1+eps))
        rows.append([kq, kx, tc, eps, n+3, c/max(n,1),
                     sigU, esigU, sigTT, esigTT, sigLT, esigLT])
rows = np.array(rows)
np.savetxt(OUT, rows, fmt="%9.4f",
           header="Q2 xB t eps npts chi2ndf sigU err sigTT err sigLT err")
print(f"{CH}: {len(rows)} (Q2,xB,t) groups -> {OUT}")
print(f"   median chi2/ndf {np.median(rows[:,5]):.2f}")
