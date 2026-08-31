"""Extract sigma_U, sigma_TT, sigma_LT on the PUBLISHED bin kinematics.

Each phi point of the published supplemental is assigned to the nearest
published (Q2,xB,t) bin centre -- no nominal t-centres, no width windows, so
the bin membership is the collaboration's, not ours.  Two extractions per bin:

    original : r = 1                  -> must reproduce the published table
    refitted : r = eta_new/eta_pub    -> our improved radiative correction

Output: sf_<ch>_both.txt with the two sets side by side.
"""
import sys
import numpy as np

MP, EBEAM = 0.938272, 5.75
PUB = "/Users/vpk/OneDrive/My_Publications"

def epsilon(q2, xb, e=EBEAM):
    nu = q2/(2*MP*xb); y = nu/e; g2 = q2/nu**2
    return (1 - y - y*y*g2/4) / (1 - y + y*y/2 + y*y*g2/4)

def wfit(ph, y, ey):
    X = np.column_stack([np.ones_like(ph), np.cos(np.radians(ph)),
                         np.cos(2*np.radians(ph))])
    W = 1/ey**2
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    p = cov @ (X.T @ (y*W))
    return p, np.sqrt(np.diag(cov)), np.sum(W*(y - X@p)**2), len(y)-3

CH, RCF, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
if CH == "pi0":
    d = np.loadtxt("/Users/vpk/exclurad-2026/refit/xs_table.txt")
    pf = f"{PUB}/2014_PRC90_025205_pi0_1405.0988/strfun_pi0.data"
else:
    d = np.loadtxt("/Users/vpk/eta-rc-2026/data/xs_corrected_eta_v020.txt")
    pf = f"{PUB}/2017_PRC95_035202_eta_1703.06982/strfun_eta.data"
q2, xb, t, phi, sig, stat, syst = (d[:, i] for i in range(7))

P = np.array([[float(z) for z in l.split()[:12]] for l in open(pf)
              if len(l.split()) >= 12])                 # published bin table

# --- RC ratio, matched point by point as in refit_phi.py -------------------
g = np.loadtxt(RCF)
r = np.ones(len(d)); nm = 0
for i in range(len(d)):
    m = ((np.abs(g[:,3]-phi[i]) < 1.0) & (np.abs(g[:,2]-t[i]) < 0.06)
         & (np.abs(g[:,0]-q2[i]) < 0.10) & (np.abs(g[:,1]-xb[i]) < 0.012))
    if m.sum():
        j = np.argmax(m)
        if np.isfinite(g[j,6]) and g[j,6] > 0:
            r[i] = g[j,6]; nm += 1
print(f"{CH}: RC ratio matched for {nm}/{len(d)} phi points")

# --- assign every phi point to a published bin, within its own (Q2,xB) group --
# The published SF table has fewer bins than the phi supplemental.  Rows of a
# (Q2,xB) group whose t bin was not published must be DROPPED, never folded
# into a neighbouring group -- that would shrink the errors below the published
# ones and quietly change the data.
groups = {}
for i in range(len(d)):
    for k in groups:
        if abs(q2[i]-k[0]) < 0.06 and abs(xb[i]-k[1]) < 0.008:
            groups[k].append(i); break
    else:
        groups[(q2[i], xb[i])] = [i]

own = np.full(len(d), -1, int)
norphan = 0
for k, idx in groups.items():
    idx = np.array(idx)
    cand = np.where((np.abs(P[:,0]-np.mean(q2[idx])) < 0.12)
                    & (np.abs(P[:,1]-np.mean(xb[idx])) < 0.015))[0]
    if len(cand) == 0:
        norphan += len(idx); continue
    for i in idx:
        j = cand[int(np.argmin(np.abs(P[cand,2]-t[i])))]
        if abs(P[j,2]-t[i]) < 0.06:
            own[i] = j
        else:
            norphan += 1
print(f"{CH}: {len(groups)} (Q2,xB) groups; "
      f"{int((own>=0).sum())} phi points assigned, {norphan} in unpublished t bins")

rows = []
for b in range(len(P)):
    s = np.where(own == b)[0]
    if len(s) < 4:                      # need at least one dof for A,B,C
        print(f"   bin {P[b,:3]} has only {len(s)} phi points -- skipped")
        continue
    eps = epsilon(np.mean(q2[s]), np.mean(xb[s]))
    fL, fT, fI = 2*np.pi, 2*np.pi/eps, 2*np.pi/np.sqrt(2*eps*(1+eps))
    out = [P[b,0], P[b,1], P[b,2], eps, len(s),
           np.mean(q2[s]), np.mean(xb[s]), np.mean(t[s])]
    # Systematics are the PUBLISHED ones, carried as a relative error.  Our own
    # propagation would have to assume how the phi points are correlated; the
    # collaboration already did that work, and on published bins their numbers
    # transfer one to one.
    for rr in (np.ones(len(s)), r[s]):
        p, e, c, n = wfit(phi[s], sig[s]/rr, stat[s]/rr)
        val = (fL*p[0], fI*p[1], fT*p[2])
        sy  = [abs(val[k]) * (P[b,5+3*k]/abs(P[b,3+3*k]) if P[b,3+3*k] else 0.0)
               for k in range(3)]
        out += [c/max(n,1),
                fL*p[0], fL*e[0], sy[0],
                fI*p[1], fI*e[1], sy[1],
                fT*p[2], fT*e[2], sy[2]]
    rows.append(out)
rows = np.array(rows)
np.savetxt(OUT, rows, fmt="%12.6f", header=(
    "Q2 xB t eps npts Q2mean xBmean tmean | chi2ndf sigU stat syst sigLT stat syst sigTT stat syst  (ORIGINAL, r=1)"
    " | chi2ndf sigU stat syst sigLT stat syst sigTT stat syst  (REFITTED, new RC)"))
print(f"{CH}: {len(rows)} published bins -> {OUT}")

# --- validation against the published numbers ------------------------------
lab = ("sigma_U", "sigma_LT", "sigma_TT")
for k, col in enumerate((9, 12, 15)):
    pv, ov = P[:len(rows), 3+3*k], rows[:, col]
    if len(pv) == len(ov):
        rel = (ov-pv)/np.maximum(np.abs(pv), 1e-9)
        print(f"   {lab[k]:9s} r=1 vs published: median {np.median(rel):+.4f}  "
              f"max |diff| {np.abs(rel).max():.4f}  ({np.sum(np.abs(rel)>0.02)} bins off >2%)")
