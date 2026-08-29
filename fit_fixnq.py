"""The production fit leaves HT_d nQ pinned on its lower bound (-6), i.e. the
d-quark H_T falling as Q^-3 relative to u, with a normalisation error of 70%.
Two repairs, both bin-averaged and with the same slope constraints:
    tied   nQ_d = nQ_u   (25 par)
    loose  nQ_d free but bounded to [-3, 3]  (26 par)
"""
import math, numpy as np
from scipy.optimize import least_squares
import fit_binned as FB, fit_slope as F
from reparam import slopes
P=np.load("fitpar_amp2026s.npy")
def run(tag, tie, lo7=-6.0):
    free=[i for i in range(27) if i!=12 and not (tie and i==7)]
    def expand(x):
        p=np.zeros(27); p[free]=x; p[12]=0.0
        if tie: p[7]=p[3]
        return p
    def resid(x):
        p=expand(x); out=FB.blocks(p,True); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        r.append((p[13]-0.54)/0.15); r+=FB.pen(p,True)
        return np.array(r)
    LO=np.array(F.LO); HI=np.array(F.HI); LO[7]=lo7; HI[7]=abs(lo7)
    best=None
    for seed in (P, np.where(np.arange(27)==7, -1.0, P), np.where(np.arange(27)==7, 1.0, P)):
        s=np.array(seed); s[7]=max(min(s[7],HI[7]),LO[7])
        r=least_squares(resid,np.clip(s[free],LO[free],HI[free]),bounds=(LO[free],HI[free]),
                        x_scale='jac',xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=FB.blocks(p,True)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    print(f"{tag:>8} {len(free):4d} {tot:9.1f} {tot/(n-len(free)):9.4f}  nQ_u {p[3]:+6.2f} nQ_d {p[7]:+6.2f}"
          f"  N_d/N_u {p[4]/p[0]:+.3f}  R_ET {p[13]:.3f} ({(p[13]-0.54)/0.15:+.1f} sig)"
          f"  pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f}",flush=True)
    np.save(f"fitpar_{tag}.npy",p); return tot
print(f"{'variant':>8} {'par':>4} {'chi2':>9} {'chi2/ndf':>9}")
run("asis",False,-6.0)
run("tied",True)
run("loose",False,-3.0)
