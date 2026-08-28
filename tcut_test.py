"""Is the d-steeper-than-u ordering acceptable if the fit is restricted to the
region where the handbag/twist-3 description is meant to hold, -t <= TCUT?"""
import math, os, numpy as np
from scipy.optimize import least_squares
import fit_slope as F
from reparam import slopes
TCUT=float(os.environ.get("TCUT","1.2"))
F.PI0=[d for d in F.PI0 if -d["t"]<=TCUT]; F.ETA=[d for d in F.ETA if -d["t"]<=TCUT]
F.BSA_PTS=[q for q in F.BSA_PTS if q[2]<=TCUT]; F.EBSA=[q for q in F.EBSA if q[2]<=TCUT]
F.C12_PTS=[q for q in F.C12_PTS if q[2]<=TCUT]
F.EG1_PTS=[q for q in F.EG1_PTS if q[3]<=TCUT]
P0=np.load("fitpar_slope.npy")
W=30.0
def pen(p,order):
    r=[]
    for x in (0.10,0.60):
        s=slopes(p,x)
        for n in ("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00"): r.append(W*min(0.0,s[n]))
        if order: r.append(W*min(0.0,s["H_T^d"]-s["H_T^u"]))
    return r
def run(tag,order):
    def resid(p):
        out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        r.append((p[13]-0.54)/0.15); r+=pen(p,order)
        return np.array(r)
    best=None
    for bd,bpd in ((None,None),(2.0,0.0),(3.0,-1.2)):
        q=np.array(P0)
        if bd is not None: q[5],q[6]=bd,bpd
        try: r=least_squares(resid,np.clip(q,F.LO,F.HI),bounds=(F.LO,F.HI),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        cc=float(np.sum(r.fun**2))
        if best is None or cc<best[0]: best=(cc,r.x)
    p=best[1]; out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    print(f"{tag:>10}  chi2 {tot:8.1f}/{n} pts  chi2/ndf {tot/(n-27):.4f}   "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f}  eta {sum(x*x for x in out['eta']):6.1f}   "
          f"slopes x=0.25: u {slopes(p,0.25)['H_T^u']:.2f} d {slopes(p,0.25)['H_T^d']:.2f}   "
          f"R_HT {p[4]/p[0]:+.3f}  R_ET {p[13]:.3f}")
    return tot,p
print(f"fit restricted to -t <= {TCUT}")
t0,_=run("slope>=0",False)
t1,p1=run("d>=u",True)
print(f"\ncost of the ordering inside the validity region: d_chi2 = {t1-t0:+.1f}"
      f"   (full t range it was +21.3)")
np.save(f"fitpar_order_tcut{TCUT}.npy",p1)
