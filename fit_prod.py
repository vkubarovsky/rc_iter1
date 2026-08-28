"""Production fit: the physical slope constraints, all of them at once.

  (a) b + b' ln xB >= 0 for every block on xB in [0.1, 0.6]     (no growing FFs)
  (b) slope(H_T^d) >= slope(H_T^u)   - the ordering of VPK's global GPD fits
  (c) Ebar_T effective slope b + b' ln xB + 2 b2 t >= 0 out to -t = TMAX,
      i.e. the b2 t^2 curvature may flatten the fall but not reverse it inside
      the generator validity window.
"""
import math, os, numpy as np
from scipy.optimize import least_squares
import fit_slope as F
from reparam import slopes

TMAX=float(os.environ.get("TMAX","2.5"))
OUTP=os.environ.get("OUTP","fitpar_prod.npy")
W=30.0; XS=(0.10,0.60)
NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")

def pen(p, order=True, curv=True):
    r=[]
    for x in XS:
        s=slopes(p,x)
        for n in NAMES: r.append(W*min(0.0, s[n]))
        if order: r.append(W*min(0.0, s["H_T^d"]-s["H_T^u"]))
        if curv:
            for n in ("Ebar_T^u","Ebar_T^d"):
                r.append(W*min(0.0, s[n]-2*p[12]*TMAX))     # slope at t = -TMAX
    return r

def run(tag, order, curv, seed):
    def resid(p):
        out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        r.append((p[13]-0.54)/0.15); r+=pen(p,order,curv)
        return np.array(r)
    best=None
    for mod in (None,(2.0,0.0),(3.0,-1.2)):
        q=np.array(seed)
        if mod: q[5],q[6]=mod
        try: r=least_squares(resid,np.clip(q,F.LO,F.HI),bounds=(F.LO,F.HI),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        cc=float(np.sum(r.fun**2))
        if best is None or cc<best[0]: best=(cc,r.x)
    p=best[1]; out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    s25=slopes(p,0.25)
    print(f"{tag:>26}  chi2 {tot:7.1f}  chi2/ndf {tot/(n-27):.4f}  pi0 {sum(x*x for x in out['pi0']):6.1f} "
          f"eta {sum(x*x for x in out['eta']):6.1f}  b2 {p[12]:+.3f}  ET turnover -t="
          f"{(slopes(p,0.60)['Ebar_T^u']/(2*p[12]) if p[12]>1e-6 else 99):5.2f}  "
          f"u {s25['H_T^u']:.2f} d {s25['H_T^d']:.2f}  R_HT {p[4]/p[0]:+.3f} R_ET {p[13]:.3f} "
          f"pen {sum(x*x for x in pen(p,order,curv)):.2f}")
    return tot,p

P0=np.load("fitpar_slope_order.npy")
o=F.blocks(P0); base=sum(sum(x*x for x in v) for v in o.values())
print(f"reference: d>=u without the curvature constraint, chi2 = {base:.1f}\n")
t1,p1=run("d>=u + curvature(-t<=2.5)",True,True,P0)
t2,p2=run("d>=u + curvature(-t<=3.5)",True,True,P0) if False else (None,None)
np.save(OUTP,p1); print(f"\nsaved {OUTP}   (cost of the curvature constraint: {t1-base:+.1f})")
