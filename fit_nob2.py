"""Is the b2 t^2 curvature of Ebar_T^u still needed?

It was introduced in the P2 GFF fit (~/pi0_eta_joint_fit) to give pi0 sigma_TT
its t-curvature without a steep d-component.  The amplitude fit has since gained
the BSA, CLAS12 and eg1-dvcs sectors.  Here b2 is frozen at 0 (26 parameters)
under the SAME production constraints, and the fit is compared to fitpar_prod.
"""
import math, numpy as np
from scipy.optimize import least_squares
import fit_slope as F
from reparam import slopes
W=30.0; XS=(0.10,0.60); TMAX=2.5
NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")
def pen(p):
    r=[]
    for x in XS:
        s=slopes(p,x)
        for n in NAMES: r.append(W*min(0.0,s[n]))
        r.append(W*min(0.0,s["H_T^d"]-s["H_T^u"]))
        for n in ("Ebar_T^u","Ebar_T^d"): r.append(W*min(0.0,s[n]-2*p[12]*TMAX))
    return r
FREE=[i for i in range(27) if i!=12]
def expand(x):
    p=np.zeros(27); p[FREE]=x; p[12]=0.0; return p
def resid(x):
    p=expand(x); out=F.blocks(p); r=[]
    for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
    r.append((p[13]-0.54)/0.15); r+=pen(p)
    return np.array(r)
P=np.load("fitpar_prod.npy")
def show(tag,p,npar):
    out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    tt_lo=tt_hi=0.
    for d in F.PI0:
        s=F.amp.structure(p,"pi0p",d["t"],d["xB"],d["Q2"])
        c=((d["TT"]-s["TT"])/d["dTT"])**2
        (tt_lo if -d["t"]<0.8 else tt_hi).__class__  # noqa
        if -d["t"]<0.8: tt_lo+=c
        else: tt_hi+=c
    print(f"{tag:>16}  {npar} par  chi2 {tot:7.1f}  chi2/ndf {tot/(n-npar):.4f}   "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f}  eta {sum(x*x for x in out['eta']):6.1f}   "
          f"pi0 sigma_TT: -t<0.8 {tt_lo:5.1f}  -t>0.8 {tt_hi:5.1f}   b2 {p[12]:+.3f}")
    return tot
c0=show("b2 free (prod)",P,27)
best=None
for s0 in (P[FREE], np.where(np.arange(26)==8,P[8]*1.1,P[FREE])):
    r=least_squares(resid,np.clip(s0,F.LO[FREE],F.HI[FREE]),bounds=(F.LO[FREE],F.HI[FREE]),
                    x_scale='jac',xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
    cc=float(np.sum(r.fun**2))
    if best is None or cc<best[0]: best=(cc,r.x)
p=expand(best[1]); c1=show("b2 = 0",p,26)
print(f"\ncost of dropping b2: dchi2 = {c1-c0:+.1f} for one parameter")
print(f"Ebar_T^u slope now: " + "  ".join(f"xB={x}: {slopes(p,x)['Ebar_T^u']:.2f}" for x in (0.10,0.25,0.60))
      + "   (no turnover at any t)")
np.save("fitpar_nob2.npy",p)
