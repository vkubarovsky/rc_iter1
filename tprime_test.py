"""Is the b2 curvature a near-tmin effect?

sigma_TT = |T01p|^2 - |U01p|^2 with T01m = U01m cancelling, and |U01p|^2 carries
the kinematic factor -t'/(8 M^2): the model vanishes linearly in t' at tmin BY
CONSTRUCTION.  Half the pi0 cross-section points sit at -t' < 0.3.  So: does b2
survive when the near-tmin region is removed from the fit?
"""
import math, os, numpy as np
from scipy.optimize import least_squares
import fit_slope as F
from reparam import slopes
import amplitudes as amp
TPCUT=float(os.environ.get("TPCUT","0.3"))
def tp_of(ch,d):
    mM=amp.Meta if ch=="etap" else amp.Mpi0
    return -(d["t"]-amp.tmin(mM,d["Q2"],d["xB"]))
n0=(len(F.PI0),len(F.ETA))
F.PI0=[d for d in F.PI0 if tp_of("pi0p",d)>TPCUT]
F.ETA=[d for d in F.ETA if tp_of("etap",d)>TPCUT]
keep=lambda pts,ch: [q for q in pts if tp_of(ch,dict(t=-q[2],xB=q[0],Q2=q[1]))>TPCUT]
F.BSA_PTS=keep(F.BSA_PTS,"pi0p"); F.EBSA=keep(F.EBSA,"etap"); F.C12_PTS=keep(F.C12_PTS,"pi0p")
F.EG1_PTS=[q for q in F.EG1_PTS if tp_of("pi0p",dict(t=-q[3],xB=q[2],Q2=q[1]))>TPCUT]
print(f"-t' > {TPCUT}: pi0 {len(F.PI0)}/{n0[0]} bins, eta {len(F.ETA)}/{n0[1]} bins, "
      f"BSA {len(F.BSA_PTS)}+{len(F.EBSA)}+{len(F.C12_PTS)}, eg1 {len(F.EG1_PTS)}")
W=30.0
def pen(p):
    r=[]
    for x in (0.10,0.60):
        s=slopes(p,x)
        for n in ("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00"): r.append(W*min(0.0,s[n]))
        r.append(W*min(0.0,s["H_T^d"]-s["H_T^u"]))
        for n in ("Ebar_T^u","Ebar_T^d"): r.append(W*min(0.0,s[n]-2*p[12]*2.5))
    return r
def fit(free,fix_b2):
    def expand(x):
        p=np.zeros(27); p[free]=x
        if fix_b2: p[12]=0.0
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        r.append((p[13]-0.54)/0.15); r+=pen(p)
        return np.array(r)
    P=np.load("fitpar_prod.npy")
    best=None
    for scale in (1.0,1.1):
        s0=np.array(P[free]); s0[0]*=scale
        r=least_squares(resid,np.clip(s0,F.LO[free],F.HI[free]),bounds=(F.LO[free],F.HI[free]),
                        x_scale='jac',xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        cc=float(np.sum(r.fun**2))
        if best is None or cc<best[0]: best=(cc,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    return tot,n,p,out
a,na,pa,oa=fit(list(range(27)),False)
b,nb,pb,ob=fit([i for i in range(27) if i!=12],True)
print(f"\n{'':10} {'chi2':>9} {'npts':>6} {'chi2/ndf':>9} {'pi0':>8} {'eta':>8}  b2")
print(f"{'b2 free':>10} {a:9.1f} {na:6d} {a/(na-27):9.4f} {sum(x*x for x in oa['pi0']):8.1f} {sum(x*x for x in oa['eta']):8.1f}  {pa[12]:+.3f}")
print(f"{'b2 = 0':>10} {b:9.1f} {nb:6d} {b/(nb-26):9.4f} {sum(x*x for x in ob['pi0']):8.1f} {sum(x*x for x in ob['eta']):8.1f}  {pb[12]:+.3f}")
print(f"\ncost of dropping b2 away from tmin: dchi2 = {b-a:+.1f}   (full data set: +71.0)")
