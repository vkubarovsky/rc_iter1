"""Night run: does an xB shape F(xB) = xB^alpha (1-xB)^n help Ebar_T?

VPK tried this years ago and concluded it works better without.  Here it is
tested inside the amplitude fit, on top of the independent Ebar_T^d block and
with no R_ET prior.  The only constraint is that every t-slope be >= 0 on
xB in [0.1, 0.6].

Ratios are quoted at -t = 0.3 (inside the data), never at t = 0: t = 0 lies
below |tmin| and is not a physical point.
"""
import math, os, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
def slopes32(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
LO=np.array(list(F.LO)+[-12.0,-5.0,-5.0,-5.0,-5.0]); HI=np.array(list(F.HI)+[12.0,5.0,8.0,5.0,8.0])
LO[3]=LO[7]=LO[11]=LO[17]=-12.0; HI[3]=HI[7]=HI[11]=HI[17]=12.0
LO[12]=-8.0; HI[12]=8.0; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.0; HI[14]=12.0
def seed():
    q=np.load("fitpar_etd_free.npy"); p=np.zeros(32); p[:28]=q; return p   # alphas, n = 0
HA=dict(Q2=1.75,xB=0.36,mt=0.27)
def report(tag,p,free,out):
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    s25=slopes32(p,0.25)
    print(f"{tag:>12} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f}  "
          f"ETd/ETu(-t=0.3) {ET[1]/ET[0]:+7.3f}  slope u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:4.2f}  "
          f"n/p {sn['TT']/sp['TT']:5.2f}  a_u {p[28]:+5.2f} n_u {p[29]:+5.2f} a_d {p[30]:+5.2f} n_d {p[31]:+5.2f}  "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f} "
          f"eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
def run(tag, free, tie=None):
    def expand(x):
        p=np.zeros(32); p[free]=x
        if tie: tie(p)
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.10,0.60):
            s=slopes32(p,xx)
            for n in s: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for jitter in (0.0,0.5,-0.5):
        q=seed()
        if 29 in free: q[29]=jitter
        if 31 in free: q[31]=jitter
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=60000)
        except Exception as e: print("  seed fail",e); continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p); report(tag,p,free,out); np.save(f"fitpar_night_{tag}.npy",p)
    return p
ALL=list(range(32))
print("all fits: no R_ET prior, independent Ebar_T^d, slope >= 0 only.  Hall-A n/p = 0.28 +- 0.07\n")
print(f"{'variant':>12} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("noF",         [i for i in ALL if i not in (28,29,30,31)])
run("F_common",    [i for i in ALL if i not in (30,31)], tie=lambda p:(p.__setitem__(30,p[28]),p.__setitem__(31,p[29])))
run("F_split",     ALL)
run("F_1mx_common",[i for i in ALL if i not in (28,30,31)], tie=lambda p:p.__setitem__(31,p[29]))
run("F_xa_common", [i for i in ALL if i not in (29,30,31)], tie=lambda p:p.__setitem__(30,p[28]))
run("F_1mx_split", [i for i in ALL if i not in (28,30)])
