"""Production fit with nQ_d = nQ_u imposed in BOTH sectors.

Why it is imposed and not fitted: at leading order the Q2 dependence factorises
out of the flavour structure -- the GPD carries no Q2 dependence (in the GK code
the log(Q2/Q0) is literally dead), and the hard subprocess is flavour-blind up to
the quark charge, a constant.  Transversity moreover does not mix with gluons, so
its evolution is non-singlet and flavour-blind as well.  The residual difference
is second order: at fixed xB the convolution samples a moving skewness (xi falls
25% over Q2 = 1.5-4.4) and the u,d x-shapes differ, which moves the d/u ratio by
4%, i.e. delta(nQ) ~ 0.04.  Our free fits were producing nQ_d - nQ_u of -7 to +3.

Everything else: no R_ET prior, b2 = 0, independent Ebar_T^d block, phi_CE = 0
(gauge), only constraint slope >= 0 on xB in [0.05, 1.0].
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
def slopes32(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.]); HI=np.array(list(F.HI)+[12.,5.,8.,5.,8.])
LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
HA=dict(Q2=1.75,xB=0.36,mt=0.27,val=0.28,err=0.07)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else float('nan')
FREE=[i for i in range(32) if i not in (7,23,27,28,29,30,31)]   # nQ_d tied, phi_CE gauge, no F
def expand(x):
    p=np.zeros(32); p[FREE]=x; p[23]=0.0
    p[7]=p[3]; p[27]=p[11]                  # nQ_d = nQ_u, both sectors
    return p
def run(tag, halla):
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        if halla: r.append((nratio(p)-HA["val"])/HA["err"])
        for xx in (0.05,0.35,1.00):     # positivity now over the WHOLE xB range
            s=slopes32(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in ("fitpar_lad_nQ.npy","fitpar_lad_free.npy","fitpar_night_noF.npy","fitpar_amp2026s.npy"):
        q=np.zeros(32); z=np.load(src); q[:len(z)]=z; q[7]=q[3]; q[27]=q[11]; q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[FREE],bounds=(LO[FREE],HI[FREE]),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2); s=slopes32(p,0.25)
    print(f"\n=== {tag} === {n} pts, {len(FREE)} par: chi2/ndf = {tot:.1f}/{n-len(FREE)} = {tot/(n-len(FREE)):.4f}")
    for k,v in out.items(): print(f"   {k:>8}: {sum(x*x for x in v):7.1f}/{len(v)}")
    print(f"   nQ: H_T {p[3]:+.3f} (u=d)   Ebar_T {p[11]:+.3f} (u=d)   T00 {p[17]:+.3f}")
    print(f"   at xB=0.25, Q2=2.2, -t=0.3:  H_T d/u {HT[1]/HT[0]:+.3f}   Ebar_T d/u {ET[1]/ET[0]:+.3f}   T00 d/u {LL[1]/LL[0]:+.3f}")
    print(f"   slopes at xB=0.25: " + "  ".join(f"{k} {v:.2f}" for k,v in s.items()))
    st=amp.structure(p,"pi0p",-0.4,0.25,1.94)
    print(f"   sigma_TT(n)/sigma_TT(p) = {nratio(p):.2f} (Hall-A 0.28+-0.07)   sigma_L/sigma_T = {st['L']/st['T']:.4f}")
    np.save(f"fitpar_{tag}.npy",p); return p
run("tieQx_data",False)
run("tieQx_halla",True)
