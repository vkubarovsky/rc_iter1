"""Ladder of u/d ties, both sectors.  Goal: find out what the data actually
demand of the d block, and whether d can be reduced to R x u.

Base: no R_ET prior, no d>=u ordering, b2 = 0, no F(xB), phi_CE = 0 (gauge),
independent Ebar_T^d block, only constraint is slope >= 0 on xB in [0.1, 0.6].

    free      everything free
    nQ        nQ_d = nQ_u in both sectors
    nQ+bx     also b'_d = b'_u
    nQ+bx+b   also b_d = b_u   -> d = R x u exactly: one normalisation per sector
    HT tied   d = R x u for H_T only
    ET tied   d = R x u for Ebar_T only
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
def slopes32(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.]); HI=np.array(list(F.HI)+[12.,5.,8.,5.,8.])
LO[3]=LO[7]=LO[11]=LO[17]=-12.; HI[3]=HI[7]=HI[11]=HI[17]=12.
LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
HA=dict(Q2=1.75,xB=0.36,mt=0.27)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else float('nan')
# tie maps: which d parameter copies which u parameter
TIE={"nQ_H":(7,3),"nQ_E":(27,11),"bx_H":(6,2),"bx_E":(12,10),"b_H":(5,1),"b_E":(14,9)}
def run(tag, ties):
    fixed={TIE[t][0] for t in ties}|{12+16,23}  # 28..31 (no F) handled below, 23 = phi_CE gauge
    free=[i for i in range(32) if i not in fixed and i not in (28,29,30,31)]
    def expand(x):
        p=np.zeros(32); p[free]=x; p[23]=0.0
        for t in ties: p[TIE[t][0]]=p[TIE[t][1]]
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.10,0.60):
            s=slopes32(p,xx)
            for n in NAMES: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in ("fitpar_night_noF.npy","fitpar_etd_free.npy","fitpar_amp2026s.npy"):
        q=np.zeros(32); z=np.load(src); q[:len(z)]=z
        for t in ties: q[TIE[t][0]]=q[TIE[t][1]]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=60000)
        except Exception as e: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2); s25=slopes32(p,0.25)
    print(f"{tag:>10} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f}  "
          f"HTd/HTu {HT[1]/HT[0]:+6.2f}  ETd/ETu {ET[1]/ET[0]:+6.2f}  "
          f"slope HT u/d {s25['H_T^u']:4.2f}/{s25['H_T^d']:5.2f}  ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f}  "
          f"nQ HT u/d {p[3]:+5.2f}/{p[7]:+5.2f} ET u/d {p[11]:+5.2f}/{p[27]:+5.2f}  n/p {nratio(p):5.2f}  "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f} eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
    np.save(f"fitpar_lad_{tag.replace('+','_')}.npy",p); return tot
print("ratios quoted at xB=0.25, Q2=2.2, -t=0.3.  Hall-A n/p = 0.28 +- 0.07\n")
print(f"{'variant':>10} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
base=run("free",[])
run("nQ",["nQ_H","nQ_E"])
run("nQ+bx",["nQ_H","nQ_E","bx_H","bx_E"])
run("nQ+bx+b",["nQ_H","nQ_E","bx_H","bx_E","b_H","b_E"])
run("HT tied",["nQ_H","bx_H","b_H"])
run("ET tied",["nQ_E","bx_E","b_E"])
