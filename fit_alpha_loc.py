"""Who resists a common alpha'?  Release one block at a time."""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]; BX={"H_T^u":2,"H_T^d":6,"Ebar_T^u":10,"Ebar_T^d":12,"T00":21}
def slopes32(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
HA=dict(Q2=1.75,xB=0.36,mt=0.27)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else float('nan')
def run(tag, freeblock):
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.]); HI=np.array(list(F.HI)+[12.,5.,8.,5.,8.])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0                      # b >= 0 hard, it is free
    tiedbx={BX[k] for k in BX if k!=freeblock and BX[k]!=2}
    tied={7,27}|tiedbx
    free=[i for i in range(32) if i not in tied|{23,28,29,30,31}]
    def expand(x):
        p=np.zeros(32); p[free]=x; p[23]=0.0; p[7]=p[3]; p[27]=p[11]
        for i in tiedbx: p[i]=p[2]
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.05,0.35,1.00):
            s=slopes32(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in ("fitpar_both.npy","fitpar_alpha1.npy","fitpar_tieQx_data.npy"):
        q=np.zeros(32); z=np.load(src); q[:len(z)]=z; q[7]=q[3]; q[27]=q[11]
        for i in tiedbx: q[i]=q[2]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    aps=" ".join(f"{-p[BX[k]]:5.2f}" for k in BL)
    bs=" ".join(f"{p[i]:5.2f}" for i in BSLOT)
    print(f"{tag:>12} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f} | alpha': {aps} | b: {bs} | n/p {nratio(p):5.2f}",flush=True)
    np.save(f"fitpar_al_{tag.replace(' ','_')}.npy",p); return tot
print("all with b >= 0 hard and nQ_d = nQ_u.  alpha' common except for the released block.")
print("order: H_T^u  H_T^d  Ebar_T^u  Ebar_T^d  T00\n")
print(f"{'released':>12} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("none",None)
for k in ("H_T^d","Ebar_T^d","T00","H_T^u"): run(k,k)
