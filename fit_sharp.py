"""b >= 0 as a hard bound, and one common alpha' for all blocks.

b is the t-slope at xB = 1, so b >= 0 forbids a form factor that grows with |t|
anywhere, by construction rather than by penalty.  alpha' = -b' is the slope of
the exchanged Regge trajectory: it belongs to the trajectory, not to the GPD, so
one value should serve every block (GK uses 0.45 throughout, both flavours, both
sectors, both versions).

    ref        current final model (b free, alpha' per block)
    b0         b >= 0 hard
    alpha1     one common alpha'
    both       b >= 0 and one common alpha'
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]; BXSLOT=[2,6,10,12,21]
def slopes32(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
HA=dict(Q2=1.75,xB=0.36,mt=0.27)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else float('nan')
def run(tag, b0, alpha1):
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.]); HI=np.array(list(F.HI)+[12.,5.,8.,5.,8.])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    if b0:
        for i in BSLOT: LO[i]=0.0
    tied={7,27}|({6,10,12,21} if alpha1 else set())
    free=[i for i in range(32) if i not in tied|{23,28,29,30,31}]
    def expand(x):
        p=np.zeros(32); p[free]=x; p[23]=0.0; p[7]=p[3]; p[27]=p[11]
        if alpha1:
            for i in (6,10,12,21): p[i]=p[2]
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.05,0.35,1.00):
            s=slopes32(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in ("fitpar_tieQx_data.npy","fitpar_tieQ_data.npy","fitpar_lad_nQ.npy"):
        q=np.zeros(32); z=np.load(src); q[:len(z)]=z; q[7]=q[3]; q[27]=q[11]
        if alpha1:
            for i in (6,10,12,21): q[i]=q[2]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception as e: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2); s=slopes32(p,0.25)
    bs=" ".join(f"{p[i]:5.2f}" for i in BSLOT); aps=" ".join(f"{-p[i]:5.2f}" for i in BXSLOT)
    print(f"{tag:>7} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f} | b: {bs} | alpha': {aps} | "
          f"HTd/u {HT[1]/HT[0]:+5.2f} ETd/u {ET[1]/ET[0]:+5.2f} n/p {nratio(p):5.2f} "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f} eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
    np.save(f"fitpar_{tag}.npy",p); return tot
print("blocks in order: H_T^u  H_T^d  Ebar_T^u  Ebar_T^d  T00        (GK: b = 0.30/0.30/0.77/0.50, alpha' = 0.45 everywhere)\n")
print(f"{'variant':>7} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("ref",False,False)
run("b0",True,False)
run("alpha1",False,True)
run("both",True,True)
