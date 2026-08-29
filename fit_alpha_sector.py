"""alpha' common between u and d WITHIN a sector, but H_T and Ebar_T free to differ.

Physics: within one sector the same trajectory is exchanged for both flavours, so
alpha' cannot depend on the quark.  Between sectors it may differ - H_T and
Ebar_T are different GPDs.  (GK uses 0.45 for everything, which is the stronger
assumption we already priced at 92 chi2.)

    ref       alpha' free in every block               (25 par)
    sector    alpha'_HT shared, alpha'_ET shared, T00 free   (23 par)
    sec+T00   T00 also tied to the Ebar_T value         (22 par)
    all       one alpha' everywhere                     (21 par, for reference)
All with b >= 0 hard, nQ_d = nQ_u, slope >= 0 on xB in [0.05, 1].
"""
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
def run(tag, ties):          # ties: dict target_slot -> source_slot
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.]); HI=np.array(list(F.HI)+[12.,5.,8.,5.,8.])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0
    tied=set(ties)|{7,27}
    free=[i for i in range(32) if i not in tied|{23,28,29,30,31}]
    def expand(x):
        p=np.zeros(32); p[free]=x; p[23]=0.0; p[7]=p[3]; p[27]=p[11]
        for t,s in ties.items(): p[t]=p[s]
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.05,0.35,1.00):
            s=slopes32(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in ("fitpar_production.npy","fitpar_al_Ebar_T_d.npy","fitpar_alpha1.npy","fitpar_tieQx_data.npy"):
        try: z=np.load(src)
        except Exception: continue
        q=np.zeros(32); q[:len(z)]=z; q[7]=q[3]; q[27]=q[11]
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    aps=" ".join(f"{-p[BX[k]]:5.2f}" for k in BL); bs=" ".join(f"{p[i]:5.2f}" for i in BSLOT)
    print(f"{tag:>8} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f} | alpha': {aps} | b: {bs} | "
          f"ETd/u {ET[1]/ET[0]:+5.2f} HTd/u {HT[1]/HT[0]:+5.2f} n/p {nratio(p):5.2f} "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f} eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
    np.save(f"fitpar_as_{tag}.npy",p); return tot
print("order of blocks: H_T^u  H_T^d  Ebar_T^u  Ebar_T^d  T00       (GK: alpha' = 0.45 for all)\n")
print(f"{'variant':>8} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("ref",{})
run("sector",{6:2, 12:10})
run("sec+T00",{6:2, 12:10, 21:10})
run("all",{6:2, 10:2, 12:2, 21:2})
