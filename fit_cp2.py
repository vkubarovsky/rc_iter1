"""Phases computed, not fitted, and what that does to the L/T separation.

    B          model B, everything real at flavour level      (23 par)
    phi_fit    relative u-d phase free                        (24 par) -> 1.28 rad
    phi_calc   relative u-d phase COMPUTED from the convolution (23 par, no new freedom)
    delta_calc phi computed AND delta(t) fixed to the computed L-T phase (21 par)
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]
def sl(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
HA=dict(Q2=1.75,xB=0.36,mt=0.27)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else float('nan')
# computed L-T phase: arg<L> - arg<T> from the same convolutions, fitted linearly in xi
XI=np.array([0.07,0.14,0.22,0.33]); DLT=np.array([0.347,0.319,0.285,0.235])
c=np.polyfit(XI,DLT,1)
def delta_calc(xB,Q2): return float(np.polyval(c,amp.ksi(xB,Q2)))
def run(tag, mode):
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.,-math.pi,-math.pi,-math.pi,-math.pi,-math.pi])
    HI=np.array(list(F.HI)+[ 12., 5., 8., 5., 8., math.pi, math.pi, math.pi, math.pi, math.pi])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0
    ties={7:3, 27:11, 6:2, 12:10}
    frozen={23,28,29,30,31,32,33,34,36}
    if mode in ("B",): frozen|={35}
    if mode=="delta_calc": frozen|={19,20}
    free=[i for i in range(37) if i not in set(ties)|frozen]
    def expand(x):
        p=np.zeros(37); p[free]=x; p[23]=0.0
        for t,s in ties.items(): p[t]=p[s]
        if mode in ("phi_calc","delta_calc"): p[35]=-99.0
        return p
    def resid(x):
        p=expand(x)
        if mode=="delta_calc":
            r=[]
            for rows,ch,tag2 in ((F.PI0,"pi0p","pi0"),(F.ETA,"etap","eta")):
                for d in rows:
                    q=np.array(p); q[19]=delta_calc(d["xB"],d["Q2"]); q[20]=0.0
                    s=amp.structure(q,ch,d["t"],d["xB"],d["Q2"])
                    if s is None: r+=[10.,10.,10.]; continue
                    e=amp.epsilon(d["xB"],d["Q2"],F.E_XS)
                    r+=[(d["U"]-(s["T"]+e*s["L"]))/d["dU"],(d["TT"]-s["TT"])/d["dTT"],(d["LT"]-s["LT"])/d["dLT"]]
            for pts,ch,E in ((F.BSA_PTS,"pi0p",F.E_BSA),(F.EBSA,"etap",F.E_BSA),(F.C12_PTS,"pi0p",F.E_C12)):
                for (xB,Q2,t,a,da) in pts:
                    q=np.array(p); q[19]=delta_calc(xB,Q2); q[20]=0.0
                    v=amp.bsa_sinphi(q,ch,-t,xB,Q2,E); r.append(((a-v)/da) if v is not None else 5.0)
            for (key,Q2,xB,mt,A,dA) in F.EG1_PTS:
                q=np.array(p); q[19]=delta_calc(xB,Q2); q[20]=0.0
                s=amp.structure(q,"pi0p",-mt,xB,Q2)
                if s is None: r.append(5.0); continue
                e=amp.epsilon(xB,Q2,F.E_EG1); s0=s["T"]+e*s["L"]
                v=dict(AULsin=math.sqrt(2*e*(1+e))*s["LTUL"]/s0, AULsin2=e*s["TTUL"]/s0,
                       ALLc=math.sqrt(1-e*e)*s["Tp"]/s0, ALLcos=math.sqrt(2*e*(1-e))*s["LTpLL"]/s0)[key]
                r.append((A-v)/dA)
        else:
            out=F.blocks(p); r=[]
            for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.05,0.35,1.00):
            s=sl(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for scale in (1.0,0.5,2.0):
        q=np.zeros(37); z=np.load("fitpar_production.npy"); q[:len(z)]=z; q[15]*=scale
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
        except Exception as e: continue
        cc=float(np.sum(r.fun**2))
        if best is None or cc<best[0]: best=(cc,r.x)
    p=expand(best[1])
    if mode=="delta_calc":
        pp=np.array(p); pp[19]=delta_calc(0.25,2.2); pp[20]=0.0
    else: pp=p
    out=F.blocks(pp); tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    s25=sl(p,0.25); st=amp.structure(pp,"pi0p",-0.4,0.25,1.94)
    print(f"{tag:>11} {len(free):4d} {tot:8.1f} | phi_du {'computed' if mode in ('phi_calc','delta_calc') else f'{p[35]:+.2f}':>9} "
          f"delta {'computed' if mode=='delta_calc' else f'{p[19]:.2f}':>8} | slope ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f} "
          f"n/p {nratio(pp):5.2f} | sigma_L/sigma_T {st['L']/st['T']:.4f} | pi0 {sum(x*x for x in out['pi0']):6.1f} "
          f"eta {sum(x*x for x in out['eta']):6.1f} bsa {sum(sum(x*x for x in out[k]) for k in ('bsa_pi0','bsa_eta','bsa_c12')):6.1f} "
          f"eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
    np.save(f"fitpar_cp_{tag}.npy",p); return tot
print(f"computed L-T phase: {delta_calc(0.25,2.2):.3f} rad at xB=0.25, Q2=2.2;  the fit had 2.00\n")
print(f"{'variant':>11} {'par':>4} {'chi2':>8}")
run("phi_calc","phi_calc")
run("delta_calc","delta_calc")
