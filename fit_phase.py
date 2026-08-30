"""Relative u-d phase, one per sector.

Until now every flavour GFF was real and the flavour sum was real; the phases in
the model (delta, phi_w, phi_CE, phi_nf) were attached AFTER the sum, i.e. one
per channel.  But each convolution is complex on its own -- Im = pi F(xi,xi,t),
Re = principal value -- and those two pieces have different t-slopes, so u and d
cannot share a phase.

This matters exactly where our problem is.  The eta channel is 2u - d, close to a
node; with real u and d the cancellation is exact and the fit can only keep eta
sigma_TT alive by making d fall steeply in t (our b_d = 5.55).  With a relative
phase the cancellation is incomplete by construction:
    |2u - d e^{i phi}|^2 >= (2|u| - |d|)^2
so eta sigma_TT survives without any steep slope.

Built on model B (sector-common alpha', 23 par).
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
def run(tag, phases):        # phases: subset of {34:'H_T', 35:'Ebar_T', 36:'T00'}
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.,-math.pi,-math.pi,-math.pi,-math.pi,-math.pi])
    HI=np.array(list(F.HI)+[ 12., 5., 8., 5., 8., math.pi, math.pi, math.pi, math.pi, math.pi])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0
    ties={7:3, 27:11, 6:2, 12:10}
    frozen={23,28,29,30,31}|({34,35,36}-set(phases))
    free=[i for i in range(37) if i not in set(ties)|frozen]
    def expand(x):
        p=np.zeros(37); p[free]=x; p[23]=0.0
        for t,s in ties.items(): p[t]=p[s]
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.05,0.35,1.00):
            s=sl(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    seeds=[]
    for ph in (0.0,0.6,-0.6,1.5,-1.5,2.5):
        q=np.zeros(37); z=np.load("fitpar_production.npy"); q[:len(z)]=z
        for i in phases: q[i]=ph
        for t,s in ties.items(): q[t]=q[s]
        seeds.append(np.clip(q,LO,HI))
    for q in seeds:
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    s25=sl(p,0.25)
    ph=" ".join(f"{p[i]:+5.2f}" for i in (34,35,36))
    print(f"{tag:>12} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f} | phases HT/ET/T00 {ph} | "
          f"slope ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f}  HT u/d {s25['H_T^u']:4.2f}/{s25['H_T^d']:4.2f}  "
          f"n/p {nratio(p):5.2f} | pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f} "
          f"bsa {sum(sum(x*x for x in out[k]) for k in ('bsa_pi0','bsa_eta','bsa_c12')):6.1f} eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
    np.save(f"fitpar_ph_{tag}.npy",p); return tot
print("model B for reference: 23 par, chi2 1133.6, slope ET u/d 1.10/5.94, n/p 1.09\n")
print(f"{'variant':>12} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("none",[])
run("ET",[35])
run("HT+ET",[34,35])
run("all three",[34,35,36])
