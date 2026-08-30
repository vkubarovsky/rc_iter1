"""Normalisation nuisance parameters: lambda (common) and rho (eta vs pi0).

The absolute normalisation of both data sets comes from ONE elastic-scattering
measurement on the SAME e1-dvcs run: delta_Norm = 0.89 quoted with 6.0% in the
pi0 paper (PRC 90 025205 Table IV, marked "overall") and 0.87 with 4.1% in the
eta paper (PRC 95 035202 Table IV).  It is one physical number estimated twice,
so it cancels in the eta/pi0 ratio.  The published per-point systematics exclude
it (eta paper Fig. 12 caption; and our points scatter from 3.8% to 24%, below the
4.1% floor a normalisation term would impose), so nothing has to be removed from
the errors -- the nuisance parameters are simply added.

    sigma_model(pi0) -> lambda * model          prior lambda = 1.00 +- 0.05
    sigma_model(eta) -> lambda * rho * model    prior rho    = 1.00 +- 0.03

BSA and the eg1 moments are ratios of cross sections, so the scale cancels there
and lambda is not applied to them.

The question this answers: how much of the steep Ebar_T^d slope is the fit
buying by sliding eta against pi0 in normalisation?
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]
def slopes32(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
HA=dict(Q2=1.75,xB=0.36,mt=0.27)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else float('nan')
def blocks_scaled(p, lam, rho):
    """cross-section residuals with the two scales applied; BSA/eg1 untouched."""
    out={"pi0":[], "eta":[]}
    for rows,ch,tag,sc in ((F.PI0,"pi0p","pi0",lam),(F.ETA,"etap","eta",lam*rho)):
        for d in rows:
            s=amp.structure(p,ch,d["t"],d["xB"],d["Q2"])
            if s is None: out[tag]+=[10.,10.,10.]; continue
            e=amp.epsilon(d["xB"],d["Q2"],F.E_XS)
            out[tag]+=[(d["U"]-sc*(s["T"]+e*s["L"]))/d["dU"],
                       (d["TT"]-sc*s["TT"])/d["dTT"],
                       (d["LT"]-sc*s["LT"])/d["dLT"]]
    o2=F.blocks(p)
    for k in ("bsa_pi0","bsa_eta","bsa_c12","eg1"): out[k]=o2[k]
    return out
def run(tag, seed, sector, use_norm, sig_lam=0.05, sig_rho=0.03):
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.,0.5,0.5]); HI=np.array(list(F.HI)+[12.,5.,8.,5.,8.,1.5,1.5])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0
    ties={7:3, 27:11}
    if sector: ties.update({6:2, 12:10})
    fixed=set(ties)|{23,28,29,30,31}
    if not use_norm: fixed|={32,33}
    free=[i for i in range(34) if i not in fixed]
    def expand(x):
        p=np.zeros(34); p[free]=x; p[23]=0.0
        if not use_norm: p[32]=p[33]=1.0
        for t,s in ties.items(): p[t]=p[s]
        return p
    def resid(x):
        p=expand(x); lam,rho=p[32],p[33]
        out=blocks_scaled(p[:32],lam,rho); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        if use_norm: r+= [(lam-1.0)/sig_lam, (rho-1.0)/sig_rho]
        for xx in (0.05,0.35,1.00):
            s=slopes32(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in (seed,"fitpar_production.npy","fitpar_as_sector.npy"):
        q=np.ones(34); z=np.load(src); q[:len(z)]=z; q[32]=q[33]=1.0
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); lam,rho=p[32],p[33]
    out=blocks_scaled(p[:32],lam,rho)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    s25=slopes32(p,0.25); HT,ET,LL=amp._flavour(p[:32],-0.3,0.25,2.2)
    print(f"{tag:>16} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f} | lam {lam:6.3f} rho {rho:6.3f} | "
          f"slope ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f}  ETd/u {ET[1]/ET[0]:+5.2f}  n/p {nratio(p[:32]):5.2f} | "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f}",flush=True)
    np.save(f"fitpar_nrm_{tag.replace(' ','_')}.npy",p); return tot
print("cross-section chi2 only in the total; BSA/eg1 are ratios and are not scaled\n")
print(f"{'variant':>16} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("prod, no norm","fitpar_production.npy",False,False)
run("prod + norm","fitpar_production.npy",False,True)
run("sector, no norm","fitpar_as_sector.npy",True,False)
run("sector + norm","fitpar_as_sector.npy",True,True)
run("sector,rho free","fitpar_as_sector.npy",True,True,0.05,0.15)
