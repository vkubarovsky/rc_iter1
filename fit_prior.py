"""The R_ET prior is applied at the wrong level -- VPK, 2026-08-28.

R_ET is the ratio of the d and u CONVOLUTIONS (GPD x hard kernel), not of the
first moments of the GPDs.  kappa_T^d/kappa_T^u = 0.54 +- 0.15 constrains the
moments; the two ratios coincide only if u and d have the same x shape, and they
do not.  So this fit compares:

    prior   R_ET = 0.54 +- 0.15          (what we have been doing - wrong level)
    none    no constraint on R_ET        (why the prior was introduced: it runs away)
    halla   sigma_TT(n)/sigma_TT(p) = 0.28 +- 0.07 at Q2=1.75, xB=0.36, -t=0.27
            (PRL 117 262001 / 118 222002) - an OBSERVABLE-level constraint that
            needs only isospin, no GPD-level assumption at all
    both    no R_ET prior + the Hall-A datum

All variants: b2 = 0, slopes constrained, nQ_d tied to nQ_u (removes the pinned
-6 boundary artefact at a cost of 9 chi2).  Point-evaluated: bin averaging was
shown to move the model by 0.6% median, so it is not worth its cost here.
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
from reparam import slopes
W=30.0
NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")
HA=dict(Q2=1.75,xB=0.36,mt=0.27,val=0.28,err=0.07)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"])
    sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else 9.9
def pen(p):
    r=[]
    for x in (0.10,0.60):
        s=slopes(p,x)
        for n in NAMES: r.append(W*min(0.0,s[n]))
        r.append(W*min(0.0,s["H_T^d"]-s["H_T^u"]))
    return r
FREE=[i for i in range(27) if i not in (7,12)]
def expand(x):
    p=np.zeros(27); p[FREE]=x; p[12]=0.0; p[7]=p[3]; return p
def run(tag,use_prior,use_halla):
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        if use_prior: r.append((p[13]-0.54)/0.15)
        if use_halla: r.append((nratio(p)-HA["val"])/HA["err"])
        r+=pen(p)
        return np.array(r)
    best=None
    for seed,mod in ((np.load("fitpar_tied.npy"),None),(np.load("fitpar_amp2026s.npy"),None),
                     (np.load("fitpar_tied.npy"),(13,0.4)),(np.load("fitpar_tied.npy"),(13,2.0))):
        q=np.array(seed); q[12]=0.0; q[7]=q[3]
        if mod: q[mod[0]]=mod[1]
        try:
            r=least_squares(resid,np.clip(q[FREE],F.LO[FREE],F.HI[FREE]),bounds=(F.LO[FREE],F.HI[FREE]),
                            x_scale='jac',xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    print(f"{tag:>7} {tot:9.1f} {tot/(n-len(FREE)):8.4f}  R_ET {p[13]:6.3f}  R_HT {p[4]/p[0]:+7.3f}  "
          f"db {p[14]:5.2f}  sigTT(n)/sigTT(p) {nratio(p):5.2f}  "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f} "
          f"bsa {sum(sum(x*x for x in out[k]) for k in ('bsa_pi0','bsa_eta','bsa_c12')):6.1f}",flush=True)
    np.save(f"fitpar_prior_{tag}.npy",p); return p
print("data chi2 only (constraint terms excluded);  Hall-A measures 0.28 +- 0.07\n")
print(f"{'variant':>7} {'chi2':>9} {'chi2/ndf':>8}")
run("prior",True,False)
run("none",False,False)
run("halla",False,True)
run("both",True,True)
