"""chi2 profile in db_ET: freeze db, refit everything else, see what R_ET does.

db free lands at 4.2 with R_ET 4.2 (a d-spike at t=0 that is dead by -t=0.5);
db = 0 lands at R_ET = 0.04 (no d at all) and costs 157 chi2.  The profile shows
whether anything sane lives in between.
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
from reparam import slopes
W=30.0; NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")
LO=np.array(F.LO); HI=np.array(F.HI)
LO[3]=LO[7]=LO[11]=LO[17]=-12.0; HI[3]=HI[7]=HI[11]=HI[17]=12.0
LO[13]=0.0; HI[13]=20.0
FREE=[i for i in range(27) if i not in (12,14)]
HA=dict(Q2=1.75,xB=0.36,mt=0.27)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else float('nan')
print(f"{'db':>5} {'chi2':>9} {'d_chi2':>8} {'R_ET':>7} {'slope ET_d':>11} {'n/p':>6} "
      f"{'pi0':>7} {'eta':>7} {'bsa':>7} {'eg1':>7}   (Hall-A n/p = 0.28+-0.07)")
base=None
for db in (0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.22,5.0):
    def expand(x,db=db):
        p=np.zeros(27); p[FREE]=x; p[12]=0.0; p[14]=db; return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.10,0.60):
            s=slopes(p,xx)
            for n in NAMES: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for f in ("fitpar_plain.npy","fitpar_nodb.npy","fitpar_amp2026s.npy"):
        q=np.clip(np.array(np.load(f)),LO,HI); q[12]=0.0
        try: r=least_squares(resid,q[FREE],bounds=(LO[FREE],HI[FREE]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values())
    if base is None: base=tot
    bsa=sum(sum(x*x for x in out[k]) for k in ("bsa_pi0","bsa_eta","bsa_c12"))
    print(f"{db:5.2f} {tot:9.1f} {tot-base:+8.1f} {p[13]:7.3f} {slopes(p,0.25)['Ebar_T^d']:11.2f} "
          f"{nratio(p):6.2f} {sum(x*x for x in out['pi0']):7.1f} {sum(x*x for x in out['eta']):7.1f} "
          f"{bsa:7.1f} {sum(x*x for x in out['eg1']):7.1f}",flush=True)
    np.save(f"fitpar_db{db:.2f}.npy",p)
