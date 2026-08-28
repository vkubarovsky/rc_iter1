"""Price of the slope constraint, block by block.

Baseline for every variant: slope = b + b' ln xB >= 0 for ALL blocks on
xB in [0.1, 0.6] (the physical requirement - form factors must not grow with |t|).
On top of that:
  scan   floor on H_T^d alone: how hard do the data resist a FALLING d-quark H_T?
  order  slope(H_T^d) >= slope(H_T^u) at both endpoints - the ordering of VPK's
         global GPD fits (d steeper than u), which amp2026 inverts.
"""
import math, numpy as np
from scipy.optimize import least_squares
import fit_slope as F
from reparam import slopes

XLO,XHI,W=0.10,0.60,30.0
NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")

def pen(p, floor_d=0.0, order=False):
    r=[]
    for x in (XLO,XHI):
        s=slopes(p,x)
        for n in NAMES:
            fl = floor_d if n=="H_T^d" else 0.0
            r.append(W*min(0.0, s[n]-fl))
        if order: r.append(W*min(0.0, s["H_T^d"]-s["H_T^u"]))
    return r

def run(tag, floor_d=0.0, order=False, seeds=((None,None),)):
    def resid(p):
        out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        r.append((p[13]-0.54)/0.15)
        r+=pen(p,floor_d,order)
        return np.array(r)
    starts=[np.array(P0)]
    for bd,bpd in ((floor_d+1.0,0.0),(floor_d+2.0,-0.5),(3.0,-1.2)):
        q=np.array(P0); q[5],q[6]=bd,bpd; starts.append(q)
    best=None
    for sd in starts:
        try:
            r=least_squares(resid,np.clip(sd,F.LO,F.HI),bounds=(F.LO,F.HI),x_scale='jac',
                            xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        cc=float(np.sum(r.fun**2))
        if best is None or cc<best[0]: best=(cc,r.x)
    p=best[1]; out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    s=F.amp.structure(p,"pi0p",-0.4,0.25,1.94)
    su,sd_=slopes(p,0.25)["H_T^u"],slopes(p,0.25)["H_T^d"]
    print(f"{tag:>12} {tot:9.1f} {tot-BASE:+8.1f} {tot/(n-27):9.4f} "
          f"{sum(x*x for x in out['pi0']):8.1f} {sum(x*x for x in out['eta']):8.1f} "
          f"{p[4]/p[0]:+7.3f} {p[13]:7.3f} {(p[13]-0.54)/0.15:+6.1f} "
          f"{su:7.2f} {sd_:7.2f} {s['L']/s['T']:7.4f}  pen {sum(x*x for x in pen(p,floor_d,order)):.2f}")
    return tot,p

P0=np.load("fitpar_amp2026_lx.npy")
out0=F.blocks(P0); BASE=sum(sum(x*x for x in v) for v in out0.values())
print(f"unconstrained amp2026: chi2 = {BASE:.1f}\n")
print(f"{'variant':>12} {'chi2':>9} {'d_chi2':>8} {'chi2/ndf':>9} {'pi0':>8} {'eta':>8} "
      f"{'R_HT':>7} {'R_ET':>7} {'pull':>6} {'sl_u.25':>7} {'sl_d.25':>7} {'sL/sT':>7}")
for fd in (0.0,0.25,0.5,0.75,1.0,1.5,2.0,2.76):
    tot,p=run(f"d>={fd:.2f}",floor_d=fd)
    np.save(f"fitpar_slope_d{fd:.2f}.npy",p)
tot,p=run("d>=u",order=True); np.save("fitpar_slope_order.npy",p)
