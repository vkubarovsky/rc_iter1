"""Final production fit: b2 removed entirely, 26 parameters.

Constraints kept (VPK's call): every block falls, b + b' ln xB >= 0 on
xB in [0.1, 0.6], and slope(H_T^d) >= slope(H_T^u).  With b2 = 0 the Ebar_T
slope is t-independent, so the large-|t| turnover cannot happen by construction
and no curvature constraint is needed.
"""
import math, os, numpy as np
from scipy.optimize import least_squares
import fit_slope as F
from reparam import slopes
OUTP=os.environ.get("OUTP","fitpar_amp2026s.npy")
W=30.0; XS=(0.10,0.60)
NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")
FREE=[i for i in range(27) if i!=12]
def expand(x):
    p=np.zeros(27); p[FREE]=x; p[12]=0.0; return p
def pen(p):
    r=[]
    for x in XS:
        s=slopes(p,x)
        for n in NAMES: r.append(W*min(0.0,s[n]))
        r.append(W*min(0.0,s["H_T^d"]-s["H_T^u"]))
    return r
def resid(x):
    p=expand(x); out=F.blocks(p); r=[]
    for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
    r.append((p[13]-0.54)/0.15); r+=pen(p)
    return np.array(r)
seeds=[]
for f in ("fitpar_nob2.npy","fitpar_prod.npy","fitpar_slope_order.npy","fitpar_slope.npy","fitpar_amp2026_lx.npy"):
    q=np.array(np.load(f)); q[12]=0.0; seeds.append((f,q))
for tag,mod in (("ET_N+20%",(8,1.2)),("ET_N-20%",(8,0.8)),("HTd flat",(5,0.0))):
    q=np.array(np.load("fitpar_nob2.npy"))
    q[mod[0]]=q[mod[0]]*mod[1] if mod[0]==8 else 1.5
    seeds.append((tag,q))
best=None
for tag,q in seeds:
    try:
        r=least_squares(resid,np.clip(q[FREE],F.LO[FREE],F.HI[FREE]),bounds=(F.LO[FREE],F.HI[FREE]),
                        x_scale='jac',xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
    except Exception as e:
        print(f"  seed {tag:22s} failed {e}"); continue
    cc=float(np.sum(r.fun**2)); print(f"  seed {tag:22s} chi2+pen {cc:9.2f}  nfev {r.nfev}")
    if best is None or cc<best[0]: best=(cc,r.x)
p=expand(best[1]); out=F.blocks(p)
tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
print(f"\n=== amp2026s (b2 removed) === {n} pts, 26 par: chi2/ndf = {tot:.1f}/{n-26} = {tot/(n-26):.4f}"
      f"   penalty {sum(x*x for x in pen(p)):.3f}  R_ET prior pull {(p[13]-0.54)/0.15:+.1f}")
for k,v in out.items(): print(f"   {k:>8}: {sum(x*x for x in v):7.1f}/{len(v)}")
print(f"\n   {'block':>10} {'b':>8} {'b_prime':>8} |" + "".join(f"{x:8.2f}" for x in (0.10,0.15,0.25,0.40,0.60)))
for name in NAMES:
    s={x:slopes(p,x)[name] for x in (0.10,0.15,0.25,0.40,0.60)}
    from reparam import BLOCKS
    ib,ibp=BLOCKS.get(name,(9,10))
    b=p[ib]+(p[14] if name=="Ebar_T^d" else 0.0); bp=p[ibp] if ibp<len(p) else 0.0
    print(f"   {name:>10} {b:8.3f} {bp:8.3f} |"+"".join(f"{s[x]:8.2f}" for x in (0.10,0.15,0.25,0.40,0.60)))
s=F.amp.structure(p,"pi0p",-0.4,0.25,1.94)
print(f"\n   R_HT(t=0) {p[4]/p[0]:+.3f}   R_ET {p[13]:.3f}   db_ET {p[14]:.3f}   R_L {p[18]:+.3f}"
      f"   delta(t) {p[19]:.2f}{p[20]:+.2f}t   sigma_L/sigma_T {s['L']/s['T']:.4f}")
np.save(OUTP,p); print(f"saved {OUTP}")
