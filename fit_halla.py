"""Plain fit PLUS the Hall-A neutron ratio (observable level, isospin only).

Originally:  The only thing imposed is that
form factors fall -- slope b + b' ln xB >= 0 on xB in [0.1, 0.6].  b2 = 0.

Bounds widened where the previous fits were pinned: R_ET was hitting 3.0 and
nQ_d was hitting -6, and a bound that the minimum sits on is a prior in disguise.
Point-evaluated (bin averaging moves the model by 0.6% median - not worth it).
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
from reparam import slopes
W=30.0
NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")
LO=np.array(F.LO); HI=np.array(F.HI)
LO[3]=LO[7]=LO[11]=LO[17]=-12.0; HI[3]=HI[7]=HI[11]=HI[17]=12.0   # nQ of every block
LO[13]=0.0; HI[13]=20.0                                           # R_ET
FREE=[i for i in range(27) if i!=12]
def expand(x):
    p=np.zeros(27); p[FREE]=x; p[12]=0.0; return p
def pen(p):
    r=[]
    for x in (0.10,0.60):
        s=slopes(p,x)
        for n in NAMES: r.append(W*min(0.0,s[n]))
    return r
HA=dict(Q2=1.75,xB=0.36,mt=0.27,val=0.28,err=0.07)   # PRL 117 262001 / 118 222002
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"])
    sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else 9.9
def resid(x):
    p=expand(x); out=F.blocks(p); r=[]
    for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
    r.append((nratio(p)-HA["val"])/HA["err"])
    r+=pen(p)
    return np.array(r)
best=None
seeds=[]
for f in ("fitpar_amp2026s.npy","fitpar_prior_none.npy","fitpar_prior_halla.npy","fitpar_tied.npy","fitpar_amp2026_lx.npy"):
    try: q=np.array(np.load(f)); q[12]=0.0; seeds.append((f,q))
    except Exception: pass
for tag,q in seeds:
    q=np.clip(q,LO,HI)
    try:
        r=least_squares(resid,q[FREE],bounds=(LO[FREE],HI[FREE]),x_scale='jac',
                        xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
    except Exception as e:
        print("seed",tag,"failed",e); continue
    c=float(np.sum(r.fun**2)); print(f"  seed {tag:26s} chi2+pen {c:9.2f}",flush=True)
    if best is None or c<best[0]: best=(c,r.x)
p=expand(best[1]); out=F.blocks(p)
tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
print(f"\n=== PLAIN FIT + Hall-A neutron ratio === {n} pts, 26 par: "
      f"chi2/ndf = {tot:.1f}/{n-26} = {tot/(n-26):.4f}   penalty {sum(x*x for x in pen(p)):.3f}")
for k,v in out.items(): print(f"   {k:>8}: {sum(x*x for x in v):7.1f}/{len(v)}")
NAME=["HT_u N","HT_u b","HT_u b'","HT_u nQ","HT_d N","HT_d b","HT_d b'","HT_d nQ","ET_u N","ET_u b",
      "ET_u b'","ET_u nQ","ET_u b2","R_ET","db_ET","T00 N","T00 b","T00 nQ","R_L","delta0","delta1",
      "T00 b'","rho_CE","phi_CE","phi_w","rho_nf","phi_nf"]
print("\n   parameter        value      on a bound?")
for i,nm in enumerate(NAME):
    flag=""
    if i in FREE and min(abs(p[i]-LO[i]),abs(p[i]-HI[i]))/max(HI[i]-LO[i],1e-9)<0.005: flag="  <-- ON THE BOUND"
    print(f"   {nm:12s} {p[i]:10.4f}{flag}")
print("\n   slopes:", {k:round(v,2) for k,v in slopes(p,0.25).items()})
sp=amp.structure(p,"pi0p",-0.27,0.36,1.75); sn=amp.structure(p,"pi0n",-0.27,0.36,1.75)
print(f"   R_HT {p[4]/p[0]:+.3f}   R_ET {p[13]:.3f}   sigma_TT(n)/sigma_TT(p) = {sn['TT']/sp['TT']:.2f}"
      f"   (Hall-A: 0.28 +- 0.07)")
st=amp.structure(p,"pi0p",-0.4,0.25,1.94); print(f"   sigma_L/sigma_T = {st['L']/st['T']:.4f}")
np.save("fitpar_halla.npy",p)
