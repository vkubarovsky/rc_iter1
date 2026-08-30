"""chi2 profile in the relative u-d phase of the Ebar_T sector.

Fix phi, refit everything else, and see whether the data determine it.
Built on model B (sector-common alpha', 23 par); phi is the 24th.
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
    a=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); b=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return b["TT"]/a["TT"] if a and b else float('nan')
LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.,-math.pi,-math.pi,-math.pi,-math.pi,-math.pi])
HI=np.array(list(F.HI)+[ 12., 5., 8., 5., 8., math.pi, math.pi, math.pi, math.pi, math.pi])
LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
for i in BSLOT: LO[i]=0.0
ties={7:3, 27:11, 6:2, 12:10}
free=[i for i in range(37) if i not in set(ties)|{23,28,29,30,31,32,33,34,35,36}]
print(f"{'phi [rad]':>9} {'deg':>5} {'chi2':>8} {'d_chi2':>8} {'slope ET d':>11} {'ETd/ETu':>8} {'n/p':>6} "
      f"{'pi0':>7} {'eta':>7} {'bsa':>7} {'eg1':>6}")
base=None
for phi in (0.0,0.2,0.4,0.6,0.8,1.0,1.28,1.6,2.0,2.4,2.8):
    def expand(x,phi=phi):
        p=np.zeros(37); p[free]=x; p[23]=0.0; p[35]=phi
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
    for src in ("fitpar_production.npy","fitpar_ph_ET.npy"):
        q=np.zeros(37); z=np.load(src); q[:len(z)]=z
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values())
    if base is None: base=tot
    s25=sl(p,0.25); HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    print(f"{phi:9.2f} {math.degrees(phi):5.0f} {tot:8.1f} {tot-base:+8.1f} {s25['Ebar_T^d']:11.2f} "
          f"{abs(ET[1])/ET[0]:8.2f} {nratio(p):6.2f} {sum(x*x for x in out['pi0']):7.1f} "
          f"{sum(x*x for x in out['eta']):7.1f} {sum(sum(x*x for x in out[k]) for k in ('bsa_pi0','bsa_eta','bsa_c12')):7.1f} "
          f"{sum(x*x for x in out['eg1']):6.1f}",flush=True)
    np.save(f"fitpar_scanphi_{phi:.2f}.npy",p)
