"""Night series 2, on top of the best xB-shape variant (F_split).

  halla    + the Hall-A neutron ratio 0.28 +- 0.07 (observable level, isospin only)
  bsa2x    De Masi BSA errors inflated x2 -- justified: its point-to-point scatter
           inside a bin is twice the quoted (figure-extracted) errors, and its
           pull rms is 1.95
  phase0   phi_CE frozen at 0 -- the three phases have an exactly flat direction
           (errors ~4e4), so one of them is pure gauge
  both     bsa2x + phase0 + Hall-A
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
def slopes32(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
LO=np.array(list(F.LO)+[-12.0,-5.0,-5.0,-5.0,-5.0]); HI=np.array(list(F.HI)+[12.0,5.0,8.0,5.0,8.0])
LO[3]=LO[7]=LO[11]=LO[17]=-12.0; HI[3]=HI[7]=HI[11]=HI[17]=12.0
LO[12]=-8.0; HI[12]=8.0; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.0; HI[14]=12.0
HA=dict(Q2=1.75,xB=0.36,mt=0.27,val=0.28,err=0.07)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else 9.9
def run(tag, halla=False, bsa2x=False, phase0=False):
    free=[i for i in range(32)]
    if phase0: free.remove(23)
    def expand(x):
        p=np.zeros(32); p[free]=x
        if phase0: p[23]=0.0
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        r+=out["pi0"]+out["eta"]
        r+=[v/2.0 for v in out["bsa_pi0"]] if bsa2x else out["bsa_pi0"]
        r+=out["bsa_eta"]+out["bsa_c12"]+out["eg1"]
        if halla: r.append((nratio(p)-HA["val"])/HA["err"])
        for xx in (0.10,0.60):
            s=slopes32(p,xx)
            for n in s: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for f in ("fitpar_night_F_split.npy","fitpar_night_F_common.npy","fitpar_etd_free.npy"):
        q=np.zeros(32); z=np.load(f); q[:len(z)]=z; q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=60000)
        except Exception as e: print("  seed fail",e); continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    # chi2 always reported with the ORIGINAL errors, so the numbers stay comparable
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2); s25=slopes32(p,0.25)
    print(f"{tag:>8} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f}  ETd/ETu(-t=.3) {ET[1]/ET[0]:+6.2f} "
          f" HTd/HTu(-t=.3) {HT[1]/HT[0]:+6.2f}  slopes u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:4.2f}  "
          f"n/p {nratio(p):5.2f}  pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f} "
          f"bsa {sum(x*x for x in out['bsa_pi0']):6.1f} eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
    np.save(f"fitpar_n2_{tag}.npy",p)
print("chi2 always quoted with the ORIGINAL De Masi errors.  Hall-A n/p = 0.28 +- 0.07\n")
print(f"{'variant':>8} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("base")
run("halla",halla=True)
run("bsa2x",bsa2x=True)
run("phase0",phase0=True)
run("all",halla=True,bsa2x=True,phase0=True)
