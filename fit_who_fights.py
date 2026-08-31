"""Which data set fights the neutron?

The Hall-A neutron is the only measurement of u + 2d anywhere.  With everything
in the fit it is described at chi2 = 13.2 on its 4 sigma_TT points, all four pulls
positive, while P2 -- an earlier fit that never saw it -- gets 1.2.  So something
we added since P2 pushes the d component up.  Each row drops one block and asks
where the neutron then lands.
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0; E_HALLA=5.55; SYST=0.10
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]
def sl(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
NEUT=[dict(zip(("Q2","xB","mt","s","ds","LT","dLT","TT","dTT"),map(float,l.split()[1:])))
      for l in open("data/halla_pi0.data") if l.startswith("n_U")]
def ha_res(p):
    r=[]
    for d in NEUT:
        s=amp.structure(p,"pi0n",-d["mt"],d["xB"],d["Q2"]); e=amp.epsilon(d["xB"],d["Q2"],E_HALLA)
        for val,err,mod in ((d["s"],d["ds"],s["T"]+e*s["L"]),(d["LT"],d["dLT"],s["LT"]),(d["TT"],d["dTT"],s["TT"])):
            r.append((val-mod)/math.hypot(err,SYST*abs(val)))
    return r
def tt_chi(p):
    c=0
    for d in NEUT:
        s=amp.structure(p,"pi0n",-d["mt"],d["xB"],d["Q2"])
        c+=((d["TT"]-s["TT"])/math.hypot(d["dTT"],SYST*abs(d["TT"])))**2
    return c
ALL=("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1")
def run(tag, drop):
    keep=[k for k in ALL if k not in drop]
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.,-math.pi,-math.pi,-math.pi,0.5,0.5])
    HI=np.array(list(F.HI)+[ 12., 5., 8., 5., 8., math.pi, math.pi, math.pi, 2.0, 2.0])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0
    ties={7:3, 27:11, 6:2, 12:10}
    free=[i for i in range(37) if i not in set(ties)|{23,28,29,30,31,32,33}]
    def expand(x):
        p=np.zeros(37); p[free]=x; p[23]=0.0
        for t,s in ties.items(): p[t]=p[s]
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in keep: r+=out[k]
        r+=ha_res(p)
        for xx in (0.05,0.35,1.00):
            s=sl(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in ("fitpar_production.npy","fitpar_as_sector.npy"):
        q=np.zeros(37); z=np.load(src); q[:len(z)]=z
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    s25=sl(p,0.25); HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    sp=amp.structure(p,"pi0p",-0.27,0.36,1.75); sn=amp.structure(p,"pi0n",-0.27,0.36,1.75)
    print(f"{tag:>16} | neutron sigma_TT chi2 {tt_chi(p):6.1f}/4   n/p {sn['TT']/sp['TT']:5.2f} | "
          f"|ETd/ETu| {abs(ET[1])/ET[0]:5.2f}  slope ET d {s25['Ebar_T^d']:5.2f} | "
          + "  ".join(f"{k} {sum(x*x for x in out[k]):6.1f}" for k in ALL),flush=True)
    np.save(f"fitpar_wf_{tag.replace(' ','_')}.npy",p)
print("Hall-A neutron is IN every fit below; the listed block is removed.")
print("P2, for reference: neutron sigma_TT chi2 = 1.2, n/p = 0.50.  Measured n/p = 0.28 +- 0.07.\n")
run("nothing dropped",())
run("no eg1",("eg1",))
run("no BSA",("bsa_pi0","bsa_eta","bsa_c12"))
run("no CLAS12",("bsa_c12",))
run("no eg1, no BSA",("eg1","bsa_pi0","bsa_eta","bsa_c12"))
run("only xsec",("eg1","bsa_pi0","bsa_eta","bsa_c12"))
