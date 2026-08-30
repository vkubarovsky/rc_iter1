"""Production fit with the Hall-A neutron in it.

Neutron: PRL 118 222002, 4 points at Q2 = 1.75, xB = 0.36, sigma_U, sigma_LT,
sigma_TT.  The xlsx has no systematics, so 10% of each value is added in
quadrature to the statistical error.  It is the only measurement anywhere in the
fit of the u + 2d combination, and the blind prediction of our previous model
missed it by 3 to 9 sigma per point.

All three flavour phases free (H_T, Ebar_T, T00).

    n_only    neutron added
    n_and_p   neutron and the Hall-A proton (PRL 117 262001, L/T separated)
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0; E_HALLA=5.55; SYST=0.10
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]
def sl(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
HALLA=[]
for line in open("data/halla_pi0.data"):
    if line.startswith("#") or not line.strip(): continue
    v=line.split()
    HALLA.append(dict(tgt=v[0],Q2=float(v[1]),xB=float(v[2]),mt=float(v[3]),
                      s=float(v[4]),ds=float(v[5]),LT=float(v[6]),dLT=float(v[7]),
                      TT=float(v[8]),dTT=float(v[9])))
def ha_res(p, which):
    r=[]
    for d in HALLA:
        if which=="n_only" and d["tgt"]!="n_U": continue
        ch="pi0p" if d["tgt"]=="p_T" else "pi0n"
        s=amp.structure(p,ch,-d["mt"],d["xB"],d["Q2"])
        if s is None: r+=[10.,10.,10.]; continue
        e=amp.epsilon(d["xB"],d["Q2"],E_HALLA)
        pred = s["T"] if d["tgt"]=="p_T" else s["T"]+e*s["L"]
        for val,err,mod in ((d["s"],d["ds"],pred),(d["LT"],d["dLT"],s["LT"]),(d["TT"],d["dTT"],s["TT"])):
            r.append((val-mod)/math.hypot(err,SYST*abs(val)))
    return r
def run(tag, which):
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
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        r+=ha_res(p,which)
        for xx in (0.05,0.35,1.00):
            s=sl(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src,ph in (("fitpar_han_syst_10_added.npy",None),("fitpar_production.npy",None),
                   ("fitpar_production.npy",0.3),("fitpar_production.npy",0.0)):
        try: z=np.load(src)
        except Exception: continue
        q=np.zeros(37); q[:len(z)]=z
        if ph is not None: q[34]=q[35]=q[36]=ph
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    ch=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    hl=sum(x*x for x in ha_res(p,which)); nh=len(ha_res(p,which))
    s25=sl(p,0.25); HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    sp=amp.structure(p,"pi0p",-0.27,0.36,1.75); sn=amp.structure(p,"pi0n",-0.27,0.36,1.75)
    st=amp.structure(p,"pi0p",-0.4,0.25,1.94)
    print(f"\n=== {tag} === {n+nh} pts, {len(free)} par: chi2/ndf = {ch+hl:.1f}/{n+nh-len(free)} = {(ch+hl)/(n+nh-len(free)):.4f}")
    for k,v in out.items(): print(f"   {k:>8}: {sum(x*x for x in v):7.1f}/{len(v)}")
    print(f"   {'HallA':>8}: {hl:7.1f}/{nh}")
    print(f"   phases  H_T {p[34]:+.3f}   Ebar_T {p[35]:+.3f}   T00 {p[36]:+.3f}  rad")
    print(f"   slopes at xB=0.25: " + "  ".join(f"{k} {v:.2f}" for k,v in s25.items()))
    print(f"   d/u at -t=0.3:  H_T {abs(HT[1])/HT[0]*np.sign(HT[1].real if hasattr(HT[1],'real') else HT[1]):+.3f}"
          f"   |Ebar_T| {abs(ET[1])/ET[0]:.3f}   n/p {sn['TT']/sp['TT']:.3f}   sigma_L/sigma_T {st['L']/st['T']:.4f}")
    np.save(f"fitpar_n_{tag}.npy",p); return p
run("n_only","n_only")
run("n_and_p","all")
