"""Hall-A in the fit with a normalisation nuisance.

The xlsx carries statistical errors only -- the syst columns are all zero -- so
the Hall-A block is currently weighted as if it had no systematics at all.  Two
repairs, run separately so their effects can be told apart:
  lam        one scale on the whole Hall-A block, Gaussian prior of the stated width
  syst10     10% of each value added in quadrature to its error, no scale
Hall-A and CLAS are different experiments with independent normalisations, so a
relative scale between them is legitimate -- unlike the eta/pi0 case, where both
channels came from the same run and the same elastic calibration.
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0; E_HALLA=5.55
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
def halla_res(p,lam=1.0,syst=0.0):
    r=[]
    for d in HALLA:
        ch="pi0p" if d["tgt"]=="p_T" else "pi0n"
        s=amp.structure(p,ch,-d["mt"],d["xB"],d["Q2"])
        if s is None: r+=[10.,10.,10.]; continue
        e=amp.epsilon(d["xB"],d["Q2"],E_HALLA)
        pred = s["T"] if d["tgt"]=="p_T" else s["T"]+e*s["L"]
        for val,err,mod in ((d["s"],d["ds"],pred),(d["LT"],d["dLT"],s["LT"]),(d["TT"],d["dTT"],s["TT"])):
            E=math.hypot(err,syst*abs(val))
            r.append((val-lam*mod)/E)
    return r
def run(tag, siglam, syst):
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.,-math.pi,-math.pi,-math.pi,0.5,0.5])
    HI=np.array(list(F.HI)+[ 12., 5., 8., 5., 8., math.pi, math.pi, math.pi, 2.0, 2.0])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0
    ties={7:3, 27:11, 6:2, 12:10}
    frozen={23,28,29,30,31,33}|(set() if siglam else {32})
    free=[i for i in range(37) if i not in set(ties)|set(frozen)]
    def expand(x):
        p=np.zeros(37); p[free]=x; p[23]=0.0
        if not siglam: p[32]=1.0
        for t,s in ties.items(): p[t]=p[s]
        return p
    def resid(x):
        p=expand(x); lam=p[32] if siglam else 1.0
        out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        r+=halla_res(p,lam,syst)
        if siglam: r.append((lam-1.0)/siglam)
        for xx in (0.05,0.35,1.00):
            s=sl(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in ("fitpar_all_3_phases_+_HallA.npy","fitpar_production.npy"):
        try: z=np.load(src)
        except Exception: continue
        q=np.ones(37); q[:len(z)]=z; q[32]=1.0
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); lam=p[32] if siglam else 1.0
    out=F.blocks(p); ch=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    hl=sum(x*x for x in halla_res(p,lam,syst)); nh=len(HALLA)*3
    s25=sl(p,0.25)
    sp=amp.structure(p,"pi0p",-0.27,0.36,1.75); sn=amp.structure(p,"pi0n",-0.27,0.36,1.75)
    st=amp.structure(p,"pi0p",-0.4,0.25,1.94)
    print(f"{tag:>18} {len(free):4d} {ch+hl:8.1f} {(ch+hl)/(n+nh-len(free)):8.4f} | lam {lam:5.3f} | "
          f"phases {p[34]:+5.2f} {p[35]:+5.2f} {p[36]:+5.2f} | ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f} "
          f"n/p {sn['TT']/sp['TT']:5.2f} L/T {st['L']/st['T']:.3f} | CLAS {ch:7.1f} HallA {hl:6.1f}/{nh}",flush=True)
    np.save(f"fitpar_han_{tag.replace(' ','_').replace('%','')}.npy",p); return p
print(f"{'variant':>18} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("no nuisance",None,0.0)
run("lam prior 10%",0.10,0.0)
run("lam prior 5%",0.05,0.0)
run("lam free",0.50,0.0)
run("syst 10% added",None,0.10)
