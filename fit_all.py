"""Everything free that should be free, and the Hall-A data in the fit.

New in this fit:
  * all three flavour phases free (H_T, Ebar_T, T00): the flavour convolutions are
    complex, and the data determine the phases (the Ebar_T one to 1.28 +- 0.09,
    4.8 sigma from zero).  Computing them instead was wrong: the computation needs
    an x shape, and that shape already fails for delta, where the data give 2.05
    model-independently against 0.29 computed.
  * Hall-A pi0, from HallA_pi0_p_n_2016_2017.xlsx:
      proton  PRL 117 262001, L/T SEPARATED -> the sigma column is sigma_T itself
      neutron PRL 118 222002, unseparated   -> sigma_U
    10 points x (sigma, sigma_LT, sigma_TT).  The neutron is the only direct
    measurement of the u + 2d combination anywhere in the fit.
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
print(f"Hall-A: {sum(1 for d in HALLA if d['tgt']=='p_T')} proton (sigma_T) + "
      f"{sum(1 for d in HALLA if d['tgt']=='n_U')} neutron (sigma_U) points\n")
def halla_res(p):
    r=[]
    for d in HALLA:
        ch="pi0p" if d["tgt"]=="p_T" else "pi0n"
        s=amp.structure(p,ch,-d["mt"],d["xB"],d["Q2"])
        if s is None: r+=[10.,10.,10.]; continue
        e=amp.epsilon(d["xB"],d["Q2"],E_HALLA)
        pred = s["T"] if d["tgt"]=="p_T" else s["T"]+e*s["L"]
        r+=[(d["s"]-pred)/d["ds"], (d["LT"]-s["LT"])/d["dLT"], (d["TT"]-s["TT"])/d["dTT"]]
    return r
def run(tag, phases, halla):
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.,-math.pi,-math.pi,-math.pi,-math.pi,-math.pi])
    HI=np.array(list(F.HI)+[ 12., 5., 8., 5., 8., math.pi, math.pi, math.pi, math.pi, math.pi])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0
    ties={7:3, 27:11, 6:2, 12:10}
    frozen={23,28,29,30,31,32,33}|({34,35,36}-set(phases))
    free=[i for i in range(37) if i not in set(ties)|frozen]
    def expand(x):
        p=np.zeros(37); p[free]=x; p[23]=0.0
        for t,s in ties.items(): p[t]=p[s]
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        if halla: r+=halla_res(p)
        for xx in (0.05,0.35,1.00):
            s=sl(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src,ph in (("fitpar_production.npy",None),("fitpar_ph_all three.npy",None),("fitpar_production.npy",0.6)):
        try: z=np.load(src)
        except Exception: continue
        q=np.zeros(37); q[:len(z)]=z
        if ph is not None:
            for i in phases: q[i]=ph
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    ch=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    hl=sum(x*x for x in halla_res(p)); nh=len(halla_res(p))
    tot=ch+(hl if halla else 0.0); N=n+(nh if halla else 0)
    s25=sl(p,0.25); HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    sp=amp.structure(p,"pi0p",-0.27,0.36,1.75); sn=amp.structure(p,"pi0n",-0.27,0.36,1.75)
    st=amp.structure(p,"pi0p",-0.4,0.25,1.94)
    print(f"{tag:>16} {len(free):4d} {tot:8.1f} {tot/(N-len(free)):8.4f} | phases {p[34]:+5.2f} {p[35]:+5.2f} {p[36]:+5.2f} | "
          f"slope ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f}  n/p {sn['TT']/sp['TT']:5.2f}  L/T {st['L']/st['T']:.3f} | "
          f"CLAS {ch:7.1f} HallA {hl:6.1f}/{nh}",flush=True)
    np.save(f"fitpar_all_{tag.replace(' ','_').replace(',','')}.npy",p); return p
print(f"{'variant':>16} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("1 phase, no HA",[35],False)
run("3 phases, no HA",[34,35,36],False)
run("1 phase + HallA",[35],True)
run("3 phases + HallA",[34,35,36],True)
