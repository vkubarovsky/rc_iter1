"""Hall-A properly: both targets, both L/T separated, sigma column = sigma_T.

Correction found 2026-08-31: the neutron paper (PRL 118 222002) is Rosenbluth
separated exactly like the proton one (PRL 117 262001) -- its four structure
functions are dsigma_T, dsigma_L, dsigma_TL, dsigma_TT.  We had been comparing
its sigma column against sigma_T + eps sigma_L.  Verified against their Fig. 5:
the four values 0.26, 0.50, 0.26, 0.41 ub/GeV^2 are the blue squares of dsigma_T.

Their sigma_L is published in the figure but is not in the spreadsheet, so it is
still missing from the fit -- it would be the only direct constraint on the
longitudinal sector anywhere in this analysis.

    n_only     neutron sigma_T, sigma_LT, sigma_TT
    n_and_p    plus the Hall-A proton, same three
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0; SYST=0.10
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]
def sl(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
HA=[]
for line in open("data/halla_pi0.data"):
    if line.startswith("#") or not line.strip(): continue
    v=line.split()
    HA.append(dict(tgt=v[0],Q2=float(v[1]),xB=float(v[2]),mt=float(v[3]),
                   s=float(v[4]),ds=float(v[5]),LT=float(v[6]),dLT=float(v[7]),
                   TT=float(v[8]),dTT=float(v[9])))
def ha_res(p,which):
    r=[]
    for d in HA:
        if which=="n_only" and d["tgt"]!="n_T": continue
        ch="pi0p" if d["tgt"]=="p_T" else "pi0n"
        s=amp.structure(p,ch,-d["mt"],d["xB"],d["Q2"])
        if s is None: r+=[10.,10.,10.]; continue
        for val,err,mod in ((d["s"],d["ds"],s["T"]),          # BOTH targets: sigma_T
                            (d["LT"],d["dLT"],s["LT"]),(d["TT"],d["dTT"],s["TT"])):
            r.append((val-mod)/math.hypot(err,SYST*abs(val)))
    return r
def run(tag,which):
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
    for src in ("fitpar_production.npy","fitpar_as_sector.npy"):
        q=np.zeros(37); z=np.load(src); q[:len(z)]=z
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
    print(f"{tag:>9} {len(free):4d} {ch+hl:8.1f} {(ch+hl)/(n+nh-len(free)):8.4f} | phase ET {p[35]:+.3f} | "
          f"ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f} |ETd/u| {abs(ET[1])/ET[0]:5.2f} "
          f"n/p {sn['TT']/sp['TT']:5.2f} L/T {st['L']/st['T']:.3f} | CLAS {ch:7.1f} HallA {hl:6.1f}/{nh} "
          f"eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
    np.save(f"fitpar_hf_{tag}.npy",p); return p
print(f"{'variant':>9} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("n_only","n_only")
run("n_and_p","all")
