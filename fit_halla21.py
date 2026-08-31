"""Add the Hall-A 12-GeV data (HallA_y21) and the 6-GeV set (HallA_y11).

From hepgen_mac/data/All_experiment.xlsx.  HallA_y21 reaches Q2 = 8.31 and
xB = 0.60 -- our CLAS-based fit stops at Q2 = 4.2 and xB = 0.54, so this doubles
the Q2 lever arm and is the first real test of the Q^nQ behaviour outside the
region where it was fitted.  These sets are UNSEPARATED: the sigma column is
sigma_U = sigma_T + eps sigma_L, with eps from the per-row beam energy
(4.49 to 10.99 GeV).

    base      as now: CLAS + Hall-A neutron
    +y21      add the 36 Hall-A 12-GeV points
    +y21+y11  add the 14 Hall-A 6-GeV points as well
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0; SYST_N=0.10
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]
def sl(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
NEUT=[]
for line in open("data/halla_pi0.data"):
    if line.startswith("#") or not line.strip(): continue
    v=line.split()
    if v[0]=="n_T":
        NEUT.append(dict(Q2=float(v[1]),xB=float(v[2]),mt=float(v[3]),s=float(v[4]),ds=float(v[5]),
                         LT=float(v[6]),dLT=float(v[7]),TT=float(v[8]),dTT=float(v[9])))
MORE=[]
for line in open("data/halla_more.data"):
    if line.startswith("#") or not line.strip(): continue
    v=line.split()
    MORE.append(dict(exp=v[0],Q2=float(v[3]),xB=float(v[4]),mt=float(v[5]),
                     U=float(v[6]),sU=float(v[7]),yU=float(v[8]),
                     LT=float(v[9]),sLT=float(v[10]),yLT=float(v[11]),
                     TT=float(v[12]),sTT=float(v[13]),yTT=float(v[14]),E=float(v[15])))
def n_res(p):
    r=[]
    for d in NEUT:
        s=amp.structure(p,"pi0n",-d["mt"],d["xB"],d["Q2"])
        if s is None: r+=[10.]*3; continue
        for val,err,mod in ((d["s"],d["ds"],s["T"]),(d["LT"],d["dLT"],s["LT"]),(d["TT"],d["dTT"],s["TT"])):
            r.append((val-mod)/math.hypot(err,SYST_N*abs(val)))
    return r
def more_res(p,which):
    r=[]
    for d in MORE:
        if d["exp"] not in which: continue
        s=amp.structure(p,"pi0p",-d["mt"],d["xB"],d["Q2"])
        if s is None: r+=[10.]*3; continue
        e=amp.epsilon(d["xB"],d["Q2"],d["E"])
        for val,st,sy,mod in ((d["U"],d["sU"],d["yU"],s["T"]+e*s["L"]),
                              (d["LT"],d["sLT"],d["yLT"],s["LT"]),
                              (d["TT"],d["sTT"],d["yTT"],s["TT"])):
            r.append((val-mod)/math.hypot(st,sy))
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
        r+=n_res(p); r+=more_res(p,which)
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
    cn=sum(x*x for x in n_res(p)); nn=len(n_res(p))
    c21=sum(x*x for x in more_res(p,("HallA_y21",))); n21=len(more_res(p,("HallA_y21",)))
    c11=sum(x*x for x in more_res(p,("HallA_y11",))); n11=len(more_res(p,("HallA_y11",)))
    used=cn+(c21 if "HallA_y21" in which else 0)+(c11 if "HallA_y11" in which else 0)
    nused=nn+(n21 if "HallA_y21" in which else 0)+(n11 if "HallA_y11" in which else 0)
    s25=sl(p,0.25); HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    print(f"{tag:>10} {len(free):4d} {ch+used:8.1f} {(ch+used)/(n+nused-len(free)):8.4f} | "
          f"nQ H_T {p[3]:+5.2f} ET {p[11]:+5.2f} | ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f} | "
          f"CLAS {ch:7.1f} n {cn:5.1f}/{nn} y21 {c21:7.1f}/{n21} y11 {c11:6.1f}/{n11}",flush=True)
    np.save(f"fitpar_h21_{tag}.npy",p); return p
print("y21 and y11 chi2 are always shown, whether or not they are in the fit\n")
print(f"{'variant':>10} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("base",())
run("+y21",("HallA_y21",))
run("+y21+y11",("HallA_y21","HallA_y11"))
