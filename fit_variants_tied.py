"""Variants A/B/C: how much do the data actually require an independent d shape?
A  27 par  nothing tied
B  25 par  H_T^d: b', c shared with u; free R_HT and db_H
C  23 par  GK-like: both d blocks are pure normalisations, db = 0
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import least_squares
import amplitudes as amp

E_XS,E_BSA,E_EG1,E_C12=5.75,5.776,5.9,10.6
KIN={"1.94":0.25,"2.83":0.40}
def load_sf(fn):
    rows=[]
    for line in open(fn):
        v=line.split()
        if len(v)<12: continue
        v=[float(x) for x in v[:12]]
        rows.append(dict(Q2=v[0],xB=v[1],t=-v[2],U=v[3],dU=math.hypot(v[4],v[5]),
                         LT=v[6],dLT=math.hypot(v[7],v[8]),TT=v[9],dTT=math.hypot(v[10],v[11])))
    return rows
PI0=load_sf("data/strfun_pi0.data"); ETA=load_sf("data/strfun_eta.data")
BSA=json.load(open("data/clas6_demasi_alu.json"))
BSA_PTS=[(b['xB'],b['Q2'],q['t'],q['alpha'],q['err']) for b in BSA for q in b['pts'] if q['err']>0]
EBSA=[(d['xB'],d['Q2'],d['t'],d['alpha'],math.hypot(d['stat'],d['syst'])) for d in json.load(open("data/clas6_zhao_eta_alu.json"))]
C12=json.load(open("data/clas12_kim_alu.json"))
C12_PTS=[(r['xB'],r['Q2'],q['t'],q['A'],q['dA']) for r in C12 for q in r['pts'] if q['dA']>0]
EG1=json.load(open("data/eg1dvcs_pi0_target_asym.json"))
RECS={"E154M5":"AULsin","E154M6":"AULsin","E154M7":"AULsin2","E154M8":"AULsin2",
      "E154M9":"ALLc","E154M10":"ALLc","E154M11":"ALLcos","E154M12":"ALLcos"}
EG1_PTS=[]
for rec,key in RECS.items():
    for row in EG1[rec]['rows']:
        Q2,mt,A,dA=row[0],row[1],row[2],row[3]
        ds=row[4] if len(row)>4 else 0.0
        EG1_PTS.append((key,Q2,KIN[f"{Q2:.2f}"],mt,A,math.hypot(dA,ds)))
def moment(s,key,e):
    s0=s["T"]+e*s["L"]
    return dict(AULsin=math.sqrt(2*e*(1+e))*s["LTUL"]/s0, AULsin2=e*s["TTUL"]/s0,
                ALLc=math.sqrt(1-e*e)*s["Tp"]/s0, ALLcos=math.sqrt(2*e*(1-e))*s["LTpLL"]/s0)[key]

FULL=np.load("fitpar_i2.npy")
# variant -> (free index list, expand function)
def make(variant):
    if variant=="A":
        free=list(range(27))
        def exp_(x):
            p=np.array(FULL); p[free]=x; return p
    elif variant=="B":
        free=[i for i in range(27) if i not in (6,7)]   # HTd_bx, HTd_c tied to u
        def exp_(x):
            p=np.array(FULL); p[free]=x; p[6]=p[2]; p[7]=p[3]; return p
    else:                                               # C: GK-like
        free=[i for i in range(27) if i not in (6,7,14)]  # + db(ET)=0
        def exp_(x):
            p=np.array(FULL); p[free]=x; p[6]=p[2]; p[7]=p[3]; p[14]=0.0
            p[5]=p[1]                                    # db_H = 0 as well
            return p
    if variant=="C": free=[i for i in free if i!=5]
    return free, exp_

def resid(p):
    r=[]
    for rows,ch in ((PI0,"pi0p"),(ETA,"etap")):
        for d in rows:
            s=amp.structure(p,ch,d["t"],d["xB"],d["Q2"])
            if s is None: r+=[10.,10.,10.]; continue
            e=amp.epsilon(d["xB"],d["Q2"],E_XS)
            r+=[(d["U"]-(s["T"]+e*s["L"]))/d["dU"],(d["TT"]-s["TT"])/d["dTT"],(d["LT"]-s["LT"])/d["dLT"]]
    for pts,ch,E in ((BSA_PTS,"pi0p",E_BSA),(EBSA,"etap",E_BSA),(C12_PTS,"pi0p",E_C12)):
        for (xB,Q2,t,a,da) in pts:
            v=amp.bsa_sinphi(p,ch,-t,xB,Q2,E); r.append(((a-v)/da) if v is not None else 5.0)
    for (key,Q2,xB,mt,A,dA) in EG1_PTS:
        s=amp.structure(p,"pi0p",-mt,xB,Q2)
        r.append(((A-moment(s,key,amp.epsilon(xB,Q2,E_EG1)))/dA) if s else 5.0)
    r.append((p[13]-0.54)/0.15)          # forward-limit prior on R_ET
    return np.array(r)

LO=np.array([-1e4,-2.,-8.,-6.]*2+[-1e4,-2.,-8.,-6.,-3.]+[0.05,-2.0]+[0.,-2.,-6.,-2.,-math.pi,-3.]+[-8.]+[0.,-math.pi,-math.pi,0.,-math.pi])
HI=np.array([ 1e4,12., 8., 6.]*2+[ 1e4,12., 8., 6., 3.]+[3.00, 8.0]+[1e4, 8., 6., 2., math.pi, 3.]+[ 8.]+[1.5, math.pi, math.pi, 3., math.pi])
print(f"{'variant':>8} {'npar':>5} {'chi2':>9} {'ndf':>5} {'chi2/ndf':>9} {'R_HT(t=0)':>10} {'R_ET(t=0)':>10} {'db_ET':>7}")
print("-"*76)
res_store={}
for V in ("A","B","C"):
    free,exp_=make(V)
    x0=FULL[free]
    def f(x): return resid(exp_(x))
    r=least_squares(f,x0,bounds=(LO[free],HI[free]),x_scale='jac',
                    xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
    p=exp_(r.x); c=float(np.sum(r.fun**2)); n=len(r.fun)-len(free)
    RHT=p[4]/p[0]; RET=p[13]
    print(f"{V:>8} {len(free):5d} {c:9.1f} {n:5d} {c/n:9.4f} {RHT:10.3f} {RET:10.3f} {p[14]:7.3f}  pull {(p[13]-0.54)/0.15:+5.1f}")
    res_store[V]=p
    np.save(f"fitpar_{V}_prior.npy",p)
print()
a,b,cc=res_store["A"],res_store["B"],res_store["C"]
print("H_T^d shape:  b       b'      c")
for V,p in (("A",a),("B",b),("C",cc)):
    print(f"    {V}:  {p[5]:7.3f} {p[6]:7.3f} {p[7]:7.3f}")
