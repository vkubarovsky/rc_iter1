"""Slope-positivity constraint: b + b' ln xB > 0 for xB in [0.1, 0.6].

The parameterisation is now exp[(b + b' ln xB) t] (reparam.py).  A negative
effective slope means a form factor that GROWS with |t| - amp2026 does that for
H_T^d over the whole measured range.  The slope is linear in ln xB, so imposing
the two endpoints xB = 0.1 and xB = 0.6 constrains the whole interval.

    SMIN=0.0 OUTP=fitpar_slope.npy ~/.venv/bin/python3 fit_slope.py

Seed and data are the amp2026 fixed point (fitpar_amp2026_lx.npy, data/ = i2).
"""
import json, math, os
import numpy as np
from scipy.optimize import least_squares
import amplitudes as amp
from reparam import BLOCKS, slopes

SEED=os.environ.get("SEED","fitpar_amp2026_lx.npy")
OUTP=os.environ.get("OUTP","fitpar_slope.npy")
SMIN=float(os.environ.get("SMIN","0.0"))     # required minimum slope, GeV^-2
WPEN=float(os.environ.get("WPEN","30.0"))    # penalty stiffness
XLO,XHI=0.10,0.60

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

def blocks(p):
    out={}
    for rows,ch,tag in ((PI0,"pi0p","pi0"),(ETA,"etap","eta")):
        rr=[]
        for d in rows:
            s=amp.structure(p,ch,d["t"],d["xB"],d["Q2"])
            if s is None: rr+=[10.,10.,10.]; continue
            e=amp.epsilon(d["xB"],d["Q2"],E_XS)
            rr+=[(d["U"]-(s["T"]+e*s["L"]))/d["dU"],(d["TT"]-s["TT"])/d["dTT"],(d["LT"]-s["LT"])/d["dLT"]]
        out[tag]=rr
    for pts,ch,E,tag in ((BSA_PTS,"pi0p",E_BSA,"bsa_pi0"),(EBSA,"etap",E_BSA,"bsa_eta"),
                         (C12_PTS,"pi0p",E_C12,"bsa_c12")):
        rr=[]
        for (xB,Q2,t,a,da) in pts:
            v=amp.bsa_sinphi(p,ch,-t,xB,Q2,E)
            rr.append(((a-v)/da) if v is not None else 5.0)
        out[tag]=rr
    rr=[]
    for (key,Q2,xB,mt,A,dA) in EG1_PTS:
        s=amp.structure(p,"pi0p",-mt,xB,Q2)
        rr.append(((A-moment(s,key,amp.epsilon(xB,Q2,E_EG1)))/dA) if s else 5.0)
    out["eg1"]=rr
    return out

def slope_pen(p, smin=None, w=None):
    """One residual per (block, endpoint): w*min(0, slope - smin)."""
    smin = SMIN if smin is None else smin
    w = WPEN if w is None else w
    r=[]
    for x in (XLO,XHI):
        s=slopes(p,x)
        for name in ("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00"):
            r.append(w*min(0.0, s[name]-smin))
    return r

def resid(p):
    out=blocks(p)
    r=[]
    for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
    r.append((p[13]-0.54)/0.15)            # forward-limit prior on R_ET
    r+=slope_pen(p)
    return np.array(r)

LO=np.array([-1e4,-2.,-8.,-6.]*2+[-1e4,-2.,-8.,-6.,-3.]+[0.05,-2.0]+[0.,-2.,-6.,-2.,-math.pi,-3.]+[-8.]+[0.,-math.pi,-math.pi,0.,-math.pi])
HI=np.array([ 1e4,12., 8., 6.]*2+[ 1e4,12., 8., 6., 3.]+[3.00, 8.0]+[1e4, 8., 6., 2., math.pi, 3.]+[ 8.]+[1.5, math.pi, math.pi, 3., math.pi])

def report(tag,p,p0=None):
    out=blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    pen=sum(x*x for x in slope_pen(p)); pri=((p[13]-0.54)/0.15)**2
    print(f"\n=== {tag} ===  {n} pts, 27 par: chi2/ndf = {tot:.1f}/{n-27} = {tot/(n-27):.4f}"
          f"   (R_ET prior {pri:.2f}, slope penalty {pen:.3f})")
    for k,v in out.items(): print(f"   {k:>8}: {sum(x*x for x in v):7.1f}/{len(v)}")
    XS=(0.10,0.15,0.25,0.40,0.60)
    tab={x:slopes(p,x) for x in XS}
    print(f"   {'block':>10} {'b':>8} {'b_prime':>8} |" + "".join(f"{x:8.2f}" for x in XS))
    for name in ("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00"):
        ib,ibp=BLOCKS.get(name,(9,10))
        b=p[ib]+(p[14] if name=="Ebar_T^d" else 0.0); bp=p[ibp] if ibp<len(p) else 0.0
        print(f"   {name:>10} {b:8.3f} {bp:8.3f} |" + "".join(f"{tab[x][name]:8.2f}" for x in XS))
    s=amp.structure(p,"pi0p",-0.4,0.25,1.94)
    print(f"   R_HT(t=0) = {p[4]/p[0]:+.3f}   R_ET = {p[13]:.3f} (pull {(p[13]-0.54)/0.15:+.1f})"
          f"   db_ET = {p[14]:.3f}   sigma_L/sigma_T(1.94,0.25,0.4) = {s['L']/s['T']:.4f}")
    return tot,n

if __name__=="__main__":
    p0=np.load(SEED)
    print(f"seed {SEED}, data/ = iteration-2, SMIN={SMIN}, WPEN={WPEN}")
    report("SEED (amp2026, unconstrained)",p0)
    # seeds: as-is, and two with H_T^d forced into the allowed region
    starts=[np.array(p0)]
    for bd,bpd in ((1.0,0.0),(2.0,-0.5)):
        q=np.array(p0); q[5],q[6]=bd,bpd; starts.append(q)
    best=None
    for i,sd in enumerate(starts):
        try:
            r=least_squares(resid,np.clip(sd,LO,HI),bounds=(LO,HI),x_scale='jac',
                            xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception as e:
            print("start",i,"failed:",e); continue
        cc=float(np.sum(r.fun**2))
        print(f"   start {i}: chi2+pen = {cc:.2f}  nfev={r.nfev}")
        if best is None or cc<best[0]: best=(cc,r.x)
    _,p=best
    report(f"CONSTRAINED (slope > {SMIN})",p)
    np.save(OUTP,p); print(f"\nsaved {OUTP}")
