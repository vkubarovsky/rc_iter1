"""Add CLAS12 BSA (PLB 849 138459, E=10.6 GeV, eps=0.83-0.87) to the global fit.
This is the second-epsilon lever on sigma_L/sigma_T."""
import json, math
import numpy as np
from scipy.optimize import least_squares
import amplitudes as amp
import os
# seeds must be in the ln xB slope convention (reparam.py); the pre-2026-08-28
# files (fitpar_i1/i2, fitpar_amp2021_published, ~/pi0_eta_amplitude_model) are NOT.
SEED=os.environ.get("SEED","/Users/vpk/rc_iter1/fitpar_slope.npy")
OUTP=os.environ.get("OUTP","fitpar_i1.npy")

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
    r=[]
    for rows,ch,tag in ((PI0,"pi0p","pi0"),(ETA,"etap","eta")):
        rr=[]
        for d in rows:
            s=amp.structure(p,ch,d["t"],d["xB"],d["Q2"])
            if s is None: rr+=[10.,10.,10.]; continue
            e=amp.epsilon(d["xB"],d["Q2"],E_XS)
            rr+=[(d["U"]-(s["T"]+e*s["L"]))/d["dU"],(d["TT"]-s["TT"])/d["dTT"],(d["LT"]-s["LT"])/d["dLT"]]
        out[tag]=rr; r+=rr
    for pts,ch,E,tag in ((BSA_PTS,"pi0p",E_BSA,"bsa_pi0"),(EBSA,"etap",E_BSA,"bsa_eta"),
                         (C12_PTS,"pi0p",E_C12,"bsa_c12")):
        rr=[]
        for (xB,Q2,t,a,da) in pts:
            v=amp.bsa_sinphi(p,ch,-t,xB,Q2,E)
            rr.append(((a-v)/da) if v is not None else 5.0)
        out[tag]=rr; r+=rr
    rr=[]
    for (key,Q2,xB,mt,A,dA) in EG1_PTS:
        s=amp.structure(p,"pi0p",-mt,xB,Q2)
        rr.append(((A-moment(s,key,amp.epsilon(xB,Q2,E_EG1)))/dA) if s else 5.0)
    out["eg1"]=rr; r+=rr
    return out, np.array(r)

def residuals(p, use_c12=True):
    out,_=blocks(p)
    keys=["pi0","eta","bsa_pi0","bsa_eta","eg1"]+(["bsa_c12"] if use_c12 else [])
    r=[]
    for k in keys: r+=out[k]
    r.append((p[13]-0.54)/0.15)
    return np.array(r)

p0=np.load(SEED)
print("=== PREDICTION for CLAS12 BSA (not in the fit) ===")
out,_=blocks(p0)
c=sum(x*x for x in out["bsa_c12"]); print(f"chi2 = {c:.1f}/{len(out['bsa_c12'])}")
for (xB,Q2,t,a,da),res in zip(C12_PTS,out["bsa_c12"]):
    v=amp.bsa_sinphi(p0,"pi0p",-t,xB,Q2,E_C12)
    print(f"  Q2={Q2:4.2f} xB={xB:4.2f} t={t:4.2f}: data {a:+.3f}+-{da:.3f}  model {v:+.3f}  pull {res:+.1f}")

LO=[-1e4,-2.,-8.,-6.]*2+[-1e4,-2.,-8.,-6.,-3.]+[0.05,-2.0]+[0.,-2.,-6.,-2.,-math.pi,-3.]+[-8.]+[0.,-math.pi,-math.pi,0.,-math.pi]
HI=[ 1e4,12., 8., 6.]*2+[ 1e4,12., 8., 6., 3.]+[3.00, 8.0]+[1e4, 8., 6., 2., math.pi, 3.]+[ 8.]+[1.5, math.pi, math.pi, 3., math.pi]
best=None
for sd in (list(p0), list(p0[:15])+[p0[15]*0.5]+list(p0[16:])):
    try:
        res=least_squares(residuals,sd,bounds=(LO,HI),x_scale='jac',xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
    except Exception: continue
    cc=float(np.sum(res.fun**2))
    if best is None or cc<best[0]: best=(cc,res.x)
_,p=best
out,_=blocks(p)
tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
print(f"\n=== REFIT WITH CLAS12 ===  {n} pts, 27 par: chi2/ndf = {tot:.1f}/{n-27} = {tot/(n-27):.3f}")
for k,v in out.items(): print(f"   {k:>8}: {sum(x*x for x in v):7.1f}/{len(v)}")
np.save(OUTP,p)
s=amp.structure(p,"pi0p",-0.4,0.25,1.94); s0=amp.structure(p0,"pi0p",-0.4,0.25,1.94)
print(f"\nsigma_L/sigma_T at (1.94,0.25,0.4): {s0['L']/s0['T']:.4f} -> {s['L']/s['T']:.4f}")
print(f"L_N {p0[15]:.2f}->{p[15]:.2f}   rho_nf {p0[25]:.3f}->{p[25]:.3f}   delta0 {p0[19]:.3f}->{p[19]:.3f}")
