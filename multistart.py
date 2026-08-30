"""Is the production minimum global?  Refit from many scattered starts.

Each start perturbs the solution: normalisations x [0.5, 2], slopes +- 0.5,
Q2 exponents +- 0.6, all phases drawn uniformly in [-pi, pi].  If the fit is
well posed, they should come back to the same chi2 and the same parameters.
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0; E_HALLA=5.55; SYST=0.10
rng=np.random.default_rng(20260830)
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]
def sl(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
HALLA=[l.split() for l in open("data/halla_pi0.data") if l.strip() and not l.startswith("#")]
NEUT=[dict(Q2=float(v[1]),xB=float(v[2]),mt=float(v[3]),s=float(v[4]),ds=float(v[5]),
           LT=float(v[6]),dLT=float(v[7]),TT=float(v[8]),dTT=float(v[9])) for v in HALLA if v[0]=="n_U"]
def ha_res(p):
    r=[]
    for d in NEUT:
        s=amp.structure(p,"pi0n",-d["mt"],d["xB"],d["Q2"])
        if s is None: r+=[10.,10.,10.]; continue
        e=amp.epsilon(d["xB"],d["Q2"],E_HALLA)
        for val,err,mod in ((d["s"],d["ds"],s["T"]+e*s["L"]),(d["LT"],d["dLT"],s["LT"]),(d["TT"],d["dTT"],s["TT"])):
            r.append((val-mod)/math.hypot(err,SYST*abs(val)))
    return r
LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.,-math.pi,-math.pi,-math.pi,0.5,0.5])
HI=np.array(list(F.HI)+[ 12., 5., 8., 5., 8., math.pi, math.pi, math.pi, 2.0, 2.0])
LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
for i in BSLOT: LO[i]=0.0
ties={7:3, 27:11, 6:2, 12:10}
FREE=[i for i in range(37) if i not in set(ties)|{23,28,29,30,31,32,33}]
def expand(x):
    p=np.zeros(37); p[FREE]=x; p[23]=0.0
    for t,s in ties.items(): p[t]=p[s]
    return p
def resid(x):
    p=expand(x); out=F.blocks(p); r=[]
    for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
    r+=ha_res(p)
    for xx in (0.05,0.35,1.00):
        s=sl(p,xx)
        for n in BL: r.append(W*min(0.0,s[n]))
    return np.array(r)
P0=np.load("fitpar_production.npy")
c0=float(np.sum(resid(P0[FREE])**2))
print(f"production solution: chi2 = {c0:.2f}\n")
print(f"{'start':>5} {'chi2':>9} {'d_chi2':>8} {'ET_d slope':>11} {'phi_ET':>7} {'n/p':>6} {'|ETd/ETu|':>10}")
res=[]
for k in range(14):
    q=np.array(P0)
    for i in (0,4,8,13,15): q[i]*= rng.uniform(0.5,2.0)
    for i in BSLOT: q[i]=max(0.0,q[i]+rng.uniform(-0.5,0.5))
    for i in (2,6,10,12,21): q[i]+=rng.uniform(-0.4,0.4)
    for i in (3,11,17): q[i]+=rng.uniform(-0.6,0.6)
    for i in (34,35,36): q[i]=rng.uniform(-math.pi,math.pi)
    for t,s in ties.items(): q[t]=q[s]
    q=np.clip(q,LO,HI)
    try:
        r=least_squares(resid,q[FREE],bounds=(LO[FREE],HI[FREE]),x_scale='jac',
                        xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
    except Exception as e:
        print(f"{k:5d}  failed: {e}"); continue
    p=expand(r.x); c=float(np.sum(r.fun**2))
    s25=sl(p,0.25); HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    sp=amp.structure(p,"pi0p",-0.27,0.36,1.75); sn=amp.structure(p,"pi0n",-0.27,0.36,1.75)
    print(f"{k:5d} {c:9.2f} {c-c0:+8.2f} {s25['Ebar_T^d']:11.2f} {p[35]:+7.3f} {sn['TT']/sp['TT']:6.3f} {abs(ET[1])/ET[0]:10.3f}",flush=True)
    res.append((c,p))
res.sort(key=lambda z:z[0])
print(f"\nbest of {len(res)} starts: chi2 = {res[0][0]:.2f}  (production {c0:.2f})")
if res[0][0] < c0-0.5:
    np.save("fitpar_multistart_best.npy",res[0][1]); print("A BETTER minimum was found -> fitpar_multistart_best.npy")
else:
    print("no better minimum found; the production solution stands")
n_same=sum(1 for c,_ in res if abs(c-c0)<1.0)
print(f"{n_same} of {len(res)} starts returned to within 1.0 of it")
