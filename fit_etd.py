"""Ebar_T^d as an independent block, same form as u (VPK's proposal).

    ET_u = N_u exp[(b_u + b'_u ln xB) t] Q^nQ_u          (b2 removed)
    ET_d = N_d exp[(b_d + b'_d ln xB) t] Q^nQ_d

28-parameter vector: the old R_ET slot p[13] is now N_d, db slot p[14] is b_d,
the freed b2 slot p[12] is b'_d, and p[27] is nQ_d.

    free    all four free
    tieQ    nQ_d = nQ_u
    nobx    b'_d = 0  (no ln xB term for d)
    tieQ+nobx
No R_ET prior anywhere.  Only constraint: every slope b + b' ln xB >= 0 on
xB in [0.1, 0.6].
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F
W=30.0
BL={"H_T^u":(1,2,None),"H_T^d":(5,6,None),"Ebar_T^u":(9,10,None),"Ebar_T^d":(14,12,None),"T00":(16,21,None)}
def slopes28(p,x):
    L=math.log(x)
    return {k:(p[b]+p[bp]*L) for k,(b,bp,_) in BL.items()}
LO=np.array(list(F.LO)+[-12.0]); HI=np.array(list(F.HI)+[12.0])
LO[3]=LO[7]=LO[11]=LO[17]=-12.0; HI[3]=HI[7]=HI[11]=HI[17]=12.0
LO[12]=-8.0; HI[12]=8.0            # b'_d
LO[13]=-1e4; HI[13]=1e4            # N_d is now an absolute normalisation
LO[14]=-2.0; HI[14]=12.0           # b_d
def seed():
    q=np.array(np.load("fitpar_plain.npy"))
    p=np.zeros(28); p[:27]=q
    p[13]=q[13]*q[8]      # N_d  = R_ET * N_u
    p[14]=q[9]+q[14]      # b_d  = b_u + db
    p[12]=q[10]           # b'_d = b'_u
    p[27]=q[11]           # nQ_d = nQ_u
    return p
def run(tag, tieQ, nobx):
    free=[i for i in range(28)]
    if tieQ: free.remove(27)
    if nobx: free.remove(12)
    def expand(x):
        p=np.zeros(28); p[free]=x
        if tieQ: p[27]=p[11]
        if nobx: p[12]=0.0
        return p
    def resid(x):
        p=expand(x); out=F.blocks(p); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        for xx in (0.10,0.60):
            s=slopes28(p,xx)
            for n in s: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for scale in (1.0,0.3,3.0):
        q=seed(); q[13]*=scale; q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=60000)
        except Exception as e: print(" seed fail",e); continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); out=F.blocks(p)
    tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
    sp=amp.structure(p,"pi0p",-0.27,0.36,1.75); sn=amp.structure(p,"pi0n",-0.27,0.36,1.75)
    s25=slopes28(p,0.25)
    print(f"{tag:>10} {len(free):4d} {tot:8.1f} {tot/(n-len(free)):8.4f}  N_d/N_u {p[13]/p[8]:+7.3f}  "
          f"b_d {p[14]:6.2f} b'_d {p[12]:+6.2f} nQ_d {p[27]:+6.2f}  slope_d(0.25) {s25['Ebar_T^d']:5.2f} "
          f"(u {s25['Ebar_T^u']:4.2f})  n/p {sn['TT']/sp['TT']:5.2f}  "
          f"pi0 {sum(x*x for x in out['pi0']):6.1f} eta {sum(x*x for x in out['eta']):6.1f} "
          f"bsa {sum(sum(x*x for x in out[k]) for k in ('bsa_pi0','bsa_eta','bsa_c12')):6.1f} "
          f"eg1 {sum(x*x for x in out['eg1']):5.1f}",flush=True)
    np.save(f"fitpar_etd_{tag}.npy",p)
print("Hall-A neutron ratio for reference: 0.28 +- 0.07\n")
print(f"{'variant':>10} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("free",False,False)
run("tieQ",True,False)
run("nobx",False,True)
run("tieQ_nobx",True,True)
