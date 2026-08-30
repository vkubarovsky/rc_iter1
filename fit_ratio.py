"""VPK's proposal: put eta into the fit through the RATIO eta/pi0, with
statistical errors only.

Both channels come from the same e1-dvcs analysis: same beam, same luminosity,
same electron and proton cuts, same acceptance code.  Those systematics are
common and cancel in R = sigma_U(eta)/sigma_U(pi0) taken at the SAME bin.  What
does not cancel is meson-specific: the M(gamma gamma) window and the background
subtraction (25% under the eta peak against 3-5% under the pi0).

This is a REPLACEMENT, not an addition: for every matched bin the eta sigma_U is
removed from the fit and R enters instead, so no information is counted twice.
All 76 eta bins have a pi0 partner when matched by bin index.

    ratio_stat   sigma_R from statistical errors only (VPK's version)
    ratio_syst   sigma_R with the eta-specific systematic added, as a check
    absolute     the current fit, for comparison
"""
import math, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, bins, fit_slope as F
W=30.0
BL={"H_T^u":(1,2),"H_T^d":(5,6),"Ebar_T^u":(9,10),"Ebar_T^d":(14,12),"T00":(16,21)}
BSLOT=[1,5,9,14,16]
def slopes32(p,x):
    L=math.log(x); return {k:(p[b]+p[bp]*L) for k,(b,bp) in BL.items()}
def load(fn):
    rows=[]
    for line in open(fn):
        v=line.split()
        if len(v)<12: continue
        v=[float(x) for x in v[:12]]
        rows.append(dict(Q2=v[0],xB=v[1],t=-v[2],U=v[3],sU=v[4],yU=v[5],
                         LT=v[6],sLT=v[7],yLT=v[8],TT=v[9],sTT=v[10],yTT=v[11]))
    return rows
PI0=load("data/strfun_pi0.data"); ETA=load("data/strfun_eta.data")
bp={}; be={}
for d in PI0:
    k=bins.bin_of(d["Q2"],d["xB"],-d["t"])
    if None not in k: bp[k]=d
for d in ETA:
    k=bins.bin_of(d["Q2"],d["xB"],-d["t"])
    if None not in k: be[k]=d
PAIRS=[(bp[k],be[k]) for k in sorted(set(bp)&set(be))]
print(f"{len(PAIRS)} matched (pi0, eta) bins\n")
HA=dict(Q2=1.75,xB=0.36,mt=0.27)
def nratio(p):
    sp=amp.structure(p,"pi0p",-HA["mt"],HA["xB"],HA["Q2"]); sn=amp.structure(p,"pi0n",-HA["mt"],HA["xB"],HA["Q2"])
    return sn["TT"]/sp["TT"] if sp and sn else float('nan')
def run(tag, mode, sector=True):
    LO=np.array(list(F.LO)+[-12.,-5.,-5.,-5.,-5.]); HI=np.array(list(F.HI)+[12.,5.,8.,5.,8.])
    LO[3]=LO[11]=LO[17]=-12.; HI[3]=HI[11]=HI[17]=12.
    LO[12]=-8.; HI[12]=8.; LO[13]=-1e4; HI[13]=1e4; LO[14]=-2.; HI[14]=12.
    for i in BSLOT: LO[i]=0.0
    ties={7:3, 27:11}
    if sector: ties.update({6:2, 12:10})
    free=[i for i in range(32) if i not in set(ties)|{23,28,29,30,31}]
    def expand(x):
        p=np.zeros(32); p[free]=x; p[23]=0.0
        for t,s in ties.items(): p[t]=p[s]
        return p
    def parts(p):
        out={"pi0":[], "eta":[], "R":[]}
        for d in PI0:
            s=amp.structure(p,"pi0p",d["t"],d["xB"],d["Q2"]); e=amp.epsilon(d["xB"],d["Q2"],F.E_XS)
            eU=math.hypot(d["sU"],d["yU"]); eTT=math.hypot(d["sTT"],d["yTT"]); eLT=math.hypot(d["sLT"],d["yLT"])
            out["pi0"]+=[(d["U"]-(s["T"]+e*s["L"]))/eU,(d["TT"]-s["TT"])/eTT,(d["LT"]-s["LT"])/eLT]
        for d in ETA:
            s=amp.structure(p,"etap",d["t"],d["xB"],d["Q2"]); e=amp.epsilon(d["xB"],d["Q2"],F.E_XS)
            eTT=math.hypot(d["sTT"],d["yTT"]); eLT=math.hypot(d["sLT"],d["yLT"])
            out["eta"]+=[(d["TT"]-s["TT"])/eTT,(d["LT"]-s["LT"])/eLT]
            if mode=="absolute":
                out["eta"].append((d["U"]-(s["T"]+e*s["L"]))/math.hypot(d["sU"],d["yU"]))
        if mode!="absolute":
            for dp,de in PAIRS:
                R=de["U"]/dp["U"]
                rel=(de["sU"]/de["U"])**2+(dp["sU"]/dp["U"])**2
                if mode=="ratio_syst": rel+=(0.09)**2          # eta-specific: M(gg) window + background
                sR=R*math.sqrt(rel)
                ep=amp.epsilon(dp["xB"],dp["Q2"],F.E_XS); ee=amp.epsilon(de["xB"],de["Q2"],F.E_XS)
                sp_=amp.structure(p,"pi0p",dp["t"],dp["xB"],dp["Q2"]); se=amp.structure(p,"etap",de["t"],de["xB"],de["Q2"])
                Rm=(se["T"]+ee*se["L"])/(sp_["T"]+ep*sp_["L"])
                out["R"].append((R-Rm)/sR)
        return out
    def resid(x):
        p=expand(x); o=parts(p); r=o["pi0"]+o["eta"]+o["R"]
        o2=F.blocks(p)
        for k in ("bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=o2[k]
        for xx in (0.05,0.35,1.00):
            s=slopes32(p,xx)
            for n in BL: r.append(W*min(0.0,s[n]))
        return np.array(r)
    best=None
    for src in ("fitpar_as_sector.npy","fitpar_production.npy","fitpar_tieQx_data.npy"):
        q=np.zeros(32); z=np.load(src); q[:len(z)]=z
        for t,s in ties.items(): q[t]=q[s]
        q=np.clip(q,LO,HI)
        try: r=least_squares(resid,q[free],bounds=(LO[free],HI[free]),x_scale='jac',
                             xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=60000)
        except Exception: continue
        c=float(np.sum(r.fun**2))
        if best is None or c<best[0]: best=(c,r.x)
    p=expand(best[1]); o=parts(p); o2=F.blocks(p)
    cp=sum(x*x for x in o["pi0"]); ce=sum(x*x for x in o["eta"]); cR=sum(x*x for x in o["R"])
    cb=sum(sum(x*x for x in o2[k]) for k in ("bsa_pi0","bsa_eta","bsa_c12")); cg=sum(x*x for x in o2["eg1"])
    n=len(o["pi0"])+len(o["eta"])+len(o["R"])+sum(len(o2[k]) for k in ("bsa_pi0","bsa_eta","bsa_c12","eg1"))
    s25=slopes32(p,0.25); HT,ET,LL=amp._flavour(p,-0.3,0.25,2.2)
    print(f"{tag:>11} {len(free):4d} {cp+ce+cR+cb+cg:8.1f} {(cp+ce+cR+cb+cg)/(n-len(free)):8.4f} | "
          f"pi0 {cp:6.1f}/{len(o['pi0'])} eta {ce:6.1f}/{len(o['eta'])} R {cR:6.1f}/{len(o['R'])} "
          f"bsa {cb:6.1f} eg1 {cg:5.1f} | slope ET u/d {s25['Ebar_T^u']:4.2f}/{s25['Ebar_T^d']:5.2f} "
          f"ETd/u {ET[1]/ET[0]:+5.2f} n/p {nratio(p):5.2f}",flush=True)
    np.save(f"fitpar_R_{tag}.npy",p); return p
print(f"{'variant':>11} {'par':>4} {'chi2':>8} {'chi2/ndf':>8}")
run("absolute","absolute")
run("ratio_stat","ratio_stat")
run("ratio_syst","ratio_syst")
