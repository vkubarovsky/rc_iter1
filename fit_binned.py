"""Refit with the model AVERAGED OVER THE ACCEPTED BIN VOLUME.

Until now the model was evaluated at the published (mean) point of each bin while
the data are averages over the bin.  Near tmin that is a real bias: sigma_TT ~ t'
vanishes at threshold and moves by tens of percent across the first t bins, and
the tmin boundary itself sweeps across the bin as Q2 and xB vary inside it.  The
b2 t^2 curvature of Ebar_T was doing most of its work exactly there (85% of its
chi2 significance sits at -t' < 0.3), so the honest test of whether b2 is needed
is this fit, not the point-evaluated one.

Bins: PRC 90 025205 / PRC 95 035202 Tables I-III, accepted volume from bins.py.
The BSA and eg1-dvcs sectors stay point-evaluated: their binning is not published
in these two papers.

Variants:  free   27 par, no slope constraints (does the H_T^d pathology survive
                  proper binning at all?)
           cons   27 par, slope >= 0 and d >= u
           nob2   26 par, same constraints, b2 = 0
"""
import math, os, numpy as np
from scipy.optimize import least_squares
import amplitudes as amp
import bins
import fit_slope as F
from reparam import slopes

W=30.0; XS=(0.10,0.60)
NAMES=("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")

# ---- precompute the quadrature once ---------------------------------------
NODES=[]
for rows,ch,tag in ((F.PI0,"pi0p","pi0"),(F.ETA,"etap","eta")):
    for d in rows:
        iq,ix,it=bins.bin_of(d["Q2"],d["xB"],-d["t"])
        nd=bins.nodes(iq,ix,it,ch)
        eps=np.array([amp.epsilon(x,q,F.E_XS) for q,x,_,_ in nd])
        NODES.append((ch,nd,eps,d))
print(f"quadrature: {sum(len(n[1]) for n in NODES)} nodes over {len(NODES)} bins")

def binned_xsec(p):
    """chi2 residuals of the bin-averaged sigma_U, sigma_TT, sigma_LT."""
    out={"pi0":[], "eta":[]}
    for ch,nd,eps,d in NODES:
        sU=sTT=sLT=0.0
        for (Q2,xB,mt,w),e in zip(nd,eps):
            s=amp.structure(p,ch,-mt,xB,Q2)
            if s is None: continue
            sU+=w*(s["T"]+e*s["L"]); sTT+=w*s["TT"]; sLT+=w*s["LT"]
        tag="pi0" if ch=="pi0p" else "eta"
        out[tag]+= [(d["U"]-sU)/d["dU"], (d["TT"]-sTT)/d["dTT"], (d["LT"]-sLT)/d["dLT"]]
    return out

def pen(p, cons):
    if not cons: return []
    r=[]
    for x in XS:
        s=slopes(p,x)
        for n in NAMES: r.append(W*min(0.0,s[n]))
        r.append(W*min(0.0,s["H_T^d"]-s["H_T^u"]))
    return r

def blocks(p, cons):
    out=binned_xsec(p)
    o2=F.blocks(p)                        # BSA + eg1, point-evaluated
    for k in ("bsa_pi0","bsa_eta","bsa_c12","eg1"): out[k]=o2[k]
    return out

def make(variant):
    cons = variant!="free"
    free = [i for i in range(27) if not (variant=="nob2" and i==12)]
    def expand(x):
        p=np.zeros(27); p[free]=x
        if variant=="nob2": p[12]=0.0
        return p
    def resid(x):
        p=expand(x); out=blocks(p,cons); r=[]
        for k in ("pi0","eta","bsa_pi0","bsa_eta","bsa_c12","eg1"): r+=out[k]
        r.append((p[13]-0.54)/0.15); r+=pen(p,cons)
        return np.array(r)
    return free, expand, resid, cons

if __name__=="__main__":
    print(f"{'variant':>8} {'par':>4} {'chi2':>9} {'chi2/ndf':>9} {'pi0':>8} {'eta':>8} {'bsa':>8} {'eg1':>7}"
          f" {'b2':>7} {'R_ET':>7} {'pull':>6} {'sl_u':>6} {'sl_d':>6}")
    for variant in ("free","cons","nob2"):
        free,expand,resid,cons=make(variant)
        best=None
        for src in ("fitpar_prod.npy","fitpar_amp2026s.npy","fitpar_amp2026_lx.npy"):
            q=np.load(src)
            if variant=="nob2": q=np.array(q); q[12]=0.0
            try:
                r=least_squares(resid,np.clip(q[free],F.LO[free],F.HI[free]),
                                bounds=(F.LO[free],F.HI[free]),x_scale='jac',
                                xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=40000)
            except Exception as e:
                print("  seed fail",src,e); continue
            cc=float(np.sum(r.fun**2))
            if best is None or cc<best[0]: best=(cc,r.x)
        p=expand(best[1]); out=blocks(p,cons)
        tot=sum(sum(x*x for x in v) for v in out.values()); n=sum(len(v) for v in out.values())
        npar=len(free); s25=slopes(p,0.25)
        bsa=sum(sum(x*x for x in out[k]) for k in ("bsa_pi0","bsa_eta","bsa_c12"))
        print(f"{variant:>8} {npar:4d} {tot:9.1f} {tot/(n-npar):9.4f} {sum(x*x for x in out['pi0']):8.1f}"
              f" {sum(x*x for x in out['eta']):8.1f} {bsa:8.1f} {sum(x*x for x in out['eg1']):7.1f}"
              f" {p[12]:7.3f} {p[13]:7.3f} {(p[13]-0.54)/0.15:+6.1f} {s25['H_T^u']:6.2f} {s25['H_T^d']:6.2f}")
        np.save(f"fitpar_binned_{variant}.npy",p)
