"""Would a new RC iteration move anything?

r = eta(model) / eta(published) is what corrects the data.  Our current strfun
files carry r computed with amp2021 (= the iteration-1 amplitude model).  The
question is the DOUBLE ratio

    R = eta(today's model) / eta(amp2021)

which is exactly the factor by which the corrected cross sections would move if
the iteration were repeated with today's model.  R = 1 means nothing to redo.

For scale: switching from the sigma_L-dominant vpk2013 to the sigma_T-dominant
amp2021 -- a change of regime -- gave r(pi0) = 1.019 and r(eta) = 1.007.
"""
import math, sys, numpy as np
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,"/Users/vpk/exclurad_py")
from exclurad_py import models
from exclurad_py.core.rc import rc_factor
from exclurad_py.core.constants import M_P
CH=sys.argv[1] if len(sys.argv)>1 else "pi0"
CFG=dict(pi0=dict(old="pi0.amp2021", new="pi0.amp2026n", vcut=0.18, m=0.1349768),
         eta=dict(old="eta.amp2021", new="eta.amp2026n", vcut=0.094, m=0.547862))[CH]
mx2=CFG["m"]**2+CFG["vcut"]
PTS=[]
for Q2,xB in ((1.15,0.13),(1.38,0.17),(1.75,0.22),(2.21,0.28),(2.71,0.34),(3.22,0.43)):
    for mt in (0.20,0.40,0.80):
        for ph in (0.,90.,180.):
            PTS.append((Q2,xB,mt,ph))
def one(i):
    Q2,xB,mt,ph=PTS[i]
    W2=M_P*M_P+Q2*(1-xB)/xB
    out=[]
    for key in (CFG["old"],CFG["new"]):
        try:
            m=models.get(key,t_nucl=-mt,Ebeam=5.75)
            r=rc_factor(5.75,Q2,W2,-mt,math.radians(ph),m,h=+1,mx2_cut=mx2)
            out.append(float(r["eta"]))
        except Exception as e:
            out.append(float('nan'))
    return (Q2,xB,mt,ph,out[0],out[1])
if __name__=="__main__":
    with ProcessPoolExecutor(max_workers=10) as ex:
        res=list(ex.map(one,range(len(PTS))))
    R=np.array([r[5]/r[4] for r in res if np.isfinite(r[4]) and np.isfinite(r[5]) and r[4]>0])
    print(f"{CH}: {len(R)} of {len(PTS)} points computed\n")
    print(f"{'Q2':>5} {'xB':>5} {'-t':>5} {'phi':>5} | {'eta(amp2021)':>13} {'eta(today)':>11} {'ratio':>8}")
    for r in res[:12]:
        if not (np.isfinite(r[4]) and np.isfinite(r[5])): continue
        print(f"{r[0]:5.2f} {r[1]:5.2f} {r[2]:5.2f} {r[3]:5.0f} | {r[4]:13.4f} {r[5]:11.4f} {r[5]/r[4]:8.4f}")
    print(f"\ndouble ratio R = eta(today)/eta(amp2021):")
    print(f"   median {np.median(R):.4f}   mean {R.mean():.4f}   16-84% [{np.percentile(R,16):.4f},{np.percentile(R,84):.4f}]")
    print(f"   min {R.min():.4f}   max {R.max():.4f}")
    print(f"   => the corrected cross sections would move by {abs(np.median(R)-1)*100:.2f}% (median), "
          f"{max(abs(R.min()-1),abs(R.max()-1))*100:.2f}% at worst")
