"""GK_2024 (hepgen, the real code) against the CLAS6 data and against our fit.

Uses the precomputed convolution tables that already sit in hepgen's
preparation/ directory (exact tmin, exact xi, vpk mixing angle, mu_eta = 1.76),
so nothing has to be recomputed.  Work is done in ~/gk_work, a copy: the
original ~/hepgen_mac is not touched.
"""
import os, re, glob, math, subprocess, json
import numpy as np
import amplitudes as amp, fit_slope as F

GK="/Users/vpk/gk_work"
BIN=f"{GK}/build/bin/get_SF"
PAT=f"{GK}/preparation/model_GK_2024_CLAS6_y12_*_Ebeam_5.75_mu_eta_1.76_tmin_exact_xi_exact_mixangle_vpk.dat"
OUR=np.load("fitpar_n2_final.npy")

tables=[]
for f in sorted(glob.glob(PAT)):
    m=re.search(r"y12_(pi0|eta)_(p|n)_Q2_([0-9.]+)_xB_([0-9.]+)_",os.path.basename(f))
    if m: tables.append(dict(ch=m.group(1),tgt=m.group(2),Q2=float(m.group(3)),xB=float(m.group(4)),f=f))
print(f"{len(tables)} GK_2024 tables\n",flush=True)

def gk(tab,mt):
    cmd=[BIN,"-pi0" if tab["ch"]=="pi0" else "-eta","-proton","-Q2",f"{tab['Q2']}","-xB",f"{tab['xB']}",
         "-t",f"{mt}","-phi","0","-Ebeam","5.75","-model","GK_2024","-prep",tab["f"],"--fmt=csv"]
    try:
        out=subprocess.run(cmd,capture_output=True,text=True,timeout=300).stdout.strip().splitlines()[-1]
        s0,sT,sL,sTT,sLT=[float(x) for x in out.split(",")]
        return dict(T=sT,L=sL,TT=sTT,LT=sLT)
    except Exception as e:
        return None

rows=[]
for tab in tables:
    data=F.PI0 if tab["ch"]=="pi0" else F.ETA
    pts=[d for d in data if abs(d["Q2"]-tab["Q2"])<0.06 and abs(d["xB"]-tab["xB"])<0.012]
    if not pts: continue
    for d in pts:
        g=gk(tab,-d["t"])
        if g is None: continue
        e=amp.epsilon(d["xB"],d["Q2"],F.E_XS)
        s=amp.structure(OUR,"pi0p" if tab["ch"]=="pi0" else "etap",d["t"],d["xB"],d["Q2"])
        rows.append(dict(ch=tab["ch"],Q2=d["Q2"],xB=d["xB"],mt=-d["t"],
                         dU=d["U"],edU=d["dU"],dTT=d["TT"],edTT=d["dTT"],dLT=d["LT"],edLT=d["dLT"],
                         gU=g["T"]+e*g["L"],gTT=g["TT"],gLT=g["LT"],gL=g["L"],gT=g["T"],
                         oU=s["T"]+e*s["L"],oTT=s["TT"],oLT=s["LT"],oL=s["L"],oT=s["T"]))
    print(f"  {tab['ch']} Q2={tab['Q2']} xB={tab['xB']}: {len(pts)} points",flush=True)
json.dump(rows,open("gk2024_rows.json","w"))
print(f"\n{len(rows)} matched points\n")
def chi(rows,key,mkey):
    return sum(((r[key]-r[mkey])/r["e"+key])**2 for r in rows), len(rows)
for ch in ("pi0","eta"):
    rr=[r for r in rows if r["ch"]==ch]
    if not rr: continue
    print(f"--- {ch}: {len(rr)} bins ---")
    print(f"{'observable':>10} {'chi2 GK_2024':>14} {'chi2 our fit':>14}")
    for key,g,o in (("dU","gU","oU"),("dTT","gTT","oTT"),("dLT","gLT","oLT")):
        cg=sum(((r[key]-r[g])/r["e"+key])**2 for r in rr)
        co=sum(((r[key]-r[o])/r["e"+key])**2 for r in rr)
        print(f"{key[1:]:>10} {cg:14.1f} {co:14.1f}")
    print(f"{'sigma_L/sigma_T at the first bin':>10}: GK {rr[0]['gL']/rr[0]['gT']:.3f}, ours {rr[0]['oL']/rr[0]['oT']:.3f}\n")
