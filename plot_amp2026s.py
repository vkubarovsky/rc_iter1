"""All comparison figures for amp2026s (26 par, b2 removed, bin-averaged fit).

Solid curve  = amp2026s, the new production model, evaluated pointwise.
Dashed curve = amp2026, the previous 27-par model (b2 kept, slopes unconstrained).
Open squares = amp2026s AVERAGED OVER THE ACCEPTED BIN VOLUME - that is what the
fit actually compares to the data, and it differs from the curve near tmin.
"""
import json, math, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as amp, bins
import fit_slope as F
from reparam import slopes

os.makedirs("figures", exist_ok=True)
NEW=np.load("fitpar_amp2026s.npy"); OLD=np.load("fitpar_amp2026_lx.npy")
E_XS,E_BSA,E_EG1,E_C12=F.E_XS,F.E_BSA,F.E_EG1,F.E_C12
CN,CO='#1f4e9c','#c23b22'

def groups(rows):
    gs=[]
    for d in rows:
        for g in gs:
            if abs(g['Q2']-d['Q2'])<0.07 and abs(g['xB']-d['xB'])<0.012: g['pts'].append(d); break
        else: gs.append(dict(Q2=d['Q2'],xB=d['xB'],pts=[d]))
    for g in gs:
        g['Q2']=np.mean([x['Q2'] for x in g['pts']]); g['xB']=np.mean([x['xB'] for x in g['pts']])
    return sorted(gs,key=lambda g:(round(g['xB'],2),g['Q2']))

def binavg(p,ch,d):
    iq,ix,it=bins.bin_of(d["Q2"],d["xB"],-d["t"])
    nd=bins.nodes(iq,ix,it,ch)
    if nd is None: return None
    U=TT=LT=0.
    for Q2,xB,mt,w in nd:
        s=amp.structure(p,ch,-mt,xB,Q2)
        if s is None: continue
        U+=w*(s["T"]+amp.epsilon(xB,Q2,E_XS)*s["L"]); TT+=w*s["TT"]; LT+=w*s["LT"]
    return U,TT,LT

# ---------- 1-2: cross sections, all bins ----------------------------------
for rows,ch,meson,fname in ((F.PI0,"pi0p",r"\pi^0","xsec_pi0"),(F.ETA,"etap",r"\eta","xsec_eta")):
    gs=groups(rows); n=len(gs); nc=5; nr=math.ceil(n/nc)
    fig,axs=plt.subplots(nr,nc,figsize=(3.3*nc,2.9*nr),squeeze=False)
    for i,g in enumerate(gs):
        a=axs[i//nc][i%nc]; Q2,xB=g['Q2'],g['xB']; e=amp.epsilon(xB,Q2,E_XS)
        T=np.linspace(0.05,2.0,200)
        for p,ls,lw in ((OLD,'--',1.3),(NEW,'-',2.0)):
            U=[];TT=[];LT=[]
            for t in T:
                s=amp.structure(p,ch,-t,xB,Q2)
                if s is None: U.append(np.nan);TT.append(np.nan);LT.append(np.nan)
                else: U.append(s["T"]+e*s["L"]);TT.append(s["TT"]);LT.append(s["LT"])
            a.plot(T,U,ls,color='k',lw=lw); a.plot(T,TT,ls,color=CN,lw=lw); a.plot(T,LT,ls,color='seagreen',lw=lw)
        mt=[-d['t'] for d in g['pts']]
        a.errorbar(mt,[d['U'] for d in g['pts']],yerr=[d['dU'] for d in g['pts']],fmt='ko',ms=4,capsize=2)
        a.errorbar(mt,[d['TT'] for d in g['pts']],yerr=[d['dTT'] for d in g['pts']],fmt='^',color=CN,ms=4,capsize=2)
        a.errorbar(mt,[d['LT'] for d in g['pts']],yerr=[d['dLT'] for d in g['pts']],fmt='s',color='seagreen',ms=4,capsize=2)
        for d in g['pts']:
            v=binavg(NEW,ch,d)
            if v is None: continue
            a.plot(-d['t'],v[0],'s',mfc='none',mec='k',ms=8,mew=1.2)
            a.plot(-d['t'],v[1],'s',mfc='none',mec=CN,ms=8,mew=1.2)
            a.plot(-d['t'],v[2],'s',mfc='none',mec='seagreen',ms=8,mew=1.2)
        a.axhline(0,color='0.6',lw=0.7); a.set_xlim(0,2.0)
        a.set_title(f"$Q^2$={Q2:.2f}, $x_B$={xB:.3f}",fontsize=9.5); a.tick_params(labelsize=8); a.grid(alpha=0.2)
        if i//nc==nr-1: a.set_xlabel(r"$-t$ [GeV$^2$]",fontsize=9)
        if i%nc==0: a.set_ylabel(r"nb/GeV$^2$",fontsize=9)
    for j in range(n,nr*nc): axs[j//nc][j%nc].axis('off')
    h=[plt.Line2D([],[],color='k',lw=2,label=r'$\sigma_U$'),
       plt.Line2D([],[],color=CN,lw=2,label=r'$\sigma_{TT}$'),
       plt.Line2D([],[],color='seagreen',lw=2,label=r'$\sigma_{LT}$'),
       plt.Line2D([],[],color='k',lw=2,ls='-',label='amp2026s (26 par, production)'),
       plt.Line2D([],[],color='k',lw=1.3,ls='--',label='amp2026 (27 par, previous)'),
       plt.Line2D([],[],color='k',marker='s',mfc='none',ls='none',ms=8,label='amp2026s averaged over the bin')]
    fig.legend(handles=h,loc='lower right',fontsize=11,ncol=3,frameon=False)
    fig.suptitle(rf'CLAS6 ${meson}$ structure functions (RC-iterated data) vs the amplitude fits',fontsize=14)
    fig.tight_layout(rect=[0,0.03,1,0.965]); fig.savefig(f"figures/{fname}.png",dpi=115); plt.close(fig)
    print(fname,"done",flush=True)

# ---------- 3: CLAS6 pi0 BSA ------------------------------------------------
BSA=json.load(open("data/clas6_demasi_alu.json"))
n=len(BSA); nc=5; nr=math.ceil(n/nc)
fig,axs=plt.subplots(nr,nc,figsize=(3.3*nc,2.9*nr),squeeze=False)
for i,b in enumerate(BSA):
    a=axs[i//nc][i%nc]; xB,Q2=b['xB'],b['Q2']
    a.errorbar([q['t'] for q in b['pts']],[q['alpha'] for q in b['pts']],
               yerr=[q['err'] for q in b['pts']],fmt='ks',ms=4.5,capsize=2.5)
    T=np.linspace(0.06,1.7,200)
    for p,ls,lw,c in ((OLD,'--',1.3,CO),(NEW,'-',2.0,CN)):
        a.plot(T,[amp.bsa_sinphi(p,"pi0p",-t,xB,Q2,E_BSA) or np.nan for t in T],ls,color=c,lw=lw)
    a.axhline(0,color='0.6',lw=0.7); a.set_xlim(0,1.7); a.set_ylim(-0.05,0.22)
    a.set_title(f"$Q^2$={Q2:.2f}, $x_B$={xB:.2f}",fontsize=9.5); a.tick_params(labelsize=8); a.grid(alpha=0.2)
    if i//nc==nr-1: a.set_xlabel(r"$-t$ [GeV$^2$]",fontsize=9)
    if i%nc==0: a.set_ylabel(r"$A_{LU}^{\sin\phi}$",fontsize=10)
for j in range(n,nr*nc): axs[j//nc][j%nc].axis('off')
fig.legend(handles=[plt.Line2D([],[],color=CN,lw=2,label='amp2026s'),
                    plt.Line2D([],[],color=CO,lw=1.3,ls='--',label='amp2026')],
           loc='lower right',fontsize=12,frameon=False)
fig.suptitle(r'CLAS6 $\pi^0$ beam-spin asymmetry (De Masi, PRC 77 042201), all 13 bins',fontsize=14)
fig.tight_layout(rect=[0,0.02,1,0.955]); fig.savefig("figures/bsa_pi0_clas6.png",dpi=115); plt.close(fig)
print("bsa clas6 done",flush=True)

# ---------- 4: CLAS12 pi0 BSA + CLAS6 eta BSA -------------------------------
C12=json.load(open("data/clas12_kim_alu.json"))
ETB=json.load(open("data/clas6_zhao_eta_alu.json"))
ebins=sorted({(d['xB'],d['Q2']) for d in ETB})
fig,axs=plt.subplots(2,5,figsize=(17,6.6),squeeze=False)
for i,r in enumerate(C12):
    a=axs[0][i]; xB,Q2=r['xB'],r['Q2']
    a.errorbar([q['t'] for q in r['pts']],[q['A'] for q in r['pts']],yerr=[q['dA'] for q in r['pts']],
               fmt='ko',ms=5,capsize=2.5)
    T=np.linspace(0.15,1.8,200)
    for p,ls,lw,c in ((OLD,'--',1.3,CO),(NEW,'-',2.0,CN)):
        a.plot(T,[amp.bsa_sinphi(p,"pi0p",-t,xB,Q2,E_C12) or np.nan for t in T],ls,color=c,lw=lw)
    a.axhline(0,color='0.6',lw=0.7); a.set_xlim(0,1.8); a.set_ylim(-0.05,0.25); a.grid(alpha=0.2)
    a.set_title(rf"CLAS12 $\pi^0$: $Q^2$={Q2:.2f}, $x_B$={xB:.2f}",fontsize=10)
    a.set_xlabel(r"$-t$ [GeV$^2$]",fontsize=9)
    if i==0: a.set_ylabel(r"$A_{LU}^{\sin\phi}$",fontsize=10)
for i,(xB,Q2) in enumerate(ebins):
    a=axs[1][i]
    pts=[d for d in ETB if (d['xB'],d['Q2'])==(xB,Q2)]
    a.errorbar([d['t'] for d in pts],[d['alpha'] for d in pts],
               yerr=[math.hypot(d['stat'],d['syst']) for d in pts],fmt='k^',ms=6,capsize=2.5)
    T=np.linspace(0.1,1.5,200)
    for p,ls,lw,c in ((OLD,'--',1.3,CO),(NEW,'-',2.0,CN)):
        a.plot(T,[amp.bsa_sinphi(p,"etap",-t,xB,Q2,E_BSA) or np.nan for t in T],ls,color=c,lw=lw)
    a.axhline(0,color='0.6',lw=0.7); a.set_xlim(0,1.5); a.set_ylim(-0.1,0.45); a.grid(alpha=0.2)
    a.set_title(rf"CLAS6 $\eta$: $Q^2$={Q2:.2f}, $x_B$={xB:.2f}",fontsize=10)
    a.set_xlabel(r"$-t$ [GeV$^2$]",fontsize=9)
    if i==0: a.set_ylabel(r"$A_{LU}^{\sin\phi}$",fontsize=10)
for j in range(len(ebins),5): axs[1][j].axis('off')
fig.legend(handles=[plt.Line2D([],[],color=CN,lw=2,label='amp2026s'),
                    plt.Line2D([],[],color=CO,lw=1.3,ls='--',label='amp2026')],
           loc='lower right',fontsize=12,frameon=False)
fig.suptitle(r'Beam-spin asymmetries: CLAS12 $\pi^0$ (Kim, PLB 849 138459) and CLAS6 $\eta$ (Zhao, PLB 789 426)',fontsize=13)
fig.tight_layout(rect=[0,0.02,1,0.94]); fig.savefig("figures/bsa_clas12_eta.png",dpi=115); plt.close(fig)
print("bsa clas12+eta done",flush=True)

# ---------- 5: eg1-dvcs target/double-spin moments --------------------------
EG1=json.load(open("data/eg1dvcs_pi0_target_asym.json"))
KIN={"1.94":0.25,"2.83":0.40}
SPEC=[("E154M5","E154M6","AULsin",r"$A_{UL}^{\sin\phi}$"),
      ("E154M7","E154M8","AULsin2",r"$A_{UL}^{\sin2\phi}$"),
      ("E154M9","E154M10","ALLc",r"$A_{LL}^{\rm const}$"),
      ("E154M11","E154M12","ALLcos",r"$A_{LL}^{\cos\phi}$")]
def moment(p,key,Q2,xB,t):
    s=amp.structure(p,"pi0p",-t,xB,Q2)
    if s is None: return np.nan
    e=amp.epsilon(xB,Q2,E_EG1); s0=s["T"]+e*s["L"]
    return dict(AULsin=math.sqrt(2*e*(1+e))*s["LTUL"]/s0, AULsin2=e*s["TTUL"]/s0,
                ALLc=math.sqrt(1-e*e)*s["Tp"]/s0, ALLcos=math.sqrt(2*e*(1-e))*s["LTpLL"]/s0)[key]
fig,axs=plt.subplots(2,4,figsize=(16,7.2),squeeze=False)
for col,(r1,r2,key,lab) in enumerate(SPEC):
    for row,rec in enumerate((r1,r2)):
        a=axs[row][col]; rows=EG1[rec]['rows']; Q2=rows[0][0]; xB=KIN[f"{Q2:.2f}"]
        a.errorbar([r[1] for r in rows],[r[2] for r in rows],
                   yerr=[math.hypot(r[3],r[4] if len(r)>4 else 0) for r in rows],
                   fmt='ko',ms=6,capsize=3)
        T=np.linspace(0.06,1.7,200)
        for p,ls,lw,c in ((OLD,'--',1.3,CO),(NEW,'-',2.0,CN)):
            a.plot(T,[moment(p,key,Q2,xB,t) for t in T],ls,color=c,lw=lw)
        a.axhline(0,color='k',lw=0.8); a.set_xlim(0,1.75); a.grid(alpha=0.25)
        a.set_title(rf"{lab}   $Q^2$={Q2}, $x_B$={xB}",fontsize=11)
        if row==1: a.set_xlabel(r'$-t$ [GeV$^2$]',fontsize=11)
fig.legend(handles=[plt.Line2D([],[],color=CN,lw=2,label='amp2026s'),
                    plt.Line2D([],[],color=CO,lw=1.3,ls='--',label='amp2026')],
           loc='lower right',fontsize=12,frameon=False)
fig.suptitle('eg1-dvcs longitudinal-target and double-spin moments (in the fit)',fontsize=13.5)
fig.tight_layout(rect=[0,0.02,1,0.94]); fig.savefig("figures/eg1_targets.png",dpi=120); plt.close(fig)
print("eg1 done",flush=True)

# ---------- 6: form factors, slopes, sigma_L/sigma_T -------------------------
fig,axs=plt.subplots(2,3,figsize=(15,8),squeeze=False)
T=np.linspace(0.02,2.4,200); Q2r,xBr=2.2,0.25
for j,(name,idx) in enumerate((("H_T",0),("\\bar E_T",1),("T_{00}",2))):
    a=axs[0][j]
    for p,ls,lw,c in ((OLD,'--',1.3,CO),(NEW,'-',2.0,CN)):
        u=[];d_=[]
        for t in T:
            HT,ET,LL=amp._flavour(p,-t,xBr,Q2r)
            pair=(HT,ET,LL)[idx]; u.append(pair[0]); d_.append(pair[1])
        a.plot(T,np.abs(u),ls,color=c,lw=lw)
        a.plot(T,np.abs(d_),ls,color=c,lw=lw,alpha=0.45)
    a.set_yscale('log'); a.set_xlabel(r'$-t$ [GeV$^2$]'); a.grid(alpha=0.25)
    a.set_title(rf'$|{name}^{{u}}|$ (dark) and $|{name}^{{d}}|$ (pale),  $Q^2$=2.2, $x_B$=0.25',fontsize=10.5)
a=axs[1][0]
X=np.linspace(0.10,0.60,80)
for nm,c in (("H_T^u",'k'),("H_T^d",'crimson'),("Ebar_T^u",'royalblue'),("Ebar_T^d",'navy'),("T00",'seagreen')):
    a.plot(X,[slopes(NEW,x)[nm] for x in X],'-',color=c,lw=2,label=nm.replace("Ebar",r"$\bar E$"))
    a.plot(X,[slopes(OLD,x)[nm] for x in X],'--',color=c,lw=1.2)
a.axhline(0,color='k',lw=1); a.set_xlabel(r'$x_B$'); a.set_ylabel(r'slope $b+b^\prime\ln x_B$ [GeV$^{-2}$]')
a.set_title('t-slopes: solid amp2026s, dashed amp2026\n(amp2026 $H_T^d$ is NEGATIVE everywhere)',fontsize=10.5)
a.legend(fontsize=8,ncol=2); a.grid(alpha=0.25)
a=axs[1][1]
for p,ls,lw,c,lab in ((OLD,'--',1.3,CO,'amp2026'),(NEW,'-',2.0,CN,'amp2026s')):
    def rat(t,p=p):
        s_=amp.structure(p,"pi0p",-t,0.25,2.2)
        return np.nan if s_ is None else s_["L"]/s_["T"]
    a.plot(T,[rat(t) for t in T],ls,color=c,lw=lw,label=lab)
a.set_xlabel(r'$-t$ [GeV$^2$]'); a.set_ylabel(r'$\sigma_L/\sigma_T$'); a.grid(alpha=0.25); a.legend(fontsize=9)
a.set_title(r'$\sigma_L/\sigma_T$ for $\pi^0$, $Q^2$=2.2, $x_B$=0.25',fontsize=10.5)
a=axs[1][2]
for ch,c,lab in (("pi0p",'k',r'$\pi^0 p$'),("etap",'royalblue',r'$\eta p$'),("pi0n",'crimson',r'$\pi^0 n$')):
    def rtt(t,ch=ch):
        s_=amp.structure(NEW,ch,-t,0.25,2.2)
        return np.nan if s_ is None else abs(s_["TT"]/s_["T"])
    a.plot(T,[rtt(t) for t in T],'-',color=c,lw=2,label=lab)
a.set_xlabel(r'$-t$ [GeV$^2$]'); a.set_ylabel(r'$|\sigma_{TT}/\sigma_T|$'); a.grid(alpha=0.25); a.legend(fontsize=10)
a.set_title(r'amp2026s prediction incl. the neutron, $Q^2$=2.2, $x_B$=0.25',fontsize=10.5)
fig.suptitle('amp2026s: flavour form factors, t-slopes, and ratios',fontsize=14)
fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig("figures/model_amp2026s.png",dpi=120); plt.close(fig)
print("model figure done",flush=True)
