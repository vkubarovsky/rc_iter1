"""Two published tests: the slope B(xB) (eta paper Fig. 15, pi0 paper Fig. 17)
and the ratio R = sigma_U(eta)/sigma_U(pi0) (eta paper Fig. 14).

B is defined exactly as in the papers: dsigma_U/dt = A exp(B t), fitted over the
measured t range of each (Q2, xB) bin.  The same fit is applied to the data and
to the model, so the comparison is like for like.
"""
import math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as amp, fit_slope as F
P=np.load("fitpar_as_sector.npy")      # the sector-common alpha' version VPK likes
OLD=np.load("fitpar_amp2026_lx.npy")
def groups(rows):
    gs=[]
    for d in rows:
        for g in gs:
            if abs(g['Q2']-d['Q2'])<0.07 and abs(g['xB']-d['xB'])<0.012: g['pts'].append(d); break
        else: gs.append(dict(Q2=d['Q2'],xB=d['xB'],pts=[d]))
    for g in gs:
        g['Q2']=np.mean([x['Q2'] for x in g['pts']]); g['xB']=np.mean([x['xB'] for x in g['pts']])
    return sorted(gs,key=lambda g:g['xB'])
def expfit(t,y,dy=None):
    m=[i for i in range(len(y)) if y[i]>0]
    if len(m)<3: return float('nan')
    x=np.array([t[i] for i in m]); ly=np.log(np.array([y[i] for i in m]))
    w=np.ones(len(m)) if dy is None else np.array([y[m[i]]/dy[m[i]] for i in range(len(m))])
    A=np.vstack([np.ones(len(x)),-x]).T
    W=np.diag(w*w)
    c=np.linalg.lstsq(A.T@W@A, A.T@W@ly, rcond=None)[0]
    return c[1]
out={}
for rows,ch,tag in ((F.PI0,"pi0p","pi0"),(F.ETA,"etap","eta")):
    res=[]
    for g in groups(rows):
        t=[-d['t'] for d in g['pts']]; y=[d['U'] for d in g['pts']]; dy=[d['dU'] for d in g['pts']]
        Bd=expfit(t,y,dy)
        e=amp.epsilon(g['xB'],g['Q2'],F.E_XS)
        ym=[]; yo=[]
        for d in g['pts']:
            s=amp.structure(P,ch,d['t'],d['xB'],d['Q2']); so=amp.structure(OLD,ch,d['t'],d['xB'],d['Q2'])
            ym.append(s["T"]+e*s["L"]); yo.append(so["T"]+e*so["L"])
        res.append((g['xB'],g['Q2'],Bd,expfit(t,ym),expfit(t,yo),len(t)))
    out[tag]=res
    print(f"\n--- {tag}:  B [GeV^-2] from dsigma_U/dt = A exp(Bt) ---")
    print(f"{'xB':>6} {'Q2':>6} {'n':>3} {'B data':>8} {'B model':>8} {'B amp2026':>10}")
    for xB,Q2,Bd,Bm,Bo,n in res:
        print(f"{xB:6.3f} {Q2:6.2f} {n:3d} {Bd:8.2f} {Bm:8.2f} {Bo:10.2f}")
fig,axs=plt.subplots(1,2,figsize=(13,5))
for a,tag,lab in ((axs[0],"pi0",r'$\pi^0$'),(axs[1],"eta",r'$\eta$')):
    r=out[tag]
    a.errorbar([x[0] for x in r],[x[2] for x in r],fmt='ko',ms=7,label='data (this RC-iterated set)')
    a.plot([x[0] for x in r],[x[3] for x in r],'s',color='#1f4e9c',ms=7,label='model (sector-common alpha\')')
    a.plot([x[0] for x in r],[x[4] for x in r],'^',color='#c23b22',ms=6,alpha=0.7,label='amp2026')
    a.set_xlabel(r'$x_B$'); a.set_ylabel(r'$B$ [GeV$^{-2}$]'); a.grid(alpha=0.25)
    a.set_title(rf'{lab}: slope of $d\sigma_U/dt = A e^{{Bt}}$'); a.legend(fontsize=9)
    a.set_ylim(0,3)
fig.suptitle('Reproduction of the published slope-vs-$x_B$ test (eta paper Fig. 15, pi0 paper Fig. 17)',fontsize=13)
fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig("figures/slope_vs_xB.png",dpi=130); plt.close(fig)
# ---- ratio R = sigma_U(eta)/sigma_U(pi0) ----
print("\n--- R = sigma_U(eta)/sigma_U(pi0) ---")
print(f"{'Q2':>5} {'xB':>5} {'-t':>5} {'R data':>8} {'R model':>8} {'R amp2026':>10}")
pairs=[]
for de in F.ETA:
    m=[d for d in F.PI0 if abs(d['Q2']-de['Q2'])<0.25 and abs(d['xB']-de['xB'])<0.03 and abs(d['t']-de['t'])<0.03]
    if not m: continue
    dp=m[0]; e=amp.epsilon(de['xB'],de['Q2'],F.E_XS)
    def U(p,ch,d):
        s=amp.structure(p,ch,d['t'],d['xB'],d['Q2']); return s["T"]+e*s["L"]
    Rm=U(P,"etap",de)/U(P,"pi0p",dp); Ro=U(OLD,"etap",de)/U(OLD,"pi0p",dp)
    pairs.append((de['Q2'],de['xB'],-de['t'],de['U']/dp['U'],Rm,Ro))
for q in pairs[:14]: print(f"{q[0]:5.2f} {q[1]:5.2f} {q[2]:5.2f} {q[3]:8.2f} {q[4]:8.2f} {q[5]:10.2f}")
a=np.array([[p[3],p[4],p[5]] for p in pairs])
print(f"\nover {len(pairs)} matched points:  R data {np.median(a[:,0]):.2f}, model {np.median(a[:,1]):.2f}, amp2026 {np.median(a[:,2]):.2f} (medians)")
fig,ax=plt.subplots(figsize=(7,5))
ax.plot([p[2] for p in pairs],[p[3] for p in pairs],'ko',ms=6,label='data')
ax.plot([p[2] for p in pairs],[p[4] for p in pairs],'s',color='#1f4e9c',ms=6,label="model (sector-common alpha')")
ax.plot([p[2] for p in pairs],[p[5] for p in pairs],'^',color='#c23b22',ms=5,alpha=0.7,label='amp2026')
ax.axhline(1.0,color='0.5',lw=1)
ax.set_xlabel(r'$-t$ [GeV$^2$]'); ax.set_ylabel(r'$R=\sigma_U(\eta)/\sigma_U(\pi^0)$')
ax.set_title('Published ratio test (eta paper Fig. 14): all matched bins'); ax.grid(alpha=0.25); ax.legend()
fig.tight_layout(); fig.savefig("figures/ratio_eta_pi0.png",dpi=130); plt.close(fig)
print("figures/slope_vs_xB.png and figures/ratio_eta_pi0.png written")
