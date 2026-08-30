"""Hall-A sigma_TT off proton and neutron against the models.

The point of the figure: the CLAS-trained model has never seen these data.  The
solid curves are its blind prediction; the dashed ones are the same model after
the Hall-A points were added to the fit; P2 is VPK's earlier GFF fit, also
CLAS-trained.

Data: HallA_pi0_p_n_2016_2017.xlsx -- proton PRL 117 262001, neutron PRL 118 222002.
Statistical errors only; the xlsx carries no systematics.
"""
import sys, math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0,"/Users/vpk/pi0_eta_joint_fit")
import amplitudes as amp
import p2
Q2,xB=1.75,0.36
PROT=[(0.18288,-116.776,19.8532),(0.23190,-243.039,23.1091),(0.28097,-321.289,28.2996)]
NEUT=[(0.18273,-47.397,47.0435),(0.23219,-80.554,42.9743),
      (0.28173,-78.214,45.5613),(0.33133,-152.733,49.5934)]
BLIND=np.load("fitpar_production.npy")            # never saw Hall-A
WITH =np.load("fitpar_han_syst_10_added.npy")     # Hall-A included, 10% syst
tmn=-amp.tmin(amp.Mpi0,Q2,xB)
T=np.linspace(tmn+0.005,0.45,240)
def tt(p,ch): return [amp.structure(p,ch,-t,xB,Q2)["TT"] for t in T]
def tt_p2(ch): return [p2.structure(ch,-t,xB,Q2)[1] for t in T]
fig,axs=plt.subplots(1,2,figsize=(13.5,5.6),sharex=True)
for a,ch,dat,lab in ((axs[0],"pi0p",PROT,r'$\pi^0$ off the proton'),
                     (axs[1],"pi0n",NEUT,r'$\pi^0$ off the neutron')):
    a.plot(T,tt(BLIND,ch),'-',color='#1f4e9c',lw=2.4,label='our fit, blind (Hall-A not in it)')
    a.plot(T,tt(WITH ,ch),'--',color='#1f4e9c',lw=1.8,alpha=0.85,label='our fit, Hall-A included')
    a.plot(T,tt_p2(ch),'-',color='#2E6B4F',lw=1.8,alpha=0.9,label='P2 (earlier GFF fit, CLAS-trained)')
    a.errorbar([d[0] for d in dat],[d[1] for d in dat],yerr=[d[2] for d in dat],
               fmt='o',color='k',ms=8,capsize=4,mfc='white',mew=1.8,
               label='Hall-A (stat. errors only)')
    a.axvline(tmn,color='0.7',lw=1,ls=':'); a.text(tmn+0.004,0.02,r'$|t_{\rm min}|$',color='0.5',fontsize=9,
                                                   transform=a.get_xaxis_transform())
    a.axhline(0,color='k',lw=0.8)
    a.set_xlabel(r'$-t$ [GeV$^2$]'); a.set_title(rf'{lab},  $Q^2$=1.75, $x_B$=0.36',fontsize=12)
    a.grid(alpha=0.25); a.set_xlim(0.1,0.45)
axs[0].set_ylabel(r'$d\sigma_{TT}/dt$ [nb/GeV$^2$]')
axs[0].set_ylim(-460,40); axs[1].set_ylim(-460,40)
axs[0].legend(fontsize=9,loc='lower left')
fig.suptitle('Hall-A $\\sigma_{TT}$ against models trained on CLAS only',fontsize=13.5)
fig.tight_layout(rect=[0,0,1,0.94]); fig.savefig("figures/halla_sigmaTT.png",dpi=135)
print("figures/halla_sigmaTT.png written\n")
print(f"{'target':>8} {'-t':>7} {'data':>9} {'+-':>6} | {'blind':>8} {'pull':>6} | {'with HA':>8} {'pull':>6} | {'P2':>8} {'pull':>6}")
for ch,dat,nm in (("pi0p",PROT,"proton"),("pi0n",NEUT,"neutron")):
    for mt,v,dv in dat:
        b=amp.structure(BLIND,ch,-mt,xB,Q2)["TT"]; w=amp.structure(WITH,ch,-mt,xB,Q2)["TT"]
        q=p2.structure(ch,-mt,xB,Q2)[1]
        print(f"{nm:>8} {mt:7.3f} {v:9.1f} {dv:6.1f} | {b:8.1f} {(v-b)/dv:+6.1f} | {w:8.1f} {(v-w)/dv:+6.1f} | {q:8.1f} {(v-q)/dv:+6.1f}")
for ch,dat,nm in (("pi0p",PROT,"proton"),("pi0n",NEUT,"neutron")):
    cb=sum(((v-amp.structure(BLIND,ch,-mt,xB,Q2)["TT"])/dv)**2 for mt,v,dv in dat)
    cw=sum(((v-amp.structure(WITH ,ch,-mt,xB,Q2)["TT"])/dv)**2 for mt,v,dv in dat)
    cp=sum(((v-p2.structure(ch,-mt,xB,Q2)[1])/dv)**2 for mt,v,dv in dat)
    print(f"  chi2 {nm:>8}: blind {cb:7.1f}   with Hall-A {cw:6.1f}   P2 {cp:6.1f}   ({len(dat)} points)")
