"""One workbook with every pi0/eta data set we have.

Sheets
  datasets        one row per set: what it is, kinematic reach, whether it is fitted
  cross_sections  one row per (Q2,xB,t) bin, ORIGINAL and RC-REFITTED side by side
  asymmetries     asymmetries and phi moments

Conventions
  t is negative everywhere (the physical Mandelstam t).
  sigma_meaning distinguishes sigma_U = sigma_T + eps sigma_L from the
    Rosenbluth-separated sigma_T of the two Hall-A E07-007 sets.
  The CLAS6 pi0 and eta rows are on the PUBLISHED bin kinematics: every phi
    point of the published supplemental is assigned to its own published bin,
    and points belonging to t bins that were never published are dropped.
    The original columns reproduce the published tables (worst bin 0.06 sigma
    for pi0, 0.27 sigma for eta); the *_rc columns are the same extraction
    after our improved radiative correction.
"""
import json, math, os
import numpy as np, pandas as pd

M, MPI, META = 0.9382720813, 0.1349768, 0.547862
def tmin_of(Q2, xB, meson="pi0"):
    """|t| at theta*=0, from Q2 and xB.  For Hall-A y21 this must be evaluated at
    the xB LABEL of the setting, not at <xB>: the published |t| column was built
    as tmin + t', with tmin computed that way, and feeding <xB> instead would
    shift t' by up to a factor three."""
    m = MPI if meson == "pi0" else META
    W2 = M*M + Q2*(1-xB)/xB; W = math.sqrt(W2)
    Eg = (W2 - Q2 - M*M)/(2*W); pg = math.sqrt(Eg*Eg + Q2)
    Ep = (W2 + m*m - M*M)/(2*W); pp = math.sqrt(max(Ep*Ep - m*m, 0.0))
    return -(m*m - Q2 - 2*(Eg*Ep - pg*pp))

# Table I of PRL 127 152301: (Q2, Ebeam) -> <xB>, W2, eps as published
HA21 = {(3.11,7.38):(0.36,6.51,0.61), (3.57,8.52):(0.36,7.29,0.62), (4.44,10.59):(0.36,8.79,0.63),
        (2.67,4.49):(0.48,3.81,0.51), (4.06,8.85):(0.45,5.62,0.71), (5.16,8.85):(0.46,6.67,0.55),
        (6.56,10.99):(0.46,8.32,0.52), (5.49,8.52):(0.59,4.58,0.66), (8.31,10.59):(0.60,6.46,0.50)}
PUB = "/Users/vpk/OneDrive/My_Publications"
ALL = ("/Users/vpk/Library/CloudStorage/OneDrive-JeffersonLab/Jlab_OneDrive/"
       "arXive/pi0_eta_papers/All_experiment.xlsx")
SF  = ["s_u","stat_U","sys_U","s_LT","stat_LT","sys_LT","s_TT","stat_TT","sys_TT"]
SFP = ["s_LTp","stat_LTp","sys_LTp"]        # only Hall-A y21 publishes sigma_LT'
rows_x, rows_a = [], []

def blank_rc():
    return {c+"_rc": np.nan for c in SF}
def blank_ltp():
    return {c: np.nan for c in SFP} | {"tp_low": np.nan, "tp_up": np.nan}
def blank_proj():
    return {"proj": None, "bin_lo": np.nan, "bin_up": np.nan}

# --- CLAS6, both variants side by side --------------------------------------
for meson, src, ref in (("pi0","sf_pi0_both.txt","PRC 90 025205 (2014)"),
                        ("eta","sf_eta_both.txt","PRC 95 035202 (2017)")):
    for v in np.loadtxt(src):
        # Q2, xB, t are the means over the phi points of the bin.  The published
        # bin label is rounded to two decimals, which is not enough for t: it
        # sits in the exponent of the slope, and for the steep Ebar_T^d that
        # rounding is worth up to 5% on the amplitude.
        r = dict(exp="CLAS6_y12", meson=meson, target="p",
                 Q2=v[6], xB=v[7], t=-v[8],
                 Q2_bin=v[0], xB_bin=v[1], t_bin=-v[2], group=f"{meson}{int(v[5]):02d}",
                 tmin=tmin_of(v[6], v[7], meson), tprime=v[8]-tmin_of(v[6], v[7], meson),
                 xB_mean=v[7],
                 eps=v[3], npts_phi=int(v[4]), Ebeam=5.75, sigma_meaning="sigma_U",
                 chi2ndf_phi=v[9], chi2ndf_phi_rc=v[19], in_fit="yes",
                 source=ref+" supplemental, phi table", **blank_proj())
        for k, c in zip(SF, range(10, 19)):  r[k]      = v[c]
        for k, c in zip(SF, range(20, 29)):  r[k+"_rc"] = v[c]
        rows_x.append(r)

# --- Hall-A E07-007, Rosenbluth separated -----------------------------------
EB = {1.50: "3.355 + 5.55", 1.75: "4.455 + 5.55", 2.00: "4.455 + 5.55"}
AE = pd.read_excel(ALL)
for _, q in AE[AE.exp.isin(["HallA_y16", "HallA_y17"])].iterrows():
    p = q.exp == "HallA_y16"; Q2 = float(q.Q2); xB = float(q.xB); mt = abs(float(q.t))
    r = dict(exp=q.exp, meson="pi0", target="p" if p else "n", Q2=Q2, xB=xB, t=-mt,
             tmin=tmin_of(Q2, xB), tprime=mt-tmin_of(Q2, xB),
             xB_mean=xB, Q2_bin=Q2, xB_bin=xB, t_bin=np.nan,
             group=("ha16_%.2f" % Q2) if p else "ha17_1.75",
             eps=np.nan, npts_phi=np.nan,
             Ebeam=EB[Q2] if p else "E07-007 pair", sigma_meaning="sigma_T",
             chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan, in_fit="no" if p else "yes",
             source=("PRL 117 262001 (2016), figure 4" if p else "PRL 118 222002 (2017), figure 5"),
             **blank_rc(), **blank_ltp(), **blank_proj())
    for k, c in zip(SF, ("s_u","stat_U","sys_U","s_LT","stat_LT","sys_LT","s_TT","stat_TT","sys_TT")):
        r[k] = float(getattr(q, c))
    rows_x.append(r)

for line in []:
    v = line.split(); p = True; Q2 = 0.0
    r = dict(exp="HallA_y16" if p else "HallA_y17", meson="pi0", target="p" if p else "n",
             Q2=Q2, xB=float(v[2]), t=-float(v[3]),
             Q2_bin=Q2, xB_bin=float(v[2]), t_bin=np.nan,
             group=("ha16_%.2f" % Q2) if p else "ha17_1.75",
             tmin=tmin_of(Q2, float(v[2])), tprime=float(v[3])-tmin_of(Q2, float(v[2])),
             xB_mean=float(v[2]), eps=np.nan, npts_phi=np.nan,
             Ebeam=EB[Q2] if p else "E07-007 pair", sigma_meaning="sigma_T",
             chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan,
             in_fit="yes" if not p else "no",
             source=("PRL 117 262001 (2016), table I" if p else "PRL 118 222002 (2017)"),
             **blank_rc(), **blank_ltp(), **blank_proj())
    for k, c in zip(SF, (4,5,None,6,7,None,8,9,None)):
        r[k] = 0.0 if c is None else float(v[c])
    rows_x.append(r)

# --- Hall-A unseparated ------------------------------------------------------
NAME11 = {(2.350,0.369):"Kin3", (1.941,0.368):"Kin2",
          (2.155,0.335):"KinX3", (2.073,0.394):"KinX2"}
for v in np.loadtxt("data/halla_y11.data"):
    Q2, xB, eps, tm, tp = v[0], v[1], v[2], v[3], v[4]
    kin = NAME11[min(NAME11, key=lambda k: (k[0]-Q2)**2 + (k[1]-xB)**2)]
    r = dict(exp="HallA_y11", meson="pi0", target="p", Q2=Q2, xB=xB, t=-(tm+tp),
             tmin=tm, tprime=tp, tp_low=np.nan, tp_up=np.nan, xB_mean=xB,
             Q2_bin=Q2, xB_bin=xB, t_bin=np.nan, group=f"ha11_{kin}",
             eps=eps, npts_phi=np.nan, Ebeam=5.752, sigma_meaning="sigma_U",
             chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan, in_fit="no",
             source="PRC 83 025201 (2011), tables IV-VII", **blank_rc(), **blank_proj())
    for k, c in zip(SF, range(5, 14)): r[k] = v[c]
    for k, c in zip(SFP, range(14, 17)): r[k] = v[c]
    rows_x.append(r)

for line in open("data/halla_more.data"):
    if line.startswith("#") or not line.strip(): continue
    v = line.split()
    if v[0] in ("HallA_y21", "HallA_y11"): continue   # rebuilt from the papers above
    Q2, xB, mt, E = float(v[3]), float(v[4]), abs(float(v[5])), float(v[15])
    key = None
    if v[0] == "HallA_y21":
        key = min(HA21, key=lambda k: (k[0]-Q2)**2 + 0.01*(k[1]-E)**2)
        if abs(key[0]-Q2) > 0.02 or abs(key[1]-E) > 0.02: key = None
    r = dict(exp=v[0], meson="pi0", target="p", Q2=Q2, xB=xB, t=-mt,
             Q2_bin=Q2, xB_bin=xB, t_bin=np.nan,
             group=(f"ha21_{sorted(HA21).index(key):02d}" if key else f"ha11_{Q2:.2f}"),
             tmin=tmin_of(Q2, xB), tprime=mt-tmin_of(Q2, xB),
             xB_mean=(HA21[key][0] if key else xB),
             eps=(HA21[key][2] if key else np.nan), npts_phi=np.nan, Ebeam=E,
             sigma_meaning="sigma_U", chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan,
             in_fit="no",
             source="PRL 127 152301 (2021)" if v[0]=="HallA_y21" else "PRC 83 025201 (2011)",
             **blank_rc(), **blank_ltp(), **blank_proj())
    for k, c in zip(SF, range(6, 15)): r[k] = float(v[c])
    rows_x.append(r)

# --- Hall-A y21 from the supplemental table II, with sigma_LT' -----------------
EOF_Y21 = {k[0]: k[1] for k in HA21}          # Q2 -> Ebeam
for v in np.loadtxt("data/halla_y21.data"):
    Q2, xB, tlo, tup, tp = v[0], v[1], v[2], v[3], v[4]
    key = min(HA21, key=lambda k: abs(k[0]-Q2))
    tm = tmin_of(Q2, xB)
    r = dict(exp="HallA_y21", meson="pi0", target="p", Q2=Q2, xB=xB, t=-(tm+tp),
             tmin=tm, tprime=tp, tp_low=tlo, tp_up=tup, xB_mean=HA21[key][0],
             Q2_bin=Q2, xB_bin=xB, t_bin=np.nan,
             group=f"ha21_{sorted(HA21).index(key):02d}",
             eps=HA21[key][2], npts_phi=np.nan, Ebeam=key[1], sigma_meaning="sigma_U",
             chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan, in_fit="no",
             source="PRL 127 152301 (2021), supplemental table II", **blank_rc(),
             **blank_proj())
    for k, c in zip(SF, range(5, 14)): r[k] = v[c]
    for k, c in zip(SFP, range(14, 17)): r[k] = v[c]
    rows_x.append(r)

# --- CLAS12 2025, preliminary -------------------------------------------------
C12 = AE[AE.exp == "CLAS12_y25_v0"]
c12g = {}
for _, q in C12.iterrows():
    key = min(c12g, key=lambda k: (k[0]-q.Q2)**2 + 100*(k[1]-q.xB)**2, default=None)
    if key is None or abs(key[0]-q.Q2) > 0.06 or abs(key[1]-q.xB) > 0.008:
        key = (round(float(q.Q2), 2), round(float(q.xB), 3)); c12g[key] = len(c12g)
    mt = abs(float(q.t))
    r = dict(exp="CLAS12_y25", meson="pi0", target="p", Q2=float(q.Q2), xB=float(q.xB), t=-mt,
             tmin=tmin_of(float(q.Q2), float(q.xB)), tprime=mt-tmin_of(float(q.Q2), float(q.xB)),
             xB_mean=float(q.xB),
             Q2_bin=key[0], xB_bin=key[1], t_bin=np.nan, group=f"c12_{c12g[key]:02d}",
             eps=np.nan, npts_phi=np.nan, Ebeam=float(q.Ebeam), sigma_meaning="sigma_U",
             chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan, in_fit="no",
             source="CLAS12 2025, preliminary (All_experiment.xlsx, provenance to be traced)",
             **blank_rc(), **blank_ltp(), **blank_proj())
    for k, c in zip(SF, ("s_u","stat_U","sys_U","s_LT","stat_LT","sys_LT","s_TT","stat_TT","sys_TT")):
        r[k] = float(getattr(q, c))
    rows_x.append(r)

# --- COMPASS 2025, three projections of the same events ----------------------
# Phys. Lett. B 870 (2025) 139832, tables 7+8 (|t|), 9+11 (Q2), 10+11 (nu).
# This supersedes the 2020 measurement, whose five rows in All_experiment.xlsx
# carried the bin lower edge as t and a single (Q2, xB) for all five bins.
# The three projections are the SAME data binned three ways: at most one may
# enter a fit.  The published systematics are asymmetric; the value stored is
# the mean of the up and down excursions.
for line in open("data/compass_y25.data"):
    if line.startswith("#") or not line.strip(): continue
    v = line.split()
    proj, Q2, mt, xB, eps = v[0], float(v[3]), float(v[5]), float(v[7]), float(v[8])
    sU, esU, syU = float(v[9]), float(v[10]), 0.5*(float(v[11])+float(v[12]))
    sT, esT, syT = float(v[13]), float(v[14]), 0.5*(float(v[15])+float(v[16]))
    # sigma_LT is quoted only for the restricted domain of table 6
    sL, esL, syL = ((float(v[17]), float(v[18]), 0.5*(float(v[19])+float(v[20])))
                    if len(v) > 20 else (np.nan, np.nan, np.nan))
    r = dict(exp="COMPASS_y25", meson="pi0", target="p", Q2=Q2, xB=xB, t=-mt,
             proj=proj, bin_lo=float(v[1]), bin_up=float(v[2]),
             tmin=tmin_of(Q2, xB), tprime=mt-tmin_of(Q2, xB),
             xB_mean=xB, Q2_bin=Q2, xB_bin=xB, t_bin=np.nan, group=f"compass_{proj}",
             eps=eps, npts_phi=np.nan, Ebeam=160.0, sigma_meaning="sigma_U",
             chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan, in_fit="no",
             source=f"COMPASS, PLB 870 139832 (2025), {proj} projection",
             **blank_rc())
    for k, val in zip(SF, (sU, esU, syU, sL, esL, syL, sT, esT, syT)): r[k] = val
    for k in SFP: r[k] = np.nan
    rows_x.append(r)

# --- COMPASS 2020 as it stood in All_experiment.xlsx, kept for the record -----
for _, q in AE[AE.exp == "COMPASS_y20"].iterrows():
    mt = abs(float(q.t))
    r = dict(exp="COMPASS_y20", meson="pi0", target="p", Q2=float(q.Q2), xB=float(q.xB),
             t=-mt, tmin=tmin_of(float(q.Q2), float(q.xB)),
             tprime=mt-tmin_of(float(q.Q2), float(q.xB)),
             xB_mean=float(q.xB), Q2_bin=float(q.Q2), xB_bin=float(q.xB), t_bin=np.nan,
             group="compass_y20", eps=np.nan, npts_phi=np.nan, Ebeam=160.0,
             sigma_meaning="sigma_U", chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan,
             in_fit="no",
             source=("superseded by PLB 870 139832 (2025); t here is the bin lower edge "
                     "and Q2, xB are nominal, not the bin averages"),
             **blank_rc(), **blank_ltp(), **blank_proj())
    for k, c in zip(SF, ("s_u","stat_U","sys_U","s_LT","stat_LT","sys_LT","s_TT","stat_TT","sys_TT")):
        r[k] = float(getattr(q, c))
    rows_x.append(r)

# --- asymmetries -------------------------------------------------------------
def A(**kw): rows_a.append(kw)
for g in json.load(open("data/clas6_demasi_alu.json")):
    for q in g["pts"]:
        A(exp="CLAS6_demasi", meson="pi0", target="p", observable="A_LU^sinphi",
          Q2=g["Q2"], xB=g["xB"], t=-abs(q["t"]), value=q["alpha"], stat=q["err"],
          syst=np.nan, Ebeam=5.776, in_fit="yes",
          source="De Masi et al., PRC 77 042201(R) (2008), figure extraction")
# Zhao: read from the vector figure of the published PDF, not by eye.  The paper
# never went to arXiv, so the PDF is all there is, but its figure 5 is vector and
# the numbers come out exactly.  The earlier digitisation had the values right to
# 0.007 but the statistical errors low by 1.3 -- they had been measured from the
# edge of the marker rather than its centre -- and a flat systematic of 0.087
# where the figure draws a step function of t running 0.022 to 0.064.
for _v in np.loadtxt("data/zhao_vector.data"):
    A(exp="CLAS6_zhao", meson="eta", target="p", observable="A_LU^sinphi",
      Q2=_v[0], xB=_v[1], t=-abs(_v[2]), value=_v[3], stat=_v[4], syst=_v[5],
      Ebeam=5.776, in_fit="yes",
      source="Zhao et al., PLB 789 426 (2019), figure 5 read from the vector PDF")
# CLAS12: the published supplemental gives sigma_LT'/sigma_0 directly, with the
# exact bin averages.  Our earlier figure extraction stored those values divided
# by sqrt(2 eps (1-eps)) instead of left alone, so they were too large by about
# a factor four once the fit multiplied by that factor again.
for _l in open("data/clas12_kim_sigLTp.data"):
    if _l.startswith("#") or not _l.strip(): continue
    _v = [float(_z) for _z in _l.split()]
    A(exp="CLAS12_y24", meson="pi0", target="p", observable="sigma_LT'/sigma_0",
      Q2=_v[0], xB=_v[1], t=-_v[2], value=_v[3], stat=_v[4], syst=_v[5],
      Ebeam=10.6, in_fit="yes",
      source="Kim et al., PLB 849 138459 (2024), supplemental tables III and IV")
EG = json.load(open("data/eg1dvcs_pi0_target_asym.json"))
KIN = {"1.94": 0.25, "2.83": 0.40}
REC = {"E154M5":"A_UL^sinphi","E154M6":"A_UL^sinphi","E154M7":"A_UL^sin2phi",
       "E154M8":"A_UL^sin2phi","E154M9":"A_LL^const","E154M10":"A_LL^const",
       "E154M11":"A_LL^cosphi","E154M12":"A_LL^cosphi"}
for rec, obs in REC.items():
    for row in EG[rec]["rows"]:
        A(exp="eg1dvcs", meson="pi0", target="p (long. polarised)", observable=obs,
          Q2=row[0], xB=KIN[f"{row[0]:.2f}"], t=-abs(row[1]), value=row[2], stat=row[3],
          syst=row[4] if len(row) > 4 else np.nan, Ebeam=5.9, in_fit="yes",
          source=f"Kim et al., PLB 768 168 (2017), CLAS DB {rec}")

# --- Kroll's calculation, kept apart: it is theory, every error is exactly zero
T = AE[AE.exp == "Peter"].copy()
T["target"] = T.target.astype(str).str.strip()
T = T.rename(columns={"exp": "source_tag"})
T["source_tag"] = "Kroll, calculation at one setting (All_experiment.xlsx)"
T = T[["source_tag","meson","target","Q2","xB","t","s_u","s_LT","s_TT","Kin_bin"]]
T = T.rename(columns={"s_u": "sigma_T", "s_LT": "sigma_LT", "s_TT": "sigma_TT",
                      "Kin_bin": "variant"})

# ---- the phi distributions, a sheet of their own ---------------------------
# 60 bins of 12 points from the CLAS database export.  These supersede the 62
# digitised amplitudes: the values agreed to 0.007 but the errors were 1.5 times
# too small, so that block used to read 2.44 per point instead of 0.97.
import pickle as _pk
rows_p = []
for _nm, _d in _pk.load(open("data/demasi_phi.pkl", "rb")):
    for _r in _d:
        if _r[5] <= 0: continue
        rows_p.append(dict(exp="CLAS6_demasi", meson="pi0", target="p", bin=_nm,
                           observable="A_LU(phi)", Q2=float(_r[1]), xB=float(_r[0]),
                           t=-float(_r[2]), phi=float(_r[3]), value=float(_r[4]),
                           stat=float(_r[5]), syst=np.nan, Ebeam=5.776, in_fit="yes",
                           source="De Masi et al., PRC 77 042201(R) (2008), "
                                  "CLAS physics database export"))
P = pd.DataFrame(rows_p)

X, Aa = pd.DataFrame(rows_x), pd.DataFrame(rows_a)
# the digitised De Masi amplitudes stay in the book as a record, out of the fit
Aa.loc[Aa.exp == "CLAS6_demasi", "in_fit"] = "no"
Aa.loc[Aa.exp == "CLAS6_demasi", "source"] = (
    "superseded by the phi distributions on the bsa_phi sheet; these were "
    "digitised from figure 5 and their errors are 1.5 times too small")
sets = []
for k, g in X.groupby(["exp","meson","target"], sort=False):
    sets.append(dict(kind="cross section", dataset=" ".join(k),
        observable=g.sigma_meaning.iloc[0] + ", sigma_LT, sigma_TT"
                   + (", sigma_LT'" if g.s_LTp.notna().any() else ""), rows=len(g),
        Q2=f"{g.Q2.min():.2f}-{g.Q2.max():.2f}", xB=f"{g.xB.min():.3f}-{g.xB.max():.3f}",
        t=f"{g.t.min():.4f}..{g.t.max():.4f}",
        Ebeam=", ".join(sorted({str(v) for v in g.Ebeam})),
        RC_refit="yes" if g.s_u_rc.notna().any() else "no",
        in_fit=g.in_fit.iloc[0], source=g.source.iloc[0]))
sets.append(dict(kind="asymmetry", dataset="CLAS6_demasi", observable="A_LU(phi)",
    rows=len(P), Q2=f"{P.Q2.min():.2f}-{P.Q2.max():.2f}",
    xB=f"{P.xB.min():.3f}-{P.xB.max():.3f}",
    t=f"{P.t.min():.4f}..{P.t.max():.4f}", Ebeam="5.776", RC_refit="no",
    in_fit="yes", source=P.source.iloc[0]))
for k, g in Aa.groupby(["exp","observable"], sort=False):
    sets.append(dict(kind="asymmetry", dataset=k[0], observable=k[1], rows=len(g),
        Q2=f"{g.Q2.min():.2f}-{g.Q2.max():.2f}", xB=f"{g.xB.min():.3f}-{g.xB.max():.3f}",
        t=f"{g.t.min():.4f}..{g.t.max():.4f}", Ebeam=str(g.Ebeam.iloc[0]),
        RC_refit="no", in_fit=g.in_fit.iloc[0], source=g.source.iloc[0]))
S = pd.DataFrame(sets)
cols = ["exp","meson","target","group","proj","bin_lo","bin_up",
        "Q2","xB","t","tmin","tprime","xB_mean",
        "tp_low","tp_up","Q2_bin","xB_bin","t_bin","eps","npts_phi",
        "Ebeam","sigma_meaning","chi2ndf_phi","chi2ndf_phi_rc","in_fit","source"] + SF + SFP + [c+"_rc" for c in SF]
X = X[[c for c in cols if c in X.columns]]
with pd.ExcelWriter("pi0_eta_database.xlsx", engine="openpyxl") as w:
    S.to_excel(w, sheet_name="datasets", index=False)
    X.to_excel(w, sheet_name="cross_sections", index=False)
    Aa.to_excel(w, sheet_name="asymmetries", index=False)
    P.to_excel(w, sheet_name="bsa_phi", index=False)
    T.to_excel(w, sheet_name="theory", index=False)
    # t and xB carry more digits than Excel shows by default
    sh = w.sheets["cross_sections"]
    fmt = {"Q2":"0.000000","xB":"0.000000","t":"0.000000","eps":"0.000000",
           "tmin":"0.000000","tprime":"0.000000","xB_mean":"0.000000",
           "Q2_bin":"0.00","xB_bin":"0.000","t_bin":"0.00",
           "chi2ndf_phi":"0.0000","chi2ndf_phi_rc":"0.0000"}
    for c in SF + SFP + [c+"_rc" for c in SF]: fmt[c] = "0.0000"
    for c in ("tp_low","tp_up"): fmt[c] = "0.00"
    for j, c in enumerate(X.columns, start=1):
        if c in fmt:
            for i in range(2, len(X)+2): sh.cell(row=i, column=j).number_format = fmt[c]
    sa = w.sheets["asymmetries"]
    for j, c in enumerate(Aa.columns, start=1):
        if c in ("t","xB","Q2"): f = "0.000000"
        elif c in ("value","stat","syst"): f = "0.00000"
        else: continue
        for i in range(2, len(Aa)+2): sa.cell(row=i, column=j).number_format = f
print(f"cross_sections {len(X)}   asymmetries {len(Aa)}   bsa_phi {len(P)}   "
      f"theory {len(T)}   datasets {len(S)}\n")
print(S.to_string(index=False))
