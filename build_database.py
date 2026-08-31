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
import json, os
import numpy as np, pandas as pd
PUB = "/Users/vpk/OneDrive/My_Publications"
SF  = ["s_u","stat_U","sys_U","s_LT","stat_LT","sys_LT","s_TT","stat_TT","sys_TT"]
rows_x, rows_a = [], []

def blank_rc():
    return {c+"_rc": np.nan for c in SF}

# --- CLAS6, both variants side by side --------------------------------------
for meson, src, ref in (("pi0","sf_pi0_both.txt","PRC 90 025205 (2014)"),
                        ("eta","sf_eta_both.txt","PRC 95 035202 (2017)")):
    for v in np.loadtxt(src):
        r = dict(exp="CLAS6_y12", meson=meson, target="p", Q2=v[0], xB=v[1], t=-v[2],
                 eps=v[3], npts_phi=int(v[4]), Ebeam=5.75, sigma_meaning="sigma_U",
                 chi2ndf_phi=v[5], chi2ndf_phi_rc=v[15], in_fit="yes",
                 source=ref+" supplemental, phi table")
        for k, c in zip(SF, range(6, 15)):   r[k]      = v[c]
        for k, c in zip(SF, range(16, 25)):  r[k+"_rc"] = v[c]
        rows_x.append(r)

# --- Hall-A E07-007, Rosenbluth separated -----------------------------------
EB = {1.50: "3.355 + 5.55", 1.75: "4.455 + 5.55", 2.00: "4.455 + 5.55"}
for line in open("data/halla_pi0.data"):
    if line.startswith("#") or not line.strip(): continue
    v = line.split(); p = v[0] == "p_T"
    Q2 = float(v[1])
    r = dict(exp="HallA_y16" if p else "HallA_y17", meson="pi0", target="p" if p else "n",
             Q2=Q2, xB=float(v[2]), t=-float(v[3]), eps=np.nan, npts_phi=np.nan,
             Ebeam=EB[Q2] if p else "E07-007 pair", sigma_meaning="sigma_T",
             chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan,
             in_fit="yes" if not p else "no",
             source=("PRL 117 262001 (2016), table I" if p else "PRL 118 222002 (2017)"),
             **blank_rc())
    for k, c in zip(SF, (4,5,None,6,7,None,8,9,None)):
        r[k] = 0.0 if c is None else float(v[c])
    rows_x.append(r)

# --- Hall-A unseparated ------------------------------------------------------
for line in open("data/halla_more.data"):
    if line.startswith("#") or not line.strip(): continue
    v = line.split()
    r = dict(exp=v[0], meson="pi0", target="p", Q2=float(v[3]), xB=float(v[4]),
             t=-abs(float(v[5])), eps=np.nan, npts_phi=np.nan, Ebeam=float(v[15]),
             sigma_meaning="sigma_U", chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan,
             in_fit="no",
             source="PRL 127 152301 (2021)" if v[0]=="HallA_y21" else "PRC 83 025201 (2011)",
             **blank_rc())
    for k, c in zip(SF, range(6, 15)): r[k] = float(v[c])
    rows_x.append(r)

# --- COMPASS ----------------------------------------------------------------
d = pd.read_excel("/Users/vpk/hepgen_mac/data/All_experiment.xlsx")
for _, q in d[d.exp == "COMPASS_y20"].iterrows():
    r = dict(exp="COMPASS_y20", meson="pi0", target="p", Q2=q.Q2, xB=q.xB, t=-abs(q.t),
             eps=np.nan, npts_phi=np.nan, Ebeam=q.Ebeam, sigma_meaning="sigma_U",
             chi2ndf_phi=np.nan, chi2ndf_phi_rc=np.nan, in_fit="no",
             source="COMPASS, Phys. Lett. B 805 135454 (2020)", **blank_rc())
    for k in SF: r[k] = getattr(q, k, 0.0) if k in d.columns else 0.0
    rows_x.append(r)

# --- asymmetries -------------------------------------------------------------
def A(**kw): rows_a.append(kw)
for g in json.load(open("data/clas6_demasi_alu.json")):
    for q in g["pts"]:
        A(exp="CLAS6_demasi", meson="pi0", target="p", observable="A_LU^sinphi",
          Q2=g["Q2"], xB=g["xB"], t=-abs(q["t"]), value=q["alpha"], stat=q["err"],
          syst=np.nan, Ebeam=5.776, in_fit="yes",
          source="De Masi et al., PRC 77 042201(R) (2008), figure extraction")
for q in json.load(open("data/clas6_zhao_eta_alu.json")):
    A(exp="CLAS6_zhao", meson="eta", target="p", observable="A_LU^sinphi",
      Q2=q["Q2"], xB=q["xB"], t=-abs(q["t"]), value=q["alpha"], stat=q["stat"],
      syst=q["syst"], Ebeam=5.776, in_fit="yes",
      source="Zhao et al., PLB 789 426 (2019), figure extraction")
for g in json.load(open("data/clas12_kim_alu.json")):
    for q in g["pts"]:
        A(exp="CLAS12_y24", meson="pi0", target="p", observable="A_LU^sinphi",
          Q2=g["Q2"], xB=g["xB"], t=-abs(q["t"]), value=q["A"], stat=q["dA"],
          syst=np.nan, Ebeam=10.6, in_fit="yes",
          source="Kim et al., PLB 849 138459 (2024), figure extraction")
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

X, Aa = pd.DataFrame(rows_x), pd.DataFrame(rows_a)
sets = []
for k, g in X.groupby(["exp","meson","target"], sort=False):
    sets.append(dict(kind="cross section", dataset=" ".join(k),
        observable=g.sigma_meaning.iloc[0] + ", sigma_LT, sigma_TT", rows=len(g),
        Q2=f"{g.Q2.min():.2f}-{g.Q2.max():.2f}", xB=f"{g.xB.min():.3f}-{g.xB.max():.3f}",
        t=f"{g.t.min():.2f}..{g.t.max():.2f}",
        Ebeam=", ".join(sorted({str(v) for v in g.Ebeam})),
        RC_refit="yes" if g.s_u_rc.notna().any() else "no",
        in_fit=g.in_fit.iloc[0], source=g.source.iloc[0]))
for k, g in Aa.groupby(["exp","observable"], sort=False):
    sets.append(dict(kind="asymmetry", dataset=k[0], observable=k[1], rows=len(g),
        Q2=f"{g.Q2.min():.2f}-{g.Q2.max():.2f}", xB=f"{g.xB.min():.3f}-{g.xB.max():.3f}",
        t=f"{g.t.min():.2f}..{g.t.max():.2f}", Ebeam=str(g.Ebeam.iloc[0]),
        RC_refit="no", in_fit=g.in_fit.iloc[0], source=g.source.iloc[0]))
S = pd.DataFrame(sets)
with pd.ExcelWriter("pi0_eta_database.xlsx", engine="openpyxl") as w:
    S.to_excel(w, sheet_name="datasets", index=False)
    X.to_excel(w, sheet_name="cross_sections", index=False)
    Aa.to_excel(w, sheet_name="asymmetries", index=False)
print(f"cross_sections {len(X)}   asymmetries {len(Aa)}   datasets {len(S)}\n")
print(S.to_string(index=False))
