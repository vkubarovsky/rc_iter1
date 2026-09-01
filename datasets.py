"""Every data set, loaded from the workbook, with its model prediction.

One place decides what each measurement means -- sigma_T against sigma_U, a
ratio against an asymmetry moment, which flavour combination, which beam energy
-- so that a fit and a plot can never disagree about it.

Each entry of SETS is a dict with
    rows   list of (Q2, xB, mt, obs, value, err, eps)
    label  what to print
    kind   "xs" or "asym"
and a predict(p, row) that returns the model number for that row.
"""
import json, math
import numpy as np
import pandas as pd
import amplitudes as amp

XS = pd.read_csv("db_csv/cross_sections.csv")
AS = pd.read_csv("db_csv/asymmetries.csv")
SYST_N = 0.10          # the neutron xlsx carries no systematics

def _rows_xs(sel, meaning, ebeam=None, obs=("U", "LT", "TT"), rc=False, addsyst=0.0):
    out = []
    for _, r in sel.iterrows():
        e = r.eps if np.isfinite(r.eps) else amp.epsilon(r.xB, r.Q2,
                ebeam if ebeam is not None else float(r.Ebeam))
        for o in obs:
            col = {"U": "s_u", "LT": "s_LT", "TT": "s_TT", "LTp": "s_LTp"}[o]
            sfx = "_rc" if rc else ""
            v = r[col + sfx] if o != "LTp" else r[col]
            st = r[{"U": "stat_U", "LT": "stat_LT", "TT": "stat_TT",
                    "LTp": "stat_LTp"}[o] + (sfx if o != "LTp" else "")]
            sy = r[{"U": "sys_U", "LT": "sys_LT", "TT": "sys_TT",
                    "LTp": "sys_LTp"}[o] + (sfx if o != "LTp" else "")]
            if not np.isfinite(v): continue
            err = math.hypot(st if np.isfinite(st) else 0.0,
                             sy if np.isfinite(sy) else 0.0)
            if addsyst: err = math.hypot(err, addsyst*abs(v))
            if err <= 0: err = SYST_N * abs(v)
            if err <= 0: continue
            out.append((float(r.Q2), float(r.xB), abs(float(r.t)),
                        meaning if o == "U" else o, float(v), err, float(e),
                        str(r.group)))
    return out

def _pred_xs(p, row, ch):
    Q2, xB, mt, obs, _, _, eps, _ = row
    s = amp.structure(p, ch, -mt, xB, Q2)
    if s is None: return None
    if obs == "U":   return s["T"] + eps*s["L"]
    if obs == "T":   return s["T"]
    if obs == "LT":  return s["LT"]
    if obs == "TT":  return s["TT"]
    if obs == "LTp": return s["LTp"]
    raise ValueError(obs)

def _mk_xs(key, label, sel, meaning, ch, ebeam=None, obs=("U","LT","TT"), rc=False,
           addsyst=0.0):
    rows = _rows_xs(sel, meaning, ebeam, obs, rc, addsyst)
    return dict(key=key, label=label, kind="xs", ch=ch, rows=rows,
                predict=lambda p, r, ch=ch: _pred_xs(p, r, ch))

C6 = XS[XS.exp == "CLAS6_y12"]
SETS = [
  _mk_xs("clas6_pi0", "CLAS6 $\\pi^0$ structure functions",
         C6[C6.meson == "pi0"], "U", "pi0p", 5.75, rc=True),
  _mk_xs("clas6_eta", "CLAS6 $\\eta$ structure functions",
         C6[C6.meson == "eta"], "U", "etap", 5.75, rc=True),
  _mk_xs("halla_n", "Hall-A neutron, $\\sigma_T$ separated (2017)",
         XS[XS.exp == "HallA_y17"], "T", "pi0n", 5.55, addsyst=SYST_N),
  _mk_xs("halla_y16", "Hall-A proton, $\\sigma_T$ separated (2016)",
         XS[XS.exp == "HallA_y16"], "T", "pi0p", 5.55, addsyst=SYST_N),
  _mk_xs("halla_y11", "Hall-A proton 6 GeV (2011)",
         XS[XS.exp == "HallA_y11"], "U", "pi0p", 5.752, obs=("U","LT","TT","LTp")),
  _mk_xs("halla_y21", "Hall-A proton 12 GeV (2021)",
         XS[XS.exp == "HallA_y21"], "U", "pi0p", None, obs=("U","LT","TT","LTp")),
  _mk_xs("clas12_xs", "CLAS12 cross sections, preliminary",
         XS[XS.exp == "CLAS12_y25"], "U", "pi0p", 10.604),
  _mk_xs("compass", "COMPASS 2025, $|t|$ projection",
         XS[(XS.exp == "COMPASS_y25") & (XS.group == "compass_t")], "U", "pi0p",
         160.0, obs=("U","TT")),
]

# ---- asymmetries -------------------------------------------------------------
def _asym_rows(sel, E):
    out = []
    for _, r in sel.iterrows():
        err = math.hypot(r.stat if np.isfinite(r.stat) else 0.0,
                         r.syst if np.isfinite(r.syst) else 0.0)
        if err <= 0: continue
        out.append((float(r.Q2), float(r.xB), abs(float(r.t)), r.observable,
                    float(r.value), err, amp.epsilon(r.xB, r.Q2, E), "all"))
    return out

def _pred_bsa(p, row, ch, E):
    Q2, xB, mt, obs, _, _, _, _ = row
    if obs == "A_LU^sinphi":
        return amp.bsa_sinphi(p, ch, -mt, xB, Q2, E)
    s = amp.structure(p, ch, -mt, xB, Q2)
    if s is None: return None
    e = amp.epsilon(xB, Q2, E)
    return s["LTp"]/(s["T"] + e*s["L"])

MP, MPI = 0.9382720813, 0.1349768
def above_threshold(Q2, xB, mt, m=MPI):
    W2 = MP*MP + Q2*(1-xB)/xB; W = math.sqrt(W2)
    Eg = (W2 - Q2 - MP*MP)/(2*W); pg = math.sqrt(Eg*Eg + Q2)
    Ep = (W2 + m*m - MP*MP)/(2*W); pp = math.sqrt(max(Ep*Ep - m*m, 0.0))
    return mt > -(m*m - Q2 - 2*(Eg*Ep - pg*pp))

for key, label, exp, ch, E in (
        ("bsa_demasi", "CLAS6 $\\pi^0$ beam-spin asymmetry", "CLAS6_demasi", "pi0p", 5.776),
        ("bsa_zhao",   "CLAS6 $\\eta$ beam-spin asymmetry",  "CLAS6_zhao",   "etap", 5.776),
        ("bsa_clas12", "CLAS12 $\\sigma_{LT'}/\\sigma_0$",   "CLAS12_y24",   "pi0p", 10.6)):
    m = MPI if ch == "pi0p" else 0.547862
    rows = [r for r in _asym_rows(AS[AS.exp == exp], E)
            if above_threshold(r[0], r[1], r[2], m)]
    SETS.append(dict(key=key, label=label, kind="asym", ch=ch, rows=rows,
                     predict=lambda p, r, ch=ch, E=E: _pred_bsa(p, r, ch, E)))

EG = json.load(open("data/eg1dvcs_pi0_target_asym.json"))
KIN = {"1.94": 0.25, "2.83": 0.40}
REC = {"E154M5": "AULsin", "E154M6": "AULsin", "E154M7": "AULsin2", "E154M8": "AULsin2",
       "E154M9": "ALLc", "E154M10": "ALLc", "E154M11": "ALLcos", "E154M12": "ALLcos"}
_eg = []
for rec, key in REC.items():
    for row in EG[rec]["rows"]:
        Q2, mt, A, dA = row[0], row[1], row[2], row[3]
        ds = row[4] if len(row) > 4 else 0.0
        _eg.append((Q2, KIN[f"{Q2:.2f}"], mt, key, A, math.hypot(dA, ds),
                    amp.epsilon(KIN[f"{Q2:.2f}"], Q2, 5.9), key))
def _pred_eg1(p, row):
    Q2, xB, mt, key, _, _, e, _ = row
    s = amp.structure(p, "pi0p", -mt, xB, Q2)
    if s is None: return None
    s0 = s["T"] + e*s["L"]
    return dict(AULsin=math.sqrt(2*e*(1+e))*s["LTUL"]/s0,
                AULsin2=e*s["TTUL"]/s0,
                ALLc=math.sqrt(1-e*e)*s["Tp"]/s0,
                ALLcos=math.sqrt(2*e*(1-e))*s["LTpLL"]/s0)[key]
SETS.append(dict(key="eg1", label="eg1-dvcs, polarised target", kind="asym",
                 ch="pi0p", rows=_eg, predict=lambda p, r: _pred_eg1(p, r)))

BY = {s["key"]: s for s in SETS}
ALL = [s["key"] for s in SETS]

def chi2(p, key):
    s = BY[key]; c = 0.0; n = 0
    for r in s["rows"]:
        v = s["predict"](p, r)
        if v is None: continue
        c += ((r[4] - v)/r[5])**2; n += 1
    return c, n

if __name__ == "__main__":
    print(f"{'set':>12} {'kind':>5} {'points':>7}   label")
    for s in SETS:
        print(f"{s['key']:>12} {s['kind']:>5} {len(s['rows']):7d}   {s['label']}")
    print(f"{'':>12} {'':>5} {sum(len(s['rows']) for s in SETS):7d}   total")
