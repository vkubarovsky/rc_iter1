"""Export the workbook to db_csv/, the form datasets.py reads.

The workbook is the only source of data.  db_csv is a cache of it, nothing
more: read_excel is slow and datasets.py is imported by every fit and every
plotting script.  Run this after any edit to pi0_eta_database.xlsx, then
check_db_sync() below will pass.
"""
import sys
import pandas as pd

XLSX = "pi0_eta_database.xlsx"
SHEETS = ("datasets", "cross_sections", "asymmetries", "bsa_phi", "theory")

def export():
    for sh in SHEETS:
        pd.read_excel(XLSX, sh).to_csv(f"db_csv/{sh}.csv", index=False)
        print(f"  {sh}")

EPS_TOL = 0.01

def check_eps():
    """The published epsilon is a closure test on (xB, Q2, Ebeam), and it works.

    It is what caught HallA_y21 carrying the group label 0.36/0.48/0.60 of
    table I of PRL 127 152301 instead of the per-setting <xB>: every other set
    agreed to 0.003 while that one was off by 0.05, and the anomaly sat in the
    file unexamined for two days.  Run it, do not eyeball it.
    """
    import sys as _s
    _s.path.insert(0, ".")
    import numpy as np, amplitudes as amp
    X = pd.read_csv("db_csv/cross_sections.csv")
    X = X[X.eps.notna() & X.Ebeam.notna()]
    d = np.array([abs(float(r.eps) - amp.epsilon(float(r.xB), float(r.Q2),
                                                 float(r.Ebeam)))
                  for _, r in X.iterrows()])
    ok = True
    for e in sorted(X.exp.unique()):
        m = (X.exp == e).values
        worst = d[m].max()
        flag = "" if worst <= EPS_TOL else "   <-- FAILS"
        if worst > EPS_TOL: ok = False
        print(f"  eps  {e:14s} worst |published - computed| = {worst:.4f}{flag}")
    return ok

def check():
    import numpy as np
    ok = True
    for sh in SHEETS:
        a = pd.read_excel(XLSX, sh); b = pd.read_csv(f"db_csv/{sh}.csv")
        same = a.shape == b.shape and list(a.columns) == list(b.columns)
        if same:
            na = a.select_dtypes("number")
            same = na.empty or np.nanmax(np.abs(na.values - b[na.columns].values)) < 1e-9
        print(f"  {sh:16s} {'ok' if same else 'DIFFERS'}")
        ok &= same
    return ok

if __name__ == "__main__":
    if "--check" in sys.argv:
        sys.exit(0 if (check() & check_eps()) else 1)
    export()
