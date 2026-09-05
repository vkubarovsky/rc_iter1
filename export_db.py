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
        sys.exit(0 if check() else 1)
    export()
