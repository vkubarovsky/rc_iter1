"""One table over every run in runs/."""
import glob, json, os
import datasets as D
recs = []
for f in sorted(glob.glob("runs/*/summary.json")):
    recs.append(json.load(open(f)))
print(f"{'run':>20}{'chi2/ndf':>10}{'chi2':>10}{'ndf':>6}   " +
      "".join(f"{k[:10]:>11}" for k in D.ALL))
print("-"*(46+11*len(D.ALL)))
for r in recs:
    row = f"{r['tag']:>20}{r['chi2_ndf']:10.4f}{r['chi2_fitted']:10.1f}{r['ndf']:6d}   "
    for k in D.ALL:
        s = r["sets"][k]; x = s["chi2"]/max(s["n"], 1)
        row += f"{('*' if s['fitted'] else ' ')}{x:10.2f}"
    print(row)
print("\n* = the set was in that fit;  the number is chi2 per point.")
print("\nparameters:")
print(f"{'run':>20}" + "".join(f"{k:>11}" for k in
      ("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")) +
      f"{'d/u H_T':>10}{'|d/u| ET':>10}{'n/p':>8}{'sigL/sigT':>11}")
for r in recs:
    s = r["slopes_xB025"]; q = r["ratios"]
    print(f"{r['tag']:>20}" + "".join(f"{s[k]:11.2f}" for k in
          ("H_T^u","H_T^d","Ebar_T^u","Ebar_T^d","T00")) +
          f"{q['du_HT']:10.3f}{q['du_ET']:10.3f}{q['n_over_p']:8.3f}{q['sigL_over_sigT']:11.4f}")
