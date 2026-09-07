"""The tied-slope test, started from a grid instead of from the free solution.

Tying b_d to b_u in the Ebar_T sector moves one parameter from 5.28 to 0.67 in
a single step, and least_squares does not recover from it: started at the free
minimum the fit stops with all three normalisations still exactly at 1.0000 and
CLAS6 pi0 at 23 per point, which is not a minimum but a stall.  Start from
several values of the common slope and keep the best.
"""
import json, os, sys
import numpy as np
import fitrun as R, datasets as D

os.environ["TIE_BET"] = "1"
tag = "T3_tie_bET"; keys = sys.argv[1].split(",")
base = np.load("runs/C1_with_compass/fitpar.npy")
os.makedirs("seeds_tmp", exist_ok=True)
seeds = []
for b in (0.7, 1.5, 2.5, 3.5, 4.5):
    q = base.copy(); q[9] = q[14] = b
    # the normalisation of the d block absorbs a slope change to first order
    q[13] = base[13]*np.exp((base[14] - b)*(-0.3))
    f = f"seeds_tmp/tie_{b:.1f}.npy"; np.save(f, q); seeds.append(f)
p, lam = R.fit(keys, seeds=tuple(seeds))
out = f"runs/{tag}"; os.makedirs(out, exist_ok=True)
np.save(f"{out}/fitpar.npy", p)
rec = R.summarise(p, keys, tag, lam)
json.dump(rec, open(f"{out}/summary.json", "w"), indent=1)
print(f"=== {tag}: chi2/ndf = {rec['chi2_ndf']:.4f} ({rec['chi2_fitted']:.1f}/{rec['ndf']})")
print(f"    b(Ebar_T) common = {p[9]:.3f}   n/p = {rec['ratios']['n_over_p']:.3f}")
for k, v in rec["norms"].items(): print(f"    norm {k:>16} = {v:.4f}")
for k in D.ALL:
    s = rec["sets"][k]
    print(f"   {'fit ' if s['fitted'] else '    '}{k:>14}: {s['chi2']:9.1f}/{s['n']:<5d} "
          f"{s['chi2']/max(s['n'],1):6.2f}")
