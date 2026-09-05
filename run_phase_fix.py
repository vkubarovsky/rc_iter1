"""Refit with the phase bounds repaired, from a grid of phase starting points.

p[35] and p[36] sat at 0.5 only because that was the floor of a window left over
from a previous meaning of those slots.  Reopening the window is not enough: the
seeds still start there, so the fit can stay in the same local minimum.  Start it
from a grid instead and keep the best.
"""
import json, os, sys
import numpy as np
import fitrun as R, datasets as D

tag = sys.argv[1]; keys = sys.argv[2].split(",")
base = np.load(sys.argv[3])
# GRID=full explores the phase plane from scratch; GRID=base (the default once
# the phases are known) seeds from the given solution and lets fitrun jitter.
os.makedirs("seeds_tmp", exist_ok=True)
seeds = []
G = (-2.0, 0.0, 2.0) if os.environ.get("GRID") == "full" else (base[35],)
H = (-2.0, 0.0, 2.0) if os.environ.get("GRID") == "full" else (base[36],)
for a in G:
    for b in H:
        q = base.copy(); q[35] = a; q[36] = b
        f = f"seeds_tmp/s_{a:+.1f}_{b:+.1f}.npy"; np.save(f, q); seeds.append(f)
p, lam = R.fit(keys, seeds=tuple(seeds))
out = f"runs/{tag}"; os.makedirs(out, exist_ok=True)
np.save(f"{out}/fitpar.npy", p)
rec = R.summarise(p, keys, tag, lam)
json.dump(rec, open(f"{out}/summary.json", "w"), indent=1)
print(f"=== {tag}: chi2/ndf = {rec['chi2_ndf']:.4f} ({rec['chi2_fitted']:.1f}/{rec['ndf']})")
print(f"    phases  H_T {p[34]:+.3f}   Ebar_T {p[35]:+.3f}   T00 {p[36]:+.3f}")
for k, v in rec["norms"].items():
    print(f"    norm {k:>16} = {v:.4f}")
for k in D.ALL:
    s = rec["sets"][k]
    print(f"   {'fit ' if s['fitted'] else '    '}{k:>14}: {s['chi2']:9.1f}/{s['n']:<5d} "
          f"{s['chi2']/max(s['n'],1):6.2f}")
