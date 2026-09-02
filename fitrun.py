"""Fit the amplitude model to a chosen list of data sets.

    python3 fitrun.py <tag> <set1,set2,...|all>

Writes runs/<tag>/fitpar.npy and runs/<tag>/summary.json, the latter holding the
chi2 of EVERY data set, whether or not it was fitted, so that a run doubles as a
blind prediction for the sets left out.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F, datasets as D

W = 30.0
BMAX  = float(os.environ.get("BMAX", "12"))    # ceiling on the slope at xB = 1
BDMAX = os.environ.get("BDMAX")                # separate ceiling for b of H_T^d
BL = {"H_T^u": (1, 2), "H_T^d": (5, 6), "Ebar_T^u": (9, 10), "Ebar_T^d": (14, 12),
      "T00": (16, 21)}
BSLOT = [1, 5, 9, 14, 16]
TIES = {7: 3, 27: 11, 6: 2, 12: 10}
# b of H_T^d tied to b of H_T^u.  Left free, the fit runs it to whatever ceiling
# it is given -- 12, 18, 25 -- and buys 2 units of chi2 over 896 points, while
# the normalisation N_d compensates by a factor sixteen.  It is a flat valley,
# not a measurement.
if os.environ.get("TIE_BD", "1") == "1": TIES[5] = 1
FROZEN = {23, 28, 29, 30, 31, 32, 33}
FREE = [i for i in range(37) if i not in set(TIES) | FROZEN]

def slopes(p, x):
    L = math.log(x)
    return {k: (p[b] + p[bp]*L) for k, (b, bp) in BL.items()}

def bounds():
    LO = np.array(list(F.LO) + [-12., -5., -5., -5., -5., -math.pi, -math.pi, -math.pi, 0.5, 0.5])
    HI = np.array(list(F.HI) + [12., 5., 8., 5., 8., math.pi, math.pi, math.pi, 2.0, 2.0])
    LO[3] = LO[11] = LO[17] = -12.; HI[3] = HI[11] = HI[17] = 12.
    LO[12] = -8.; HI[12] = 8.; LO[13] = -1e4; HI[13] = 1e4; LO[14] = -2.; HI[14] = 12.
    for i in BSLOT:
        LO[i] = 0.0
        HI[i] = BMAX          # 12 was inherited from a generic bounds array, not chosen
    if BDMAX is not None: HI[5] = float(BDMAX)
    return LO, HI

def expand(x):
    p = np.zeros(37); p[FREE] = x; p[23] = 0.0
    for t, s in TIES.items(): p[t] = p[s]
    return p

def fit(keys, seeds=("fitpar_production_pub.npy", "fitpar_n_n_and_p.npy")):
    LO, HI = bounds()
    def resid(x):
        p = expand(x); r = []
        for k in keys:
            s = D.BY[k]
            for row in s["rows"]:
                v = s["predict"](p, row)
                r.append((row[4] - v)/row[5] if v is not None else 5.0)
        for xx in (0.05, 0.35, 1.00):
            sl = slopes(p, xx)
            for n in BL: r.append(W*min(0.0, sl[n]))
        return np.array(r)
    best = None
    for src in seeds:
        if not os.path.exists(src): continue
        z = np.load(src)
        for jitter in (0.0, 0.05):
            q = np.zeros(37); q[:len(z)] = z
            if jitter:
                rng = np.random.default_rng(len(keys))
                q = q*(1 + jitter*rng.standard_normal(37))
            for t, s in TIES.items(): q[t] = q[s]
            q = np.clip(q, LO, HI)
            try:
                r = least_squares(resid, q[FREE], bounds=(LO[FREE], HI[FREE]),
                                  x_scale='jac', xtol=1e-13, ftol=1e-13, gtol=1e-13,
                                  max_nfev=60000)
            except Exception:
                continue
            c = float(np.sum(r.fun**2))
            if best is None or c < best[0]: best = (c, r.x)
    return expand(best[1])

def summarise(p, keys, tag):
    used = set(keys)
    rec = dict(tag=tag, fitted=sorted(used), npar=len(FREE), sets={})
    cf = nf = 0.0, 0
    cf, nf = 0.0, 0
    for k in D.ALL:
        c, n = D.chi2(p, k)
        rec["sets"][k] = dict(chi2=c, n=n, fitted=k in used, label=D.BY[k]["label"])
        if k in used: cf += c; nf += n
    rec["chi2_fitted"] = cf; rec["n_fitted"] = nf
    rec["ndf"] = nf - len(FREE)
    rec["chi2_ndf"] = cf/max(nf - len(FREE), 1)
    s25 = slopes(p, 0.25)
    HT, ET, LL = amp._flavour(p, -0.3, 0.25, 2.2)
    sp = amp.structure(p, "pi0p", -0.27, 0.36, 1.75)
    sn = amp.structure(p, "pi0n", -0.27, 0.36, 1.75)
    st = amp.structure(p, "pi0p", -0.4, 0.25, 1.94)
    rec["slopes_xB025"] = {k: float(v) for k, v in s25.items()}
    rec["phases"] = dict(H_T=float(p[34]), Ebar_T=float(p[35]), T00=float(p[36]))
    rec["ratios"] = dict(du_HT=float(abs(HT[1])/HT[0]*np.sign(np.real(HT[1]))),
                         du_ET=float(abs(ET[1])/ET[0]),
                         n_over_p=float(sn["TT"]/sp["TT"]),
                         sigL_over_sigT=float(st["L"]/st["T"]))
    return rec

if __name__ == "__main__":
    tag = sys.argv[1]
    keys = D.ALL if sys.argv[2] == "all" else sys.argv[2].split(",")
    out = f"runs/{tag}"; os.makedirs(out, exist_ok=True)
    p = fit(keys)
    np.save(f"{out}/fitpar.npy", p)
    rec = summarise(p, keys, tag)
    json.dump(rec, open(f"{out}/summary.json", "w"), indent=1)
    print(f"=== {tag}: chi2/ndf = {rec['chi2_ndf']:.4f} "
          f"({rec['chi2_fitted']:.1f}/{rec['ndf']})")
    for k in D.ALL:
        s = rec["sets"][k]
        mark = "fit " if s["fitted"] else "    "
        print(f"   {mark}{k:>12}: {s['chi2']:9.1f}/{s['n']:<5d} {s['chi2']/max(s['n'],1):6.2f}")
