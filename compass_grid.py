"""Average the model over the COMPASS grid, the way the measurement is defined.

PLB 870 139832 builds the cross section on a four-dimensional grid in
(phi, |t|, Q2, nu) -- table 1 -- and equation (16) combines the cells with the
weight of the cell volume, flat in Q2 and nu, not weighted by events.  Our model
lives in (Q2, xB, t), and a cell rectangular in (Q2, nu) is not rectangular in
xB = Q2/(2 M nu): the full grid spans xB from 0.013 to 0.666.  Evaluating the
model at the quoted <Q2>, <xB> therefore has nothing to do with what was
measured -- it overshoots by a factor five to eight.  Averaged properly it
agrees to within a factor 1.0 to 1.7.
"""
import math
import numpy as np
import amplitudes as amp

M = 0.9382720813
Q2E = [1.0, 1.5, 2.24, 3.34, 5.0, 8.0]
NUE = [6.4, 8.5, 11.45, 15.43, 20.78, 26.0, 40.0]
TE  = [0.08, 0.15, 0.22, 0.36, 0.50, 0.64]
GL = np.array([-0.7745966692, 0.0, 0.7745966692]); GW = np.array([5/9, 8/9, 5/9])

def _sig(p, Q2, nu, mt, which):
    xb = Q2/(2*M*nu)
    s = amp.structure(p, "pi0p", -mt, xb, Q2)
    if s is None: return None
    if which == "U":  return s["T"] + amp.epsilon(xb, Q2, 160.)*s["L"]
    if which == "TT": return s["TT"]
    if which == "LT": return s["LT"]
    raise ValueError(which)

def _cell(p, q0, q1, n0, n1, t0, t1, which):
    tot = w = 0.0
    for gq, wq in zip(GL, GW):
        Q2 = 0.5*(q0+q1) + 0.5*(q1-q0)*gq
        for gn, wn in zip(GL, GW):
            nu = 0.5*(n0+n1) + 0.5*(n1-n0)*gn
            for gt, wt in zip(GL, GW):
                mt = 0.5*(t0+t1) + 0.5*(t1-t0)*gt
                v = _sig(p, Q2, nu, mt, which)
                tot += wq*wn*wt*(v if v is not None else 0.0); w += wq*wn*wt
    return tot/w

def average(p, which, q2rng=None, nurng=None, trng=None):
    """Model averaged over the grid cells inside the given ranges, cell-volume
    weighted -- exactly equation (16)."""
    # tables 9 and 10 quote their ranges rounded, so they do not line up with the
    # grid edges of table 1 (1.5-2.1 against the cell 1.5-2.24).  Take every cell
    # whose centre falls inside the quoted range.
    def pick(edges, rng):
        cs = list(zip(edges, edges[1:]))
        if rng is None: return cs
        return [(a, b) for a, b in cs if rng[0] <= 0.5*(a+b) <= rng[1]]
    qs, ns, ts = pick(Q2E, q2rng), pick(NUE, nurng), pick(TE, trng)
    num = den = 0.0
    for q0, q1 in qs:
        for n0, n1 in ns:
            for t0, t1 in ts:
                wgt = (q1-q0)*(n1-n0)*(t1-t0)
                num += _cell(p, q0, q1, n0, n1, t0, t1, which)*wgt
                den += wgt
    return num/den

if __name__ == "__main__":
    p = np.load("fitpar_production_pub.npy")
    d = [l.split() for l in open("data/compass_y25.data")
         if not l.startswith("#") and l.strip()]
    print(f"{'proj':>6}{'bin':>12} | {'sigma_U':>26} | {'sigma_TT':>26}")
    print(f"{'':>6}{'':>12} | {'data':>8}{'point':>8}{'grid':>9}"
          f" | {'data':>8}{'point':>8}{'grid':>9}")
    chi = n = 0.0, 0
    chi, n = 0.0, 0
    for r in d:
        pr = r[0]
        q2r = (float(r[1]), float(r[2])) if pr == "Q2" else None
        nur = (float(r[1]), float(r[2])) if pr == "nu" else None
        tr  = (float(r[1]), float(r[2])) if pr in ("t", "ref27") else None
        Q2, xb, mt, e = float(r[3]), float(r[7]), float(r[5]), float(r[8])
        s = amp.structure(p, "pi0p", -mt, xb, Q2)
        pU = s["T"] + e*s["L"]; pT = s["TT"]
        gU = average(p, "U", q2r, nur, tr); gT = average(p, "TT", q2r, nur, tr)
        dU, eU = float(r[9]), math.hypot(float(r[10]), 0.5*(float(r[11])+float(r[12])))
        dT, eT = float(r[13]), math.hypot(float(r[14]), 0.5*(float(r[15])+float(r[16])))
        chi += ((dU-gU)/eU)**2 + ((dT-gT)/eT)**2; n += 2
        print(f"{pr:>6}{r[1]+'-'+r[2]:>12} | {dU:8.1f}{pU:8.1f}{gU:9.1f}"
              f" | {dT:8.1f}{pT:8.1f}{gT:9.1f}")
    print(f"\nchi2 with the grid-averaged model: {chi:.1f}/{n} = {chi/n:.2f} per point")
