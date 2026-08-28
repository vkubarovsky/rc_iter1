"""CLAS6 pi0/eta binning (PRC 90 025205 Tables I-III; PRC 95 035202 Tables I-III,
identical grids) and averaging of the model over the ACCEPTED bin volume.

The published Q2, xB, t of a bin are the MEAN values over the accepted volume
assuming constant event density -- not bin centres.  The accepted volume is the
box cut by W > 2 GeV, E' > 0.8 GeV, 21 deg < theta_e < 45 deg, and, in t, by the
physical boundary |t| >= |tmin(Q2, xB)|.  Near tmin that boundary removes a large
part of the first t bins, which is exactly where sigma_TT ~ t' varies fastest.
"""
import math
import numpy as np
import amplitudes as amp

E_BEAM = 5.75
Mp = amp.Mp
Q2_EDGES = [1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.6]
XB_EDGES = [0.10,0.15,0.20,0.25,0.30,0.38,0.48,0.58]
MT_EDGES = [0.09,0.15,0.20,0.30,0.40,0.60,1.00,1.50,2.00]

def accepted(Q2, xB, E=E_BEAM):
    """W > 2, E' > 0.8, 21 < theta_e < 45 deg."""
    nu = Q2/(2*Mp*xB)
    Ep = E - nu
    if Ep <= 0.8: return False
    W2 = Mp*Mp + 2*Mp*nu - Q2
    if W2 <= 4.0: return False
    s = Q2/(4*E*Ep)
    if s <= 0 or s >= 1: return False
    th = 2*math.asin(math.sqrt(s))*180/math.pi
    return 21.0 <= th <= 45.0

def _edge_index(v, edges):
    for i in range(len(edges)-1):
        if edges[i] <= v <= edges[i+1]: return i
    return None

def bin_of(Q2, xB, mt):
    """Bin indices of a published (mean) point; None if it falls outside."""
    return (_edge_index(Q2,Q2_EDGES), _edge_index(xB,XB_EDGES), _edge_index(mt,MT_EDGES))

def nodes(iq, ix, it, ch, nq=5, nx=5, nt=6):
    """Quadrature nodes over the accepted part of a bin, with weights summing to 1.

    (Q2, xB) uses a Gauss-Legendre product rule masked by the acceptance cuts;
    the t integration is done over the PHYSICAL sub-interval
    [max(mt_lo, |tmin|), mt_hi] of each (Q2, xB) node, so the tmin boundary is
    treated exactly rather than by masking.
    """
    mM = amp.Meta if ch == "etap" else amp.Mpi0
    qa,qb = Q2_EDGES[iq], Q2_EDGES[iq+1]
    xa,xb = XB_EDGES[ix], XB_EDGES[ix+1]
    ta,tb = MT_EDGES[it], MT_EDGES[it+1]
    gq,wq = np.polynomial.legendre.leggauss(nq)
    gx,wx = np.polynomial.legendre.leggauss(nx)
    gt,wt = np.polynomial.legendre.leggauss(nt)
    out=[]
    for a,wa in zip(gq,wq):
        Q2 = 0.5*(qb-qa)*a + 0.5*(qb+qa); WQ = wa*0.5*(qb-qa)
        for b,wb in zip(gx,wx):
            xB = 0.5*(xb-xa)*b + 0.5*(xb+xa); WX = wb*0.5*(xb-xa)
            if not accepted(Q2,xB): continue
            mtmin = -amp.tmin(mM,Q2,xB)          # tmin() returns the negative t
            lo = max(ta, mtmin)
            if lo >= tb: continue                 # bin entirely below threshold
            for c,wc in zip(gt,wt):
                mt = 0.5*(tb-lo)*c + 0.5*(tb+lo); WT = wc*0.5*(tb-lo)
                out.append((Q2,xB,mt,WQ*WX*WT))
    if not out: return None
    a=np.array(out); a[:,3] /= a[:,3].sum()
    return a
