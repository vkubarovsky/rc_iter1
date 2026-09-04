# Asymmetry data for the hepgen branch

Four independent measurements, 785 points.  Everything below is stored in the
form it was **published**: no conversions have been applied, because every
conversion we found in the old files turned out to be wrong in one direction or
another.  Convert inside the fit instead, using the formulas given here.

Beam energies: CLAS6 5.776 GeV, CLAS12 10.6 GeV, eg1-dvcs 5.9 GeV.

`eps` everywhere means

    y  = Q2 / (2 M xB E)
    g2 = 4 M^2 xB^2 / Q2
    eps = (1 - y - g2 y^2/4) / (1 - y + y^2/2 + g2 y^2/4)

and `sigma_0 = sigma_T + eps sigma_L`.

---

## 1. CLAS6 pi0, De Masi 2008 — 703 points, phi level

    file    data/demasi_phi.pkl        (pickle: list of (name, Nx6 array))
    raw     2024_CLAS6_pi0_BSA.txt     export from the CLAS physics database
    columns xB, Q2, -t, phi[deg], A, dA
    ref     PRC 77 042201(R) (2008), arXiv:0711.4736

60 bins, 12 phi points each.  **Fit the phi distribution directly**, do not
extract an amplitude first:

    A(phi) = sqrt(2 eps (1-eps)) (sigma_LT'/sigma_0) sin(phi)
             -----------------------------------------------------------------
             1 + sqrt(2 eps (1+eps)) (sigma_LT/sigma_0) cos(phi)
               + eps (sigma_TT/sigma_0) cos(2 phi)

Why the phi level and not the amplitude: the denominator matters.  Extracting
alpha with a plain sin fit, as the paper did, and with the full form differ by a
factor 0.853 in the median.  At the phi level nothing has to be assumed, and the
cos and cos 2phi content constrains sigma_LT and sigma_TT in the same bins.

Of the 703 residuals only 180 carry information -- three moments per bin -- and
the rest is a constant noise floor of about 439 units of chi2 that no model can
remove.  It shifts the total but cannot move the minimum.

## 2. CLAS6 eta, Zhao 2019 — 12 points

    file    data/zhao_vector.data
    columns Q2  xB  |t|  alpha  stat  syst
    ref     PLB 789 426 (2019).  NO arXiv preprint exists.

    alpha = sqrt(2 eps (1-eps)) * sigma_LT' / sigma_0

Read from the vector figure 5 of the published PDF (kept as data/zhao_fig5.svg).
The values of the earlier by-eye digitisation were right to 0.007, but its
statistical errors were 1.3 times too small -- measured from the edge of the
marker rather than its centre -- and it carried a flat systematic of 0.087 where
the figure draws a step function of t running 0.022 to 0.064.  Using the old
numbers makes this block read 0.39 per point instead of 1.67.

## 3. CLAS12 pi0, Kim 2024 — 30 points

    file    data/clas12_kim_sigLTp.data
    columns Q2  xB  |t|  sigma_LT'/sigma_0  stat  syst
    ref     PLB 849 138459 (2024), arXiv:2307.07874, supplemental tables III, IV

**This is the ratio itself, not an asymmetry.**  Compare it directly with
sigma_LT'/sigma_0 and do NOT multiply by sqrt(2 eps (1-eps)).  Our earlier file
had these values divided by that factor, and the fit then multiplied by it
again, so the model was being pushed to a ratio four times the measured one.

Two torus polarities, 15 points each, at the same five (Q2, xB) settings; they
are independent measurements and both should be used.

## 4. eg1-dvcs, polarised target, Kim 2017 — 40 points

    file    data/eg1dvcs_pi0_target_asym.json
    ref     PLB 768 168 (2017), arXiv:1511.03338, CLAS database records E154M*

Two settings: Q2 = 1.94 with xB = 0.25, and Q2 = 2.83 with xB = 0.40, five t
bins each.  Four moments, ten points each:

    record      moment          formula
    E154M5,M6   A_UL^sin(phi)   sqrt(2 eps (1+eps)) sigma_LT^UL  / sigma_0
    E154M7,M8   A_UL^sin(2phi)  eps sigma_TT^UL / sigma_0
    E154M9,M10  A_LL^const      sqrt(1 - eps^2) sigma_T' / sigma_0
    E154M11,M12 A_LL^cos(phi)   sqrt(2 eps (1-eps)) sigma_LT'^LL / sigma_0

These need the four polarised structure functions.  In our amplitude basis they
are, with T00, T01, U01 the helicity amplitudes,

    sigma_LT^UL   = Im[ -sqrt(2) <T00|U01> ]
    sigma_TT^UL   = Im[ 2 <T01|U01> ]
    sigma_T'      = Re[ 2 <T01|U01> ]
    sigma_LT'^LL  = Re[ -sqrt(2) <T00|U01> ]

If the hepgen branch cannot produce these, eg1 has to be left out -- but note
what that costs: eg1 is the only thing in the whole analysis that constrains
sigma_L.  Drop it and sigma_L/sigma_T goes from 0.089 to 0.119, and n/p from
1.02 to 0.59.  The two Rosenbluth separations do not help: their sigma_L is
measured to +-100 to +-320 nb/GeV^2 against a model value of 4 to 61, and some
points come out negative.

---

## Two things about the cross sections, since they cross over

**Hall-A 2011** (PRC 83 025201, arXiv:1003.2938) writes its phi decomposition
with eps_L where everyone else uses eps, and defines eps_L/eps = 4 M^2 xB^2/Q^2
= gamma^2 under its equation (9).  So

    sigma_LT(ours) = gamma * sigma_TL(theirs),   gamma = 2 M xB / Q

and the same for sigma_LT'.  gamma is 0.43 to 0.51 at their four settings, so
the published numbers are about twice ours.  sigma_U and sigma_TT need no
conversion.  The 2016, 2017 and 2021 Hall-A papers all use plain eps.

**COMPASS** builds its cross section on a four-dimensional grid in
(phi, |t|, Q2, nu) and combines the cells by cell volume, its equation (16).
One cell spans xB from 0.02 to 0.47.  Evaluating a model at the quoted <Q2>,
<xB> overshoots by a factor five to eight; averaged over the grid it agrees to
1.0-1.7.  See compass_grid.py.

---

## Known tension, worth keeping in view

At matched kinematics -- 32 pairs within |dxB| < 0.04, |dQ2| < 1.2, |dt| < 0.20
-- CLAS6 gives sigma_LT'/sigma_0 about 2.7 times the CLAS12 value, and with the
corrected errors the mean difference is 1.4 sigma per pair with 5 pairs beyond
2 sigma.  Part of it, a factor 0.853, is the extraction method.  The rest is
real: this is a ratio of structure functions and cannot depend on beam energy.
Dropping De Masi takes the CLAS12 block from 1.23 to 0.65 per point.
