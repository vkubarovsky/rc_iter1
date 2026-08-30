# Production amplitude model — 2026-08-30

`~/rc_iter1/fitpar_production.npy`   26 free parameters   chi2 = 1160.8 / 712 = **1.630**
on 738 points.  14 of 14 randomised restarts return to it within 1.0 in chi2.

## Data

| block | points | source |
|---|---|---|
| pi0 cross sections | 354 | CLAS6, PRC 90 025205, RC-iterated |
| eta cross sections | 228 | CLAS6, PRC 95 035202, RC-iterated |
| pi0 BSA | 62 | CLAS6, De Masi PRC 77 042201 |
| eta BSA | 12 | CLAS6, Zhao PLB 789 426 |
| pi0 BSA | 30 | CLAS12, Kim PLB 849 138459 |
| target and double-spin moments | 40 | eg1-dvcs |
| **pi0 off the NEUTRON** | **12** | **Hall-A PRL 118 222002** |

## Parameters

| block | N | b [GeV^-2] | alpha' | nQ | phase [rad] |
|---|---|---|---|---|---|
| H_T^u | 17.77 ± 1.26 | 0.000 ± 0.130 | 0.898 ± 0.120 | 1.480 ± 0.161 | 0 (reference) |
| H_T^d | -7.22 ± 5.38 | 0.000 ± 0.520 | 0.898 (tied) | 1.480 (= nQ_u) | 1.200 ± 0.281 |
| Ebar_T^u | 60.67 ± 4.73 | 0.728 ± 0.111 | 0.294 ± 0.076 | 1.199 ± 0.189 | 0 (reference) |
| Ebar_T^d | 232.31 ± 45.96 | 5.863 ± 0.638 | 0.294 (tied) | 1.199 (= nQ_u) | 0.500 ± 0.483 |
| T00^u | 27.50 ± 4.93 | 0.708 ± 0.283 | 0.983 ± 0.239 | 1.632 ± 0.287 | 0 (reference) |
| T00^d | -0.584 ± 0.305 x u | = u | = u | = u | 0.913 ± 0.442 |

| delta0 | delta1 | rho_CE | phi_w | rho_nf | phi_nf |
|---|---|---|---|---|---|
| 2.050 ± 0.145 | 0.585 ± 0.215 | 0.168 ± 0.034 | -1.027 ± 0.206 | 1.007 ± 0.257 | 1.702 ± 0.251 |

chi2 by block: pi0 550.5/354, eta 235.7/228, BSA CLAS6 pi0 246.4/62, BSA eta 3.4/12,
BSA CLAS12 40.7/30, eg1-dvcs 59.4/40, Hall-A neutron 24.7/12.

## At a physical point (xB = 0.25, Q2 = 2.2, -t = 0.3)

| quantity | value |
|---|---|
| H_T d/u | -0.406 |
| \|Ebar_T d/u\| | 0.820 |
| sigma_TT(n)/sigma_TT(p) at the Hall-A point | 0.970 (measured 0.28 ± 0.07) |
| sigma_L/sigma_T | 0.0582 |
| t-slope H_T^u at xB = 0.25 | 1.25 GeV^-2 |
| t-slope H_T^d at xB = 0.25 | 1.25 GeV^-2 |
| t-slope Ebar_T^u at xB = 0.25 | 1.14 GeV^-2 |
| t-slope Ebar_T^d at xB = 0.25 | 6.27 GeV^-2 |
| t-slope T00 at xB = 0.25 | 2.07 GeV^-2 |

## Two measured phases

| phase | from the fit | model-independent from data | from a valence Regge convolution |
|---|---|---|---|
| delta, longitudinal vs transverse | 1.99 ± 0.14 | **2.05** (from tan delta = -sigma_LT'/sigma_LT on 35 bins) | 0.29-0.47 |
| phi, Ebar_T^d vs u | 0.50 | — | 0.07-0.18 |

The first is a genuine measurement: sigma_LT is small and consistent with zero while
sigma_LT' is clearly positive, which forces the phase near 90-120 degrees.  No valence
Regge x-shape reproduces it: even an intercept of 1.3 gives only 1.14 rad.

## Open

The Ebar_T^d t-slope, 6.27 against 1.14 for u.  It is barely determined -- one restart
found 7.65 at 0.11 lower chi2 -- and it has no counterpart in GK, where the two flavours
differ by at most 0.3.  Adding the Hall-A neutron fixed the SIZE of the d component
(|d/u| 1.98 -> 0.82, into the GK range) but not its t shape.