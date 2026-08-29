# Production amplitude model — 2026-08-29

`~/rc_iter1/fitpar_production.npy`  ( = fitpar_b0.npy )   25 free parameters
chi2 = 1092.0 / 701 = **1.558** on 726 points.

## What is imposed, and why

| constraint | reason | cost in chi2 |
|---|---|---|
| nQ_d = nQ_u, both sectors | the GPD carries no Q2 dependence and the hard subprocess is flavour-blind up to the quark charge; transversity does not mix with gluons, so its evolution is flavour-blind too. Verified in GK_2024 via proton (2u+d) vs neutron (u+2d): the d/u ratio runs as Q^+0.08 (Ebar_T), Q^-0.02 (H_T) | 9.3 for 2 par |
| every t-slope >= 0 over xB in [0.05, 1] | a form factor must not grow with abs(t) anywhere, not only where there are data | 6.2 |
| b >= 0 as a hard bound | b is the slope at xB = 1; a negative b puts the zero crossing inside the physical range (it was at xB = 0.785) | 1.5 |
| phi_CE = 0 | pure gauge: the three phases have an exactly flat direction | 0.0 |
| b2 = 0 | VPK's call; without the R_ET prior the curvature is not needed | — |
| no R_ET prior | it constrains kappa_T^d/kappa_T^u, a ratio of MOMENTS; R_ET is a ratio of CONVOLUTIONS. They agree only if u and d share an x shape | — |

## Parameters

Convention: exp[(b + b' ln xB) t] Q^nQ, so b is the slope at xB = 1 and alpha' = -b'.
Same convention as GK, so the numbers compare directly.

| block | N | b [GeV^-2] | alpha' = -b' | nQ |
|---|---|---|---|---|
| H_T^u | 17.82 ± 1.06 | 0.000 ± 0.123 | 0.871 ± 0.114 | 1.494 ± 0.154 |
| H_T^d | -3.53 ± 0.69 | 1.034 ± 0.387 | -0.345 ± 0.131 | 1.494 (= nQ_u) |
| Ebar_T^u | 57.36 ± 4.67 | 0.587 ± 0.111 | 0.324 ± 0.072 | 0.973 ± 0.214 |
| Ebar_T^d | 253.87 ± 31.34 | 2.706 ± 0.816 | 1.721 ± 0.566 | 0.973 (= nQ_u) |
| T00 | 29.82 ± 5.19 | 0.709 ± 0.279 | 1.024 ± 0.244 | 1.555 ± 0.287 |

| R_L | delta0 | delta1 | rho_CE | phi_w | rho_nf | phi_nf |
|---|---|---|---|---|---|---|
| -0.292 ± 0.207 | 1.990 ± 0.134 | 0.592 ± 0.205 | 0.191 ± 0.032 | -0.825 ± 0.157 | 0.959 ± 0.241 | 1.833 ± 0.213 |

No parameter has a runaway error.  chi2 by block: pi0 502.8/354, eta 235.1/228,
BSA CLAS6 pi0 245.8/62, BSA eta 3.9/12, BSA CLAS12 40.4/30, eg1-dvcs 63.6/40.

## At a physical point (xB = 0.25, Q2 = 2.2, -t = 0.3)

t = 0 lies below abs(tmin) and is not physical; ratios quoted there are meaningless.

| quantity | value |
|---|---|
| H_T^d / H_T^u | -0.241 |
| Ebar_T^d / Ebar_T^u | +1.311 |
| T00^d / T00^u | -0.292 |
| t-slope H_T^u at xB = 0.25 | 1.21 GeV^-2 |
| t-slope H_T^d at xB = 0.25 | 0.56 GeV^-2 |
| t-slope Ebar_T^u at xB = 0.25 | 1.04 GeV^-2 |
| t-slope Ebar_T^d at xB = 0.25 | 5.09 GeV^-2 |
| t-slope T00 at xB = 0.25 | 2.13 GeV^-2 |

## The common-alpha' test (not adopted, but a result)

alpha' is the slope of the exchanged trajectory and should be one number for
every block.  Forcing that: **alpha' = 0.56**, chi2 1092 -> 1182, i.e. 92 for 4
parameters.  GK uses 0.45 for every block, both flavours, both versions.
So the single value the data pick sits next to the Regge expectation; the
scatter we see when each block is free (-0.35 to +1.72) is freedom, not physics.

Releasing one block at a time from the common alpha' (all else tied):

| released | chi2 | recovered | the alpha' it wants |
|---|---|---|---|
| none | 1182.2 | — | 0.56 |
| Ebar_T^d | 1147.2 | 35.0 | 3.16 (and its b drops 5.05 -> 1.04) |
| H_T^d | 1163.8 | 18.4 | -0.57 |
| T00 | 1172.4 | 9.8 | 1.13 |

Two thirds of the resistance is the two d blocks.  A common alpha' is compatible
with everything except them.

## GK for comparison (libGKPi0.cpp:89-98, same slope convention)

| | alpha0 | alpha' | b | N |
|---|---|---|---|---|
| GK_2011 H_T^u / H_T^d | -0.17 | 0.45 | 0.30 / 0.30 | +0.830 / -0.052 |
| GK_2011 Ebar_T^u / ^d | +0.30 | 0.45 | 0.50 / 0.50 | 2.075 / 1.345 |
| GK_2024 H_T^u / H_T^d | -0.17 | 0.45 | 0.30 / 0.30 | +0.830 / -0.052 |
| GK_2024 Ebar_T^u / ^d | -0.10 | 0.45 | 0.77 / 0.50 | 3.351 / 2.031 |

## The one open problem

Ebar_T^d.  Its t-slope is 5.09 +- 0.40 against 1.04 +- 0.04 for u; in GK the two
flavours differ by at most 0.3, and even identical GPD slopes give convolution
slopes differing by 0.29 through the x-shapes alone.  Every constraint we have
tried this week is paid for by this one block, and it always pays in pi0 sigma_TT
below -t = 0.4.  It also drives the Hall-A conflict: this model gives
sigma_TT(n)/sigma_TT(p) = 1.41 against the measured 0.28 +- 0.07, while the model
with d = R x u reaches 0.21 on its own and costs 150 chi2.

Figures: solid = this model, dashed = amp2026, open squares = this model averaged
over the accepted bin volume.