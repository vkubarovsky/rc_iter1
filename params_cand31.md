# Candidate model, night of 2026-08-28/29

31 free parameters.  No R_ET prior of any kind.  Independent Ebar_T^d block,
b2 gone, F(xB) = xB^alpha (1-xB)^n on both Ebar_T blocks, phi_CE fixed at 0
(pure gauge: the three phases have an exactly flat direction), Hall-A neutron
ratio included as a datum.  Only imposed inequality: every t-slope >= 0 on
xB in [0.1, 0.6].   chi2 = 1069.8 / 699 = 1.539.

## Anchored at a PHYSICAL point: xB = 0.25, Q2 = 2.2, -t = 0.3

Never quote anything at t = 0: that point lies below |tmin| and is not physical.
At t = 0 the normalisations carry 50-260% errors and the d/u ratio reads 4-5;
at -t = 0.3 the same model gives the numbers below.

| block | value | error |
|---|---|---|
| H_T^u | 21.919 | 0.494 |
| H_T^d | -3.433 | 0.712 |
| Ebar_T^u | 60.804 | 3.925 |
| Ebar_T^d | 87.795 | 11.786 |
| T00^u | 27.385 | 3.710 |
| T00^d | -4.657 | 7.031 |

| ratio at that point | value |
|---|---|
| Ebar_T^d / Ebar_T^u | +1.444 ± 0.255 |
| H_T^d / H_T^u | -0.157 ± 0.032 |
| T00^d / T00^u | -0.170 ± 0.237 |

| t-slope at xB = 0.25 | value [GeV^-2] |
|---|---|
| H_T^u | 1.17 ± 0.07 |
| H_T^d | 0.47 ± 0.20 |
| Ebar_T^u | 1.09 ± 0.04 |
| Ebar_T^d | 5.11 ± 0.50 |
| T00 | 2.13 ± 0.14 |

## Raw parameters (normalisations here are at t = 0 and are strongly correlated
with alpha and n -- use the anchored table above for anything physical)

| block | N | b | b' | nQ | alpha | n |
|---|---|---|---|---|---|---|
| H_T^u | 17.77 ± 1.15 | -0.112 ± 0.130 | -0.928 ± 0.127 | 1.426 ± 0.160 | - | - |
| H_T^d | -11.64 ± 3.01 | -0.275 ± 0.132 | -0.534 ± 0.233 | -2.742 ± 0.822 | - | - |
| Ebar_T^u | 269.5 ± 131.7 | 1.189 ± 0.155 | 0.073 ± 0.111 | 0.122 ± 0.352 | 0.851 ± 0.239 | 0.105 ± 0.514 |
| Ebar_T^d | 852.9 ± 2224.8 | 4.544 ± 2.271 | -0.407 ± 1.383 | 3.886 ± 1.286 | 0.091 ± 1.100 | 7.463 ± 3.128 |
| T00 | 28.72 ± 5.22 | 0.788 ± 0.287 | -0.965 ± 0.250 | 1.497 ± 0.307 | - | - |

| R_L | delta0 | delta1 | rho_CE | phi_CE | phi_w | rho_nf | phi_nf |
|---|---|---|---|---|---|---|---|
| -0.170 ± 0.237 | 2.042 ± 0.136 | 0.597 ± 0.204 | 0.183 ± 0.033 | 0 (gauge) | -0.978 ± 0.174 | 0.940 ± 0.252 | 1.648 ± 0.230 |

## chi2 by block

| pi0 xsec | eta xsec | BSA CLAS6 pi0 | BSA CLAS6 eta | BSA CLAS12 | eg1-dvcs | Hall-A n/p |
|---|---|---|---|---|---|---|
| 491.3 / 354 | 228.1 / 228 | 246.7 / 62 | ~4 / 12 | ~39 / 30 | 59.7 / 40 | 0.45 vs 0.28 ± 0.07 |
