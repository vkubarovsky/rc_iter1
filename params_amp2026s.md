# amp2026s parameter table

26 free parameters, b2 frozen at 0.  chi2/ndf = 1177.2/713 = 1.651.
Form factor of each block:   F(t, xB, Q2) = N * exp[(b + b' ln xB) t + b2 t^2] * Q^nQ

## Block parameters

| block | N | b [GeV^-2] | b' [GeV^-2] | nQ | b2 [GeV^-4] |
|---|---|---|---|---|---|
| H_T^u | 22.17 ± 1.35 | -0.485 ± 0.139 | -1.439 ± 0.164 | 0.948 ± 0.151 | - |
| H_T^d | -7.12 ± 5.11 | -0.485 ± 0.150 | -1.439 ± 0.168 | -6.00 ± 5.47 | - |
| Ebar_T^u | 69.12 ± 4.95 | 0.888 ± 0.115 | -0.240 ± 0.077 | 1.232 ± 0.188 | 0 (frozen) |
| Ebar_T^d | 81.63 (= R_ET x N_u) | 4.593 (= b_u + db) | -0.240 (= b'_u) | 1.232 (= nQ_u) | 0 (frozen) |
| T00 | 28.61 ± 5.02 | 0.800 ± 0.289 | -0.758 ± 0.258 | 1.745 ± 0.287 | - |

## Flavour ratios and the longitudinal sector

| R_HT = N_d/N_u (H_T) | R_ET (Ebar_T) | db_ET [GeV^-2] | R_L = L_d/L_u |
|---|---|---|---|
| -0.321 | 1.181 ± 0.176 | 3.705 ± 0.706 | -0.530 ± 0.159 |

## Phases and the level-3 amplitudes

| delta0 [rad] | delta1 [rad/GeV^2] | rho_CE | phi_CE [rad] | phi_w [rad] | rho_nf | phi_nf [rad] |
|---|---|---|---|---|---|---|
| 2.015 ± 0.132 | 0.666 ± 0.206 | 0.184 ± 0.054 | 0.073 ± 4e4 | -0.751 ± 4e4 | 0.926 ± 0.590 | 1.936 ± 4e4 |

## The same parameters in the previous model, amp2026 (27 par, b2 kept)

| block | N | b | b' | nQ | b2 |
|---|---|---|---|---|---|
| H_T^u | 18.06 | -0.344 | -1.134 | 1.237 | - |
| H_T^d | -1.79 | -0.883 | -0.189 | -0.840 | - |
| Ebar_T^u | 94.82 | 1.762 | -0.219 | 1.353 | 0.408 |
| Ebar_T^d | 79.50 | 5.111 | -0.219 | 1.353 | 0.408 |
| T00 | 25.87 | 0.906 | -0.820 | 1.795 | - |

R_HT -0.099, R_ET 0.838, db_ET 3.349, R_L -0.232, delta0 2.003, delta1 0.604, rho_CE 0.194, rho_nf 0.924

## Derived

| t-slope b + b' ln xB at xB = | 0.10 | 0.15 | 0.25 | 0.40 | 0.60 |
|---|---|---|---|---|---|
| H_T^u | 2.83 | 2.25 | 1.51 | 0.83 | 0.25 |
| H_T^d | 2.83 | 2.24 | 1.51 | 0.83 | 0.25 |
| Ebar_T^u | 1.44 | 1.34 | 1.22 | 1.11 | 1.01 |
| Ebar_T^d | 5.15 | 5.05 | 4.93 | 4.81 | 4.72 |
| T00 | 2.55 | 2.24 | 1.85 | 1.49 | 1.19 |

## chi2 by block

| pi0 xsec | eta xsec | BSA CLAS6 pi0 | BSA CLAS6 eta | BSA CLAS12 pi0 | eg1-dvcs |
|---|---|---|---|---|---|
| 580.1 / 354 | 222.7 / 228 | 243.6 / 62 | 4.1 / 12 | 38.6 / 30 | 69.8 / 40 |

## Three defects

1. nQ of H_T^d is PINNED on its lower bound (-6): the d-quark H_T then falls as
   Q^-3 relative to u (ratio -0.32 at Q2=1, -0.02 at 2.2, -0.001 at 5), and its
   normalisation carries a 70% error.  A fit-boundary artefact, not physics.
2. R_ET = 1.18 pulls the forward-limit prior (0.54 +- 0.15) by +4.3 sigma.
3. phi_CE, phi_w, phi_nf have errors ~4e4: an exactly flat direction, only a
   combination of the three phases is determined.
