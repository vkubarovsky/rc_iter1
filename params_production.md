# Production amplitude model &mdash; 2026-08-29

`~/rc_iter1/fitpar_production.npy` ( = fitpar_as_sector.npy)   23 free parameters
chi2 = 1133.6 / 703 = **1.613** on 726 points.

## What is imposed

| constraint | reason | cost |
|---|---|---|
| nQ_d = nQ_u, both sectors | the GPD carries no Q2 dependence and the hard subprocess is flavour-blind up to the quark charge; transversity does not mix with gluons. Verified inside GK_2024 through proton (2u+d) vs neutron (u+2d): the d/u ratio runs as Q^+0.08 and Q^-0.02 | 9.3 |
| one alpha' per sector, shared by u and d | alpha' is the slope of the exchanged trajectory, which cannot depend on the quark. Per-block fitting gave alpha'(H_T^d) = -0.345, which has no Regge meaning | 41.6 |
| slope >= 0 over all xB (0.05 to 1) | a form factor must not grow with abs(t) anywhere | 6.2 |
| b >= 0 as a hard bound | b is the slope at xB = 1 | 1.5 |
| phi_CE = 0 | pure gauge, the three phases have a flat direction | 0.0 |
| b2 = 0, no R_ET prior | the prior constrained moments, ours is a ratio of convolutions | &mdash; |

## Parameters

Convention `exp[(b + b' ln xB) t] Q^nQ`: b is the slope at xB = 1, alpha' = -b'. Same as GK.

| block | N | b [GeV^-2] | alpha' | nQ |
|---|---|---|---|---|
| H_T^u | 18.23 ± 1.09 | 0.000 ± 0.134 | 0.875 ± 0.119 | 1.456 ± 0.159 |
| H_T^d | -3.55 ± 0.89 | 0.000 ± 0.378 | 0.875 (= H_T^u) | 1.456 (= nQ_u) |
| Ebar_T^u | 56.02 ± 4.38 | 0.716 ± 0.111 | 0.276 ± 0.075 | 1.285 ± 0.181 |
| Ebar_T^d | 235.27 ± 31.14 | 5.554 ± 0.519 | 0.276 (= Ebar_T^u) | 1.285 (= nQ_u) |
| T00 | 27.71 ± 4.92 | 0.802 ± 0.278 | 0.948 ± 0.238 | 1.661 ± 0.288 |

| R_L | delta0 | delta1 | rho_CE | phi_w | rho_nf | phi_nf |
|---|---|---|---|---|---|---|
| -0.234 ± 0.234 | 2.001 ± 0.137 | 0.559 ± 0.210 | 0.168 ± 0.033 | -0.938 ± 0.190 | 0.939 ± 0.245 | 1.730 ± 0.246 |

chi2 by block: pi0 537.4/354, eta 240.5/228, BSA CLAS6 pi0 245.9/62, BSA eta 4.2/12,
BSA CLAS12 41.6/30, eg1-dvcs 64.0/40.

In this model H_T^d is simply H_T^u times -0.195: same b (both at the bound, 0.000),
same alpha', same nQ. Its flavour ratio is constant in t and in xB.

## At a physical point (xB = 0.25, Q2 = 2.2, -t = 0.3)

| quantity | value |
|---|---|
| H_T d/u | -0.195 |
| Ebar_T d/u | +0.984 |
| T00 d/u | -0.234 |
| t-slope H_T^u at xB = 0.25 | 1.21 GeV^-2 |
| t-slope H_T^d at xB = 0.25 | 1.21 GeV^-2 |
| t-slope Ebar_T^u at xB = 0.25 | 1.10 GeV^-2 |
| t-slope Ebar_T^d at xB = 0.25 | 5.94 GeV^-2 |
| t-slope T00 at xB = 0.25 | 2.12 GeV^-2 |
| sigma_TT(n)/sigma_TT(p) at the Hall-A point | 1.09 (measured 0.28 ± 0.07) |

## Rejected on the way, and why

| version | why it is out |
|---|---|
| per-block alpha' (25 par, chi2/ndf 1.558) | alpha'(H_T^d) = -0.345: a negative trajectory slope |
| eta entered through the ratio eta/pi0 | sharper (the ratio is known to 5-8% instead of 15-25%) but the model fails it at 2.2 per point, and tying the d slope to cure that drives Ebar_T d/u to -36 and the neutron to 4.4 |
| normalisation nuisance parameters | lambda and rho both settle at 1.000 with no change in chi2: the fit was not buying the steep d slope by sliding the two channels |
| F(xB) = xB^a (1-x)^n on Ebar_T | (1-x)^n does nothing, x^a buys 7 chi2 for one parameter, and the normalisations become degenerate with a |
| b2 t^2 curvature | not needed once the R_ET prior is gone |

## The open problem

Ebar_T^d: t-slope 5.94 +- 0.52 against 1.10 +- 0.11 for u. In GK the two flavours differ
by at most 0.3, and even identical GPD slopes give convolution slopes differing by 0.29
through the x-shapes alone. Four independent things point at it: the published slope test
(the model runs too steep in eta at middle xB), the published eta/pi0 ratio (the model sits
below the data and increasingly so with abs(t)), the Hall-A neutron ratio, and the fact that
every constraint tried is paid for in pi0 sigma_TT below -t = 0.4.

Figures: solid = this model, dashed = amp2026, open squares = this model averaged over the
accepted bin volume. slope_vs_xB.png and ratio_eta_pi0.png reproduce the published tests.