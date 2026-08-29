# Production amplitude model, 2026-08-29 (final of the day)

`~/rc_iter1/fitpar_tieQx_data.npy`   25 free parameters   chi2 = 1090.5 / 701 = **1.556**

Two constraints imposed, both for stated reasons rather than to improve chi2:

1. **nQ_d = nQ_u in both sectors.**  A GFF is (GPD) x (hard subprocess); the GPD
   carries no Q2 dependence (in the GK code the log(Q2/Q0) is computed and never
   used) and the subprocess is flavour-blind up to the quark charge, a constant.
   Transversity does not mix with gluons, so its evolution is non-singlet and
   flavour-blind as well.  Verified in GK_2024 through the proton (2u+d) and
   neutron (u+2d) channels: the d/u ratio runs as Q^+0.079 (Ebar_T) and Q^-0.020
   (H_T), consistent with the 4% expected from the skewness moving at fixed xB.
   Cost: 9.3 chi2 for 2 parameters.  Without it the fitted nQ_d - nQ_u came out
   -7.0 and +2.8 in successive fits, with flavour ratios running by an order of
   magnitude across Q2.
2. **Every t-slope >= 0 over the whole xB range** (0.05 to 1), not just where
   there are data.  Cost: 6.2 chi2.  Before it, H_T^u had b = -0.26: its slope
   crossed zero at xB = 0.785 and the form factor grew with |t| above that.

No R_ET prior of any kind (it constrains a ratio of moments, ours is a ratio of
convolutions), b2 = 0, Ebar_T^d has its own N, b, b', phi_CE = 0 (pure gauge).

## Parameters

Slope convention: exp[(b + b' ln xB) t] Q^nQ, so b is the slope at xB = 1 and
alpha' = -b'.  Same convention as GK.

| block | N | b [GeV^-2] | b' | alpha' = -b' | nQ |
|---|---|---|---|---|---|
| H_T^u | 18.00 ± 0.78 | -0.029 ± 0.039 | -0.894 ± 0.056 | 0.894 | 1.465 ± 0.100 |
| H_T^d | -3.52 ± 0.69 | 0.998 ± 0.369 | 0.333 ± 0.125 | -0.333 | 1.465 (= nQ_u) |
| Ebar_T^u | 57.18 ± 4.65 | 0.602 ± 0.091 | -0.315 ± 0.062 | 0.315 | 0.991 ± 0.198 |
| Ebar_T^d | 251.26 ± 29.50 | 2.731 ± 0.803 | -1.703 ± 0.560 | 1.703 | 0.991 (= nQ_u) |
| T00 | 29.60 ± 5.10 | 0.723 ± 0.275 | -1.012 ± 0.240 | 1.012 | 1.568 ± 0.283 |

| R_L | delta0 | delta1 | rho_CE | phi_w | rho_nf | phi_nf |
|---|---|---|---|---|---|---|
| -0.286 ± 0.208 | 1.991 ± 0.134 | 0.591 ± 0.205 | 0.191 ± 0.032 | -0.824 ± 0.157 | 0.955 ± 0.240 | 1.833 ± 0.213 |

No parameter has a runaway error.

## At a PHYSICAL point (xB = 0.25, Q2 = 2.2, -t = 0.3)

t = 0 is below |tmin| and is not physical -- ratios quoted there are meaningless.

| quantity | value |
|---|---|
| H_T^d / H_T^u | -0.240 |
| Ebar_T^d / Ebar_T^u | +1.302 |
| T00^d / T00^u | -0.286 |
| t-slope H_T^u at xB=0.25 | 1.21 GeV^-2 |
| t-slope H_T^d at xB=0.25 | 0.54 GeV^-2 |
| t-slope Ebar_T^u at xB=0.25 | 1.04 GeV^-2 |
| t-slope Ebar_T^d at xB=0.25 | 5.09 GeV^-2 |
| t-slope T00 at xB=0.25 | 2.13 GeV^-2 |

## chi2 by block

| pi0 xsec | eta xsec | BSA CLAS6 pi0 | BSA eta | BSA CLAS12 | eg1-dvcs |
|---|---|---|---|---|---|
| 503.0 / 354 | 233.7 / 228 | 245.8 / 62 | 3.9 / 12 | 40.4 / 30 | 63.6 / 40 |

## The GK parameters, for comparison (libGKPi0.cpp:89-98, same slope convention)

| | alpha0 | alpha' | b | N |
|---|---|---|---|---|
| GK_2011 H_T^u | -0.17 | 0.45 | 0.30 | +0.830 |
| GK_2011 H_T^d | -0.17 | 0.45 | 0.30 | -0.052 |
| GK_2011 Ebar_T^u | +0.30 | 0.45 | 0.50 | 2.075 |
| GK_2011 Ebar_T^d | +0.30 | 0.45 | 0.50 | 1.345 |
| GK_2024 H_T^u/d | -0.17 | 0.45 | 0.30 | +0.830 / -0.052 |
| GK_2024 Ebar_T^u | -0.10 | 0.45 | **0.77** | 3.351 |
| GK_2024 Ebar_T^d | -0.10 | 0.45 | **0.50** | 2.031 |

alpha' = 0.45 for every block, both flavours, both versions.  Ours ranges from
-0.33 to +1.70, which is the next thing to look at: alpha' is the slope of one
exchanged trajectory and has no business differing between blocks.

## Open problems, in order

1. Ebar_T t-slopes: 1.04 +- 0.04 (u) against 5.09 +- 0.40 (d).  In GK the two
   flavours differ by at most 0.3, and even identical GPD slopes give convolution
   slopes differing by 0.29 through the x-shapes alone.  Tying ours costs 150
   chi2, paid almost entirely by pi0 sigma_TT below -t = 0.4.
2. The Hall-A neutron ratio: this fit gives sigma_TT(n)/sigma_TT(p) = 1.41
   against the measured 0.28 +- 0.07.  Forcing it costs 55 chi2 and lands at
   0.67 (fitpar_tieQx_halla.npy).  The fully tied d = R x u model reaches 0.21
   on its own but costs 150.
3. alpha' scattered from -0.33 to +1.70 across blocks (see above).

Figures: solid = this model, dashed = amp2026, open squares = this model
averaged over the accepted bin volume.