# RC iteration: self-consistent radiative corrections with the amplitude model

State as of 2026-08-28, ~14:00. **Read this first after any context loss.**

## What this is

Iterating the loop: amplitude model -> RC -> re-corrected CLAS6 data -> refit
-> new amplitude model -> ... to a fixed point.  Continues
`~/pi0eta-rc-2026-paper` (July 2026), whose finding was that EXCLURAD's RC
depends significantly on the assumed sigma_T/sigma_L.  There the sigma_T-dominant
model was a guess; here sigma_L is measured (amplitude fit, sigma_L/sigma_T = 0.056).

## The correction

    r = eta_py(new model) / eta_py(published model)      sigma_new = sigma_pub / r

BOTH etas are computed with exclurad_py so the CODE cancels and only the MODEL
difference survives.  Published models: pi0.vpk2013, eta.tableB1.
v_cut: pi0 0.18, eta 0.094 (from the original input.dat files, NOT 0.17/0.20).

VALIDATED: eta_py(vpk2013)/delta13_fortran = 1.0021, 16-84% [0.998,1.004] on
1802 matched points.  The python and Fortran RC agree to 0.2%.

## Results so far

| step | pi0 sigma_U shift | eta sigma_U shift | chi2/ndf |
|---|---|---|---|
| published data (amp2021 fit) | -- | -- | 1.479 |
| iteration 1 (RC from amp2021) | -2.9% | -0.8% | 1.548 |
| iteration 2 (RC from i1) | NOT DONE -- see below | | |

r(pi0) median 1.019, r(eta) median 1.007.  The channel-DEPENDENT shift is the
physically important part: it moves the eta/pi0 ratio and hence the flavour
ratios.  Compare July: pi0 -6%..-1.5% with a strong xB slope; ours is nearly
flat, i.e. the amplitude model already sits close to the fixed point in shape.

## Where iteration 2 stopped

Launched (`./iterate.sh i2 fitpar_i1.npy`) and killed mid-way on purpose: the
context was tightening, and while it ran I restored
`exclurad_py/models/amp2021_par.npy` to the published-data parameters, which
swapped the model under the running workers.  Its partial output was deleted.
Nothing downstream used it.  To resume, from a clean tree:

    cd ~/rc_iter1 && ./iterate.sh i2 fitpar_i1.npy      # ~25 min
    # then compare fitpar_i2.npy with fitpar_i1.npy; if the drift is small,
    # freeze as amp2026 and install into exclurad_py.

NOTE the coupling: `iterate.sh` writes the seed into the exclurad_py model
parameter file.  Do not touch that file while an iteration runs.

## Files

    compute_rc.py <ch> <jobs>   RC with both models -> rc_<ch>.txt
    refit_phi.py <ch> <rc> <out>  sigma(phi)=A+Bcos+Ccos2 refit -> sf_<ch>_<it>.txt
    to_strfun.py <ch> <sf> <out>  -> strfun format for the amplitude fit
    fit_clas12.py               amplitude fit (env SEED, OUTP)
    iterate.sh <label> <seed>   one full turn of the loop
    fitpar_i1.npy               iteration-1 parameters
    data/                       CURRENT iteration's data (overwritten each turn!)

## Traps found the hard way

* `rc_factor(h=0)` is NOT helicity-averaged -- it breaks the normalisation
  (eta ~ 1e6).  Use h=+1.
* amp2021 validity floor must be W2 >= 3.3, not 4.0: a 4.0 floor costs up to
  5.6% of eta in the lowest-W2 bins (production/tail_reach.py in exclurad_py).
* `iterate.sh` OVERWRITES exclurad_py/models/amp2021_par.npy.  The frozen
  published-data parameters are ~/pi0_eta_amplitude_model/fitpar_amp27c.npy.

## Variant study (2026-08-28, with the R_ET prior restored)

| variant | par | chi2/ndf | R_HT(t=0) | R_ET(t=0) | db_ET |
|---|---|---|---|---|---|
| A free (= amp2026) | 27 | 1.554 | -0.099 | 0.838 | 3.35 |
| B H_T^d shape tied to u | 25 | 1.561 | -0.026 | 0.781 | 4.29 |
| C GK-like, db=0 both | 23 | 1.651 | +0.249 | 0.050 | 0 |

C is rejected (dchi2 = 75 for 4 par, prior pulled -3.3 sigma, R_HT comes out
POSITIVE against lattice/JAM/GK).  B costs only dchi2 = 8 for 2 par and cures
the H_T^d pathology (b_d -0.525 -> +0.121, c_d -0.840 -> +1.230).

WARNING: the first run of fit_variants_tied.py DROPPED the R_ET prior.  Without
it variant A runs to R_ET = 2.4 -- the exotic d>2u branch -- and gives a BETTER
chi2/ndf (1.529) than amp2026 (1.554).  So part of the chi2 rise 1.479 -> 1.554
is the prior tension, not the RC iteration.  Separate the two before quoting.

## Open problem: t-slope ordering disagrees with VPK's own global GPD fits

VPK's global fits (3/2/26, full t-range, FULL tmin/xi, 194 points, his preferred
version with HT bu = HT alphaStr D = 0) give, as effective slopes b - alpha' ln x
at x = 0.15:

    EbarT^u 1.78   EbarT^d 0.70     (d FLATTER than u)
    H_T^u   0.50   H_T^d   2.76     (d STEEPER than u)

ours (amp2026): EbarT^u 2.18, EbarT^d 5.53 (d steeper); H_T^u 1.81, H_T^d 0.12
(d flatter).  BOTH sectors are inverted.  Robust against his t-range choice
(the |t|<1 fit has the same ordering).

Normalisation ratios are also far apart: his n_d/n_u = 0.17 (EbarT), -0.78 (HT);
ours 0.84 and -0.10.  But these are NOT directly comparable: his x-shapes differ
strongly between flavours (Delta alpha0 = 0.40 for EbarT, 1.06 for HT), which
decouples the convolution ratio from the moment ratio.

NEXT STEP: push his GPDs through our hard kernel and compute the convolutions,
then compare convolution-to-convolution.  Only that decides whether the two
analyses agree.

Also learned from his tables: a negative d-slope is a LIMITED-t-RANGE artefact --
his ET bd was -2.54 for |t|<1 and +0.57 on the full range.  And Peter's
approximate tmin/xi vs the exact ones shifts his parameters by 20-25%; our
amplitude fit uses exact tmin and xi, matching his preferred column.
