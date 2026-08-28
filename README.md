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

## alpha' correction and the reparameterisation (end of 2026-08-28)

We DO have alpha': it is minus our b'.  Matching the GK form
exp[t(b_GK - alpha' ln x)] to ours exp[t(b + b'(ln xB - ln 0.15))] gives

    b' = -alpha'          b_ours = b_GK + 1.897 alpha'

so our b IS the effective slope at xB = 0.15, by construction of the offset.
Our alpha' are 3-5x larger than VPK's global-fit values (0.22 vs 0.04-0.07 for
EbarT; 1.13 vs 0.265 for H_T) -- likely the source of the inverted slope
ordering, since a large alpha' makes the effective slope strongly x-dependent
and we compared at a single x = 0.15.

VPK's proposal (adopt): drop the ln 0.15 offset entirely, b_new = b + 1.897 b'.
EXACT reparameterisation, no refit.  Then slope(xB) = b_new + b' ln xB and the
physical requirement slope > 0 is trivial to check.  Values for amp2026:

    block       N     b_new     b'   | slope at xB = 0.10 0.15 0.25 0.40 0.60
    H_T^u    18.06   -0.344 -1.134   |   2.27  1.81  1.23  0.70  0.24
    H_T^d    -1.79   -0.883 -0.189   |  -0.45 -0.53 -0.62 -0.71 -0.79  NEGATIVE
    EbarT^u  94.82    1.762 -0.219   |   2.27  2.18  2.07  1.96  1.87
    EbarT^d  79.50    5.111 -0.219   |   5.62  5.53  5.41  5.31  5.22
    T00^u/d  25.87/-6.00 0.906 -0.820|   2.79  2.46  2.04  1.66  1.33

H_T^d has a NEGATIVE t-slope over the WHOLE measured xB range -- the form factor
grows with |t| everywhere, not just at an edge.  H_T^u reaches only 0.24 at
xB = 0.6, i.e. the whole H_T sector misbehaves at high xB.

## Step 1 DONE (2026-08-28 evening): reparameterisation + slope constraint

Convention changed everywhere in this repo: `exp[(b + b' ln xB) t]`, i.e. b is
the slope at xB = 1 and -b' is alpha'.  `reparam.py` converts (b_new = b + 1.8971 b',
`to_old` is the inverse) and the change was verified EXACT: 4.99e-14 max relative
difference on 2619 structure-function values over a (channel, xB, Q2, t) grid.

TRAP: parameter files are now in TWO conventions.  New (ln xB): `fitpar_amp2026_lx.npy`,
`fitpar_slope*.npy`.  Old (offset): `fitpar_i1/i2`, `fitpar_A/B/C_prior`,
`fitpar_amp2021_published`, everything in `~/pi0_eta_amplitude_model`, and BOTH
.npy files inside exclurad_py (whose `models/_amplitude_fit.py` still has the
offset).  `iterate.sh` now converts with `reparam.to_old` before installing;
`fit_clas12.py`'s default seed is the new-convention `fitpar_slope.npy`.

Constrained refit (`fit_slope.py`, data/ = iteration 2, R_ET prior kept):
slope >= 0 on xB in [0.1, 0.6] for all five blocks, penalty at the two endpoints
(the slope is linear in ln xB).  Cost: **chi2 1083.7 -> 1086.9, +3.1 for zero new
parameters**, chi2/ndf 1.5504 -> 1.5549.  `fitpar_slope.npy`.

H_T^d lands exactly ON the boundary (slope 0.00, i.e. a t-INDEPENDENT d-quark
H_T; b -0.883 -> -0.006, b' -0.189 -> -0.001, N -1.79 -> -2.87).  So the data
really do pull that slope negative; the constraint is active, not decorative.
R_HT(t=0) -0.099 -> -0.157, R_ET 0.838 -> 0.920 (prior pull +2.0 -> +2.5),
db_ET 3.35 -> 3.04.  In the measured region the model barely moves: sigma_U(pi0)
within 0.1-2.4%, sigma_TT within 3%, sigma_U(eta) within 4%, sigma_L/sigma_T
0.055 -> 0.056.  **The RC fixed point therefore survives - no new RC iteration
is needed for this change.**

### How much does a FALLING H_T^d cost? (`scan_slope.py`, logs/scan_slope2.log)

All blocks kept >= 0, floor scanned on H_T^d alone:

| floor on slope(H_T^d) | 0 | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 | 2.0 | 2.76 | >= slope(H_T^u) |
|---|---|---|---|---|---|---|---|---|---|
| d_chi2 vs amp2026 | +3.1 | +7.2 | +12.5 | +16.7 | +18.4 | +21.4 | +23.7 | +26.4 | +24.5 |

The floor is always active (H_T^d sits on it), but the curve is FLAT: forcing the
d-quark H_T to be as steep as VPK's global fit (2.76) costs only dchi2 = 26 for
zero parameters, and imposing his ORDERING (d steeper than u) costs 24.5.  The
H_T^d slope is therefore only weakly determined - the sign of the disagreement in
the open problem below is real, but its chi2 significance is modest.

WARNING about the first scan (logs/scan_slope.log, superseded): it applied the
same floor to ALL blocks, so at floor >= 1 it was also squeezing H_T^u (whose own
slope is only 0.24 at xB = 0.6) and T00, giving a spurious dchi2 = +285.  The
expensive constraint is on H_T^u at high xB, not on H_T^d.

### Left over: the Ebar_T curvature

The b2 t^2 term is untouched by this constraint and still turns the Ebar_T^u
effective slope (b + b' ln xB + 2 b2 t) negative at -t = 2.7 (xB = 0.15) and 2.45
(xB = 0.40) - INSIDE the declared generator validity window
`_AMP_VALIDITY mt = (0.0, 2.5)` in exclurad_py.  Beyond the turnover Ebar_T grows
with |t|.  Either tighten that window to -t <= 2.0 or constrain b2; the fit
region only reaches -t = 1.75, so nothing in the fit decides it.

## NEXT STEPS, in order

1. Decide what to install: `fitpar_slope.npy` (slope >= 0, +3.1) is the minimal
   physical fix; `fitpar_slope_d2.76.npy` / `fitpar_slope_order.npy` (+26 / +24.5)
   are the variants that agree with the global-fit slope ordering.  Whatever is
   chosen has to be converted with `reparam.to_old` on the way into exclurad_py,
   and the amp2026 generator baseline in OneDrive
   Work/2026_pi0_amplitudes/generator/ re-run against it.
2. Fix the Ebar_T large-|t| turnover (constrain b2, or cap the validity window).
3. Push VPK's global-fit GPDs through our hard kernel and compare convolution to
   convolution -- the only way to settle the inverted slope ordering.  Now sharper:
   we know the chi2 cost of adopting his ordering outright is only ~25.
4. Rewrite ~/pi0eta-rc-2026-paper with the self-consistent result.

## Operational lesson

After pkill on a script name, ALSO kill "spawn_main": multiprocessing workers on
macOS carry only `from multiprocessing.spawn import spawn_main` in their command
line, survive the parent, and keep burning CPU.  One such orphan ran 6h42m at
99% on phallbvpk-mac after a mistaken run there was "stopped".
