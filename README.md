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

### Is the d-steeper-than-u ordering acceptable?  YES (tcut_test.py)

Imposing slope(H_T^d) >= slope(H_T^u) costs dchi2 = +21.3 over the full range, and
ALL of it sits in one place: eta sigma_U at -t > 1.2 (+20.6 on 26 points).  pi0
sigma_U improves by 5 and pi0 sigma_LT by 2.6; sigma_TT and the BSA blocks do not
move.  Restricting the fit to -t <= 1.2, where the handbag/twist-3 description is
meant to hold at all, the cost drops to **+6.1 chi2 on 513 points** - noise.

Independent support: forcing his slope ordering also drags the H_T normalisation
ratio to his value without being asked.  R_HT = -0.157 (slope >= 0 only) ->
-0.516 (d >= u, full range) -> -0.880 (d >= u, -t <= 1.2), against VPK's global
n_d/n_u = -0.78.  Two separate pieces of his global fit come out consistent at
once, so this is not a fit to a single number.

Caveat: the constraint SATURATES - the fit puts d exactly equal to u (1.30-1.37
at xB = 0.25), never strictly steeper.  A genuinely steeper d has to be imposed
by hand (fitpar_slope_d2.76.npy, +26.4 on the full range).

### Left over: the Ebar_T curvature

The b2 t^2 term is untouched by this constraint and still turns the Ebar_T^u
effective slope (b + b' ln xB + 2 b2 t) negative at -t = 2.7 (xB = 0.15) and 2.45
(xB = 0.40) - INSIDE the declared generator validity window
`_AMP_VALIDITY mt = (0.0, 2.5)` in exclurad_py.  Beyond the turnover Ebar_T grows
with |t|.  Either tighten that window to -t <= 2.0 or constrain b2; the fit
region only reaches -t = 1.75, so nothing in the fit decides it.

## PRODUCTION MODEL: amp2026s (2026-08-28, VPK's call)

`fit_prod.py` -> `fitpar_prod.npy`, installed by `install_amp2026s.sh` as
exclurad_py models `pi0.amp2026s` / `eta.amp2026s` / `pi0n.amp2026s`.
Constraints, all at zero new parameters:
  b + b' ln xB >= 0 per block on xB in [0.1, 0.6];
  slope(H_T^d) >= slope(H_T^u)  (the global-fit ordering);
  Ebar_T slope including 2 b2 t >= 0 out to -t = 2.5 -- FREE (dchi2 = -0.1);
    it only moves b2 0.391 -> 0.365 and pushes the turnover to -t = 2.49.
chi2 = 1108.1, chi2/ndf 1.5853 (amp2026: 1083.7 / 1.5504).  H_T slopes at
xB = 0.25: u = d = 1.39 (the ordering constraint saturates).  R_HT = -0.495,
R_ET = 0.957.  Verified end to end: exclurad_py amp2026s vs fitpar_prod agrees
to 2e-14, so the to_old conversion on install is right.

amp2026 is left registered and frozen - the OneDrive generator baseline
(Work/2026_pi0_amplitudes/generator/) was produced with it.

WHAT ACTUALLY MOVED (Q2 = 2.2, xB = 0.25): pi0 sigma_T within 5% out to -t = 1.75
(-15% at 2.4, the curvature constraint biting), but **eta sigma_T falls by 6.5%
at -t = 0.2, 17% at 1.2, 26% at 1.75, 36% at 2.4**.  That is the ordering being
paid for, and it is large enough that the RC fixed point should be re-checked:
`./iterate.sh i3 fitpar_prod.npy` (~25 min, overwrites data/ and the exclurad_py
amp2021 parameters).  Generator runs must be redone against amp2026s.

### Is b2 needed at all?  YES, strongly (fit_nob2.py)

Freezing b2 = 0 under the same production constraints: chi2 1108.1 -> 1179.1,
**dchi2 = +71 for one parameter**.  The damage is exactly where the term was
invented for: pi0 sigma_TT at -t < 0.8 goes 123.2 -> 155.2, while -t > 0.8 is
unchanged (41.4 -> 40.8) and eta does not move (227.4 -> 229.8).  Without the
curvature the exponential has to flatten to reach the large-|t| points
(Ebar_T^u slope at xB = 0.25: 2.07 -> 1.22) and then misses at small |t|.

So the protez is load-bearing: the data demand curvature in Ebar_T^u at small
|t|.  Its unphysical large-|t| tail is contained by the -t <= 2.5 constraint at
zero chi2 cost.  If a shape that is monotone by construction is wanted for the
paper (dipole-like exp(bt)/(1 - t/L^2)^n instead of exp(bt + b2 t^2)), that is a
form change to be tested separately - not needed for production.

## FINAL MODEL: amp2026s with b2 REMOVED (2026-08-28, VPK's call)

"We don't need exact data description.  We have what we have."  b2 was mostly a
near-tmin fudge (see below), so it is gone: 26 parameters, every form factor a
pure falling exponential, no slope reversal possible at any t.
`fit_final.py` -> `fitpar_amp2026s.npy`, installed as exclurad_py `*.amp2026s`
(verified end to end, 1.7e-14).  All EIGHT seeds converge to the same minimum.

| model | par | chi2 | chi2/ndf |
|---|---|---|---|
| amp2026 (unconstrained) | 27 | 1083.7 | 1.5504 |
| slope >= 0 only | 27 | 1086.9 | 1.5549 |
| + d >= u ordering, b2 kept | 27 | 1108.1 | 1.5853 |
| **+ b2 removed = amp2026s** | **26** | **1179.1** | **1.6845** |

Blocks: pi0 587.0/354, eta 229.8/228, bsa_pi0 246.4/62, bsa_eta 6.0/12,
bsa_c12 38.2/30, eg1 71.8/40.  Slopes (b + b' ln xB), all positive everywhere:

    block        b     b'   |  xB = 0.10  0.15  0.25  0.40  0.60
    H_T^u    -0.468 -1.448  |      2.86  2.28  1.54  0.86  0.27
    H_T^d    -0.468 -1.447  |      2.86  2.28  1.54  0.86  0.27   (constraint saturated)
    Ebar_T^u  0.854 -0.262  |      1.46  1.35  1.22  1.09  0.99
    Ebar_T^d  4.548 -0.262  |      5.15  5.05  4.91  4.79  4.68
    T00       0.751 -0.723  |      2.42  2.12  1.75  1.41  1.12

R_HT = -0.357, db_ET = 3.694, R_L = -0.657, delta(t) = 1.98 + 0.68 t,
sigma_L/sigma_T = 0.0560 (unmoved through all of this).

**WATCH THIS: R_ET = 1.214, a +4.5 sigma pull on the forward-limit prior
(0.54 +- 0.15).**  Removing the curvature is paid for in the Ebar_T d/u ratio -
the fit needs the steep d exponential to bend the pi0 sigma_TT shape that b2
used to bend.  It was +2.0 in amp2026 and +2.7 with b2 kept.  Everything that
depends on R_ET (neutron predictions, the flavour interpretation) now rests on a
number the prior disagrees with.

Generator-level change vs amp2026 at Q2 = 2.2, xB = 0.25 (this is what the
overnight sweep will show): pi0 sigma_T +1% at -t = 0.15, -6% at 0.3, +7% at 1.0,
**-42% at 2.0**; pi0 sigma_TT **-30% at 0.15**, -20% at 0.3, +16% at 1.0, -45% at
2.0; eta sigma_T -12% to -28% over 0.3-1.5 and -59% at 2.0.  The large-|t| drop
is the removed b2 tail; the small-|t| sigma_TT drop is the removed curvature.

### Overnight sweep launched 2026-08-28 18:36, restarted 18:39 and split 18:56

The first launch hit the generator's own warning: with sigma_LT' nonzero and
EXACT_ACCEPT off, the eta table is built at a single helicity and the sample's
beam-spin asymmetry is biased.  Our models have sigma_LT' by construction and the
morning baseline recorded no environment, so the baseline cannot serve as the
comparison set.  Killed, cleaned, and relaunched with the frozen production
exports, generating BOTH models under identical conditions.

Split across the two machines by CHANNEL (VPK freed phallbvpk-mac), so that every
model-to-model comparison stays inside one machine:

    pi0   vpkmacmini      /Volumes/wd_14tb/mc/ampgen_s_run/run_pi0.sh  --jobs 8
    eta   phallbvpk-mac   ~/ampgen_cmp_run/run_eta.sh                  --jobs 10

phallbvpk-mac runs an rsync of the mini working tree (~/exclurad_mini): its git
clone was at 5be26dd with no amplitude model files, and generating on older
generator code is exactly the mismatch farm_switch.sh warns about.  caffeinate
holds it awake.  Both scripts honour a STOP file between jobs.  Passes: 200k then
1M per mode; eta results live on phallbvpk-mac and must be rsynced back to
/Volumes/wd_14tb/mc/ampgen_cmp/ before plotting.  Full note in that run dir's README.

## BIN-AVERAGED REFIT -- the final amp2026s (2026-08-28 night)

VPK supplied the two publications; both give the same grid (PRC 90 025205 and
PRC 95 035202, Tables I-III): 7 Q2 bins 1.0-4.6, 7 xB bins 0.10-0.58, 8 |t| bins
0.09-2.00, and -- the point -- the quoted Q2, xB, t of a bin are the MEAN over
the accepted volume, not bin centres.  `bins.py` rebuilds that volume (box cut by
W > 2, E' > 0.8, 21-45 deg, and by the |t| >= |tmin(Q2,xB)| boundary that sweeps
across the low-t bins).  VALIDATION: the reconstructed volume means reproduce the
published means -- eta rms 0.0012 in -t, 0.009 in Q2; pi0 rms 0.017 in -t,
0.014 in Q2, all 194 bins matched.  `fit_binned.py` then averages the model over
that volume inside the fit (15522 quadrature nodes; BSA and eg1 stay
point-evaluated, their binning is not in these papers).

| variant | par | chi2 | chi2/ndf | b2 | R_ET | slope H_T^d at xB=0.25 |
|---|---|---|---|---|---|---|
| free (no slope constraints) | 27 | 1080.1 | 1.5452 | 0.399 | 0.817 (+1.8) | **-0.60** |
| + slope >= 0 and d >= u | 27 | 1090.7 | 1.5604 | 0.382 | 0.917 (+2.5) | 1.37 |
| **+ b2 removed = amp2026s** | **26** | **1159.0** | **1.6557** | 0 | 1.181 (+4.3) | 1.51 |

TWO NEGATIVE RESULTS, both worth keeping:
* Bin averaging does NOT remove b2.  Dropping it still costs +68 (it was +71
  point-evaluated).  The near-tmin leverage is real, but the curvature the data
  want is not a bin-centring artefact.
* Bin averaging does NOT remove the H_T^d pathology.  The free fit still runs to
  a slope of -0.60 at xB = 0.25 (-0.62 point-evaluated).  The growing d-quark
  form factor is what the data prefer; only the constraint stops it.
The bias itself is small: chi2 1083.7 -> 1080.1 free, 1108.1 -> 1090.7 constrained.

PRODUCTION MODEL (VPK's call: remove b2, refit, install): the 26-parameter
bin-averaged constrained fit, `fitpar_amp2026s.npy`, installed as exclurad_py
`pi0.amp2026s` / `eta.amp2026s` / `pi0n.amp2026s`, verified to 8e-15.
R_HT = -0.321, R_ET = 1.181, db_ET = 3.705, slopes u = d = 1.51 at xB = 0.25.
The +4.3 sigma R_ET pull is the standing caveat on this model.

## Production runs launched 2026-08-28 19:52

`pi0.amp2026s` on vpkmacmini, `eta.amp2026s` on phallbvpk-mac, both from GitHub
branch `amp2026s` commit 1dfc6ef (the laptop got a fresh clone in
~/exclurad_amp2026s; its old rsync tree and the earlier ampgen_cmp attempts were
deleted).  200k born + 200k rad + 1M born + 1M rad per channel, frozen production
conventions with EXACT_ACCEPT=1.  200k born done in under a minute per channel
(200000 events, 8/10 parts).  Run dir and full note:
/Volumes/wd_14tb/mc/ampgen2026s_run/README.md.

## NIGHT OF 2026-08-28/29: the prior comes out, Ebar_T^d gets its own block

VPK's corrections, in the order they landed, each one right:

1. **The R_ET prior was applied at the wrong level.**  R_ET is a ratio of
   CONVOLUTIONS (GPD x hard kernel), kappa_T^d/kappa_T^u = 0.54 +- 0.15 is a
   ratio of first MOMENTS.  They coincide only if u and d share an x shape.
   The prior entered as a literal extra residual, `(p[13]-0.54)/0.15`, in every
   fit since P2.  Removed.
2. **With the prior gone, b2 is not needed at all.**  b2 = 0 and one parameter
   fewer now gives chi2 1079.9 vs 1083.7 for amp2026 WITH b2.  Every earlier
   statement that "dropping b2 costs 71" was measuring the prior's pull on
   R_ET, not the data.
3. **db = 4 was not a parameterisation artefact.**  Giving Ebar_T^d its own free
   N, b, b', nQ (VPK's proposal) changes nothing: the fit still wants a steep d,
   slope 5.1 against 1.05 for u, and gains only 5 chi2.  It is the data.
4. **R_ET at t = 0 is meaningless** - t = 0 lies below |tmin|.  At the physical
   point -t = 0.3 the same fits give Ebar_T^d/Ebar_T^u = 1.2-1.4, not 4-5.
   All ratios are now quoted at xB = 0.25, Q2 = 2.2, -t = 0.3.
5. **F(xB) = xB^alpha (1-xB)^n**: the (1-xB)^n factor does nothing (n -> 0), the
   xB^alpha factor buys 7 chi2 for one parameter.  Split u/d shapes buy 14 for
   four.  It also makes the normalisations degenerate with alpha - exactly what
   VPK remembered from the earlier attempts.

| variant | par | chi2 | chi2/ndf | ET_d/ET_u at -t=0.3 | sigTT(n)/sigTT(p) |
|---|---|---|---|---|---|
| amp2026 (b2 + prior) | 27 | 1083.7 | 1.5504 | - | 0.59 |
| amp2026s as installed (prior, d>=u) | 26 | 1159.0 | 1.6557 | - | 0.9 |
| no prior, tied ET_d, b2=0 | 26 | 1079.9 | 1.5427 | 1.20 | 1.21 |
| + independent ET_d block | 28 | 1075.0 | 1.5401 | 1.20 | 1.42 |
| + F(xB) split | 32 | 1060.7 | 1.5284 | 1.42 | 0.85 |
| **+ Hall-A + phase gauge (CANDIDATE)** | **31** | **1069.8** | **1.5392** | **1.44** | **0.45** |
| no F, + Hall-A + gauge | 27 | 1117.6 | 1.5989 | 0.68 | 0.60 |
| no F, gauge only | 27 | 1075.0 | 1.5379 | 1.20 | 1.42 |

Two more findings:
* **phi_CE is pure gauge.**  Fixing it at 0 costs exactly 0.0 chi2 and removes the
  4e4 errors on the three phases: only combinations of them are observable.
* **eg1-dvcs is what forces a large d.**  Between the "no d" and "d spike"
  branches the eg1 chi2 goes 120.9 -> 62.2, and it is the DOUBLE-SPIN moments
  that do it: A_LL^const 39.9 -> 18.6 and A_LL^cos 33.2 -> 9.7.  A_LL^const
  ~ Re(T01* U01) is the direct H_T x Ebar_T interference.
  The Hall-A neutron ratio pulls the opposite way.  With F(xB) free, satisfying
  Hall-A costs only 9 chi2; without it, 43.
* **The De Masi BSA carries 245 of the 1070.**  Its pull rms is 1.95 and its
  point-to-point scatter inside a bin is twice the quoted (figure-extracted)
  errors, so most of that is underestimated errors, not model failure.  Doubling
  those errors changes no ratio by more than 0.03.

CANDIDATE for production: `fitpar_n2_final.npy` (31 par).  Figures and the
anchored parameter table: OneDrive Work/2026_pi0_amplitudes/candidate_2026_08_29/,
repo figures/*_cand31.png, params_cand31.md.

STILL OPEN, in order of how much they bother me:
1. Ebar_T^d slope 5.11 +- 0.50 against 1.09 +- 0.04 for u - eight sigma apart,
   and it means the d distribution is twice as wide in impact parameter.  Nothing
   in the data forbids it and nothing in the model explains it.
2. sigma_TT(n)/sigma_TT(p) = 0.45 against Hall-A 0.28 +- 0.07 even WITH that
   datum in the fit: still 2.4 sigma, and the model cannot go lower without
   breaking eg1.
3. The normalisations at t = 0 are degenerate with alpha; only the anchored
   values mean anything.  A proper fix is to fit in anchored variables.

## NEXT STEPS, in order

1. Re-check the RC fixed point with amp2026s (`./iterate.sh i3 fitpar_amp2026s.npy`).
   The model moved a lot in eta at large |t|, so the fixed point is not obviously
   still there.  NOTE iterate.sh calls the point-evaluated fit_clas12.py, not the
   bin-averaged fit_binned.py - decide which before running it.
2. R_ET is pulled +4.3 sigma.  Either the forward-limit prior is wrong for this
   parameterisation, or Ebar_T needs the shape freedom that b2 used to provide in
   a form that cannot reverse (dipole).  This is the first thing to settle.
2. Push VPK's global-fit GPDs through our hard kernel and compare convolution to
   convolution -- the only way to settle the inverted slope ordering.  Now sharper:
   we know the chi2 cost of adopting his ordering outright is only ~25.
4. Rewrite ~/pi0eta-rc-2026-paper with the self-consistent result.

## Operational lesson

After pkill on a script name, ALSO kill "spawn_main": multiprocessing workers on
macOS carry only `from multiprocessing.spawn import spawn_main` in their command
line, survive the parent, and keep burning CPU.  One such orphan ran 6h42m at
99% on phallbvpk-mac after a mistaken run there was "stopped".
