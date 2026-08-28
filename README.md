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
