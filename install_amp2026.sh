#!/bin/bash
# Freeze the converged iteration as amp2026 and register it in exclurad_py.
set -e
PAR=${1:-/Users/vpk/rc_iter1/fitpar_i2.npy}
cp $PAR ~/exclurad_py/exclurad_py/models/amp2026_par.npy
cd ~/exclurad_py
~/.venv/bin/python3 - <<'PY'
p = "exclurad_py/models/pseudoscalar.py"
s = open(p).read()
if "amp2026" in s:
    print("already present"); raise SystemExit
block = '''

# --------------------------------------------------------------------------- #
#  amp2026: the self-consistent fixed point of the RC iteration
# --------------------------------------------------------------------------- #
#
# Same amplitude parameterisation as amp2021, but fitted to CLAS6 data whose
# radiative correction was recomputed WITH the amplitude model itself, iterated
# to a fixed point (~/rc_iter1).  amp2021 is fitted to the published data, whose
# RC assumed a sigma_L-dominant model; the July 2026 study
# (~/pi0eta-rc-2026-paper) showed the RC depends significantly on that
# assumption.  amp2026 removes the assumption: sigma_L is measured, not guessed.
#
# Use amp2026 for production.  Keep amp2021 to reproduce the published data.

_AMP26_PAR = _np.load(_os.path.join(_os.path.dirname(__file__), "amp2026_par.npy"))
_AMP26_PROV = (
    "Self-consistent fixed point of the RC iteration: RC recomputed with the "
    "amplitude model over the full published grids, cross sections re-corrected "
    "(sigma_new = sigma_pub / r, r = eta(new)/eta(published), both from "
    "exclurad_py so the code cancels), sigma(phi) refitted, amplitudes refitted. "
    "See ~/rc_iter1/README.md. ")


def _amp26_sfs(channel, t_nucl, xB, Q2):
    t = _np.atleast_1d(_np.asarray(t_nucl, dtype=float))
    x = _np.atleast_1d(_np.asarray(xB, dtype=float))
    q = _np.atleast_1d(_np.asarray(Q2, dtype=float))
    t, x, q = _np.broadcast_arrays(t, x, q)
    out = _np.zeros((5,) + t.shape)
    o = out.reshape(5, -1)
    for i, (ti, xi, qi) in enumerate(zip(t.ravel(), x.ravel(), q.ravel())):
        s = _AF.structure(_AMP26_PAR, channel, float(ti), float(xi), float(qi))
        if s is None:
            continue
        o[0, i], o[1, i] = s["T"], s["L"]
        o[2, i], o[3, i], o[4, i] = s["TT"], s["LT"], s["LTp"]
    if out.shape[1:] == (1,):
        return tuple(float(v) for v in out[:, 0])
    return tuple(out)


@register
class Pi0_Amp2026(ProtonDetectedModel):
    name = "pi0.amp2026"
    channel = "pi0"
    version = "amp2026"
    topology = _PI0_TOPO
    validity = _AMP_VALIDITY
    provenance = _AMP26_PROV + "pi0 off the proton."
    reference = _AMP_REF

    def _sfs_tn(self, t_nucl, xB, Q2):
        return _amp26_sfs("pi0p", t_nucl, xB, Q2)


@register
class Eta_Amp2026(ProtonDetectedModel):
    name = "eta.amp2026"
    channel = "eta"
    version = "amp2026"
    topology = _ETA_TOPO
    validity = _AMP_VALIDITY
    provenance = _AMP26_PROV + "eta off the proton, flavour-rotated."
    reference = _AMP_REF

    def _sfs_tn(self, t_nucl, xB, Q2):
        return _amp26_sfs("etap", t_nucl, xB, Q2)


@register
class Pi0n_Amp2026(ProtonDetectedModel):
    name = "pi0n.amp2026"
    channel = "pi0"
    version = "amp2026-n"
    topology = _PI0_TOPO
    validity = _AMP_VALIDITY
    provenance = _AMP26_PROV + "pi0 off the NEUTRON. Prediction."
    reference = _AMP_REF

    def _sfs_tn(self, t_nucl, xB, Q2):
        return _amp26_sfs("pi0n", t_nucl, xB, Q2)


@register
class Etan_Amp2026(ProtonDetectedModel):
    name = "etan.amp2026"
    channel = "eta"
    version = "amp2026-n"
    topology = _ETA_TOPO
    validity = _AMP_VALIDITY
    provenance = _AMP26_PROV + "eta off the NEUTRON. Prediction."
    reference = _AMP_REF

    def _sfs_tn(self, t_nucl, xB, Q2):
        return _amp26_sfs("etan", t_nucl, xB, Q2)
'''
open(p, "a").write(block)
print("amp2026 models appended")
PY
PYTHONPATH=~/exclurad_py ~/.venv/bin/python3 -c "
from exclurad_py.models.registry import names
print('registered:', sorted(n for n in names() if 'amp2026' in n))"
