#!/bin/bash
# Register the slope-constrained amplitude model as exclurad_py's amp2026s.
# amp2026 itself is left frozen: the OneDrive generator baseline was produced
# with it.  amp2026s is the production model from now on.
set -e
PAR=${1:-/Users/vpk/rc_iter1/fitpar_prod.npy}
PY=~/.venv/bin/python3
# exclurad_py/models/_amplitude_fit.py still uses the ln 0.15 offset -> convert
$PY -c "import numpy as np,sys; sys.path.insert(0,'$HOME/rc_iter1'); import reparam; \
np.save('$HOME/exclurad_py/exclurad_py/models/amp2026s_par.npy', reparam.to_old(np.load('$PAR')))"
cd ~/exclurad_py
$PY - <<'PY'
p = "exclurad_py/models/pseudoscalar.py"
s = open(p).read()
if "amp2026s" in s:
    print("already present"); raise SystemExit
block = '''

# --------------------------------------------------------------------------- #
#  amp2026s: amp2026 with the physical t-slope constraints imposed
# --------------------------------------------------------------------------- #
#
# Same data and same RC fixed point as amp2026; what changed is that the fit is
# no longer allowed to produce form factors that GROW with |t|.  amp2026 had a
# negative H_T^d slope over the whole measured xB range.  Imposed here:
#   b + b' ln xB >= 0 for every block on xB in [0.1, 0.6];
#   slope(H_T^d) >= slope(H_T^u), the ordering of the global GPD fits;
#   Ebar_T slope including the b2 t^2 curvature >= 0 out to -t = 2.5, so the
#   curvature cannot reverse the fall inside the validity window.
# Cost: chi2 1083.7 -> 1108.1 (chi2/ndf 1.5504 -> 1.5853) for zero new parameters,
# and all of it sits in eta sigma_U at -t > 1.2; inside -t <= 1.2 the ordering
# costs 6 chi2 on 513 points.  See ~/rc_iter1/README.md.
#
# Use amp2026s for production.  amp2026 is kept frozen because the OneDrive
# generator baseline (Work/2026_pi0_amplitudes/generator/) was produced with it.

_AMP26S_PAR = _np.load(_os.path.join(_os.path.dirname(__file__), "amp2026s_par.npy"))
_AMP26S_PROV = (
    "Slope-constrained refit of the amp2026 fixed point: form factors required "
    "to fall with |t| over xB in [0.1, 0.6], H_T^d at least as steep as H_T^u, "
    "Ebar_T curvature not reversing the fall inside -t <= 2.5. "
    "See ~/rc_iter1/README.md. ")


def _amp26s_sfs(channel, t_nucl, xB, Q2):
    t = _np.atleast_1d(_np.asarray(t_nucl, dtype=float))
    x = _np.atleast_1d(_np.asarray(xB, dtype=float))
    q = _np.atleast_1d(_np.asarray(Q2, dtype=float))
    t, x, q = _np.broadcast_arrays(t, x, q)
    out = _np.zeros((5,) + t.shape)
    o = out.reshape(5, -1)
    for i, (ti, xi, qi) in enumerate(zip(t.ravel(), x.ravel(), q.ravel())):
        s = _AF.structure(_AMP26S_PAR, channel, float(ti), float(xi), float(qi))
        if s is None:
            continue
        o[0, i], o[1, i] = s["T"], s["L"]
        o[2, i], o[3, i], o[4, i] = s["TT"], s["LT"], s["LTp"]
    if out.shape[1:] == (1,):
        return tuple(float(v) for v in out[:, 0])
    return tuple(out)


@register
class Pi0_Amp2026s(ProtonDetectedModel):
    name = "pi0.amp2026s"
    channel = "pi0"
    version = "amp2026s"
    topology = _PI0_TOPO
    validity = _AMP_VALIDITY
    provenance = _AMP26S_PROV + "pi0 off the proton."
    reference = _AMP_REF

    def _sfs_tn(self, t_nucl, xB, Q2):
        return _amp26s_sfs("pi0p", t_nucl, xB, Q2)


@register
class Eta_Amp2026s(ProtonDetectedModel):
    name = "eta.amp2026s"
    channel = "eta"
    version = "amp2026s"
    topology = _ETA_TOPO
    validity = _AMP_VALIDITY
    provenance = _AMP26S_PROV + "eta off the proton, flavour-rotated."
    reference = _AMP_REF

    def _sfs_tn(self, t_nucl, xB, Q2):
        return _amp26s_sfs("etap", t_nucl, xB, Q2)


@register
class Pi0n_Amp2026s(ProtonDetectedModel):
    name = "pi0n.amp2026s"
    channel = "pi0"
    version = "amp2026s-n"
    topology = _PI0_TOPO
    validity = _AMP_VALIDITY
    provenance = _AMP26S_PROV + "pi0 off the NEUTRON. Prediction."
    reference = _AMP_REF

    def _sfs_tn(self, t_nucl, xB, Q2):
        return _amp26s_sfs("pi0n", t_nucl, xB, Q2)
'''
open(p, "a").write(block)
print("amp2026s registered in", p)
PY
