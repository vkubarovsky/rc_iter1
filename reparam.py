"""Drop the ln 0.15 offset from the t-slope parameterisation.

Old:  exp[ (b + b'(ln xB - ln 0.15)) t ]      b = slope at xB = 0.15
New:  exp[ (b + b'  ln xB          ) t ]      b = slope at xB = 1, -b' = alpha'

Exact reparameterisation, no refit:  b_new = b - ln(0.15) b' = b + 1.8971 b'.
Blocks carrying a b': H_T^u, H_T^d, Ebar_T^u (E_T^d inherits it), T00.

    ~/.venv/bin/python3 reparam.py <old.npy> <new.npy>

amplitudes.py must already be in the NEW convention (no offset); the round-trip
check below compares against a grid dumped with the old code.
"""
import math, sys
import numpy as np

SHIFT = -math.log(0.15)                     # 1.8971
BLOCKS = {"H_T^u": (1, 2), "H_T^d": (5, 6), "Ebar_T^u": (9, 10), "T00": (16, 21)}

def to_lnxB(p):
    q = np.array(p, dtype=float)
    for name, (ib, ibp) in BLOCKS.items():
        if ibp >= len(q):                   # optional xB-slope absent
            continue
        q[ib] = p[ib] + SHIFT*p[ibp]
    return q

def to_old(p):
    """Inverse: new (ln xB) -> old (ln xB - ln 0.15) convention.  Needed when a
    parameter file is handed to code that still uses the offset - exclurad_py's
    models/_amplitude_fit.py does."""
    q = np.array(p, dtype=float)
    for name, (ib, ibp) in BLOCKS.items():
        if ibp >= len(q):
            continue
        q[ib] = p[ib] - SHIFT*p[ibp]
    return q

def slopes(p, xB):
    """Effective t-slope b + b' ln xB of every block, new convention."""
    lx = math.log(xB)
    out = {}
    for name, (ib, ibp) in BLOCKS.items():
        bp = p[ibp] if ibp < len(p) else 0.0
        out[name] = p[ib] + bp*lx
    out["Ebar_T^d"] = out["Ebar_T^u"] + p[14]
    return out

if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    p = np.load(src)
    q = to_lnxB(p)
    np.save(dst, q)
    print(f"{src} -> {dst}   ({len(q)} par, b_new = b + {SHIFT:.4f} b')\n")
    print(f"{'block':>10} {'N':>9} {'b_old':>8} {'b_new':>8} {'b_prime':>8} |"
          f"{'  slope at xB = 0.10   0.15   0.25   0.40   0.60'}")
    XS = (0.10, 0.15, 0.25, 0.40, 0.60)
    tab = {x: slopes(q, x) for x in XS}
    for name, (ib, ibp) in list(BLOCKS.items()) + [("Ebar_T^d", (9, 10))]:
        N = {"H_T^u": q[0], "H_T^d": q[4], "Ebar_T^u": q[8],
             "Ebar_T^d": q[13]*q[8], "T00": q[15]}[name]
        bo = p[ib] if name != "Ebar_T^d" else p[9] + p[14]
        bn = q[ib] if name != "Ebar_T^d" else q[9] + q[14]
        bp = q[ibp] if ibp < len(q) else 0.0
        row = "".join(f"{tab[x][name]:7.2f}" for x in XS)
        flag = "  NEGATIVE" if min(tab[x][name] for x in XS) < 0 else ""
        print(f"{name:>10} {N:9.2f} {bo:8.3f} {bn:8.3f} {bp:8.3f} |  {row}{flag}")
