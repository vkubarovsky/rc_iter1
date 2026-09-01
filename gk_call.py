"""Thin wrapper around libGKPi0.  Runs under python3.13 only -- the module is
built for that ABI, while our venv is 3.14.

    /opt/homebrew/bin/python3.13 gk_call.py <model> <Q2> <xB> <Ebeam> <t1,t2,...>

prints one line per t:  t sigma0 sigmaT sigmaL sigmaTT sigmaLT
"""
import os, sys, time
sys.path.insert(0, "/Users/vpk/gk_work/build")
import gkpi0

PREP = os.path.expanduser("~/rc_iter1/gkprep")

def switches():
    gkpi0.set_reaction(0)          # pi0 on the proton
    gkpi0.set_tmin_peter(False)
    gkpi0.set_xi_peter(False)
    gkpi0.set_eta_mixing(False)
    gkpi0.set_mu_eta(0.0)

def sf(model, Q2, xB, E, ts):
    switches()
    stem = f"{model}_pi0p_Q{Q2:.4f}_x{xB:.5f}_E{E:.3f}.dat"
    gkpi0.prepare_or_load(float(Q2), float(xB), float(E), model,
                          os.path.join(PREP, stem))
    return gkpi0.get_structure_functions_grid(float(Q2), float(xB),
              [-abs(t) for t in ts], 0.0, float(E), model)

if __name__ == "__main__":
    model, Q2, xB, E = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    ts = [float(z) for z in sys.argv[5].split(",")]
    t0 = time.time()
    r = sf(model, Q2, xB, E, ts)
    keys = list(r.keys())
    print("# keys:", keys, file=sys.stderr)
    print(f"# elapsed {time.time()-t0:.1f} s", file=sys.stderr)
    for i, t in enumerate(ts):
        print(t, " ".join(f"{r[k][i]:.6g}" for k in keys))
