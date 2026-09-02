"""All plots for one fit run.

    python3 plots.py <tag>

Every figure says in its title whether the set was in the fit, so a run and its
blind predictions can be told apart at a glance.
"""
import json, math, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as amp, datasets as D

C_IN, C_OUT, C_MOD = "#123a5e", "#9a4f2b", "#0d6a72"
OBSLAB = {"U": r"$\sigma_U=\sigma_T+\epsilon\sigma_L$", "T": r"$\sigma_T$",
          "LT": r"$\sigma_{LT}$", "TT": r"$\sigma_{TT}$", "LTp": r"$\sigma_{LT'}$",
          "A_LU^sinphi": r"$A_{LU}^{\sin\phi}$",
          "sigma_LT'/sigma_0": r"$\sigma_{LT'}/\sigma_0$",
          "AULsin": r"$A_{UL}^{\sin\phi}$", "AULsin2": r"$A_{UL}^{\sin 2\phi}$",
          "ALLc": r"$A_{LL}^{\rm const}$", "ALLcos": r"$A_{LL}^{\cos\phi}$"}

def cluster(rows, dq=0.25, dx=0.03):
    """Group rows into kinematic settings.

    The cross-section sets carry the group label worked out when the database
    was built -- one label per (Q2, xB) setting -- so use it.  Only the
    asymmetries, which have no label, fall back on clustering by proximity."""
    lab = {r[7] for r in rows}
    if lab and "all" not in lab:
        return [ [r for r in rows if r[7] == k]
                 for k in sorted(lab, key=lambda k: ([r for r in rows if r[7]==k][0][1],
                                                     [r for r in rows if r[7]==k][0][0])) ]
    g = []
    for r in rows:
        for k in g:
            if abs(k[0][0]-r[0]) < dq and abs(k[0][1]-r[1]) < dx:
                k.append(r); break
        else:
            g.append([r])
    return sorted(g, key=lambda k: (k[0][1], k[0][0]))

def grid(n):
    c = min(5, max(1, int(math.ceil(math.sqrt(n)))))
    if n > 12: c = 5
    if n > 25: c = 6
    return int(math.ceil(n/c)), c

def panels(key, obs, p, fitted, out, ylog=False):
    s = D.BY[key]
    rows = [r for r in s["rows"] if r[3] == obs]
    if not rows: return None
    gs = cluster(rows)
    nr, nc = grid(len(gs))
    fig, ax = plt.subplots(nr, nc, figsize=(3.2*nc, 2.5*nr), squeeze=False)
    for a in ax.flat: a.set_visible(False)
    col = C_IN if fitted else C_OUT
    tot = ntot = 0
    for i, g in enumerate(gs):
        a = ax.flat[i]; a.set_visible(True)
        # a projection must be drawn against the variable it scans: the COMPASS
        # Q2 and nu projections all sit at nearly the same <|t|>, so plotting them
        # against |t| makes a nu dependence look like a t dependence
        gl = str(g[0][7])
        xv, axlab = ([r[2] for r in g], r"$-t$  [GeV$^2$]")
        if gl == "compass_Q2": xv, axlab = ([r[0] for r in g], r"$Q^2$  [GeV$^2$]")
        elif gl == "compass_nu":
            xv, axlab = ([r[0]/(2*0.9382720813*r[1]) for r in g], r"$\nu$  [GeV]")
        t = np.array(xv); v = np.array([r[4] for r in g])
        e = np.array([r[5] for r in g])
        o = np.argsort(t); t, v, e = t[o], v[o], e[o]
        a.errorbar(t, v, e, fmt='o', ms=4, color=col, ecolor=col, capsize=2, zorder=3)
        m = np.array([s["predict"](p, r) if s["predict"](p, r) is not None else np.nan
                      for r in [g[j] for j in o]], dtype=float)
        good = np.isfinite(m)
        if good.sum() > 1 and key != "compass":
            a.plot(t[good], m[good], '-', lw=1.8, color=C_MOD, zorder=2)
        a.plot(t[good], m[good], 's', ms=5 if key == "compass" else 3,
               color=C_MOD, zorder=2)
        c = float(np.nansum(((v-m)/e)**2)); n = int(good.sum())
        tot += c; ntot += n
        # a panel whose kinematics vary inside the group must say so, not quote
        # the first row and pretend the rest sit there too
        qs = [r[0] for r in g]; xs = [r[1] for r in g]
        qlab = (f"$Q^2$={qs[0]:.2f}" if max(qs)-min(qs) < 0.05*qs[0]
                else f"$Q^2$={min(qs):.2f}-{max(qs):.2f}")
        xlab = (f"$x_B$={xs[0]:.3f}" if max(xs)-min(xs) < 0.05*xs[0]
                else f"$x_B$={min(xs):.3f}-{max(xs):.3f}")
        a.text(0.97, 0.95, f"{qlab}\n{xlab}\n$\\chi^2$={c:.1f}/{n}",
               transform=a.transAxes, ha='right', va='top', fontsize=7.5)
        a.grid(alpha=.25, lw=.5); a.tick_params(labelsize=8)
        if ylog and (v > 0).all(): a.set_yscale('log')
        a.set_xlabel(axlab, fontsize=9)
        if i % nc == 0: a.set_ylabel(OBSLAB.get(obs, obs), fontsize=9)
    tag = "IN THE FIT" if fitted else "NOT in the fit  (blind prediction)"
    fig.suptitle(f"{s['label']}   —   {OBSLAB.get(obs,obs)}   —   {tag}\n"
                 f"$\\chi^2$ = {tot:.1f} / {ntot} = {tot/max(ntot,1):.2f} per point",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 1-0.055*(2 if nr < 3 else 1)])
    fig.savefig(out, dpi=130); plt.close(fig)
    return tot, ntot

def summary_figure(rec, out):
    ks = [k for k in D.ALL]
    v = [rec["sets"][k]["chi2"]/max(rec["sets"][k]["n"], 1) for k in ks]
    f = [rec["sets"][k]["fitted"] for k in ks]
    fig, a = plt.subplots(figsize=(10, 4.6))
    b = a.bar(range(len(ks)), v, color=[C_IN if x else C_OUT for x in f])
    for i, k in enumerate(ks):
        a.text(i, v[i], f" {rec['sets'][k]['chi2']:.0f}/{rec['sets'][k]['n']}",
               rotation=90, ha='center', va='bottom', fontsize=8)
    a.axhline(1, ls='--', lw=1, color='#666666')
    a.set_xticks(range(len(ks))); a.set_xticklabels(ks, rotation=35, ha='right', fontsize=9)
    a.set_ylabel(r"$\chi^2$ per point"); a.set_yscale('log')
    a.set_ylim(0.2, max(v)*3)
    a.grid(axis='y', alpha=.3, lw=.5)
    a.set_title(f"run {rec['tag']}   —   fitted $\\chi^2$/ndf = {rec['chi2_ndf']:.3f}"
                f"  ({rec['chi2_fitted']:.1f}/{rec['ndf']})\n"
                f"dark = in the fit,  brown = left out (blind prediction)", fontsize=11)
    fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)

PLAN = [("clas6_pi0", ("U","LT","TT")), ("clas6_eta", ("U","LT","TT")),
        ("halla_n", ("T","LT","TT")), ("halla_y16", ("T","LT","TT")),
        ("halla_y11", ("U","LT","TT","LTp")), ("halla_y21", ("U","LT","TT","LTp")),
        ("clas12_xs", ("U","LT","TT")), ("compass", ("U","TT")),
        ("bsa_demasi", ("A_LU^sinphi",)), ("bsa_zhao", ("A_LU^sinphi",)),
        ("bsa_clas12", ("sigma_LT'/sigma_0",)),
        ("eg1", ("AULsin","AULsin2","ALLc","ALLcos"))]
ORDER = {k: i for i, (k, _) in enumerate(PLAN)}

if __name__ == "__main__":
    tag = sys.argv[1]
    d = f"runs/{tag}"
    p = np.load(f"{d}/fitpar.npy"); rec = json.load(open(f"{d}/summary.json"))
    os.makedirs(f"{d}/plots", exist_ok=True)
    summary_figure(rec, f"{d}/plots/00_chi2_summary.png")
    for key, obs in PLAN:
        fitted = rec["sets"][key]["fitted"]
        mark = "in" if fitted else "out"
        for j, o in enumerate(obs):
            name = o.replace("/", "_over_").replace("'", "p")
            f = f"{d}/plots/{ORDER[key]+1:02d}_{key}_{mark}_{name}.png"
            r = panels(key, o, p, fitted, f)
            if r: print(f"   {os.path.basename(f):48s} chi2 {r[0]:8.1f}/{r[1]}")
    print(f"plots -> {d}/plots")
