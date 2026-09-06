"""Generate every number in the report from the run, so nothing is transcribed."""
import json, os, numpy as np, fitrun as R, datasets as D

TAG = "C0_with_compass"
d = f"runs/{TAG}"
p = np.load(f"{d}/fitpar.npy"); rec = json.load(open(f"{d}/summary.json"))
keys = rec["fitted"]; lamd = rec["norms"]; nn = [k for k in keys if k in R.NORMS]
nf = len(R.FREE)

def resid(x):
    q = R.expand(x[:nf]); lam = dict(zip(nn, x[nf:])); r = []
    for k in keys:
        s = D.BY[k]; L = lam.get(k, 1.0)
        for row in s["rows"]:
            v = s["predict"](q, row)
            r.append((row[4] - L*v)/row[5] if v is not None else 5.0)
    for k in nn: r.append((lam[k] - 1.0)/R.NORMS[k])
    for xx in (0.05, 0.35, 1.00):
        sl = R.slopes(q, xx)
        for n in R.BL: r.append(R.W*min(0.0, sl[n]))
    return np.array(r)

x0 = np.concatenate([p[R.FREE], [lamd[k] for k in nn]]); f0 = resid(x0)
J = np.zeros((len(f0), len(x0)))
for i in range(len(x0)):
    h = 1e-6*max(abs(x0[i]), 1e-3); xp = x0.copy(); xp[i] += h
    J[:, i] = (resid(xp) - f0)/h
err = np.sqrt(np.abs(np.diag(np.linalg.pinv(J.T @ J))))
E = np.zeros(37)
for j, i in enumerate(R.FREE): E[i] = err[j]
for t, s in R.TIES.items(): E[t] = E[s]
LO, HI = R.bounds()

def cell(i, tie=None):
    if tie: return tie
    at = "$^{\\dagger}$" if (abs(p[i]-LO[i]) < 1e-6 or abs(p[i]-HI[i]) < 1e-6) else ""
    return f"${p[i]:.3f} \\pm {E[i]:.3f}${at}"

BLK = [("$H_T^u$", 0, 1, 2, 3), ("$H_T^d$", 4, 5, 6, 7),
       ("$\\bar E_T^u$", 8, 9, 10, 11), ("$\\bar E_T^d$", 13, 14, 12, 27),
       ("$L$ (T00)", 15, 16, 21, 17)]
TIED = {5: "$=b_u$", 6: "$=b_u'$", 7: "$=n_{Q,u}$",
        12: "$=b_u'(\\bar E_T)$", 27: "$=n_{Q,u}(\\bar E_T)$"}
# Whole tabulars, not bodies: \input of a body inside a tabular puts
# \bottomrule in the middle of a row.
def table(fn, spec, head, rows):
    with open(fn, "w") as f:
        f.write("\\begin{tabular}{" + spec + "}\n\\toprule\n"
                + head + " \\\\\n\\midrule\n")
        f.write(" \\\\\n".join(rows) + " \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")

table("tex/tab_blocks.tex", "lcccc",
      "block & $N$ & $b$ & $b'$ & $n_Q$",
      [f"{nm} & " + " & ".join(cell(i, TIED.get(i)) for i in (iN, ib, ibp, inq))
       for nm, iN, ib, ibp, inq in BLK])

OTH = [("$R_L=N_d/N_u$ (T00)", 18), ("$\\delta_0$", 19), ("$\\delta_1$", 20),
       ("$\\rho_{CE}$", 22), ("$\\phi_{CE}$ (frozen)", 23), ("$\\phi_w$", 24),
       ("$\\rho_{nf}$", 25), ("$\\phi_{nf}$", 26),
       ("$\\varphi_{d/u}(H_T)$", 34), ("$\\varphi_{d/u}(\\bar E_T)$", 35),
       ("$\\varphi_{d/u}(T_{00})$", 36)]
table("tex/tab_other.tex", "lccc",
      "parameter & value & lower & upper",
      [f"{nm} & {cell(i)} & {LO[i]:.2f} & {HI[i]:.2f}" for nm, i in OTH])

_nr = []
for j, k in enumerate(nn):
    v = lamd[k]; e = err[nf+j]; sg = R.NORMS[k]
    _nr.append(f"\\texttt{{{k.replace('_', chr(92)+'_')}}} & ${100*sg:.2f}\\,\\%$ & "
               f"${v:.4f} \\pm {e:.4f}$ & ${(v-1)/sg:+.2f}\\,\\sigma$")
table("tex/tab_norms.tex", "lccc", "set & prior & $\\lambda$ & pull", _nr)

LAB = {"clas6_pi0": "CLAS6 $\\pi^0$ p", "clas6_eta": "CLAS6 $\\eta$ p",
       "halla_n": "Hall A 2017, $\\pi^0$ n", "halla_y16": "Hall A 2016, $\\pi^0$ p",
       "halla_y11": "Hall A 2011, $\\pi^0$ p", "halla_y21": "Hall A 2021, $\\pi^0$ p",
       "clas12_xs": "CLAS12 cross sections", "compass": "COMPASS 2025",
       "bsa_demasi": "De Masi, published $\\alpha$", "bsa_zhao": "Zhao, $\\eta$ BSA",
       "bsa_clas12": "CLAS12 $\\sigma_{LT'}/\\sigma_0$",
       "bsa_demasi_phi": "De Masi, $A_{LU}(\\phi)$", "eg1": "eg1-dvcs, pol.\\ target"}
_sr = []
for k in D.ALL:
    st = rec["sets"][k]
    mark = "yes" if st["fitted"] else "--"
    _sr.append(f"{LAB[k]} & \\texttt{{{k.replace('_', chr(92)+'_')}}} & {mark} & "
               f"{st['n']} & {st['chi2']:.1f} & {st['chi2']/max(st['n'],1):.2f}")
table("tex/tab_sets.tex", "llccrr",
      "set & key & fitted & points & $\\chi^2$ & per point", _sr)

with open("tex/tab_head.tex", "w") as f:
    f.write(f"\\newcommand{{\\fittag}}{{{TAG.replace('_','\\_')}}}\n"
            f"\\newcommand{{\\chindf}}{{{rec['chi2_ndf']:.4f}}}\n"
            f"\\newcommand{{\\chitot}}{{{rec['chi2_fitted']:.1f}}}\n"
            f"\\newcommand{{\\ndf}}{{{rec['ndf']}}}\n"
            f"\\newcommand{{\\npar}}{{{rec['npar']}}}\n"
            f"\\newcommand{{\\nfit}}{{{rec['n_fitted']}}}\n")
    for k, v in rec["ratios"].items():
        f.write(f"\\newcommand{{\\r{k.replace('_','')}}}{{{v:.3f}}}\n")
# \input inside a tabular: the trailing "\\" of the last row would put
# \bottomrule inside a row instead of between rows.
for _f in ("tex/tab_sets.tex", "tex/tab_blocks.tex", "tex/tab_other.tex",
           "tex/tab_norms.tex"):
    _s = open(_f).read().rstrip()
    if _s.endswith("\\\\"): _s = _s[:-2]
    open(_f, "w").write(_s + "\n")

print("tables written; npar =", rec["npar"], " ndf =", rec["ndf"])
