"""Relative u-d phase of the convolutions, computed rather than fitted.

The convolution has a pole: Im<F> = pi F(xi,xi,t), Re<F> = principal value.  The
two sample the GPD with different x weights, so the phase is set by the x shape.
Our parameterisation carries only the modulus, so an x shape has to be supplied:
GK's, F(x) = N x^-a0 (1-x)^beta exp[t(b - alpha' ln x)], with beta = 4 for u and
5 for d (hard-wired in EBarU/EBarD of libGKPi0), a0 = -0.10 for Ebar_T and -0.17
for H_T, alpha' = 0.45.

Output: phi_du(xi) = arg<F^d> - arg<F^u>, and the sector-to-sector phase that
sets the longitudinal-transverse interference.
"""
import numpy as np, math
from scipy.integrate import quad
def conv(xi,t,a0,beta,b,ap=0.45,N=1.0):
    f=lambda x: N*x**(-a0)*(1-x)**beta*np.exp(t*(b-ap*np.log(x)))
    fx=f(xi)
    pv=quad(lambda x:(f(x)-fx)/(x-xi),0,1,points=[xi],limit=300)[0]+fx*np.log((1-xi)/xi)
    reg=quad(lambda x: f(x)/(x+xi),0,1,limit=300)[0]
    return complex(pv+reg, np.pi*fx)
XI=np.linspace(0.05,0.42,20)
tab={}
for name,(a0,bu,bd) in (("ET",(-0.10,0.77,0.50)),("HT",(-0.17,0.30,0.30))):
    ph=[]
    for xi in XI:
        u=conv(xi,-0.3,a0,4,bu); d=conv(xi,-0.3,a0,5,bd)
        ph.append(np.angle(d)-np.angle(u))
    tab[name]=np.array(ph)
np.savez("phase_table.npz",xi=XI,ET=tab["ET"],HT=tab["HT"])
print("computed relative u-d phase, at -t = 0.3 (it barely moves with t)")
print(f"{'xi':>6} {'xB':>6} {'phi_du(Ebar_T)':>15} {'phi_du(H_T)':>13}")
for i in range(0,20,3):
    xi=XI[i]; xB=2*xi/(1+xi)
    print(f"{xi:6.3f} {xB:6.3f} {tab['ET'][i]:15.3f} {tab['HT'][i]:13.3f}")
print("\nand the phase of the convolution itself, which sets the L-T interference:")
print(f"{'xi':>6} | {'arg<Ebar_T^u>':>14} {'arg<H_T^u>':>11} {'arg(L) - arg(T)':>16}")
for xi in (0.07,0.14,0.22,0.33):
    et=conv(xi,-0.3,-0.10,4,0.77); ht=conv(xi,-0.3,-0.17,4,0.30)
    lo=conv(xi,-0.3, 0.30,4,0.50)     # a twist-2-like longitudinal shape, alpha0 = +0.3
    print(f"{xi:6.3f} | {np.angle(et):14.3f} {np.angle(ht):11.3f} {np.angle(lo)-np.angle(ht):16.3f}")
print("\nour fit has delta(t) = 2.001 + 0.559 t, i.e. 2.00 rad = 115 deg at t = 0")
