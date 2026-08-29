"""Does the t-dependence of a GPD survive the convolution, or is it distorted?

GPD (Regge profile, the standard GK form):
    F(x, t) = N x^-a0 (1-x)^beta exp[ t (b - alpha' ln x) ]
Convolution (LO handbag, valence combination):
    <F>(xi, t) = int_0^1 dx F(x,t) [ 1/(x - xi - i eps) + 1/(x + xi) ]

Three things are measured:
  1. the effective slope of |<F>| near t = 0 and averaged over -t < 1;
  2. whether the convolution is still a pure exponential (fit b1 t + b2 t^2);
  3. Re and Im separately -- Im samples x = xi exactly, Re averages over all x.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit

def F(x, N, a0, beta, b, ap, t):
    return N*x**(-a0)*(1-x)**beta*np.exp(t*(b - ap*np.log(x)))
def conv(xi, t, N=1.0, a0=0.30, beta=4.0, b=0.50, ap=0.45):
    f=lambda x: F(x,N,a0,beta,b,ap,t); fx=f(xi)
    pv=quad(lambda x:(f(x)-fx)/(x-xi),0,1,points=[xi],limit=300)[0]+fx*np.log((1-xi)/xi)
    reg=quad(lambda x: f(x)/(x+xi),0,1,limit=300)[0]
    return complex(pv+reg, np.pi*fx)
def slopes(xi, **kw):
    T=np.linspace(0.0,1.0,21)
    A=np.array([abs(conv(xi,-mt,**kw)) for mt in T])
    R=np.array([conv(xi,-mt,**kw).real for mt in T])
    I=np.array([conv(xi,-mt,**kw).imag for mt in T])
    g=lambda mt,c,b1,b2: c-b1*mt+b2*mt*mt
    p,_=curve_fit(g,T,np.log(A),p0=[np.log(A[0]),1.0,0.0])
    b0=(np.log(A[0])-np.log(A[1]))/(T[1]-T[0])                 # local slope at t->0
    bR=(np.log(abs(R[0]))-np.log(abs(R[1])))/(T[1]-T[0])
    bI=(np.log(I[0])-np.log(I[1]))/(T[1]-T[0])
    return b0,p[1],p[2],bR,bI

print("GPD input: x^-0.3 (1-x)^4, slope b - alpha' ln x   (alpha' = 0.45, b = 0.50)\n")
print(f"{'xB':>5} {'xi':>6} | {'b(GPD) at x=xi':>14} | {'b_eff(conv) t->0':>17} {'b1 over -t<1':>13} "
      f"{'curvature b2':>13} | {'b(Re)':>7} {'b(Im)':>7}")
for xB in (0.13,0.17,0.22,0.27,0.34,0.43,0.52):
    xi=xB/(2-xB)
    bgpd=0.50-0.45*np.log(xi)
    b0,b1,b2,bR,bI=slopes(xi)
    print(f"{xB:5.2f} {xi:6.3f} | {bgpd:14.2f} | {b0:17.2f} {b1:13.2f} {b2:13.3f} | {bR:7.2f} {bI:7.2f}")

print("\ncontrol: alpha' = 0, i.e. the GPD slope does not depend on x")
print(f"{'xB':>5} {'b(GPD)':>8} {'b_eff(conv)':>12} {'curvature b2':>13}")
for xB in (0.17,0.34,0.52):
    xi=xB/(2-xB); b0,b1,b2,_,_=slopes(xi,ap=0.0)
    print(f"{xB:5.2f} {0.50:8.2f} {b0:12.3f} {b2:13.4f}")

print("\nhow the distortion grows with alpha'  (xB = 0.27)")
xi=0.27/(2-0.27)
print(f"{'alphaP':>7} {'b(GPD) at xi':>13} {'b_eff(conv)':>12} {'shift':>8} {'b2':>8}")
for ap in (0.0,0.15,0.3,0.45,0.6,1.0):
    b0,b1,b2,_,_=slopes(xi,ap=ap)
    print(f"{ap:7.2f} {0.50-ap*np.log(xi):13.2f} {b0:12.2f} {b0-(0.50-ap*np.log(xi)):8.2f} {b2:8.3f}")
