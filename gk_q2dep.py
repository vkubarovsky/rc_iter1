"""Q2 dependence of the GK_2024 amplitudes at FIXED xB and t.

Uses the HallA tables that exist at xB = 0.360 for Q2 = 1.5 ... 4.44, plus the
pi0-off-NEUTRON table at Q2 = 1.75, which separates u from d without any eta
mixing factor.

Our parameterisation writes each flavour GFF as N Q^nQ, and the amplitude
carries sqrt(A) ~ 1/Q^2, so a structure function goes as Q^(2 nQ - 4).
The exponent measured here is converted back to nQ on that convention.
"""
import json, subprocess, sys, math
files=[("p",1.500,3.35),("p",1.750,4.46),("p",2.000,5.75),("p",3.110,7.38),("p",3.570,8.52),("p",4.440,10.59)]
