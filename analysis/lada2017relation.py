"""
The Lada+ 2017 fit to the California cloud
"""
from astropy import units as u

kappa_california = 0.36 # +/- 0.09 stars pc^-2 mag^-3.31
beta_california = 3.31 # +/- 0.23
AK0_california = 0.51 # mag
kappa_orionA = 1.18 # +/- 0.09, Lombardi+ 2014
beta_orionA = 1.99 # +/- 0.05
AK0_orionA = 0.050 # 0.028
kappa_orionB = 0.60 # +/- 0.09, Lombardi+ 2014
beta_orionB = 2.16 # +/- 0.10
AK0_orionB = 0.20 # +/- 0.10

# uses mu=1.37, alpha_K = 1.67e22 cm^-2/mag
mag_to_col = 183*u.M_sun/u.pc**2/u.mag

# H is the heaviside function; off at <0, on at >0
# Sigma_* = kappa A_K^beta H(A_K-A_K0)

def sigma_star(surfdens, AK0=AK0_california, beta=beta_california, kappa=kappa_california):
    A_K = (surfdens / mag_to_col).to(u.mag).value
    return kappa * A_K**beta * (A_K > AK0)

def sigma_star_california(*args, **kwargs):
    return sigma_star(*args, **kwargs)

def sigma_star_orionA(surfdens):
    return sigma_star(surfdens, AK0=AK0_orionA, beta=beta_orionA, kappa=kappa_orionA)

def sigma_star_orionB(surfdens):
    return sigma_star(surfdens, AK0=AK0_orionB, beta=beta_orionB, kappa=kappa_orionB)
