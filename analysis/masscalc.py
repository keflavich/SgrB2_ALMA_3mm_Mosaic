import numpy as np
from astropy import units as u
from constants import distance
from astropy import constants

# just hacked together quickly:
# BGPS center is 271.1 GHz
# Assume ALMA center freq is 100 GHz
# use beta = 1.5 for alpha=3.5
onemmto3mm = (271.1*u.GHz/100*u.GHz)**3.5

def mass_conversion_factor(TK=20, d=distance.to(u.kpc).value):
    return 14.30 * (np.exp(13.01/TK) - 1)*d**2 * onemmto3mm

def col_conversion_factor(TK=20):
    return 2.19e22 * (np.exp(13.01/TK - 1)) * onemmto3mm

def Jnu(T, nu):
    return (2*constants.h*nu**3 / constants.c**2 *
            np.exp(constants.h*nu/(constants.k_B*T)-1)**-1)
