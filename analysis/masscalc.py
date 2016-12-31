import numpy as np
from astropy import units as u
from constants import distance
from astropy import constants
from dust_emissivity import dust

# just hacked together quickly:
# BGPS center is 271.1 GHz
# Assume ALMA center freq is 100 GHz
# use beta = 1.XX for alpha=3.XX
"""
  SpwID  Name                           #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs
  0      ALMA_RB_03#BB_1#SW-01#FULL_RES   7680   TOPO   92220.885      -244.141   1875000.0  91283.5068        1  XX
  1      ALMA_RB_03#BB_2#SW-01#FULL_RES   7680   TOPO   90417.113      -244.141   1875000.0  89479.7356        2  XX
  2      ALMA_RB_03#BB_3#SW-01#FULL_RES   7680   TOPO  100429.129       244.141   1875000.0 101366.5068        3  XX
  3      ALMA_RB_03#BB_4#SW-01#FULL_RES   7680   TOPO  102289.244       244.141   1875000.0 103226.6223        4  XX

"""

band_edges = u.Quantity([
    u.Quantity((092220.885*u.MHz - 244.141*7680*u.kHz,  92220.885*u.MHz)),
    u.Quantity((090417.113*u.MHz - 244.141*7680*u.kHz,  90417.113*u.MHz)),
    u.Quantity((100429.129*u.MHz + 244.141*7680*u.kHz, 100429.129*u.MHz)),
    u.Quantity((102289.244*u.MHz + 244.141*7680*u.kHz, 102289.244*u.MHz)),
])

# because all bands are the same width
centerfreq = band_edges.mean()

beta = 1.75

def mass_conversion_factor(TK=20, d=distance.to(u.kpc), beta=beta):
    return dust.massofsnu(nu=centerfreq, snu=1*u.Jy, distance=d,
                          temperature=u.Quantity(TK, u.K), beta=beta)

def col_conversion_factor(snu, beamomega, TK=20, beta=beta):
    return dust.colofsnu(nu=centerfreq, snu=snu, beamomega=beamomega,
                         temperature=u.Quantity(TK, u.K), beta=beta)

def Jnu(T, nu):
    return (2*constants.h*nu**3 / constants.c**2 *
            np.exp(constants.h*nu/(constants.k_B*T)-1)**-1)
