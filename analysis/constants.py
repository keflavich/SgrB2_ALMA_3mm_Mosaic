from astropy import units as u
import numpy as np

distance = 8.4*u.kpc # CMZ distance
#distance = 7.8*u.kpc # Sgr B2 distance ?
frequency = 95*u.GHz # TODO: make this more accurate
# added long after publication (Apr 24, 2019):
frequency = np.mean((89.48, 91.28, 101.37, 103.23))*u.GHz
# 96.34 GHz

# from stellar_mass_estimates
mass_represented_by_a_source = 21.8*u.M_sun / 0.2277
