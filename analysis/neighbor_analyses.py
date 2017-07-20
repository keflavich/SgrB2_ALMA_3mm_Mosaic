import numpy as np
import pylab as pl

from astropy import coordinates, units as u
from astropy.table import Table

from constants import distance
import paths

cont_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD_andclusters.ipac"), format='ascii.ipac')
sgrb2_coords = coordinates.SkyCoord(cont_tbl['RA'], cont_tbl['Dec'],
                                    unit=(u.deg, u.deg), frame='fk5',)

nn = {}
for ii in range(2,15):
    nn[ii] = coordinates.match_coordinates_sky(sgrb2_coords, sgrb2_coords,
                                               nthneighbor=ii)

pl.clf()
for ii in nn:
    pl.plot(sorted((nn[ii][1]*distance).to(u.pc,
                                           u.dimensionless_angles()).value),
            np.linspace(0, 1, len(cont_tbl)),
            label=ii)
pl.xlabel("Separation (pc)")
