import numpy as np
import paths
from astropy.table import Table, Column
from astropy import units as u
from astropy import coordinates
import masscalc
import pylab as pl
pl.matplotlib.rc_file('/Users/adam/.matplotlib/pubfiguresrc')

core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

h,l,p = ax1.hist(core_phot_tbl['peak'], log=False, bins=np.logspace(-4,-1,50))
ax1.set_xscale('log')
ax1.set_xlim(l[h>0].min()/1.1, l[h>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("$S_{3 mm}$ (Jy)")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_flux_histogram.png'))

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

shift = np.log10(masscalc.mass_conversion_factor())
h,l,p = ax1.hist(core_phot_tbl['peak']*masscalc.mass_conversion_factor(),
                 log=False, bins=np.logspace(-4+shift,-1+shift,50))
ax1.set_xscale('log')
ax1.set_xlim(l[h>0].min()/1.1, l[h>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("M(T=20 K, $\\beta=1.5$)")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_mass_histogram.png'))


fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.gca()

coords = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                              frame='fk5').galactic
ax2.plot(coords.l, coords.b, '.')
ax2.set_xlim(*ax2.get_xlim()[::-1])
ax2.set_ylabel('Galactic Latitude')
ax2.set_xlabel('Galactic Longitude')
ax2.set_aspect(1)
fig2.savefig(paths.fpath('core_spatial_distribution.png'))
