import numpy as np
import paths
from astropy.table import Table, Column, join
from astropy import units as u
from astropy import coordinates
from astroquery.vizier import Vizier
import masscalc
import pylab as pl
from constants import frequency, distance
pl.matplotlib.rc_file('/Users/adam/.matplotlib/pubfiguresrc')

core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')

highconf = core_phot_tbl['color']=='green'
lowconf = core_phot_tbl['color']=='orange'

Vizier.ROW_LIMIT = 1000000
HOPS_classes = Vizier.get_catalogs('2016ApJS..224....5F')[0]
HOPS_fluxes = Vizier.get_catalogs('2016ApJS..224....5F')[1]
HOPS_tbl = join(HOPS_classes, HOPS_fluxes, keys='HOPS')

f870 = HOPS_tbl['F870']
d_hops = 415*u.pc
flux_3mm_cmz = (u.Quantity(f870,u.Jy) *
                (870*u.um/(frequency.to(u.um,u.spectral())))**3.5 *
                (d_hops/distance)**2).to(u.Jy)
OK = np.isfinite(flux_3mm_cmz)

class0 = OK & (HOPS_tbl['Class'] == b'0')
classI = OK & (HOPS_tbl['Class'] == b'I')
classII = OK & (HOPS_tbl['Class'] == b'II')
classflat = OK & (HOPS_tbl['Class'] == b'flat')

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

ax1.hist([flux_3mm_cmz[class0],
          flux_3mm_cmz[classI],
          flux_3mm_cmz[classII],
          flux_3mm_cmz[classflat],
         ], label=['HOPS Class 0', 'HOPS Class I', 'HOPS Class II', 'HOPS Flat'],
         bins=np.logspace(-7,-3.5,50), histtype='barstacked')

(hh,hl),l,p = ax1.hist([core_phot_tbl['peak'][highconf],
                        core_phot_tbl['peak'][lowconf]], log=False,
                       label=['Sgr B2 bright','Sgr B2 faint'],
                       bins=np.logspace(-4,0.2,50), histtype='barstacked')
ax1.set_xscale('log')
#ax1.set_xlim(l[:-1][hh>0].min()/1.1, l[1:][hh>0].max()*1.1)
ax1.set_xlim(8e-6, 10)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("$S_{3 mm}$ (Jy)")
ax1.set_ylabel("$N(cores)$")
pl.legend(loc='best', fontsize=14)

fig1.savefig(paths.fpath('core_peak_intensity_histogram_withHOPS.png'))
