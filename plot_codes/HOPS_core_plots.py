import numpy as np
import paths
from astropy.table import Table, Column, join
from astropy import units as u
from astropy import coordinates
from astroquery.vizier import Vizier
import masscalc
import pylab as pl
from constants import frequency, distance
pl.matplotlib.rc_file('pubfiguresrc')

pl.rcParams['figure.figsize']=(12,8)
pl.rcParams['figure.dpi']=75
pl.rcParams['savefig.dpi']=150
pl.rcParams['font.size'] = 14

core_phot_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD.ipac"),
                           format='ascii.ipac')

highconf = core_phot_tbl['color']=='green'
lowconf = core_phot_tbl['color']=='orange'
hii = core_phot_tbl['SIMBAD_OTYPE'] == 'HII'

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


hs,l,p = ax1.hist([core_phot_tbl['peak'][highconf & ~hii],
                   core_phot_tbl['peak'][lowconf & ~hii],
                   core_phot_tbl['peak'][hii],
                  ],
                  log=False,
                  label=['Sgr B2 optimistic',
                         'Sgr B2 conservative',
                         'Sgr B2 HII'],
                  color=['#d62728','#2ca02c','#17bcef'],
                  edgecolor='none',
                  rwidth=1,
                  stacked=True,
                  histtype='barstacked',
                  bins=np.logspace(-4,0.2,50))
(hh,hl,hhii) = hs
for pc, hatch in zip(p, ['/','-','x']):
    for patch in pc:
        patch.set_hatch(hatch)
        color = patch.get_facecolor()
        patch.set_edgecolor(color)
        color = (*color[:3], 0.5)
        patch.set_facecolor(color)


hs,l,p = ax1.hist([flux_3mm_cmz[class0],
                   flux_3mm_cmz[classI],
                   #flux_3mm_cmz[classII],
                   flux_3mm_cmz[classflat],
                  ], label=['HOPS Class 0', 'HOPS Class I',
                            #'HOPS Class II',
                            'HOPS Flat'],
                  edgecolor='none',
                  rwidth=1,
                  stacked=True,
                  bins=np.logspace(-7,-3.5,50), histtype='barstacked')
for pc, hatch in zip(p, ['\\','|','+']):
    for patch in pc:
        patch.set_hatch(hatch)
        color = patch.get_facecolor()
        patch.set_edgecolor(color)
        color = (*color[:3], 0.5)
        patch.set_facecolor(color)


ylim = ax1.get_ylim()
#mx = np.logspace(-4, 0.2, 50)
#ax1.plot(mx, (mx/2e-3)**-2.35 * 30, 'k--')

ax1.set_xscale('log')
#ax1.set_xlim(l[:-1][hh>0].min()/1.1, l[1:][hh>0].max()*1.1)
ax1.set_xlim(8e-6, 10)
#ax1.set_ylim(0.6, 15)
ax1.set_ylim(*ylim)
pl.setp(ax1.get_xticklabels(), rotation='horizontal', fontsize=20)
pl.setp(ax1.get_yticklabels(), rotation='vertical', fontsize=20)
ax1.set_xlabel("$S_{3 mm}$ (Jy)", fontsize=22)
ax1.set_ylabel("$N(cores)$", fontsize=22)
pl.legend(loc='best', fontsize=20)

fig1.savefig(paths.fpath('core_peak_intensity_histogram_withHOPS.png'),
             bbox_inches='tight')

# Signal-to-noise plot
fig2 = pl.figure(2)
fig2.clf()
ax2 = fig2.gca()
sn=core_phot_tbl['peak']/core_phot_tbl['bgmad']
hs,l,p = ax2.hist([sn[highconf & ~hii],
                   sn[lowconf & ~hii],
                   sn[hii],
                  ],
                  log=False,
                  label=['Sgr B2 optimistic',
                         'Sgr B2 conservative',
                         'Sgr B2 HII'],
                  color=['#d62728','#2ca02c','#17bcef'],
                  bins=np.logspace(0,2,50), histtype='barstacked')
for pc, hatch in zip(p, ['/','-','x']):
    for patch in pc:
        patch.set_hatch(hatch)
        color = patch.get_facecolor()
        patch.set_edgecolor(color)
        color = (*color[:3], 0.5)
        patch.set_facecolor(color)
ax2.set_xscale('log')
ax2.set_xlabel("S/N", fontsize=12)
ax2.set_ylabel("$N(cores)$", fontsize=12)
pl.legend(loc='best', fontsize=10)
