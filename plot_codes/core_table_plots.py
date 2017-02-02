import numpy as np
import paths
from astropy.table import Table, Column
from astropy import units as u
from astropy import coordinates
import masscalc
import pylab as pl
pl.matplotlib.rc_file('/Users/adam/.matplotlib/pubfiguresrc')

core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')

peak_brightness = core_phot_tbl['peak_K']

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

h,l,p = ax1.hist(core_phot_tbl['peak'], log=False, bins=np.logspace(-4,0.2,50))
ax1.set_xscale('log')
ax1.set_xlim(l[:-1][h>0].min()/1.1, l[1:][h>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("$S_{3 mm}$ (Jy)")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_peak_intensity_histogram.png'))

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

h,l,p = ax1.hist(peak_brightness, log=False, bins=np.logspace(-0.6, 3.0,30))
ax1.set_xscale('log')
ax1.set_xlim(l[:-1][h>0].min()/1.1, l[1:][h>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("$T_{B,3 mm}$ [K]")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_peak_brightness_histogram.png'))


fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

h,l,p = ax1.hist((core_phot_tbl['peak']/core_phot_tbl['bgmad'])[core_phot_tbl['color']=='green'],
                 log=False, bins=np.logspace(0,2.5))
ho,lo,po = ax1.hist((core_phot_tbl['peak']/core_phot_tbl['bgmad'])[core_phot_tbl['color']=='orange'],
                    log=False, bins=np.logspace(0,2.5))
ax1.set_xscale('log')
#ax1.set_xlim(l[:-1][h>0].min()/1.1, l[1:][h>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("S/N")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_SN_histogram.png'))




fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

shift = np.log10(masscalc.mass_conversion_factor().value)
h,l,p = ax1.hist(core_phot_tbl['peak']*masscalc.mass_conversion_factor(TK=40, beta=1.75),
                 log=False, bins=np.logspace(-4+shift,-1+shift,50))
ax1.set_xscale('log')
ax1.set_xlim(l[:-1][h>0].min()/1.1, l[1:][h>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("$M(T=40$ K, $\\beta=1.75$) [$M_\odot$]")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_mass_histogram_40K.png'))


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


# spectral index derivation and plotting

cont_tbl = core_phot_tbl
fluxratio = cont_tbl['peak_100GHz'] / cont_tbl['peak_90GHz']
spindx = np.log(fluxratio) / np.log(100./90.)
sigma_R = ((cont_tbl['bgmad_100GHz'] / cont_tbl['peak_100GHz'])**2 +
           (cont_tbl['bgmad_90GHz'] / cont_tbl['peak_90GHz'])**2
          )**0.5 * fluxratio
spindx_err = sigma_R / fluxratio / np.log(100./90.)

significant_mask = np.abs(spindx) - 3*spindx_err > 0

#pl.errorbar(cont_tbl['peak_90GHz'], cont_tbl['peak_100GHz'],
#            xerr=cont_tbl['bgmad_90GHz'], yerr=cont_tbl['bgmad_100GHz'],
#            marker='.', linestyle='')
pl.semilogx()
pl.errorbar(cont_tbl['peak_90GHz'][significant_mask], spindx[significant_mask],
            xerr=cont_tbl['bgmad_90GHz'][significant_mask], yerr=spindx_err[significant_mask],
            marker='s', linestyle='', color='k', zorder=10)
pl.errorbar(cont_tbl['peak_90GHz'][~significant_mask], spindx[~significant_mask],
            xerr=cont_tbl['bgmad_90GHz'][~significant_mask], yerr=spindx_err[~significant_mask],
            marker='.', linestyle='', color='b')
pl.ylim(-4,4)
