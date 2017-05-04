import numpy as np
import paths
from astropy.table import Table, Column
from astropy import units as u
from astropy import coordinates
import masscalc
import pylab as pl
pl.matplotlib.rc_file('pubfiguresrc')

cont_tbl = core_phot_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD.ipac"), format='ascii.ipac')

highconf = core_phot_tbl['color']=='green'
lowconf = core_phot_tbl['color']=='orange'
hii = core_phot_tbl['SIMBAD_OTYPE'] == 'HII'





peak_brightness = core_phot_tbl['peak_K']

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

(hh,hl),l,p = ax1.hist([core_phot_tbl['peak'][highconf],
                        core_phot_tbl['peak'][lowconf]], log=False,
                       bins=np.logspace(-4,0.2,50), histtype='barstacked')
ax1.set_xscale('log')
ax1.set_xlim(l[:-1][hh>0].min()/1.1, l[1:][hh>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("$S_{3 mm}$ (Jy)")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_peak_intensity_histogram.png'), bbox_inches='tight')

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

(hh,hl),l,p = ax1.hist([peak_brightness[highconf], peak_brightness[lowconf]],
                       log=False, bins=np.logspace(-0.6, 3.0,30),
                       histtype='barstacked')
ax1.set_xscale('log')
ax1.set_xlim(l[:-1][hh>0].min()/1.1, l[1:][hh>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("$T_{B,3 mm}$ [K]")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_peak_brightness_histogram.png'), bbox_inches='tight')


fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

h,l,p = ax1.hist([(core_phot_tbl['peak']/core_phot_tbl['bgmad'])[highconf],
                  (core_phot_tbl['peak']/core_phot_tbl['bgmad'])[lowconf]],
                 log=False, bins=np.logspace(0,2.5), histtype='barstacked')
#ho,lo,po = ax1.hist((core_phot_tbl['peak']/core_phot_tbl['bgmad'])[core_phot_tbl['color']=='orange'],
#                    log=False, bins=np.logspace(0,2.5))
ax1.set_xscale('log')
#ax1.set_xlim(l[:-1][h>0].min()/1.1, l[1:][h>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
ax1.set_xlabel("S/N")
ax1.set_ylabel("$N(cores)$")

fig1.savefig(paths.fpath('core_SN_histogram.png'), bbox_inches='tight')




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

fig1.savefig(paths.fpath('core_mass_histogram_40K.png'), bbox_inches='tight')


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
fig2.savefig(paths.fpath('core_spatial_distribution.png'), bbox_inches='tight')


# spectral index derivation and plotting

cont_tbl = core_phot_tbl
fluxratio = cont_tbl['peak_100GHz'] / cont_tbl['peak_90GHz']
spindx = np.log(fluxratio) / np.log(100./90.)
sigma_R = ((cont_tbl['bgmad_100GHz'] / cont_tbl['peak_100GHz'])**2 +
           (cont_tbl['bgmad_90GHz'] / cont_tbl['peak_90GHz'])**2
          )**0.5 * fluxratio
spindx_err = sigma_R / fluxratio / np.log(100./90.)

significant_mask = (np.abs(spindx) - 3*spindx_err > 0) | (spindx_err < 0.1)

#pl.errorbar(cont_tbl['peak_90GHz'], cont_tbl['peak_100GHz'],
#            xerr=cont_tbl['bgmad_90GHz'], yerr=cont_tbl['bgmad_100GHz'],
#            marker='.', linestyle='')
fig3 = pl.figure(3)
fig3.clf()
ax3 = fig3.gca()
ax3.semilogx()
ax3.errorbar(cont_tbl['peak_90GHz'][significant_mask],
             spindx[significant_mask],
             xerr=cont_tbl['bgmad_90GHz'][significant_mask],
             yerr=spindx_err[significant_mask], marker='s', linestyle='',
             color='k', zorder=10)
ax3.errorbar(cont_tbl['peak_90GHz'][~significant_mask],
             spindx[~significant_mask],
             xerr=cont_tbl['bgmad_90GHz'][~significant_mask],
             yerr=spindx_err[~significant_mask], marker='.', linestyle='',
             alpha=0.25,
             color='b')
ax3.set_ylim(-4,4)
ax3.set_xlabel("90 GHz peak $S_\\nu$")
ax3.set_ylabel("90-100 GHz Spectral Index")



alphaok_mask = (np.abs(cont_tbl['alpha']) > cont_tbl['alphaerror']*5) | (cont_tbl['alphaerror'] < 0.1)

fig4 = pl.figure(4)
fig4.clf()
ax4 = fig4.gca()
ax4.set_xscale('log')
ax4.errorbar(cont_tbl['peak'][alphaok_mask], cont_tbl['alpha'][alphaok_mask],
             xerr=cont_tbl['bgmad'][alphaok_mask], yerr=cont_tbl['alphaerror'][alphaok_mask],
             linewidth=0.5, alpha=0.5,
             marker='s', linestyle='', color='k', zorder=10)

fig5 = pl.figure(5)
fig5.clf()
ax5 = fig5.gca()
hs,l,p = ax5.hist([core_phot_tbl['alpha'][highconf & ~hii & alphaok_mask],
                   core_phot_tbl['alpha'][lowconf & ~hii & alphaok_mask],
                   core_phot_tbl['alpha'][hii & alphaok_mask],
                  ], log=False,
                  label=['Sgr B2 conservative','Sgr B2 aggressive',
                         'Sgr B2 HII'],
                  color=['#d62728','#2ca02c','#17bcef'],
                  bins=np.linspace(-3.0, 5, 20),
                  edgecolor='w',
                  rwidth=1,
                  stacked=True,
                  histtype='barstacked')
(hh,hl,hhii) = hs

#ax5.set_xscale('log')
ax5.set_xlim(l[:-1][hh>0].min()*1.1, l[1:][hh>0].max()*1.1)
ax5.set_xlim(-3.0, 5)
ax5.set_ylim(0, 11)
pl.setp(ax5.get_xticklabels(), rotation='horizontal', fontsize=20)
pl.setp(ax5.get_yticklabels(), rotation='vertical', fontsize=20)
ax5.set_xlabel("Spectral Index $\\alpha$", fontsize=22)
ax5.set_ylabel("$N(cores)$", fontsize=22)
pl.legend(loc='best', fontsize=20)
pl.savefig(paths.fpath("core_alpha_coloredbyclass.png"), bbox_inches='tight')


fig6 = pl.figure(6)
fig6.clf()
ax6 = fig6.gca()
hs,l,p = ax6.hist([core_phot_tbl['peak'][highconf & ~hii & alphaok_mask],
                   core_phot_tbl['peak'][lowconf & ~hii & alphaok_mask],
                   core_phot_tbl['peak'][hii & alphaok_mask],
                  ], log=False,
                  label=['Sgr B2 conservative','Sgr B2 aggressive',
                         'Sgr B2 HII'],
                  color=['#d62728','#2ca02c','#17bcef'],
                  bins=np.logspace(-3.0, 0.5, 20),
                  edgecolor='w',
                  rwidth=1,
                  stacked=True,
                  histtype='barstacked')
(hh,hl,hhii) = hs

ax6.set_xscale('log')
ax6.set_xlim(l[:-1][hh>0].min()/1.1, l[1:][hh>0].max()*1.1)
ax6.set_xlim(1e-3, 2)
ax6.set_ylim(0, 15)
pl.setp(ax6.get_xticklabels(), rotation='horizontal', fontsize=20)
pl.setp(ax6.get_yticklabels(), rotation='vertical', fontsize=20)
ax6.set_xlabel("$S_{3 mm}$ (Jy)", fontsize=22)
ax6.set_ylabel("$N(cores)$", fontsize=22)
pl.legend(loc='best', fontsize=20)
pl.savefig(paths.fpath("core_peak_withalphameasurements_coloredbyclass.png"), bbox_inches='tight')



# Separate figure just shows that the calculated-by-hand version doesn't match
# the measured version (or at least, didn't)
pl.figure(7).clf()
pl.plot(cont_tbl['alpha'],
        spindx,
        'r.', markersize=2)
pl.plot(cont_tbl['alpha'][alphaok_mask & significant_mask],
        spindx[alphaok_mask & significant_mask],
        'ko', alpha=0.5)
pl.plot([-2,4],[-2,4], 'k--')
pl.xlabel("CASA Alpha")
pl.ylabel("90-100 GHz spectral index")
pl.axis([-2,4,-2,4])

pl.draw()
pl.show()
