import numpy as np
import paths
from astropy.table import Table, Column
from astropy import units as u
from astropy import coordinates
import masscalc
import powerlaw
import plfit
import pylab as pl
pl.matplotlib.rc_file('pubfiguresrc')

pl.rcParams['figure.figsize']=(12,8)
pl.rcParams['figure.dpi']=75
pl.rcParams['savefig.dpi']=150
pl.rcParams['font.size'] = 14

core_phot_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD.ipac"), format='ascii.ipac')

highconf = core_phot_tbl['color']=='green'
lowconf = core_phot_tbl['color']=='orange'
hii = core_phot_tbl['SIMBAD_OTYPE'] == 'HII'

peak_fluxdens = core_phot_tbl['peak']
peak_brightness = core_phot_tbl['peak_K']

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

hs,l,p = ax1.hist([core_phot_tbl['peak'][highconf & ~hii],
                   core_phot_tbl['peak'][lowconf & ~hii],
                   core_phot_tbl['peak'][hii],
                  ], log=False,
                  label=['Sgr B2 conservative','Sgr B2 aggressive',
                         'Sgr B2 HII'],
                  color=['#d62728','#2ca02c','#17bcef'],
                  bins=np.logspace(-4,0.2,50),
                  edgecolor='none',
                  rwidth=1,
                  stacked=True,
                  histtype='barstacked')
(hh,hl,hhii) = hs

ax1.set_xscale('log')
ax1.set_xlim(l[:-1][hh>0].min()/1.1, l[1:][hh>0].max()*1.1)
#ax1.set_ylim(0.6, 15)
pl.setp(ax1.get_xticklabels(), rotation='horizontal', fontsize=20)
pl.setp(ax1.get_yticklabels(), rotation='vertical', fontsize=20)
ax1.set_xlabel("$S_{3 mm}$ (Jy)", fontsize=22)
ax1.set_ylabel("$N(cores)$", fontsize=22)
pl.legend(loc='best', fontsize=20)
pl.savefig(paths.fpath("core_peak_fluxdensity_coloredbyclass.png"), bbox_inches='tight')


#fig1.clf()
ax1 = fig1.gca()
plf = plfit.plfit(peak_fluxdens[~hii])
plf.plfit(discrete=False, verbose=True)
plf.plotpdf(dohist=False, fill=True, histcolor='none', plcolor='navy')

pl.setp(ax1.get_xticklabels(), rotation='horizontal')#, fontsize=10)
pl.setp(ax1.get_yticklabels(), rotation='vertical')#, fontsize=10)
ax1.set_xlabel("$S_{3 mm}$ (Jy)")#, fontsize=12)
ax1.set_ylabel("$N(cores)$")#, fontsize=12)
ax1.set_ylim(0.5,30)
ax1.set_yscale('linear')
ax1.set_xlim(0.0003, 2)

fig1.savefig(paths.fpath('core_peak_fluxdensity_powerlawfit.png'), bbox_inches='tight')
p,ksv = plf.test_pl()
print("All Data Consistent with power-law? p={0}".format(p))



plfhi = plfit.plfit(peak_fluxdens[highconf & ~hii])
plfhi.plfit(discrete=False, verbose=True)
plfhi.plotpdf(dohist=False, fill=True, histcolor='none', plcolor='k')
p,ksv = plfhi.test_pl()
print("High-confidence Consistent with power-law? p={0}".format(p))
ax1.set_xlim(0.0003, 2)

fig1.savefig(paths.fpath('core_peak_fluxdensity_powerlawfit_justhi.png'), bbox_inches='tight')


pl_alstott = powerlaw.Fit(peak_fluxdens[~hii])
print("Powerlaw from Alstott: xmin={0}, alpha={1}".format(pl_alstott.xmin, pl_alstott.alpha))
#pl.figure(2).clf()
#pl_alstott.plot_pdf()
