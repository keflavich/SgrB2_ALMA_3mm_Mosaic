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

core_phot_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD_andclusters.ipac"), format='ascii.ipac')

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
                  label=['Sgr B2 high-confidence',
                         'Sgr B2 low-confidence',
                         'Sgr B2 HII'],
                  color=['#d62728','orange','#17bcef'],
                  bins=np.logspace(-4,0.2,50),
                  edgecolor='none',
                  rwidth=1,
                  stacked=True,
                  histtype='barstacked')
(hh,hl,hhii) = hs
for pc, hatch in zip(p, ['/','\\','+']):
    for patch in pc:
        patch.set_hatch(hatch)
        patch.set_edgecolor('w')
        #color = patch.get_facecolor()
        #patch.set_edgecolor(color)
        #color = (*color[:3], 0.5)
        #patch.set_facecolor(color)

ax1.set_xscale('log')
ax1.set_xlim(l[:-1][hh>0].min()/1.1, l[1:][(hh>0)|(hhii>0)].max()*1.1)
#ax1.set_ylim(0.6, 15)
pl.setp(ax1.get_xticklabels(), rotation='horizontal', fontsize=20)
pl.setp(ax1.get_yticklabels(), rotation='vertical', fontsize=20)
ax1.set_xlabel("$S_{3 mm}$ (Jy)", fontsize=22)
ax1.set_ylabel("$N(cores)$", fontsize=22)
pl.legend(loc='best', fontsize=20)
pl.savefig(paths.fpath("core_peak_fluxdensity_coloredbyclass.pdf"), bbox_inches='tight')

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

clusterM = core_phot_tbl['Cluster'] == 'M'
clusterN = core_phot_tbl['Cluster'] == 'N'
clusterNE = core_phot_tbl['Cluster'] == 'NE'
clusterS = core_phot_tbl['Cluster'] == 'S'

#hs,l,p = ax1.hist([core_phot_tbl['peak'][clusterM],
#                   core_phot_tbl['peak'][clusterN],
#                   core_phot_tbl['peak'][clusterNE],
#                   core_phot_tbl['peak'][~(clusterNE|clusterM|clusterN)],
#                  ], log=False,
#                  label=['Sgr B2 M', 'Sgr B2 N',
#                         'Sgr B2 NE', 'No cluster'],
#                  color=['#0066DD', '#00DD66', '#00AAAA', '#DD5522', ],
#                  bins=np.logspace(-4,0.2,30),
#                  edgecolor='none',
#                  rwidth=1,
#                  stacked=True,
#                  histtype='barstacked')
#(hM,hN,hNE,hOther) = hs
#for pc, hatch in zip(p, ['/','\\','+']):
#    for patch in pc:
#        patch.set_hatch(hatch)
#        patch.set_edgecolor('w')
        #color = patch.get_facecolor()
        #patch.set_edgecolor(color)
        #color = (*color[:3], 0.5)
        #patch.set_facecolor(color)

hM,l,pM = ax1.hist(core_phot_tbl['peak'][clusterM], label='Sgr B2 M',
                   histtype='step',
                   linewidth=5,
                   zorder=10,
                   alpha=0.8,
                   edgecolor='#0066DD',
                   facecolor='none',
                   #facecolor=(0, 0x66/256., 0xDD/256., 0.2),
                   bins=np.logspace(-4,0.2,20),)
hN,l,pN = ax1.hist(core_phot_tbl['peak'][clusterN], label='Sgr B2 N',
                   histtype='stepfilled',
                   linewidth=2,
                   edgecolor='#000000',
                   facecolor=(0x00/256.,0x00/256,0x00/256,0.2),
                   bins=np.logspace(-4,0.2,20),)
hNE,l,pNE = ax1.hist(core_phot_tbl['peak'][clusterNE], label='Sgr B2 NE',
                     histtype='stepfilled',
                     linewidth=2,
                     edgecolor='#00AAAA',
                     facecolor=(0, 0xAA/256, 0xAA/256, 0.2),
                     bins=np.logspace(-4,0.2,20),)
hS,l,pS = ax1.hist(core_phot_tbl['peak'][clusterS], label='Sgr B2 S',
                     histtype='stepfilled',
                     linewidth=2,
                     edgecolor='#AAAA00',
                     facecolor=(0xAA/256, 0xAA/256, 0, 0.2),
                     bins=np.logspace(-4,0.2,20),)
hOther,l,pOther = ax1.hist(core_phot_tbl['peak'][~(clusterNE|clusterM|clusterN|clusterS)],
                           label='No cluster', edgecolor='#DD5522',
                           facecolor=(0xDD/256., 0x55/256., 0x22/256., 0.2),
                           histtype='stepfilled',
                           zorder=-5,
                           bins=np.logspace(-4,0.2,50),)
hOther,l,p = ax1.hist(core_phot_tbl['peak'][~(clusterNE|clusterM|clusterN|clusterS)],
                      edgecolor='#DD5522',
                      histtype='step',
                      label=None,
                      zorder=5,
                      bins=np.logspace(-4,0.2,50),)

for pc, hatch in zip([pM, pN, pNE, pOther], ['','\\','+']):
    for patch in pc:
        patch.set_hatch(hatch)
        #patch.set_edgecolor('w')


ax1.set_xscale('log')
ax1.set_xlim(3e-4,10**0.22)
#ax1.set_xlim(l[:-1][hOther>0].min()/1.1, l[1:][(hM>0)|(hN>0)].max()*1.1)
#ax1.set_ylim(0.6, 15)
pl.setp(ax1.get_xticklabels(), rotation='horizontal', fontsize=20)
pl.setp(ax1.get_yticklabels(), rotation='vertical', fontsize=20)
ax1.set_xlabel("$S_{3 mm}$ (Jy)", fontsize=22)
ax1.set_ylabel("$N(cores)$", fontsize=22)
pl.legend(loc='best', fontsize=20)

pl.savefig(paths.fpath("core_peak_fluxdensity_coloredbycluster.pdf"), bbox_inches='tight')



#fig1.clf()
ax1 = fig1.gca()
plf = plfit.plfit(peak_fluxdens[~hii])
plf.plfit(discrete=False, verbose=True)
plf.plotpdf(dohist=False, histcolor='none', plcolor='navy')

pl.setp(ax1.get_xticklabels(), rotation='horizontal')#, fontsize=10)
pl.setp(ax1.get_yticklabels(), rotation='vertical')#, fontsize=10)
ax1.set_xlabel("$S_{3 mm}$ (Jy)")#, fontsize=12)
ax1.set_ylabel("$N(cores)$")#, fontsize=12)
ax1.set_ylim(0.5,30)
ax1.set_yscale('linear')
ax1.set_xlim(0.0003, 2)

fig1.savefig(paths.fpath('core_peak_fluxdensity_powerlawfit.pdf'), bbox_inches='tight')
p,ksv = plf.test_pl()
print("All Data Consistent with power-law? p={0}".format(p))



plfhi = plfit.plfit(peak_fluxdens[highconf & ~hii])
plfhi.plfit(discrete=False, verbose=True)
plfhi.plotpdf(fill=True, histcolor='none', plcolor='k')
p,ksv = plfhi.test_pl()
print("High-confidence Consistent with power-law? p={0}".format(p))
ax1.set_xlim(0.0003, 2)

fig1.savefig(paths.fpath('core_peak_fluxdensity_powerlawfit_justhi.pdf'), bbox_inches='tight')


pl_alstott = powerlaw.Fit(peak_fluxdens[~hii])
print("Powerlaw from Alstott: xmin={0}, alpha={1}".format(pl_alstott.xmin, pl_alstott.alpha))
#pl.figure(2).clf()
#pl_alstott.plot_pdf()
