import numpy as np
import pylab as pl
from astropy import constants
from astropy import units as u
from sedfitter.sed import SEDCube
from astropy.table import Table

seds_S = SEDCube.read('spubsmi/flux.fits')
pars_S = Table.read('spubsmi/parameters.fits')
seds_H = SEDCube.read('spubhmi/flux.fits')
pars_H = Table.read('spubhmi/parameters.fits')

luminosities_H = (constants.sigma_sb * (pars_H['star.temperature']*u.K)**4 *
                  (4*np.pi*(pars_H['star.radius']*u.R_sun)**2)).to(u.L_sun)

flux_sgrb2_H = seds_H.val[:,4,6]*(1/8.4)**2

luminosities_S = (constants.sigma_sb * (pars_S['star.temperature']*u.K)**4 *
                  (4*np.pi*(pars_S['star.radius']*u.R_sun)**2)).to(u.L_sun)

flux_sgrb2_S = seds_S.val[:,4,6]*(1/8.4)**2

pl.figure(1).clf()
pl.plot(luminosities_H, flux_sgrb2_H, '.', label='spubhmi')
pl.plot(luminosities_S, flux_sgrb2_S, '.', label='spubsmi')
pl.ylabel("mJy @ Sgr B2")
pl.xlabel("Luminosity (Lsun)")
pl.plot([1,1e7],[0.5,0.5],'k--')
pl.legend(loc='best')
pl.loglog()
pl.axis([1e3,1e6,0.1,10])
pl.savefig("Robitaille_Models_Flux_vs_Luminosity.pdf")

pl.figure(2).clf()
pl.hist(luminosities_H[flux_sgrb2_H>0.5*u.mJy], bins=np.logspace(3,6,31))
pl.hist(luminosities_S[flux_sgrb2_S>0.5*u.mJy], bins=np.logspace(3,6,31))
pl.semilogx()
pl.savefig("Robitaille_Models_flux_histograms.pdf")

pl.hist(luminosities_H[(flux_sgrb2_H>0.5*u.mJy) & (flux_sgrb2_H<1.0*u.mJy)], bins=np.logspace(3,6,31))
pl.hist(luminosities_S[(flux_sgrb2_S>0.5*u.mJy) & (flux_sgrb2_S<1.0*u.mJy)], bins=np.logspace(3,6,31))
