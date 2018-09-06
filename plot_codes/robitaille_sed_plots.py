import os
import numpy as np
import pylab as pl
import glob
from astropy import constants
from astropy import units as u
from sedfitter.sed import SEDCube
from astropy.table import Table

d_sgrb2 = 8.4 # kpc

seds_S = SEDCube.read('spubsmi/flux.fits')
pars_S = Table.read('spubsmi/parameters.fits')
seds_H = SEDCube.read('spubhmi/flux.fits')
pars_H = Table.read('spubhmi/parameters.fits')

luminosities_H = (constants.sigma_sb * (pars_H['star.temperature']*u.K)**4 *
                  (4*np.pi*(pars_H['star.radius']*u.R_sun)**2)).to(u.L_sun)

flux_sgrb2_H = seds_H.val[:,4,6]*(1/d_sgrb2)**2

luminosities_S = (constants.sigma_sb * (pars_S['star.temperature']*u.K)**4 *
                  (4*np.pi*(pars_S['star.radius']*u.R_sun)**2)).to(u.L_sun)

flux_sgrb2_S = seds_S.val[:,4,6]*(1/d_sgrb2)**2

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


threshold = 0.5 * u.mJy

for seds_fn, pars_fn in zip(glob.glob("*/flux.fits"), glob.glob("*/parameters.fits")):
    seds = SEDCube.read(seds_fn)
    pars = Table.read(pars_fn)

    luminosities = (constants.sigma_sb * (pars['star.temperature']*u.K)**4 *
                    (4*np.pi*(pars['star.radius']*u.R_sun)**2)).to(u.L_sun)

    flux_sgrb2 = seds.val[:,4,6]*(1/d_sgrb2)**2

    match = flux_sgrb2 > threshold
    l104 = (luminosities > 1e4*u.L_sun) & (np.isfinite(flux_sgrb2))
    n_match = (l104 & match).sum()

    l105 = (luminosities > 1e5*u.L_sun) & (np.isfinite(flux_sgrb2))
    n_match_105 = (l105 & match).sum()

    print(os.path.split(seds_fn)[0])
    print("Recovered {0} of {1} luminosities >10^4 Lsun ({2:0.1f}%) at >{3}"
          .format(n_match, l104.sum(), n_match/l104.sum() * 100, threshold))

    print("Recovered {0} of {1} luminosities >10^5 Lsun ({2:0.1f}%) at >{3}"
          .format(n_match_105, l105.sum(), n_match_105/l105.sum() * 100, threshold))
