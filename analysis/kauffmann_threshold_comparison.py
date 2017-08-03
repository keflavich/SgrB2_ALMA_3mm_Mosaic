import numpy as np
from constants import distance
import radio_beam
from astropy import units as u
from astropy import wcs
from astropy.io import fits
import higal_sedfitter
import paths
import reproject
import dust_emissivity
from astropy import convolution
import pylab as pl

herschel25_colfile = paths.Fpath('other/Herschel25umcolum_regridded_match_stellar.fits')
herschel36_colfile = paths.Fpath('other/Herschel36umcolum_regridded_match_stellar.fits')
herschel25 = fits.getdata(herschel25_colfile)*u.M_sun/u.pc**2
herschel36 = fits.getdata(herschel36_colfile)*u.M_sun/u.pc**2

sharc_colfile = paths.cpath('column_maps/sharc_col_herscheltem.fits')
scuba_colfile = paths.cpath('column_maps/scuba_col_herscheltem.fits')
sharc_col = fits.getdata(sharc_colfile)*u.cm**-2
scuba_col = fits.getdata(scuba_colfile)*u.cm**-2

sharc_beam = 10*u.arcsec
scuba_beam = 10*u.arcsec
scuba_beam_pc = (scuba_beam*distance).to(u.pc, u.dimensionless_angles())
h25bm_pc = (25*u.arcsec*distance).to(u.pc, u.dimensionless_angles())
h36bm_pc = (36*u.arcsec*distance).to(u.pc, u.dimensionless_angles())
FWHM_SCALE = np.sqrt(8*np.log(2))


scuba_mass = (scuba_col * 2.8*u.Da *
              2*np.pi * (scuba_beam_pc/FWHM_SCALE)**2).to(u.M_sun)
herschel25mass = (herschel25 * 2*np.pi*(h25bm_pc/FWHM_SCALE)**2).to(u.M_sun)
herschel36mass = (herschel36 * 2*np.pi*(h36bm_pc/FWHM_SCALE)**2).to(u.M_sun)


sgrb2contfile = fits.open(paths.Fpath('merge/continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy.image.tt0.pbcor.fits'))
observed_region,_ = reproject.reproject_interp(sgrb2contfile, fits.getheader(scuba_colfile))
mask_scuba = np.isfinite(observed_region)
observed_region,_ = reproject.reproject_interp(sgrb2contfile, fits.getheader(herschel25_colfile))
mask_herschel25 = np.isfinite(observed_region)
observed_region,_ = reproject.reproject_interp(sgrb2contfile, fits.getheader(herschel36_colfile))
mask_herschel36 = np.isfinite(observed_region)

pl.figure(1).clf()
pl.violinplot([herschel25mass[np.isfinite(herschel25mass) & (herschel25mass>0) & mask_herschel25].value,
               herschel36mass[np.isfinite(herschel36mass) & (herschel36mass>0) & mask_herschel36].value,
               scuba_mass[np.isfinite(scuba_mass) & (scuba_mass>0) & mask_scuba].value],
              positions=([25,36,10]*u.arcsec*distance).to(u.pc, u.dimensionless_angles()).value)
pl.plot([0.1,2], np.array([0.1, 2])**1.33 * 870, 'k-')
pl.gca().set_yscale('log')


pl.figure(2).clf()
pl.imshow(scuba_mass.value, interpolation='nearest', origin='lower',
          cmap='gray', norm=pl.matplotlib.colors.LogNorm())
pl.contour(scuba_mass.value, levels=[870*scuba_beam_pc.to(u.pc).value**1.33], colors=['b'])
pl.contour(mask_scuba, levels=[0.5], colors=['r'])
