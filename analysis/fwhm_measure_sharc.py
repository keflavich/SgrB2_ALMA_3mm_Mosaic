import numpy as np
import gaussfitter
from astropy.io import fits
from astropy import units as u
from astropy import wcs

im = fits.getdata('SgrB2_350um_gal.fits')[970:1070, 830:930]
pars = gaussfitter.gaussfit(im, params=[1, 1300, 50, 50, 10, 10, 0])
w = wcs.WCS(fits.getheader('SgrB2_350um_gal.fits'))

fwhm = (wcs.utils.proj_plane_pixel_scales(w)*u.deg).to(u.arcsec) * pars[4:6] * np.sqrt(8*np.log(2))

print("FWHM of Sgr B2 N = {0}".format(fwhm))

im = fits.getdata('SgrB2_350um_gal.fits')[933-50:933+50, 994-50:994+50]
pars = gaussfitter.gaussfit(im, params=[1, 1300, 50, 50, 10, 10, 0])
w = wcs.WCS(fits.getheader('SgrB2_350um_gal.fits'))

fwhm = (wcs.utils.proj_plane_pixel_scales(w)*u.deg).to(u.arcsec) * pars[4:6] * np.sqrt(8*np.log(2))

print("FWHM of Sgr B2 M = {0}".format(fwhm))

im = fits.getdata('SgrB2_350um_gal.fits')[956-50:956+50, 756-50:756+50]
pars = gaussfitter.gaussfit(im, params=[250, 50, 50, 50, 14, 14, 0])
w = wcs.WCS(fits.getheader('SgrB2_350um_gal.fits'))

fwhm = (wcs.utils.proj_plane_pixel_scales(w)*u.deg).to(u.arcsec) * pars[4:6] * np.sqrt(8*np.log(2))

print("FWHM of Sgr B2 NE = {0}".format(fwhm))
