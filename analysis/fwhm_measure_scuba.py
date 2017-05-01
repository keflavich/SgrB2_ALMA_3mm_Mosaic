
import numpy as np
import gaussfitter
from astropy.io import fits
from astropy import units as u
from astropy import wcs

cutout_size=5

im = fits.getdata('gc450.fits.gz')[462-cutout_size:462+cutout_size,1392-cutout_size:1392+cutout_size]
pars = gaussfitter.gaussfit(im, params=[1, 20, cutout_size, cutout_size, 2, 2, 0])
w = wcs.WCS(fits.getheader('gc450.fits.gz'))

fwhm = (wcs.utils.proj_plane_pixel_scales(w)*u.deg).to(u.arcsec) * pars[4:6] * np.sqrt(8*np.log(2))

print("FWHM of Sgr B2 N = {0}".format(fwhm))

im = fits.getdata('gc450.fits.gz')[454-cutout_size:454+cutout_size, 1405-cutout_size:1405+cutout_size]
pars = gaussfitter.gaussfit(im, params=[1, 20, cutout_size, cutout_size, 2, 2, 0])
w = wcs.WCS(fits.getheader('gc450.fits.gz'))

fwhm = (wcs.utils.proj_plane_pixel_scales(w)*u.deg).to(u.arcsec) * pars[4:6] * np.sqrt(8*np.log(2))

print("FWHM of Sgr B2 M = {0}".format(fwhm))

im = fits.getdata('gc450.fits.gz')[456-cutout_size:456+cutout_size, 1417-cutout_size:1417+cutout_size]
pars = gaussfitter.gaussfit(im, params=[2, 10, cutout_size, cutout_size, 2, 2, 0])
w = wcs.WCS(fits.getheader('gc450.fits.gz'))

fwhm = (wcs.utils.proj_plane_pixel_scales(w)*u.deg).to(u.arcsec) * pars[4:6] * np.sqrt(8*np.log(2))

print("FWHM of Sgr B2 S = {0}".format(fwhm))


for cx,cy in [(3239,536),(1826,255),(2087,298),(2112,285),(1201,716)]:
    im = fits.getdata('gc450.fits.gz')[cy-cutout_size:cy+cutout_size, cx-cutout_size:cx+cutout_size]
    pars = gaussfitter.gaussfit(im, params=[2, 10, cutout_size, cutout_size, 2, 2, 0])
    w = wcs.WCS(fits.getheader('gc450.fits.gz'))

    fwhm = (wcs.utils.proj_plane_pixel_scales(w)*u.deg).to(u.arcsec) * pars[4:6] * np.sqrt(8*np.log(2))

    print("FWHM of some isolated source = {0}".format(fwhm))
