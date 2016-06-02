import numpy as np
from astropy.convolution import convolve_fft
from spectral_cube import SpectralCube
from astropy import units as u
import fft_psd_tools
from astropy.io import fits
from astropy import wcs
from reproject import reproject_interp

cube7m = SpectralCube.read('SgrB2_b3_7m.HC3N.fits.image.pbcor.fits')
cube12m = SpectralCube.read('SgrB2_b3_12m.HC3N.fits.image.pbcor.fits')
vcube7m12m = SpectralCube.read('SgrB2_b3_7M_12M_natural.HC3N.image.pbcor.fits')
tp_fn = 'HC3N_natural_tpcube_k_rg_65kms.fits'

jtok_7m = np.array([bm.jtok(x).value for bm,x in zip(cube7m.beams,
                                                     cube7m.spectral_axis)])
jtok_12m = np.array([bm.jtok(x).value for bm,x in zip(cube12m.beams,
                                                      cube12m.spectral_axis)])
jtok_7m12m = np.array([bm.jtok(x).value for bm,x in zip(vcube7m12m.beams,
                                                        vcube7m12m.with_spectral_unit(u.GHz).spectral_axis)])

vcube7m = cube7m.with_spectral_unit(u.km/u.s, velocity_convention='radio')
vcube12m = cube12m.with_spectral_unit(u.km/u.s, velocity_convention='radio')

bg7m = vcube7m.spectral_slab(-150*u.km/u.s, -5*u.km/u.s).median(axis=0)
bg12m = vcube12m.spectral_slab(-150*u.km/u.s, -5*u.km/u.s).median(axis=0)
bg7m12m = vcube7m12m.spectral_slab(-150*u.km/u.s, -5*u.km/u.s).median(axis=0)

ch_65kms_7m = vcube7m.closest_spectral_channel(65*u.km/u.s)
ch_65kms_12m = vcube12m.closest_spectral_channel(65*u.km/u.s)
ch_65kms_7m12m = vcube7m12m.closest_spectral_channel(65*u.km/u.s)

slc7m = cube7m[ch_65kms_7m] - bg7m
slc12m = cube12m[ch_65kms_12m] - bg12m
slc7m12m = vcube7m12m[ch_65kms_12m] - bg7m12m
slcTP = reproject_interp(fits.open(tp_fn), slc12m.header)[0]

pixscale7m = wcs.utils.proj_plane_pixel_area(cube7m.wcs.celestial)**0.5 * u.deg
pixscale12m = wcs.utils.proj_plane_pixel_area(cube12m.wcs.celestial)**0.5 * u.deg
pixscaleTP = pixscale12m
#old pixscaletp = wcs.utils.proj_plane_pixel_area(wcs.WCS(fits.getheader(tp_fn)))**0.5 * u.deg

convbeam_12mto7m = cube7m.beams[ch_65kms_7m].deconvolve(cube12m.beams[ch_65kms_12m])
kernel = convbeam_12mto7m.as_kernel(pixscale12m)
slc12m_conv_to_7m = convolve_fft(slc12m.value, kernel)

frq7m,pow7m = fft_psd_tools.psds.pspec(fft_psd_tools.psds.PSD2(slc7m.value*jtok_7m[ch_65kms_7m]), view=False)
frq12m,pow12m = fft_psd_tools.psds.pspec(fft_psd_tools.psds.PSD2(slc12m.value*jtok_12m[ch_65kms_12m]), view=False)
frq12m_conv,pow12m_conv = fft_psd_tools.psds.pspec(fft_psd_tools.psds.PSD2(slc12m_conv_to_7m*jtok_12m[ch_65kms_12m]), view=False)
frqTP,powTP = fft_psd_tools.psds.pspec(fft_psd_tools.psds.PSD2(slcTP), view=False)

# some UVcombine tests; useful to compare to the above.  Need psds still
import uvcombine
hdu7m = slc7m.hdu.copy()
hdu7m.data *= jtok_7m[ch_65kms_7m]
hdu7m.header['BUNIT'] = 'K'
hdu12m = slc12m.hdu.copy()
hdu12m.data *= jtok_12m[ch_65kms_12m]
hdu12m.header['BUNIT'] = 'K'
im7mTP = uvcombine.feather_simple(hdu7m, tp_fn)
hdu7mTP = hdu7m.copy()
hdu7mTP.data = im7mTP.real
im12m7mTP = uvcombine.feather_simple(hdu12m, hdu7mTP)
hdu12m7mTP = hdu12m.copy()
hdu12m7mTP.data = im12m7mTP.real



frq12m7mTP,pow12m7mTP = fft_psd_tools.psds.pspec(fft_psd_tools.psds.PSD2(hdu12m7mTP.data), view=False)
frq7mTP,pow7mTP = fft_psd_tools.psds.pspec(fft_psd_tools.psds.PSD2(hdu7mTP.data), view=False)

import pylab as pl

pl.figure(1).clf()
pl.loglog(pixscale7m.to(u.arcsec).value/frq7m, pow7m, label='7m')
pl.loglog(pixscale12m.to(u.arcsec).value/frq12m, pow12m, label='12m')
pl.loglog(pixscaleTP.to(u.arcsec).value/frqTP, powTP, label='TP')
pl.loglog(pixscale12m.to(u.arcsec).value/frq12m_conv, pow12m_conv, label='12m->7m')
pl.loglog(pixscale7m.to(u.arcsec).value/frq7mTP, pow7mTP, label='7m+TP')
pl.loglog(pixscale12m.to(u.arcsec).value/frq12m7mTP, pow12m7mTP, label='7m+12m+TP')
pl.legend(loc='best')
pl.xlabel("Angular Scale (arcsec)")
pl.ylabel("Power Spectrum (K$^2$)")
pl.savefig("array_pspeccomparison.png")

pl.figure(2).clf()
ax1 = pl.subplot(2,3,1)
ax1.set_title("7m")
ax1.imshow(slc7m.value*jtok_7m[ch_65kms_7m], cmap='hot', vmin=-5, vmax=50)
ax2 = pl.subplot(2,3,2)
ax2.set_title("12m")
ax2.imshow(slc12m.value*jtok_12m[ch_65kms_12m], cmap='hot', vmin=-5, vmax=50)
ax3 = pl.subplot(2,3,3)
ax3.set_title("TP")
ax3.imshow(slcTP, cmap='hot', vmin=-5, vmax=50)
ax4 = pl.subplot(2,3,4)
ax4.set_title("12m convolved to 7m")
ax4.imshow(slc12m_conv_to_7m*jtok_12m[ch_65kms_12m], cmap='hot', vmin=-5, vmax=50)
ax5 = pl.subplot(2,3,5)
ax5.set_title("7m feathered with TP")
ax5.imshow(hdu7mTP.data, cmap='hot', vmin=-5, vmax=50)
ax6 = pl.subplot(2,3,6)
ax6.set_title("12m feathered with 7mTP")
ax6.imshow(hdu12m7mTP.data, cmap='hot', vmin=-5, vmax=50)
pl.savefig("array_imagecomparison.png")


# Weighting Scale Tests
scale20as = uvcombine.feather_simple(hdu12m, hdu7mTP, lowresfwhm=20*u.arcsec)
scale15as = uvcombine.feather_simple(hdu12m, hdu7mTP, lowresfwhm=15*u.arcsec)
scale10as = uvcombine.feather_simple(hdu12m, hdu7mTP, lowresfwhm=10*u.arcsec)
pl.figure(3).clf()
ax1 = pl.subplot(2,2,1)
ax1.set_title("Default: {0:0.1f}arcsec".format(slc7m.header['BMAJ']*3600))
ax1.imshow(hdu12m7mTP.data, cmap='hot', vmin=-5, vmax=50)
ax2 = pl.subplot(2,2,2)
ax2.set_title("20arcsec")
ax2.imshow(scale20as.real, cmap='hot', vmin=-5, vmax=50)
ax3 = pl.subplot(2,2,3)
ax3.set_title("15arcsec")
ax3.imshow(scale15as.real, cmap='hot', vmin=-5, vmax=50)
ax4 = pl.subplot(2,2,4)
ax4.set_title("10arcsec")
ax4.imshow(scale10as.real, cmap='hot', vmin=-5, vmax=50)
pl.savefig("combination_weighting_scales.png")

# Weighting Scale Tests round 2
scale70as_no7m = uvcombine.feather_simple(hdu12m, tp_fn, lowresfwhm=70*u.arcsec)
scale50as_no7m = uvcombine.feather_simple(hdu12m, tp_fn, lowresfwhm=50*u.arcsec)
scale30as_no7m = uvcombine.feather_simple(hdu12m, tp_fn, lowresfwhm=30*u.arcsec)
scale10as_no7m = uvcombine.feather_simple(hdu12m, tp_fn, lowresfwhm=10*u.arcsec)
pl.figure(4).clf()
ax1 = pl.subplot(2,2,1)
ax1.set_title("70arcsec")
ax1.imshow(scale70as_no7m.real, cmap='hot', vmin=-5, vmax=50)
ax2 = pl.subplot(2,2,2)
ax2.set_title("50arcsec")
ax2.imshow(scale50as_no7m.real, cmap='hot', vmin=-5, vmax=50)
ax3 = pl.subplot(2,2,3)
ax3.set_title("30arcsec")
ax3.imshow(scale30as_no7m.real, cmap='hot', vmin=-5, vmax=50)
ax4 = pl.subplot(2,2,4)
ax4.set_title("10arcsec")
ax4.imshow(scale10as_no7m.real, cmap='hot', vmin=-5, vmax=50)
pl.savefig("combination_weighting_scales_no7m.png")



pl.draw()
pl.show()
