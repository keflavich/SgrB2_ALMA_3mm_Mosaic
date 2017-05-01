from astropy import units as u
from astropy import wcs
from astropy import coordinates
import radio_beam
import numpy as np
import uvcombine
import image_registration
from astropy.io import fits
from astropy.convolution import convolve_fft, Gaussian2DKernel
from astropy.nddata.utils import Cutout2D
import pylab as pl

pl.rcParams['figure.dpi'] = 75
pl.rcParams['figure.figsize'] = (12,8)

sharc_fn = 'SgrB2_350um_gal.fits'
herschel_fn = 'igls_l000_pmw_deglitch_hh.fits'

### FIX NANS in Herschel ###
herschel_hdu = fits.open(herschel_fn)
data350 = herschel_hdu[0].data
data350[1091,590] = 185000
kernel = Gaussian2DKernel(0.5, x_size=15).array
kernel /= kernel.sum()
sm_herschel = convolve_fft(data350, kernel, interpolate_nan=True)
badmask = np.isnan(data350) | (data350 == 0)
data350[badmask] = sm_herschel[badmask]
herschel_hdu[0].data = data350
herschel_nonan_fn = 'igls_l000_pmw_deglitch_hh_nonan.fits'
herschel_hdu.writeto(herschel_nonan_fn, overwrite=True)
herschel_sgrb2_cutout = Cutout2D(herschel_hdu[0].data, coordinates.SkyCoord.from_name("Sgr B2(M)"), 20*u.arcmin,
                                 wcs=wcs.WCS(herschel_hdu[0].header))
herschel_sgrb2_cutout_hdu = fits.PrimaryHDU(data=herschel_sgrb2_cutout.data,
                                            header=herschel_sgrb2_cutout.wcs.to_header())
herschel_sgrb2_cutout_hdu.writeto("herschel_350um_cutout.fits", overwrite=True)
new_header = herschel_sgrb2_cutout_hdu.header.copy()
new_header['CDELT1'] /= 4.
new_header['CDELT2'] /= 4.
new_header['NAXIS1'] *= 4.
new_header['NAXIS2'] *= 4.
new_header['CRPIX1'] = (new_header['CRPIX1']-1)*4.+1
new_header['CRPIX2'] = (new_header['CRPIX2']-1)*4.+1
new_header['BUNIT'] = 'MJy/sr'
herschel_upsampled_hdu = fits.PrimaryHDU(data=image_registration.fft_tools.upsample_image(herschel_sgrb2_cutout_hdu.data,
                                                                                          upsample_factor=4).real,
                                         header=new_header)
herschel_upsampled_hdu.writeto('herschel_350um_upsampled.fits', overwrite=True)


### REGISTER SHARC TO HERSCHEL ###
#result = image_registration.cross_correlation_shifts_FITS(
#    herschel_upsampled_hdu, sharc_fn, return_cropped_images=True,
#    register_method=image_registration.chi2_shift_iterzoom, verbose=True)
result = image_registration.FITS_tools.match_images.register_fits(
    herschel_upsampled_hdu, sharc_fn, return_cropped_images=True,
    return_error=False, return_shifted_image=True,
    register_method=image_registration.chi2_shift_iterzoom, verbose=True)
xoff,yoff,xoff_wcs,yoff_wcs,herschel_reproj,sharc_crop,sharc_reproj = result
sharc_hdu = herschel_upsampled_hdu.copy()
sharc_hdu.data = sharc_reproj
sharc_hdu.writeto("sharc_shifted.fits", overwrite=True)

### SHARC SCALE FACTOR ###
sharc_beam_fwhm = 9*u.arcsec # from Bally+2010
sharc_beam_fwhm = 11.5*u.arcsec # guess
sharc_beam = radio_beam.Beam(sharc_beam_fwhm)
sharc_data_Jybeam = sharc_hdu.data*u.Jy / sharc_beam
sharc_data_MJySr = sharc_data_Jybeam.to(u.MJy/u.sr)
print("Factor = {0}".format((1*u.Jy / sharc_beam).to(u.MJy/u.sr)))
sharc_hdu.data = sharc_data_MJySr.value
sharc_hdu.header['BUNIT'] = 'MJy/sr'
sharc_hdu.writeto("sharc_shifted_MJySr.fits", overwrite=True)



# scalefactor was ~1/0.0035, pretty close to the 9"-> 1/0.00215, but closer to 11"
# SHARC looks really good past a Herschel primary beam, so selected a larger beam
uvcombine.uvcombine.feather_plot('sharc_shifted_MJySr.fits',
                                 'herschel_350um_upsampled.fits',
                                 lowresfwhm=45*u.arcsec,
                                 highresscalefactor=1.0)

pl.figure(1).clf()
compare_results = uvcombine.uvcombine.feather_compare('sharc_shifted_MJySr.fits',
                                                      'herschel_350um_upsampled.fits',
                                                      SAS=1.5*30*u.arcsec,
                                                      LAS=120*u.arcsec,
                                                      min_beam_fraction=0.7,
                                                      lowresfwhm=28*u.arcsec)
pl.savefig("sharc_herschel_scaling_comparison.png")
highresscalefactor = 1./compare_results['median_sc']
print("Theoretical high-res scalefactor: {0}".format(highresscalefactor))
# actually, we force the high-res scale factor to 1, because we should
# be weighting the larger angular scales much more than the smaller.
compare_results = uvcombine.uvcombine.feather_compare('sharc_shifted_MJySr.fits',
                                                      'herschel_350um_upsampled.fits',
                                                      SAS=1.5*30*u.arcsec,
                                                      LAS=120*u.arcsec,
                                                      lowresfwhm=28*u.arcsec,
                                                      min_beam_fraction=0.7,
                                                      beam_divide_lores=False)
pl.savefig("sharc_herschel_scaling_comparison_nodeconv.png")

# straight feathering, without using replace_hires, resulted in negative bowls.
# also, use deconvSD because it's in a high S/N range and otherwise we have
# some weird smoothing effects in place
# (of course, with deconvSD, it is important that we get the lowresFWHM right)
result = uvcombine.uvcombine.feather_simple('sharc_shifted_MJySr.fits',
                                            'herschel_350um_upsampled.fits',
                                            lowresfwhm=25*u.arcsec,
                                            highresscalefactor=1.0,
                                            replace_hires=0.1,
                                            deconvSD=True,
                                           )


sharc_hdu.data = result.real
sharc_hdu.writeto("SHARC_Herschel_Feathered.fits", overwrite=True)

pl.close(2)
pl.figure(2, figsize=(18,9)).clf()
pl.subplot(2,3,1)
vmin = herschel_upsampled_hdu.data.min()
vmax = herschel_upsampled_hdu.data.max()
pl.imshow(herschel_upsampled_hdu.data, cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("Herschel 350um")

pl.subplot(2,3,4)
pl.hist(herschel_upsampled_hdu.data.ravel(), bins=50, log=True, histtype='step',
        linewidth=3)
pl.hist(fits.getdata('sharc_shifted_MJySr.fits').ravel(), bins=50, log=True, histtype='step')
pl.hist(sharc_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.subplot(2,3,2)
pl.imshow(fits.getdata('sharc_shifted_MJySr.fits'), cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("SHARC2 350um")

pl.subplot(2,3,5)
pl.hist(herschel_upsampled_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.hist(fits.getdata('sharc_shifted_MJySr.fits').ravel(), bins=50, log=True, histtype='step',
        linewidth=3)
pl.hist(sharc_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.subplot(2,3,3)
pl.imshow(sharc_hdu.data, cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("Feathered 350um")

pl.subplot(2,3,6)
pl.hist(herschel_upsampled_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.hist(fits.getdata('sharc_shifted_MJySr.fits').ravel(), bins=50, log=True, histtype='step')
pl.hist(sharc_hdu.data.ravel(), bins=50, log=True, histtype='step',
        linewidth=3)
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.savefig("SHARC_Herschel_images_feathered.png")


pl.figure(3).clf()
pl.suptitle("Full images side-by-side comparison")
pl.subplot(2,1,1)
pl.loglog(herschel_upsampled_hdu.data.ravel(), sharc_data_MJySr.value.ravel(), ',', alpha=0.25)
pl.loglog([1e3,1e6], [1e3,1e6], 'k--')
pl.xlabel('Herschel')
pl.ylabel('Sharc')
pl.axis([1e3,1e6,1e3,1e6])
pl.subplot(2,1,2)
pl.loglog(herschel_upsampled_hdu.data.ravel(), result.real.ravel(), ',', alpha=0.25)
pl.loglog([1e3,1e6], [1e3,1e6], 'k--')
pl.xlabel('Herschel')
pl.ylabel('Feathered')
pl.axis([1e3,1e6,1e3,1e6])
pl.savefig("SHARC_Herschel_before-after_flux_comparison.png")
