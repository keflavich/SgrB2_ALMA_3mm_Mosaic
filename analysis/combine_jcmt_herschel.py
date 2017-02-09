from astropy import units as u
from astropy import wcs
from astropy import coordinates
import reproject
import radio_beam
import numpy as np
import uvcombine
import image_registration
from image_registration.fft_tools import upsample_image
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
from astropy.convolution import convolve_fft, Gaussian2DKernel
import pylab as pl

scuba_fn = 'gc450.fits.gz'
herschel_fn = 'igls_l000_plw_deglitch_hh.fits'

### FIX NANS in Herschel ###
herschel_hdu = fits.open(herschel_fn)
herschel_sgrb2_cutout = Cutout2D(herschel_hdu[0].data,
                                 coordinates.SkyCoord.from_name("Sgr B2(M)"), 20*u.arcmin,
                                 wcs=wcs.WCS(herschel_hdu[0].header))
herschel_sgrb2_cutout_hdu = fits.PrimaryHDU(data=herschel_sgrb2_cutout.data,
                                            header=herschel_sgrb2_cutout.wcs.to_header())
herschel_sgrb2_cutout_hdu.writeto("herschel_500um_cutout.fits", overwrite=True)
new_header = herschel_sgrb2_cutout_hdu.header.copy()
upsample_factor = 5.
new_header['CDELT1'] /= upsample_factor
new_header['CDELT2'] /= upsample_factor
new_header['NAXIS1'] *= upsample_factor
new_header['NAXIS2'] *= upsample_factor
new_header['NAXIS1'] = int(new_header['NAXIS1'])
new_header['NAXIS2'] = int(new_header['NAXIS2'])
new_header['CRPIX1'] = (new_header['CRPIX1']-1)*upsample_factor+1
new_header['CRPIX2'] = (new_header['CRPIX2']-1)*upsample_factor+1
#upsampled_data = upsample_image(herschel_sgrb2_cutout_hdu.data,
#                                upsample_factor=upsample_factor).real
upsampled_data = reproject.reproject_interp(herschel_sgrb2_cutout_hdu,
                                            new_header)[0]
herschel_upsampled_hdu = fits.PrimaryHDU(data=upsampled_data,
                                         header=new_header)
herschel_upsampled_hdu.writeto('herschel_500um_upsampled.fits', overwrite=True)

### FIX crap in JCMT ###
jcmt_hdu = fits.open(scuba_fn)
data450 = jcmt_hdu[0].data
kernel = Gaussian2DKernel(1.5, x_size=15).array
kernel /= kernel.sum()
sm_jcmt = convolve_fft(data450, kernel, interpolate_nan=True)
data450[441:448,1377:1384] = sm_jcmt[441:448,1377:1384]
data450[432:438,1390:1398] = sm_jcmt[432:438,1390:1398]
data450[479:485,1400:1407] = sm_jcmt[479:485,1400:1407]
data450[470:476,1413:1420] = sm_jcmt[470:476,1413:1420]
jcmt_hdu[0].data = data450
jcmt_hdu.writeto('scuba450_fixed.fits', overwrite=True)



### REGISTER scuba TO HERSCHEL ###
#result = image_registration.cross_correlation_shifts_FITS(
#    herschel_upsampled_hdu, scuba_fn, return_cropped_images=True,
#    register_method=image_registration.chi2_shift_iterzoom, verbose=True)
result = image_registration.FITS_tools.match_images.register_fits(
    herschel_upsampled_hdu, jcmt_hdu, return_cropped_images=True,
    return_error=False, return_shifted_image=True,
    register_method=image_registration.chi2_shift_iterzoom, verbose=True)
xoff,yoff,xoff_wcs,yoff_wcs,herschel_reproj,scuba_crop,scuba_reproj = result
scuba_hdu = herschel_upsampled_hdu.copy()
scuba_hdu.data = scuba_reproj
scuba_hdu.writeto("scuba_shifted.fits", overwrite=True)

### scuba SCALE FACTOR ###
scuba_beam_fwhm = 8*u.arcsec # from pierce-price+
#scuba_beam_fwhm = 9*u.arcsec
#scuba_beam_fwhm = 11.5*u.arcsec
scuba_beam = radio_beam.Beam(scuba_beam_fwhm)
scuba_data_Jybeam = scuba_hdu.data*u.Jy / scuba_beam
scuba_data_MJySr = scuba_data_Jybeam.to(u.MJy/u.sr)
print("Factor = {0}".format((1*u.Jy / scuba_beam).to(u.MJy/u.sr)))
scuba_hdu.data = scuba_data_MJySr.value
scuba_hdu.header['BUNIT'] = 'MJy/sr'
scuba_hdu.writeto("scuba_shifted_MJySr.fits", overwrite=True)



# scalefactor = 3 forces a match near 100"
plot_rslt = uvcombine.uvcombine.feather_plot('scuba_shifted_MJySr.fits',
                                             'herschel_500um_upsampled.fits',
                                             lowresfwhm=45*u.arcsec,
                                             #highpassfilterSD='reconvolve',
                                             #deconvSD=True,
                                             highresscalefactor=3)

pl.figure(1).clf()
compare_results = uvcombine.uvcombine.feather_compare('scuba_shifted_MJySr.fits',
                                                      'herschel_500um_upsampled.fits',
                                                      SAS=1.5*45*u.arcsec,
                                                      LAS=120*u.arcsec,
                                                      min_beam_fraction=0.7,
                                                      lowresfwhm=45*u.arcsec)
pl.savefig("scuba_herschel_scaling_comparison.png")
highresscalefactor = 1./compare_results['median_sc']
print("Theoretical high-res scalefactor: {0}".format(highresscalefactor))
compare_results = uvcombine.uvcombine.feather_compare('scuba_shifted_MJySr.fits',
                                                      'herschel_500um_upsampled.fits',
                                                      SAS=1.5*45*u.arcsec,
                                                      LAS=120*u.arcsec,
                                                      lowresfwhm=45*u.arcsec,
                                                      min_beam_fraction=0.7,
                                                      beam_divide_lores=False)
pl.savefig("scuba_herschel_scaling_comparison_nodeconv.png")


result = uvcombine.uvcombine.feather_simple('scuba_shifted_MJySr.fits',
                                            'herschel_500um_upsampled.fits',
                                            lowresfwhm=45*u.arcsec,
                                            # scale factor b/c we think JCMT
                                            # is miscalibrated?
                                            highresscalefactor=3,
                                           )

scuba_hdu.data = result.real
scuba_hdu.writeto("scuba_Herschel_Feathered.fits", overwrite=True)


pl.close(2)
pl.figure(2, figsize=(18,9)).clf()
pl.subplot(2,3,1)
vmin = herschel_upsampled_hdu.data.min()
vmax = herschel_upsampled_hdu.data.max()
pl.imshow(herschel_upsampled_hdu.data, cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("Herschel 500um")

pl.subplot(2,3,4)
pl.hist(herschel_upsampled_hdu.data.ravel(), bins=50, log=True, histtype='step',
        linewidth=3)
pl.hist(fits.getdata('scuba_shifted_MJySr.fits').ravel(), bins=50, log=True, histtype='step')
pl.hist(scuba_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.subplot(2,3,2)
pl.imshow(fits.getdata('scuba_shifted_MJySr.fits'), cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("scuba2 500um")

pl.subplot(2,3,5)
pl.hist(herschel_upsampled_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.hist(fits.getdata('scuba_shifted_MJySr.fits').ravel(), bins=50, log=True, histtype='step',
        linewidth=3)
pl.hist(scuba_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.subplot(2,3,3)
pl.imshow(scuba_hdu.data, cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("Feathered 500um")

pl.subplot(2,3,6)
pl.hist(herschel_upsampled_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.hist(fits.getdata('scuba_shifted_MJySr.fits').ravel(), bins=50, log=True, histtype='step')
pl.hist(scuba_hdu.data.ravel(), bins=50, log=True, histtype='step',
        linewidth=3)
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.savefig("scuba_Herschel_images_feathered.png")
