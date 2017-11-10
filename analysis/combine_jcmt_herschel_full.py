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
import dust_emissivity
import pylab as pl
from astropy import log

pl.rcParams['figure.dpi'] = 75
pl.rcParams['figure.figsize'] = (12,8)

scuba_fn = 'gc450.fits.gz'
herschel_fn = 'igls_l000_plw_deglitch_hh.fits'
herschel_fn = 'gcmosaic_500um.fits'

jcmt_hdu = fits.open(scuba_fn)[0]

### FIX NANS in Herschel ###
herschel_hdu = fits.open(herschel_fn)[0]

# find JCMT corners
jcmtwcs = wcs.WCS(jcmt_hdu.header)
herschel_wcs = wcs.WCS(herschel_hdu.header)
lllb = coordinates.SkyCoord(*jcmtwcs.wcs_pix2world(0,0,0), unit=(u.deg, u.deg), frame='fk5').galactic
urlb = coordinates.SkyCoord(*jcmtwcs.wcs_pix2world(jcmt_hdu.data.shape[1],
                            jcmt_hdu.data.shape[0],
                            0), unit=(u.deg, u.deg), frame='fk5').galactic
lll, llb = lllb.l, lllb.b
url,urb = urlb.l, urlb.b

def rnd(x):
    return int(np.round(x))

urx,ury = map(rnd, herschel_wcs.wcs_world2pix(url,urb,0))
llx,lly = map(rnd, herschel_wcs.wcs_world2pix(lll,llb,0))

header = herschel_hdu.header.copy()
header.update(herschel_wcs[lly:ury, llx:urx].to_header())
herschel_cutout_hdu = fits.PrimaryHDU(data=herschel_hdu.data[lly:ury, llx:urx],
                                      header=header)


new_header = herschel_cutout_hdu.header.copy()
upsample_factor = 2.
new_header['CDELT1'] /= upsample_factor
new_header['CDELT2'] /= upsample_factor
new_header['NAXIS1'] *= upsample_factor
new_header['NAXIS2'] *= upsample_factor
new_header['NAXIS1'] = int(new_header['NAXIS1'])
new_header['NAXIS2'] = int(new_header['NAXIS2'])
new_header['CRPIX1'] = (new_header['CRPIX1']-1)*upsample_factor+1
new_header['CRPIX2'] = (new_header['CRPIX2']-1)*upsample_factor+1
new_header['BUNIT'] = 'MJy/sr'
#upsampled_data = upsample_image(herschel_sgrb2_cutout_hdu.data,
#                                upsample_factor=upsample_factor).real
upsampled_data = reproject.reproject_interp(herschel_cutout_hdu,
                                            new_header)[0]
herschel_upsampled_hdu = fits.PrimaryHDU(data=upsampled_data,
                                         header=new_header)
herschel_upsampled_hdu.writeto('herschel_500um_upsampled_fullCMZ.fits', overwrite=True)

### FIX crap in JCMT ###
data450 = jcmt_hdu.data
kernel = Gaussian2DKernel(1.5, x_size=15).array
kernel /= kernel.sum()
sm_jcmt = convolve_fft(data450, kernel, nan_treatment='interpolate')
data450[441:448,1377:1384] = sm_jcmt[441:448,1377:1384]
data450[432:438,1390:1398] = sm_jcmt[432:438,1390:1398]
data450[479:485,1400:1407] = sm_jcmt[479:485,1400:1407]
data450[470:476,1413:1420] = sm_jcmt[470:476,1413:1420]

# remove very large scales
kernel2 = Gaussian2DKernel(25)
sm_jcmt2 = convolve_fft(data450, kernel2, nan_treatment='interpolate', fft_pad=False)

jcmt_hdu.data = data450-sm_jcmt2
jcmt_hdu.writeto('scuba450_fixed_unsharp.fits', overwrite=True)



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
scuba_hdu.writeto("scuba_shifted_fullCMZ.fits", overwrite=True)

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
scuba_hdu.writeto("scuba_shifted_fullCMZ_MJySr.fits", overwrite=True)



# scalefactor = 3 forces a match near 100"
plot_rslt = uvcombine.uvcombine.feather_plot('scuba_shifted_fullCMZ_MJySr.fits',
                                             'herschel_500um_upsampled_fullCMZ.fits',
                                             lowresfwhm=45*u.arcsec,
                                             #highpassfilterSD='reconvolve',
                                             #deconvSD=True,
                                             highresscalefactor=3)

pl.figure(1).clf()
compare_results = uvcombine.uvcombine.feather_compare('scuba_shifted_fullCMZ_MJySr.fits',
                                                      'herschel_500um_upsampled_fullCMZ.fits',
                                                      SAS=1.5*45*u.arcsec,
                                                      LAS=120*u.arcsec,
                                                      min_beam_fraction=0.1,
                                                      lowresfwhm=45*u.arcsec)
pl.savefig("scuba_herschel_scaling_comparison_fullCMZ.png")
highresscalefactor = 1./compare_results['median_sc']
log.info("Measured high-res scalefactor (deconvolved): {0}".format(highresscalefactor))

pl.figure(7)
compare_results = uvcombine.uvcombine.feather_compare('scuba_shifted_fullCMZ_MJySr.fits',
                                                      'herschel_500um_upsampled_fullCMZ.fits',
                                                      SAS=1.5*45*u.arcsec,
                                                      LAS=120*u.arcsec,
                                                      lowresfwhm=45*u.arcsec,
                                                      min_beam_fraction=0.1,
                                                      beam_divide_lores=False)
highresscalefactor = 1./compare_results['median_sc']
log.info("Measured high-res scalefactor (not deconvolved): {0}".format(highresscalefactor))
pl.savefig("scuba_herschel_scaling_comparison_nodeconv_fullCMZ.png")


result = uvcombine.uvcombine.feather_simple('scuba_shifted_fullCMZ_MJySr.fits',
                                            'herschel_500um_upsampled_fullCMZ.fits',
                                            lowresfwhm=45*u.arcsec,
                                            # scale factor b/c we think JCMT
                                            # is miscalibrated?
                                            highresscalefactor=highresscalefactor,
                                           )

scuba_hdu.data = result.real
scuba_hdu.writeto("scuba_Herschel_Feathered_fullCMZ.fits", overwrite=True)
log.info("Wrote full CMZ feathered image")


herschel_tem = fits.open('gcmosaic_temp_conv25.fits')[0]
herschel_tem_regrid,_ = reproject.reproject_interp(herschel_tem, scuba_hdu.header)
herschel_tem_smoothed = convolve_fft(herschel_tem_regrid,
                                     Gaussian2DKernel(35))
herschel_tem_regrid[np.isnan(herschel_tem_regrid)] = herschel_tem_smoothed[np.isnan(herschel_tem_regrid)]

colmap_scuba_herscheltem = dust_emissivity.dust.colofsnu((500*u.um).to(u.GHz, u.spectral()),
                                                         u.Quantity(scuba_hdu.data,
                                                                    unit=u.Unit(scuba_hdu.header['BUNIT'])),
                                                         temperature=u.Quantity(herschel_tem_regrid, u.K),
                                                         beta=1.750,)

colheader = header=scuba_hdu.header.copy()
colheader['BUNIT'] = 'cm^-2'
colmap_hdu = fits.PrimaryHDU(data=colmap_scuba_herscheltem.to(u.cm**-2).value,
                             header=colheader)
colmap_hdu.writeto('SCUBA_Herschel_feather_columndensity_interpolatedHerscheltem.fits',
                   overwrite=True)


pl.close(2)
pl.figure(2, figsize=(18,9)).clf()
pl.subplot(2,3,1)
vmin = np.nanmin(herschel_upsampled_hdu.data)
vmax = np.nanmax(herschel_upsampled_hdu.data)
pl.imshow(herschel_upsampled_hdu.data[np.isfinite(herschel_upsampled_hdu.data)], cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("Herschel 500um")

pl.subplot(2,3,4)
_,_,(L,) = pl.hist(herschel_upsampled_hdu.data[np.isfinite(herschel_upsampled_hdu.data)],
                   bins=50, log=True, histtype='step', linewidth=3, alpha=0.75)
pl.hist(herschel_upsampled_hdu.data[np.isfinite(herschel_upsampled_hdu.data)],
        bins=50, log=True, histtype='step', linewidth=1, zorder=-5,
        color=L.get_edgecolor())
pl.hist(fits.getdata('scuba_shifted_MJySr.fits').ravel()*highresscalefactor,
        bins=50, log=True, histtype='step')
pl.hist(scuba_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.subplot(2,3,2)
pl.imshow(fits.getdata('scuba_shifted_MJySr.fits')*highresscalefactor,
          cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("SCUBA2 500um (scaled)")

pl.subplot(2,3,5)
pl.hist(herschel_upsampled_hdu.data.ravel(), bins=50, log=True, histtype='step')
_,_,(L,) = pl.hist(fits.getdata('scuba_shifted_MJySr.fits').ravel()*highresscalefactor,
        bins=50, log=True, histtype='step', linewidth=3, alpha=0.75, zorder=10)
pl.hist(fits.getdata('scuba_shifted_MJySr.fits').ravel()*highresscalefactor,
        bins=50, log=True, histtype='step', linewidth=1, zorder=-5,
        color=L.get_edgecolor())
pl.hist(scuba_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.subplot(2,3,3)
pl.imshow(scuba_hdu.data, cmap='gray', norm=pl.matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
pl.colorbar()
pl.title("Feathered 500um")

pl.subplot(2,3,6)
pl.hist(herschel_upsampled_hdu.data.ravel(), bins=50, log=True, histtype='step')
pl.hist(fits.getdata('scuba_shifted_MJySr.fits').ravel()*highresscalefactor,
        bins=50, log=True, histtype='step')
_,_,(L,) = pl.hist(scuba_hdu.data.ravel(), bins=50, log=True, histtype='step',
                   linewidth=3, alpha=0.75, zorder=10)
pl.hist(scuba_hdu.data.ravel(), bins=50, log=True, histtype='step',
        linewidth=1, zorder=-5, color=L.get_edgecolor())
pl.semilogx()
pl.yscale('log', nonposy='clip')

pl.savefig("scuba_Herschel_images_feathered_fullCMZ.png")
