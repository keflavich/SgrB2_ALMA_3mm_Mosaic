import os
from astropy import units as u
from astropy.io import fits
import uvcombine
import radio_beam
import pylab as pl

def pth(x):
    return os.path.join('/Volumes/passport/alma/sgrb2_b3/merge', x)

inthdu = fits.open(pth('singleframes/HC3N_natural_cubek_65kms.fits'))
sdhdu = fits.open(pth('singleframes/HC3N_natural_tpcube_k_spatialandspectralregrid_65kms.fits'))
sdbm = radio_beam.Beam.from_fits_header(sdhdu[0].header)

pl.figure(1).clf()
uvcombine.uvcombine.feather_compare(inthdu[0], sdhdu[0], SAS=60*u.arcsec,
                                    LAS=120*u.arcsec,
                                    lowresfwhm=sdbm.major)
pl.figure(2).clf()
uvcombine.uvcombine.feather_plot(inthdu[0], sdhdu[0],
                                 lowresfwhm=sdbm.major)

pl.figure(3).clf()
uvcombine.uvcombine.feather_plot(inthdu[0], sdhdu[0],
                                 lowresfwhm=sdbm.major,
                                 xaxisunit='lambda',
                                )

comb1 = uvcombine.uvcombine.feather_simple(inthdu[0],
                                           sdhdu[0],
                                           lowresfwhm=sdbm.major,
                                           highpassfilterSD=False,
                                           replace_hires=False,
                                           deconvSD=False,
                                           return_hdu=True)
comb1.writeto(pth('singleframes/HC3N_natural_feather_FFF.fits'), clobber=True)

comb2 = uvcombine.uvcombine.feather_simple(inthdu[0],
                                           sdhdu[0],
                                           lowresfwhm=sdbm.major,
                                           highpassfilterSD=True,
                                           replace_hires=False,
                                           deconvSD=False,
                                           return_hdu=True)
comb2.writeto(pth('singleframes/HC3N_natural_feather_convLoRes.fits'), clobber=True)

comb3 = uvcombine.uvcombine.feather_simple(inthdu[0],
                                           sdhdu[0],
                                           lowresfwhm=sdbm.major,
                                           highpassfilterSD=False,
                                           replace_hires=0.95,
                                           deconvSD=False,
                                           return_hdu=True)
comb3.writeto(pth('singleframes/HC3N_natural_feather_replaceHiRes.fits'), clobber=True)

comb4 = uvcombine.uvcombine.feather_simple(inthdu[0],
                                           sdhdu[0],
                                           lowresfwhm=sdbm.major,
                                           highpassfilterSD=False,
                                           replace_hires=False,
                                           deconvSD=True,
                                           return_hdu=True)
comb4.writeto(pth('singleframes/HC3N_natural_feather_deconvLoRes.fits'), clobber=True)

pl.figure(4)
pl.subplot(2,2,1).imshow(comb1.data, vmin=-5, vmax=50, cmap=pl.cm.gray)
pl.subplot(2,2,2).imshow(comb2.data, vmin=-5, vmax=50, cmap=pl.cm.gray)
pl.subplot(2,2,3).imshow(comb3.data, vmin=-5, vmax=50, cmap=pl.cm.gray)
pl.subplot(2,2,4).imshow(comb4.data, vmin=-5, vmax=50, cmap=pl.cm.gray)
