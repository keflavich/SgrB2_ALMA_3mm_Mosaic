import numpy as np
import warnings
from astropy import units as u
from astropy.modeling import models, fitting
from radio_beam import Beam
from astropy.io import fits
import regions
from astropy import wcs
from astropy.stats import mad_std
from astropy.utils.console import ProgressBar
import pylab as pl
import paths
import signal, sys, time

def signal_handler(signal, frame):
    sys.exit(0)

STDDEV_TO_FWHM = np.sqrt(8*np.log(2))

def gaussfit_catalog(fitsfile, region_list, radius=1.0*u.arcsec,
                     max_radius_in_beams=2,
                     background_estimator=np.nanmedian,
                     noise_estimator=lambda x: mad_std(x, ignore_nan=True)):

    fh = fits.open(fitsfile)
    data = fh[0].data
    header = fh[0].header
    datawcs = wcs.WCS(header)
    beam = Beam.from_fits_header(header)
    pixscale = wcs.utils.proj_plane_pixel_area(datawcs)**0.5 * u.deg
    bmmaj_px = (beam.major.to(u.deg) / pixscale).decompose()

    noise = noise_estimator(data)

    fit_data = {}

    for reg in ProgressBar(region_list):

        phot_reg = regions.CircleSkyRegion(center=reg.center, radius=radius)
        pixreg = phot_reg.to_pixel(datawcs)
        mask = pixreg.to_mask()
        cutout = mask.cutout(data) * mask.data

        smaller_phot_reg = regions.CircleSkyRegion(center=reg.center,
                                                   radius=beam.major/STDDEV_TO_FWHM)
        smaller_pixreg = smaller_phot_reg.to_pixel(datawcs)
        smaller_mask = smaller_pixreg.to_mask()
        smaller_cutout = smaller_mask.cutout(data) * smaller_mask.data

        background = np.nanmedian(cutout)

        sz = cutout.shape[0]
        mx = np.nanmax(smaller_cutout)
        ampguess = mx-background

        p_init = models.Gaussian2D(amplitude=ampguess,
                                   x_mean=sz/2,
                                   y_mean=sz/2,
                                   x_stddev=bmmaj_px/STDDEV_TO_FWHM,
                                   y_stddev=bmmaj_px/STDDEV_TO_FWHM,
                                   bounds={'x_stddev':(bmmaj_px/STDDEV_TO_FWHM*0.75,
                                                       bmmaj_px*max_radius_in_beams/STDDEV_TO_FWHM),
                                           'y_stddev':(bmmaj_px/STDDEV_TO_FWHM*0.75,
                                                       bmmaj_px*max_radius_in_beams/STDDEV_TO_FWHM),
                                           'x_mean':(sz/2-2, sz/2+2),
                                           'y_mean':(sz/2-2, sz/2+2),
                                           'amplitude':(ampguess*0.9, ampguess*1.1)
                                          }
                                  )

        result, fit_info, chi2 = gaussfit_image(image=(cutout-background)*mask.data,
                                                gaussian=p_init,
                                                weights=1/noise**2,
                                                plot=True
                                               )
        sourcename = reg.meta['text'].strip('{}')
        pl.savefig(paths.fpath('gaussfits/{0}.png'.format(sourcename)),
                   bbox_inches='tight')

        cx,cy = pixreg.bounding_box.ixmin+result.x_mean, pixreg.bounding_box.iymin+result.y_mean
        clon,clat = datawcs.wcs_pix2world(cx, cy, 0)
        fit_data[sourcename] = {'amplitude': result.amplitude,
                                'center_x': float(clon)*u.deg,
                                'center_y': float(clat)*u.deg,
                                'fwhm_x': result.x_stddev * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                                'fwhm_y': result.x_stddev * STDDEV_TO_FWHM * pixscale.to(u.arcsec),
                                'pa': (result.theta*u.rad).to(u.deg),
                               }

        signal.signal(signal.SIGINT, signal_handler)

    return fit_data


def gaussfit_image(image, gaussian, weights=None,
                   fitter=fitting.LevMarLSQFitter(), plot=False):

    yy, xx = np.mgrid[:image.shape[0], :image.shape[1]]
    
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        fitted = fitter(gaussian, xx, yy, image, weights=weights)

    fitim = fitted(xx,yy)
    residual = image-fitim
    residualsquaredsum = np.nansum(residual**2*weights)

    if plot:
        pl.clf()
        ax1 = pl.subplot(2,2,1)
        ax1.imshow(image, cmap='viridis', origin='lower',
                   interpolation='nearest')
        ax2 = pl.subplot(2,2,2)
        ax2.imshow(fitim, cmap='viridis', origin='lower',
                   interpolation='nearest')
        ax3 = pl.subplot(2,2,3)
        ax3.imshow(residual, cmap='viridis', origin='lower',
                   interpolation='nearest')
        ax4 = pl.subplot(2,2,4)
        ax4.hist(residual[(residual!=0) & (image!=0)], bins=50)
    
    return fitted, fitter.fit_info, residualsquaredsum

if __name__ == "__main__":
    regs = regions.read_ds9(paths.rpath('cores_with_names.reg'))

    from files import contfilename as contfnpath

    fit_data = gaussfit_catalog(contfnpath, regs)
