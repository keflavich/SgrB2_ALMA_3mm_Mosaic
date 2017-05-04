import numpy as np
import paths
from astropy.io import fits
import pylab as pl
from astropy import wcs
import astropy.visualization
from astropy.nddata import Cutout2D
from astropy import coordinates
from astropy import units as u
import regions

noise_regions = regions.read_ds9(paths.rpath('noise_regions.reg'))

sc5tt = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy.image.tt0.pbcor.fits'))
mywcs = wcs.WCS(sc5tt[0].header)

data = sc5tt[0].data
peak = np.nanmax(data)
print("Peak flux: {0} Jy".format(peak))

for reg in noise_regions:

    pixreg = reg.to_pixel(mywcs)
    mask = pixreg.to_mask()
    cutout = mask.cutout(data) * mask.data

    name = reg.meta['text'].strip("{}")

    print("{0}: std={1:0.2f} mJy  dr={2:0.2f}".format(name, cutout.std()*1000, peak/cutout.std()))
