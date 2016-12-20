from __future__ import print_function
from astropy import units as u
import numpy as np
from spectral_cube import SpectralCube, BooleanArrayMask
import os
import glob

#cube = SpectralCube.read('SgrB2_b3_12M.HCOp.flux.fits')
#rendermask_arr = (cube[0] > 0.75)
#rendermask = BooleanArrayMask(rendermask_arr, cube.wcs, cube.shape)
#rcube = cube.with_mask(rendermask).minimal_subcube()
#os.chdir('/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic/FITS/merge/')
#os.chdir('/Volumes/passport/alma/sgrb2_b3/merge/')
os.chdir('/lustre/aginsbur/sgrb2/2013.1.00269.S/merge')

#files = [
#         'SgrB2_b3_7M_12M.CH3CN.image.pbcor_medsub.fits',
#         #'SgrB2_b3_7M_12M.H41a.image.pbcor_medsub.fits',
#         #'SgrB2_b3_7M_12M.HC3N.image.pbcor_medsub.fits',
#         'SgrB2_b3_7M_12M.HCN.image.pbcor_medsub.fits',
#         #'SgrB2_b3_7M_12M.HCOp.image.pbcor_medsub.fits',
#         'SgrB2_b3_7M_12M.HNC.image.pbcor_medsub.fits',
#]

files = glob.glob("SgrB2_b3_7M_12M*.image.pbcor.fits")

for outdir in ('max','moment0', 'moment1', 'median'):
    if not os.path.exists(outdir):
        os.mkdir(outdir)

for fn in files:
    pfx = os.path.splitext(fn)[0]

    if (os.path.exists('max/{0}_max.fits'.format(pfx)) and
        os.path.exists('moment0/{0}_moment0_mask.fits'.format(pfx))
       ):
        pass
        #print("Skipping {0}".format(pfx))
        #continue

    cube = SpectralCube.read(fn).minimal_subcube().with_spectral_unit(u.km/u.s, velocity_convention='radio')
    cube.beam_threshold = 1000
    print(fn,cube)
    m0 = cube.moment0(axis=0)
    max = cube.max(axis=0)
    m0.hdu.writeto('moment0/{0}_moment0.fits'.format(pfx), clobber=True)
    max.hdu.writeto('max/{0}_max.fits'.format(pfx), clobber=True)

    std = cube.std(axis=0)
    mask = cube > std
    mcube = cube.with_mask(mask)
    m0mask = mcube.moment0(axis=0)
    m1mask = mcube.moment1(axis=0)
    m0mask.hdu.writeto('moment0/{0}_moment0_mask.fits'.format(pfx), clobber=True)
    m1mask.hdu.writeto('moment1/{0}_moment1_mask.fits'.format(pfx), clobber=True)

    med = cube.median(axis=0)
    med.hdu.writeto('median/{0}_median.fits'.format(pfx), clobber=True)

    cube.allow_huge_operations=True
    medsub_cube = cube - med
    medmax = medsub_cube.max(axis=0)
    medmax.hdu.writeto('max/{0}_max_medsub.fits'.format(pfx), clobber=True)
    medm0 = medsub_cube.moment0(axis=0)
    medm0.hdu.writeto('moment0/{0}_moment0_medsub.fits'.format(pfx), clobber=True)
