from __future__ import print_function
from astropy import units as u
import numpy as np
from spectral_cube import SpectralCube, BooleanArrayMask
import glob
import os

cube = SpectralCube.read('SgrB2_b3_12M.HCOp.flux.fits')
rendermask_arr = (cube[0] > 0.75)
rendermask = BooleanArrayMask(rendermask_arr, cube.wcs, cube.shape)
rcube = cube.with_mask(rendermask).minimal_subcube()

files = [
'SgrB2_b3_12M.HC3N.image.pbcor.fits',
'SgrB2_b3_12M.HNC.image.pbcor.fits',
'SgrB2_b3_12M.H41a.image.pbcor.fits',
'SgrB2_b3_12M.CH3CN.image.pbcor.fits',
'SgrB2_b3_12M.H2CO615-616.image.pbcor.fits',
'SgrB2_b3_12M.H2CS303-202.image.pbcor.fits',
'SgrB2_b3_12M.HCN.image.pbcor.fits',
'SgrB2_b3_12M.HCOp.image.pbcor.fits',
'SgrB2_a_03_12M.HCOp.image.pbcor.fits',
'SgrB2_a_03_12M.HCN.image.pbcor.fits',
]

for fn in files:
    pfx = os.path.splitext(fn)[0]
    if 'contsub' in pfx:
        continue

    if (os.path.exists('max/{0}_max_contsub.fits'.format(pfx)) and
        os.path.exists('{0}.contsub.fits'.format(pfx)) and
        os.path.exists('moment0/{0}_moment0_mask_contsub.fits'.format(pfx))
       ):
        pass
        #print("Skipping {0}".format(pfx))
        #continue

    cube = SpectralCube.read(fn).with_spectral_unit(u.km/u.s, velocity_convention='radio')
    print(fn,cube)
    m0 = cube.moment0(axis=0)
    max = cube.max(axis=0)
    med = cube.median(axis=0)
    m0.hdu.writeto('moment0/{0}_moment0.fits'.format(pfx), clobber=True)
    max.hdu.writeto('max/{0}_max.fits'.format(pfx), clobber=True)
    med.hdu.writeto('med/{0}_med.fits'.format(pfx), clobber=True)

    cube._data -= med.value
    m0 = cube.moment0(axis=0)
    m1 = cube.moment1(axis=0)
    max = cube.max(axis=0)
    m0.hdu.writeto('moment0/{0}_moment0_contsub.fits'.format(pfx), clobber=True)
    m1.hdu.writeto('moment1/{0}_moment1_contsub.fits'.format(pfx), clobber=True)
    max.hdu.writeto('max/{0}_max_contsub.fits'.format(pfx), clobber=True)
    cube.write('{0}.contsub.fits'.format(pfx), overwrite=True)

    std = cube.std(axis=0)
    mask = cube > std
    mcube = cube.with_mask(mask)
    m0mask = mcube.moment0(axis=0)
    m1mask = mcube.moment1(axis=0)
    m0mask.hdu.writeto('moment0/{0}_moment0_mask_contsub.fits'.format(pfx), clobber=True)
    m1mask.hdu.writeto('moment1/{0}_moment1_mask_contsub.fits'.format(pfx), clobber=True)

    rendermask = BooleanArrayMask(rendermask_arr, mcube.wcs, mcube.shape)
    rcube = mcube.with_mask(rendermask).minimal_subcube()
    print("Minimal rendering cube: ",fn,rcube)
    ytc = rcube.to_yt()
    ytc.quick_isocontour(level=3*np.nanmean(std.value),
                         title='Sgr B2 12mTC {0}'.format(pfx),
                         description='ALMA observations of Sgr B2.  File {0}'.format(pfx),
                         filename=pfx+".ply")


    del mask
    del cube

