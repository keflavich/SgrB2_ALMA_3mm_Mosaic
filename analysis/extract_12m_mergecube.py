import os
from astropy import log
import numpy as np
from spectral_cube import SpectralCube
import pyregion

dpath = lambda x: os.path.join("/Volumes/passport/alma/sgrb2_b3/12m_tc/",x)
rpath = lambda x: os.path.join("/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic/regions/",x)
spath = lambda x: os.path.join("/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic/data/spectra/",x)

mergecubes = [
'SgrB2_12m_spw0_lines.fits',
'SgrB2_12m_spw1_lines.fits',
'SgrB2_12m_spw2_lines.fits',
'SgrB2_12m_spw3_lines.fits',
]


regions = (pyregion.open(rpath('tc_continuum_core_extraction_regions.reg')) +
           pyregion.open(rpath('ionizationfront_circle.reg')) +
           pyregion.open(rpath('extraction_regions_n_and_m.reg')) +
           pyregion.open(rpath('ch3cn_large_cores.reg'))
          )

for cubename in mergecubes:
    for reg in regions:
        name = reg.attr[1]['text']
        fname = name.replace(" ","_").lower()
        reg = pyregion.ShapeList([reg])

        suffix = os.path.splitext(cubename)[0]
        if os.path.exists(spath("{1}_{0}.fits".format(suffix,fname))):
            continue

        cube = SpectralCube.read(dpath(cubename))
        print(cube)
        log.info(name)
        log.info(fname)
        scube = cube.subcube_from_ds9region(reg)
        print(scube)
        spectrum = scube.mean(axis=(1,2))

        spectrum.hdu.writeto(spath("{1}_{0}.fits".format(suffix,fname)),
                             clobber=True)
        print(spath("{1}_{0}.fits".format(suffix,fname)))
