import os
from astropy import log
from astropy import units as u
import numpy as np
from spectral_cube import SpectralCube
import pyregion
import radio_beam
import socket

if socket.gethostname() == 'cyg':
    dpath = lambda x: os.path.join("/Volumes/passport/alma/sgrb2_b3/12m_tc/",x)
    rpath = lambda x: os.path.join("/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic/regions/",x)
    spath = lambda x: os.path.join("/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic/data/spectra/",x)
    mergecubes = [
                  'SgrB2_12m_spw0_lines.fits',
                  'SgrB2_12m_spw1_lines.fits',
                  'SgrB2_12m_spw2_lines.fits',
                  'SgrB2_12m_spw3_lines.fits',
    ]
elif socket.gethostname() == 'shamash':
    dpath = lambda x: os.path.join("/data/SgrB2/ALMA/2013.1.00269.S/center",x)
    rpath = lambda x: os.path.join("/data/SgrB2/ALMA/2013.1.00269.S/regions",x)
    spath = lambda x: os.path.join("/data/SgrB2/ALMA/2013.1.00269.S/spectra",x)
    mergecubes = [
        'full_SgrB2_12m_r-0.5_spw0_lines.fits',
        'full_SgrB2_12m_r-0.5_spw1_lines.fits',
        'full_SgrB2_12m_r-0.5_spw2_lines.fits',
        'full_SgrB2_12m_r-0.5_spw3_lines.fits',
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
        try:
            scube = cube.subcube_from_ds9region(reg)
        except ValueError as ex:
            print("Skipping {0} because {1}".format(name, ex))
            continue
        print(cube)
        log.info("Source name: {0}  filename: {1}".format(name,fname))
        print(scube)
        spectrum = scube.mean(axis=(1,2))
        spectrum.meta['beam'] = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in spectrum.beams]),
                                                minor=np.nanmedian([bm.minor.to(u.deg).value for bm in spectrum.beams]),
                                                pa=np.nanmedian([bm.pa.to(u.deg).value for bm in spectrum.beams]),
                                               )

        hdu = spectrum.hdu
        pixel_scale = np.abs(cube.wcs.celestial.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
        hdu.header['PPBEAM'] = (spectrum.meta['beam'].sr / pixel_scale**2).decompose().value

        hdu.writeto(spath("{1}_{0}.fits".format(suffix,fname)), clobber=True)
        print(spath("{1}_{0}.fits".format(suffix,fname)))
