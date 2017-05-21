import os
from astropy import log
from astropy import units as u
import numpy as np
from spectral_cube import SpectralCube
import pyregion
import radio_beam
import socket

if 'nmpost' in socket.gethostname():
    dpath = lambda x: os.path.join("/lustre/aoc/users/aginsbur/sgrb2/2013.1.00269.S/merge/fullcubes_r0",x)
    rpath = lambda x: os.path.join("/lustre/aoc/users/aginsbur/sgrb2/2013.1.00269.S/SgrB2_ALMA_3mm_Mosaic/regions",x)
    spath = lambda x: os.path.join("/lustre/aoc/users/aginsbur/sgrb2/2013.1.00269.S/fullcube_spectra",x)
    mergecubes = [
        'full_SgrB2_TETC7m_r0_spw0_lines.fits',
        'full_SgrB2_TETC7m_r0_spw1_lines.fits',
        'full_SgrB2_TETC7m_r0_spw2_lines.fits',
        'full_SgrB2_TETC7m_r0_spw3_lines.fits',
    ]
else:
    raise ValueError("No match to socket hostname {0}.".format(socket.gethostname()))


regions = (
           pyregion.open(rpath('cores_with_names.reg'))
           #pyregion.open(rpath('tc_continuum_core_extraction_regions.reg')) +
           #pyregion.open(rpath('ionizationfront_circle.reg')) +
           #pyregion.open(rpath('extraction_regions_n_and_m.reg')) +
           #pyregion.open(rpath('ch3cn_large_cores.reg'))
          )

for cubename in mergecubes:
    for reg in regions:
        name = reg.attr[1]['text']
        fname = name.replace(" ","_").lower()
        CL = reg.coord_list
        if reg.name == 'circle':
            radius = CL[2]
            reg = pyregion.ShapeList([reg])
        else:
            radius = 0.5
            reg = pyregion.parse("fk5; circle({0},{1},0.5\")"
                                 .format(CL[0], CL[1]))

        suffix = os.path.splitext(cubename)[0]
        if os.path.exists(spath("{1}_{0}.fits".format(suffix,fname))):
            print("Skipping {0} {1} because it exists".format(suffix, fname))
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
        # I think this is a hack left over from old versions of SpectralCube
        spectrum.meta['beam'] = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in spectrum.beams]),
                                                minor=np.nanmedian([bm.minor.to(u.deg).value for bm in spectrum.beams]),
                                                pa=np.nanmedian([bm.pa.to(u.deg).value for bm in spectrum.beams]),
                                               )

        hdu = spectrum.hdu
        pixel_scale = np.abs(cube.wcs.celestial.pixel_scale_matrix.diagonal().prod())**0.5 * u.deg
        hdu.header['PPBEAM'] = (spectrum.meta['beam'].sr / pixel_scale**2).decompose().value

        hdu.header['OBJECT'] = name

        hdu.writeto(spath("{1}_{0}.fits".format(suffix,fname)), clobber=True)
        print(spath("{1}_{0}.fits".format(suffix,fname)))

        bgSL = pyregion.parse("fk5; circle({0},{1},{2}\")"
                              .format(CL[0],
                                      CL[1],
                                      2*radius))
        bgsc = cube.subcube_from_ds9region(bgSL)
        npix = np.count_nonzero(np.isfinite(bgsc[0,:,:]))
        bgspec = (bgsc.sum(axis=(1,2)) - scube.sum(axis=(1,2))) / npix
        bgspec.meta['beam'] = radio_beam.Beam(major=np.nanmedian([bm.major.to(u.deg).value for bm in spectrum.beams]),
                                              minor=np.nanmedian([bm.minor.to(u.deg).value for bm in spectrum.beams]),
                                              pa=np.nanmedian([bm.pa.to(u.deg).value for bm in spectrum.beams]),
                                             )
        bgspec.hdu.writeto(spath("{0}_{1}_background_mean{2}.fits".format(fname,
                                                                          name,
                                                                          suffix)
                                ), clobber=True)
