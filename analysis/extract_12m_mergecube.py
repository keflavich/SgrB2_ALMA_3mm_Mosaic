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
elif socket.gethostname() == 'nergal':
    dpath = lambda x: os.path.join("/scratch/aginsbur/SgrB2/ALMA/2013.1.00269.S/merge/fullcube/",x)
    rpath = lambda x: os.path.join("/data/SgrB2/ALMA/2013.1.00269.S/regions",x)
    spath = lambda x: os.path.join("/data/SgrB2/ALMA/2013.1.00269.S/spectra",x)
    mergecubes = [
        'full_SgrB2_TETC7m_r2_spw0_lines.fits',
        'full_SgrB2_TETC7m_r2_spw1_lines.fits',
        'full_SgrB2_TETC7m_r2_spw2_lines.fits',
        'full_SgrB2_TETC7m_r2_spw3_lines.fits',
    ]
elif 'nmpost' in socket.gethostname():
    dpath = lambda x: os.path.join("/lustre/aoc/users/aginsbur/sgrb2/2013.1.00269.S/merge/fits",x)
    rpath = lambda x: os.path.join("/lustre/aoc/users/aginsbur/sgrb2/2013.1.00269.S/SgrB2_ALMA_3mm_Mosaic/regions",x)
    spath = lambda x: os.path.join("/lustre/aoc/users/aginsbur/sgrb2/2013.1.00269.S/spectra",x)
    mergecubes = [
        'SgrB2_b3_7M_12M.CFp.image.pbcor.fits',
        'SgrB2_b3_7M_12M.CH3CN.image.pbcor.fits',
        'SgrB2_b3_7M_12M.CH3OH7m26-716.image.pbcor.fits',
        'SgrB2_b3_7M_12M.H15NC.image.pbcor.fits',
        'SgrB2_b3_7M_12M.H2CO615-616.image.pbcor.fits',
        'SgrB2_b3_7M_12M.H2CS303-202.image.pbcor.fits',
        'SgrB2_b3_7M_12M.H2CS313-212.image.pbcor.fits',
        'SgrB2_b3_7M_12M.H2CS321-220.image.pbcor.fits',
        'SgrB2_b3_7M_12M.H2CS322-221.image.pbcor.fits',
        'SgrB2_b3_7M_12M.H41a.image.pbcor.fits',
        'SgrB2_b3_7M_12M.HC3N.image.pbcor.fits',
        'SgrB2_b3_7M_12M.HCN.image.pbcor.fits',
        'SgrB2_b3_7M_12M.HCOp.image.pbcor.fits',
        'SgrB2_b3_7M_12M.HNC.image.pbcor.fits',
    ]
else:
    raise ValueError("No match to socket hostname {0}.".format(socket.gethostname()))


regions = (
           pyregion.open(rpath('cores_with_names.reg'))
           pyregion.open(rpath('tc_continuum_core_extraction_regions.reg')) +
           pyregion.open(rpath('ionizationfront_circle.reg')) +
           pyregion.open(rpath('extraction_regions_n_and_m.reg')) +
           pyregion.open(rpath('ch3cn_large_cores.reg')) +
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

        bgSL = pyregion.parse("fk5; circle({0},{1},{2}\")"
                              .format(CL[0],
                                      CL[1],
                                      2*radius*3600))
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
