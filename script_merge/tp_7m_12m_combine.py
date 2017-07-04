import numpy as np
import os
import glob
from spectral_cube import SpectralCube
# from singledish_combine import spectral_regrid, feather_simple, fourier_combine_cubes
from uvcombine import feather_simple, spectral_regrid, spectral_smooth_and_downsample
from astropy.io import fits
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy import log
import re
import socket

if 'nmpost' in socket.gethostname():
    dpath = lambda x: os.path.join("/lustre/aoc/users/aginsbur/sgrb2/2013.1.00269.S/merge",x)
    tppath = lambda x: os.path.join("/lustre/aoc/users/aginsbur/sgrb2/2013.1.00269.S/tp",x)
else:
    raise ValueError("No match to socket hostname {0}.".format(socket.gethostname()))



speciesre = re.compile('SgrB2_b3_7M_12M([_a-z]*).([-a-zA-Z0-9]*)(\.r[0-9\.]*)?.image.pbcor.fits')

# '../tp/tp_concat.spw17.image.fits: 90345430335.98326 Hz,92220007685.55771 Hz'
# '../tp/tp_concat.spw19.image.fits: 88541267991.9948 Hz,90415845341.33504 Hz'
# '../tp/tp_concat.spw21.image.fits: 100427226138.55649 Hz,102302047605.2752 Hz'
# '../tp/tp_concat.spw23.image.fits: 102287973625.33496 Hz,104162550974.6752 Hz'

velocity_ranges = {'HC3N': (-150,-25),
                   'CH3CN': (-100,-25),
                   'HCN': (-150,-120),
                   'HNC': (-150,-120),
                   'HCOp': (-150,-120),
                   'H41a': (-150,-25),
                  }

for interferometer_fn in (
    'SgrB2_b3_7M_12M.HCOp.r2.image.pbcor.fits',
    "SgrB2_b3_7M_12M.HCN.r2.image.pbcor.fits",
    'SgrB2_b3_7M_12M.HNC.r2.image.pbcor.fits',
    'SgrB2_b3_7M_12M.CH3CN.r2.image.pbcor.fits',
    'SgrB2_b3_7M_12M.H41a.r2.image.pbcor.fits',
    'SgrB2_b3_7M_12M.HCOp.r0.5.image.pbcor.fits',
    "SgrB2_b3_7M_12M.HCN.r0.5.image.pbcor.fits",
    'SgrB2_b3_7M_12M.HNC.r0.5.image.pbcor.fits',
    'SgrB2_b3_7M_12M.CH3CN.r0.5.image.pbcor.fits',
    'SgrB2_b3_7M_12M.H41a.r0.5.image.pbcor.fits',
    #"SgrB2_b3_7M_12M_natural.CH3CN.image.pbcor.fits",
    #'SgrB2_b3_7M_12M.HC3N.image.pbcor.fits',
    #"SgrB2_b3_7M_12M_natural.H2CS303-202.image.pbcor.fits",
    #"SgrB2_b3_7M_12M_natural.H2CS321-220.image.pbcor.fits",
    #"SgrB2_b3_7M_12M_natural.H2CS322-221.image.pbcor.fits",
    #"SgrB2_b3_7M_12M_natural.H41a.image.pbcor.fits",
    #"SgrB2_b3_7M_12M_natural.HC3N.image.pbcor.fits",
    #"SgrB2_b3_7M_12M_natural.HCN.image.pbcor.fits",
    #"SgrB2_b3_7M_12M_natural.HCOp.image.pbcor.fits",
    #"SgrB2_b3_7M_12M_natural.HNC.image.pbcor.fits",
):

    suffix,species,robust = speciesre.search(interferometer_fn).groups()
    outfilename = ('{species}{suffix}{robust}_TP_7m_12m_feather.fits'
                   .format(species=species, suffix=suffix, robust=robust))

    medsubfn = ('{species}{suffix}{robust}_7m_12m_medsub.fits'
                .format(species=species, suffix=suffix, robust=robust))
    if not os.path.exists(medsubfn):
        cube = (SpectralCube.read(dpath(interferometer_fn))
                .with_spectral_unit(u.km/u.s, velocity_convention='radio'))
        cube.beam_threshold = 100
        # try to avoid contamination... won't work universally; need to examine
        # individual cubes and have this as a parameter
        #med = cube.spectral_slab(90*u.km/u.s, 160*u.km/u.s).median(axis=0).value
        med = cube.spectral_slab(velocity_ranges[species][0]*u.km/u.s,
                                 velocity_ranges[species][1]*u.km/u.s).median(axis=0).value
        cube.write(medsubfn, overwrite=True)
        fh = fits.open(medsubfn, mode='update')
        log.info("Median subtracting")
        pb = ProgressBar(len(fh[0].data))
        for ii,imslice in enumerate(fh[0].data):
            fh[0].data[ii] = imslice - med
            fh.flush()
            pb.update()
        fh.close()

    cube = SpectralCube.read(medsubfn).with_spectral_unit(u.GHz)

    if hasattr(cube, 'beam'):
        avgbeam = cube.beam
        jtok = cube.beam.jtok(cube.spectral_axis).value
        print("Jansky/beam -> Kelvin factor = {0}".format(jtok))
    else:
        jtok = np.array([bm.jtok(x).value for bm,x in zip(cube.beams,
                                                          cube.spectral_axis)])
        avgbeam = cube.average_beams(0.1)
        print("Median Jansky/beam -> Kelvin factor = {0}".format(np.median(jtok)))
        print("Average Jansky/beam -> Kelvin factor = {0}"
              .format(avgbeam.jtok(cube.spectral_axis.mean())))

    minghz, maxghz = cube.spectral_extrema

    OK = False
    for fn in glob.glob(tppath("tp_concat*fits")):
        tpcube = SpectralCube.read(fn)
        if minghz > tpcube.spectral_extrema[0] and maxghz < tpcube.spectral_extrema[1]:
            OK = True
            break
    if not OK:
        raise ValueError("No matching TP cube")


    crop_channels = sorted((tpcube.closest_spectral_channel(minghz),
                            tpcube.closest_spectral_channel(maxghz)))
    tpcube = tpcube[crop_channels[0]-1:crop_channels[1]+1]
    log.info("Read tp freq")
    tpcube_k = tpcube.to(u.K, tpcube.beam.jtok_equiv(tpcube.spectral_axis))
    log.info("Converted TP to K")
    # determine smooth factor kw = kernel width
    kw = np.abs((cube.spectral_axis.diff().mean() /
                 tpcube_k.spectral_axis.diff().mean()).decompose().value)
    log.info("determined kernel = {0}".format(kw))

    tpcube_k_ds_hdu = spectral_smooth_and_downsample(tpcube_k, kw)

    tpdscube = SpectralCube.read(tpcube_k_ds_hdu)
    tpkrg_hdu = spectral_regrid(tpdscube, cube.spectral_axis)
    log.info("done regridding")
    #tpkrg.writeto('HC3N_tp_freq_ds_interp.fits', clobber=True)

    cube_tpkrg = SpectralCube.read(tpkrg_hdu) #'HC3N_tp_freq_ds_interp.fits')
    tpoutfn = '{species}{suffix}_TP.fits'.format(species=species, suffix=suffix)
    cube_tpkrg.with_spectral_unit(u.km/u.s,
                                  rest_value=cube.wcs.wcs.restfrq*u.Hz,
                                  velocity_convention='radio').write(tpoutfn,
                                                                     overwrite=True)
    print("Single dish beam = {0}".format(cube_tpkrg.beam))

    # intermediate work: test that a single frame has been properly combined
    frq = cube.wcs.wcs.restfrq*u.Hz * (1-65/3e5) # approximately 65 kms
    closestchan = cube.closest_spectral_channel(frq)
    imJy = cube[closestchan]
    imK = imJy.to(u.K, cube.beams[closestchan].jtok_equiv(frq))
    sdim = cube_tpkrg[cube_tpkrg.closest_spectral_channel(frq)]
    imJy.write('singleframes/{species}{suffix}{robust}_cubeJy_65kms.fits'.format(species=species, suffix=suffix, robust=robust), overwrite=True)
    imK.write('singleframes/{species}{suffix}{robust}_cubek_65kms.fits'.format(species=species, suffix=suffix, robust=robust), overwrite=True)
    sdim.write('singleframes/{species}{suffix}{robust}_tpcube_k_rg_65kms.fits'.format(species=species, suffix=suffix, robust=robust), overwrite=True)

    # for sanity checking purposes, write out the Jy version too
    sdim_Jy = sdim.to(u.Jy, tpcube.beam.jtok_equiv(frq))
    sdim_Jy._unit = u.Jy/u.beam
    sdim_Jy.write('singleframes/{species}{suffix}{robust}_tpcube_Jy_rg_65kms.fits'
                  .format(species=species, suffix=suffix, robust=robust), overwrite=True)

    imKhdu = imK.hdu
    #imKhdu.data[imKhdu.data<0] = 0 # DELETE THIS IT'S WRONG
    combohdu, hdu2 = feather_simple(imKhdu, sdim.hdu,
                                    lowresscalefactor=1.0,
                                    highresscalefactor=1.0,
                                    highpassfilterSD=False,
                                    replace_hires=False,
                                    deconvSD=False,
                                    return_regridded_lores=True,
                                    return_hdu=True)
    combohdu.header.update(cube.header)
    combohdu.header.update(avgbeam.to_header_keywords())
    combohdu.writeto('singleframes/{species}{suffix}{robust}_TP_7m_12m_feather_65kms.fits'
                     .format(species=species, suffix=suffix, robust=robust), clobber=True)
    hdu2.header['BMAJ'] = tpcube_k_ds_hdu.header['BMAJ']
    hdu2.header['BMIN'] = tpcube_k_ds_hdu.header['BMIN']
    hdu2.header['BPA'] = tpcube_k_ds_hdu.header['BPA']
    hdu2.writeto('singleframes/{species}{suffix}{robust}_tpcube_k_spatialandspectralregrid_65kms.fits'.format(species=species, suffix=suffix, robust=robust), clobber=True)

    combohdu, hdu2 = feather_simple(imK.hdu, sdim.hdu, return_regridded_lores=True, return_hdu=True, replace_hires=0.5)
    combohdu.writeto('singleframes/{species}{suffix}{robust}_TP_7m_12m_replacefeather_65kms.fits'.format(species=species, suffix=suffix, robust=robust), clobber=True)

    if 'do_full_cube' not in locals() or do_full_cube:
        # final goal
        assert os.system('cp {0} {1}'.format(interferometer_fn, outfilename)) == 0
        assert len(cube) == len(cube_tpkrg)
        if not np.sign(cube.spectral_axis.diff()[0]) == np.sign(cube_tpkrg.spectral_axis.diff()[0]):
            cube_tpkrg = cube_tpkrg[::-1]
        fh = fits.open(outfilename, mode='update')
        fh[0].header['BUNIT'] = 'K'

        # first set everything to zero to make sure it's actually working
        pb = ProgressBar(len(fh[0].data))
        for ii,layer in (enumerate(fh[0].data)):
            fh[0].data[ii,:,:] = 0
            fh.flush()
            pb.update()

        pb = ProgressBar(len(cube))

        for ii,(im, sdim, jtok_int) in enumerate(zip(cube, cube_tpkrg, jtok)):
            imhdu = im.hdu
            if im.unit == u.Jy:
                imhdu.data = imhdu.data*jtok_int
                imhdu.header['BUNIT'] = 'K'
            comb = feather_simple(imhdu, sdim.hdu,
                                  lowresscalefactor=1.0,
                                  highresscalefactor=0.1,
                                  highpassfilterSD=False,
                                  replace_hires=False,
                                  deconvSD=False,
                                  # divide by 2 because I think the kernel being used now is too harsh on the TP data
                                  #lowresfwhm=cube_tpkrg.beam.major/2.,
                                  # actually, I'm going to leave the low res fwhm at its default and do no downweighting
                                 )
            fh[0].data[ii,:,:] = comb.real
            fh.flush()
            pb.update()

        fh.close()
            
        #combhdu = fourier_combine_cubes(cube_k, cube_tpkrg, return_hdu=True,
        #                                lowresfwhm=cube_tpkrg.beam.major)
        #combcube = SpectralCube.read(combhdu).with_spectral_unit(u.km/u.s,
        #                                                         velocity_convention='radio')
        #combcube.write(outfilename,
        #               overwrite=True)
