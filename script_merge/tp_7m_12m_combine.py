import numpy as np
import os
import glob
import FITS_tools
from spectral_cube import SpectralCube
from singledish_combine import spectral_regrid, feather_simple, fourier_combine_cubes
from astropy.io import fits
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy import log
import re

speciesre = re.compile('SgrB2_b3_7M_12M.([a-zA-Z0-9]*).image.pbcor.fits')

# '../tp/tp_concat.spw17.image.fits: 90345430335.98326 Hz,92220007685.55771 Hz'
# '../tp/tp_concat.spw19.image.fits: 88541267991.9948 Hz,90415845341.33504 Hz'
# '../tp/tp_concat.spw21.image.fits: 100427226138.55649 Hz,102302047605.2752 Hz'
# '../tp/tp_concat.spw23.image.fits: 102287973625.33496 Hz,104162550974.6752 Hz'

for interferometer_fn in (
    'SgrB2_b3_7M_12M.HCOp.image.pbcor.fits',
    'SgrB2_b3_7M_12M.HNC.image.pbcor.fits',
    'SgrB2_b3_7M_12M.CH3CN.image.pbcor.fits',
    'SgrB2_b3_7M_12M.H41a.image.pbcor.fits',
    'SgrB2_b3_7M_12M.HC3N.image.pbcor.fits',
):

    species = speciesre.search(interferometer_fn).groups()[0]
    outfilename = '{species}_TP_7m_12m_feather.fits'.format(species=species)

    medsubfn = '{species}_7m_12m_medsub.fits'.format(species=species)
    if not os.path.exists(medsubfn):
        cube = SpectralCube.read(interferometer_fn).with_spectral_unit(u.km/u.s,
                                                                       velocity_convention='radio')
        cube.beam_threshold = 100
        # try to avoid contamination... won't work universally; need to examine
        # individual cubes and have this as a parameter
        med = cube.spectral_slab(90*u.km/u.s, 160*u.km/u.s).median(axis=0).value
        os.system('cp {0} {1}'.format(interferometer_fn, medsubfn))
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
        jtok = cube.beam.jtok(cube.spectral_axis).value
    else:
        jtok = np.array([bm.jtok(x).value for bm,x in zip(cube.beams,
                                                          cube.spectral_axis)])

    minghz, maxghz = cube.spectral_extrema

    OK = False
    for fn in glob.glob("../tp/tp_concat*fits"):
        tpcube = SpectralCube.read(fn)
        if minghz > tpcube.spectral_extrema[0] and maxghz < tpcube.spectral_extrema[1]:
            OK = True
            break
    if not OK:
        raise ValueError("No matching TP cube")


    crop_channels = sorted((tpcube.closest_spectral_channel(minghz),
                            tpcube.closest_spectral_channel(maxghz)))
    tpcube = tpcube[crop_channels[0]-1:crop_channels[1]+1]
    log.debug("Read tp freq")
    tpcube_k = tpcube.to(u.K, tpcube.beam.jtok_equiv(tpcube.spectral_axis[:,None,None]))
    log.debug("Converted TP to K")
    # determine smooth factor
    kw = (cube.spectral_axis.diff().mean() / tpcube_k.spectral_axis.diff().mean()).decompose().value
    log.debug("determined kernel")
    tbcube_k_smooth = FITS_tools.cube_regrid.spectral_smooth_cube(tpcube_k,
                                                                  kw/np.sqrt(8*np.log(2)))
    log.debug("completed cube smooth")

    integer_dsfactor = int(np.floor(kw))

    tpcube_k_ds = tbcube_k_smooth[::integer_dsfactor,:,:]
    log.debug("downsampled")
    print(tpcube_k)
    print(tpcube_k.hdu)
    tpcube_k.hdu
    log.debug("did nothing")
    tpcube_k_ds_hdu = tpcube_k.hdu
    log.debug("made hdu")
    tpcube_k_ds_hdu.data = tpcube_k_ds
    log.debug("put data in hdu")
    tpcube_k_ds_hdu.header['CRPIX3'] = 1
    # why min? because we're forcing CDELT3 to be positive, therefore the 0'th channel
    # must be the reference value.  Since we're using a symmetric kernel to downsample,
    # the reference channel - wherever we pick it - must stay fixed.
    tpcube_k_ds_hdu.header['CRVAL3'] = tpcube_k.spectral_axis[0].to(u.Hz).value
    tpcube_k_ds_hdu.header['CUNIT3'] = tpcube_k.spectral_axis[0].to(u.Hz).unit.to_string('FITS')
    tpcube_k_ds_hdu.header['CDELT3'] = tpcube_k.wcs.wcs.cdelt[2] * integer_dsfactor
    log.debug("completed header making")
    #tpcube_k_ds_hdu.writeto('{0}_tp_freq_ds.fits', clobber=True)
    #log.debug("wrote kelvin tp downsampled cube")

    tpdscube = SpectralCube.read(tpcube_k_ds_hdu)
    tpkrg = spectral_regrid(tpdscube,
                            cube.spectral_axis)
    log.debug("done regridding")
    #tpkrg.writeto('HC3N_tp_freq_ds_interp.fits', clobber=True)

    cube_tpkrg = SpectralCube.read(tpkrg) #'HC3N_tp_freq_ds_interp.fits')
    tpoutfn = '{species}_TP.fits'.format(species=species)
    cube_tpkrg.with_spectral_unit(u.km/u.s,
                                  rest_value=cube.wcs.wcs.restfrq*u.Hz,
                                  velocity_convention='radio').write(tpoutfn,
                                                                     overwrite=True)
    print("Single dish beam = {0}".format(cube_tpkrg.beam))

    # intermediate work: test that a single frame has been properly combined
    frq = cube.wcs.wcs.restfrq*u.Hz * (1-65/3e5) # approximately 65 kms
    closestchan = cube.closest_spectral_channel(frq)
    im = cube[closestchan] * jtok[closestchan]
    sdim = cube_tpkrg[cube_tpkrg.closest_spectral_channel(frq)]
    im.write('singleframes/{species}_cubek_65kms.fits'.format(species=species), overwrite=True)
    sdim.write('singleframes/{species}_tpcube_k_rg_65kms.fits'.format(species=species), overwrite=True)
    combohdu, hdu2 = feather_simple(im.hdu, sdim.hdu, return_regridded_lores=True, return_hdu=True)
    combohdu.writeto('singleframes/{species}_TP_7m_12m_feather_65kms.fits'.format(species=species), clobber=True)

    if do_full_cube:
        # final goal
        os.system('cp {0} {1}'.format(interferometer_fn, outfilename))
        assert len(cube) == len(cube_tpkrg)
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
            imhdu.data *= jtok_int
            comb = feather_simple(im.hdu, sdim.hdu,
                                  lowresfwhm=cube_tpkrg.beam.major)
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
