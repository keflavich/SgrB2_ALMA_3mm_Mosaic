"""
http://docs.astropy.org/en/stable/io/fits/appendix/faq.html#how-can-i-create-a-very-large-fits-file-from-scratch
"""
from astropy import log
from astropy.io import fits
from astropy import wcs
import numpy as np
import glob
import re
import os
from astropy.utils.console import ProgressBar
from spectral_cube import SpectralCube

nchans_total = {0: 7680, 1: 7680, 2: 7680, 3: 7680}
min_nchans = 7425
frange = {0: [90357.27912, 92220.858],
          1: [88552.69612, 90417.088],
          2: [100440.526,102304.10488],
          3: [102301.456,104164.21788],
         }
fstep = {0:250., # kHz
         1:250., # kHz
         2:250., # kHz
         3:250., # kHz
        }

# Extract the appropriate pixel indices from the file name.
# A more sophisticated approach is probably better, in which the individual
# cubes are inspected for their start/end frequencies.
# But, on the other hand, for this process to make any sense at all, you
# have to have done the original cube imaging right
def getinds(fn):
    inds = re.search('channels([0-9]*)to([0-9]*)', fn).groups()
    return [int(ii) for ii in inds]

def get_max_ind(globstr):
    # replace nchans_total with the correct version from the actual data on
    # disk
    files = glob.glob(globstr)
    maxind = max([max(getinds(fn)) for fn in files])
    return maxind

def make_spw_cube(spw='spw{0}', spwnum=0, fntemplate='SgrB2',
                  overwrite_existing=False, bmaj_limits=None,
                  fnsuffix="", filesuffix='image.pbcor.fits',
                  first_endchannel='*',
                  cropends=False,
                  minimize=True,
                  debug_mode=False,
                  add_beam_info=True):
    """
    Parameters
    ----------
    spw : str
        String template for the input/output name
    spwnum : int
        The spectral window number
    fntemplate : str
        Filename template (goes into the glob)
    overwrite_existing : bool
        Overwrite data in the output cube?
    cropends: bool or int
        Number of pixels to crop off the ends of an image
    minimize: bool
        Compute the spatial minimal subcube before building the cube?  Slices
        for all subsequent cubes will be computed from the first cube.
    """
    if debug_mode:
        lvl = log.getEffectiveLevel()
        log.setLevel('DEBUG')

    spw = spw.format(spwnum)

    big_filename = '{1}_{0}{2}_lines.fits'.format(spw, fntemplate, fnsuffix)

    header_fn = glob.glob('piece_of_{1}_cube{2}.{0}.channels0to{4}.{3}'
                          .format(spw, fntemplate, fnsuffix, filesuffix,
                                  first_endchannel))
    if len(header_fn) != 1:
        raise ValueError("Found too many or too few matches: {0}".format(header_fn))
    else:
        header_fn = header_fn[0]

    # First set up an empty file
    if not os.path.exists(big_filename):
        log.info("Creating large cube based on header {0}".format(header_fn))

        if minimize:
            cube0 = SpectralCube.read(header_fn)
            slices = cube0.subcube_slices_from_mask(cube0.mask,
                                                    spatial_only=True)
            # use the calculated 3rd dimension, plus the difference of the
            # x and y slices
            #header['NAXIS2'] = slices[1].stop-slices[1].start
            #header['NAXIS1'] = slices[2].stop-slices[2].start
            header = cube0[slices].header
        else:
            header = fits.getheader(header_fn)

        # Make an arbitrary, small data before prepping the header
        data = np.zeros((100, 100), dtype=np.float32)
        hdu = fits.PrimaryHDU(data=data, header=header)
        cdelt_sign = np.sign(hdu.header['CDELT3'])
        # Set the appropriate output size (this can be extracted from the LISTOBS)
        naxis3_in = header['NAXIS3']
        header['NAXIS3'] = nchans_total[spwnum]
        header_wcs = wcs.WCS(fits.getheader(header_fn))
        header_specwcs = header_wcs.sub([wcs.WCSSUB_SPECTRAL])
        if cdelt_sign == -1:
            ind0, ind1 = getinds(header_fn)
            #5/20/2017: redoing some of this, and the text below is frightening but no longer relevant
            # a +1 was on the next line before an edit on 4/10/2017
            # it may have been rendered irrelevant when I included +1
            # channel in each cube?  Not clear - the arithmetic no longer
            # makes sense but is empirically necessary.
            assert ind0 == 0
            header['CRPIX3'] = nchans_total[spwnum]
            header['CRVAL3'] = header_specwcs.wcs_pix2world([naxis3_in],1)[0][0]
            assert wcs.WCS(header).sub([wcs.WCSSUB_SPECTRAL]).wcs_pix2world([header['CRPIX3']], 1)[0][0] == header['CRVAL3']

        shape = (header['NAXIS3'], header['NAXIS2'], header['NAXIS1'])



        # Write to disk
        header.tofile(big_filename)
        # Using the 'append' io method, update the *header*
        with open(big_filename, 'rb+') as fobj:
            # Seek past the length of the header, plus the length of the
            # data we want to write.
            # The -1 is to account for the final byte that we are about to
            # write:
            # 'seek' works on bytes, so divide #bits / (bytes/bit)
            fobj.seek(len(header.tostring()) + (shape[0] *
                                                shape[1] *
                                                shape[2] *
                                                int(np.abs(header['BITPIX'])/8)) -
                      1)
            fobj.write(b'\0')

        big_cube = SpectralCube.read(big_filename)
        header_cube = SpectralCube.read(header_fn)
        # in both cases, SpectralCube sorts the extrema
        assert big_cube.spectral_extrema[0] == header_cube.spectral_extrema[0]
        assert np.all(big_cube.wcs.wcs.cdelt == header_cube.wcs.wcs.cdelt)

        log.info("Cube creation completed.  Now moving on to populating it.")


    # Find the appropriate files (this is NOT a good way to do this!  Better to
    # provide a list.  But wildcards are quick & easy...
    files = glob.glob("piece_of_{1}_cube{2}.{0}.chan*{3}".format(spw,fntemplate,fnsuffix,filesuffix))
    log.info("Files to be merged: ")
    log.info(str(files))

    # open the file in update mode (it should have the right dims now)
    hdul = fits.open(big_filename, mode='update')
    main_wcs = wcs.WCS(hdul[0].header)

    if add_beam_info:
        shape = hdul[0].data.shape[0]
        if len(hdul) > 1 and isinstance(hdul[1], fits.BinTableHDU):
            pass
        else:
            hdul.append(fits.BinTableHDU(np.recarray(shape,
                                                     names=['BMAJ','BMIN','BPA','CHAN','POL'],
                                                     formats=['f4','f4','f4','i4','i4'])))

    # sorted so that we deal with zero first, since it has potential to be a problem.
    for fn in ProgressBar(sorted(files)):
        log.info("{0} {1}".format(getinds(fn), fn))
        ind0,ind1 = getinds(fn)

        if 'slices' not in locals():
            if minimize:
                log.info("Determining slices")
                cube0 = SpectralCube.read(header_fn)
                slices = cube0.subcube_slices_from_mask(cube0.mask,
                                                        spatial_only=True)
            else:
                slices = (slice(None),)*3

        cdelt = fits.getheader(fn)['CDELT3']
        if 'cdelt_sign' not in locals():
            cdelt_sign = np.sign(cdelt)
            log.warn("cdelt_sign was not defined: overwriting a"
                     " previously-existing file.  "
                     "This may not be what you want; the data could be going "
                     "opposite the parent cube.  Check that the original "
                     "header is OK. sign(CDELT) is now {0}, "
                     "while for the big header it is {1}"
                     .format(cdelt_sign,
                             np.sign(fits.getheader(big_filename)['CDELT3'])))

        if cropends:
            # don't crop 1st or last pixel in full cube
            if ind0 > 0 or cdelt_sign == -1:
                ind0 = ind0 + cropends
                dataind0 = cropends
                extra = 0
            else:
                dataind0 = 0
                extra = 0 # was an outdated correction; no longer used

            if ((ind1 < nchans_total[spwnum] - 1) and cdelt_sign == 1) or (ind0 > 0 and cdelt_sign == -1):
                ind1 = ind1 - cropends
                dataind1 = - cropends - extra
            elif (ind0 == 0) and cdelt_sign == -1:
                dataind1 = None
            else:
                dataind1 = None
        else:
            dataind0 = 0
            dataind1 = None

        if cdelt_sign == -1:
            log.debug("Reversing indices from {0} {1} to ".format(ind0,ind1))
            ind1, ind0 = (nchans_total[spwnum] - ind0,
                          nchans_total[spwnum] - ind1)
            log.debug("{0} {1}".format(ind0, ind1))
            if ind0 < 0:
                ind0 = 0

        plane = hdul[0].data[ind0]
        if np.all(plane == 0) or overwrite_existing:
            log.info("Replacing indices {0}->{2} {1}"
                     .format(getinds(fn), fn, (ind0,ind1)))

            data = fits.getdata(fn)
            dwcs = wcs.WCS(fits.getheader(fn))

            dwcs0 = dwcs.sub([wcs.WCSSUB_SPECTRAL]).wcs_pix2world([dataind0], 0)[0][0]
            dwcs1 = dwcs.sub([wcs.WCSSUB_SPECTRAL]).wcs_pix2world([data.shape[0]+(dataind1 or -1)], 0)[0][0]
            hwcs0 = main_wcs.sub([wcs.WCSSUB_SPECTRAL]).wcs_pix2world([ind0], 0)[0][0]
            hwcs1 = main_wcs.sub([wcs.WCSSUB_SPECTRAL]).wcs_pix2world([ind1], 0)[0][0]
            
            if dwcs0 != hwcs0:
                log.error("current data, big cube indices: {0},{1} and {2},{3}"
                          .format(dataind0,dataind1,ind0,ind1))
                raise ValueError("World coordinates of first pixels do not match: {0} - {1} = {2} ({3} cdelt)"
                                 .format(dwcs0,hwcs0,dwcs0-hwcs0,(dwcs0-hwcs0)/cdelt))
            if hwcs1 != dwcs1:
                log.error("current data, big cube indices: {0},{1} and {2},{3}"
                          .format(dataind0,dataind1,ind0,ind1))
                raise ValueError("World coordinates of last pixels do not match: {0} - {1} = {2} ({3} cdelt)"
                                 .format(dwcs1,hwcs1,dwcs1-hwcs1,(dwcs1-hwcs1)/cdelt))


            if bmaj_limits is not None:
                beamtable = fits.open(fn)[1]
                ok_beam = ((beamtable.data['BMAJ'] > bmaj_limits[0]) &
                           (beamtable.data['BMAJ'] < bmaj_limits[1]))
                data[~ok_beam] = np.nan

            if not debug_mode:
                if add_beam_info:
                    beamtable = fits.open(fn)[1]
                    hdul[1].data[ind0:ind1] = beamtable.data[dataind0:dataind1]


                hdul[0].data[ind0:ind1,:,:] = data[dataind0:dataind1, slices[1], slices[2]]
                hdul.flush()

    if debug_mode:
        log.setLevel(lvl)

if __name__ == "__main__":
    for robust in (0,2):
        for spw in (0,1,2,3):

            mxind = get_max_ind('piece_of_full_SgrB2_TETC7m_r{1}_cube.spw{0}.channels*fits'.format(spw, robust))
            if mxind < min_nchans:
                log.critical("Skipping {0}:{1} b/c only {2} chans".format(robust, spw, mxind))
                continue
            nchans_total[spw] = mxind
            log.info("nchans_total[{0},{1}] = {2}".format(spw, robust, mxind))

            if os.path.exists('piece_of_full_SgrB2_TETC7m_r{1}_cube.spw{0}.channels0to75.image.pbcor.fits'.format(spw, robust)):

                make_spw_cube(spw='spw{0}', spwnum=spw,
                              fntemplate='full_SgrB2_TETC7m_r{0}'.format(robust),
                              overwrite_existing=False, bmaj_limits=None,
                              fnsuffix="", filesuffix='image.pbcor.fits',
                              first_endchannel=75,
                              debug_mode=True,
                              cropends=1, minimize=True, add_beam_info=True)
            elif os.path.exists('piece_of_full_SgrB2_TETC7m_r{1}_cube.spw{0}.channels0to373.image.pbcor.fits'.format(spw, robust)):

                make_spw_cube(spw='spw{0}', spwnum=spw,
                              fntemplate='full_SgrB2_TETC7m_r{0}'.format(robust),
                              overwrite_existing=False, bmaj_limits=None,
                              fnsuffix="", filesuffix='image.pbcor.fits',
                              cropends=1, minimize=True, add_beam_info=True)
