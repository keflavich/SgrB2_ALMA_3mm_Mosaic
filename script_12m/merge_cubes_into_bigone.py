"""
http://docs.astropy.org/en/stable/io/fits/appendix/faq.html#how-can-i-create-a-very-large-fits-file-from-scratch
"""
from astropy import log
from astropy.io import fits
import numpy as np
import glob
import re
import os
from astropy.utils.console import ProgressBar

def make_spw_cube(spw='spw0'):

    # First set up an empty file
    if not os.path.exists('SgrB2_12m_{0}_lines.fits'.format(spw)):
        header = fits.getheader('calibrated.ms.line.{0}.channels0to384.image.fits'.format(spw))
        # Make an arbitrary, small data before prepping the header
        data = np.zeros((100, 100), dtype=np.float32)
        hdu = fits.PrimaryHDU(data=data, header=header)
        # Set the appropriate output size (this can be extracted from the LISTOBS)
        header['NAXIS3'] = 7680
        # Write to disk
        header.tofile('SgrB2_12m_{0}_lines.fits'.format(spw))
        # Using the 'append' io method, update the *header*
        with open('SgrB2_12m_{0}_lines.fits'.format(spw), 'rb+') as fobj:
             # Seek past the length of the header, plus the length of the
             # data we want to write.
             # The -1 is to account for the final byte that we are about to
             # write:
             # 'seek' works on bytes, so divide #bits / (bytes/bit)
             fobj.seek(len(header.tostring()) + (header['NAXIS1'] *
                                                 header['NAXIS2'] *
                                                 header['NAXIS3'] *
                                                 np.abs(header['BITPIX'])/8) -
                       1)
             fobj.write('\0')

    # Find the appropriate files (this is NOT a good way to do this!  Better to
    # provide a list.  But wildcards are quick & easy...
    files = glob.glob("*{0}.chan*fits".format(spw))
    log.info(str(files))

    # Extract the appropriate pixel indices from the file name.
    # A more sophisticated approach is probably better, in which the individual
    # cubes are inspected for their start/end frequencies.
    # But, on the other hand, for this process to make any sense at all, you
    # have to have done the original cube imaging right
    def getinds(fn):
        inds = re.search('channels([0-9]*)to([0-9]*)', fn).groups()
        return [int(ii) for ii in inds]

    # open the file in update mode (it should have the right dims now)
    hdul = fits.open('SgrB2_12m_{0}_lines.fits'.format(spw), mode='update')
    for fn in ProgressBar(files):
        log.info("{0} {1}".format(getinds(fn), fn))
        ind0,ind1 = getinds(fn)
        plane = hdul[0].data[ind0]
        if np.all(plane == 0):
            log.info("Replacing indices {0} {1}".format(getinds(fn), fn))
            hdul[0].data[ind0:ind1,:,:] = fits.getdata(fn)
            hdul.flush()
