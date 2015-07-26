"""
http://docs.astropy.org/en/stable/io/fits/appendix/faq.html#how-can-i-create-a-very-large-fits-file-from-scratch
"""
from astropy.io import fits
import numpy as np
import glob
import re
import os
from astropy.utils.console import ProgressBar

def make_spw_cube(spw='spw0'):
    if not os.path.exists('SgrB2_12m_{0}_lines.fits'.format(spw)):
        header = fits.getheader('calibrated.ms.line.{0}.channels7296to7680.image.fits'.format(spw))
        data = np.zeros((100, 100), dtype=np.float32)
        hdu = fits.PrimaryHDU(data=data, header=header)
        header['NAXIS3'] = 7680
        header.tofile('SgrB2_12m_{0}_lines.fits'.format(spw))
        with open('SgrB2_12m_{0}_lines.fits'.format(spw), 'rb+') as fobj:
             # Seek past the length of the header, plus the length of the
             # data we want to write.
             # The -1 is to account for the final byte taht we are about to
             # write:
             fobj.seek(len(header.tostring()) + (header['NAXIS1'] * header['NAXIS2'] * header['NAXIS3'] * np.abs(header['BITPIX'])/8) - 1)
             fobj.write('\0')

    files = glob.glob("*{0}.chan*fits".format(spw))
    print(files)

    def getinds(fn):
        inds = re.search('channels([0-9]*)to([0-9]*)', fn).groups()
        return [int(ii) for ii in inds]

    hdul = fits.open('SgrB2_12m_{0}_lines.fits'.format(spw), mode='update')
    for fn in ProgressBar(files):
        print(getinds(fn), fn)
        ind0,ind1 = getinds(fn)
        hdul[0].data[ind0:ind1,:,:] = fits.getdata(fn)
        hdul.flush()
