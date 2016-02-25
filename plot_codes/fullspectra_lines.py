import numpy as np
from astropy import units as u
from astropy import constants
from astropy.table import Table
import paths
import pyspeckit
import glob
import matplotlib
matplotlib.use('Qt4Agg')
import os
import pylab as pl

tbl = Table.read(paths.mpath('diffuse_lines_mainonly.ipac'),
                 format='ascii.ipac')
rfreq = tbl['FreqGHz']
rfreq[rfreq.mask] = tbl['MeasFreqGHz'][rfreq.mask]

freq = rfreq * (1-65*u.km/u.s/constants.c)

for fn in glob.glob(paths.spath("*.fits")):
    sp = pyspeckit.Spectrum(fn)
    sp.xarr.convert_to_unit('GHz')
    sp.plotter(figure=pl.figure(1), clear=True)
    sp.plotter.line_ids(tbl['Species'], freq*u.GHz, unit='GHz')
    outname = os.path.splitext(os.path.split(fn)[1])[0]
    sp.plotter.savefig(paths.sppath("{0}.png".format(outname)))


"""
import paths
import pyspeckit
import numpy as np
import glob
sp = pyspeckit.Spectrum('ionization_front_circle_SgrB2_12m_spw0_lines.fits')
sp.plotter()
tbl = Table.read('/Users/adam/work/sgrb2/molecules/diffuse_lines_mainonly.ipac',format='ascii.ipac')
freq = tbl['FreqGHz']
freq[freq.mask] = tbl['MeasFreqGHz'][freq.mask]
sp.plotter.line_ids(tbl['Species'], freq*u.GHz, unit='GHz')
"""
