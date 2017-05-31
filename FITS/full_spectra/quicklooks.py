import numpy as np
import pylab as pl
import pyspeckit
import glob
from astropy import log

fig = pl.figure(figsize=(12,16))

for fn in glob.glob("*_full_SgrB2_TETC7m_r0_spw0_lines.fits"):

    for spw in (0,1,2,3):
        fn_ = fn.replace("spw0","spw{0}".format(spw))

        ax = pl.subplot(4,1,spw+1)

        sp = pyspeckit.Spectrum(fn_)
        sp.xarr.convert_to_unit('GHz')
        if not any(np.isfinite(sp.data)):
            log.warn("{0} is bad".format(fn_))
        else:
            log.info("{0} is good".format(fn_))
            sp.plotter(figure=fig, axis=ax, clear=True)

        if spw != 0:
            ax.set_title("")
        if spw == 3:
            sp.plotter.savefig('quicklooks/{0}.png'.format(fn.replace("_spw0","")[:-5]))
