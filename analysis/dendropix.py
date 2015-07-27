"""
Try to identify independent velocity peaks along each LOS in a cube
"""

from __future__ import print_function
import pylab as pl
from astropy.io import fits
import numpy as np
import pyspeckit
from spectral_cube import SpectralCube
import os
import glob
from astropy.table import Table,Column
from astropy import units as u
import pyregion
import time
from astropy import stats
import json
import astrodendro
from astropy.utils.console import ProgressBar

dpath = lambda x: os.path.join("/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic/FITS/",x)
tpath = lambda x: os.path.join("/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic/dendropix/",x)

def dendropix(fileprefix='SgrB2_b3_12M.HC3N'):
    cube = SpectralCube.read(dpath('{0}.image.pbcor.contsub.fits'.format(fileprefix))).minimal_subcube()
    noise = cube.spectral_slab(-200*u.km/u.s, -100*u.km/u.s).std(axis=0)
    keep_mask = cube.max(axis=0) > noise

    tblfile = tpath('{0}.dendrotable.ecsv'.format(fileprefix))
    if os.path.exists(tblfile):
        table = Table.read(tblfile, format='ascii.ecsv')
    else:
        table = Table([Column(name='xpix'),
                       Column(name='ypix'),
                       Column(name='zpix'),
                       Column(name='lon'),
                       Column(name='lat'),
                       Column(name='velo'),
                       Column(name='peakval'),])

    xpyp_done = set(zip(table['xpix'], table['ypix']))
    all_keepers = zip(*np.where(keep_mask))
    xpyp = [x for x in all_keepers if x not in xpyp_done]
    print(len(xpyp), len(all_keepers), len(xpyp_done))
    if len(xpyp_done) > 0:
        assert len(xpyp) < len(all_keepers)


    for ii,(ypix,xpix) in enumerate(ProgressBar(xpyp)):

        data = cube[:,ypix,xpix].value

        error = noise[ypix,xpix].value
        # alternative:
        #error = stats.sigma_clipped_stats(data)[2]

        D = astrodendro.Dendrogram.compute(data, min_value=0,
                                           min_delta=2*error, min_npix=7,
                                           is_independent=astrodendro.pruning.min_peak(5*error))
        if not D.leaves:
            table.add_row([xpix,ypix,]+[np.nan]*5)
            del D
            continue

        #peaks = [S.get_peak()[0][0] for S in D]
        #peak_vals = [S.get_peak()[1] for S in D]
        #peaks = [cube.spectral_axis[S.get_peak()[0][0]].to(u.km/u.s).value for S in D]

        for S in D:
            (peak_pix,),peak_val = S.get_peak()
            velo,lat,lon = cube.world[peak_pix, ypix, xpix]
            table.add_row([xpix,ypix,peak_pix,lon,lat,velo,peak_val])

        if ii % 100 == 0:
            table.write(tblfile, format='ascii.ecsv')

        del D
        del S
        del data
