import os
import copy

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import table
from astropy import coordinates
from astropy import wcs
from astropy.convolution import convolve_fft,Gaussian2DKernel
from astropy.utils.console import ProgressBar
import astropy.visualization
import reproject
import pyregion
import radio_beam
import pylab as pl
import uvcombine
from scipy import ndimage
from visualization import hide_labels_nonwcs

from constants import distance

import paths

datapath = '/Users/adam/work/sgrb2/continuumdata/scuba_herschel_merge'

scubafile = fits.open(os.path.join(datapath,
                                   'scuba_shifted_MJySr.fits'))[0]

sgrb2contfile = fits.open(paths.Fpath('merge/continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy.image.tt0.pbcor.fits'))

tbl = table.Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac',)
coords = coordinates.SkyCoord(tbl['RA'], tbl['Dec'], frame='fk5', unit=(u.deg, u.deg))
mywcs = wcs.WCS(scubafile.header)
xpix,ypix = mywcs.wcs_world2pix(np.transpose([coords.ra.deg, coords.dec.deg]), 0).T
xpix = np.round(xpix).astype('int')
ypix = np.round(ypix).astype('int')

observed_region,_ = reproject.reproject_interp(sgrb2contfile, scubafile.header)



lo = 1e3
hi = 3e5
bins = np.logspace(np.log10(lo), np.log10(hi), 100)
mx = 0

observed_mask = np.isfinite(observed_region)
fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

fig2 = pl.figure(2, figsize=(20,6))
fig2.clf()


slices = ndimage.find_objects(np.broadcast_arrays(observed_mask,
                                                  scubafile.data)[0])[0]

for ii,scalefactor in enumerate([0.5,1,3,5,8,10]):

    color = 'k' if scalefactor == 3 else None

    result2 = uvcombine.uvcombine.feather_simple(scubafile,
                                                 os.path.join(datapath, 'herschel_500um_upsampled.fits'),
                                                 lowresfwhm=45*u.arcsec,
                                                 # scale factor b/c we think JCMT
                                                 # is miscalibrated?
                                                 highresscalefactor=scalefactor,
                                                )


    H,L,P = ax1.hist(result2[observed_mask], bins=bins, log=True, alpha=0.75,
                     linewidth=1, label=str(scalefactor),
                     color=color,
                     histtype='step')
    mx = max([mx, H.max()])

    ax2 = fig2.add_subplot(2,6,ii+1)
    ax2.imshow(result2.real[slices], cmap='gray', vmin=lo, vmax=hi,
               norm=astropy.visualization.simple_norm(result2.real[slices],
                                                      stretch='asinh',
                                                      min_cut=lo,
                                                      max_cut=hi,
                                                      asinh_a=0.001),
               origin='lower',
              )
    hide_labels_nonwcs(ax2)
    ax2.set_title(scalefactor)

    ax3 = fig2.add_subplot(2,6,ii+1+6)
    ax3.hist(result2.real[observed_mask], bins=bins, log=True, alpha=0.75,
             linewidth=2, label=str(scalefactor), histtype='step',
             color=P[0].get_edgecolor())

    ax3.set_xscale('log')
    ax3.set_xlim(1e3,3e5)
    ax3.set_ylim(2, mx*1.1)
    ax3.set_yscale('log', nonposy='clip')

    ax4 = ax3.twinx()
    ax4.plot(sorted(result2.real[ypix,xpix]),
             np.arange(len(xpix))/len(xpix), linestyle='--',
             color=P[0].get_edgecolor(), alpha=0.5, linewidth=1)
    ax3.set_xlim(1e3,3e5)

    if ii == 0:
        ax3.set_ylabel("Number of pixels", fontsize=16)
    else:
        for tt in ax3.yaxis.get_ticklabels():
            tt.set_visible(False)
    for tt in ax4.yaxis.get_ticklabels():
        tt.set_visible(False)


    if scalefactor == 3:
        ax3.set_xlabel("$S_{450 \mu m}$ [MJy sr$^{-1}$]", fontsize=16)
        for jj in range(6):
            ax3 = fig2.add_subplot(2,6,jj+1+6)
            ax3.hist(result2.real[observed_mask], bins=bins, log=True, alpha=0.75,
                     linewidth=2, histtype='step',
                     color='k', zorder=-5)
            ax4 = ax3.twinx()
            ax4.plot(sorted(result2.real[ypix,xpix]),
                     np.arange(len(xpix))/len(xpix), linestyle='--',
                     color='k', alpha=0.5, linewidth=1, zorder=-5)
            for tt in ax4.yaxis.get_ticklabels():
                tt.set_visible(False)
            ax3.set_xlim(1e3,3e5)
            ax4.set_xlim(1e3,3e5)



ax1.set_xscale('log')
ax1.set_xlim(np.min([L.min(), L.min()]), L.max())
ax1.set_xlim(lo,hi)
ax1.set_ylim(2, mx*1.1)
ax1.set_yscale('log', nonposy='clip')
ax1.set_xlabel("$S_{450 \mu m}$ [MJy sr$^{-1}$]", fontsize=16)
ax1.set_ylabel("Number of pixels", fontsize=16)

pl.figure(1)
pl.legend(loc='best')

fig2.subplots_adjust(wspace=0)

fig1.savefig(paths.fpath("scuba_feather_scalefactor_comparison_histograms.png"),
             bbox_inches='tight')
fig2.savefig(paths.fpath("scuba_feather_scalefactor_comparison.pdf"),
             bbox_inches='tight')
