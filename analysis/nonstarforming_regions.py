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
import reproject
import pyregion
import radio_beam

from constants import distance

import paths

datapath = '/Users/adam/work/sgrb2/alma/FITS/continuumdata'
colfilename = datapath+'/column_maps/scuba_col_herscheltem.fits'

colfile = fits.open(colfilename)[0]

sgrb2contfile = fits.open(paths.Fpath('merge/continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy.image.tt0.pbcor.fits'))

tbl = table.Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac',)

# Set the beam to be approximately the measured beam size
beam = radio_beam.Beam(11*u.arcsec)
beam_rad = beam.major.to(u.deg).value

observed_region,_ = reproject.reproject_interp(sgrb2contfile, colfile.header)

# make the WCS grid
yy,xx = np.indices(colfile.data.shape)
inds = np.transpose((xx.ravel(), yy.ravel()))
ra,dec = wcs.WCS(colfile.header).celestial.wcs_pix2world(inds, 0).T
ra = ra.reshape(xx.shape)
dec = dec.reshape(yy.shape)

mask = np.isfinite(observed_region)

for row in ProgressBar(tbl):
    source_dist = ((ra-row['RA'])**2 + (dec-row['Dec'])**2)**0.5
    mask[source_dist < beam_rad] = False

lo = 5e22
hi = 3e25
bins = np.logspace(np.log10(lo)-0.5, np.log10(hi), 100)

import pylab as pl
observed_mask = np.isfinite(observed_region)
pl.figure(1).clf()
H,L,P = pl.hist(colfile.data[observed_mask], bins=bins, log=True,
                alpha=0.5, color='k',
                linewidth=2,
                #normed=True,
                histtype='step')

H2,L2,P2 = pl.hist(colfile.data[mask], bins=bins, log=True,
                   alpha=0.5, color='b', linewidth=2,
                   zorder=10,
                   #normed=True,
                   histtype='step')

H3,L3,P3 = pl.hist(colfile.data[observed_mask & ~mask], bins=bins, log=True,
                   alpha=0.8, color='r',
                   zorder=5,
                   linestyle='--',
                   #normed=True,
                   histtype='step')
pl.xscale('log')

pl.xlim(np.min([L.min(), L.min()]), L.max())
pl.xlim(lo,hi)
pl.ylim(0.5,np.max([H.max(), H.max()])*1.1)
pl.semilogx()
pl.yscale('log', nonposy='clip')
pl.xlabel("Column Density N(H$_2$) [cm$^{-2}$]")
pl.ylabel("Number of pixels")

ax1 = pl.gca()
ax2 = ax1.twiny()
print("ax1 xlims: {0}".format(ax1.get_xlim()))
pl.draw()
ax2.set_xlim(np.array(ax1.get_xlim())*(2.8*u.Da).to(u.g).value)
print("ax2 xlims: {0}".format(ax2.get_xlim()))
ax2.set_xscale('log')
ax2.set_xlabel("Column Density [g cm$^{-2}$]")
pl.draw()

pl.savefig(paths.fpath("column_density_distribution_with_and_without_SF.png"),
           bbox_inches='tight')
