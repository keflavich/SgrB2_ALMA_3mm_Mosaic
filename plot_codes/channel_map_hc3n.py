import numpy as np
from astropy import units as u
import os
import pylab as pl
import aplpy
import paths
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import wcs
import paths
from spectral_cube import SpectralCube
from astropy.visualization import SqrtStretch,AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize

pl.style.use('classic')

annotation_fontsize = 10

hc3nfn = paths.Fpath('merge/lines/HC3N_TP_7m_12m_feather.fits')
#hc3nfn = paths.Fpath('merge/lines/HC3N_TP_7m_12m_feather_r05.fits')

sourcename = 'SgrB2'
species = 'HC3N'
cubefn = hc3nfn

dx = 5
slabs = [(x, x+dx) for x in range(5,105,5)]*u.km/u.s

cube = SpectralCube.read(cubefn)
cube._unit = u.K
scube = cube
scube.allow_huge_operations=True
#cutout = Cutout2D(cube[0,:,:], source, radius, wcs=cube.wcs.celestial)

#scube = cube_cutout = cube[(slice(None),)+cutout.slices_original]
ghzaxis = scube.with_spectral_unit(u.GHz).spectral_axis
if hasattr(scube, 'beam'):
    scube = scube.to(u.K,
                     scube.beam.jtok_equiv(ghzaxis))
else:
    beam = scube.average_beam(0.1)
    scube = scube.to(u.K,
                     beam.jtok_equiv(ghzaxis))


# Begin channel map code here
Nrows = 4
Ncols = 5
figsize = ((Ncols/Nrows)*12,12)
fig3 = pl.figure(3)
if any(fig3.get_size_inches() != figsize):
    pl.close(3)
pl.figure(3, figsize=figsize).clf()
fig, axes = pl.subplots(Nrows, Ncols,
                        sharex=True,
                        sharey=True, num=3)

layers = [scube.spectral_slab(v0,v1).moment0()
          if len(scube.spectral_slab(v0,v1)) > 1
          else np.zeros(scube.shape[1:])*u.K
          for v0,v1 in slabs
         ]
# Determine the maximum value to display
mx = np.max([np.nanmax(x).value for x in layers])
mn = 0
mn = np.max([mn, np.min([np.nanmin(x).value for x in layers])])


for ii,(v1,v2) in enumerate(slabs):
    #pl.subplot(4, 4, ii+1)
    layer = layers[ii]
    ax = axes[int(ii / Ncols), int(ii % Ncols)]
    im = ax.imshow(layer.value, norm=ImageNormalize(vmin=mn, vmax=mx,
                                                    stretch=AsinhStretch(),),
                   cmap=pl.cm.gray_r, interpolation='nearest', origin='lower')
    axlims = ax.axis()
    #ax.plot(refvec_pix[0], refvec_pix[1], 'r--', alpha=0.5)
    ax.annotate("${0:d} < v < {1:d}$".format(int(v1.value),
                                             int(v2.value)), (0.1,
                                                              0.8),
                xycoords='axes fraction', color='k', fontsize=annotation_fontsize)
    ax.axis(axlims)

#fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.92, 0.05, 0.04, 0.9])
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_label("K km s$^{-1}$")

pl.subplots_adjust(hspace=0,
                   wspace=0)

for i in range(Nrows):
    for j in range(Ncols):
        if i == 0:
            axes[i,j].xaxis.set_ticks_position('top')
            pl.setp(axes[i,j].get_xticklabels(), visible=False)
            axes[i,j].xaxis.set_ticklabels([])
        elif i == Nrows-1:
            axes[i,j].xaxis.set_ticks_position('bottom')
            pl.setp(axes[i,j].get_xticklabels(), visible=True)
        else:
            axes[i,j].xaxis.set_ticks_position('none')
            pl.setp(axes[i,j].get_xticklabels(), visible=False)
            axes[i,j].xaxis.set_ticklabels([])

        if j == 0:
            axes[i,j].yaxis.set_ticks_position('left')
        elif j == Ncols-1:
            axes[i,j].yaxis.set_ticks_position('right')
            pl.setp(axes[i,j].get_yticklabels(), visible=False)
            axes[i,j].yaxis.set_ticklabels([])
        else:
            axes[i,j].yaxis.set_ticks_position('none')
            pl.setp(axes[i,j].get_yticklabels(), visible=False)
            axes[i,j].yaxis.set_ticklabels([])

pl.subplots_adjust(hspace=0,
                   wspace=0)
pl.savefig(paths.fpath('channelmaps/{0}_{1}_channelmaps.png'.format(sourcename,species)),
           bbox_inches='tight', dpi=150)
pl.savefig(paths.fpath('channelmaps/{0}_{1}_channelmaps.pdf'.format(sourcename,species)),
           bbox_inches='tight')
