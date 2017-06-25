import numpy as np
import regions
import paths
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import pylab as pl
from astropy.io import fits
from astropy import wcs
#from wcsaxes import WCS as WCSaxes
from astropy.convolution import convolve, Gaussian2DKernel
from mpl_plot_templates import asinh_norm
import matplotlib
from visualization import make_scalebar, hide_scalebar
from constants import distance
from overlay_common import core_phot_tbl, plotcores, cores

pl.matplotlib.rc_file('pubfiguresrc')

pl.rcParams['image.interpolation'] = 'nearest'
pl.rcParams['figure.figsize'] = (12,8)
pl.rcParams['figure.dpi'] = 75.
pl.rcParams['savefig.dpi'] = 300.
pl.rcParams['axes.labelsize'] = 13
pl.rcParams['xtick.labelsize'] = 12
pl.rcParams['ytick.labelsize'] = 12
tick_fontsize = 12
if matplotlib.__version__[0] == '1':
    markersize = 6
elif matplotlib.__version__[0] == '2':
    markersize = 0.5


linenames = {'HC3N': '\mathrm{HC}_3\mathrm{N}',
             'HCN': '\mathrm{HCN}',
             'HNC': '\mathrm{HNC}',
             'HCOp': '\mathrm{HCO}^+',
            }

line = 'HC3N'
mylims_fk5 = ((266.8854477, -28.44946601), (266.7838311, -28.33589021))

bubbles = regions.read_ds9(paths.rpath('bubble_ellipses.reg'))


hdu_line_r05 = fits.open(paths.Fpath('merge/max/SgrB2_b3_7M_12M.{0}.image.pbcor_max_medsub.fits'.format(line)))[0]
hdu_line_lores = fits.open(paths.Fpath('merge/max/SgrB2_b3_7M_12M.{0}.image.pbcor_medsub_max.fits'.format(line)))[0]

for hdu_line, suffix, (vmin,vmax) in [
        (hdu_line_r05, 'r05', (-0.1,150)),
        (hdu_line_lores, 'lores', (-0.5, 750))]:

    mywcs = wcs.WCS(hdu_line.header).celestial
    wcsaxes = mywcs # WCSaxes(mywcs.to_header())

    fig3 = pl.figure(3)
    fig3.clf()
    ax = fig3.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)

    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    ra.set_axislabel("RA (J2000)", fontsize=pl.rcParams['axes.labelsize'])
    dec.set_axislabel("Dec (J2000)", fontsize=pl.rcParams['axes.labelsize'], minpad=0.0)
    ra.ticklabels.set_fontsize(tick_fontsize)
    ra.set_ticks(exclude_overlapping=True)
    dec.ticklabels.set_fontsize(tick_fontsize)
    dec.set_ticks(exclude_overlapping=True)

    im = ax.imshow(hdu_line.data.squeeze()*1e3,
                   transform=ax.get_transform(wcs.WCS(hdu_line.header).celestial),
                   vmin=vmin, vmax=vmax, cmap=pl.cm.gray_r,
                   origin='lower', norm=asinh_norm.AsinhNorm())
    tr_fk5 = ax.get_transform("fk5")
    #(x1,y1),(x2,y2) = (680,350),(2720,3150)
    (x1,y1),(x2,y2) = mywcs.wcs_world2pix(mylims_fk5, 0)
    ax.axis([x1,x2,y1,y2])

    fig3.canvas.draw()
    assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

    cax = fig3.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                        ax.bbox._bbox.y1-ax.bbox._bbox.y0])
    cb = fig3.colorbar(mappable=im, cax=cax)
    cb.set_label("$S_{{{0}}}$ [mJy beam$^{{-1}}$]".format(linenames[line]))

    scalebarpos = coordinates.SkyCoord("17:47:13", "-28:26:15.0",
                                       unit=(u.h, u.deg), frame='fk5')
    sb = make_scalebar(ax, scalebarpos,
                       length=(2.0*u.pc / distance).to(u.arcsec,
                                                       u.dimensionless_angles()),
                       color='k',
                       label='2 pc',
                       text_offset=1.0*u.arcsec,
                      )

    for b in bubbles:
        ell = b.to_pixel(mywcs).as_patch()
        ell.set_facecolor('none')
        ell.set_edgecolor('r')
        ax.add_artist(ell)

    fig3.savefig(paths.fpath("bubbles_on_{0}_peak{1}.png".format(line, suffix)), bbox_inches='tight')
