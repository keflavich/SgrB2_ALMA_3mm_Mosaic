import numpy as np
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
from collections import defaultdict
from files import contfilename
import warnings
from visualization import make_scalebar
from constants import distance
from overlay_common import core_phot_tbl, plotcores

warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)


if int(matplotlib.__version__[0]) >= 2:
    pl.rcParams['figure.dpi'] = 75.
    pl.rcParams['savefig.dpi'] = 300.
    pl.rcParams['axes.labelsize'] = 8
    pl.rcParams['axes.titlesize'] = 8
    pl.rcParams['xtick.labelsize'] = 6
    pl.rcParams['ytick.labelsize'] = 6
    pl.rcParams['axes.linewidth'] = 0.15
    tick_fontsize = 6
    markersize = 3
else:
    pl.rcParams['savefig.dpi'] = 300.
    tick_fontsize = 10
    markersize = 8


for regionname in ('DeepSouth', 'MandN'):

    vmax_hi = defaultdict(lambda: 0.25*1e3)
    vmax_hi['continuum'] = 0.01*1e3
    vmin_hi = defaultdict(lambda: -0.0001*1e3)
    vmin_hi['continuum'] = -0.0005*1e3
    vmax_lo = defaultdict(lambda: 0.1*1e3)
    vmax_lo['continuum'] = 0.0025*1e3
    vmin_lo = defaultdict(lambda: -0.001*1e3)
    vmin_lo['continuum'] = -0.0002*1e3

    if regionname == 'DeepSouth':
        # Deep South
        bottomleft = coordinates.SkyCoord("17:47:24.199", "-28:26:02.565", unit=(u.h, u.deg), frame='fk5')
        topright = coordinates.SkyCoord("17:47:17.666", "-28:23:30.722", unit=(u.h, u.deg), frame='fk5')
        scalebarpos = coordinates.SkyCoord("17:47:23.7", "-28:23:45.0", unit=(u.h, u.deg), frame='fk5')
    elif regionname == 'MandN':
        bottomleft = coordinates.SkyCoord("17:47:24.199", "-28:23:30.722", unit=(u.h, u.deg), frame='fk5')
        topright = coordinates.SkyCoord("17:47:14.666", "-28:21:04.980", unit=(u.h, u.deg), frame='fk5')
        scalebarpos = coordinates.SkyCoord("17:47:23.7", "-28:23:25.0", unit=(u.h, u.deg), frame='fk5')
    else:
        raise Exception

    for line in ("continuum","HC3N","HCN","HNC","HCOp"):
        if line == 'continuum':
            hdu_line = fits.open(contfilename)[0]
        else:
            hdu_line = fits.open(paths.Fpath('merge/max/SgrB2_b3_7M_12M.{0}.image.pbcor_max_medsub.fits'.format(line)))[0]
        mywcs = wcs.WCS(hdu_line.header).celestial
        #wcsaxes = WCSaxes(mywcs.to_header())
        wcsaxes = mywcs # hack for newer version of wcsaxes that is included in astropy

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
                       transform=ax.get_transform(mywcs),
                       vmin=vmin_hi[line], vmax=vmax_hi[line], cmap=pl.cm.gray_r,
                       origin='lower', norm=asinh_norm.AsinhNorm())
        tr_fk5 = ax.get_transform("fk5")
        #(x1,y1),(x2,y2) = (1200,434),(2142,1743)
        # wrong (x1,y1),(x2,y2) = tr_fk5.transform_point([bottomleft.ra.deg, bottomleft.dec.deg]),tr_fk5.transform_point([topright.ra.deg, topright.dec.deg])
        (x1,y1),(x2,y2) = (mywcs.wcs_world2pix([[bottomleft.ra.deg,
                                                 bottomleft.dec.deg]],0)[0],
                           mywcs.wcs_world2pix([[topright.ra.deg,
                                                 topright.dec.deg]],0)[0]
                          )

        make_scalebar(ax, scalebarpos,
                      length=(0.5*u.pc / distance).to(u.arcsec,
                                                      u.dimensionless_angles()),
                      color='k',
                      label='0.5 pc',
                      text_offset=1.0*u.arcsec,
                     )


        ax.set_aspect(1)
        ax.axis([x1,x2,y1,y2])

        fig3.canvas.draw()
        assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

        cax = fig3.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                             ax.bbox._bbox.y1-ax.bbox._bbox.y0])
        cb = fig3.colorbar(mappable=im, cax=cax)
        cb.set_label("$S_{3 mm}$ [mJy beam$^{-1}$]")

        fig3.savefig(paths.fpath("core_overlays/{0}_peak_{1}.png".format(line,regionname)), bbox_inches='tight')
        fig3.savefig(paths.fpath("core_overlays/{0}_peak_{1}.pdf".format(line,regionname)), bbox_inches='tight')

        coredots = plotcores(ax, alpha=0.5, transform=tr_fk5,
                             markersize=markersize, zorder=50)
        ax.axis([x1,x2,y1,y2])
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}.png".format(line,regionname)), bbox_inches='tight')
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}.pdf".format(line,regionname)), bbox_inches='tight')

        im.set_visible(False)

        im = ax.imshow(hdu_line.data.squeeze()*1e3,
                       transform=ax.get_transform(mywcs),
                       vmin=vmin_lo[line], vmax=vmax_lo[line], cmap=pl.cm.gray_r,
                       origin='lower', norm=asinh_norm.AsinhNorm())
        cb.on_mappable_changed(mappable=im)
        fig3.canvas.draw()

        ax.axis([x1,x2,y1,y2])
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_saturated.png".format(line,regionname)), bbox_inches='tight')
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_saturated.pdf".format(line,regionname)), bbox_inches='tight')

        coredots.set_visible(False)
        fig3.savefig(paths.fpath("core_overlays/{0}_peak_{1}_saturated.png".format(line,regionname)), bbox_inches='tight')
        fig3.savefig(paths.fpath("core_overlays/{0}_peak_{1}_saturated.pdf".format(line,regionname)), bbox_inches='tight')
