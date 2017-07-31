import re
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


scalebarpos = coordinates.SkyCoord("17:47:32.", "-28:26:33.0", unit=(u.h, u.deg), frame='fk5')

contfnr2 = paths.Fpath('merge/continuum/SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_taper1.5as_r2_mask5mJy.image.tt0.pbcor.fits')
contfnr0 = paths.Fpath('merge/continuum/SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_deeper_mask1.5mJy.image.tt0.pbcor.fits')
for contfn,name,(vmin,vmax) in [(contfnr2, "taper", (-1,50)), (contfnr0, "hires",(-1,5))]:
    hdu = fits.open(contfn)[0]

    mywcs = wcs.WCS(hdu.header).celestial
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

    # manual normalization...
    data = hdu.data.squeeze()*1e3
    #data[data<1.5] = 0
    #data[data>4000] = 4000
    #ldata = np.log10(data)

    im = ax.imshow(data,
                   transform=ax.get_transform(mywcs),
                   cmap=pl.cm.gray_r,
                   origin='lower',
                   interpolation='nearest',
                   vmin=vmin, vmax=vmax, norm=asinh_norm.AsinhNorm(),
                   #matplotlib.colors.LogNorm(),
                  )
    tr_fk5 = ax.get_transform("fk5")
    #(x1,y1),(x2,y2) = (1200,434),(2142,1743)
    # wrong (x1,y1),(x2,y2) = tr_fk5.transform_point([bottomleft.ra.deg, bottomleft.dec.deg]),tr_fk5.transform_point([topright.ra.deg, topright.dec.deg])
    bottomleft = coordinates.SkyCoord("17:47:36.2", "-28:27:00.0", unit=(u.h, u.deg), frame='fk5')
    topright = coordinates.SkyCoord("17:47:03.1", "-28:20:00.0", unit=(u.h, u.deg), frame='fk5')
    (x1,y1),(x2,y2) = (mywcs.wcs_world2pix([[bottomleft.ra.deg,
                                             bottomleft.dec.deg]],0)[0],
                       mywcs.wcs_world2pix([[topright.ra.deg,
                                             topright.dec.deg]],0)[0]
                      )

    make_scalebar(ax, scalebarpos,
                  length=(2.5*u.pc / distance).to(u.arcsec,
                                                  u.dimensionless_angles()),
                  color='k',
                  label='2.5 pc',
                  text_offset=5.0*u.arcsec,
                 )

    ax.set_aspect(1)
    ax.axis([x1,x2,y1,y2])
    cb = pl.colorbar(mappable=im)
    cb.set_label("$S_{3 \mathrm{mm}}$ [mJy beam$^{-1}$]")

    fig3.canvas.draw()

    fig3.savefig(paths.fpath("overview_figure_{0}.pdf".format(name)), bbox_inches='tight', dpi=600)

    with open(paths.rpath('labels.reg'),'r') as fh:
        for line in fh.readlines():
            if 'text' in line:
                text = re.compile("text={(.*)}").search(line).groups()[0]
                coordinate_str = re.compile("text\(([0-9,:\.-]*)\)").search(line).groups()[0]
                coord = coordinates.SkyCoord(*coordinate_str.split(","),
                                             unit=(u.hour, u.deg), frame='fk5')
                ax.text(coord.ra.deg, coord.dec.deg, s=text, transform=tr_fk5,
                        ha='left',
                        color='k')


    ax.axis([x1,x2,y1,y2])

    fig3.canvas.draw()

    fig3.savefig(paths.fpath("overview_figure_{0}_labeled.pdf".format(name)), bbox_inches='tight', dpi=600)
