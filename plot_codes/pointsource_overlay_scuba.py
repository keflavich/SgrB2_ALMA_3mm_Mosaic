import numpy as np
import paths
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import pylab as pl
import matplotlib
from astropy.io import fits
from astropy import wcs
from mpl_plot_templates import asinh_norm
from overlay_common import core_phot_tbl, plotcores

pl.matplotlib.rc_file('pubfiguresrc')

pl.rcParams['figure.figsize'] = (12,8)
pl.rcParams['figure.dpi'] = 75.
pl.rcParams['savefig.dpi'] = 300.
pl.rcParams['axes.labelsize'] = 9
pl.rcParams['xtick.labelsize'] = 8
pl.rcParams['ytick.labelsize'] = 8
tick_fontsize = 6
if matplotlib.__version__[0] == '1':
    markersize = 10
elif matplotlib.__version__[0] == '2':
    markersize = 0.5


hdu = fits.open(paths.Fpath('./other/scuba_Herschel_Feathered.fits'))[0]
mywcs = wcs.WCS(hdu.header).celestial
wcsaxes = mywcs # WCSaxes(mywcs.to_header())

fig3 = pl.figure(3)
fig3.clf()
ax = fig3.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)

ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
ra.set_axislabel("RA (J2000)", fontsize=pl.rcParams['axes.labelsize'])
dec.set_axislabel("Dec (J2000)", fontsize=pl.rcParams['axes.labelsize'], minpad=0)
ra.ticklabels.set_fontsize(tick_fontsize)
ra.set_ticks(exclude_overlapping=True)
dec.ticklabels.set_fontsize(tick_fontsize)
dec.set_ticks(exclude_overlapping=True)

im = ax.imshow(hdu.data.squeeze(),
               transform=ax.get_transform(wcs.WCS(hdu.header).celestial),
               vmin=2300, vmax=55000, cmap=pl.cm.gray_r, origin='lower',
               norm=asinh_norm.AsinhNorm())
tr_fk5 = ax.get_transform("fk5")
(x1,y1),(x2,y2) = (206,174),(348,333)
ax.axis([x1,x2,y1,y2])

coredots = plotcores(ax, alpha=0.5,
                     transform=tr_fk5,
                     markersize=markersize, zorder=50)
cb = pl.colorbar(mappable=im)
cb.set_label("$S_{450 \mu m}$ [MJy sr$^{-1}$]")
fig3.savefig(paths.fpath("cores_on_SCUBA_feathered.png"), bbox_inches='tight')

im = ax.imshow(hdu.data.squeeze(),
               transform=ax.get_transform(wcs.WCS(hdu.header).celestial),
               vmin=2300, vmax=25000, cmap=pl.cm.gray_r, origin='lower',
               norm=asinh_norm.AsinhNorm())
cb.set_clim(im.get_clim())
cb.draw_all()

fig3.savefig(paths.fpath("cores_on_SCUBA_feathered_saturated.png"), bbox_inches='tight')
