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


hdu = fits.open('/Users/adam/work/sgrb2/continuumdata/SGRB2_1.3CM_fix.fits')[0]
hdu = fits.open('/Users/adam/work/sgrb2/continuumdata/SGRB2_1.3CM_fix_gal.fits')[0]

mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
#wcsaxes = WCSaxes(mywcs.to_header())
wcsaxes = mywcs

fig = pl.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)
im = ax.imshow(hdu.data.squeeze()*1e3, cmap=pl.cm.gray_r, origin='lower', vmin=-0.03e3, vmax=0.05e3)
lon = ax.coords['glon']
#lon.set_major_formatter('hh:mm:ss.s')
lat = ax.coords['glat']
lon.set_axislabel("Galactic Longitude", fontsize=pl.rcParams['axes.labelsize'])
lat.set_axislabel("Galactic Latitude", fontsize=pl.rcParams['axes.labelsize'])
#ra = ax.coords['ra']
#ra.set_major_formatter('hh:mm:ss.s')
#dec = ax.coords['dec']
#ra.set_axislabel("RA (J2000)")
#dec.set_axislabel("Dec (J2000)")

lims = ax.axis()
mylims_gal = ((0.6504856, -0.021613), (0.68863811, -0.065048))
mylims_gal = ((0.6297901196, -0.09914333816), (0.6856479113, 0.02787655912))
mylims_fk5 = ((266.81808,-28.362349), (266.86288, -28.40042))
mylims_fk5 = ((266.8744477, -28.44946601), (266.7838311, -28.33589021))
(x1,y1),(x2,y2) = mywcs.wcs_world2pix(mylims_gal, 0)
x2,x1 = min([x1,x2]), max([x1,x2])
y1,y2 = min([y1,y2]), max([y1,y2])
# fk5 (x1,y1),(x2,y2) = (500,-1500),(2500,2500)
#(x1,y1),(x2,y2) = (000,-100),(900,800)


clims = mywcs.wcs_pix2world([[lims[0],lims[2]], [lims[1],lims[3]]], 0)
bins = [np.linspace(clims[1,0], clims[0,0], 200),
        np.linspace(clims[0,1], clims[1,1], 200)]

tr_fk5 = ax.get_transform("fk5")

#coredots, = ax.plot(cores.ra, cores.dec, 'r.', transform=tr_fk5, markersize=markersize,
#                    zorder=50)
coredots = plotcores(ax, transform=tr_fk5, markersize=markersize, zorder=50)

ax.axis([x1,x2,y1,y2])

lon = ax.coords['glon']
lat = ax.coords['glat']
lon.set_major_formatter('d.ddd')
lat.set_major_formatter('d.ddd')
lon.ticklabels.set_fontsize(tick_fontsize)
lon.set_ticks(exclude_overlapping=True)
lat.ticklabels.set_fontsize(tick_fontsize)
lat.set_ticks(exclude_overlapping=True)


fig.savefig(paths.fpath("cores_on_1.3cm_continuum.png"), bbox_inches='tight')

H,bx,by = np.histogram2d(cores.ra, cores.dec, bins=bins)
H2 = convolve(H, Gaussian2DKernel(2))
cx = (bx[1:]+bx[:-1])/2.
cy = (by[1:]+by[:-1])/2.
con = ax.contour(cx,
                 cy,
                 H2.T, transform=tr_fk5,
                 levels=[0.12,  0.18,  0.24,  0.30,  0.36,  0.42],
                 colors=['k']*10,
                 linewidth=0.5,
                 alpha=0.5,
                 interpolation='bicubic',
                )
#ax.axis(lims)
ax.axis([x1,x2,y1,y2])


fig.canvas.draw()
assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

cax = fig.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                    ax.bbox._bbox.y1-ax.bbox._bbox.y0])
cb = fig.colorbar(mappable=im, cax=cax)
cb.set_label("$S_{1.3 cm}$ [mJy beam$^{-1}$]")

fig.savefig(paths.fpath("coredensity_on_1.3cm_continuum_withdots.png"), bbox_inches='tight')
#coredots.set_visible(False)
for cd in coredots:
    cd.set_visible(False)
fig.savefig(paths.fpath("coredensity_on_1.3cm_continuum.png"), bbox_inches='tight')

hdu2 = fits.open('/Users/adam/work/sgrb2/continuumdata/sgrb2_20cm_12as.fits')[0]
mywcs = wcs.WCS(hdu2.header).celestial
wcsaxes = mywcs # WCSaxes(mywcs.to_header())

fig = pl.figure(2)
fig.clf()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)

ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
ra.set_axislabel("RA (J2000)", fontsize=pl.rcParams['axes.labelsize'],)
dec.set_axislabel("Dec (J2000)", fontsize=pl.rcParams['axes.labelsize'], minpad=0)
ra.ticklabels.set_fontsize(tick_fontsize)
ra.set_ticks(exclude_overlapping=True)
dec.ticklabels.set_fontsize(tick_fontsize)
dec.set_ticks(exclude_overlapping=True)

im = ax.imshow(hdu2.data.squeeze()*1e3,
               transform=ax.get_transform(wcs.WCS(hdu2.header).celestial),
               vmax=0.45*1e3, cmap=pl.cm.gray_r, origin='lower', )
tr_fk5 = ax.get_transform("fk5")
#(x1,y1),(x2,y2) = (1807,2100),(2221,2697)
(x1,y1),(x2,y2) = mywcs.wcs_world2pix(mylims_fk5, 0)
ax.axis([x1,x2,y1,y2])

scalebarpos = coordinates.SkyCoord("17:47:27", "-28:26:15.0",
                                   unit=(u.h, u.deg), frame='fk5')
make_scalebar(ax, scalebarpos,
              length=(2.0*u.pc / distance).to(u.arcsec,
                                              u.dimensionless_angles()),
              color='k',
              label='2 pc',
              text_offset=1.0*u.arcsec,
             )

#coredots, = ax.plot(cores.ra, cores.dec, 'r.', transform=tr_fk5, markersize=markersize,
#                    alpha=0.5,
#                    zorder=50, )
coredots = plotcores(ax, transform=tr_fk5, markersize=markersize, zorder=50,
                     alpha=0.5)

fig.canvas.draw()
assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

cax = fig.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                    ax.bbox._bbox.y1-ax.bbox._bbox.y0])
cb = fig.colorbar(mappable=im, cax=cax)
cb.set_label("$S_{20 \mathrm{cm}}$ [mJy beam$^{-1}$]")

fig.savefig(paths.fpath("cores_on_20cm_continuum.png"), bbox_inches='tight')


hdu_h41a = fits.open(paths.Fpath('merge/max/SgrB2_b3_7M_12M.H41a.image.pbcor_max_medsub.fits'))[0]
mywcs = wcs.WCS(hdu_h41a.header).celestial
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

ax.imshow(hdu_h41a.data.squeeze(), transform=ax.get_transform(wcs.WCS(hdu_h41a.header).celestial),
          vmin=-0.0001, vmax=0.25, cmap=pl.cm.gray_r, origin='lower', norm=asinh_norm.AsinhNorm())
tr_fk5 = ax.get_transform("fk5")
#(x1,y1),(x2,y2) = (680,350),(2720,3150)
(x1,y1),(x2,y2) = mywcs.wcs_world2pix(mylims_fk5, 0)
ax.axis([x1,x2,y1,y2])

#coredots, = ax.plot(cores.ra, cores.dec, 'r.', transform=tr_fk5,
# markersize=markersize, alpha=0.5,
#                    zorder=50, )
coredots = plotcores(ax, transform=tr_fk5, markersize=markersize, zorder=50,
                     alpha=0.5)
fig3.savefig(paths.fpath("cores_on_h41a_peak.png"), bbox_inches='tight')

ax.imshow(hdu_h41a.data.squeeze(), transform=ax.get_transform(wcs.WCS(hdu_h41a.header).celestial),
          vmin=-0.0001, vmax=0.05, cmap=pl.cm.gray_r, origin='lower', norm=asinh_norm.AsinhNorm())
fig3.savefig(paths.fpath("cores_on_h41a_peak_saturated.png"), bbox_inches='tight')




linenames = {'HC3N': '\mathrm{HC}_3\mathrm{N}',
             'HCN': '\mathrm{HCN}',
             'HNC': '\mathrm{HNC}',
             'HCOp': '\mathrm{HCO}^+',
            }

for line in ("HC3N","HCN","HNC","HCOp"):
    hdu_line = fits.open(paths.Fpath('merge/max/SgrB2_b3_7M_12M.{0}.image.pbcor_max_medsub.fits'.format(line)))[0]
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
                   vmin=-0.0001*1e3, vmax=0.25*1e3, cmap=pl.cm.gray_r,
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

    scalebarpos = coordinates.SkyCoord("17:47:27", "-28:26:15.0",
                                       unit=(u.h, u.deg), frame='fk5')
    sb = make_scalebar(ax, scalebarpos,
                       length=(2.0*u.pc / distance).to(u.arcsec,
                                                       u.dimensionless_angles()),
                       color='k',
                       label='2 pc',
                       text_offset=1.0*u.arcsec,
                      )

    if matplotlib.__version__[0] == '1':
        markersize = 6
    elif matplotlib.__version__[0] == '2':
        markersize = 0.5
    #coredots, = ax.plot(cores.ra, cores.dec, 'r.', transform=tr_fk5,
    #                    markersize=markersize, alpha=0.5, zorder=50, )
    fig3.savefig(paths.fpath("{0}_peak.png".format(line)), bbox_inches='tight')
    coredots = plotcores(ax, transform=tr_fk5, markersize=markersize, zorder=50,
                         alpha=0.5)
    fig3.savefig(paths.fpath("cores_on_{0}_peak.png".format(line)), bbox_inches='tight')

    im = ax.imshow(hdu_line.data.squeeze()*1e3,
                   transform=ax.get_transform(wcs.WCS(hdu_line.header).celestial),
                   vmin=-0.001*1e3, vmax=0.1*1e3, cmap=pl.cm.gray_r,
                   origin='lower', norm=asinh_norm.AsinhNorm())
    fig3.savefig(paths.fpath("cores_on_{0}_peak_saturated.png".format(line)), bbox_inches='tight')
    cb.on_mappable_changed(mappable=im)


    # Deep South
    if matplotlib.__version__[0] == '1':
        markersize = 12
    elif matplotlib.__version__[0] == '2':
        markersize = 2
    bottomleft = coordinates.SkyCoord("17:47:24.199", "-28:26:02.565", unit=(u.h, u.deg), frame='fk5')
    topright = coordinates.SkyCoord("17:47:17.666", "-28:23:30.722", unit=(u.h, u.deg), frame='fk5')
    im = ax.imshow(hdu_line.data.squeeze()*1e3,
                   transform=ax.get_transform(wcs.WCS(hdu_line.header).celestial),
                   vmin=-0.0001*1e3, vmax=0.25*1e3, cmap=pl.cm.gray_r, origin='lower',
                   norm=asinh_norm.AsinhNorm())
    tr_fk5 = ax.get_transform("fk5")
    (x1,y1),(x2,y2) = (1200,434),(2142,1743)
    # wrong (x1,y1),(x2,y2) = tr_fk5.transform_point([bottomleft.ra.deg, bottomleft.dec.deg]),tr_fk5.transform_point([topright.ra.deg, topright.dec.deg])
    (x1,y1),(x2,y2) = mywcs.wcs_world2pix([[bottomleft.ra.deg, bottomleft.dec.deg]],0)[0],mywcs.wcs_world2pix([[topright.ra.deg, topright.dec.deg]],0)[0]
    ax.set_aspect(1)
    ax.axis([x1,x2,y1,y2])

    cb.on_mappable_changed(mappable=im)

    hide_scalebar(sb)
    scalebarpos = coordinates.SkyCoord("17:47:23.7", "-28:23:45.0",
                                       unit=(u.h, u.deg), frame='fk5')
    sb = make_scalebar(ax, scalebarpos,
                       length=(2.0*u.pc / distance).to(u.arcsec,
                                                       u.dimensionless_angles()),
                       color='k',
                       label='2 pc',
                       text_offset=1.0*u.arcsec,
                      )


    coredots = plotcores(ax, transform=tr_fk5, markersize=markersize, zorder=50,
                         alpha=0.5)
    #coredots, = ax.plot(cores.ra, cores.dec, 'r.', transform=tr_fk5,
    #                    markersize=markersize, alpha=0.5, zorder=50, )
    ax.axis([x1,x2,y1,y2])
    fig3.savefig(paths.fpath("cores_on_{0}_peak_DeepSouth.png".format(line)), bbox_inches='tight')
    fig3.savefig(paths.fpath("cores_on_{0}_peak_DeepSouth.pdf".format(line)), bbox_inches='tight')

    im = ax.imshow(hdu_line.data.squeeze()*1e3,
                   transform=ax.get_transform(wcs.WCS(hdu_line.header).celestial),
                   vmin=-0.001*1e3, vmax=0.1*1e3, cmap=pl.cm.gray_r,
                   origin='lower', norm=asinh_norm.AsinhNorm())
    cb.on_mappable_changed(mappable=im)
    ax.axis([x1,x2,y1,y2])
    fig3.savefig(paths.fpath("cores_on_{0}_peak_DeepSouth_saturated.png".format(line)), bbox_inches='tight')
    fig3.savefig(paths.fpath("cores_on_{0}_peak_DeepSouth_saturated.pdf".format(line)), bbox_inches='tight')
