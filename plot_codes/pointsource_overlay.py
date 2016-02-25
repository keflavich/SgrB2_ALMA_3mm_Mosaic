import numpy as np
import paths
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import pylab as pl
from astropy.io import fits
from astropy import wcs
from wcsaxes import WCS as WCSaxes
from astropy.convolution import convolve, Gaussian2DKernel

core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
cores = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                             frame='fk5')

hdu = fits.open('/Users/adam/work/sgrb2/continuumdata/SGRB2_1.3CM_fix.fits')[0]
hdu = fits.open('/Users/adam/work/sgrb2/continuumdata/SGRB2_1.3CM_fix_gal.fits')[0]

mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
wcsaxes = WCSaxes(mywcs.to_header())

fig = pl.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)
im = ax.imshow(hdu.data.squeeze()*1e3, cmap=pl.cm.gray_r, origin='lower', vmin=-0.03e3, vmax=0.05e3)
lon = ax.coords['glon']
#lon.set_major_formatter('hh:mm:ss.s')
lat = ax.coords['glat']
lon.set_axislabel("Galactic Longitude")
lat.set_axislabel("Galactic Latitude")
#ra = ax.coords['ra']
#ra.set_major_formatter('hh:mm:ss.s')
#dec = ax.coords['dec']
#ra.set_axislabel("RA (J2000)")
#dec.set_axislabel("Dec (J2000)")

lims = ax.axis()
mylims = ((0.6504856, -0.021613), (0.68863811, -0.065048))
mylims = ((266.81808,-28.362349), (266.86288, -28.40042))
(x1,y1),(x2,y2) = mywcs.wcs_world2pix(mylims, 0)
# fk5 (x1,y1),(x2,y2) = (500,-1500),(2500,2500)
(x1,y1),(x2,y2) = (000,-100),(900,800)


clims = mywcs.wcs_pix2world([[lims[0],lims[2]], [lims[1],lims[3]]], 0)
bins = [np.linspace(clims[1,0], clims[0,0], 200),
        np.linspace(clims[0,1], clims[1,1], 200)]

tr_fk5 = ax.get_transform("fk5")

coredots, = ax.plot(cores.ra, cores.dec, 'r.', transform=tr_fk5, markersize=2,
                    zorder=50)

ax.axis([x1,x2,y1,y2])

fig.savefig(paths.fpath("cores_on_1.3cm_continuum.png"))

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
fig.savefig(paths.fpath("coredensity_on_1.3cm_continuum_withdots.png"))
coredots.set_visible(False)
fig.savefig(paths.fpath("coredensity_on_1.3cm_continuum.png"))

hdu2 = fits.open('/Users/adam/work/sgrb2/continuumdata/sgrb2_20cm_12as.fits')[0]
mywcs = wcs.WCS(hdu2.header).celestial
wcsaxes = WCSaxes(mywcs.to_header())

fig = pl.figure(2)
fig.clf()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)

ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
ra.set_axislabel("RA (J2000)")
dec.set_axislabel("Dec (J2000)")

ax.imshow(hdu2.data.squeeze(), transform=ax.get_transform(wcs.WCS(hdu2.header).celestial),
          vmax=0.45, cmap=pl.cm.gray_r, origin='lower', )
tr_fk5 = ax.get_transform("fk5")
(x1,y1),(x2,y2) = (1807,2100),(2221,2697)
ax.axis([x1,x2,y1,y2])

coredots, = ax.plot(cores.ra, cores.dec, 'r.', transform=tr_fk5, markersize=2,
                    zorder=50, )
fig.savefig(paths.fpath("cores_on_20cm_continuum.png"))
pl.draw()
pl.show()
