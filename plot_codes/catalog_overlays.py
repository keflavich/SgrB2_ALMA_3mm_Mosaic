from astroquery.vizier import Vizier
from astropy import coordinates
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.table import Table
import matplotlib
import pylab as pl
import paths
from files import contfilename as fncont

tick_fontsize = 6

Chambers2014 = Vizier.get_catalogs('J/A+A/563/A68')
Chambers2014_H2O = Chambers2014[0]
h2o_coords = coordinates.SkyCoord(Chambers2014_H2O['GLON'],
                                  Chambers2014_H2O['GLAT'], frame='galactic',
                                  unit=(u.deg, u.deg))
cont_tbl = Table.read(paths.tpath('continuum_photometry.ipac'), format='ascii.ipac')
sgrb2_coords = coordinates.SkyCoord(cont_tbl['RA'], cont_tbl['Dec'],
                                    unit=(u.deg, u.deg), frame='fk5',)
caswell_maser_results = Vizier.query_region(sgrb2_coords.fk5, radius=2*u.arcsec, catalog='VIII/96/catalog')['VIII/96/catalog']
caswell_coords = coordinates.SkyCoord(caswell_maser_results['RAJ2000'].astype('str'), caswell_maser_results['DEJ2000'].astype('str'), frame='fk5', unit=(u.hour, u.deg))

muno_xray_results = Vizier.query_region(sgrb2_coords.fk5, radius=2*u.arcsec,
                                        catalog='J/ApJS/165/173')['J/ApJS/165/173/table2']
muno_xray_coords = coordinates.SkyCoord(muno_xray_results['RAJ2000'],
                                        muno_xray_results['DEJ2000'],
                                        frame='fk5', unit=(u.deg, u.deg))


# old fncont = paths.Fpath('merge/SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase.image.pbcor.fits')
hdu = fits.open(fncont)[0]

mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
#wcsaxes = WCSaxes(mywcs.to_header())
wcsaxes = mywcs

fig = pl.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)
im = ax.imshow(hdu.data.squeeze()*1e3 + 31, cmap=pl.cm.gray_r, origin='lower', vmin=5, vmax=90,
               norm=matplotlib.colors.LogNorm())

axlims = ax.axis()

ax.plot(h2o_coords.fk5.ra.deg, h2o_coords.fk5.dec.deg, color='orange',
        marker='.', linestyle='none', alpha=0.25, transform=ax.get_transform('world'))
ax.plot(caswell_coords.fk5.ra.deg, caswell_coords.fk5.dec.deg, 'r.',
        alpha=0.25, transform=ax.get_transform('world'))
ax.plot(muno_xray_coords.fk5.ra.deg, muno_xray_coords.fk5.dec.deg, 'b.',
        alpha=0.25, transform=ax.get_transform('world'))

ax.axis(axlims)

ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
ra.set_axislabel("RA (J2000)", fontsize=pl.rcParams['axes.labelsize'])
dec.set_axislabel("Dec (J2000)", fontsize=pl.rcParams['axes.labelsize'], minpad=0.0)
ra.ticklabels.set_fontsize(tick_fontsize)
ra.set_ticks(exclude_overlapping=True)
dec.ticklabels.set_fontsize(tick_fontsize)
dec.set_ticks(exclude_overlapping=True)
