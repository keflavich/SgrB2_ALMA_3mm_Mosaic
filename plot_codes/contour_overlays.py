import paths
from astropy.io import fits
from astropy.table import Table
from astropy import coordinates
from astropy import units as u
import pylab as pl
from astropy import wcs

from visualization import make_scalebar, hide_scalebar
from constants import distance

import matplotlib
matplotlib.use('Qt5Agg')
assert matplotlib.get_backend() == 'Qt5Agg'

pl.matplotlib.rc_file('pubfiguresrc')

textsize = 14
pl.rcParams['figure.figsize'] = (12,8)
pl.rcParams['figure.dpi'] = 75.
pl.rcParams['savefig.dpi'] = 300.
pl.rcParams['axes.labelsize'] = textsize
pl.rcParams['xtick.labelsize'] = textsize
pl.rcParams['ytick.labelsize'] = textsize
if matplotlib.__version__[0] == '1':
    markersize = 10
    tick_fontsize = textsize
elif matplotlib.__version__[0] == '2':
    markersize = 0.5
    tick_fontsize = 6


fnscubacolmap = paths.cpath('column_maps/scuba_col_herscheltem.fits')
fnhc3n = paths.Fpath('merge/max/SgrB2_b3_7M_12M.HC3N.image.pbcor_medsub_max_K.fits')

scubafh = fits.open(fnscubacolmap)[0]
hc3nfh = fits.open(fnhc3n)[0]

scubawcs = wcs.WCS(scubafh.header)

hc3nslice = slice(100,-100),slice(100,-100)
hc3ndata = hc3nfh.data[hc3nslice]
hc3nwcs = wcs.WCS(hc3nfh.header)[hc3nslice]

fig = pl.figure(1)
fig.clf()
ax = fig.add_subplot(111, projection=scubawcs)
im = ax.imshow(scubafh.data, vmin=5e22, vmax=5e24,
               norm=matplotlib.colors.LogNorm(),
               cmap='gray_r',
               interpolation='nearest',
               origin='lower',
               )
cb = fig.colorbar(im)
cb.ax.tick_params(labelsize=textsize)
cb.set_label("$N(H_2)$ [cm$^{-2}$]")

ax.contour(hc3ndata, colors=['red']*7,
           levels=[3,7,11,15,19,23],
           linewidths=1,
           alpha=0.8, transform=ax.get_transform(hc3nwcs))
ax.axis([195,340,180,320])

ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
ra.set_axislabel("RA (J2000)", fontsize=pl.rcParams['axes.labelsize'])
dec.set_axislabel("Dec (J2000)", fontsize=pl.rcParams['axes.labelsize'], minpad=0)
ra.ticklabels.set_fontsize(tick_fontsize)
ra.set_ticks(exclude_overlapping=True)
dec.ticklabels.set_fontsize(tick_fontsize)
dec.set_ticks(exclude_overlapping=True)

scalebarpos = coordinates.SkyCoord("17:47:27.5", "-28:26:00.0",
                                   unit=(u.h, u.deg), frame='fk5')
make_scalebar(ax, scalebarpos,
              length=(2.0*u.pc / distance).to(u.arcsec,
                                              u.dimensionless_angles()),
              color='k',
              label='2 pc',
              text_offset=1.0*u.arcsec,
             )

fig.savefig(paths.fpath("HC3N_contours_on_SCUBA_column.png"), bbox_inches='tight')

core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
cores = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                             frame='fk5')
markersize = 6
tr_fk5 = ax.get_transform("fk5")
coredots, = ax.plot(cores.ra, cores.dec, '.', color='lime', transform=tr_fk5,
                    markersize=markersize, zorder=50)
ax.axis([195,340,180,320])

fig.savefig(paths.fpath("HC3N_contours_on_SCUBA_column_withcores.png"),
            bbox_inches='tight')
