import paths
from astropy.io import fits
import pylab as pl
from astropy import wcs
from astropy import visualization
import matplotlib
matplotlib.use('Qt5Agg')
assert matplotlib.get_backend() == 'Qt5Agg'

pl.matplotlib.rc_file('pubfiguresrc')

pl.rcParams['figure.figsize'] = (12,8)
pl.rcParams['figure.dpi'] = 75.
pl.rcParams['savefig.dpi'] = 300.
pl.rcParams['axes.labelsize'] = 12
pl.rcParams['xtick.labelsize'] = 12
pl.rcParams['ytick.labelsize'] = 12
if matplotlib.__version__[0] == '1':
    markersize = 10
    tick_fontsize = 12
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
cb.ax.tick_params(labelsize=12)

ax.contour(hc3ndata, colors=['orange']*5,
           levels=[3,6,9,12,15,18],
           linewidths=1,
           alpha=0.5, transform=ax.get_transform(hc3nwcs))
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

fig.savefig(paths.fpath("HC3N_contours_on_SCUBA_column.png"), bbox_inches='tight')
