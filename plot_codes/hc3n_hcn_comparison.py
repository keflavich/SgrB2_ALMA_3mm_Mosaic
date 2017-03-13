import aplpy
import paths
import numpy as np
from astropy import wcs
from wcsaxes import WCS as WCSaxes
from astropy.io import fits
import reproject

import pylab as pl

pl.matplotlib.rc_file('pubfiguresrc')

pl.rcParams['figure.figsize'] = (12,8)
pl.rcParams['figure.dpi'] = 300.
pl.rcParams['savefig.dpi'] = 300.
pl.rcParams['axes.labelsize'] = 9
pl.rcParams['xtick.labelsize'] = 8
pl.rcParams['ytick.labelsize'] = 8
pl.rcParams['font.size'] = 6
tick_fontsize = 6


fnhc3n = paths.Fpath('merge/max/SgrB2_b3_7M_12M.HC3N.image.pbcor_max_medsub.fits')
fnhcn = paths.Fpath('merge/max/SgrB2_b3_7M_12M.HCN.image.pbcor_max_medsub.fits')



# FF2 = aplpy.FITSFigure(fnhc3n)
# FF2.show_grayscale(vmax=0.06, invert=True)
# FF2.show_contour(fnhcn, levels=[0.02,0.03,0.04,0.05,0.06,0.1,0.2,0.3,0.4],
#                  colors=[(1,0,0,x) for x in np.linspace(0,1,10)],
#                  filled=True, )
# FF2.save(paths.fpath("HC3N_grayscale_with_HCN_red_filled_contours.png"), dpi=150)
# FF2.recenter(266.8318615, -28.3940598, width=0.062, height=0.105)
# FF2.save(paths.fpath("HC3N_grayscale_with_HCN_red_filled_contours_zoom.png"), dpi=150)

hdu = fits.open(fnhc3n)[0]
mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
wcsaxes = WCSaxes(mywcs.to_header())

outhdr = fits.getheader(fnhc3n)
hcn = reproject.reproject_interp(fnhcn, outhdr, order=1)[0]

fig = pl.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)
im = ax.imshow(hdu.data.squeeze(), cmap=pl.cm.gray_r, origin='lower', vmax=0.06)
con = ax.contourf(hcn, levels=[0.02,0.03,0.04,0.05,0.06,0.1,0.2,0.3,0.4],
                  colors=[(1,0,0,x) for x in np.linspace(0.2,0.7,10)], origin='lower')
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
ra.set_axislabel("RA (J2000)", fontsize=pl.rcParams['axes.labelsize'],)
dec.set_axislabel("Dec (J2000)", fontsize=pl.rcParams['axes.labelsize'], minpad=0)
ra.ticklabels.set_fontsize(tick_fontsize)
ra.set_ticks(exclude_overlapping=True)
dec.ticklabels.set_fontsize(tick_fontsize)
dec.set_ticks(exclude_overlapping=True)


fig.savefig(paths.fpath("HC3N_grayscale_with_HCN_red_filled_contours.png"), dpi=300, bbox_inches='tight')

tr_fk5 = ax.get_transform("fk5")
tr_pix = ax.get_transform("pixel")
xmin,xmax,ymin,ymax = [266.8318615-0.062/2, 266.8318615+0.062/2, -28.3940598-0.105/2, -28.3940598+0.105/2]
xmin_,ymin_ = mywcs.wcs_world2pix(xmin, ymin, 0) #tr_pix.transform([xmin,ymin])
xmax_,ymax_ = mywcs.wcs_world2pix(xmax, ymax, 0) #tr_pix.transform([xmax,ymax])
print((xmin_,xmax_,ymin_,ymax_))
for blah in (xmin_,xmax_,ymin_,ymax_):
    assert blah > 0 and blah < 4096
ax.axis([xmax_,xmin_,ymin_,ymax_])
fig.savefig(paths.fpath("HC3N_grayscale_with_HCN_red_filled_contours_zoom.png"), dpi=300, bbox_inches='tight')
