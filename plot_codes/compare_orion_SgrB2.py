from astropy.io import fits
import paths
import os
from astropy.nddata import Cutout2D
from astropy import coordinates
from astropy import units as u
from astropy import wcs
import pylab as pl
from visualization import make_scalebar, hide_labels
from constants import distance as d_sgr
from files import contfilename

pl.rcParams['figure.figsize'] = (12,8)
pl.rcParams['figure.dpi'] = 75.
pl.rcParams['savefig.dpi'] = 300.
pl.rcParams['axes.labelsize'] = 9
pl.rcParams['xtick.labelsize'] = 8
pl.rcParams['ytick.labelsize'] = 8

d_ori = 415*u.pc

fh_o = fits.open(os.path.join(paths.root,'other_continuum/Orion_MUSTANG_at_SgrB2.fits'))
fh_o = fits.open(os.path.join(paths.root,'other_continuum/OrionClean9.0-24nov08.fits.gz'))
fh_s = fits.open(contfilename)

fh_o[0].data *= 1000 # Jy->mJy
fh_s[0].data *= 1000 # Jy->mJy

sgrb2_hii_T_coord = coordinates.SkyCoord(266.8642831*u.deg, -28.35005004*u.deg, frame='fk5')
sgrb2_hii_T_cutout = Cutout2D(fh_s[0].data, sgrb2_hii_T_coord, 12*u.arcsec, wcs=wcs.WCS(fh_s[0].header))
sgrb2_hii_L_coord = coordinates.SkyCoord(266.8442039*u.deg, -28.36499239*u.deg, frame='fk5')
sgrb2_hii_L_cutout = Cutout2D(fh_s[0].data, sgrb2_hii_L_coord, 12*u.arcsec, wcs=wcs.WCS(fh_s[0].header))

orion_center = coordinates.SkyCoord(83.82187403*u.deg, -5.387010177*u.deg, frame='fk5')
orion_cutout = Cutout2D(fh_o[0].data.squeeze(), orion_center,
                        12*u.arcsec*d_sgr/d_ori,
                        wcs=wcs.WCS(fh_o[0].header).celestial)

fig = pl.figure(1)
fig.clf()
ax = pl.subplot(1,3,1, projection=wcs.WCS(fh_s[0].header))
pl.imshow(sgrb2_hii_T_cutout.data, vmin=-1e-1, vmax=4e-0, cmap='viridis',
          origin='lower', interpolation='nearest')
pl.title("Sgr B2 HII T")
hide_labels(ax)
make_scalebar(ax, left_side=coordinates.SkyCoord('266d54m39s', '-28d27m48.5s',
                                                 unit=(u.deg, u.deg), frame='fk5'),
              length=(0.1*u.pc / d_sgr).to(u.arcsec, u.dimensionless_angles()),
              label='0.1 pc',
              fontsize=8,
             )

pl.subplot(1,3,2)
pl.imshow(orion_cutout.data*(d_ori/d_sgr).decompose()**2, vmin=-1e-1,
          vmax=4e-0, cmap='viridis', origin='lower', interpolation='nearest')
pl.title("M42")
pl.gca().get_xaxis().set_ticks([])
pl.gca().get_yaxis().set_ticks([])

ax = pl.subplot(1,3,3, projection=wcs.WCS(fh_s[0].header))
im = ax.imshow(sgrb2_hii_L_cutout.data, vmin=-1e-1, vmax=4e-0, cmap='viridis',
               origin='lower', interpolation='nearest')
ax.set_title("Sgr B2 HII L")
hide_labels(ax)
make_scalebar(ax, left_side=coordinates.SkyCoord('266d54m30s', '-28d27m48.5s',
                                                 unit=(u.deg, u.deg), frame='fk5'),
              length=(0.1*u.pc / d_sgr).to(u.arcsec, u.dimensionless_angles()),
              label='0.1 pc',
              fontsize=8,
             )
fig.canvas.draw()
print(ax.bbox._bbox.x1 - 0.9)
print(ax.bbox._bbox.y0 - 0.32404411764705876)
x1 = ax.bbox.x1 / (fig.bbox.x1-fig.bbox.x0)
y0 = ax.bbox.y0 / (fig.bbox.y1-fig.bbox.y0)
height = (ax.bbox.y1-ax.bbox.y0) / (fig.bbox.y1-fig.bbox.y0)
cax_bbox = [x1 + 0.02, y0, 0.02, height]
print(cax_bbox)
cb = pl.colorbar(mappable=im, cax=pl.axes(cax_bbox))
cb.set_label('mJy/beam')
pl.savefig(paths.fpath("Orion_SgrB2HII_side_by_side.pdf"), bbox_inches='tight')
