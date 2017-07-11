import os
import copy

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.table import Table
from astropy import coordinates
from astropy import wcs
from astropy.convolution import convolve_fft,Gaussian2DKernel
import reproject
import pyregion

from constants import distance, mbar
import pylab as pl

import paths


# define a grid with 1 pc pixels covering the observed area
center = coordinates.SkyCoord('17:47:19.335 -28:23:31.993', frame='fk5', unit=(u.hour, u.deg))
size = 6*u.arcmin
cell_size = 0.25*u.pc
size_pc = (size*distance).to(u.pc, u.dimensionless_angles())

grid_size = int(np.ceil(size_pc/cell_size).value)

grid = np.zeros([grid_size, grid_size])

header = fits.Header()
header['NAXIS'] = 2
header['NAXIS1'] = header['NAXIS2'] = grid_size
header['CRVAL1'] = center.ra.deg
header['CRVAL2'] = center.dec.deg
header['CRPIX1'] = (grid_size+1)/2.
header['CRPIX2'] = (grid_size+1)/2.
header['CDELT1'] = -(cell_size/distance).to(u.deg, u.dimensionless_angles()).value
header['CDELT2'] = (cell_size/distance).to(u.deg, u.dimensionless_angles()).value
header['CUNIT1'] = 'deg'
header['CUNIT2'] = 'deg'
header['CTYPE1'] = 'RA---SIN'
header['CTYPE2'] = 'DEC--SIN'

mywcs = wcs.WCS(header)
x_edges,_ = mywcs.wcs_pix2world(np.arange(0, grid_size+1)-0.5, np.zeros(grid_size+1), 0)
_,y_edges = mywcs.wcs_pix2world(np.zeros(grid_size+1), np.arange(0, grid_size+1)-0.5, 0)

yy,xx = np.indices([grid_size, grid_size])
x_centers,y_centers = mywcs.wcs_pix2world(xx.flat, yy.flat, 0)
grid_coords = coordinates.SkyCoord(x_centers*u.deg, y_centers*u.deg, frame='fk5')


cont_tbl = Table.read(paths.tpath('continuum_photometry.ipac'), format='ascii.ipac')
sgrb2_coords = coordinates.SkyCoord(cont_tbl['RA'], cont_tbl['Dec'],
                                    unit=(u.deg, u.deg), frame='fk5',)

gridded_stars = np.histogram2d(sgrb2_coords.ra.deg, sgrb2_coords.dec.deg, bins=[x_edges[::-1], y_edges])[0][::-1,:].T

hdu = fits.PrimaryHDU(data=gridded_stars, header=header)

hdu.writeto(paths.Fpath('stellar_density_grid.fits'), overwrite=True)

_,nn11_grid,_ = coordinates.match_coordinates_sky(grid_coords, sgrb2_coords, nthneighbor=11)
nn11_grid = nn11_grid.reshape([grid_size, grid_size])
nn11_grid_pc = (nn11_grid * distance).to(u.pc, u.dimensionless_angles())


datapath = '/Users/adam/work/sgrb2/alma/FITS/continuumdata'
herschel25 = fits.open(os.path.join(datapath, 'gcmosaic_column_conv25.fits'))
herschel36 = fits.open(os.path.join(datapath, 'gcmosaic_column_conv36.fits'))

herschel25reproj,_ = reproject.reproject_interp(herschel25, header)
herschel36reproj,_ = reproject.reproject_interp(herschel36, header)

fits.PrimaryHDU(data=herschel25reproj, header=header).writeto(paths.Fpath('other/Herschel25umcolum_regridded_match_stellar.fits'), overwrite=True)
fits.PrimaryHDU(data=herschel36reproj, header=header).writeto(paths.Fpath('other/Herschel36umcolum_regridded_match_stellar.fits'), overwrite=True)

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

ok = np.isfinite(herschel25reproj) & (gridded_stars > 0)
gridded_star_massdensity = (gridded_stars * mbar / (cell_size**2)).to(u.M_sun/u.pc**2)
gas_massdensity25 = (herschel25reproj * 1e22*u.cm**-2 * 2.8*u.Da).to(u.M_sun/u.pc**2)
ax1.loglog(gas_massdensity25[ok], gridded_star_massdensity[ok], 'k.', alpha=0.7, markeredgecolor=(0,0,0,0.5))

logas = (~np.isfinite(herschel25reproj)) & (gridded_stars > 0)
ax1.plot(np.nanmax(gas_massdensity25) * np.ones(logas.sum()),
         gridded_star_massdensity[logas],
         '>')

lostars = (np.isfinite(herschel25reproj)) & (gridded_stars == 0)
ax1.plot(gas_massdensity25[lostars],
         np.nanmin(gridded_star_massdensity[ok])*0.5*np.ones(lostars.sum()),
         'v')
ax1.loglog([1e3,1e6], [1e0, 1e5], 'k--')

ax1.fill_between([0.1, 1e5],
                 np.array([0.1, 1e5])**2.67/(100**2.67),
                 np.array([0.1, 1e5])**2.67/(100**2.67)*10,
                 alpha=0.5,
                 color='orange',
                 label='Mon R2')
ax1.fill_between([0.1, 1e5],
                 np.array([0.1, 1e5])**1.87/(100**1.87),
                 np.array([0.1, 1e5])**1.87/(100**1.87)*10,
                 color='blue',
                 alpha=0.5,
                 label='Ophiucus')
ax1.plot([0.1, 1e5], np.array([0.1, 1e5])**1.87/(1e4**1.87)*(1e4**(5/3.)/1e5), 'b:', linewidth=3, alpha=0.5)

ax1.set_ylabel("Gridded NN11 Stellar Surface Density\n[M$_\odot$ pc$^{-2}$]")
ax1.set_xlabel("Herschel-derived Surface Density [M$_\odot$ pc$^{-2}$]")
ax1.plot([0.1, 1e5], np.array([0.1, 1e5])*1e-2, 'r-', linewidth=3, alpha=0.5, zorder=-10)
ax1.axis([1e3,1e5,1e0,1e5])
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel.png"), bbox_inches='tight')
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel.pdf"), bbox_inches='tight')

nn = 11
nn11_msunpersqpc = ((nn-1) * mbar / (np.pi*nn11_grid_pc)**2).to(u.M_sun/u.pc**2)

fig2 = pl.figure(2, figsize=(10,10))
fig2.clf()
ax2 = fig2.gca()

ax2.loglog(gas_massdensity25.ravel().value, nn11_msunpersqpc.ravel().value, 'k.', alpha=0.7, markeredgecolor=(0,0,0,0.5))
lims = ax2.axis()
ax2.loglog([1e3,1e6], [1e0, 1e5], 'k--')

ax2.fill_between([0.1, 1e5],
                 np.array([0.1, 1e5])**2.67/(100**2.67),
                 np.array([0.1, 1e5])**2.67/(100**2.67)*10,
                 alpha=0.5,
                 color='green',
                 label='Mon R2')
ax2.fill_between([0.1, 1e5],
                 np.array([0.1, 1e5])**1.87/(100**1.87),
                 np.array([0.1, 1e5])**1.87/(100**1.87)*10,
                 color='blue',
                 alpha=0.5,
                 label='Ophiucus')
ax2.plot([0.1, 1e5], np.array([0.1, 1e5])**1.87/(1e4**1.87)*(1e4**(5/3.)/1e5), 'b:', linewidth=3, alpha=0.5)
ax2.plot([0.1, 1e5], np.array([0.1, 1e5])*1e-2, 'r-', linewidth=3, alpha=0.5, zorder=-10)

ax2.set_ylabel("Gridded NN11 Stellar Surface Density\n[M$_\odot$ pc$^{-2}$]")
ax2.set_xlabel("Herschel-derived Surface Density [M$_\odot$ pc$^{-2}$]")
ax2.axis(lims)
ax2.axis([1e3,1e5,1e0,1e5])
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel.png"), bbox_inches='tight')
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel.pdf"), bbox_inches='tight')
