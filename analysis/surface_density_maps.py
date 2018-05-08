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

from constants import distance, mass_represented_by_a_source

from gutermuth2011_law import gas_depletion_law, sigma_gas_of_t
import lada2017relation
import elmegreen2018
import label_lines

import paths

import pylab as pl

pl.close('all')
pl.rcParams['figure.figsize'] = (12,8)
pl.rcParams['figure.dpi'] = 75
pl.rcParams['savefig.dpi'] = 150
pl.rcParams['axes.titlesize']= 40
pl.rcParams['axes.labelsize']= 24
pl.rcParams['xtick.labelsize']= 20
pl.rcParams['ytick.labelsize']= 20

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

gridded_stars = np.histogram2d(sgrb2_coords.ra.deg, sgrb2_coords.dec.deg,
                               bins=[x_edges[::-1], y_edges])[0][::-1,:].T

hdu = fits.PrimaryHDU(data=gridded_stars, header=header)

hdu.writeto(paths.Fpath('stellar_density_grid.fits'), overwrite=True)
hdu.data *= mass_represented_by_a_source.value
hdu.writeto(paths.Fpath('stellar_mass_density_grid.fits'), overwrite=True)

_,nn11_grid,_ = coordinates.match_coordinates_sky(grid_coords, sgrb2_coords,
                                                  nthneighbor=11)
nn11_grid = nn11_grid.reshape([grid_size, grid_size])
nn11_grid_pc = (nn11_grid * distance).to(u.pc, u.dimensionless_angles())



herschel25 = fits.open(paths.cpath('gcmosaic_column_conv25.fits'))
herschel36 = fits.open(paths.cpath('gcmosaic_column_conv36.fits'))

herschel25reproj,_ = reproject.reproject_interp(herschel25, header)
herschel36reproj,_ = reproject.reproject_interp(herschel36, header)

gas_massdensity25 = (herschel25reproj * 1e22*u.cm**-2 * 2.8*u.Da).to(u.M_sun/u.pc**2)
gas_massdensity36 = (herschel36reproj * 1e22*u.cm**-2 * 2.8*u.Da).to(u.M_sun/u.pc**2)

header['BUNIT'] = 'Msun / pc^2'
fits.PrimaryHDU(data=gas_massdensity25.value,
                header=header).writeto(paths.Fpath('other/Herschel25umcolum_regridded_match_stellar.fits'),
                                       overwrite=True)
fits.PrimaryHDU(data=gas_massdensity36.value,
                header=header).writeto(paths.Fpath('other/Herschel36umcolum_regridded_match_stellar.fits'),
                                       overwrite=True)

fig1 = pl.figure(1)
fig1.clf()
ax1 = fig1.gca()

ok = np.isfinite(herschel25reproj) & (gridded_stars > 0)
gridded_star_massdensity = (gridded_stars * mass_represented_by_a_source / (cell_size**2)).to(u.M_sun/u.pc**2)

pixdots, = ax1.loglog(gas_massdensity25[ok], gridded_star_massdensity[ok],
                      'k.', alpha=0.7, markeredgecolor=(0,0,0,0.5))

logas = (~np.isfinite(herschel25reproj)) & (gridded_stars > 0)
ax1.plot(np.nanmax(gas_massdensity25) * np.ones(logas.sum()),
         gridded_star_massdensity[logas],
         '>')

lostars = (np.isfinite(herschel25reproj)) & (gridded_stars == 0)
lostarvs = ax1.plot(gas_massdensity25[lostars],
                    np.nanmin(gridded_star_massdensity[ok])*0.5*np.ones(lostars.sum()),
                    '^')
# 5/3 slope
#ax1.loglog([1e3,1e6], [3e0, 3e5], 'k--')


ax1.axis([1e3,1e5,1e0,1e5])
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_nomodel_nolocal.png"), bbox_inches='tight')
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_nomodel_nolocal.pdf"), bbox_inches='tight')

monr2_fill = ax1.fill_between([0.1, 1e5],
                              np.array([0.1, 1e5])**2.67/(100**2.67),
                              np.array([0.1, 1e5])**2.67/(100**2.67)*10,
                              alpha=0.5,
                              color='orange',
                              label='Mon R2')
oph_fill = ax1.fill_between([0.1, 1e5],
                            np.array([0.1, 1e5])**1.87/(100**1.87),
                            np.array([0.1, 1e5])**1.87/(100**1.87)*10,
                            color='blue',
                            alpha=0.5,
                            label='Ophiuchus')
local_plotobjs = [monr2_fill, oph_fill]

ax1.plot([550,500,15],[300,30,0.4],'bo')
ax1.plot([300,300,15],[500,54,0.22],'go')
ax1.plot([0.1, 1e5], np.array([0.1, 1e5])**1.87/(1e4**1.87)*(1e4**(5/3.)/1e5),
         'b:', linewidth=3, alpha=0.5)

ax1.axis([1e3,1e5,1e0,1e5])
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_nomodel.png"), bbox_inches='tight')
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_nomodel.pdf"), bbox_inches='tight')

sigma_gas = np.logspace(1,6) * u.M_sun / u.pc**2

model_plotobjs = []
model_labels = []
for time in (0.01, 0.1, 0.74)*u.Myr:
    model_plotobjs += ax1.loglog(sigma_gas_of_t(sigma_gas, time, alpha=1,
                                                k=0.1/u.Myr),
                                 gas_depletion_law(sigma_gas, time, alpha=1,
                                                   k=0.1/u.Myr), label=time,
                                 color='r', linewidth=3, alpha=0.5,
                                 zorder=-10,)
    model_plotobjs += ax1.loglog(sigma_gas_of_t(sigma_gas, time),
                                 gas_depletion_law(sigma_gas, time),
                                 label=time, color='orange', linewidth=3,
                                 alpha=0.5, zorder=-10,)
    model_labels.append(time)

ax1.set_ylabel("Gridded NN11 Stellar Surface Density\n$\Sigma_*$ [M$_\odot$ pc$^{-2}$]", fontsize=24)
ax1.set_xlabel("Gas Surface Density $\Sigma_{gas}$ [M$_\odot$ pc$^{-2}$]", fontsize=24)
#ax1.plot([0.1, 1e5], np.array([0.1, 1e5])*4e-2, 'r-', linewidth=3, alpha=0.5, zorder=-10)
ax1.axis([1e3,1e5,1e0,1e5])
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel.png"), bbox_inches='tight')
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel.pdf"), bbox_inches='tight')
ax1.axis([1e0,1e5,1e-1,1e5])
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_full.png"), bbox_inches='tight')


for obj in local_plotobjs:
    obj.set_visible(False)


ax1.axis([1e3,1e5,1e0,1e5])
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_nolocal.pdf"), bbox_inches='tight')

for obj in local_plotobjs:
    obj.set_visible(True)


ax1.axis([1e3,1e5,1e0,1e5])
ax1.loglog(np.logspace(3,5),
           lada2017relation.sigma_star_california(np.logspace(3,5)*u.M_sun/u.pc**2),
           linewidth=3, color='m', label='Lada2017_cali')
ax1.loglog(np.logspace(3,5),
           lada2017relation.sigma_star_orionA(np.logspace(3,5)*u.M_sun/u.pc**2),
           linewidth=3, linestyle='--', color='m', label='Lada2017_orionA')
ax1.loglog(np.logspace(3,5),
           lada2017relation.sigma_star_orionB(np.logspace(3,5)*u.M_sun/u.pc**2),
           linewidth=3, linestyle=':', color='m', label='Lada2017_orionB')

fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_withLada2017.png"), bbox_inches='tight')
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_withLada2017.pdf"), bbox_inches='tight')

for obj in local_plotobjs+model_plotobjs:
    obj.set_visible(False)

fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_nomodel_nolocal_withLada2017.png"), bbox_inches='tight')
fig1.savefig(paths.fpath("stellar_vs_gas_column_density_gridded_herschel_nomodel_nolocal_withLada2017.pdf"), bbox_inches='tight')


nn = 11
nn11_msunpersqpc = ((nn-1) * mass_represented_by_a_source /
                    (np.pi*(nn11_grid_pc)**2)).to(u.M_sun/u.pc**2)

hdu = fits.PrimaryHDU(data=nn11_msunpersqpc.value, header=header)
hdu.writeto(paths.Fpath('nn11_stellar_massdensity_grid.fits'), overwrite=True)

fig2 = pl.figure(2, figsize=(10,10))
fig2.clf()
ax2 = fig2.gca()

ax2.loglog(gas_massdensity25.ravel().value,
           nn11_msunpersqpc.ravel().value,
           'k.', alpha=0.5, markeredgecolor=(0,0,0,0.5))
lims = ax2.axis()
# 5/3 slope
#ax2.loglog([1e3,1e6], [3e0, 3e5], 'k--')

# Handle lower limits
logas = (~np.isfinite(gas_massdensity25)) & (nn11_msunpersqpc > 0)
ax2.plot(np.nanmax(gas_massdensity25) * np.ones(logas.sum()),
         nn11_msunpersqpc[logas],
         markeredgecolor='k',
         markerfacecolor='none',
         linestyle='none',
         marker='>')

# No limits are possible here.
# lostars = (np.isfinite(herschel25reproj)) & (nn11_msunpersqpc == 0)
# ax2.plot(gas_massdensity25[lostars],
#          np.nanmin(nn11_msunpersqpc[ok])*0.5*np.ones(lostars.sum()),
#          'v')

ax2.set_ylabel("Gridded NN11 Stellar Surface Density\n$\Sigma_*$ [M$_\odot$ pc$^{-2}$]", fontsize=24)
ax2.set_xlabel("Gas Surface Density $\Sigma_{gas}$ [M$_\odot$ pc$^{-2}$]", fontsize=24)

ax2.axis([1e3,1e5,1e0,1e5])
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_nomodel_nolocal.png"), bbox_inches='tight')
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_nomodel_nolocal.pdf"), bbox_inches='tight')

monr2_lowerline =np.array([0.1, 1e5])**2.67/(100**2.67) * 2.5
monr2_fill = ax2.fill_between([0.1, 1e5],
                              monr2_lowerline,
                              monr2_lowerline*10,
                              alpha=0.5,
                              color='green',
                              label='Mon R2')
oph_lowerline = np.array([0.1, 1e5])**1.87/(100**1.87) * 1.5
oph_fill = ax2.fill_between([0.1, 1e5],
                            oph_lowerline,
                            oph_lowerline*10,
                            color='blue',
                            alpha=0.5,
                            label='Ophiuchus')
oph_scalefactor = 50.
oph_lowerline_line = ax2.plot([0.1, 1e5], oph_lowerline/oph_scalefactor, 'b:', linewidth=3, alpha=0.5)
#ax2.plot([0.1, 1e5], np.array([0.1, 1e5])*4e-2, 'r-', linewidth=3, alpha=0.5, zorder=-10)

monr2_text = ax2.text(2.3e3, 9e3, "Mon R2", color='k', rotation=50,
                      fontsize=18, verticalalignment='bottom',
                      horizontalalignment='center')
oph_text = ax2.text(2.5e3, 4.5e2, "Ophiuchus", color='k', rotation=38,
                    fontsize=18, verticalalignment='bottom',
                    horizontalalignment='center')

local_plotobjs = [monr2_text, oph_text, monr2_fill, oph_fill,] + oph_lowerline_line


ax2.axis([1e3,1e5,1e0,1e5])
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_nomodel.png"), bbox_inches='tight')
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_nomodel.pdf"), bbox_inches='tight')

model_labels = []
model_plotobjs = []
for time in (0.01, 0.1, 0.74)*u.Myr:
    model_plotobjs += ax2.loglog(sigma_gas_of_t(sigma_gas, time, alpha=1,
                                                k=0.1/u.Myr),
                                 gas_depletion_law(sigma_gas, time, alpha=1,
                                                   k=0.1/u.Myr), label=time,
                                 color='r', linewidth=3, alpha=0.5,
                                 zorder=-10,)
    model_plotobjs += ax2.loglog(sigma_gas_of_t(sigma_gas, time),
                                 gas_depletion_law(sigma_gas, time),
                                 label=time, color='orange', linewidth=3,
                                 alpha=0.5, zorder=-10,)

    model_labels.append(time)


ax2.axis(lims)

# Arrows showing the shift if you subtract off the "most aggressive plausible"
# uniform foreground value
bg_5e22 = (5e22*u.cm**-2*2.8*u.Da).to(u.M_sun/u.pc**2)
ax2.arrow(3e3, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(3e3-bg_5e22.value), color='k')
ax2.arrow(1e4, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(1e4-bg_5e22.value), color='k')
ax2.arrow(3e4, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(3e4-bg_5e22.value), color='k')
#ax2.arrow(3e5, 1.1, -(3e5-bg_5e22.value), 0, head_width=6, shape='left')
ax2.axis([1e3,1e5,1e0,1e5])
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel.png"), bbox_inches='tight')
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel.pdf"), bbox_inches='tight')

ax2.plot([550,500,15],[300,30,0.4],'bo')
ax2.plot([300,300,15],[500,54,0.22],'go')
ax2.axis([1e0,1e5,1e-1,1e5])
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_full.png"), bbox_inches='tight')

for obj in local_plotobjs:
    obj.set_visible(False)

line_labels = label_lines.label_lines(model_plotobjs,
                                      xvals=[7e4, 5e4, 7e4, 6e3, 7e4, 1.2e3],
                                      backgroundcolor=(1,1,1,0.8),
                                      fontsize=14)
line_labels.append(ax2.text(3e4, 4e0, '$\\alpha=1$', color='r', fontsize=20, alpha=0.9))
line_labels.append(ax2.text(1.1e4, 4e4, '$\\alpha=2$', color='orange', fontsize=20, alpha=0.9))

ax2.axis([1e3,1e5,1e0,1e5])
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_nolocal.pdf"), bbox_inches='tight')

for obj in line_labels:
    obj.set_visible(False)

for obj in local_plotobjs:
    obj.set_visible(True)


ax2.axis([1e3,1e5,1e0,1e5])
lada_cali = ax2.loglog(np.logspace(3,5),
                       lada2017relation.sigma_star_california(np.logspace(3,5)*u.M_sun/u.pc**2),
                       linewidth=3, color='m', label='Lada2017_cali')
lada_oria = ax2.loglog(np.logspace(3,5),
                       lada2017relation.sigma_star_orionA(np.logspace(3,5)*u.M_sun/u.pc**2),
                       linewidth=3, linestyle='--', color='m',
                       label='Lada2017_orionA')
lada_orib = ax2.loglog(np.logspace(3,5),
                       lada2017relation.sigma_star_orionB(np.logspace(3,5)*u.M_sun/u.pc**2),
                       linewidth=3, linestyle=':', color='m',
                       label='Lada2017_orionB')

fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_withLada2017.png"), bbox_inches='tight')
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_withLada2017.pdf"), bbox_inches='tight')

for obj in local_plotobjs+model_plotobjs:
    obj.set_visible(False)

fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_nomodel_nolocal_withLada2017.png"), bbox_inches='tight')
fig2.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_nomodel_nolocal_withLada2017.pdf"), bbox_inches='tight')


# for obj in [lada_cali, lada_orib, lada_oria]:
#     obj.set_visible(False)
# 
# 
#can't plot this: don't have a prediction for sigma_star, just sigma_sfr
# ax2.loglog(np.logspace(3,5),
#            elmegreen2018.
# 





total_gas_mass = np.nansum(cell_size**2 * gas_massdensity25)
total_gas_mass_bgsub = np.nansum(cell_size**2 * (gas_massdensity25-bg_5e22))
total_stellar_mass = len(cont_tbl) * mass_represented_by_a_source
total_mass = (total_stellar_mass+total_gas_mass)
total_mass_bgsub = (total_stellar_mass+total_gas_mass_bgsub)

print("Total SFE = {0} / {1} = {2}".format(total_stellar_mass, total_mass,
                                           total_stellar_mass/total_mass,
                                          ))

print("Total SFE BGsub = {0} / {1} = {2}"
      .format(total_stellar_mass, total_mass_bgsub,
              total_stellar_mass/total_mass_bgsub,))




fig3 = pl.figure(3, figsize=(10,10))
fig3.clf()
ax3 = fig3.gca()

timescale = elmegreen2018.tff(2 * 200*u.M_sun/u.pc**2 / (10*u.pc))
# Diederik says  use higher density
timescale = elmegreen2018.tff(2 * 1000*u.M_sun/u.pc**2 / (10*u.pc))
#timescale = 1.8*u.Myr

ax3.loglog(gas_massdensity25.ravel().value,
           (nn11_msunpersqpc.ravel()/timescale).to(u.Msun/u.pc**2/u.Myr).value,
           'k.', alpha=0.5, markeredgecolor=(0,0,0,0.5))
lims = ax3.axis()
# 5/3 slope
#ax2.loglog([1e3,1e6], [3e0, 3e5], 'k--')

# Handle lower limits
logas = (~np.isfinite(gas_massdensity25)) & (nn11_msunpersqpc > 0)
ax3.plot(np.nanmax(gas_massdensity25) * np.ones(logas.sum()),
         (nn11_msunpersqpc[logas].ravel()/u.Myr).to(u.Msun/u.pc**2/u.Myr).value,
         markeredgecolor='k',
         markerfacecolor='none',
         linestyle='none',
         marker='>')

ax3.plot(np.logspace(3,5),
         elmegreen2018.Sigma_sfr_eqn21(np.logspace(3,5)*u.M_sun/u.pc**2, epsilon_ff=0.001).to(u.Msun/u.pc**2/u.Myr),
         label='$\epsilon_{ff}=0.001$',
         linestyle='--',
        )
ax3.plot(np.logspace(3,5),
         elmegreen2018.Sigma_sfr_eqn21(np.logspace(3,5)*u.M_sun/u.pc**2, epsilon_ff=0.01).to(u.Msun/u.pc**2/u.Myr),
         label='$\epsilon_{ff}=0.01$',
         linestyle=':',
        )
ax3.plot(np.logspace(3,5),
         elmegreen2018.Sigma_sfr_eqn21(np.logspace(3,5)*u.M_sun/u.pc**2, epsilon_ff=0.1).to(u.Msun/u.pc**2/u.Myr),
         label='$\epsilon_{ff}=0.1$',
         linewidth=2, alpha=0.5,
        )
ax3.set_xlabel("Column Density $\Sigma_{gas}$ [M$_\odot$ pc$^{-2}$]")
ax3.set_ylabel("SFR Surface Density $\Sigma_{SFR}$ [M$_\odot$ pc$^{-2}$ Myr$^{-1}$]"
               "\nassuming $t_{sf} = t_{ff}(\\rho_{midplane}) = "+
               "{0:0.1f}".format(timescale.value)+"$ Myr")
ax3.axis(lims)

pl.legend(loc='best', fontsize=18)
fig3.savefig(paths.fpath("stellar_vs_gas_column_density_gridNN11_herschel_Elmegreen2018_SFR_models.pdf"), bbox_inches='tight')
