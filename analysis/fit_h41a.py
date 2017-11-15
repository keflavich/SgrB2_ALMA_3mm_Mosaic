import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u, table, coordinates
import pyspeckit
import pyregion
import paths
import pylab as pl

regs = pyregion.open(paths.rpath('SgrB2_1.3cm_hiiRegions_shiftedtoALMA.reg'))

cube = SpectralCube.read(paths.Fpath('M_cutouts/SgrB2_b3_12M_TE.H41a.image.pbcor_M_cutout.medsub.fits'))

fit_values = {}
coords = {}

for reg in regs:
    if 'text' not in reg.attr[1]:
        continue
    name = reg.attr[1]['text'].strip("{}")
    print(name)

    try:
        cutout = cube.subcube_from_ds9region(pyregion.ShapeList([reg])).to(u.K,
                                                                           cube.beam.jtok_equiv(cube.with_spectral_unit(u.GHz).spectral_axis))
    except ValueError as ex:
        print(ex)
        continue

    spectrum = cutout.mean(axis=(1,2))

    sp = pyspeckit.Spectrum.from_hdu(spectrum.hdu)
    sp.specname = name

    sp.plotter(figure=pl.figure(1))

    sp.baseline(order=0)
    sp.specfit(fittype='gaussian')
    sp.baseline(excludefit=True, order=0)
    sp.specfit(fittype='gaussian')

    sp.plotter.savefig(paths.fpath("spectra/h41afits/{0}_h41a_fit.png".format(name)))

    if sp.specfit.parinfo['WIDTH0'].value > 4 and sp.specfit.parinfo['AMPLITUDE0'].value > 0:
        fit_values[name] = sp.specfit.parinfo
        coords[name] = coordinates.SkyCoord(reg.coord_list[0],
                                            reg.coord_list[1],
                                            unit=(u.deg, u.deg),
                                            frame=reg.coord_format)

names = sorted(fit_values.keys())
columns = [table.Column(name=par, data=[fit_values[key][par].value for key in names])
           for par in sp.specfit.parinfo.keys()]
errcolumns = [table.Column(name='e_'+par, data=[fit_values[key][par].error for key in names])
              for par in sp.specfit.parinfo.keys()]
ra, dec = (table.Column(name='RA', data=u.Quantity([coords[key].ra for key in names])),
           table.Column(name='Dec', data=u.Quantity([coords[key].dec for key in names])))

allcols = [table.Column(data=names, name='Name'), ra, dec] + [z for x,y in zip(columns, errcolumns) for z in (x,y)]
tbl = table.Table(allcols)

ok = tbl['AMPLITUDE0'] > tbl['e_AMPLITUDE0']*3

tbl[ok].write(paths.tpath('H41a_fits.ipac'), format='ascii.ipac', overwrite=True)

center = coordinates.SkyCoord('17:47:20.174', '-28:23:04.233', frame='fk5', unit=(u.hour, u.deg))

depree = table.Table.read(paths.tpath('DePree2011.txt'), format='ascii.fixed_width', delimiter='|')

tbl.rename_column('Name','Source')
tbl.rename_column('SHIFT0','VLSR41')
tbl.rename_column('e_SHIFT0','eVLSR41')
tbl.add_column(table.Column(data=tbl['WIDTH0']*np.sqrt(8*np.log(2)), name='VWFHM41'))
tbl.add_column(table.Column(data=tbl['e_WIDTH0']*np.sqrt(8*np.log(2)), name='eVWFHM41'))
depree_merged = table.join(depree, tbl, keys='Source', join_type='outer').filled(np.nan)
depree_merged.write(paths.tpath('DePree2011_plus_H41afits.ipac'), format='ascii.ipac')


ok41 = depree_merged['AMPLITUDE0'] > depree_merged['e_AMPLITUDE0']*3
std66 = np.nanstd(depree_merged['VLSR66'][ok41])
std52 = np.nanstd(depree_merged['VLSR52'][ok41])
std41 = np.nanstd(depree_merged['VLSR41'][ok41])
print("1D Velocity Dispersion of 41a: {0}, 52a: {1}, 66a: {2}".format(std41, std52, std66))
vlsrmean66 = np.nanmean(depree_merged['VLSR66'][ok41])
vlsrmean52 = np.nanmean(depree_merged['VLSR52'][ok41])
vlsrmean41 = np.nanmean(depree_merged['VLSR41'][ok41])
print("Mean VLSR of 41a: {0}, 52a: {1}, 66a: {2}".format(vlsrmean41, vlsrmean52, vlsrmean66))

# some plots
fig = pl.figure(4)
fig.clf()
ax = fig.gca()
sc = ax.scatter(depree_merged['RA'], depree_merged['Dec'], c=depree_merged['VLSR41'], vmin=35, vmax=80, s=150)
pl.colorbar(mappable=sc)

ra, dec = depree_merged['RA'], depree_merged['Dec']
dist = coordinates.SkyCoord(ra, dec, frame='fk5').separation(center)

fig = pl.figure(2)
fig.clf()
ax = fig.gca()
ax.plot((dist*8500*u.pc).to(u.pc, u.dimensionless_angles())[ok41], depree_merged['VLSR41'][ok41], 'o')
ax.set_xlabel("Distance from 'center' (pc)")
ax.set_ylabel("$V_{\mathrm{LSR}}(\mathrm{H}41\\alpha)$ [km s$^{-1}$]")

fig = pl.figure(1)
fig.clf()
pl.hist(depree_merged['VLSR41'][ok41])
