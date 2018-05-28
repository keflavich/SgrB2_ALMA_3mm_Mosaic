import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u, table, coordinates
import pyspeckit
import pyregion
import paths
import pylab as pl

from latex_info import (latexdict, format_float, round_to_n,
                        strip_trailing_zeros, exp_to_tex)
latexdict = latexdict.copy()

regs = pyregion.open(paths.rpath('SgrB2_1.3cm_hiiRegions_shiftedtoALMA.reg'))

cube = SpectralCube.read(paths.Fpath('M_cutouts/SgrB2_b3_12M_TE.H41a.image.pbcor_M_cutout.medsub.fits'))

if 'fit_values' not in locals():
    fit_values = {}
    coords = {}

    for reg in regs:
        if 'text' not in reg.attr[1]:
            continue
        name = reg.attr[1]['text'].strip("{}")

        try:
            cutout = cube.subcube_from_ds9region(pyregion.ShapeList([reg])).to(u.K,
                                                                               cube.beam.jtok_equiv(cube.with_spectral_unit(u.GHz).spectral_axis))
            print(name)
        except ValueError as ex:
            print("{1} - An error was raised: {0}".format(ex, name))
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

snr = tbl['AMPLITUDE0'] / tbl['e_AMPLITUDE0']
tbl.add_column(table.Column(name='SNR41', data=snr))
ok = snr > 3

tbl[ok].write(paths.tpath('H41a_fits.ipac'), format='ascii.ipac', overwrite=True)
for cn in [xx.name for xx in columns]+[xx.name for xx in errcolumns]:
    tbl[~ok][cn] = np.nan
    tbl[cn][~ok] = np.nan

print(tbl['Name','SNR41'][~ok])

center = coordinates.SkyCoord('17:47:20.174', '-28:23:04.233', frame='fk5', unit=(u.hour, u.deg))

depree = table.Table.read(paths.tpath('DePree2011.txt'), format='ascii.fixed_width', delimiter='|')

tbl.rename_column('Name','Source')
tbl.rename_column('SHIFT0','VLSR41')
tbl.rename_column('e_SHIFT0','eVLSR41')
tbl['VLSR41'].unit = u.km/u.s
tbl['eVLSR41'].unit = u.km/u.s
tbl.add_column(table.Column(data=tbl['WIDTH0']*np.sqrt(8*np.log(2)), unit=u.km/u.s, name='VFWHM41'))
tbl.add_column(table.Column(data=tbl['e_WIDTH0']*np.sqrt(8*np.log(2)), unit=u.km/u.s, name='eVFWHM41'))
tbl.remove_column('WIDTH0')
tbl.remove_column('e_WIDTH0')
depree_merged = table.join(depree, tbl, keys='Source', join_type='outer').filled(np.nan)
depree_merged.write(paths.tpath('DePree2011_plus_H41afits.ipac'), format='ascii.ipac')

depree_merged_tex = depree_merged.copy()

latexdict['header_start'] = '\label{tab:h41afits}'#\n\\footnotesize'
latexdict['preamble'] = '\caption{H41$\\alpha$ Line Fits}\n'
latexdict['col_align'] = 'l'*len(depree_merged.columns)
latexdict['tabletype'] = 'table'
latexdict['tablefoot'] = ("\par\nFits to the H41$\\alpha$ line from \citet{Ginsburg2018a}.")

for key in depree_merged_tex.colnames:
    newname = key
    if key[0].lower() == 'e':
        del depree_merged_tex[key]
        continue
    #if 'e_' in key:
    #    newname = key.replace("_","")
    if 'LSR' in key:
        newname = key.split("LSR")[0] + "$_\mathrm{LSR}$(" + key.split("LSR")[1] +")"
    if 'FWHM' in key:
        newname = "$\mathrm{FWHM}$(" + key.split("FWHM")[1] + ")"
    if 'AMPLITUDE' in key:
        #newname = key.replace("AMPLITUDE","Peak")
        del depree_merged_tex[key]
        continue
    depree_merged_tex.rename_column(key, newname)

merged_tex_ok = (np.isfinite(depree_merged_tex['RA']) &
                 (np.isfinite(depree_merged_tex['V$_\\mathrm{LSR}$(41)']) | 
                  np.isfinite(depree_merged_tex['V$_\\mathrm{LSR}$(52)']) | 
                  np.isfinite(depree_merged_tex['V$_\\mathrm{LSR}$(66)']))
                )
depree_merged_tex = depree_merged_tex[merged_tex_ok]

formats = {key: lambda x: (strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2)))
                           if not np.isnan(x) else '-')
           for key in depree_merged_tex.colnames}
del formats['Source']
del formats['RA']
del formats['Dec']

coords_column = coordinates.SkyCoord(depree_merged_tex['RA'],
                                     depree_merged_tex['Dec'], frame='fk5',
                                     unit=(u.deg, u.deg))
depree_merged_tex.remove_column('RA')
depree_merged_tex.remove_column('Dec')
depree_merged_tex.add_column(table.Column(name="Coordinates", data=coords_column))
def cformat(x):
    try:
        return x.to_string('hmsdms', sep=":")
    except:
        return '-'
formats['Coordinates'] = cformat

column_order = ['Source',
 'Coordinates',
 'V$_\\mathrm{LSR}$(41)',
 '$\\mathrm{FWHM}$(41)',
 'V$_\\mathrm{LSR}$(52)',
 '$\\mathrm{FWHM}$(52)',
 'V$_\\mathrm{LSR}$(66)',
 '$\\mathrm{FWHM}$(66)',
               ]
depree_merged_tex = depree_merged_tex[column_order]

depree_merged_tex.write(paths.cfepath('h41afits.tex'), format='ascii.latex',
                        overwrite=True, latexdict=latexdict, formats=formats,)


ok41 = depree_merged['AMPLITUDE0'] > depree_merged['e_AMPLITUDE0']*3
std66 = np.nanstd(depree_merged['VLSR66'][ok41])
std52 = np.nanstd(depree_merged['VLSR52'][ok41])
std41 = np.nanstd(depree_merged['VLSR41'][ok41])
print("1D Velocity Dispersion of 41a: {0}, 52a: {1}, 66a: {2}".format(std41, std52, std66))
vlsrmean66 = np.nanmean(depree_merged['VLSR66'][ok41])
vlsrmean52 = np.nanmean(depree_merged['VLSR52'][ok41])
vlsrmean41 = np.nanmean(depree_merged['VLSR41'][ok41])
print("Mean VLSR of 41a: {0}, 52a: {1}, 66a: {2}".format(vlsrmean41, vlsrmean52, vlsrmean66))

print()
print("Same, excluding the 40 km/s group:")
exclude_low = depree_merged['VLSR41'] > 45
std66 = np.nanstd(depree_merged['VLSR66'][ok41&exclude_low])
std52 = np.nanstd(depree_merged['VLSR52'][ok41&exclude_low])
std41 = np.nanstd(depree_merged['VLSR41'][ok41&exclude_low])
print("1D Velocity Dispersion of 41a: {0}, 52a: {1}, 66a: {2}".format(std41, std52, std66))
vlsrmean66 = np.nanmean(depree_merged['VLSR66'][ok41&exclude_low])
vlsrmean52 = np.nanmean(depree_merged['VLSR52'][ok41&exclude_low])
vlsrmean41 = np.nanmean(depree_merged['VLSR41'][ok41&exclude_low])
print("Mean VLSR of 41a: {0}, 52a: {1}, 66a: {2}".format(vlsrmean41, vlsrmean52, vlsrmean66))

# some plots
fig = pl.figure(4)
fig.clf()
ax = fig.gca()
sc = ax.scatter(depree_merged['RA'], depree_merged['Dec'],
                c=depree_merged['VLSR41'], vmin=35, vmax=80, s=150)
cb = pl.colorbar(mappable=sc)
cb.set_label("VLSR(H41$\\alpha$)")
pl.title("Location of sources colored by VLSR")

ra, dec = depree_merged['RA'], depree_merged['Dec']
dist = coordinates.SkyCoord(ra, dec, frame='fk5').separation(center)

fig = pl.figure(2)
fig.clf()
ax = fig.gca()
ax.plot((dist*8500*u.pc).to(u.pc, u.dimensionless_angles())[ok41],
        depree_merged['VLSR41'][ok41], 'o')
ax.set_xlabel("Distance from 'center' (pc)")
ax.set_ylabel("$V_{\mathrm{LSR}}(\mathrm{H}41\\alpha)$ [km s$^{-1}$]")

fig = pl.figure(3)
fig.clf()
pl.hist(depree_merged['VLSR41'][ok41])
pl.xlabel("VLSR(H41$\\alpha$)")
