import numpy as np
from astropy import units as u
from astropy import log
import paths
import pyregion
from astropy.io import fits
from astropy.table import Table,Column
import masscalc
import radio_beam

regions = pyregion.open(paths.rpath('sgrb2_cores_TE.reg'))

contfnpath = paths.tmpath('te/SgrB2_selfcal_full_TE_selfcal4_ampphase.image.pbcor.fits')
contfile = fits.open(contfnpath)
data = contfile[0].data
beam = radio_beam.Beam.from_fits_header(contfnpath)

results = {}
units = {'peak':u.Jy/u.beam,
         'sum':u.Jy/u.beam,
         'npix':u.dimensionless_unscaled,
         'beam_area':u.sr,
         'peak_mass':u.M_sun,
         'peak_col':u.cm**-2,
         'RA': u.deg,
         'Dec': u.deg,
        }

for ii,reg in enumerate(regions):
    if 'text' not in reg.attr[1]:
        name = str(ii)
    else:
        name = reg.attr[1]['text']

    # all regions are points: convert them to 0.5" circles
    reg.name = 'circle'
    reg.coord_list.append(0.5/3600.)
    reg.params.append(pyregion.region_numbers.AngularDistance('0.5"'))

    shreg = pyregion.ShapeList([reg])
    log.info(name)

    mask = shreg.get_mask(hdu=contfile[0])

    results[name] = {'peak': data[mask].max(),
                     'sum': data[mask].sum(),
                     'npix': mask.sum(),
                     'beam_area': beam.sr,
                     'RA': reg.coord_list[0],
                     'Dec': reg.coord_list[1],
                    }
    results[name]['peak_mass'] = masscalc.mass_conversion_factor()*results[name]['peak']*u.M_sun
    results[name]['peak_col'] = masscalc.col_conversion_factor()*results[name]['peak']*u.cm**-2

# invert the table to make it parseable by astropy...
# (this shouldn't be necessary....)
results_inv = {'name':{}}
columns = {'name':[]}
for k,v in results.items():
    results_inv['name'][k] = k
    columns['name'].append(k)
    for kk,vv in v.items():
        if kk in results_inv:
            results_inv[kk][k] = vv
            columns[kk].append(vv)
        else:
            results_inv[kk] = {k:vv}
            columns[kk] = [vv]

for c in columns:
    if c in units:
        columns[c] = columns[c] * units[c]

tbl = Table([Column(data=columns[k],
                    name=k)
             for k in ['name','RA','Dec','peak','sum','npix','beam_area','peak_mass','peak_col']])

tbl.sort('peak_mass')
tbl.write(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')
