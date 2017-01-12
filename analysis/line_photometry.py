import numpy as np

import regions
import radio_beam
from astropy import units as u
from astropy import log
from astropy import wcs
from astropy.stats import mad_std
from astropy.io import fits
from astropy.table import Table,Column

import paths
from core_photometry import photometry

lines = ["CFp", "CH3OH7m26-716", "H15NC", "H2CO615-616", "H2CS303-202",
         "H2CS313-212", "H2CS322-221", "H41a", "HC3N", "HCN", "HNC",]

if __name__ == "__main__":
    regs = regions.read_ds9(paths.rpath('sgrb2_cores_TE.reg'))

    all_results = {}

    for line in lines:
        fn = paths.mergepath("max/SgrB2_b3_7M_12M.{0}.image.pbcor_max_medsub.fits".format(line))
        ffile = fits.open(fn)
        data = ffile[0].data
        header = ffile[0].header
        beam = radio_beam.Beam.from_fits_header(header)
        mywcs = wcs.WCS(header)

        units = {'peak':u.Jy/u.beam,
                 'sum':u.Jy/u.beam,
                 'npix':u.dimensionless_unscaled,
                 'beam_area':u.sr,
                 'peak_brightness': u.K,
                 'freq': u.GHz,
                 'RA': u.deg,
                 'Dec': u.deg,
                }

        results = photometry(data, mywcs, regs, beam)

        for name in results:
            rslt_dict = {line+"_"+key: value
                         for key,value in results[name].items()}
            if name in all_results:
                all_results[name].update(rslt_dict)
            else:
                all_results[name] = rslt_dict

            restfrq = header['RESTFRQ'] * u.Hz
            tbmax = (results[name]['peak']*u.Jy).to(u.K,
                                                    u.brightness_temperature(results[name]['beam_area'],
                                                                             restfrq))
            all_results[name][line+'_peak_brightness'] = tbmax

            all_results[name][line+"_freq"] = restfrq

    # invert the table to make it parseable by astropy...
    # (this shouldn't be necessary....)
    results_inv = {'name':{}}
    columns = {'name':[]}
    colnames = {'name'} # set
    for k,v in all_results.items():
        results_inv['name'][k] = k
        columns['name'].append(k)
        for kk,vv in v.items():
            colnames.add(kk)
            if kk in results_inv:
                results_inv[kk][k] = vv
                columns[kk].append(vv)
            else:
                results_inv[kk] = {k:vv}
                columns[kk] = [vv]

    for c in columns:
        umatch = [k for k in units if k in c]
        if len(umatch) >= 1:
            cn = max(umatch, key=len)
            columns[c] = u.Quantity(columns[c], units[cn])

    tbl = Table([Column(data=columns[k],
                        name=k)
                 for k in sorted(colnames)])

    tbl.write(paths.tpath("line_photometry.csv"), format='ascii.csv',
              overwrite=True)
    tbl.write(paths.tpath("line_photometry.tab"), format='ascii.fixed_width',
              delimiter='|',
              overwrite=True)
    tbl.write(paths.tpath("line_photometry.ipac"), format='ascii.ipac',
              overwrite=True)
