import regions
import paths
from gaussfit_catalog import gaussfit_catalog, gaussfit_image
from astropy.table import Table, Column

def data_to_table(fit_data):
    names = fit_data.keys()
    numnames = [int(nm[5:].split("_")[0]) for nm in names]
    stripnames = [nm[5:] for nm in names]
    stripnames = [fullname for nnm,fullname in sorted(zip(numnames,stripnames))]
    names = [fullname for nnm,fullname in sorted(zip(numnames,names))]
    namecol = Column(name='Name', data=stripnames)
    colnames = ['amplitude', 'center_x', 'center_y', 'fwhm_x', 'fwhm_y', 'pa',
                'chi2', 'chi2/n', 'e_amplitude', 'e_center_x', 'e_center_y',
                'e_fwhm_x', 'e_fwhm_y', 'e_pa', 'success',]
    columns = [Column(name=k, data=[fit_data[entry][k].value
                                    if hasattr(fit_data[entry][k],'value')
                                    else fit_data[entry][k]
                                    for entry in names],
                      unit=(fit_data[names[0]][k].unit
                            if hasattr(fit_data[names[0]][k], 'unit')
                            else None))
               for k in colnames]

    return Table([namecol]+columns)

if __name__ == "__main__":
    regs = regions.read_ds9(paths.rpath('cores_with_names.reg'))

    from files import contfilename as contfnpath

    fit_data = gaussfit_catalog(contfnpath, regs, savepath=paths.fpath('gaussfits'))

    tbl = data_to_table(fit_data)

    tbl.rename_column("chi2/n", "chi2_n")
    tbl.write(paths.tpath("gaussian_fit_table.ipac"), format='ascii.ipac',
              overwrite=True)
