import regions
import paths
from gaussfit_catalog import gaussfit_catalog, gaussfit_image, data_to_table

if __name__ == "__main__":
    regs = regions.read_ds9(paths.rpath('cores_with_names.reg'))

    from files import contfilename as contfnpath

    fit_data = gaussfit_catalog(contfnpath, regs, savepath=paths.fpath('gaussfits'))

    tbl = data_to_table(fit_data)

    tbl.rename_column("chi2/n", "chi2_n")
    tbl.write(paths.tpath("gaussian_fit_table.ipac"), format='ascii.ipac',
              overwrite=True)
