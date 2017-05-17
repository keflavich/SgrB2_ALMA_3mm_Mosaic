from astropy.table import Table

tbl = Table.read('../tables/continuum_photometry.ipac', format='ascii.ipac')

with open('ncores.tex', 'w') as fh:
    fh.write("\\newcommand{{\\ncores}}{{{0}\\xspace}}".format(len(tbl)))
