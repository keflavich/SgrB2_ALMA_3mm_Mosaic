import numpy as np
from astropy.table import Table

tbl = Table.read('../tables/continuum_photometry_withSIMBAD_andclusters.ipac', format='ascii.ipac')

OK = tbl['peak'] != 0
nhii = sum([('HII' in x) and k for x,k in zip(tbl['Classification'], OK)])
xray = np.array(['X' in row['Classification'] for row in tbl], dtype='bool')
maser = np.array(['M' in row['Classification'] for row in tbl], dtype='bool')
ncores = OK.sum()


with open('ncores.tex', 'w') as fh:
    fh.write("\\newcommand{{\\ncores}}{{{0}\\xspace}}\n".format(ncores))
    fh.write("\\newcommand{{\\nhii}}{{{0}\\xspace}}\n".format(nhii))
    fh.write("\\newcommand{{\\nmasermatch}}{{{0}\\xspace}}\n".format(maser.sum()))
    fh.write("\\newcommand{{\\nxraymatch}}{{{0}\\xspace}}\n".format(xray.sum()))
