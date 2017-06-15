from astropy.table import Table

tbl = Table.read('../tables/continuum_photometry_withSIMBAD_andclusters.ipac', format='ascii.ipac')

nhii = sum(['HII' in x for x in tbl['Classification']])

with open('ncores.tex', 'w') as fh:
    fh.write("\\newcommand{{\\ncores}}{{{0}\\xspace}}\n".format(len(tbl)))
    fh.write("\\newcommand{{\\nhii}}{{{0}\\xspace}}\n".format(nhii))
