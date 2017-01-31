import pylab as pl
from astropy.table import Table

import paths

lines = ["CFp", "CH3OH7m26-716", "H15NC", "H2CO615-616", "H2CS303-202",
         "H2CS313-212", "H2CS322-221", "H41a", "HC3N", "HCN", "HNC",]

cont_tbl = Table.read(paths.tpath('continuum_photometry.ipac'), format='ascii.ipac')
line_tbl = Table.read(paths.tpath('line_photometry.tab'),
                      format='ascii.fixed_width', delimiter='|')

pl.figure(1).clf()

for ii,linename in enumerate(lines):

    ax = pl.subplot(4,3,ii+1)

    mask = line_tbl['{0}_peak'.format(linename)] / line_tbl['{0}_bgrms'.format(linename)] > 3

    ax.loglog(cont_tbl['peak'][mask], line_tbl['{0}_peak'.format(linename)][mask], '.')
    ax.set_title(linename)
