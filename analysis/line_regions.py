from astropy import coordinates
from astropy import units as u
from astropy.table import Table

import paths

lines = ["CFp", "CH3OH7m26-716", "H15NC", "H2CO615-616", "H2CS303-202",
         "H2CS313-212", "H2CS322-221", "H41a", "HC3N", "HCN", "HNC",]

line_tbl = Table.read(paths.tpath('line_photometry.tab'),
                      format='ascii.fixed_width', delimiter='|')

for ii,linename in enumerate(lines):

    mask = (line_tbl['{0}_peak'.format(linename)] /
            line_tbl['{0}_bgrms'.format(linename)] > 4)

    print("For line {0}, found {1} detections out of {2} ({3:0.1f}%)"
          .format(linename, mask.sum(), mask.size, mask.sum()/mask.size*100))

    coords = coordinates.SkyCoord(line_tbl['{0}_RA'.format(linename)][mask],
                                  line_tbl['{0}_Dec'.format(linename)][mask],
                                  unit=(u.deg, u.deg),
                                  frame='fk5',)

    with open(paths.rpath('cores_{0}.reg'.format(linename)), 'w') as fh:
        fh.write('fk5\n')
        for crd in coords:
            fh.write("point({0},{1}) # point=x\n".format(crd.ra.deg, crd.dec.deg))
