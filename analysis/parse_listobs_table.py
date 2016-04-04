from astropy.table import Table
import re

path='script_12m_te/listobs/field_table.tsv'

def parse_table(path):
    data = []

    with open(path,'r') as fh:
        lines = fh.readlines()
        names = lines[0].split()
        for row in lines[1:]:
            data.append({name:val for name, val in zip(names, row.split())})

    tbl = Table(data)

    return tbl

def tbl_to_reg(tbl, outname, radius=None):
    dotstocolons = re.compile("(-?[0-9]?[0-9]?)\.([0-9]?[0-9]?)\.([0-9]?[0-9]?)\.([0-9]*)")
    with open(outname, 'w') as fh:
        fh.write('fk5\n')
        for row in tbl:
            rowdict = dict(zip(row.colnames, row.data))
            rowdict['Decl'] = "{0}:{1}:{2}.{3}".format(*dotstocolons.search(rowdict['Decl']).groups())
            if radius is None:
                fh.write("point({RA}, {Decl}) # text={{{ID}}}\n".format(**rowdict))
            else:
                fh.write("circle({RA}, {Decl}, {radius}) # text={{{ID}}}\n".format(radius=radius, **rowdict))
