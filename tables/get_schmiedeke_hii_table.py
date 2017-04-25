import requests
import paths
from bs4 import BeautifulSoup
from astropy.table import Table,Column
from astropy import coordinates
from astropy import units as u
import regions

rslt = requests.get('http://www.aanda.org/articles/aa/full_html/2016/04/aa27311-15/T4.html')
soup = BeautifulSoup(rslt.text, 'html5lib')
htmltable = soup.findAll('table')[-1]

allrows = (htmltable.findAll('tr'))

colnames = [x.text.strip("0123456789\n \t") for x in allrows[0].findAll('td')]
units = [x.text.strip() for x in allrows[1].findAll('td')]

data = {colname: [] for colname in colnames}

for row in allrows[3:]:

    for colname, val in zip(colnames, row.findAll('td')):

        data[colname].append(val.text)

for ii,xx in enumerate(data['Dec']):
    data['Dec'][ii] = xx.replace('–','-')


tbl = Table(data=[Column(data=data[colname], name=colname, unit=None or unit)
                  for colname, unit in zip(colnames, units)]
           )
tbl.rename_column('log()','log(Q_lyc)')

tbl.write(paths.tpath("Schmiedeke2016_HIIregions_tableB1.txt"),
          format='ascii.fixed_width', overwrite=True)


reglist = [regions.CircleSkyRegion(coordinates.SkyCoord(row['RA'], row['Dec'], frame='fk5', unit=(u.hour, u.deg)),
                                   radius=(u.Quantity(float(row['robs'])*1000, u.au) / (8.5*u.kpc)).to(u.arcsec, u.dimensionless_angles()),
                                   meta={'name': row['ID']},
                                   visual={'name': row['ID']},
                                  )
           for row in tbl
           if row['ID']
          ]

regions.write_ds9(reglist, paths.rpath('Schmiedeke2016_HIIregions_tableB1.reg'))
