import numpy as np
from imf import imf
import paths
from astropy.table import Table

core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')

kroupa = imf.Kroupa()

over8fraction = kroupa.integrate(8, 200)[0]
over8median = 2**(1/(1.7-1)) * 8

x = np.linspace(8,200,50000)
y = kroupa(x)
over8mean = (x*y).sum()/y.sum()

nsources = len(core_phot_tbl)

print("Mass estimate if all sources are 8 msun = {0}".format(nsources*8/over8fraction))
print("Mass estimate if Mbar={1} = {0}".format(nsources*over8mean/over8fraction, over8mean))
