import numpy as np
from imf import imf
import paths
from astropy.table import Table

core_phot_tbl = Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac')

core_powerlaw_index = 1.94

kroupa = imf.Kroupa()

o_mmin = 8
mmax = 200
over8fraction = (kroupa.m_integrate(o_mmin, mmax)[0] /
                 kroupa.m_integrate(kroupa.mmin, mmax)[0])
# from wikipedia, median of power-law is well-defined for more slopes than the
# mean, but we're interested in the mean so I just compute that numerically
# below
over8median = 2**(1/(1.94-1)) * o_mmin

x = np.linspace(o_mmin,mmax,50000)
y = kroupa(x)
over8mean = (x*y).sum()/y.sum()

nsources = len(core_phot_tbl)

print("Mass fraction M>8 = {0}".format(over8fraction))
print("Mass of observed sources, assuming all are 8 msun = {0}".format(nsources*8))
print("Total Mass estimate if all sources are 8 msun = {0}".format(nsources*8/over8fraction))
print("Total Mass estimate if Mbar={1} = {0}".format(nsources*over8mean/over8fraction, over8mean))
