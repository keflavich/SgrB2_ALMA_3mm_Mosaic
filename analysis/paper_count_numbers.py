import numpy as np
import paths
from astropy.table import Table, Column
from astropy import units as u
import imf

core_phot_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD_andclusters.ipac"), format='ascii.ipac')
hii = core_phot_tbl['SIMBAD_OTYPE'] == 'HII'

# Some numbers quoted in the paper
nfaint = ((u.Quantity(core_phot_tbl['peak'], u.Jy/u.beam) < 1.5*u.mJy/u.beam)).sum()
print("number of faint sources < 1.5 mJy: {0}".format(nfaint))
# should be 100

# number not hii:
n_not_hii_gt_1pt5 = ((~hii) & (u.Quantity(core_phot_tbl['peak'], u.Jy/u.beam) > 1.5*u.mJy/u.beam)).sum()
# should be 140
print("Number of not-HII regions at flux >1.5 mJy: {0}".format(n_not_hii_gt_1pt5))
n_not_hii_mid = ((~hii) & (u.Quantity(core_phot_tbl['peak'], u.Jy/u.beam) > 1.5*u.mJy/u.beam) & ((u.Quantity(core_phot_tbl['peak'], u.Jy/u.beam) < 10*u.mJy/u.beam))).sum()
# should be 119
print("Number of not-HII regions at 10 mJy > flux >1.5 mJy: {0}".format(n_not_hii_mid))
nbright = (u.Quantity(core_phot_tbl['peak'], u.Jy/u.beam) > 10*u.mJy/u.beam).sum()
print("nbright = {0}".format(nbright))
nbright_hii = (hii & (u.Quantity(core_phot_tbl['peak'], u.Jy/u.beam) > 10*u.mJy/u.beam)).sum()
print("nbright,hii = {0}".format(nbright_hii))


# in order to produce a >10 mJy flux, a star must produce > 2e47 photons, which
# implies a star B1.5V or earlier, >10 Msun.
kroupa = imf.Kroupa()

mmax = 200
cutoff1 = 8
cutoff2 = 10
x = np.linspace(cutoff2,mmax,50000)
y = kroupa(x)
over10mean = (x*y).sum()/y.sum()
x = np.linspace(cutoff1,cutoff2,50000)
y = kroupa(x)
eighttotenmean = (x*y).sum()/y.sum()

over8fraction = (kroupa.m_integrate(cutoff1, mmax)[0] /
                 kroupa.m_integrate(kroupa.mmin, mmax)[0])
over10fraction = (kroupa.m_integrate(cutoff2, mmax)[0] /
                  kroupa.m_integrate(kroupa.mmin, mmax)[0])
mass8to10fraction = over8fraction - over10fraction
print("Mass 8-10 msun fraction: {0}".format(mass8to10fraction))
predicted_mass_from_gt10 = nbright * over10mean / over10fraction
print("Mean mass of >10 Msun stars and total massive mass predicted assuming all are average stars = {0}, {1}".format(over10mean, nbright*over10mean))
print("Expected mass of stars from the 50 massive ones: {0}".format(predicted_mass_from_gt10))
print("Expected mass of 8-10 Msun stars: {0}".format(predicted_mass_from_gt10 * mass8to10fraction))
print("Expected # of 8-10 msun stars: {0}".format(predicted_mass_from_gt10 * mass8to10fraction / eighttotenmean))
