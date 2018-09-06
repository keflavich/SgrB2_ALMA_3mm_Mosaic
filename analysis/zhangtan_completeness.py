import os
import glob
import numpy as np
from astropy import table
from astropy import units as u
from astropy import convolution
from astropy.utils.console import ProgressBar

# parse ZT2018 models

datapath = '/Volumes/external/zhang_sedmodels/yczhang04-sedfit-959e99c/Model_SEDs/model_info/'
sedpath = '/Volumes/external/zhang_sedmodels/yczhang04-sedfit-959e99c/Model_SEDs/seds/'

fitstable = '/Volumes/external/zhang_sedmodels/model_info_table.fits'

if os.path.exists(fitstable):
    tbl = table.Table.read(fitstable)
else:
    alldata = []

    for fn in ProgressBar(glob.glob('{0}/*.dat'.format(datapath))):
        x,y,z = map(int, os.path.split(fn)[-1][:8].split("_"))
        with open(fn, 'r') as fh:
            lines = fh.readlines()
            headings = lines[0].split(',')
            data = np.array(lines[1].split(), dtype='float')

        sed_data_fn = fn.replace("model_info", "sed").replace(".dat","_{0:02d}.dat")
        sed_vals = []
        for viewangle in range(1,21):
            with open(sed_data_fn.format(viewangle), 'r') as fh:
                wave, flux = np.loadtxt(fh).T
                flux[flux == 0] = np.nan
                fluxcorr = convolution.convolve(flux, convolution.Gaussian1DKernel(3))
                # 3mm = 3e3 microns
                lum_3mm = fluxcorr[np.nanargmin(np.abs((wave-3e3)))]
                sed_vals.append(lum_3mm)


        alldata.append(np.array([x,y,z] + data.tolist() + sed_vals))

    unit_mappings = {'msun': 'Msun', 'rsun': 'Rsun', 'lsun': 'Lsun',
                     'msun/yr': 'Msun/yr'}
    def remap_units(x):
        if x in unit_mappings:
            return unit_mappings[x]
        else:
            return x


    units=[None, None, None] + [u.Unit(remap_units(x.split()[1].strip("()"))) for x in headings] + [u.L_sun]*20

    datacols = [table.Column(data=xx, unit=unit)
                for xx, unit in zip(np.array(alldata).T, units)]


    tbl = table.Table(data=datacols,
                      names=['xx', 'yy', 'zz'] + [x.split()[0] for x in headings] +
                      ["th{0:0.3f}".format(th) for th in
                       np.arccos(np.linspace(0.975, 0.025, 20))*180/np.pi],
                     )

    tbl.write(fitstable, overwrite=True)

# luminosity factor:
lfac = ((1*u.L_sun) / (4*np.pi*(8.4*u.kpc)**2) / (3.3*u.mm).to(u.GHz, u.spectral())).to(u.mJy)

flux_headers = tbl.colnames[-20:]

fluxes = np.array(tbl[flux_headers].as_array().tolist()) * lfac.value

import pylab as pl
pl.figure(1).clf()
pl.xlabel("$L_{tot}$ [L$_{\odot}$]")
pl.ylabel(r"$S_{3 \mathrm{mm}}$ [mJy]")
for fh in flux_headers:
    pl.loglog(tbl['ltot'], tbl[fh]*lfac, 'o', alpha=0.25)


pl.figure(2).clf()
pl.xlabel("$L_{*}$ [L$_{\odot}$]")
pl.ylabel(r"$S_{3\mathrm{mm}}$ [mJy]")
for fh in flux_headers:
    pl.loglog(tbl['lstar'], tbl[fh]*lfac, 'o', alpha=0.25)



pl.figure(3).clf()
pl.xlabel("$M_{star}$ [M$_{\odot}$]")
pl.ylabel(r"$S_{3 \mathrm{mm}}$ [mJy]")
for fh in flux_headers:
    pl.loglog(tbl['mstar'], tbl[fh]*lfac, 'o', alpha=0.25)
