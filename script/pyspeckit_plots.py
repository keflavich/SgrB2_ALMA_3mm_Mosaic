from spectral_cube import SpectralCube, BooleanArrayMask
from astropy import units as u
import numpy as np
import matplotlib as mpl
import pyspeckit
import os
import pylab as pl

molecules = ['CFp', 'CH3CN_5-4_3', 'H2CO615-616', 'H2CS', 'H2CS303-202',
             'H2CS321-220', 'H41a', 'HC3N', 'HCN', 'HCOp', 'HNC',]


cubes = {mol: SpectralCube.read('SgrB2_a_03_7M.{0}.image.pbcor.fits'.format(mol))
                          .with_spectral_unit(u.km/u.s,
                                              velocity_convention='radio')
         for mol in molecules}

pcubes = {mol: pyspeckit.Cube(cube=cubes[mol]) for mol in molecules}


fig1 = pl.figure(1)
fig1.clf()

for pcube in pcubes.values():
    pcube.plot_spectrum(101,90, clear=False, figure=fig1)
