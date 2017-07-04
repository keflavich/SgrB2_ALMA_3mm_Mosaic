import numpy as np
import os
import glob
import pyspeckit
from astropy import units as u
from spectral_cube import SpectralCube
import paths

from restfreqs import restfreqs

import pylab as pl
pl.close('all')



velo = 62.5*u.km/u.s

species_names = list(restfreqs.keys())
frequencies = u.Quantity([restfreqs[key] for key in species_names],
                         unit=u.GHz)

for fn in glob.glob(paths.Fpath('tp/tp_concat*fits')):
    basefn = os.path.basename(fn)
    outf = paths.Fpath('tp/avgspectra/avg_{0}'.format(basefn))
    if not os.path.exists(outf):
        cube = SpectralCube.read(fn)
        cube.allow_huge_operations=True
        mad = cube.mad_std(axis=0)
        sum = (cube*mad).sum(axis=(1,2))
        mean = sum/np.nansum(mad)
        mean.write(paths.Fpath('tp/avgspectra/weighted_avg_{0}'.format(basefn)), overwrite=True)
        mn = cube.mean(axis=(1,2))
        mn.write(outf, overwrite=True)

for fn in glob.glob(paths.Fpath('tp/tp_concat*fits')):
    basefn = os.path.basename(fn)
    #spw = pyspeckit.Spectrum('avspec/weighted_avg_{0}'.format(fn))
    sp = pyspeckit.Spectrum(paths.Fpath('tp/avgspectra/avg_{0}'.format(basefn)))
    sp.plotter()
    #spw.plotter(axis=sp.plotter.axis, clear=False, color='b')
    sp.plotter.axis.set_ylim(-0.5, 4.0)
    sp.plotter.savefig(paths.fpath('tpavgspectra/avg_{0}'
                       .format(basefn.replace(".fits",".png"))))


    plot_kwargs = {'color':'r', 'linestyle':'--'}
    annotate_kwargs = {'color': 'r'}
    sp.plotter.line_ids(species_names, u.Quantity(frequencies),
                        velocity_offset=velo, plot_kwargs=plot_kwargs,
                        annotate_kwargs=annotate_kwargs)

    sp.plotter.savefig(paths.fpath('tpavgspectra/lineid_avg_{0}'
                                   .format(basefn.replace(".fits",".png"))))
