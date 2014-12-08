import numpy as np
import pvextractor
from astropy import units as u
from spectral_cube import SpectralCube
import aplpy
import pylab as pl
import glob

pl.ioff()

paths = pvextractor.pvregions.paths_from_regfile('loop_segment_30kms.reg')

for mol in ("HCN","HNC","HCOp","HC3N","H2CS303-202","H2CS303202"):
    for pathno,path in enumerate(paths):
        for fn in (glob.glob("*{0}*image.pbcor.fits".format(mol)) +
                   glob.glob("feathered/Feathered*{0}*.fits".format(mol)) +
                   glob.glob("feathered/Regridded*{0}*fits".format(mol))):

            if any(x in fn for x in ('mask','max','sum','mom','pvext')):
                continue

            rcube = (SpectralCube.read(fn)
                     .with_spectral_unit(u.km/u.s, velocity_convention='radio')
                     .spectral_slab(0.0*u.km/u.s, 160*u.km/u.s))

            slice = pvextractor.extract_pv_slice(rcube, path)
            slice.writeto(fn.replace(".fits",".pvextraction_path{0}.fits".format(pathno)),
                          clobber=True)

            F = aplpy.FITSFigure(slice)
            F.show_grayscale()
            F.save(fn.replace(".fits",".pvextraction_path{0}.png".format(pathno)))
            print fn.replace(".fits",".pvextraction_path{0}.png".format(pathno))
            F.show_contour('SgrB2_a_03_7M.HNC.image.pbcor.pvextraction_path0.fits')
            F.save(fn.replace(".fits",".pvextraction_path{0}_contourHNCACA.png".format(pathno)))
            F.remove_layer('contour_set_1')
            F.show_contour('feathered/Feathered_HNC.pvextraction_path0.fits')
            F.save(fn.replace(".fits",".pvextraction_path{0}_contourHNCFeath.png".format(pathno)))

            for ii in pl.get_fignums():
                pl.close(ii)

F = aplpy.FITSFigure('SgrB2_a_03_7M.HNC.image.pbcor.max.fits')
F.show_grayscale()
#F.show_lines([np.array([path._coords.ra.value, path._coords.dec.value])])
xy = np.array(path.sample_points(1, wcs=rcube.wcs))
F._ax1.plot(xy[0,:], xy[1,:], color='g', linewidth=5, alpha=0.3, zorder=10)
F._ax1.plot(xy[0,:], xy[1,:], color='g', linewidth=2, alpha=0.5, zorder=100)
pl.draw()
F.save('path_on_HNCmax.png')
