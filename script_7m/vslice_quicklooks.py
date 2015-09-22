from astropy import units as u
from spectral_cube import SpectralCube
import pylab as pl
import glob
import os

pl.ioff()

def pdir(fn):
    pathlist = os.path.split(fn)
    pathlist = pathlist[:-1] + ("projections",) + pathlist[-1:]
    return os.path.join(*pathlist)

def pngdir(fn):
    pathlist = os.path.split(fn)
    pathlist = pathlist[:-1] + ("projections","pngs") + pathlist[-1:]
    return os.path.join(*pathlist)

for mol in ("HCN","HNC","HCOp","HC3N","H2CS303-202","H2CS303202"):
    for vrange in ((30,51),(52,93),(94,134),(55,66)):
        for fn in (glob.glob("12m/*{0}*image.pbcor.fits".format(mol)) +
                   glob.glob("7m/*{0}*image.pbcor.fits".format(mol))
                   #glob.glob("feathered/Feathered*{0}*.fits".format(mol)) +
                   #glob.glob("feathered/Regridded*{0}*fits".format(mol))
                  ):

            if any(x in fn for x in ('mask','max','sum','mom','pvext')):
                continue

            rcube = (SpectralCube.read(fn)[:,32:-32,32:-32]
                     .with_spectral_unit(u.km/u.s, velocity_convention='radio')
                     .spectral_slab(vrange[0]*u.km/u.s, vrange[1]*u.km/u.s))
            mask = rcube > 0.5*rcube.unit
            cube = rcube.with_mask(mask)

            cubemax = rcube.max(axis=0)
            cubemax.quicklook()
            cubemax.FITSFigure.save(filename=pngdir(fn.replace(".fits",".v{0}to{1}.max.png".format(*vrange))))
            cubemax.hdu.writeto(pdir(fn.replace(".fits",".v{0}to{1}.max.fits".format(*vrange))),
                                clobber=True)
            cubemom0 = rcube.moment0(axis=0)
            cubemom0.quicklook()
            cubemom0.FITSFigure.save(filename=pngdir(fn.replace(".fits",".v{0}to{1}.mom0.png".format(*vrange))))
            cubemom0.hdu.writeto(pdir(fn.replace(".fits",".v{0}to{1}.mom0.fits".format(*vrange))),
                                 clobber=True)


            cubemom1 = cube.moment1(axis=0)
            cubemom1.quicklook()
            cubemom1.FITSFigure.show_colorscale(vmin=vrange[0], vmax=vrange[1])
            cubemom1.FITSFigure.save(pngdir(fn.replace(".fits",".v{0}to{1}.mom1.png".format(*vrange))))
            cubemom1.hdu.writeto(pdir(fn.replace(".fits",".v{0}to{1}.mom1.fits".format(*vrange))),
                                 clobber=True)

            cubemom2 = cube.moment2(axis=0)
            cubemom2.quicklook()
            cubemom2.FITSFigure.show_colorscale(vmin=00, vmax=50)
            cubemom2.FITSFigure.save(pngdir(fn.replace(".fits",".v{0}to{1}.mom2.png".format(*vrange))))
            cubemom2.hdu.writeto(pdir(fn.replace(".fits",".v{0}to{1}.mom2.fits".format(*vrange))),
                                 clobber=True)

            for ii in pl.get_fignums():
                pl.close(ii)

