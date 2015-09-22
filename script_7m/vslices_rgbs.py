import aplpy
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

for globstr,zranges in ((("12m/*{0}*image.pbcor.fits"),     ((-3,35),(-5,60),(-5,40),(-5,20))),
                        (("7m/*{0}*image.pbcor.fits"),      ((-3,35),(-5,60),(-5,40),(-5,20))),
                        (("feathered/Feathered*{0}*.fits"), ((-3,25),(-5,60),(-5,40),(-5,20))),
                        (("feathered/Regridded*{0}*.fits"), ((-3,35), (-5,60), (-5,40), (-5,20))),
                       ):
    for (vmin,vmax),(zmin,zmax) in zip(((30,51),(52,93),(94,134),(55,66)),
                                       zranges):
        rgbfn = [glob.glob(pdir(globstr.format(mol)
                                .replace(".fits",".v{0}to{1}.mom0.fits".format(vmin,vmax))))[0]
                 for mol in ("HCN","HNC","HCO*")
                 ]
        if 'Regridded' in globstr:
            txt = 'single'
        else:
            txt = ''
        out = pngdir(os.path.join(os.path.split(globstr)[0],
                                  "rgb_hcn_hnc_hcop.{2}v{0}to{1}.mom0.png".format(vmin,vmax, txt)))
        print rgbfn,out

        aplpy.rgb.make_rgb_image(rgbfn, out, vmin_r=zmin, vmin_g=zmin,
                                 vmin_b=zmin, vmax_r=zmax, vmax_g=zmax,
                                 vmax_b=zmax)
