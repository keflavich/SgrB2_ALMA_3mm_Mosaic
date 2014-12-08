import pylab as pl
from astropy import units as u
from spectral_cube import SpectralCube
import glob

for fn in glob.glob("*image.pbcor.fits"):
    rcube = (SpectralCube.read(fn).with_spectral_unit(u.km/u.s,
                                                     velocity_convention='radio')
             .spectral_slab(13*u.km/u.s, 136*u.km/u.s))
    mask = rcube > 0.5
    cube = rcube.with_mask(mask)

    cubesum = cube.sum(axis=0)
    cubesum.quicklook(filename=fn.replace(".fits",".mask.sum.png"))
    cubesum.hdu.writeto(fn.replace(".fits",".mask.sum.fits"),clobber=True)
    cubemax = cube.max(axis=0)
    cubemax.quicklook(filename=fn.replace(".fits",".mask.max.png"))
    cubemax.hdu.writeto(fn.replace(".fits",".mask.max.fits"),clobber=True)
    cubemom0 = cube.moment0(axis=0)
    cubemom0.quicklook(filename=fn.replace(".fits",".mask.mom0.png"))
    cubemom0.hdu.writeto(fn.replace(".fits",".mask.mom0.fits"),clobber=True)
    cubemom1 = cube.moment1(axis=0)
    cubemom1.quicklook(filename=fn.replace(".fits",".mask.mom1.png"))
    cubemom1.hdu.writeto(fn.replace(".fits",".mask.mom1.fits"),clobber=True)
    cubemom2 = cube.moment2(axis=0)
    cubemom2.quicklook(filename=fn.replace(".fits",".mask.mom2.png"))
    cubemom2.hdu.writeto(fn.replace(".fits",".mask.mom2.fits"),clobber=True)
    print('masked',fn)

    #ytc = cube.to_yt()
    #name = fn.replace(".image.pbcor.fits","")
    #ytc.quick_render_movie('movies/'+name, size=512, nframes=60)
    #ytc.quick_isocontour(title=name,
    #                     description='Sgr B2 7m ACA/ALMA data at 3mm: Isocontours of {0}'.format(name))
    

for fn in glob.glob("*image.pbcor.fits"):
    cube = (SpectralCube.read(fn).with_spectral_unit(u.km/u.s,
                                                     velocity_convention='radio')
            .spectral_slab(13*u.km/u.s, 136*u.km/u.s))
    cubesum = cube.sum(axis=0)
    cubesum.quicklook(filename=fn.replace(".fits",".sum.png"))
    cubesum.hdu.writeto(fn.replace(".fits",".sum.fits"),clobber=True)
    cubemax = cube.max(axis=0)
    cubemax.quicklook(filename=fn.replace(".fits",".max.png"))
    cubemax.hdu.writeto(fn.replace(".fits",".max.fits"),clobber=True)
    cubemom0 = cube.moment0(axis=0)
    cubemom0.quicklook(filename=fn.replace(".fits",".mom0.png"))
    cubemom0.hdu.writeto(fn.replace(".fits",".mom0.fits"),clobber=True)
    cubemom1 = cube.moment1(axis=0)
    cubemom1.quicklook(filename=fn.replace(".fits",".mom1.png"))
    cubemom1.hdu.writeto(fn.replace(".fits",".mom1.fits"),clobber=True)
    cubemom2 = cube.moment2(axis=0)
    cubemom2.quicklook(filename=fn.replace(".fits",".mom2.png"))
    cubemom2.hdu.writeto(fn.replace(".fits",".mom2.fits"),clobber=True)
    print(fn)

for fn in glob.glob("feathered/Feathered*.fits")+glob.glob("feathered/Regridded*fits"):
    if any(x in fn for x in ('max','mom','sum')):
       continue
    rcube = (SpectralCube.read(fn)[:,32:-32,32:-32]
             .with_spectral_unit(u.km/u.s, velocity_convention='radio')
             .spectral_slab(13*u.km/u.s, 136*u.km/u.s))
    mask = rcube > 0.5
    cube = rcube.with_mask(mask)

    cubesum = cube.sum(axis=0)
    cubesum.quicklook(filename=fn.replace(".fits",".mask.sum.png"))
    cubesum.hdu.writeto(fn.replace(".fits",".mask.sum.fits"),clobber=True)
    cubemax = cube.max(axis=0)
    cubemax.quicklook(filename=fn.replace(".fits",".mask.max.png"))
    cubemax.hdu.writeto(fn.replace(".fits",".mask.max.fits"),clobber=True)
    cubemom0 = cube.moment0(axis=0)
    cubemom0.quicklook(filename=fn.replace(".fits",".mask.mom0.png"))
    cubemom0.hdu.writeto(fn.replace(".fits",".mask.mom0.fits"),clobber=True)
    cubemom1 = cube.moment1(axis=0)
    cubemom1.quicklook(filename=fn.replace(".fits",".mask.mom1.png"))
    cubemom1.hdu.writeto(fn.replace(".fits",".mask.mom1.fits"),clobber=True)
    cubemom2 = cube.moment2(axis=0)
    cubemom2.quicklook(filename=fn.replace(".fits",".mask.mom2.png"))
    cubemom2.hdu.writeto(fn.replace(".fits",".mask.mom2.fits"),clobber=True)
    print('masked',fn)

    #ytc = cube.to_yt()
    #name = fn.replace("Feathered_","")
    #ytc.quick_render_movie('movies_feathered/'+name, size=512, nframes=60)
    #ytc.quick_isocontour(title=name,
    #                     description='Sgr B2 7m ACA/ALMA data at 3mm combined with MOPRA data: Isocontours of {0}'.format(name))
