from spectral_cube import SpectralCube

for fn in glob.glob("*image.pbcor.fits"):
    rcube=SpectralCube.read(fn)
    mask = rcube > 0.5
    cube = rcube.with_mask(mask)
    
    cubesum = cube.sum(axis=0)
    cubesum.hdu.writeto(fn.replace(".fits",".mask.sum.fits"),clobber=True)
    cubemax = cube.max(axis=0)
    cubemax.hdu.writeto(fn.replace(".fits",".mask.max.fits"),clobber=True)
    cubemom0 = cube.moment0(axis=0)
    cubemom0.hdu.writeto(fn.replace(".fits",".mask.mom0.fits"),clobber=True)
    cubemom1 = cube.moment1(axis=0)
    cubemom1.hdu.writeto(fn.replace(".fits",".mask.mom1.fits"),clobber=True)
    print('masked',fn)

for fn in glob.glob("*image.pbcor.fits"):
    cube=SpectralCube.read(fn)
    cubesum = cube.sum(axis=0)
    cubesum.hdu.writeto(fn.replace(".fits",".sum.fits"),clobber=True)
    cubemax = cube.max(axis=0)
    cubemax.hdu.writeto(fn.replace(".fits",".max.fits"),clobber=True)
    cubemom0 = cube.moment0(axis=0)
    cubemom0.hdu.writeto(fn.replace(".fits",".mom0.fits"),clobber=True)
    cubemom1 = cube.moment1(axis=0)
    cubemom1.hdu.writeto(fn.replace(".fits",".mom1.fits"),clobber=True)
    print(fn)
