from ch3cn_fits import SpectralCube, pyspeckit, fits, u

cubefn = '../FITS/merge/SgrB2_b3_7M_12M.CH3CN.image.pbcor_medsub.fits'
cube = SpectralCube.read(cubefn)
cubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                              cube.wcs.wcs.restfrq*u.Hz))
err = cubeK[:30].std(axis=0)
peak = fits.getdata('../FITS/merge/max/SgrB2_b3_7M_12M.CH3CN.image.pbcor_medsub_max.fits')
mask = peak > 0.5

pcube = pyspeckit.Cube(cube=cubeK)

pcube.fiteach(fittype='ch3cn', guesses=[0,5,150,1e15], integral=False,
              verbose_level=3, start_from_point=(470,429),
              use_neighbor_as_guess=True, position_order=1/peak,
              maskmap=mask,
              errmap=err.value, multicore=4)
