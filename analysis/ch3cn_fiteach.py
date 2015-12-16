from ch3cn_fits import SpectralCube, pyspeckit, fits, u, np
import os
T=True
F=False

cubefn = '../FITS/merge/SgrB2_b3_7M_12M.CH3CN.image.pbcor_medsub.fits'
if not os.path.exists(cubefn):
    cube = SpectralCube.read('../FITS/merge/SgrB2_b3_7M_12M.CH3CN.image.pbcor.fits').minimal_subcube()
    med = cube.percentile(25, axis=0)
    cube.allow_huge_operations=True
    cube = cube - med
    cube.write(cubefn)
else:
    cube = SpectralCube.read(cubefn).minimal_subcube()
cubeK = cube.to(u.K, u.brightness_temperature(cube.beam,
                                              cube.wcs.wcs.restfrq*u.Hz))
err = cubeK[:30].std(axis=0)
peak = cubeK.max(axis=0)
mask = (peak > 1*u.K) & (peak > 6*err)

subcube = cubeK.spectral_slab(-50*u.km/u.s, 100*u.km/u.s)

pcube = pyspeckit.Cube(cube=subcube)

vguesses = subcube.spectral_axis[subcube.argmax(axis=0)]
colguesses = np.ones_like(mask)*1e15
temguesses = np.ones_like(mask)*150.
widths = np.ones_like(mask)*5.0
guesses = np.array([vguesses.value, widths, temguesses, colguesses])

# For laptop
#mask &= (peak>10*u.K)

start_point = np.unravel_index(np.nanargmax(peak*mask), peak.shape)

position_order = 1./peak.value
position_order[np.isnan(peak)] = np.inf

pcube.fiteach(fittype='ch3cn', guesses=guesses, integral=False,
              verbose_level=3, start_from_point=start_point,
              use_neighbor_as_guess=True, position_order=position_order,
              limitedmax=[T,T,T,T],
              maxpars=[100,15,500,1e17],
              maskmap=mask,
              errmap=err.value, multicore=4)
