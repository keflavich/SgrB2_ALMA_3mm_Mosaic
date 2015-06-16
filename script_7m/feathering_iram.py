from __future__ import print_function
# TODO: put this in a package
execfile("/Users/adam/repos/w51evlareductionscripts/singledish_combine/singledish_combine.py")
from spectral_cube import SpectralCube


pairs = [ 
         ('CH3CN', 'CH3CN54_hires_lmv.fits', 'SgrB2_a_03_7M.CH3CN_5-4_3.image.pbcor.fits', 91.98705*1e9), #91971465000.0),
         ('H41a',  'H41a_lmv.fits', 'SgrB2_a_03_7M.H41a.image.pbcor.fits', 92034434000.0),
         ]
 
for species,iram,aca,restfrq in pairs:

    print(species)
    iramcube = SpectralCube.read(iram)
    acacube = SpectralCube.read(aca)
    mc = iramcube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.Hz,
                                      velocity_convention='radio')
    acac = acacube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.Hz,
                                      velocity_convention='radio')

    h1 = acacube.header
    jyk = (1*u.Jy).to(u.K, u.brightness_temperature(h1['BMAJ']*u.deg *
                                                    h1['BMIN']*u.deg * 2 *
                                                    np.pi, h1['RESTFRQ']*u.Hz))

    combcube,f2 = fourier_combine_cubes(acac, mc, highresscalefactor=jyk.value,
                                        lowresfwhm=0.30*u.arcmin,
                                        return_regridded_cube2=True)
    f2.writeto("Regridded_"+iram, clobber=True)
    f1 = fits.open(aca)
    f1[0].data = combcube
    f1[0].header = acacube.header
    f1.writeto('Feathered_{0}.fits'.format(species), clobber=True)
    print()
