from __future__ import print_function
# TODO: put this in a package
execfile("/Users/adam/repos/w51evlareductionscripts/singledish_combine/singledish_combine.py")
from spectral_cube import SpectralCube


pairs = [('HCOp',  'CMZ_3mm_HCO+.fits', 'SgrB2_a_03_7M.HCOp.image.pbcor.fits', 89188526000.0),
         ('HCN',   'CMZ_3mm_HCN.fits', 'SgrB2_a_03_7M.HCN.image.pbcor.fits', 88633936000.0),
         ('HNC',   'CMZ_3mm_HNC.fits', 'SgrB2_a_03_7M.HNC.image.pbcor.fits', 90663564000.0),
         ('CH3CN', 'CMZ_3mm_CH3CN.fits', 'SgrB2_a_03_7M.CH3CN_5-4_3.image.pbcor.fits', 91.98705*1e9), #91971465000.0),
         ('HC3N',  'CMZ_3mm_HC3N.fits', 'SgrB2_a_03_7M.HC3N.image.pbcor.fits', 92034434000.0),
         # Wrong name in MOPRA files
         ('H41a',  'CMZ_3mm_RRL_H42a.fits', 'SgrB2_a_03_7M.H41a.image.pbcor.fits', 92034434000.0),
         ]
 
for species,mopra,aca,restfrq in pairs:

    print(species)
    mopracube = SpectralCube.read(mopra)
    acacube = SpectralCube.read(aca)
    mc = mopracube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.Hz,
                                      velocity_convention='radio')
    acac = acacube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.Hz,
                                      velocity_convention='radio')

    h1 = acacube.header
    jyk = (1*u.Jy).to(u.K, u.brightness_temperature(h1['BMAJ']*u.deg *
                                                    h1['BMIN']*u.deg * 2 *
                                                    np.pi, h1['RESTFRQ']*u.Hz))

    combcube = fourier_combine_cubes(acac, mc, highresscalefactor=jyk.value,
                                     lowresfwhm=0.65*u.arcmin)
    f1 = fits.open(aca)
    f1[0].data = combcube
    f1.writeto('Feathered_{0}.fits'.format(species), clobber=True)
    print()
