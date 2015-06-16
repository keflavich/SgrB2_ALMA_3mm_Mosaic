from __future__ import print_function
import os
# TODO: put this in a package
execfile("/Users/adam/repos/w51evlareductionscripts/singledish_combine/singledish_combine.py")
from spectral_cube import SpectralCube

mopradir = '/Volumes/passport/alma/sgrb2_b3/mopra'
acadir = '/Volumes/passport/alma/sgrb2_b3/FITS'
outdir = '/Volumes/passport/alma/sgrb2_b3/FITS/feathered'

pairs = [('HCOp',  'CMZ_3mm_HCO+.fits', 'SgrB2_a_03_7M.HCOp.image.pbcor.fits', 89188526000.0),
         ('HCN',   'CMZ_3mm_HCN.fits', 'SgrB2_a_03_7M.HCN.image.pbcor.fits',  88.6316024e9),
         ('HNC',   'CMZ_3mm_HNC.fits', 'SgrB2_a_03_7M.HNC.image.pbcor.fits', 90663564000.0),
         ('CH3CN', 'CMZ_3mm_CH3CN.fits', 'SgrB2_a_03_7M.CH3CN_5-4_3.image.pbcor.fits', 91.98705*1e9), #91971465000.0),
         ('HC3N',  'CMZ_3mm_HC3N.fits', 'SgrB2_a_03_7M.HC3N.dirty.image.fits', 92034434000.0),
         # Wrong name in MOPRA files
         ('H41a',  'CMZ_3mm_RRL_H42a.fits', 'SgrB2_a_03_7M.H41a.image.pbcor.fits', 92034434000.0),
         ('H2CS303202',  'SgrB2_103_H2CS_nb_cube.fits', 'SgrB2_a_03_7M.H2CS303-202.image.pbcor.fits', 104.61711e9),
         #('13CS',  'CMZ_3mm_13CS.fits', 'SgrB2_a_03_7M.13CS.image.pbcor.fits', 92.49431e9), # BEYOND OUR RANGE =(
         ]
pairs = [ ('HCN',   'CMZ_3mm_HCN.fits', 'SgrB2_a_03_7M.HCN.image.pbcor.fits',
           88.6316024e9),]
 
for species,mopra,aca,restfrq in pairs:

    print(species)
    mopracube = SpectralCube.read(os.path.join(mopradir, mopra))
    acacube = SpectralCube.read(os.path.join(acadir, aca))
    mc = mopracube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.Hz,
                                      velocity_convention='radio')
    acac = acacube.with_spectral_unit(u.km/u.s, rest_value=restfrq*u.Hz,
                                      velocity_convention='radio')

    h1 = acacube.header
    jyk = (1*u.Jy).to(u.K, u.brightness_temperature(h1['BMAJ']*u.deg *
                                                    h1['BMIN']*u.deg * 2 *
                                                    np.pi, h1['RESTFRQ']*u.Hz))

    combcube,f2 = fourier_combine_cubes(acac, mc, highresscalefactor=jyk.value,
                                        lowresfwhm=0.65*u.arcmin,
                                        return_regridded_cube2=True)
    f2.writeto(os.path.join(outdir, "Regridded_"+mopra), clobber=True)
    f1 = fits.open(os.path.join(acadir, aca))
    f1[0].data = combcube
    f1[0].header = acacube.header
    f1.writeto(os.path.join(outdir, 'Feathered_{0}.fits'.format(species)),
               clobber=True)
    print()
