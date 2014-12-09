from astropy import units as u
import ds9
from spectral_cube import SpectralCube

files = [
    "Feathered_CH3CN.fits",
    "Feathered_H2CS303202.fits",
    "Feathered_H41a.fits",
    "Feathered_HC3N.fits",
    "Feathered_HCN.fits",
    "Feathered_HCOp.fits",
    "Feathered_HNC.fits",
]

dd = ds9.ds9()

vrange = [-95, 135]*u.km/u.s

for fn in files:
    species = 'feathered '+fn.split(".")[0].split("_")[1]
    if 'CH3CN' in species:
        rest_value = 91987.094*u.MHz
    else:
        rest_value = None
    cube = SpectralCube.read(fn).with_spectral_unit(u.km/u.s, rest_value=rest_value,
                                                    velocity_convention='radio').spectral_slab(vrange[0],
                                                                                               vrange[1])
    cube.to_ds9(dd.id, newframe=True)
    dd.set('wcs append', "OBJECT  = '{0}'".format(species))
    print cube

dd.set('tile yes')
dd.set('frame frameno 1')
dd.set('frame delete')
dd.set('scale limits -0.5 2')
dd.set('wcs fk5')
dd.set('lock frame wcs')
dd.set('lock slice image')
dd.set('lock scale yes')
dd.set('lock color yes')


#SgrB2_a_03_7M.lowres.spw0.fits
#SgrB2_a_03_7M.lowres.spw1.fits
#SgrB2_a_03_7M.lowres.spw2.fits
#SgrB2_a_03_7M.lowres.spw3.fits
