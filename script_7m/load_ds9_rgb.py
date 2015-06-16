import ds9
from spectral_cube import SpectralCube

dd = ds9.ds9()
dd.set('rgb')

cubehcn = SpectralCube.read('Feathered_HCN.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=88.6316024*u.GHz).spectral_slab(-95*u.km/u.s, 140*u.km/u.s)
cubehnc = SpectralCube.read('Feathered_HNC.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_slab(-95*u.km/u.s, 140*u.km/u.s)
cubehcop = SpectralCube.read('Feathered_HCOp.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_slab(-95*u.km/u.s, 140*u.km/u.s)

dd.set('rgb channel red')
cubehcn.to_ds9(dd.id)
dd.set('scale limits -0.1 3')
dd.set('rgb channel green')
cubehnc.to_ds9(dd.id)
dd.set('scale limits -0.1 3')
dd.set('rgb channel blue')
cubehcop.to_ds9(dd.id)
dd.set('scale limits -0.1 3')

dd.set('rgb lock slice yes')
