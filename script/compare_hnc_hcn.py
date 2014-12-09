from spectral_cube import SpectralCube
import copy

#cubehcn = SpectralCube.read('Feathered_HCN.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=88.6316024*u.GHz).spectral_slab(-95*u.km/u.s, 140*u.km/u.s)
#cubehnc = SpectralCube.read('Feathered_HNC.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_slab(-95*u.km/u.s, 140*u.km/u.s)
cubehcn = SpectralCube.read('SgrB2_a_03_7M.HCN.image.pbcor.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=88.6316024*u.GHz).spectral_slab(-90*u.km/u.s, 130*u.km/u.s)
cubehnc = SpectralCube.read('SgrB2_a_03_7M.HNC.image.pbcor.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_slab(-90*u.km/u.s, 130*u.km/u.s)

cubehnc._wcs = cubehcn.wcs
cubehnc.mask._wcs = cubehcn.wcs

mask  = (cubehcn > 0.5) & (cubehnc.mask)
mask2 = (cubehcn > 0.5) & (cubehnc.mask)
mask2._wcs = cubehnc.wcs
hcn_flat = cubehcn.with_mask(mask).flattened()
hnc_flat = cubehnc.with_mask(mask2).flattened()

mask3 = cubehnc > 0.5
mask4 = cubehnc > 0.5
mask4._wcs = cubehcn.wcs
hnc_flat2 = cubehnc.with_mask(mask3).flattened()
hcn_flat2 = cubehcn.with_mask(mask4).flattened()

import pylab as pl
pl.clf()
pl.plot(hnc_flat, hcn_flat, '.', alpha=0.05)
pl.plot(hnc_flat2, hcn_flat2, '.', alpha=0.05)
pl.plot([0,5],[0,5],'k--')
pl.axis([0,5,0,5])
