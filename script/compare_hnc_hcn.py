from spectral_cube import SpectralCube, LazyMask, BooleanArrayMask
import copy
from scipy import stats

#cubehcn = SpectralCube.read('Feathered_HCN.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=88.6316024*u.GHz).spectral_slab(-95*u.km/u.s, 140*u.km/u.s)
#cubehnc = SpectralCube.read('Feathered_HNC.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_slab(-95*u.km/u.s, 140*u.km/u.s)
cubehcn = SpectralCube.read('SgrB2_a_03_7M.HCN.image.pbcor.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=88.6316024*u.GHz).spectral_slab(-90*u.km/u.s, 130*u.km/u.s)
cubehcn = SpectralCube.read('SgrB2_a_03_7M.HCN.image.pbcor.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=88.63042*u.GHz).spectral_slab(-90*u.km/u.s, 130*u.km/u.s)
cubehnc = SpectralCube.read('SgrB2_a_03_7M.HNC.image.pbcor.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_slab(-90*u.km/u.s, 130*u.km/u.s)
#hcn_freqs = [88.63042, 88.6316, 88.63185, 88.63394, 88.6316024]*u.GHz
#for f in hcn_freqs:
#    vrange=[-16,120]*u.km/u.s
#    cubehnc = SpectralCube.read('SgrB2_a_03_7M.HNC.image.pbcor.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio').spectral_slab(vrange[0],vrange[1])
#    cubehcn = SpectralCube.read('SgrB2_a_03_7M.HCN.image.pbcor.fits').with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=f).spectral_slab(vrange[0],vrange[1])
#    cubehcn=cubehcn.with_mask(BooleanArrayMask(cubehnc.mask.include(), cubehcn.wcs))
#    print f, stats.pearsonr(cubehnc.flattened(), cubehcn.flattened())
#    88.63042 GHz (0.52182746, 0.0)
#    88.6316 GHz (0.49761587, 0.0)
#    88.63185 GHz (0.48257208, 0.0)
#    88.63394 GHz (0.3509686, 0.0)
#    88.6316024 GHz (0.49761587, 0.0)

from FITS_tools import regrid_cube
hnc_regrid = regrid_cube(cubehnc.filled_data[:], cubehnc.header, cubehcn.header)
cubehnc = SpectralCube(data=hnc_regrid, wcs=cubehcn.wcs)

#cubehnc._wcs = cubehcn.wcs
#cubehnc.mask._wcs = cubehcn.wcs

mask  = (cubehcn > 0.5)
mask2 = (cubehcn > 0.5)
#mask2._wcs = cubehnc.wcs
hcn_flat = cubehcn.with_mask(mask).flattened()
hnc_flat = cubehnc.with_mask(mask2).flattened()

mask3 = cubehnc > 0.5
mask4 = cubehnc > 0.5
#mask4._wcs = cubehcn.wcs
hnc_flat2 = cubehnc.with_mask(mask3).flattened()
hcn_flat2 = cubehcn.with_mask(mask4).flattened()

import pylab as pl
pl.clf()
pl.plot(hnc_flat, hcn_flat, ',', alpha=0.5)
pl.plot(hnc_flat2, hcn_flat2, ',', alpha=0.5)
pl.plot([0,5],[0,5],'k--')
pl.axis([0,5,0,5])

dd = cubehcn.to_ds9()
cubehnc.to_ds9(dd.id, newframe=True)
dd.set('tile yes')
dd.set('scale limits -0.5 2')
dd.set('wcs fk5')
dd.set('lock frame wcs')
dd.set('lock slice image')
dd.set('lock scale yes')
dd.set('lock color yes')

