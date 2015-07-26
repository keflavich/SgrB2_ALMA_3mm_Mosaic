import numpy as np
from astropy import units as u
import ds9
from spectral_cube import SpectralCube


if 'dd9' not in locals():
    dd9 = ds9.ds9()

niters = (1,5,100,1000,2000,2250,2500,2750,3000,3250, 3500, 5000,8000)
outcube = np.empty([216,216,len(niters)])

for ii,niter in enumerate(niters):
    fn = 'SgrB2_a_03_7M.HC3N.clean.niter{0}.image'.format(niter)
    cube = spectral_cube.SpectralCube.read(fn)[140:166]
    cube.to_ds9(dd9.id, newframe=True)
    dd9.set('wcs append', "OBJECT  = '{0}'".format(fn[:-6]))
    outcube[:,:,ii] = cube[10,:,:]

dd9.set('tile yes')
dd9.set('frame frameno 1')
dd9.set('frame delete')
dd9.set('scale limits -0.5 2')
dd9.set('wcs fk5')
dd9.set('lock frame wcs')
dd9.set('lock slice image')
dd9.set('lock scale yes')
dd9.set('lock color yes')
