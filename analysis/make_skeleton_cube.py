import numpy as np
import paths
from astropy.io import fits
from astropy.utils.console import ProgressBar
import fil_finder
from spectral_cube import SpectralCube

cube = SpectralCube.read(paths.Fpath('12m/SgrB2_b3_12M.HC3N.image.pbcor.contsub.fits'))
subcube = cube[100:155].minimal_subcube()

newcube = np.empty_like(subcube, dtype='int16')

for ii,imslice in enumerate(ProgressBar(subcube)):
    img = imslice.value
    hdr = header = imslice.header
    try:
        ff = fil_finder.fil_finder_2D(img, hdr, beamwidth=3.5, distance=8500,
                                      skel_thresh=20, branch_thresh=10,
                                      glob_thresh=65, adapt_thresh=10,
                                      pad_size=0,
                                     )
    except ValueError:
        print("FilFinder creation failed on iter {0}".format(ii))
        continue
    try:
        mask = ff.create_mask(verbose=True)
    except ValueError:
        print("Mask creation failed on iter {0}".format(ii))
        continue
    try:
        skels = ff.medskel(verbose=True)
    except ValueError:
        print("Skeleton creation failed on iter {0}".format(ii))
        continue

    try:
        anaskel = ff.analyze_skeletons(verbose=True, skel_thresh=40, branch_thresh=40)
    except ValueError:
        print("Skeleton analysis failed on iter {0}".format(ii))
        continue

    newcube[ii,:,:] = skels.skeleton

newhdu = fits.PrimaryHDU(data=newcube.astype('int16'), header=subcube.header)
newhdu.writeto(paths.Fpath('12m/SgrB2_b3_12M.HC3N.image.pbcor.skeletons.fits'),
               clobber=True)
