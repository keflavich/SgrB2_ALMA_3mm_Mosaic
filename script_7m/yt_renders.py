from spectral_cube import SpectralCube, BooleanArrayMask
from astropy import units as u
import numpy as np
import matplotlib as mpl
import yt

red = mpl.colors.LinearSegmentedColormap('red',
                                         {'red':[(0,1,1),(1,1,1)],
                                          'green':[(0,0,0),(1,1,1)],
                                          'blue':[(0,0,0),(1,1,1)]})
green = mpl.colors.LinearSegmentedColormap('green',
                                           {'red':[(0,0,0),(1,1,1)],
                                            'green':[(0,1,1),(1,1,1)],
                                            'blue':[(0,0,0),(1,1,1)]})
blue = mpl.colors.LinearSegmentedColormap('blue', {'red':[(0,0,0),(1,1,1)],
                                                   'green':[(0,0,0),(1,1,1)],
                                                   'blue':[(0,1,1),(1,1,1)]})

brown = mpl.colors.LinearSegmentedColormap('brown', {'red':[(0,0.5,0.5),(1,1,1)],
                                                     'green':[(0,0.25,0.25),(1,1,1)],
                                                     'blue':[(0,0,0),(1,1,1)]})

molecules = ('HNC','HC3N','HCOp','HCN',)
molecules = ('HC3N',)

colors = {'HNC': red,
          'HCN': green,
          'HC3N': brown,
          'HCOp': blue,
         }

cubes = {mol: SpectralCube.read('Feathered_{0}.fits'.format(mol))
                          .with_spectral_unit(u.km/u.s,
                                              velocity_convention='radio')
         for mol in molecules}
hc3nmean = cubes['HC3N'].mean(axis=(1,2))
bad = np.abs(hc3nmean)>1
mh3cn = cubes['HC3N'].with_mask(BooleanArrayMask(bad[:,None,None],
                                                 cubes['HC3N'].wcs))
cubes['HC3N'] = mh3cn

ytcubes = {mol: cubes[mol].to_yt() for mol in molecules}

try:
    os.mkdir("ChemMovie")
except:
    pass

def ramp(vals, minval, maxval):
    v = (vals)**0.25
    return (v - v.min())/(v.max() - v.min())

nframes = 60
size = 512

all_images = {}
for mol in molecules:
    print "MAKING MOLECULE ",mol
    transfer_function = ytcubes[mol].auto_transfer_function([1.0, 4.0],
                                                            colormap=colors[mol])
    #transfer_function.tf.add_layers(10, w=0.01, colormap=colors[mol],
    #                                alpha=np.linspace(0.3,1,10))
    transfer_function.tf.map_to_colormap(1, 4, colormap=colors[mol],
                                         scale_func=ramp)
    transfer_function.tf.show()
    print "ROTATION SET 1",mol
    images1 = ytcubes[mol].quick_render_movie('ChemMovie/',
                                    size=size,
                                    nframes=nframes,
                                    start_index=0,
                                    image_prefix=mol,
                                    transfer_function=transfer_function.tf,
                                   )
    print "ROTATION SET 2",mol
    images2 = ytcubes[mol].quick_render_movie('ChemMovie/',
                                    size=size,
                                    nframes=nframes,
                                    start_index=nframes,
                                    camera_angle=(1,0,0),
                                    rot_vector=(0,1,0),
                                    image_prefix=mol,
                                    transfer_function=transfer_function.tf,
                                   )
    print "ROTATION SET 3",mol
    images3 = ytcubes[mol].quick_render_movie('ChemMovie/',
                                    size=size,
                                    nframes=nframes,
                                    start_index=nframes*2,
                                    camera_angle=(0,0,1),
                                    rot_vector=(0.5,0.5,0),
                                    image_prefix=mol,
                                    transfer_function=transfer_function.tf,
                                   )
    all_images[mol] = images1+images2+images3

molecules = ('HNC','HCOp','HCN','HC3N')
from matplotlib import image
from yt.data_objects.image_array import ImageArray
combined = [np.sum([image.imread("ChemMovie/{mol}{ii:04d}.png".format(mol=mol,
                                                                      ii=ii))[:,:,:3]
                    for mol in molecules], axis=0)
            for ii in range(nframes*3)]
cmax = max(np.percentile(img[:, :, :3].sum(axis=2), 99.5) for img in combined)
for ii, img in enumerate(combined):
    img = ImageArray(img)
    img = img.rescale(cmax=cmax).swapaxes(0,1)
    img.write_png("ChemMovie/combined_%04i.png" % (ii), rescale=False)

import spectral_cube.ytcube
spectral_cube.ytcube._make_movie('ChemMovie', prefix='combined_',
                                 filename='combined_%i.mp4' % size)
