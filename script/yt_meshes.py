from spectral_cube import SpectralCube, BooleanArrayMask
from astropy import units as u
import numpy as np
import matplotlib as mpl
import yt
import os
from restfreqs import restfreqs

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
from yt.visualization._colormap_data import color_map_luts
color_map_luts['red'] = red(np.linspace(0,1,256)).T
color_map_luts['green'] = green(np.linspace(0,1,256)).T
color_map_luts['blue'] = blue(np.linspace(0,1,256)).T
color_map_luts['brown'] = brown(np.linspace(0,1,256)).T

molecules = ('HNC','HC3N','HCN',)

colors = {'HNC': red,
          'HCN': green,
          'HC3N': blue,
          #'HCOp': blue,
         }
colornames = {'HNC': 'red',
              'HCN': 'green',
              'HC3N': 'blue',
              #'HCOp': 'blue',
             }
levels = {'HNC':(0.5,1.0,1.5),
          'HCN':(1.2,1.5,2.2),
          'HC3N':(0.4,1,1.2),
         }

#cubes = {mol: SpectralCube.read('SgrB2_a_03_7M.{0}.image.pbcor.fits'.format(mol))
#                          .with_spectral_unit(u.km/u.s,
#                                              velocity_convention='radio',
#                                              rest_value=restfreqs[mol])
#         for mol in molecules}
cubes = {mol: SpectralCube.read('feathered/Feathered_{0}.fits'.format(mol))
                          .with_spectral_unit(u.km/u.s,
                                              velocity_convention='radio',
                                              rest_value=restfreqs[mol])
                          .spectral_slab(-20*u.km/u.s, 137*u.km/u.s)
         for mol in molecules}
hc3nmean = cubes['HC3N'].mean(axis=(1,2))
bad = np.abs(hc3nmean)>1
mhc3n = cubes['HC3N'].with_mask(BooleanArrayMask(~bad[:,None,None],
                                                 cubes['HC3N'].wcs))
cubes['HC3N'] = mhc3n
for cube in cubes.values():
    cube._data -= cube.mean(axis=0).value

ytcubes = {mol: cubes[mol][::2,::2,::2].to_yt() for mol in molecules}

for ytc in ytcubes.values():
    ytc.dataset.periodicity = (True,)*3

try:
    os.mkdir("IsoSurfs")
except:
    pass


#for ii,(level,transparency) in enumerate(zip((2, 3, 4), (0.1,0.3,0.6))):
for ii,(level,transparency) in enumerate(zip((0,1,2), (0.1, 0.2, 0.6))):
    surfaces = {mol:
                ytcubes[mol].dataset.surface(ytcubes[mol].dataset.all_data(),
                                             ytcubes[mol].dataset.field_list[0],
                                             levels[mol][level])
                for mol in molecules}

    for jj,(mol,surf) in enumerate(surfaces.iteritems()):
        #filename = os.path.join('IsoSurfs',
        #                        '{mol}_surf'.format(mol=mol, level=level,
        #                                            transparency=transparency))
        #surf.export_obj(filename, transparency=transparency,
        #                color_field=ytcubes[mol].dataset.field_list[0],
        #                #color_map=colors[mol],
        #                plot_index=ii)
        filename = os.path.join('IsoSurfs',
                                'all_surfs')
        surf.export_obj(filename, transparency=transparency,
                        color_field='ones', #ytcubes[mol].dataset.field_list[0],
                        color_map=colornames[mol],
                        plot_index=jj+ii*len(surfaces))

import zipfile
zfn = 'IsoSurfs/all_surfs12.zip'
zf = zipfile.ZipFile(zfn, mode='w')
zf.write('IsoSurfs/all_surfs.obj')
zf.write('IsoSurfs/all_surfs.mtl')
zf.close()

from yt.config import ytcfg
api_key = ytcfg.get("yt","sketchfab_api_key")
import requests
import os
data = {'title': 'Sgr B2 meshes colored (try 5)',
        'token': api_key,
        'fileModel': zfn,
        'filenameModel': os.path.basename(zfn)}
response = requests.post('https://api.sketchfab.com/v1/models', files=data)
response.raise_for_status()
