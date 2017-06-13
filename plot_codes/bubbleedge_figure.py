import paths
import pyregion
import aplpy
import numpy as np
import pylab as pl

regions = pyregion.open(paths.rpath('bubble_edges.reg'))

fig1 = pl.figure(1)
fig1.clf()
F = aplpy.FITSFigure(paths.tmpath('max/SgrB2_b3_12M.HC3N.image.pbcor.contsub_max.fits'),
                     figure=fig1)

F.show_grayscale()

cm = pl.cm.rainbow
norm = pl.Normalize(vmin=18, vmax=56)

for reg in regions:
    seg = np.array(list(zip(reg.coord_list[::2], reg.coord_list[1::2]))).T
    vel = float(reg.attr[1]['text'])

    color = cm(norm(vel))

    F.show_lines([seg], color=color)

F.recenter(266.83188, -28.391511, radius=0.0585)
pl.draw()
pl.show()
