import numpy as np
import paths
from astropy.table import Table
from astropy import units as u
from astropy import coordinates
import pylab as pl
from astropy.io import fits
from astropy import wcs
import astropy.visualization
#from wcsaxes import WCS as WCSaxes
from astropy.convolution import convolve, Gaussian2DKernel
from mpl_plot_templates import asinh_norm
import matplotlib
from collections import defaultdict
from files import contfilename
import warnings
from visualization import make_scalebar
from constants import distance
from overlay_common import core_phot_tbl, plotcores
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)


pl.rcParams['image.interpolation'] = 'nearest'
pl.rcParams['image.origin'] = 'lower'
if int(matplotlib.__version__[0]) >= 2:
    pl.rcParams['figure.dpi'] = 75.
    pl.rcParams['savefig.dpi'] = 300.
    pl.rcParams['axes.labelsize'] = 8
    pl.rcParams['axes.titlesize'] = 8
    pl.rcParams['xtick.labelsize'] = 6
    pl.rcParams['ytick.labelsize'] = 6
    pl.rcParams['axes.linewidth'] = 0.15
    tick_fontsize = 6
    markersize = 3
else:
    pl.rcParams['savefig.dpi'] = 300.
    tick_fontsize = 10
    markersize = 3
pl.rcParams['xtick.direction'] = 'in'
pl.rcParams['ytick.direction'] = 'in'

zoomregions = {'SouthOfSouth':
               {'bottomleft': coordinates.SkyCoord("17:47:21.352",
                                                   "-28:24:06.985",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:19.587",
                                                 "-28:23:53.629",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'DeepSouth',
                'bbox':[-0.075,1],
                'loc': 2,
                'l1':3,
                'l2':1,
                'min': -0.2,
                'max': 2.5,
                'zoom': 3,
               },
               'MidDS':
               {'bottomleft': coordinates.SkyCoord("17:47:22.272",
                                                   "-28:24:44.405",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:20.835",
                                                 "-28:24:16.756",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'DeepSouth',
                'bbox':[-0.070,0.7],
                'loc': 2,
                'l1':1,
                'l2':4,
                'min': -0.2,
                'max': 2.5,
                'zoom': 3,
               },
               'LowerDS':
               {'bottomleft': coordinates.SkyCoord("17:47:23.673",
                                                   "-28:25:54.266",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:20.271",
                                                 "-28:25:40.432",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'DeepSouth',
                'bbox':[-0.075,0.05],
                'loc': 2,
                'l1':2,
                'l2':4,
                'min': -0.2,
                'max': 2.5,
                'zoom': 3,
               },
               'M':
               {'bottomleft': coordinates.SkyCoord("17:47:20.929",
                                                   "-28:23:12.813",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:18.469",
                                                 "-28:22:49.771",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'MandN',
                'bbox': [-0.25,0.45],
                #'bbox':[-0.024,0.425],
                'loc': 2,
                'l1':1,
                'l2':4,
                'min': -1,
                'max': 50,
                'zoom': 3,
               },
               'Mfull':
               {'bottomleft': coordinates.SkyCoord("17:47:20.929",
                                                   "-28:23:12.9",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:18.469",
                                                 "-28:22:49.5",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'full',
                'bbox': [0.0,0.6],
                #'bbox':[-0.024,0.425],
                'loc': 2,
                'l1':1,
                'l2':4,
                'min': -1,
                'max': 50,
                'zoom': 5,
               },
               'N':
               {'bottomleft': coordinates.SkyCoord("17:47:21.006",
                                                   "-28:22:24.555",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:18.206",
                                                 "-28:22:02.383",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'MandN',
                'bbox':[-0.4,0.95],
                'loc': 2,
                'l1':1,
                'l2':3,
                'min': -1,
                'max': 50,
                'zoom': 3,
               },
               'M_inner':
               {'bottomleft': coordinates.SkyCoord("17:47:20.364",
                                                   "-28:23:08.792",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:19.886",
                                                 "-28:22:59.814",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'MandN',
                'bbox':[-0.43,0.425],
                'loc': 2,
                'l1':1,
                'l2':4,
                'min': -1,
                'max': 300,
                'zoom': 2.25,
                'inset_axes': 'M',
               },
              }
zoomregions_order = ['Mfull', 'M', 'N', 'M_inner', 'SouthOfSouth', 'MidDS', 'LowerDS']


filenames = {'continuum': contfilename,
             '1.3cm': '/Users/adam/work/sgrb2/continuumdata/SGRB2_1.3CM_J2000_realigned.fits',
            }


for regionname in ('full', ):#'MandN', 'DeepSouth', ):

    legloc = 'upper right' if regionname == 'MandN' else 'upper left'
    leg_bbox = [0.5, -0.10] if regionname == 'DeepSouth' else (1,1)

    vmax_hi = defaultdict(lambda: 0.25*1e3)
    vmax_hi['continuum'] = 0.01*1e3
    vmax_hi['1.3cm'] = 0.02*1e3
    vmin_hi = defaultdict(lambda: -0.0001*1e3)
    vmin_hi['continuum'] = -0.0005*1e3
    vmax_lo = defaultdict(lambda: 0.1*1e3)
    vmax_lo['continuum'] = 0.0025*1e3
    vmin_lo = defaultdict(lambda: -0.001*1e3)
    vmin_lo['continuum'] = -0.0002*1e3

    for line in ("continuum", "1.3cm"):#"HC3N",):
    #for line in ("1.3cm",):#"HC3N",):
    #for line in ("continuum",):#"HC3N",):


        if regionname == 'DeepSouth':
            # Deep South
            bottomleft = coordinates.SkyCoord("17:47:24.199", "-28:26:02.565", unit=(u.h, u.deg), frame='fk5')
            topright = coordinates.SkyCoord("17:47:17.666", "-28:23:30.722", unit=(u.h, u.deg), frame='fk5')
            scalebarpos = coordinates.SkyCoord("17:47:23.7", "-28:23:45.0", unit=(u.h, u.deg), frame='fk5')
        elif regionname == 'MandN':
            bottomleft = coordinates.SkyCoord("17:47:24.199", "-28:23:30.722", unit=(u.h, u.deg), frame='fk5')
            #if line == '1.3cm':
            #    topright = coordinates.SkyCoord("17:47:16.25", "-28:21:36", unit=(u.h, u.deg), frame='fk5')
            #    scalebarpos = coordinates.SkyCoord("17:47:18.0", "-28:23:25.0", unit=(u.h, u.deg), frame='fk5')
            #else:
            topright = coordinates.SkyCoord("17:47:14.666", "-28:21:04.980", unit=(u.h, u.deg), frame='fk5')
            scalebarpos = coordinates.SkyCoord("17:47:17.0", "-28:23:25.0", unit=(u.h, u.deg), frame='fk5')
        elif regionname == 'full':
            mylims_fk5 = ((266.8744477, -28.44946601), (266.7838311, -28.33589021))
            topright = coordinates.SkyCoord(*mylims_fk5[1],
                                            unit=(u.deg, u.deg),
                                            frame='fk5')
            bottomleft = coordinates.SkyCoord(*mylims_fk5[0],
                                              unit=(u.deg, u.deg),
                                              frame='fk5')
            scalebarpos = coordinates.SkyCoord("17:47:27", "-28:26:15.0",
                                               unit=(u.h, u.deg), frame='fk5')
        else:
            raise Exception


        if line in filenames:
            hdu_line = fits.open(filenames[line])[0]
        else:
            hdu_line = fits.open(paths.Fpath('merge/max/SgrB2_b3_7M_12M.{0}.image.pbcor_max_medsub.fits'.format(line)))[0]
        mywcs = wcs.WCS(hdu_line.header).celestial
        #wcsaxes = WCSaxes(mywcs.to_header())
        wcsaxes = mywcs # hack for newer version of wcsaxes that is included in astropy

        fig3 = pl.figure(3)
        fig3.clf()
        ax = fig3.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)

        ra = ax.coords['ra']
        ra.set_major_formatter('hh:mm:ss.s')
        dec = ax.coords['dec']
        ra.set_axislabel("RA (J2000)", fontsize=pl.rcParams['axes.labelsize'])
        dec.set_axislabel("Dec (J2000)", fontsize=pl.rcParams['axes.labelsize'], minpad=0.0)
        ra.ticklabels.set_fontsize(tick_fontsize)
        ra.set_ticks(exclude_overlapping=True)
        dec.ticklabels.set_fontsize(tick_fontsize)
        dec.set_ticks(exclude_overlapping=True)


        im = ax.imshow(hdu_line.data.squeeze()*1e3,
                       transform=ax.get_transform(mywcs),
                       vmin=vmin_hi[line], vmax=vmax_hi[line], cmap=pl.cm.gray_r,
                       interpolation='nearest',
                       origin='lower', norm=asinh_norm.AsinhNorm())
        tr_fk5 = ax.get_transform("fk5")
        #(x1,y1),(x2,y2) = (1200,434),(2142,1743)
        # wrong (x1,y1),(x2,y2) = tr_fk5.transform_point([bottomleft.ra.deg, bottomleft.dec.deg]),tr_fk5.transform_point([topright.ra.deg, topright.dec.deg])
        (x1,y1),(x2,y2) = (mywcs.wcs_world2pix([[bottomleft.ra.deg,
                                                 bottomleft.dec.deg]],0)[0],
                           mywcs.wcs_world2pix([[topright.ra.deg,
                                                 topright.dec.deg]],0)[0]
                          )

        make_scalebar(ax, scalebarpos,
                      length=(0.5*u.pc / distance).to(u.arcsec,
                                                      u.dimensionless_angles()),
                      color='k',
                      label='0.5 pc',
                      text_offset=1.0*u.arcsec,
                     )

        if regionname == 'full':

            markersize = 1
            coredots = plotcores(ax, alpha=1,
                                 transform=ax.get_transform('fk5'),
                                 markerfacecolor='none',
                                 markersize=markersize, zorder=50)

        ax.set_aspect(1)
        ax.axis([x1,x2,y1,y2])


        for zoomregion in zoomregions_order:

            ZR = zoomregions[zoomregion]
            if ZR['inregion'] != regionname:
                continue

            parent_ax = zoomregions[ZR['inset_axes']]['axins'] if 'inset_axes' in ZR else ax

            bl, tr = ZR['bottomleft'],ZR['topright'],
            (zx1,zy1),(zx2,zy2) = (mywcs.wcs_world2pix([[bl.ra.deg,
                                                         bl.dec.deg]],0)[0],
                                   mywcs.wcs_world2pix([[tr.ra.deg,
                                                         tr.dec.deg]],0)[0]
                                  )

            axins = zoomed_inset_axes(parent_ax, zoom=ZR['zoom'], loc=ZR['loc'],
                                      bbox_to_anchor=ZR['bbox'],
                                      bbox_transform=fig3.transFigure,
                                      axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                                      axes_kwargs=dict(wcs=wcsaxes))
            ZR['axins'] = axins
            imz = axins.imshow(hdu_line.data.squeeze()*1e3,
                               #transform=parent_ax.get_transform(mywcs),
                               vmin=ZR['min'], vmax=ZR['max'], cmap=pl.cm.gray_r,
                               interpolation='nearest',
                               origin='lower', norm=asinh_norm.AsinhNorm())


            coredots = plotcores(axins, alpha=1,
                                 transform=axins.get_transform('fk5'),
                                 markerfacecolor='none',
                                 markersize=markersize, zorder=50)
            ax.axis([x1,x2,y1,y2])
            axins.axis([zx1,zx2,zy1,zy2])

            axins.set_xticklabels([])
            axins.set_yticklabels([])

            lon = axins.coords['ra']
            lat = axins.coords['dec']
            lon.set_ticklabel_visible(False)
            lat.set_ticklabel_visible(False)

            # draw a bbox of the region of the inset axes in the parent axes and
            # connecting lines between the bbox and the inset axes area
            mark_inset(parent_axes=parent_ax, inset_axes=axins,
                       loc1=ZR['l1'], loc2=ZR['l2'], fc="none", ec="0.5")


            fig3.canvas.draw()
            assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

            cax = fig3.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                                 ax.bbox._bbox.y1-ax.bbox._bbox.y0])
            cb = fig3.colorbar(mappable=im, cax=cax)
            cb.set_label("$S_{3 mm}$ [mJy beam$^{-1}$]")


            print(("core_overlays/cores_on_{0}_peak_{1}_zoomin: {2}"
                   .format(line,regionname,zoomregion)))
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_zoomin.png"
                                 .format(line,regionname)),
                     bbox_inches='tight',
                     dpi=200
                    )
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_zoomin.pdf"
                                 .format(line,regionname)),
                     bbox_inches='tight')

        leg = ax.legend(handles=coredots, labels=[cd.get_label() for cd in
                                                  coredots],
                        bbox_to_anchor=leg_bbox,
                        loc=legloc,
                        fontsize=8)

        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_zoomin_legend.png"
                                 .format(line,regionname)),
                     bbox_inches='tight')
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_zoomin_legend.pdf"
                                 .format(line,regionname)),
                     bbox_inches='tight')
