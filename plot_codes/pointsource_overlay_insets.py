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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import TransformedBbox, BboxPatch, BboxConnector
from matplotlib.transforms import Bbox

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

filenames = {'continuum': contfilename,
             '1.3cm': '/Users/adam/work/sgrb2/continuumdata/SGRB2_1.3CM_J2000_realigned.fits',
            }


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
                                                   "-28:23:14.0",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:18.8",
                                                 "-28:22:51.5",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'full',
                'bbox': [-0.15,0.75],
                #'bbox':[-0.024,0.425],
                'loc': 2,
                'l1':1,
                'l2':4,
                'min': -1,
                'max': 50,
                'zoom': 9,
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
               'Nhires':
               {'bottomleft': coordinates.SkyCoord("17:47:20.072",
                                                   "-28:22:22.0",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("17:47:19.706",
                                                 "-28:22:15.6",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'fullN',
                'fitsfile':paths.lbpath('sgr_b2m.N.B6.allspw.continuum.r0.5.clean1000.image.tt0.pbcor.fits'),
                'bbox':[-0.15,0.95],
                'loc': 2,
                'l1':1,
                'l2':4,
                'min': -1, # mJy
                'max': 25,
                'width': 2.5,
                'height': 4,
                'scalebarpos': coordinates.SkyCoord("17:47:20.004", "-28:22:20.9",
                                               unit=(u.h, u.deg), frame='icrs')

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
               'LowerDS_1mm':
               {'bottomleft': coordinates.SkyCoord("17:47:23.63",
                                                   "-28:25:58",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:19.772",
                                                 "-28:25:35.949",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'DeepestSouth',
                'fitsfile':paths.dspath('downsampled_minimal_member.uid___A001_X1290_X46.Sgr_B2_DS_sci.spw25_27_29_31.mfs.I.manual.image.pbcor.fits'),
                'bbox':[0.1,0.33],
                'loc': 2,
                'l1':1,
                'l2':2,
                'min': -0.2,
                'max': 5.0,
                'zoom': 4,
                'width': 6,
                'height': 4,
               },
               'LowerDS_1mm_hires':
               {'bottomleft': coordinates.SkyCoord("17:47:21.848",
                                                   "-28:25:57",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:20.623",
                                                 "-28:25:40.7",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'DeepestSouth',
                'inset_axes': 'LowerDS_1mm',
                'fitsfile':paths.dspath('downsampled_minimal_member.uid___A001_X1290_X46.Sgr_B2_DS_sci.spw25_27_29_31.mfs.I.manual.image.pbcor.fits'),
                'bbox':[0.1,-1.25],
                'loc': 3,
                'l1':1,
                'l2':2,
                'min': -0.2,
                'max': 1.0,
                'zoom': 4,
                'width': 4,
                'height': 4,
               },
               'LowerDS_1mm_higherres':
               {'bottomleft': coordinates.SkyCoord("17:47:21.358",
                                                   "-28:25:55.262",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:21.008",
                                                 "-28:25:51.605",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'DeepestSouth',
                'inset_axes': 'LowerDS_1mm_hires',
                'fitsfile':paths.dspath('downsampled_minimal_member.uid___A001_X1290_X46.Sgr_B2_DS_sci.spw25_27_29_31.mfs.I.manual.image.pbcor.fits'),
                'bbox':[0.75,-0.55],
                'contour':{'fn':filenames['continuum'],
                           'kwargs': {'levels': np.array([0.25, 0.5, 1, 1.5])/1e3,
                                      'colors': ['g']*5},
                          },
                'loc': 2,
                'l1':2,
                'l2':3,
                'min': -0.2,
                'max': 2.5,
                'zoom': 2,
                'width': 4,
                'height': 3,
               },
               'LowerDS_1mm_infullwide':
               {'bottomleft': coordinates.SkyCoord("17:47:24.274",
                                                   "-28:26:05.011",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:18.026",
                                                 "-28:23:27.725",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'inregion': 'fullwide',
                #'inset_axes': 'LowerDS_1mm_hires',
                'fitsfile':paths.dspath('downsampled_minimal_member.uid___A001_X1290_X46.Sgr_B2_DS_sci.spw25_27_29_31.mfs.I.manual.image.pbcor.fits'),
                'bbox':[0.00,0.68],
                #'contour':{'fn':filenames['continuum'],
                #           'kwargs': {'levels': np.array([0.25, 0.5, 1, 1.5])/1e3,
                #                      'colors': ['g']*5},
                #          },
                'loc': 2,
                'l1':1,
                'l2':4,
                'min': -0.2,
                'max': 2.5,
                'zoom': 2,
                'width': 5,
                'height': 8,
                'hide_axes': True,
                'show_cores': False,
               },
               'Nhires_fullwide':
               {'bottomleft': coordinates.SkyCoord("17:47:20.072",
                                                   "-28:22:22.0",
                                                   unit=(u.h, u.deg),
                                                   frame='icrs'),
                'topright': coordinates.SkyCoord("17:47:19.706",
                                                 "-28:22:15.6",
                                                 unit=(u.h, u.deg),
                                                 frame='icrs'),
                'inregion': 'fullwide',
                'fitsfile':paths.lbpath('sgr_b2m.N.B6.allspw.continuum.r0.5.clean1000.image.tt0.pbcor.fits'),
                'bbox':[0.75,0.75],
                'loc': 2,
                'l1':2,
                'l2':3,
                'min': -1, # mJy
                'max': 25,
                'width': 4,
                'height': 4,
                'hide_axes': True,
                'show_cores': False,
                #'scalebarpos': coordinates.SkyCoord("17:47:20.004", "-28:22:20.9",
                #                               unit=(u.h, u.deg), frame='icrs')

               },
               'M_inner_fullwide':
               {'bottomleft': coordinates.SkyCoord("17:47:20.364",
                                                   "-28:23:08.792",
                                                   unit=(u.h, u.deg),
                                                   frame='fk5'),
                'topright': coordinates.SkyCoord("17:47:19.886",
                                                 "-28:22:59.814",
                                                 unit=(u.h, u.deg),
                                                 frame='fk5'),
                'fitsfile': paths.lbpath('sgr_b2m.M.B6.allspw.continuum.r0.5.clean1000.image.tt0.pbcor.fits'),
                'inregion': 'fullwide',
                'bbox':[0.75,0.43],
                'loc': 2,
                'l1':2,
                'l2':3,
                'min': -1,
                'max': 100,
                'zoom': 5,
                'width': 4,
                'height': 4,
                'hide_axes': True,
                'show_cores': False,
                #'inset_axes': 'M',
               },
               #'LowerDS_1mm_higherres':
               #{'bottomleft': coordinates.SkyCoord("17:47:21.78",
               #                                    "-28:25:45.871",
               #                                    unit=(u.h, u.deg),
               #                                    frame='fk5'),
               # 'topright': coordinates.SkyCoord("17:47:21.442",
               #                                  "-28:25:42.000",
               #                                  unit=(u.h, u.deg),
               #                                  frame='fk5'),
               # 'inregion': 'DeepestSouth',
               # 'inset_axes': 'LowerDS_1mm_hires',
               # 'fitsfile':paths.dspath('downsampled_minimal_member.uid___A001_X1290_X46.Sgr_B2_DS_sci.spw0_1_2_3.mfs.I.manual.image.pbcor.fits'),
               # 'bbox':[0.7,-0.55],
               # 'loc': 2,
               # 'l1':1,
               # 'l2':3,
               # 'min': -0.2,
               # 'max': 2.5,
               # 'zoom': 2,
               # 'width': 4,
               # 'height': 3,
               #},
              }


zoomregions_order = ['LowerDS_1mm_infullwide', 'Nhires_fullwide', 'M_inner_fullwide',
                     'LowerDS_1mm', 'LowerDS_1mm_hires', 'LowerDS_1mm_higherres', 'Nhires', 'Mfull',
                     'M', 'N', 'M_inner', 'SouthOfSouth',
                     'MidDS', 'LowerDS']
#zoomregions_order = zoomregions_order[:2]

bigregion_parameters = {'fullwide': {'hide_axes': True, 'dpi':600, 'figsize':(12,17)},
                        'DeepestSouth': {'dpi':200, 'figsize':None},
                        'fullN': {'dpi':200, 'figsize':None},
                        'full': {'dpi':200, 'figsize':None},
                        'MandN': {'dpi':200, 'figsize':None},
                        'DeepSouth': {'dpi':200, 'figsize':None},
                       }


for regionname in ('fullwide',):#('DeepestSouth',): #'fullN', ):#'full', ):#'MandN', 'DeepSouth', ):
#for regionname in ('MandN',): #'fullN', ):#'full', ):#'MandN', 'DeepSouth', ):

    legloc = 'upper right' if regionname == 'MandN' else 'upper left'
    leg_bbox = [0.5, -0.10] if regionname in ('DeepSouth', 'DeepestSouth') else (1,1)

    vmax_hi = defaultdict(lambda: 0.25*1e3)
    vmax_hi['continuum'] = 0.01*1e3
    vmax_hi['1.3cm'] = 0.02*1e3
    vmin_hi = defaultdict(lambda: -0.0001*1e3)
    vmin_hi['continuum'] = -0.0005*1e3
    vmax_lo = defaultdict(lambda: 0.1*1e3)
    vmax_lo['continuum'] = 0.0025*1e3
    vmin_lo = defaultdict(lambda: -0.001*1e3)
    vmin_lo['continuum'] = -0.0002*1e3

    for line in ("continuum", ):#"1.3cm"):#"HC3N",):
    #for line in ("1.3cm",):#"HC3N",):
    #for line in ("continuum",):#"HC3N",):


        if regionname == 'DeepSouth':
            # Deep South
            bottomleft = coordinates.SkyCoord("17:47:24.199", "-28:26:02.565", unit=(u.h, u.deg), frame='fk5')
            topright = coordinates.SkyCoord("17:47:17.666", "-28:23:30.722", unit=(u.h, u.deg), frame='fk5')
            scalebarpos = coordinates.SkyCoord("17:47:23.7", "-28:23:45.0", unit=(u.h, u.deg), frame='fk5')
        elif regionname == 'DeepestSouth':
            # Deep South
            bottomleft = coordinates.SkyCoord("17:47:24.199", "-28:26:02.565", unit=(u.h, u.deg), frame='fk5')
            topright = coordinates.SkyCoord("17:47:17.666", "-28:25:20.722", unit=(u.h, u.deg), frame='fk5')
            scalebarpos = coordinates.SkyCoord("17:47:18.5", "-28:26:00.0", unit=(u.h, u.deg), frame='fk5')
        elif regionname == 'MandN':
            bottomleft = coordinates.SkyCoord("17:47:24.199", "-28:23:30.722", unit=(u.h, u.deg), frame='fk5')
            #if line == '1.3cm':
            #    topright = coordinates.SkyCoord("17:47:16.25", "-28:21:36", unit=(u.h, u.deg), frame='fk5')
            #    scalebarpos = coordinates.SkyCoord("17:47:18.0", "-28:23:25.0", unit=(u.h, u.deg), frame='fk5')
            #else:
            topright = coordinates.SkyCoord("17:47:14.666", "-28:21:04.980", unit=(u.h, u.deg), frame='fk5')
            scalebarpos = coordinates.SkyCoord("17:47:17.0", "-28:23:25.0", unit=(u.h, u.deg), frame='fk5')
        elif regionname in ('fullwide',):
            mylims_fk5 = ((266.8989023, -28.44894567), (266.7627384, -28.33548893))
            topright = coordinates.SkyCoord(*mylims_fk5[1],
                                            unit=(u.deg, u.deg),
                                            frame='fk5')
            bottomleft = coordinates.SkyCoord(*mylims_fk5[0],
                                              unit=(u.deg, u.deg),
                                              frame='fk5')
            scalebarpos = coordinates.SkyCoord("17:47:11", "-28:26:15.0",
                                               unit=(u.h, u.deg), frame='fk5')
        elif regionname in ('full', 'fullN'):
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
        toplevel_wcs = wcs.WCS(hdu_line.header).celestial
        #wcsaxes = WCSaxes(mywcs.to_header())
        #toplevel_wcs = mywcs # hack for newer version of toplevel_wcs that is included in astropy

        fig3 = pl.figure(3, figsize=bigregion_parameters[regionname]['figsize'])
        fig3.clf()
        ax = fig3.add_axes([0.15, 0.1, 0.8, 0.8], projection=toplevel_wcs)

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
                       transform=ax.get_transform(toplevel_wcs),
                       vmin=vmin_hi[line], vmax=vmax_hi[line], cmap=pl.cm.gray_r,
                       interpolation='nearest',
                       origin='lower', norm=asinh_norm.AsinhNorm())
        tr_fk5 = ax.get_transform("fk5")
        #(x1,y1),(x2,y2) = (1200,434),(2142,1743)
        # wrong (x1,y1),(x2,y2) = tr_fk5.transform_point([bottomleft.ra.deg, bottomleft.dec.deg]),tr_fk5.transform_point([topright.ra.deg, topright.dec.deg])
        (x1,y1),(x2,y2) = (toplevel_wcs.wcs_world2pix([[bottomleft.ra.deg,
                                                 bottomleft.dec.deg]],0)[0],
                           toplevel_wcs.wcs_world2pix([[topright.ra.deg,
                                                 topright.dec.deg]],0)[0]
                          )

        make_scalebar(ax, scalebarpos,
                      length=(0.5*u.pc / distance).to(u.arcsec,
                                                      u.dimensionless_angles()),
                      color='k',
                      label='0.5 pc',
                      text_offset=1.0*u.arcsec,
                     )

        if regionname in ('full', 'fullN'):

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
                print(f"{zoomregion} is not in {regionname}, it's in {ZR['inregion']}")
                continue

            if 'fitsfile' in ZR:
                inset_hdu = fits.open(ZR['fitsfile'])[0]
                ff_zr_wcs = wcs.WCS(inset_hdu.header).celestial
                # this was commented out, but it might be needed?
                if ff_zr_wcs.wcs.radesys != toplevel_wcs.wcs.radesys:
                    cntr = coordinates.SkyCoord(ff_zr_wcs.wcs.crval[0]*u.deg,
                                                ff_zr_wcs.wcs.crval[1]*u.deg,
                                                frame=ff_zr_wcs.wcs.radesys.lower())
                    cntr_new = cntr.transform_to(toplevel_wcs.wcs.radesys.lower())
                    ff_zr_wcs.wcs.radesys = toplevel_wcs.wcs.radesys
                    ff_zr_wcs.wcs.crval = [cntr_new.ra.deg,
                                       cntr_new.dec.deg]

                zr_wcs = ff_zr_wcs
            else:
                inset_hdu = hdu_line
                zr_wcs = toplevel_wcs

            parent_ax = zoomregions[ZR['inset_axes']]['axins'] if 'inset_axes' in ZR else ax

            bl, tr = ZR['bottomleft'],ZR['topright'],
            (zx1,zy1),(zx2,zy2) = (zr_wcs.wcs_world2pix([[bl.ra.deg,
                                                         bl.dec.deg]],0)[0],
                                   zr_wcs.wcs_world2pix([[tr.ra.deg,
                                                         tr.dec.deg]],0)[0]
                                  )

            # cut out the inset image before plotting to reduce time
            data = inset_hdu.data.squeeze()[int(zy1)-5:int(zy2)+5, int(zx1)-5:int(zx2)+5]
            zr_wcs = zr_wcs[int(zy1)-5:int(zy2)+5, int(zx1)-5:int(zx2)+5]

            if 'fitsfile' in ZR:
                axins = inset_axes(parent_ax,
                                   loc=ZR['loc'],
                                   width=ZR['width'],
                                   height=ZR['height'],
                                   bbox_to_anchor=ZR['bbox'],
                                   bbox_transform=fig3.transFigure,
                                   axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                                   axes_kwargs=dict(wcs=zr_wcs))
            else:
                axins = zoomed_inset_axes(parent_ax, zoom=ZR['zoom'], loc=ZR['loc'],
                                          bbox_to_anchor=ZR['bbox'],
                                          bbox_transform=fig3.transFigure,
                                          axes_class=astropy.visualization.wcsaxes.core.WCSAxes,
                                          axes_kwargs=dict(wcs=toplevel_wcs))

            ZR['axins'] = axins
            imz = axins.imshow(data*1e3,
                               #transform=parent_ax.get_transform(mywcs),
                               vmin=ZR['min'], vmax=ZR['max'], cmap=pl.cm.gray_r,
                               interpolation='nearest',
                               origin='lower', norm=asinh_norm.AsinhNorm())


            if 'show_cores' not in ZR or ZR['show_cores']:
                coredots = plotcores(axins, alpha=1,
                                     transform=axins.get_transform('fk5'),
                                     markerfacecolor='none',
                                     markersize=markersize, zorder=50)
            ax.axis([x1,x2,y1,y2])
            #axins.axis([zx1,zx2,zy1,zy2])
            # no need for this; should already be zoomed
            #axins.axis([0,zx2-zx1,0,zy2-zy1])

            axins.set_xticklabels([])
            axins.set_yticklabels([])

            lon = axins.coords['ra']
            lat = axins.coords['dec']
            lon.set_ticklabel_visible(False)
            lat.set_ticklabel_visible(False)

            if 'hide_axes' in ZR and ZR['hide_axes']:
                axins.xaxis.set_visible(False)
                axins.yaxis.set_visible(False)
                axins.xaxis.set_ticks([])
                axins.yaxis.set_ticks([])
                axins.set_facecolor('none')
                axins.axis('off')
                for spine in axins.spines.values():
                    spine.set_visible(False)
                for crd in axins.coords:
                    crd.set_ticklabel_visible(False)
                    crd.set_axislabel('')
                    crd.set_ticks_visible(False)
                    #crd.frame.set_visible(False)


            if 'fitsfile' not in ZR:
                # draw a bbox of the region of the inset axes in the parent axes and
                # connecting lines between the bbox and the inset axes area
                mark_inset(parent_axes=parent_ax, inset_axes=axins,
                           loc1=ZR['l1'], loc2=ZR['l2'], fc="none", ec="0.5")
            else: # if 'fitsfile' in ZR
                blt = bl.transform_to(parent_ax.wcs.wcs.radesys.lower())
                trt = tr.transform_to(parent_ax.wcs.wcs.radesys.lower())
                (rx1,ry1),(rx2,ry2) = (parent_ax.wcs.wcs_world2pix([[blt.ra.deg,
                                                               blt.dec.deg]],0)[0],
                                       parent_ax.wcs.wcs_world2pix([[trt.ra.deg,
                                                               trt.dec.deg]],0)[0]
                                      )
                bbox = Bbox(np.array([(rx1,ry1),(rx2,ry2)]))
                rect = TransformedBbox(bbox, parent_ax.transData)

                markinkwargs = dict(fc='none', ec='0.5')

                pp = BboxPatch(rect, fill=False, **markinkwargs)
                parent_ax.add_patch(pp)

                p1 = BboxConnector(axins.bbox, rect, loc1=ZR['l1'], **markinkwargs)
                axins.add_patch(p1)
                p1.set_clip_on(False)
                p2 = BboxConnector(axins.bbox, rect, loc1=ZR['l2'], **markinkwargs)
                axins.add_patch(p2)
                p2.set_clip_on(False)

                if 'scalebarpos' in ZR:
                    make_scalebar(axins, ZR['scalebarpos'],
                                  length=(5000*u.au / distance).to(u.arcsec,
                                                                   u.dimensionless_angles()),
                                  color='k',
                                  label='5000 AU',
                                  text_offset=0.08*u.arcsec,
                                 )

            if 'contour' in ZR:
                cfh = fits.open(ZR['contour']['fn'])
                cont_wcs = wcs.WCS(cfh[0].header)


                blt = bl.transform_to(cont_wcs.wcs.radesys.lower())
                trt = tr.transform_to(cont_wcs.wcs.radesys.lower())
                (rx1,ry1),(rx2,ry2) = (cont_wcs.wcs_world2pix([[blt.ra.deg,
                                                               blt.dec.deg]],0)[0],
                                       cont_wcs.wcs_world2pix([[trt.ra.deg,
                                                               trt.dec.deg]],0)[0]
                                      )

                cutout_dat = cfh[0].data[int(ry1)-5:int(ry2)+5, int(rx1)-5:int(rx2)+5]
                cutout_wcs = cont_wcs[int(ry1)-5:int(ry2)+5, int(rx1)-5:int(rx2)+5]

                cont_transform = axins.get_transform(cutout_wcs)

                axins.contour(cutout_dat, transform=cont_transform,
                              **ZR['contour']['kwargs'])


            fig3.canvas.draw()

            #assert np.abs(ax.bbox._bbox.x1 - 0.95) > 1e-4

            if 'hide_axes' in bigregion_parameters[regionname] and bigregion_parameters[regionname]['hide_axes']:
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.xaxis.set_visible(False)
                ax.yaxis.set_visible(False)
                ra.set_ticks_visible(False)
                dec.set_ticks_visible(False)
                ra.set_axislabel('')
                dec.set_axislabel('')
                ra.ticklabels.set_visible(False)
                dec.ticklabels.set_visible(False)
            else:
                cax = fig3.add_axes([ax.bbox._bbox.x1+0.01, ax.bbox._bbox.y0, 0.02,
                                     ax.bbox._bbox.y1-ax.bbox._bbox.y0])
                cb = fig3.colorbar(mappable=im, cax=cax)
                cb.set_label("$S_{3 mm}$ [mJy beam$^{-1}$]")


            print(("core_overlays/cores_on_{0}_peak_{1}_zoomin: {2}"
                   .format(line,regionname,zoomregion)))
        print("Saving figures")
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_zoomin.png"
                                 .format(line, regionname, )),
                     bbox_inches='tight',
                     dpi=bigregion_parameters[regionname]['dpi']
                    )
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_zoomin.pdf"
                                 .format(line, regionname, )),
                     bbox_inches='tight')

        if 'coredots' in locals():
            leg = ax.legend(handles=coredots, labels=[cd.get_label() for cd in
                                                      coredots],
                            bbox_to_anchor=leg_bbox,
                            loc=legloc,
                            fontsize=8)

        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_zoomin_legend.png"
                                 .format(line,regionname, )),
                     bbox_inches='tight')
        fig3.savefig(paths.fpath("core_overlays/cores_on_{0}_peak_{1}_zoomin_legend.pdf"
                                 .format(line,regionname, )),
                     bbox_inches='tight')
