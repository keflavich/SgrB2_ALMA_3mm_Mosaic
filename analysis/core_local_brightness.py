import os
import copy

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import table
from astropy import coordinates
from astropy import wcs
from astropy.convolution import convolve_fft,Gaussian2DKernel
import reproject
import pyregion

from constants import distance, mass_represented_by_a_source
from gutermuth2011_law import gas_depletion_law, sigma_gas_of_t
import lada2017relation

import paths

datapath = '/Users/adam/work/sgrb2/alma/FITS/continuumdata'

cara_higal_fit_scaling = 1e22


names = {'HerschelColumn25': 'Herschel (25\")',
         'HerschelColumn36': 'Herschel (36\")',
         'Sharc20Column': 'SHARC (20 K)',
         'Sharc50Column': 'SHARC (50 K)',
         'Scuba20Column': 'SCUBA (20 K)',
         'Scuba50Column': 'SCUBA (50 K)',
         'ScubaHTemColumn': 'SCUBA ($T=T_{Herschel}$)',
         'SharcHTemColumn': 'SHARC ($T=T_{Herschel}$)',
        }

files = {"SHARC": {'wavelength': 350*u.um, 'bmarea':9.55e-10*u.sr, 'bunit':u.Jy, 'filename':'SgrB2_350um_gal.fits',},
         "ATLASGAL": {'wavelength': 870*u.um, 'bmarea':4.91e-9*u.sr, 'bunit':u.Jy, 'filename':'AG-Laboca-Planck.1.5.fits',},
         "SCUBA450": {'wavelength': 450*u.um, 'bmarea':7.59e-10*u.sr, 'bunit':u.Jy, 'filename':'gc450_gal_gal_zoom.fits',},
         "SCUBA850": {'wavelength': 850*u.um, 'bmarea':2.71e-9*u.sr, 'bunit':u.Jy, 'filename':'gc850_gal_gal_zoom.fits',},
         #"Herschel70": {'wavelength': 350*u.um, 'bmarea':10*u.arcsec**2, 'bunit':u.Jy, 'filename':'igls_l000_blue_cnr_hh_gal_zoom.fits',},
         #"Herschel160": {'wavelength': 350*u.um, 'bmarea':10*u.arcsec**2, 'bunit':u.Jy, 'filename':'igls_l000_red_cnr_hh_gal_zoom.fits',},
         "Herschel500": {'wavelength': 500*u.um, 'bmarea':1.72e-8*u.sr, 'bunit':u.MJy/u.sr, 'filename':'igls_l000_plw_deglitch_hh_gal_zoom.fits',},
         "Herschel350": {'wavelength': 350*u.um, 'bmarea':8.43e-9*u.sr, 'bunit':u.MJy/u.sr, 'filename':'igls_l000_pmw_deglitch_hh_gal_zoom.fits',},
         "Herschel250": {'wavelength': 250*u.um, 'bmarea':4.30*u.sr, 'bunit':u.MJy/u.sr, 'filename':'igls_l000_psw_deglitch_hh_gal_zoom.fits',},
         "BGPS": {'wavelength': 1100*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.Jy, 'filename':'v2.1_ds2_l001_13pca_map20.fits',},
         "HerschelColumn25": {'wavelength': 0*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.cm**-2, 'filename':'gcmosaic_column_conv25.fits',},
         "HerschelColumn36": {'wavelength': 0*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.cm**-2, 'filename':'gcmosaic_column_conv36.fits',},
         "Sharc20Column": {'wavelength': 350*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.cm**-2, 'filename':'column_maps/sharc_col_20K.fits',},
         "Sharc50Column": {'wavelength': 350*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.cm**-2, 'filename':'column_maps/sharc_col_50K.fits',},
         "Scuba20Column": {'wavelength': 450*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.cm**-2, 'filename':'column_maps/scuba_col_20K.fits',},
         "Scuba50Column": {'wavelength': 450*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.cm**-2, 'filename':'column_maps/scuba_col_50K.fits',},
         "SharcHTemColumn": {'wavelength': 350*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.cm**-2, 'filename':'column_maps/sharc_col_herscheltem.fits',},
         "ScubaHTemColumn": {'wavelength': 450*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.cm**-2, 'filename':'column_maps/scuba_col_herscheltem.fits',},
        }

"""
Goal:
    For each core in the core photometry catalog, determine the corresponding
    brightness in the above images, then put the core CDF on top of the
    brightness ("column") PDF.
"""

alma_hdr = fits.Header(dict(NAXIS=2,
                            NAXIS1=450,
                            NAXIS2=450,
                            CTYPE1='RA---SIN',
                            CRVAL1=2.668301750000E+02,
                            CDELT1=-0.0002777777777777778,
                            CRPIX1=225.0,
                            CUNIT1='deg     ',
                            CTYPE2='DEC--SIN',
                            CRVAL2=-2.839256111111E+01,
                            CDELT2=0.0002777777777777778,
                            CRPIX2=225.0,
                            CUNIT2='deg     ',)
                      )


tbl = table.Table.read(paths.tpath("continuum_photometry_withSIMBAD.ipac"), format='ascii.ipac',)

regfile = pyregion.open(paths.rpath('coverage.reg'))

for imname in files:
    print("Loading {0}".format(imname))
    tbl.add_column(table.Column(name=imname, data=np.zeros(len(tbl),
                                                           dtype='float'),
                                unit=files[imname]['bunit']))

    ffile = fits.open(os.path.join(datapath, files[imname]['filename']))
    new_data,_ = reproject.reproject_interp(ffile, alma_hdr)
    new_hdu = fits.PrimaryHDU(data=new_data, header=alma_hdr)

    files[imname]['file'] = new_hdu

    ww = wcs.WCS(files[imname]['file'].header)
    files[imname]['wcs'] = ww

    files[imname]['pixarea'] = (ww.wcs.cdelt[1]*u.deg)**2

    files[imname]['mask'] = regfile.get_mask(hdu=new_hdu)

    assert ww.naxis == 2

# Fix column around Sgr B2
col = files['HerschelColumn36']['file'].data
col_conv = convolve_fft(col, Gaussian2DKernel(5), nan_treatment='interpolate',
                        normalize_kernel=True)
files['HerschelColumn36']['file'].data[np.isnan(col)] = col_conv[np.isnan(col)]
if np.nanmax(col) > 1e22:
    raise ValueError("The Herschel column file being used appears inconsistent with this code.")
files['HerschelColumn36']['file'].data *= cara_higal_fit_scaling
files['HerschelColumn36']['bunit'] = u.cm**-2

col = files['HerschelColumn25']['file'].data
col_conv = convolve_fft(col, Gaussian2DKernel(5), nan_treatment='interpolate',
                        normalize_kernel=True)
files['HerschelColumn25']['file'].data[np.isnan(col)] = col_conv[np.isnan(col)]
if np.nanmax(col) > 1e22:
    raise ValueError("The Herschel column file being used appears inconsistent with this code.")
files['HerschelColumn25']['file'].data *= cara_higal_fit_scaling
files['HerschelColumn25']['bunit'] = u.cm**-2


for row in tbl:

    coord = coordinates.SkyCoord(row['RA'], row['Dec'], frame='fk5', unit=(u.deg, u.deg))

    for imname,fmeta in files.items():
        if 'GLON' in fmeta['wcs'].wcs.ctype[0]:
            xcrd, ycrd = fmeta['wcs'].wcs_world2pix(coord.galactic.l, coord.galactic.b, 0)
        elif fmeta['wcs'].wcs.equinox == 1950:
            xcrd, ycrd = fmeta['wcs'].wcs_world2pix(coord.fk4.ra, coord.fk4.dec, 0)
        else: # should always be this one now...
            xcrd, ycrd = fmeta['wcs'].wcs_world2pix(coord.ra, coord.dec, 0)

        # TODO: convert to K
        # (todo is no longer necessary: just use column now)
        row[imname] = (fmeta['file'].data[int(ycrd), int(xcrd)] * fmeta['bunit']).value

tbl.write(paths.tpath("continuum_photometry_plusbackground.ipac"), format='ascii.ipac',
          overwrite=True)


brick_files = copy.deepcopy(files)
brick_hdr = fits.Header(dict(NAXIS=2,
                             NAXIS1=450,
                             NAXIS2=450,
                             CTYPE1='RA---SIN',
                             CRVAL1=266.5335253,
                             CDELT1=-0.0002777777777777778,
                             CRPIX1=225.0,
                             CUNIT1='deg     ',
                             CTYPE2='DEC--SIN',
                             CRVAL2=-28.70780832,
                             CDELT2=0.0002777777777777778,
                             CRPIX2=225.0,
                             CUNIT2='deg     ',)
                       )
brick_files['SHARC']['filename'] = 'SHARC_350_Dowell.fits'
brick_files['SCUBA450']['filename'] = 'gc450.fits.gz'
brick_files['SCUBA850']['filename'] = 'gc850.fits.gz'
brick_files['ATLASGAL']['filename'] = 'AG-Laboca-Planck.0.0.fits'
brick_files["Herschel500"]['filename'] = 'igls_l000_plw_deglitch_hh.fits'
brick_files["Herschel350"]['filename'] = 'igls_l000_pmw_deglitch_hh.fits'
brick_files["Herschel250"]['filename'] = 'igls_l000_psw_deglitch_hh.fits'
brick_files["BGPS"]['filename'] = 'v2.1_ds2_l000_13pca_map20_crop.fits'
brick_files['ScubaHTemColumn']['filename'] = 'column_maps/brick_scuba_col_herscheltem.fits'

for imname in brick_files:

    ffile = fits.open(os.path.join(datapath, brick_files[imname]['filename']))
    new_data,_ = reproject.reproject_interp(ffile, brick_hdr)
    new_hdu = fits.PrimaryHDU(data=new_data, header=brick_hdr)

    brick_files[imname]['file'] = new_hdu

    ww = wcs.WCS(brick_files[imname]['file'].header)
    brick_files[imname]['wcs'] = ww

brick_files['HerschelColumn25']['file'].data *= cara_higal_fit_scaling
if np.nanmax(brick_files['HerschelColumn25']['file'].data) > 1e30:
    raise ValueError("Inconsistent column density")
brick_files['HerschelColumn25']['bunit'] = u.cm**-2
brick_files['HerschelColumn36']['file'].data *= cara_higal_fit_scaling
if np.nanmax(brick_files['HerschelColumn36']['file'].data) > 1e30:
    raise ValueError("Inconsistent column density")
brick_files['HerschelColumn36']['bunit'] = u.cm**-2


def nan_to_val(x, val):
    y = x.copy()
    y[np.isnan(x)] = val
    return y


def plotit():
    import astropy.stats
    import pylab as pl
    pl.rcParams['axes.prop_cycle'] = pl.cycler('color', ['#1f77b4', '#ff7f0e',
                                                         '#2ca02c', '#d62728',
                                                         '#9467bd', '#8c564b',
                                                         '#e377c2', '#7f7f7f',
                                                         '#bcbd22', '#17becf'])
    pl.rcParams['figure.figsize'] = (12,8)
    pl.rcParams['figure.dpi'] = 75
    pl.rcParams['savefig.dpi'] = 150
    pl.rcParams['axes.titlesize']= 40
    pl.rcParams['axes.labelsize']= 24
    pl.rcParams['xtick.labelsize']= 20
    pl.rcParams['ytick.labelsize']= 20

    pl.figure(1).clf()
    pl.figure(2).clf()
    pl.figure(3).clf()
    pl.figure(4).clf()
    pl.figure(5, figsize=(12,8), dpi=75).clf()
    for ii,imname in enumerate(f for f in files if 'column' in f.lower()):


        print("Plotting {0}".format(imname))
        data = files[imname]['file'].data
        mask = files[imname]['mask']
        brickdata = brick_files[imname]['file'].data
        # https://github.com/astropy/astropy/pull/5232 for ignore_nan
        lo = astropy.stats.mad_std(data, ignore_nan=True)
        hi = np.nanmax(data)
        bins = np.logspace(np.log10(lo), np.log10(hi), 100)

        pl.figure(1)
        pl.subplot(2,4,ii+1)
        brickweights = np.ones(np.isfinite(brickdata).sum(),
                               dtype='float')/np.isfinite(brickdata).sum()*np.isfinite(data).sum()
        bH,bL,bP = pl.hist(brickdata[np.isfinite(brickdata)], bins=bins,
                           log=True, alpha=0.5, histtype='step',
                           #weights=brickweights,
                           bottom=1.1,
                           color='b', zorder=-1)
        #weights = np.ones(np.isfinite(data).sum(), dtype='float')/np.isfinite(data).sum()
        H,L,P = pl.hist(data[np.isfinite(data) & mask], bins=bins, log=True,
                        alpha=0.5, color='k',
                        #normed=True,
                        histtype='step')
        # Lada threshold, approximately (116 Msun/pc^2)
        pl.vlines(5e21, 1.1, H.max(), label='Lada+ 2010 Threshold')
        #pl.hist(tbl[imname], bins=bins, log=True, alpha=0.5)
        pl.xlim(np.min([bL.min(), L.min()]), L.max())
        pl.semilogx()
        pl.yscale('log', nonposy='clip')
        ax2 = pl.gca().twinx()
        ax2.plot(np.sort(tbl[imname]), np.arange(len(tbl),
                                                 dtype='float')/len(tbl),
                 'k-', linewidth=3, alpha=0.5, zorder=10)
        ax2.set_xlim(L.min(), L.max())
        ax2.set_ylim(0,1)
        pl.title(imname, fontsize=12)

        pl.figure(2)
        pl.subplot(2,4,ii+1)
        pl.imshow(data, cmap='gray_r')
        pl.contour(data, levels=np.percentile(tbl[imname],[5,10,25,50,75,90,95]))
        pl.contour(mask, levels=[0.5])
        pl.title(imname)

        sorted_col = u.Quantity(np.sort(data[mask & (data>0)]), u.cm**-2)
        cumul = np.cumsum(sorted_col)[::-1]
        pixarea_cm2 = (files[imname]['pixarea']*distance**2)
        cum_mass = (cumul * pixarea_cm2 * 2.8*u.Da).to(u.M_sun, u.dimensionless_angles())
        pl.figure(3)
        pl.plot(sorted_col, cum_mass, label=names[imname])
        pl.legend(loc='best')
        pl.xlim(1e21,2e25)
        pl.gca().set_xscale('log')
        pl.gca().set_yscale('log')

        pl.figure(4)
        assert cum_mass.min() > 0.1*u.M_sun
        pl.plot(sorted_col-5e22*u.cm**-2, cum_mass-(5e22*u.cm**-2*pixarea_cm2*2.8*u.Da).to(u.M_sun, u.dimensionless_angles()), label=names[imname])
        pl.legend(loc='best')
        pl.xlim(1e21,2e25)
        pl.gca().set_xscale('log')
        pl.gca().set_yscale('log')


    pl.figure(1)
    pl.tight_layout()
    pl.savefig(paths.fpath("flux_histograms_with_core_location_CDF.png"), bbox_inches='tight')


    pl.figure(3)
    pl.vlines(5e21, 0.5, 1e6, color='k', linestyle='--', linewidth=1,
              label='Lada+ 2010 Threshold')
    pl.tight_layout()
    pl.legend(loc='best')
    pl.savefig(paths.fpath("mass_cdf_histograms.png"), bbox_inches='tight')
    pl.figure(4)
    pl.vlines(5e21, 0.5, 1e6, color='k', linestyle='--', linewidth=1,
              label='Lada+ 2010 Threshold')
    pl.legend(loc='best')
    pl.tight_layout()
    pl.savefig(paths.fpath("mass_cdf_histograms_bgsubd.png"), bbox_inches='tight')

    pl.figure(5).clf()
    for imname,color in [('HerschelColumn25', 'k'), ('HerschelColumn36','g'),
                         ('SharcHTemColumn','b'), ('ScubaHTemColumn','r')]:
        pl.plot(np.sort(tbl[imname]), np.arange(len(tbl),
                                                dtype='float')/len(tbl),
                linestyle='-', linewidth=4, alpha=1, zorder=10,
                color=color,
                label=names[imname],
               )
        pl.xscale('log')
        pl.xlim(1e21,2e25)
        pl.ylim(0,1)
    fb1 = pl.fill_betweenx(y=np.arange(len(tbl), dtype='float')/len(tbl),
                           x1=np.sort(nan_to_val(tbl['Sharc20Column'], 1e26)),
                           x2=np.sort(nan_to_val(tbl['Sharc50Column'], 1e26)),
                           alpha=0.5,
                           label='SHARC 20-50 K',
                           edgecolor='b',
                           facecolor='none',
                           zorder=-10,
                           hatch='//',
                           linewidth=2,
                          )
    fb2 = pl.fill_betweenx(y=np.arange(len(tbl), dtype='float')/len(tbl),
                           x1=np.sort(nan_to_val(tbl['Scuba20Column'], 1e26)),
                           x2=np.sort(nan_to_val(tbl['Scuba50Column'], 1e26)),
                           alpha=0.5,
                           label='SCUBA 20-50 K',
                           edgecolor='r',
                           facecolor='none',
                           zorder=-10,
                           hatch='\\\\',
                           linewidth=2,
                          )


    pl.figure(5)
    pl.vlines(5e21, 0, 1, color='k', linestyle='--', linewidth=1,
              label='Lada+ 2010 Threshold')
    pl.vlines(2e23, 0, 1,
              color='k',
              linestyle=':',
              label='Krumholz+ 2008 Threshold')
    pl.legend(loc='best', fontsize=20)
    pl.tight_layout()
    pl.xlabel("Column Density [$N(\mathrm{H}_2)$ cm$^{-2}$]", fontsize=24)
    pl.ylabel("Cumulative fraction\nof cores at column $<N$", fontsize=24)
    ax1 = pl.gca()
    pl.draw()
    ax2 = ax1.twiny()
    print("ax1 xlims: {0}".format(ax1.get_xlim()))
    pl.draw()
    ax2.set_xlim(np.array(ax1.get_xlim())*(2.8*u.Da).to(u.g).value)
    print("ax2 xlims: {0}".format(ax2.get_xlim()))
    ax2.set_xscale('log')
    ax2.set_xlabel("Column Density [g cm$^{-2}$]")
    pl.draw()
    pl.savefig(paths.fpath("core_background_column_cdf.png"), bbox_inches='tight')

    #pl.figure(3)
    #pl.tight_layout()
    #pl.savefig(paths.fpath("cumulative_mass_histograms.png"), bbox_inches='tight')

    pl.figure(6).clf()
    ax1 = pl.gca()
    imname = 'ScubaHTemColumn'
    print("fig6: Plotting {0}".format(imname))
    data = files[imname]['file'].data
    mask = files[imname]['mask']
    brickdata = brick_files[imname]['file'].data
    lo = astropy.stats.mad_std(brickdata, ignore_nan=True)
    hi = np.nanmax(data)
    bins = np.logspace(np.log10(lo)-0.5, np.log10(hi), 100)
    bH,bL,bP = ax1.hist(brickdata[np.isfinite(brickdata)], bins=bins,
                        log=True, alpha=0.5, histtype='step',
                        bottom=0.1,
                        linewidth=2,
                        color='b', zorder=-1)
    #weights = np.ones(np.isfinite(data).sum(), dtype='float')/np.isfinite(data).sum()
    H,L,P = ax1.hist(data[np.isfinite(data) & mask], bins=bins, log=True,
                     alpha=0.5, color='k',
                     linewidth=2,
                     #normed=True,
                     histtype='step')
    # Lada threshold, approximately (116 Msun/pc^2)
    ax1.vlines(5e21, 1.1, H.max(),
               label='Lada+ 2010 Threshold')
    ax1.vlines(2e23, 0.1, H.max()*2,
               color='k',
               linestyle=':',
               label='Krumholz+ 2008 Threshold')
    #ax1.hist(tbl[imname], bins=bins, log=True, alpha=0.5)
    ax1.set_xlim(np.min([bL.min(), L.min()]), L.max())
    ax1.set_ylim(0.5,np.max([bH.max(), H.max()])*1.1)
    ax1.semilogx()
    ax1.set_yscale('log', nonposy='clip')
    ax1.set_xlabel("Column Density [N(H$_2$) cm$^{-2}$]")
    ax1.set_ylabel("Number of pixels")

    ax3 = ax1.twiny()
    print("ax1 xlims: {0}".format(ax1.get_xlim()))
    pl.draw()
    ax3.set_xlim(np.array(ax1.get_xlim())*(2.8*u.Da).to(u.g).value)
    ax3lims = ax3.get_xlim()
    print("ax2 xlims: {0}".format(ax2.get_xlim()))
    ax3.set_xscale('log')
    ax3.set_xlabel("Column Density [g cm$^{-2}$]")
    pl.draw()

    pl.savefig(paths.fpath("compare_brick_sgrb2_colPDF_nofractions.png"), bbox_inches='tight')

    ax2 = ax1.twinx()
    ax2.plot(np.sort(tbl[imname]), np.arange(len(tbl),
                                             dtype='float')/len(tbl),
             'k-', linewidth=3, alpha=0.5, zorder=10)
    ax2.set_xlim(L.min(), L.max())
    ax2.set_ylim(0,1)

    ax2.set_ylabel("Fraction of point sources below N(H$_2$)")
    ax3.set_xlim(ax3lims)
    ax2.set_xlim(L.min(), L.max())

    pl.savefig(paths.fpath("compare_brick_sgrb2_colPDF.png"), bbox_inches='tight')


    nn11_pc = (u.Quantity(tbl['nn11'], u.arcsec) * distance).to(u.pc,
                                                                u.dimensionless_angles())
    nn = 11
    nn11_msunpersqpc_starcentered = ((nn-1) * mass_represented_by_a_source /
                                     (np.pi*(nn11_pc)**2)).to(u.M_sun/u.pc**2)

    herschelsurfdens = (u.Quantity(tbl['HerschelColumn25']).to(u.cm**-2) *
                        2.8*u.Da).to(u.M_sun/u.pc**2)

    assert herschelsurfdens.min() < 10**5*u.M_sun/u.pc**2
    assert herschelsurfdens.max() > 10**3*u.M_sun/u.pc**2

    fig7 = pl.figure(7, figsize=(10,10))
    fig7.clf()
    ax7 = fig7.gca()
    ax7.loglog(
               herschelsurfdens,
               nn11_msunpersqpc_starcentered,
               'k.', alpha=0.7, markeredgecolor=(0,0,0,0.5))
    lims = ax7.axis()
    # 5/3 slope
    #ax7.loglog([1e3,1e6], [3e0, 3e5], 'k--')


    ax7.set_ylabel("Source-centered NN11 Surface Density\n$\Sigma_*$ [M$_\odot$ pc$^{-2}$]")
    ax7.set_xlabel("Gas Surface Density $\Sigma_{gas}$ [M$_\odot$ pc$^{-2}$]")
    ax7.axis(lims)

    # Arrows showing the shift if you subtract off the "most aggressive plausible"
    # uniform foreground value
    bg_5e22 = (5e22*u.cm**-2*2.8*u.Da).to(u.M_sun/u.pc**2)
    ax7.arrow(3e3, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(3e3-bg_5e22.value), color='k')
    ax7.arrow(1e4, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(1e4-bg_5e22.value), color='k')
    ax7.arrow(3e4, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(3e4-bg_5e22.value), color='k')

    ax7.axis([1e3,1e5,1e0,1e5])
    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel_nomodels_nolocal.png"), bbox_inches='tight')
    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel_nomodels_nolocal.pdf"), bbox_inches='tight')


    monr2_lowerline = np.array([0.1, 1e5])**2.67/(100**2.67) * 2.5
    monr2_fill = ax7.fill_between([0.1, 1e5],
                                  monr2_lowerline,
                                  monr2_lowerline*10,
                                  alpha=0.5,
                                  color='green',
                                  label='Mon R2')
    oph_lowerline = np.array([0.1, 1e5])**1.87/(100**1.87) * 1.5
    oph_fill = ax7.fill_between([0.1, 1e5],
                                oph_lowerline,
                                oph_lowerline*10,
                                color='blue',
                                alpha=0.5,
                                label='Ophiucus')

    local_plotobjs = [monr2_fill, oph_fill]


    T = ax7.text(2.3e3, 9e3, "Mon R2", color='k', rotation=50,
                 fontsize=18,
                 verticalalignment='bottom', horizontalalignment='center')
    local_plotobjs.append(T)
    T = ax7.text(2.5e3, 4.5e2, "Ophiucus", color='k', rotation=38,
                 fontsize=18,
                 verticalalignment='bottom', horizontalalignment='center')
    local_plotobjs.append(T)

    #ax7.plot([0.1, 1e5], np.array([0.1, 1e5])**1.87/(1e4**1.87)*(1e4**(5/3.)/1e5), 'b:', linewidth=3, alpha=0.5)
    oph_scalefactor = 50.
    L = ax7.plot([0.1, 1e5], oph_lowerline/oph_scalefactor, 'b:', linewidth=3, alpha=0.5)
    local_plotobjs += L


    ax7.axis([1e3,1e5,1e0,1e5])
    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel_nomodels.png"), bbox_inches='tight')
    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel_nomodels.pdf"), bbox_inches='tight')

    sigma_gas = np.logspace(1,6) * u.M_sun / u.pc**2

    mdlplots = {}

    for time in (0.01, 0.1, 0.74)*u.Myr:
        mdlplots[(1, time)] = ax7.loglog(sigma_gas_of_t(sigma_gas, time, alpha=1, k=0.1/u.Myr),
                                         gas_depletion_law(sigma_gas, time, alpha=1, k=0.1/u.Myr), label=time,
                                         color='r', linewidth=3, alpha=0.5, zorder=-10,)
        mdlplots[(2,time)] = ax7.loglog(sigma_gas_of_t(sigma_gas, time),
                                        gas_depletion_law(sigma_gas, time), label=time,
                                        color='orange', linewidth=3, zorder=-5, alpha=0.5)


    ax7.axis([1e3,1e5,1e0,1e5])
    #ax7.plot([0.1, 1e5], np.array([0.1, 1e5])*4e-2, 'r-', linewidth=3, alpha=0.5, zorder=-10)
    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel.png"), bbox_inches='tight')
    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel.pdf"), bbox_inches='tight')

    ax7.loglog(np.logspace(3,5),
               lada2017relation.sigma_star_california(np.logspace(3,5)*u.M_sun/u.pc**2),
               linewidth=3, color='m', label='Lada2017_cali')
    ax7.loglog(np.logspace(3,5),
               lada2017relation.sigma_star_orionA(np.logspace(3,5)*u.M_sun/u.pc**2),
               linewidth=3, linestyle='--', color='m', label='Lada2017_orionA')
    ax7.loglog(np.logspace(3,5),
               lada2017relation.sigma_star_orionB(np.logspace(3,5)*u.M_sun/u.pc**2),
               linewidth=3, linestyle=':', color='m', label='Lada2017_orionB')

    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel_withLada2017.png"), bbox_inches='tight')
    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel_withLada2017.pdf"), bbox_inches='tight')

    for line in mdlplots.values():
        try:
            for ll in line:
                ll.set_visible(False)
        except TypeError:
            line.set_visible(False)
    for obj in local_plotobjs:
        obj.set_visible(False)

    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel_nomodels_withLada2017.png"), bbox_inches='tight')
    fig7.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_herschel_nomodels_withLada2017.pdf"), bbox_inches='tight')


    scubasurfdens = (u.Quantity(tbl['ScubaHTemColumn'], u.cm**-2) * 2.8*u.Da).to(u.M_sun/u.pc**2)


    fig8 = pl.figure(8)
    fig8.clf()
    ax8 = fig8.gca()
    ax8.loglog(
               scubasurfdens,
               nn11_msunpersqpc_starcentered,
               'k.', alpha=0.7, markeredgecolor=(0,0,0,0.5))
    lims = ax8.axis()
    # 5/3 slope
    #ax8.loglog([1e3,1e6], [3e0, 3e5], 'k--')

    #ax8.plot([0.1, 1e5], np.array([0.1, 1e5])**1.87/(1e4**1.87)*(1e4**(5/3.)/1e5), 'b:', linewidth=3, alpha=0.5)
    ax8.plot([0.1, 1e5], oph_lowerline/oph_scalefactor, 'b:', linewidth=3, alpha=0.5)
    ax8.axis(lims)
    ax8.set_ylabel("Source-centered NN11 Surface Density\n$\Sigma_*$ [M$_\odot$ pc$^{-2}$]")
    ax8.set_xlabel("SCUBA-derived Surface Density [M$_\odot$ pc$^{-2}$]")

    # arrows showing the shift if you subtract off the "most aggressive plausible"
    # uniform foreground value
    bg_5e22 = (5e22*u.cm**-2*2.8*u.Da).to(u.M_sun/u.pc**2)
    ax8.arrow(3e3, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(3e3-bg_5e22.value), color='k')
    ax8.arrow(1e4, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(1e4-bg_5e22.value), color='k')
    ax8.arrow(3e4, 1.1, -bg_5e22.value, 0, head_width=0.1, head_length=0.033*(3e4-bg_5e22.value), color='k')

    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba_nomodels_nolocal.png"), bbox_inches='tight')
    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba_nomodels_nolocal.pdf"), bbox_inches='tight')

    monr2_lowerline = np.array([0.1, 1e5])**2.67/(100**2.67) * 2.5
    monr2_fill = ax8.fill_between([0.1, 1e5],
                                  monr2_lowerline,
                                  monr2_lowerline*10,
                                  alpha=0.5,
                                  color='green',
                                  label='Mon R2')
    oph_lowerline = np.array([0.1, 1e5])**1.87/(100**1.87) * 1.5
    oph_fill = ax8.fill_between([0.1, 1e5],
                                oph_lowerline,
                                oph_lowerline*10,
                                color='blue',
                                alpha=0.5,
                                label='Ophiucus')

    local_plotobjs = [monr2_fill, oph_fill]

    monr2txt = ax8.text(2.3e3, 9e3, "Mon R2", color='k', rotation=50,
                        fontsize=18, verticalalignment='bottom',
                        horizontalalignment='center')
    ophtxt = ax8.text(2.5e3, 4.5e2, "Ophiucus", color='k', rotation=38,
                      fontsize=18, verticalalignment='bottom',
                      horizontalalignment='center')
    local_plotobjs.append(monr2txt)
    local_plotobjs.append(ophtxt)

    ax8.axis([1e3,1e5,1e0,1e5])
    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba_nomodels.png"), bbox_inches='tight')
    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba_nomodels.pdf"), bbox_inches='tight')


    mdlplots = {}

    for time in (0.01, 0.1, 0.74)*u.Myr:
        mdlplots[(1, time)] = ax8.loglog(sigma_gas_of_t(sigma_gas, time, alpha=1, k=0.1/u.Myr),
                                         gas_depletion_law(sigma_gas, time, alpha=1, k=0.1/u.Myr), label=time,
                                         color='r', linewidth=3, alpha=0.5, zorder=-10,)
        mdlplots[(2,time)] = ax8.loglog(sigma_gas_of_t(sigma_gas, time),
                                        gas_depletion_law(sigma_gas, time), label=time,
                                        color='orange', linewidth=3, zorder=-5, alpha=0.5)


    #ax8.axis(lims)
    ax8.axis([1e3,1e5,1e0,1e5])
    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba.png"), bbox_inches='tight')
    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba.pdf"), bbox_inches='tight')

    ax8.loglog(np.logspace(3,5),
               lada2017relation.sigma_star(np.logspace(3,5)*u.M_sun/u.pc**2),
               linewidth=3, color='m', label='Lada2017')
    ax8.loglog(np.logspace(3,5),
               lada2017relation.sigma_star_orionA(np.logspace(3,5)*u.M_sun/u.pc**2),
               linewidth=3, linestyle='--', color='m', label='Lada2017_orionA')
    ax8.loglog(np.logspace(3,5),
               lada2017relation.sigma_star_orionB(np.logspace(3,5)*u.M_sun/u.pc**2),
               linewidth=3, linestyle=':', color='m', label='Lada2017_orionB')

    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba_withLada2017.png"), bbox_inches='tight')
    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba_withLada2017.pdf"), bbox_inches='tight')

    for line in mdlplots.values():
        try:
            for ll in line:
                ll.set_visible(False)
        except TypeError:
            line.set_visible(False)
    for obj in local_plotobjs:
        obj.set_visible(False)

    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba_nomodels_withLada2017.png"), bbox_inches='tight')
    fig8.savefig(paths.fpath("stellar_vs_gas_column_density_starcentered_scuba_nomodels_withLada2017.pdf"), bbox_inches='tight')

    # TODO: plot the same (?) histograms for The Brick
    # DONE!

if __name__ == "__main__":
    import matplotlib
    matplotlib.rcParams['figure.figsize'] = (12,8)
    matplotlib.rcParams['figure.dpi'] = 75.
    matplotlib.rcParams['savefig.dpi'] = 300.
    matplotlib.rcParams['axes.labelsize'] = 9
    matplotlib.rcParams['xtick.labelsize'] = 8
    matplotlib.rcParams['ytick.labelsize'] = 8
    matplotlib.use('Qt5Agg')
    import pylab as pl
    pl.ion()

    plotit()
