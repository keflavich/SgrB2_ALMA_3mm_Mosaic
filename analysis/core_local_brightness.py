import os

import numpy as np 
from astropy.io import fits
from astropy import units as u
from astropy import table
from astropy import coordinates
from astropy import wcs
import reproject

import paths

datapath = '/Users/adam/work/sgrb2/continuumdata'


files = {"SHARC": {'wavelength': 350*u.um, 'bmarea':9.55e-10*u.sr, 'bunit':u.Jy, 'filename':'SgrB2_350um_gal.fits',},
         "ATLASGAL": {'wavelength': 870*u.um, 'bmarea':4.91e-9*u.sr, 'bunit':u.Jy, 'filename':'ATLASGAL.1.5_reprojectoverlap9.fits',},
         "SCUBA450": {'wavelength': 450*u.um, 'bmarea':7.59e-10*u.sr, 'bunit':u.Jy, 'filename':'gc450_gal_gal_zoom.fits',},
         "SCUBA850": {'wavelength': 850*u.um, 'bmarea':2.71e-9*u.sr, 'bunit':u.Jy, 'filename':'gc850_gal_gal_zoom.fits',},
         #"Herschel70": {'wavelength': 350*u.um, 'bmarea':10*u.arcsec**2, 'bunit':u.Jy, 'filename':'igls_l000_blue_cnr_hh_gal_zoom.fits',},
         #"Herschel160": {'wavelength': 350*u.um, 'bmarea':10*u.arcsec**2, 'bunit':u.Jy, 'filename':'igls_l000_red_cnr_hh_gal_zoom.fits',},
         "Herschel500": {'wavelength': 500*u.um, 'bmarea':1.72e-8*u.sr, 'bunit':u.MJy/u.sr, 'filename':'igls_l000_plw_deglitch_hh_gal_zoom.fits',},
         "Herschel350": {'wavelength': 350*u.um, 'bmarea':8.43e-9*u.sr, 'bunit':u.MJy/u.sr, 'filename':'igls_l000_pmw_deglitch_hh_gal_zoom.fits',},
         "Herschel250": {'wavelength': 250*u.um, 'bmarea':4.30*u.sr, 'bunit':u.MJy/u.sr, 'filename':'igls_l000_psw_deglitch_hh_gal_zoom.fits',},
         "BGPS": {'wavelength': 1100*u.um, 'bmarea':2.90e-8*u.sr, 'bunit':u.Jy, 'filename':'igls_l000_psw_deglitch_hh_gal_zoom.fits',},
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

tbl = table.Table.read(paths.tpath("continuum_photometry.ipac"), format='ascii.ipac',)

for imname in files:
    tbl.add_column(table.Column(name=imname, data=np.zeros(len(tbl),
                                                           dtype='float'),
                                unit=files[imname]['bunit']))

    ffile = fits.open(os.path.join(datapath, files[imname]['filename']))
    new_data,_ = reproject.reproject_interp(ffile, alma_hdr)
    new_hdu = fits.PrimaryHDU(data=new_data, header=alma_hdr)

    files[imname]['file'] = new_hdu

    ww = wcs.WCS(files[imname]['file'].header)
    files[imname]['wcs'] = ww

    assert ww.naxis == 2

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
        row[imname] = (fmeta['file'].data[int(ycrd), int(xcrd)] * fmeta['bunit']).value

tbl.write(paths.tpath("continuum_photometry_plusbackground.ipac"), format='ascii.ipac',
          overwrite=True)

def plotit():
    import astropy.stats
    import pylab as pl

    pl.figure(1).clf()
    pl.figure(2).clf()
    for ii,imname in enumerate(files):

        data = files[imname]['file'].data
        # https://github.com/astropy/astropy/pull/5232 for ignore_nan
        lo = astropy.stats.mad_std(data, ignore_nan=True)
        hi = np.nanmax(data)
        bins = np.logspace(np.log10(lo), np.log10(hi), 100)

        pl.figure(1)
        pl.subplot(3,3,ii+1)
        H,L,P = pl.hist(data[np.isfinite(data)], bins=bins, log=True, alpha=0.5, normed=True)
        #pl.hist(tbl[imname], bins=bins, log=True, alpha=0.5)
        pl.xlim(L.min(), L.max())
        pl.semilogx()
        ax2 = pl.gca().twinx()
        ax2.plot(sorted(tbl[imname]), np.arange(len(tbl),
                                                dtype='float')/len(tbl)*H.max(),
                 'k-', linewidth=3, alpha=0.5, zorder=10)
        ax2.set_xlim(L.min(), L.max())
        pl.title(imname)

        pl.figure(2)
        pl.subplot(3,3,ii+1)
        pl.imshow(data, cmap='gray_r')
        pl.contour(data, levels=np.percentile(tbl[imname],[5,10,25,50,75,90,95]))
        pl.title(imname)

    # TODO: plot the same (?) histograms for The Brick
