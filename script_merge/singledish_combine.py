import image_registration
from image_registration.fft_tools.zoom import zoom_on_pixel
from FITS_tools.cube_regrid import regrid_fits_cube,regrid_cube_hdu
from FITS_tools.hcongrid import hcongrid,hcongrid_hdu
import FITS_tools
import fft_psd_tools
from astropy import wcs
from astropy.io import fits
from astropy import coordinates
from astropy import units as u
from astropy import log
from astropy.utils.console import ProgressBar
import numpy as np
from spectral_cube import SpectralCube

def fourier_combine_cubes(cube1, cube2, highresextnum=0,
                          highresscalefactor=1.0,
                          lowresscalefactor=1.0, lowresfwhm=1*u.arcmin,
                          return_regridded_cube2=False,
                          return_hdu=False,
                         ):
    """
    Fourier combine two data cubes

    Parameters
    ----------
    cube1 : SpectralCube
    highresfitsfile : str
        The high-resolution FITS file
    cube2 : SpectralCube
    lowresfitsfile : str
        The low-resolution (single-dish) FITS file
    highresextnum : int
        The extension number to use from the high-res FITS file
    highresscalefactor : float
    lowresscalefactor : float
        A factor to multiply the high- or low-resolution data by to match the
        low- or high-resolution data
    lowresfwhm : `astropy.units.Quantity`
        The full-width-half-max of the single-dish (low-resolution) beam;
        or the scale at which you want to try to match the low/high resolution
        data
    return_hdu : bool
        Return an HDU instead of just a cube.  It will contain two image
        planes, one for the real and one for the imaginary data.
    return_regridded_cube2 : bool
        Return the 2nd cube regridded into the pixel space of the first?
    """
    if isinstance(cube1, str):
        cube1 = SpectralCube.read(cube1)
    if isinstance(cube2, str):
        cube2 = SpectralCube.read(cube2)
    #cube1 = spectral_cube.io.fits.load_fits_cube(highresfitsfile,
    #                                             hdu=highresextnum)
    im1 = cube1._data # want the raw data for this
    hd1 = cube1.header
    assert hd1['NAXIS'] == im1.ndim == 3
    w1 = cube1.wcs
    pixscale = np.abs(w1.wcs.get_cdelt()[0]) # REPLACE EVENTUALLY...

    cube2 = cube2.to(cube1.unit)

    assert cube1.unit == cube2.unit, 'Cubes must have same or equivalent unit'
    assert cube1.unit.is_equivalent(u.Jy/u.beam) or cube1.unit.is_equivalent(u.K), "Cubes must have brightness units."

    #f2 = regrid_fits_cube(lowresfitsfile, hd1)
    f2 = regrid_cube_hdu(cube2.hdu, hd1)
    w2 = wcs.WCS(f2.header)

    nax1,nax2,nax3 = (hd1['NAXIS1'],
                      hd1['NAXIS2'],
                      hd1['NAXIS3'])

    dcube1 = im1 * highresscalefactor
    dcube2 = f2.data * lowresscalefactor
    outcube = np.empty_like(dcube1)

    xgrid,ygrid = (np.indices([nax2,nax1])-np.array([(nax2-1.)/2,(nax1-1.)/2.])[:,None,None])
    fwhm = np.sqrt(8*np.log(2))
    # sigma in pixels
    sigma = ((lowresfwhm/fwhm/(pixscale*u.deg)).decompose().value)
    #sigma_fftspace = (1/(4*np.pi**2*sigma**2))**0.5
    sigma_fftspace = (2*np.pi*sigma)**-1
    log.debug('sigma = {0}, sigma_fftspace={1}'.format(sigma, sigma_fftspace))

    kernel = np.fft.fftshift(np.exp(-(xgrid**2+ygrid**2)/(2*sigma**2)))
    # convert the kernel, which is just a gaussian in image space,
    # to its corresponding kernel in fourier space
    kfft = np.abs(np.fft.fft2(kernel)) # should be mostly real
    # normalize the kernel
    kfft/=kfft.max()
    ikfft = 1-kfft

    pb = ProgressBar(dcube1.shape[0])

    for ii,(im1,im2) in enumerate(zip(dcube1, dcube2)):

        fft1 = np.fft.fft2(np.nan_to_num(im1))
        fft2 = np.fft.fft2(np.nan_to_num(im2))

        fftsum = kfft*fft2 + ikfft*fft1

        combo = np.fft.ifft2(fftsum)
        outcube[ii,:,:] = combo.real

        pb.update(ii+1)

    if return_regridded_cube2:
        return outcube, f2
    elif return_hdu:
        return fits.PrimaryHDU(data=outcube, header=w1.to_header())
    else:
        return outcube

def fourier_combine(highresfitsfile, lowresfitsfile,
                    matching_scale=60*u.arcsec, scale=False,
                    return_hdu=False,
                   ):
    """
    Simple reimplementation of 'feather' for 2D images
    """
    f1 = fits.open(highresfitsfile)
    w1 = wcs.WCS(f1[0].header)
    f2 = fits.open(lowresfitsfile)
    w2 = wcs.WCS(f2[0].header)

    nax1,nax2 = f1[0].header['NAXIS1'], f1[0].header['NAXIS2']
    # We take care of zooming later...
    #if not(nax1 == f2[0].header['NAXIS1'] and nax2 == f2[0].header['NAXIS2']):
    #    raise ValueError("Images are not in the same pixel space; reproject "
    #                     "them to common pixel space first.")

    pixscale1 = w1.wcs.get_cdelt()[1]
    pixscale2 = w2.wcs.get_cdelt()[1]

    center = w1.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([nax1/2.],
                                                          [nax2/2.],
                                                          1)
    frame = 'icrs' if w1.celestial.wcs.ctype[0][:2] == 'RA' else 'galactic'
    if w2.celestial.wcs.ctype[0][:2] == 'RA':
        center = coordinates.SkyCoord(*(center*u.deg), frame=frame).fk5
        cxy = center.ra.deg,center.dec.deg
    elif w2.celestial.wcs.ctype[0][:4] == 'GLON':
        center = coordinates.SkyCoord(*(center*u.deg), frame=frame).galactic
        cxy = center.l.deg,center.b.deg

    im1 = f1[0].data.squeeze()
    im1[np.isnan(im1)] = 0
    shape = im1.shape
    im2raw = f2[0].data.squeeze()
    im2raw[np.isnan(im2raw)] = 0
    if len(shape) != im2raw.ndim:
        raise ValueError("Different # of dimensions in the interferometer and "
                         "single-dish images")
    if len(shape) == 3:
        if shape[0] != im2raw.shape[0]:
            raise ValueError("Spectral dimensions of cubes do not match.")

    center_pixel = w2.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix(cxy[0], cxy[1],
                                                                0)[::-1]

    zoomed = zoom_on_pixel(np.nan_to_num(im2raw),
                           center_pixel,
                           usfac=np.abs(pixscale2/pixscale1),
                           outshape=shape)

    im2 = zoomed

    xax,psd1 = fft_psd_tools.PSD2(im1, oned=True)
    xax,psd2 = fft_psd_tools.PSD2(im2, oned=True)

    xax_as = (pixscale1/xax*u.deg).to(u.arcsec)

    if scale:
        closest_point = np.argmin(np.abs(xax_as-matching_scale))

        scale_2to1 = (psd1[closest_point] / psd2[closest_point])**0.5
    else:
        scale_2to1 = 1

    fft1 = np.fft.fft2(im1)
    fft2 = np.fft.fft2(im2) * scale_2to1

    xgrid,ygrid = (np.indices(shape)-np.array([(shape[0]-1.)/2,(shape[1]-1.)/2.])[:,None,None])
    
    sigma = np.abs(shape[0]/((matching_scale/(pixscale1*u.deg)).decompose().value)) / np.sqrt(8*np.log(2))
    kernel = np.fft.fftshift(np.exp(-(xgrid**2+ygrid**2)/(2*sigma**2)))
    kernel/=kernel.max()

    fftsum = kernel*fft2 + (1-kernel)*fft1

    combo = np.fft.ifft2(fftsum)

    if not return_hdu:
        return combo
    elif return_hdu:
        combo_hdu = fits.PrimaryHDU(data=np.abs(combo), header=w1.to_header())
        return combo_hdu


def feather_simple(hires, lores,
                   highresextnum=0,
                   lowresextnum=0,
                   highresscalefactor=1.0,
                   lowresscalefactor=1.0, lowresfwhm=1*u.arcmin,
                   return_hdu=False,
                   return_regridded_lores=False):
    """
    Fourier combine two data cubes

    Parameters
    ----------
    highresfitsfile : str
        The high-resolution FITS file
    lowresfitsfile : str
        The low-resolution (single-dish) FITS file
    highresextnum : int
        The extension number to use from the high-res FITS file
    highresscalefactor : float
    lowresscalefactor : float
        A factor to multiply the high- or low-resolution data by to match the
        low- or high-resolution data
    lowresfwhm : `astropy.units.Quantity`
        The full-width-half-max of the single-dish (low-resolution) beam;
        or the scale at which you want to try to match the low/high resolution
        data
    return_hdu : bool
        Return an HDU instead of just an image.  It will contain two image
        planes, one for the real and one for the imaginary data.
    return_regridded_cube2 : bool
        Return the 2nd cube regridded into the pixel space of the first?
    """
    if isinstance(hires, (fits.ImageHDU, fits.PrimaryHDU)):
        hdu1 = hires
    else:
        hdu1 = fits.open(hires)[highresextnum]
    if isinstance(lores, (fits.ImageHDU, fits.PrimaryHDU)):
        hdu2 = lores
    else:
        hdu2 = fits.open(lores)[lowresextnum]
    im1 = hdu1.data.squeeze()
    hd1 = FITS_tools.strip_headers.flatten_header(hdu1.header)
    assert hd1['NAXIS'] == im1.ndim == 2
    pixscale = FITS_tools.header_tools.header_to_platescale(hd1)
    log.debug('pixscale = {0}'.format(pixscale))

    im2raw = hdu2.data.squeeze()
    hd2 = FITS_tools.strip_headers.flatten_header(hdu2.header)
    assert hd2['NAXIS'] == im2raw.ndim == 2
    hdu2 = fits.PrimaryHDU(data=im2raw, header=hd2)

    hdu2 = hcongrid_hdu(hdu2, hd1)
    im2 = hdu2.data.squeeze()

    nax1,nax2 = (hd1['NAXIS1'],
                 hd1['NAXIS2'],
                )

    ygrid,xgrid = (np.indices([nax2,nax1])-np.array([(nax2-1.)/2,(nax1-1.)/2.])[:,None,None])
    fwhm = np.sqrt(8*np.log(2))
    # sigma in pixels
    sigma = ((lowresfwhm/fwhm/(pixscale*u.deg)).decompose().value)
    #sigma_fftspace = (1/(4*np.pi**2*sigma**2))**0.5
    sigma_fftspace = (2*np.pi*sigma)**-1
    log.debug('sigma = {0}, sigma_fftspace={1}'.format(sigma, sigma_fftspace))

    kernel = np.fft.fftshift(np.exp(-(xgrid**2+ygrid**2)/(2*sigma**2)))
    # convert the kernel, which is just a gaussian in image space,
    # to its corresponding kernel in fourier space
    kfft = np.abs(np.fft.fft2(kernel)) # should be mostly real
    # normalize the kernel
    kfft/=kfft.max()
    ikfft = 1-kfft

    fft1 = np.fft.fft2(np.nan_to_num(im1*highresscalefactor))
    fft2 = np.fft.fft2(np.nan_to_num(im2*lowresscalefactor))

    fftsum = kfft*fft2 + ikfft*fft1

    combo = np.fft.ifft2(fftsum)

    if return_hdu:
        combo_hdu = fits.PrimaryHDU(data=combo.real, header=hdu1.header)
        combo = combo_hdu

    if return_regridded_lores:
        return combo, hdu2
    else:
        return combo

def spectral_regrid(cube, outgrid):
    """
    Spectrally regrid a cube onto a new spectral output grid

    (this is apparently redundant with regrid_cube_hdu)
    """

    assert isinstance(cube, SpectralCube)

    inaxis = cube.spectral_axis.to(outgrid.unit)

    indiff = np.mean(np.diff(inaxis))
    outdiff = np.mean(np.diff(outgrid))
    if outdiff < 0:
        outgrid=outgrid[::-1]
        outdiff = np.mean(np.diff(outgrid))
    if indiff < 0:
        cubedata = cube.filled_data[::-1]
        inaxis = cube.spectral_axis.to(outgrid.unit)[::-1]
        indiff = np.mean(np.diff(inaxis))
    else:
        cubedata = cube.filled_data[:]
    if indiff < 0 or outdiff < 0:
        raise ValueError("impossible.")

    assert np.all(np.diff(outgrid) > 0)
    assert np.all(np.diff(inaxis) > 0)

    np.testing.assert_allclose(np.diff(outgrid), outdiff,
                               err_msg="Output grid must be linear")

    if outdiff > 2 * indiff:
        raise ValueError("Input grid has too small a spacing.  It needs to be "
                         "smoothed prior to resampling.")

    newcube = np.empty([outgrid.size, cube.shape[1], cube.shape[2]])

    yy,xx = np.indices(cube.shape[1:])

    pb = ProgressBar(xx.size)
    for ix, iy in (zip(xx.flat, yy.flat)):
        newcube[:,iy,ix] = np.interp(outgrid.value, inaxis.value,
                                     cubedata[:,iy,ix].value)
        pb.update()

    newheader = cube.header
    newheader['CRPIX3'] = 1
    newheader['CRVAL3'] = outgrid[0].value
    newheader['CDELT3'] = outdiff.value
    newheader['CUNIT3'] = outgrid.unit.to_string('FITS')

    return fits.PrimaryHDU(data=newcube, header=newheader)

