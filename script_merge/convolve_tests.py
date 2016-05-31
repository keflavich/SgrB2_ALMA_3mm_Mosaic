import numpy as np
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from astropy.io import fits
import radio_beam
from astropy import units as u
from astropy import wcs

tpfn = 'singleframes/HC3N_natural_tpcube_k_rg_65kms.fits'
featherfn = 'singleframes/HC3N_natural_TP_7m_12m_feather_65kms.fits'
intfn = 'singleframes/HC3N_natural_cubek_65kms.fits'

intfh = fits.open(intfn)
featherfh = fits.open(featherfn)
tpfh = fits.open(tpfn)

beam_int = radio_beam.Beam.from_fits_header(intfh[0].header)
beam_feath = radio_beam.Beam.from_fits_header(featherfh[0].header)
beam_tp = radio_beam.Beam.from_fits_header(tpfh[0].header)


pixscale = wcs.utils.proj_plane_pixel_area(wcs.WCS(intfh[0].header))**0.5 * u.deg
kernel = Gaussian2DKernel(beam_tp.major / (8*np.log(2))**0.5 / pixscale)
int_conv = convolve_fft(intfh[0].data, kernel)
feather_conv = convolve_fft(featherfh[0].data, kernel)

featherfh[0].data = feather_conv
featherfh.writeto("singleframes/_test_HC3N_feather_convolved.fits", clobber=True)
intfh[0].data = int_conv
intfh.writeto("singleframes/_test_HC3N_int_convolved.fits", clobber=True)
