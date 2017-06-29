"""
This is obsolete, even though it was used.  tp_7m_12m_combine.py is an improved
and more efficient version of this.
"""
import numpy as np
from spectral_cube import SpectralCube
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
from uvcombine import spectral_regrid, feather_simple, fourier_combine_cubes
from astropy import units as u
from astropy import log

log.setLevel("DEBUG")
log.debug('Reading...')
cube = SpectralCube.read('SgrB2_b3_7M_12M.HC3N.image.pbcor_medsub.fits').with_spectral_unit(u.GHz)
cube.allow_huge_operations=True
log.debug("Read file")
cube_k = cube.to(u.K, cube.beam.jtok_equiv(cube.spectral_axis))
log.debug("Converted first file to K")


tpcube = SpectralCube.read('../tp/tp_concat.spw17.image.fits').with_spectral_unit(u.GHz)
crop_freqs = cube_k.spectral_axis[0], cube_k.spectral_axis[-1]
crop_channels = sorted((tpcube.closest_spectral_channel(crop_freqs[0]), tpcube.closest_spectral_channel(crop_freqs[1])))
tpcube = tpcube[crop_channels[0]-1:crop_channels[1]+1]
log.debug("Read tp freq")
tpcube_k = tpcube.to(u.K, tpcube.beam.jtok_equiv(tpcube.spectral_axis))
log.debug("Converted TP to K")
# determine smooth factor
kw = cube_k.spectral_axis.diff().mean() / tpcube_k.spectral_axis.diff().mean()
log.debug("determined kernel")
tpcube_k_smooth = tpcube_k.spectral_smooth(Gaussian1DKernel(kw/(8*np.log(2))**0.5))
#tpcube_k_smooth = FITS_tools.cube_regrid.spectral_smooth_cube(tpcube_k,
#                                                              kw.value/np.sqrt(8*np.log(2)))
log.debug("completed cube smooth")

tpcube_k_ds = tpcube_k_smooth[::2,:,:]
tpcube_k_ds.hdu.writeto('HC3N_tp_freq_ds.fits', clobber=True)
log.debug("wrote kelvin tp cube")

tpcube_k_ds_rg = tpcube_k_ds.spectral_interpolate(cube_k.spectral_axis)
#tpkrg = spectral_regrid(tpcube_k_ds, cube_k.spectral_axis)
log.debug("done regridding")
#tpkrg.writeto('HC3N_tp_freq_ds_interp.fits', clobber=True)
tpcube_k_ds_rg.write('HC3N_tp_freq_ds_interp.fits', overwrite=True)
#
#cube_tpkrg = SpectralCube.read('HC3N_tp_freq_ds_interp.fits')

frq = 90.9662*u.GHz
im = cube_k[cube_k.closest_spectral_channel(frq)]
sdim = tpcube_k_ds_rg[tpcube_k_ds_rg.closest_spectral_channel(frq)]
im.write('cubek_ch126.fits', overwrite=True)
sdim.write('tpcube_k_rg_ch126.fits', overwrite=True)
combohdu, hdu2 = feather_simple(im.hdu, sdim.hdu, return_regridded_lores=True, return_hdu=True)
combohdu.writeto('HC3N_TP_7m_12m_feather_126.fits', clobber=True)

# final goal
combhdu = fourier_combine_cubes(cube_k, tpcube_k_ds_rg, return_hdu=True,
                                lowresfwhm=tpcube_k_ds_rg.beam.major,
                                maximum_cube_size=3e8)
combhdu.header.update(cube_k.header)
combhdu.header.update(cube_k.beam.to_header_keywords())
combcube = SpectralCube.read(combhdu).with_spectral_unit(u.km/u.s,
                                                         velocity_convention='radio')
combcube.write('HC3N_TP_7m_12m_feather.fits', overwrite=True)
