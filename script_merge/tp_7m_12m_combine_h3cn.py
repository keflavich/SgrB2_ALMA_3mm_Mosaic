import numpy as np
from spectral_cube import SpectralCube
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel
from singledish_combine import spectral_regrid, feather_simple, fourier_combine_cubes
from astropy import units as u
from astropy import log

log.debug('Reading...')
cube = SpectralCube.read('SgrB2_b3_7M_12M.HC3N.image.pbcor_medsub.fits').with_spectral_unit(u.GHz)
cube.allow_huge_operations=True
log.debug("Read file")
cube_k = cube.to(u.K, cube.beam.jtok_equiv(cube.spectral_axis[:,None,None]))
log.debug("Converted first file to K")


tpcube = SpectralCube.read('../tp/tp_concat.spw17.image.fits').with_spectral_unit(u.GHz)
crop_freqs = cube_k.spectral_axis[0], cube_k.spectral_axis[-1]
crop_channels = sorted((tpcube.closest_spectral_channel(crop_freqs[0]), tpcube.closest_spectral_channel(crop_freqs[1])))
tpcube = tpcube[crop_channels[0]-1:crop_channels[1]+1]
log.debug("Read tp freq")
tpcube_k = tpcube.to(u.K, tpcube.beam.jtok_equiv(tpcube.spectral_axis[:,None,None]))
log.debug("Converted TP to K")
# determine smooth factor
kw = cube_k.spectral_axis.diff().mean() / tpcube_k.spectral_axis.diff().mean()
log.debug("determined kernel")
tbcube_k_smooth = tpcube_k.spectral_smooth(Gaussian1DKernel(kw/(8*np.log(2))**0.5))
#tbcube_k_smooth = FITS_tools.cube_regrid.spectral_smooth_cube(tpcube_k,
#                                                              kw.value/np.sqrt(8*np.log(2)))
log.debug("completed cube smooth")

tpcube_k_ds = tbcube_k_smooth[::2,:,:]
log.debug("downsampled")
log.debug("time to sleep")
print(tpcube_k)
print(tpcube_k.hdu)
tpcube_k.hdu
log.debug("did nothing")
tpcube_k_ds_hdu = tpcube_k.hdu
log.debug("made hdu")
tpcube_k_ds_hdu.data = tpcube_k_ds
log.debug("put data in hdu")
tpcube_k_ds_hdu.header['CRPIX3'] = 1
# why min? because we're forcing CDELT3 to be positive, therefore the 0'th channel
# must be the reference value.  Since we're using a symmetric kernel to downsample,
# the reference channel - wherever we pick it - must stay fixed.
tpcube_k_ds_hdu.header['CRVAL3'] = tpcube_k.spectral_axis[0].to(u.Hz).value
tpcube_k_ds_hdu.header['CUNIT3'] = tpcube_k.spectral_axis[0].to(u.Hz).unit.to_string('FITS')
tpcube_k_ds_hdu.header['CDELT3'] = tpcube_k.wcs.wcs.cdelt[2] * 2.
log.debug("completed header making")
tpcube_k_ds_hdu.writeto('HC3N_tp_freq_ds.fits', clobber=True)
log.debug("wrote kelvin tp cube")

tpdscube = SpectralCube.read('HC3N_tp_freq_ds.fits')
tpkrg = spectral_regrid(tpdscube,
                        cube_k.spectral_axis)
log.debug("done regridding")
tpkrg.writeto('HC3N_tp_freq_ds_interp.fits', clobber=True)

cube_tpkrg = SpectralCube.read('HC3N_tp_freq_ds_interp.fits')

frq = 90.9662*u.GHz
im = cube_k[cube_k.closest_spectral_channel(frq)]
sdim = cube_tpkrg[cube_tpkrg.closest_spectral_channel(frq)]
im.write('cubek_ch126.fits', overwrite=True)
sdim.write('tpcube_k_rg_ch126.fits', overwrite=True)
combohdu, hdu2 = feather_simple(im.hdu, sdim.hdu, return_regridded_lores=True, return_hdu=True)
combohdu.writeto('HC3N_TP_7m_12m_feather_126.fits', clobber=True)

# final goal
combhdu = fourier_combine_cubes(cube_k, cube_tpkrg, return_hdu=True,
                                lowresfwhm=cube_tpkrg.beam.major)
combcube = SpectralCube.read(combhdu).with_spectral_unit(u.km/u.s, velocity_convention='radio')
combcube.write('HC3N_TP_7m_12m_feather.fits', overwrite=True)
