import numpy as np
import radio_beam
from astropy import units as u
from astropy import wcs
from astropy.io import fits
import higal_sedfitter
import paths
import reproject
import dust_emissivity
from astropy import convolution
import pylab as pl

sharc_fn = paths.Fpath('other/SHARC_Herschel_Feathered.fits')
scuba_fn = paths.Fpath('other/scuba_Herschel_Feathered.fits')
scuba_hdr = fits.getheader(scuba_fn)

pixscale = wcs.WCS(scuba_hdr).celestial.pixel_scale_matrix[1,1]*u.deg

sharc_data,_ = reproject.reproject_interp(sharc_fn, scuba_hdr)
scuba_data = fits.getdata(scuba_fn)
ratio_map_sharc_to_scuba = sharc_data/scuba_data

# Column density is independent of pixel scale; it is related to surface
# brightness
imagecube = u.Quantity([sharc_data,scuba_data], u.MJy/u.sr)


if False:
    # don't need to fit it myself
    pixelfitter = higal_sedfitter.fit.PixelFitter(bfixed=True)

    tem,beta,col = higal_sedfitter.fit.fit_modified_blackbody_to_imagecube(
        imagecube, outheader=scuba_hdr, wavelengths=[350,450],
        pixelfitter=pixelfitter,
        out_prefix=paths.cpath('column_maps/sharcscuba_'),
    )

temperatures = np.linspace(10,200)*u.K
frequencies = ([350,500]*u.um).to(u.GHz, u.spectral())
fluxes = [dust_emissivity.blackbody.modified_blackbody(frequencies, temperature=T,
                                                       beta=1.75,
                                                       column=1e19*u.cm**-2)
          for T in temperatures]
fluxes = u.Quantity(fluxes)
ratios = fluxes[:,0] / fluxes[:,1]

temperature_map = np.interp(ratio_map_sharc_to_scuba, ratios, temperatures)
pl.figure(2).clf()
pl.imshow(temperature_map, cmap='viridis')
pl.colorbar()
pl.contour(imagecube[0,:,:], levels=[5e4,1e5,2e5], colors=['k']*3)
pl.savefig(paths.fpath("SHARC_SCUBA_temperature_map_experiment.png"))

# not used (and if it is, it needs to be syned with combine_sharc)
# sharc_beam = radio_beam.Beam(11.5*u.arcsec)

scuba_hdr['BUNIT'] = ('cm^-2', 'N(H2), column density of H2')
#sharc_hdr['BUNIT'] = ('cm^-2', 'N(H2), column density of H2')

colmap_sharc_20 = dust_emissivity.dust.colofsnu(frequencies[0],
                                                imagecube[0,:,:],
                                                temperature=20*u.K,
                                                beta=1.750,
                                               )
colmap_sharc_50 = dust_emissivity.dust.colofsnu(frequencies[0],
                                                imagecube[0,:,:],
                                                temperature=50*u.K,
                                                beta=1.750,
                                               )
fits.writeto(filename=paths.cpath('column_maps/sharc_col_50K.fits'), data=colmap_sharc_50.value, header=scuba_hdr, overwrite=True)
fits.writeto(filename=paths.cpath('column_maps/sharc_col_20K.fits'), data=colmap_sharc_20.value, header=scuba_hdr, overwrite=True)

scuba_beam = radio_beam.Beam(8*u.arcsec)

colmap_scuba_20 = dust_emissivity.dust.colofsnu(frequencies[1],
                                                imagecube[1,:,:],
                                                beta=1.750,
                                                temperature=20*u.K)
colmap_scuba_50 = dust_emissivity.dust.colofsnu(frequencies[1],
                                                imagecube[1,:,:],
                                                beta=1.750,
                                                temperature=50*u.K)

fits.writeto(filename=paths.cpath('column_maps/scuba_col_50K.fits'), data=colmap_scuba_50.value, header=scuba_hdr, overwrite=True)
fits.writeto(filename=paths.cpath('column_maps/scuba_col_20K.fits'), data=colmap_scuba_20.value, header=scuba_hdr, overwrite=True)

pl.figure(1).clf()
pl.hist(colmap_sharc_50[np.isfinite(colmap_sharc_50)], bins=np.logspace(19,24),alpha=0.5, log=True)
pl.hist(colmap_scuba_50[np.isfinite(colmap_scuba_50)], bins=np.logspace(19,24),alpha=0.5, log=True)
pl.semilogx()
pl.draw()
pl.show()

herschel_tem = fits.open(paths.Fpath('other/gcmosaic_temp_conv25.fits'))[0]
herschel_tem_sharcgrid,_ = reproject.reproject_interp(herschel_tem, scuba_hdr)
herschel_tem_smoothed = convolution.convolve_fft(herschel_tem_sharcgrid,
                                                 convolution.Gaussian2DKernel(35))
herschel_tem_sharcgrid[np.isnan(herschel_tem_sharcgrid)] = herschel_tem_smoothed[np.isnan(herschel_tem_sharcgrid)]
colmap_sharc_herscheltem = dust_emissivity.dust.colofsnu(frequencies[0],
                                                         imagecube[0,:,:],
                                                         temperature=herschel_tem_sharcgrid*u.K,
                                                         beta=1.750,)

colmap_scuba_herscheltem = dust_emissivity.dust.colofsnu(frequencies[1],
                                                         imagecube[1,:,:],
                                                         temperature=herschel_tem_sharcgrid*u.K,
                                                         beta=1.750,)

fits.writeto(filename=paths.cpath('column_maps/sharc_col_herscheltem.fits'), data=colmap_sharc_herscheltem.value, header=scuba_hdr, overwrite=True)
fits.writeto(filename=paths.cpath('column_maps/scuba_col_herscheltem.fits'), data=colmap_scuba_herscheltem.value, header=scuba_hdr, overwrite=True)

brick_scuba = fits.open(paths.Fpath('continuumdata/scuba_herschel_merge/scuba_Herschel_Feathered_brick.fits'))[0]
herschel_tem_brickgrid,_ = reproject.reproject_interp(herschel_tem, brick_scuba.header)

colmap_scuba_herscheltem_brick = dust_emissivity.dust.colofsnu(frequencies[1],
                                                               u.Quantity(brick_scuba.data, u.Unit(brick_scuba.header['BUNIT'])),
                                                               temperature=herschel_tem_brickgrid*u.K,
                                                               beta=1.750,)

fits.writeto(filename=paths.cpath('column_maps/brick_scuba_col_herscheltem.fits'),
             data=colmap_scuba_herscheltem_brick.value,
             header=brick_scuba.header, overwrite=True)
