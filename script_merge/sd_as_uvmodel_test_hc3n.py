
# cvel(vis='HC3N_SgrB2_a_03_cvel_merge.ms',
#      outputvis='HC3N_65kms_merge_cvel.ms',
#      mode='velocity',
#      nchan=1,
#      width='2km/s',
#      start='65km/s',
#      restfreq='90979020000.0Hz',
#     )
# Results in error:
# 2016-06-01 07:56:41     SEVERE  ms::cvel        Exception Reported: Table DataManager error: Invalid operation: TSM: no array in row 0 of column SIGMA_SPECTRUM in /scratch/aginsbur/SgrB2/ALMA/2013.1.00269.S/merge/HC3N_65kms_merge_cvel.ms/table.f20
# *** Error ***  Table DataManager error: Invalid operation: TSM: no array in row 0 of column SIGMA_SPECTRUM in /scratch/aginsbur/SgrB2/ALMA/2013.1.00269.S/merge/HC3N_65kms_merge_cvel.ms/table.f20
# 2016-06-01 07:56:41     SEVERE  cvel::::        An error occurred running task cvel.


# extract a single channel of HC3N at 65 km/s from the UV data
os.system('rm -rf HC3N_65kms_merge_cvel.ms')
split(vis='HC3N_SgrB2_a_03_cvel_merge.ms',
      outputvis='HC3N_65kms_merge_cvel.ms',
      spw='0:40',
      datacolumn='data',
     )

# might have to run this section from a different python session, or just make
# sure CASA has spectral-cube and astropy installed
from spectral_cube import SpectralCube
from astropy import units as u
cube = SpectralCube.read('../HC3N/HC3N_tp_freq_ds_interp.fits')
restfreq = 90979020000.0*u.Hz
chan = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
                               rest_value=restfreq).closest_spectral_channel(65*u.km/u.s)
im_k = cube[chan,:,:]
im_jy = im_k.to(u.Jy, cube.beam.jtok_equiv(cube.spectral_axis[chan]))
im_jy._unit = u.Jy/u.beam
hdu = im_jy.hdu
# add on extra irrelevant dimensions for CASA
hdu.data = hdu.data[None,None,:,:]
hdu.header['CRPIX3'] = 1
hdu.header['CRVAL3'] = restfreq.to(u.Hz).value
hdu.header['CDELT3'] = (cube.spectral_axis[1]-cube.spectral_axis[0]).to(u.Hz).value
hdu.header['CDELT3'] = (1e9, 'This is a hack')
hdu.header['CTYPE3'] = 'FREQ'
hdu.header['CUNIT3'] = 'Hz'
hdu.header['CRPIX4'] = 1
hdu.header['CRVAL4'] = 1
hdu.header['CDELT4'] = 1
hdu.header['CTYPE4'] = 'STOKES'
hdu.header['CUNIT4'] = ''
hdu.writeto('HC3N_65kms_SingleDish_JyBeam.fits', clobber=True)

importfits('HC3N_65kms_SingleDish_JyBeam.fits',
           'HC3N_65kms_SingleDish_JyBeam.image', overwrite=True)

#os.system('cp -r HC3N_65kms_SingleDish_JyBeam.image HC3N_65kms_SingleDish_JyBeam.model')

os.system('rm -rf HC3N_65kms_with_sdmodel_tclean.*/')
tclean(vis='HC3N_65kms_merge_cvel.ms',
       datacolumn='data',
       imagename='HC3N_65kms_with_sdmodel_tclean',
       imsize=[1024,1024],
       cell='0.5arcsec',
       startmodel='HC3N_65kms_SingleDish_JyBeam.image',
       gridder='mosaic',
       weighting='briggs',
       robust=2,
       niter=10000,
       threshold='100mJy',
       phasecenter='J2000 17:47:19.242 -28.23.33.22',
      )

grid = imregrid(imagename='HC3N_65kms_with_sdmodel_tclean.image', template='get')

os.system('rm -rf HC3N_65kms_sdimage_regrid.image')
imregrid(imagename='HC3N_65kms_SingleDish_JyBeam.image', template=grid,
         output='HC3N_65kms_sdimage_regrid.image', overwrite=True)

os.system('rm -rf HC3N_65kms_with_sdmodel_clean.*/')
clean(vis='HC3N_65kms_merge_cvel.ms',
      imagename='HC3N_65kms_with_sdmodel_clean',
      imsize=[1024,1024],
      cell='0.5arcsec',
      modelimage='HC3N_65kms_sdimage_regrid.image',
      imagermode='mosaic',
      weighting='briggs',
      robust=2,
      niter=10000,
      threshold='100mJy',
      #phasecenter='J2000 17:47:19.242 -28.23.33.22',
      phasecenter='J2000 17h47m19.242s -28d23m33.2199s'
     )
myimagebase = 'HC3N_65kms_with_sdmodel_clean'
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux',
        outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor',
           fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True,
           dropdeg=True)

