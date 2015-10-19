split(vis='calibrated.ms', outputvis='calibrated.spw1.cont.ms', datacolumn='data', field='SgrB2', spw='1;5;9;13:4000~5000', width=100, )
listobs(vis='calibrated.spw1.cont.ms', listfile='spw1.listobs.txt', overwrite=True)
clean(vis='calibrated.spw1.cont.ms', imagename='spw1_continuum_channels.dirty',
      imsize=[216,216], field='SgrB2', spw='0~3', interactive=F, cell='2.5arcsec',
      phasecenter = 'J2000 17h47m19.4 -28d23m29', weighting = 'briggs', niter = 0,
      threshold = '300mJy', robust = 0.5, usescratch = True, mode='mfs',
      outframe='lsrk', imagermode='mosaic')
exportfits(imagename='spw1_continuum_channels.dirty.image', fitsimage='spw1_continuum_channels.dirty.fits', dropdeg=True, overwrite=True)

clean(vis='calibrated.spw1.cont.ms', imagename='spw1_continuum_channels.clean',
      imsize=[216,216], field='SgrB2', spw='0~3', interactive=F, cell='2.5arcsec',
      phasecenter = 'J2000 17h47m19.4 -28d23m29', weighting = 'briggs', niter = 1000,
      threshold = '300mJy', robust = 0.5, usescratch = True, mode='mfs',
      outframe='lsrk', imagermode='mosaic')
exportfits(imagename='spw1_continuum_channels.clean.image', fitsimage='spw1_continuum_channels.clean.fits', dropdeg=True, overwrite=True)


# do an experiment with 10 different NITER's
for niter in (100,1000,1e4,1e5,1e6,1e7):
    clean(vis='calibrated.spw1.cont.ms', imagename='spw1_continuum_channels.clean.niter{0}'.format(niter),
          imsize=[216,216], field='SgrB2', spw='0~3', interactive=F, cell='2.5arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29', weighting = 'briggs', niter = niter,
          threshold = '30mJy', robust = 0.5, usescratch = True, mode='mfs',
          outframe='lsrk', imagermode='mosaic')
    exportfits(imagename='spw1_continuum_channels.clean.niter{0}.image'.format(niter),
               fitsimage='spw1_continuum_channels.clean.niter{0}.fits'.format(niter),
               dropdeg=True, overwrite=True)
