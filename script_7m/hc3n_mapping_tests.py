
niters = (1,5,100,1000,2000,2250,2500,2750,3000,3250, 3500, 5000,8000,10000,20000)
for niter in niters:
    os.system('rm -rf SgrB2_a_03_7M.HC3N.clean.hogbom.niter{0}*'.format(niter))
    clean(vis = 'SgrB2_a_03_7M.cal.contsub',
          imagename = 'SgrB2_a_03_7M.HC3N.clean.hogbom.niter{0}'.format(niter),
          field = '0~52', # SgrB2
          spw = '0,4,8,12',
          mode = 'velocity', restfreq='90979.02MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
          psfmode='hogbom',
          ftmachine='mosaic',
          imagermode = 'mosaic',
          interactive = F,
          imsize = [216, 216],
          cell = '2.5arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          #mask = maskname,
          weighting = 'briggs',
          niter = niter, threshold = '300mJy',
          robust = 0.5, usescratch = True,
          minpb = 0.4, # attempt to mitigate 'inf' artifacts
          mask = True,
           )
