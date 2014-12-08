vis = 'SgrB2_a_03_7M.lowres.cal.cont'

mask = ['circle[[17h47m20.150s,-28d22m17.35s], 10.2971arcsec]',
        'circle[[17h47m20.067s,-28d23m04.65s], 8.59101arcsec]',]

for spw,freq in zip(('0,4,8,12', '1,5,9,13', '2,6,10,14', '3,7,11,15'),
                    ('92GHz', '90GHz', '100GHz', '102GHz')):
    imagename = 'SgrB2_a_03_7M.{0}.briggs05.continuum'
    clean(vis = vis,
          imagename = imagename.format(freq)+".dirty",
          field = '0~52', # SgrB2
          spw = spw,
          mode = 'mfs', width = 1, outframe = 'lsrk',
          psfmode = 'clark',
          imagermode = 'mosaic',
          gain = 0.05,
          interactive = F,
          imsize = [216, 216],
          cell = '2.5arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          weighting = 'briggs',
          negcomponent=1,
          niter = 0, threshold = '30mJy',
          robust = 0.5, usescratch = True)

    ia.open(imagename.format(freq)+".dirty.image")
    maskname = imagename.format(freq)+'.mask'
    ia.calcmask(imagename.format(freq)+'.dirty.image > 0.2', name=maskname)
    ia.done()

    clean(vis = vis,
          imagename = imagename.format(freq),
          field = '0~52', # SgrB2
          spw = spw,
          mode = 'mfs', width = 1, outframe = 'lsrk',
          psfmode = 'clark',
          imagermode = 'mosaic',
          gain = 0.05,
          interactive = F,
          imsize = [216, 216],
          cell = '2.5arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          mask = maskname,
          weighting = 'briggs',
          negcomponent=1,
          niter = 10000, threshold = '30mJy',
          robust = 0.5, usescratch = True)

    exportfits(imagename.format(freq)+".image",
               imagename.format(freq)+".image.fits", dropdeg=True,
               overwrite=True)

    imagename = 'SgrB2_a_03_7M.{0}.uniform.continuum'
    clean(vis = vis,
          imagename = imagename.format(freq)+".dirty",
          field = '0~52', # SgrB2
          spw = spw,
          mode = 'mfs', width = 1, outframe = 'lsrk',
          psfmode = 'clark',
          imagermode = 'mosaic',
          interactive = F,
          imsize = [512, 512],
          cell = '1.1arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          weighting = 'uniform',
          niter = 0, threshold = '10mJy',
          gain = 0.05,
          negcomponent=1,
          usescratch = True)

    ia.open(imagename.format(freq)+".dirty.image")
    maskname = imagename.format(freq)+'.mask'
    ia.calcmask(imagename.format(freq)+'.dirty.image > 0.2', name=maskname)
    ia.done()

    clean(vis = vis,
          imagename = imagename.format(freq),
          field = '0~52', # SgrB2
          spw = spw,
          mode = 'mfs', width = 1, outframe = 'lsrk',
          psfmode = 'clark',
          imagermode = 'mosaic',
          interactive = F,
          imsize = [512, 512],
          cell = '1.1arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          mask = maskname,
          weighting = 'uniform',
          niter = 50000, threshold = '10mJy',
          gain = 0.05,
          negcomponent=1,
          usescratch = True)


    exportfits(imagename.format(freq)+".image",
               imagename.format(freq)+".image.fits", dropdeg=True,
               overwrite=True)
