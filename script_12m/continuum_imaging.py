vis = 'SgrB2_a_03_TC.calibrated.ms'

mask = ['circle[[17h47m20.150s,-28d22m17.35s], 10.2971arcsec]',
        'circle[[17h47m20.067s,-28d23m04.65s], 8.59101arcsec]',]
mask = None

for spw,freq in zip(('0,4,8,12', '1,5,9,13', '2,6,10,14', '3,7,11,15'),
                    ('92GHz', '90GHz', '100GHz', '102GHz')):
    imagename = 'SgrB2_a_03_12M.{0}.briggs05.continuum'
    clean(vis = vis,
          imagename = imagename.format(freq)+".dirty",
          field = 'SgrB2', # SgrB2
          spw = spw,
          mode = 'mfs', width = 1, outframe = 'lsrk',
          psfmode = 'clark',
          imagermode = 'mosaic',
          gain = 0.05,
          interactive = F,
          imsize = [1296, 1296],
          cell = '0.45 arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          weighting = 'briggs',
          negcomponent=1,
          niter = 0, threshold = '30mJy',
          robust = 0.5, usescratch = True)

    #maskname = imagename.format(freq)+'.mask'
    #os.system('rm -rf {0}'.format(maskname))
    #immath(imagename = imagename.format(freq)+".dirty.image",
    #       outfile = maskname,
    #       expr = 'iif(IM0 > 0.4, 1.0, 0.0)')
    exportfits(imagename.format(freq)+".dirty.image",
               imagename.format(freq)+".dirty.image.fits", dropdeg=True,
               overwrite=True)
    #exportfits(maskname,
    #           maskname+".fits", dropdeg=True,
    #           overwrite=True)
    

    clean(vis = vis,
          imagename = imagename.format(freq),
          field = 'SgrB2', # SgrB2
          spw = spw,
          mode = 'mfs', width = 1, outframe = 'lsrk',
          psfmode = 'clark',
          imagermode = 'mosaic',
          gain = 0.05,
          interactive = F,
          imsize = [1296, 1296],
          cell = '0.45 arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          mask = mask,
          weighting = 'briggs',
          negcomponent=1,
          niter = 10000, threshold = '30mJy',
          robust = 0.5, usescratch = True)

    exportfits(imagename.format(freq)+".image",
               imagename.format(freq)+".image.fits", dropdeg=True,
               overwrite=True)

    imagename = 'SgrB2_a_03_12M.{0}.uniform.continuum'
    clean(vis = vis,
          imagename = imagename.format(freq)+".dirty",
          field = 'SgrB2', # SgrB2
          spw = spw,
          mode = 'mfs', width = 1, outframe = 'lsrk',
          psfmode = 'clark',
          imagermode = 'mosaic',
          interactive = F,
          imsize = [1296, 1296],
          cell = '0.45 arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          weighting = 'uniform',
          niter = 0, threshold = '10mJy',
          gain = 0.05,
          negcomponent=1,
          usescratch = True)

    #ia.open(imagename.format(freq)+".dirty.image")
    #maskname = imagename.format(freq)+'.mask'
    #ia.calcmask(imagename.format(freq)+'.dirty.image > 0.2', name=maskname)
    #ia.done()
    #makemask(mode='copy', inpimage=imagename.format(freq)+".dirty.image",
    #         inpmask="{imn}:{maskn}".format(imn=imagename.format(freq)+".dirty.image",
    #                                        maskn=maskname),
    #         output=maskname,
    #         overwrite=True
    #        )
    exportfits(imagename.format(freq)+".dirty.image",
               imagename.format(freq)+".dirty.image.fits", dropdeg=True,
               overwrite=True)

    clean(vis = vis,
          imagename = imagename.format(freq),
          field = 'SgrB2', # SgrB2
          spw = spw,
          mode = 'mfs', width = 1, outframe = 'lsrk',
          psfmode = 'clark',
          imagermode = 'mosaic',
          interactive = F,
          imsize = [1296, 1296],
          cell = '0.45 arcsec',
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          mask = mask,
          weighting = 'uniform',
          niter = 50000, threshold = '10mJy',
          gain = 0.05,
          negcomponent=1,
          usescratch = True)


    exportfits(imagename.format(freq)+".image",
               imagename.format(freq)+".image.fits", dropdeg=True,
               overwrite=True)
