vis = 'SgrB2_a_03_7M.lowres.cal.cont'

mask = ['circle[[17h47m20.150s,-28d22m17.35s], 10.2971arcsec]',
        'circle[[17h47m20.067s,-28d23m04.65s], 8.59101arcsec]',]

spw = '0,4,8,12'
freq = '92GHz'

for briggs, briggsname in zip((-5,0.5,5),
                              ('m5','05','p5')):
    for niter in (0, 1000, 10000):
        for threshold in ('10mJy','30mJy'):
            for negcomponent,negname in zip((-1,1),
                                            ('negm1','negp1')):

                imagename = 'SgrB2_a_03_7M.{0}.{1}.{2}.{3}.{4}.continuum'.format(freq,
                                                                             briggsname,
                                                                             niter,
                                                                             threshold,
                                                                                negname)

                os.system('rm -rf {0}.*'.format(imagename))

                clean(vis = vis,
                      imagename = imagename,
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
                      niter = niter, threshold = threshold,
                      robust = briggs, usescratch = True)

                exportfits(imagename+".image",
                           imagename+".image.fits", dropdeg=True,
                           overwrite=True)
