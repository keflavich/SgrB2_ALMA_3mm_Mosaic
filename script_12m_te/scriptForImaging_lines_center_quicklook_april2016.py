import os
phasecenter='J2000 17:47:20.048 -28.23.35.72'
cell='0.15arcsec' # cell size for imaging (0.5" resoln)
imsize = [1024,1024] # size of image in pixels.

inputvis = 'SgrB2_a_03_TE_center.ms'
for line, restfreq in (
                       ('Heplus63a','102786.9830702896MHz'),
                       ('SO2_220-313','100.87811GHz'),
                       ('SO2_313-202','104.02942GHz'),
                       # this isn't in band; I don't know how I found it('SO2_835-928','92.66036GHz'),
):
    output = 'SgrB2_b3_12M_TE_center.{0}'.format(line)
    os.system('rm -rf ' + output + '*')
    clean(vis = inputvis,
          imagename = output,
          field = 'SgrB2', # SgrB2
          spw = '',
          imagermode = 'mosaic',
          mode = 'velocity',
          width = '2km/s',
          start = '-200km/s',
          nchan = 200,
          restfreq = restfreq,
          outframe = 'LSRK',
          interactive = F,
          niter = 2000,
          imsize = imsize,
          cell = cell,
          weighting = 'briggs',
          phasecenter = phasecenter,
          robust = -0.5,
          threshold = '34.1mJy',
          pbcor = F,
          usescratch= T)
    myimagebase = output
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
    exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)
