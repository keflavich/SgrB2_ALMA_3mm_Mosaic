vis = 'calibrated.ms'
phasecenter = 'J2000 17:47:19.438000 -28.23.29.78000'

imsize = [4096,4096]
cell = '0.125arcsec'

print "# running clean on spw0 continuum"
# imaging continuum in spw0 - subtracting stongest lines:
clean(vis = vis,
  imagename = 'SgrB2_TE.contWOlines.spw0',
  field = 'SgrB2', # SgrB2
  spw = '0:1~700;900~1000;1100~2700;2900~3650;3750~4200;4300~5100;5250~6400;6550~6850;7000~7679,'
        '5:1~700;900~1000;1100~2700;2900~3650;3750~4200;4300~5100;5250~6400;6550~6850;7000~7679,'
        '9:1~700;900~1000;1100~2700;2900~3650;3750~4200;4300~5100;5250~6400;6550~6850;7000~7679,'
        '13:1~700;900~1000;1100~2700;2900~3650;3750~4200;4300~5100;5250~6400;6550~6850;7000~7679',
  mode = 'mfs',
  outframe='LSRK',
  imagermode = 'mosaic',
  imsize = imsize,
  cell = cell,
  phasecenter = phasecenter,
  interactive = F,
  niter=20000,
  threshold='0.6mJy')

myimagebase = 'SgrB2_TE.contWOlines.spw0'

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


print "# running clean on spw3 continuum"
# imaging continuum in spw3 - subtracting stongest lines:

clean(vis = vis,
  imagename = 'SgrB2_TE.contWOlines.spw3',
  field = 'SgrB2',
  spw = '3:5050~5150;6100~6200;7000~7200,'
        '7:5050~5150;6100~6200;7000~7200,'
        '11:5050~5150;6100~6200;7000~7200,'
        '15:5050~5150;6100~6200;7000~7200', # to avoid strong lines
  mode = 'mfs',
  outframe='LSRK',
  imagermode = 'mosaic',
  imsize = imsize,
  cell = cell,
  phasecenter = phasecenter,
  interactive = F,
  niter=20000,
  threshold='0.6mJy',
  weighting = 'briggs',
  robust = 0.5,
  usescratch=True)


myimagebase = 'SgrB2_TE.contWOlines.spw3'
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True, dropdeg=True)
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True, dropdeg=True)
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)


inputvis = vis
for line, restfreq in (
                       ('HC3N','90979.02MHz'),
                       ('HNC','90.663574GHz'),
                       ('H41a','92034.43MHz'),
                       #('SiO','86.8469850GHz'),
                       ('CH3CN','91971.465MHz'),
                       ('HCOp','89.18853GHz'),
                       ('HCN','88.631847GHz'),
                       ('C3H2440-515','91.51601GHz'),
                       ('H2CS303-202','103040.548MHz'),
                       ('H2CO615-616','101332.993MHz'),
                       ('H2CS322-221','103039.93MHz'),
                       ('H2CS321-220','103051.867MHz'),
                       ('H2CS313212','101477.885MHz'),
                       ('CFp','102587.476MHz'),
                      ):
    output = 'SgrB2_b3_12M_TE.{0}'.format(line)
    #---------------------------------------------------
    # LINE IMAGING (MOSAIC MODE)
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
          robust = 0.5,
          threshold = '34.1mJy',
          pbcor = F,
          usescratch= T)
    myimagebase = output
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
    exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)


nchans_total = 7680
ncubes_per_window = 20
nchans_per_cube = nchans_total/ncubes_per_window

for spw in '0123':
    print "# running clean on all lines in spw{0}".format(spw)
    inputvis = vis
    for ii in range(ncubes_per_window):
        start = nchans_per_cube*ii
        end = nchans_per_cube*(ii+1)
        output = 'piece_of_SgrB2_TE.line.spw{0}.channels{1}to{2}'.format(spw, start, end)
        #---------------------------------------------------
        # LINE IMAGING (MOSAIC MODE)
        if not os.path.exists(output+".image"):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            clean(vis = inputvis,
                   imagename = output,
                   field = '3~151', # SgrB2
                   spw = spw,
                   imagermode = 'mosaic',
                   mode = 'channel',
                   width = 1,
                   start = start,
                   nchan = nchans_per_cube,
                   chaniter = True,
                   outframe = 'LSRK',
                   interactive = F,
                   niter = 2000,
                   imsize = imsize,
                   cell = cell,
                   weighting = 'briggs',
                   phasecenter = phasecenter,
                   robust = 0.5,
                   threshold = '34.1mJy',
                   pbcor = F,
                   usescratch= T)
          
        if not os.path.exists(output+".image.pbcor"):
            myimagebase = output
            impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
            exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
            exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)
