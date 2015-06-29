#####################
#Start imaging script
#####################

import re

if re.search('^4.2.2', casadef.casa_version) == None:
   sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.2.2')

vis = 'SgrB2_a_03_TC.calibrated.ms'


print "# running clean on spw0 continuum"
# imaging continuum in spw0 - subtracting stongest lines:
clean(vis = vis,
  imagename = 'calibrated.ms.contWOlines.spw0',
  field = '3~151', # SgrB2
#  spw = '0',
  spw = 
'0:1~700;900~1000;1100~2700;2900~3650;3750~4200;4300~5100;5250~6400;6550~6850;7000~7679', # to avoid strong lines
  mode = 'mfs',
  outframe='LSRK',
  imagermode = 'mosaic',
  imsize = [1296, 1296],
  cell = '0.45arcsec',
  phasecenter = 3,
  interactive = F,
  niter=2000,
  threshold='0.6mJy')

# Restoring Beam   : 4.685 arcsec, 2.529 arcsec, 71.9 deg
# RMS 38.8 mJy

myimagebase = 'calibrated.ms.contWOlines.spw0'

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits')
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits')


print "# running clean on spw3 continuum"
# imaging continuum in spw3 - subtracting stongest lines:

clean(vis = vis,
  imagename = 'calibrated.ms.contWOlines.spw3',
  field = '3~151', # SgrB2
#  spw = '3',
  spw = '3:5050~5150;6100~6200;7000~7200', # to avoid strong lines
  mode = 'mfs',
  outframe='LSRK',
  imagermode = 'mosaic',
  imsize = [1296, 1296],
  cell = '0.45arcsec',
  phasecenter = 3,
  interactive = F,
  niter=2000,
  threshold='0.6mJy',
  weighting = 'briggs',
  robust = 0.5,
  usescratch=True)

# Restoring Beam   : 3.56797 arcsec, 1.82822 arcsec, 73.8483 deg
# RMS 0.0283 Jy/beam

myimagebase = 'calibrated.ms.contWOlines.spw3'
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)



for spw in '0123':
    print "# running clean on all lines in spw{0}".format(spw)
    inputvis = vis
    output = 'calibrated.ms.line.spw{0}'.format(spw)
    #---------------------------------------------------
    # LINE IMAGING (MOSAIC MODE)
    os.system('rm -rf ' + output + '*')
    clean(vis = inputvis,
           imagename = output,
           field = '3~151', # SgrB2
           spw = spw,
           imagermode = 'mosaic',
           mode = 'channel',
           width = 1,
           outframe = 'LSRK',
           interactive = F,
           niter = 200,
           imsize = [1296,1296],
           cell = '0.45arcsec',
           weighting = 'briggs',
           phasecenter = 3,
           robust = 0.5,
           threshold = '34.1mJy',
           pbcor = F,
           usescratch= T)

    # RMS 0.0386 Jy/beam
    # Beam    4.098 arcsec x    1.963 arcsec pa= 74.4 deg


    # RMS 0.0391 Jy/beam
    # Beam    3.57 arcsec x    1.75 arcsec pa= 73.59 deg


    # RMS 0.03569205 Jy/beam
    # Beam    4.17 arcsec x    1.99 arcsec pa= 74.54 deg


    # RMS 0.03676751 Jy/beam
    # Beam    3.67 arcsec x    1.75 arcsec pa= 74.32 deg
      
    myimagebase = 'calibrated.ms.line.spw{0}'.format(spw)
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
    exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)

inputvis = vis
output = 'calibrated.ms.line.hcop'
#---------------------------------------------------
# LINE IMAGING (MOSAIC MODE)
os.system('rm -rf ' + output + '*')
clean(vis = inputvis,
      imagename = output,
      field = '3~151', # SgrB2
      spw = spw,
      imagermode = 'mosaic',
      mode = 'velocity',
      width = '2km/s',
      start = '-100km/s',
      nchan = 200,
      restfreq = '89.18853GHz',
      outframe = 'LSRK',
      interactive = F,
      niter = 200,
      imsize = [1296,1296],
      cell = '0.45arcsec',
      weighting = 'briggs',
      phasecenter = 3,
      robust = 0.5,
      threshold = '34.1mJy',
      pbcor = F,
      usescratch= T)
myimagebase = output
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)

output = 'calibrated.ms.line.hcn'
#---------------------------------------------------
# LINE IMAGING (MOSAIC MODE)
os.system('rm -rf ' + output + '*')
clean(vis = inputvis,
      imagename = output,
      field = '3~151', # SgrB2
      spw = spw,
      imagermode = 'mosaic',
      mode = 'velocity',
      width = '2km/s',
      start = '-100km/s',
      nchan = 200,
      restfreq = '88633.9360MHz',
      outframe = 'LSRK',
      interactive = F,
      niter = 200,
      imsize = [1296,1296],
      cell = '0.45arcsec',
      weighting = 'briggs',
      phasecenter = 3,
      robust = 0.5,
      threshold = '34.1mJy',
      pbcor = F,
      usescratch= T)
myimagebase = output
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)
