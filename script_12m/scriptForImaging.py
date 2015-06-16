#####################
#Start imaging script
#####################

import re

if re.search('^4.2.2', casadef.casa_version) == None:
   sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.2.2')

# time-averaging to speed up clean
split(vis='SgrB2_a_03_TC.calibrated.ms', outputvis='calibrated.ms.av60s.split',
      datacolumn='data', timebin='60s')


print "# running clean on spw0 continuum"
# imaging continuum in spw0 - subtracting stongest lines:
clean(vis = 'calibrated.ms.av60s.split',
  imagename = 'calibrated.ms.av60s.split.contWOlines.spw0',
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
  interactive = T,
  niter=200,
  threshold='0.6mJy')

# Restoring Beam   : 4.685 arcsec, 2.529 arcsec, 71.9 deg
# RMS 38.8 mJy

myimagebase = 'calibrated.ms.av60s.split.contWOlines.spw0'

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits')
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits')


print "# running clean on spw3 continuum"
# imaging continuum in spw3 - subtracting stongest lines:

clean(vis = 'calibrated.ms.av60s.split',
  imagename = 'calibrated.ms.av60s.split.contWOlines.spw3',
  field = '3~151', # SgrB2
#  spw = '3',
  spw = '3:5050~5150;6100~6200;7000~7200', # to avoid strong lines
  mode = 'mfs',
  outframe='LSRK',
  imagermode = 'mosaic',
  imsize = [1296, 1296],
  cell = '0.45arcsec',
  phasecenter = 3,
  interactive = T,
  niter=1000,
  threshold='0.6mJy',
  weighting = 'briggs',
  robust = 0.5,
  usescratch=True)

# Restoring Beam   : 3.56797 arcsec, 1.82822 arcsec, 73.8483 deg
# RMS 0.0283 Jy/beam

myimagebase = 'calibrated.ms.av60s.split.contWOlines.spw3'
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits')
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits')



print "# running clean on a line in spw0"
input = 'calibrated.ms.av60s.split'
output = 'calibrated.ms.av60s.split.line.spw0'
#---------------------------------------------------
# LINE IMAGING (MOSAIC MODE)
os.system('rm -rf ' + output + '*')
clean(vis = input,
       imagename = output,
       field = '3~151', # SgrB2
       spw = '0',
       imagermode = 'mosaic',
       mode = 'velocity',
       nchan = 70,
       start = '-70km/s',
       width = '2.0km/s',
       restfreq = '90.96GHz',
       outframe = 'LSRK',
       interactive = T,
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

print "#creating Fits files"
imagename = output
impbcor(imagename=imagename+'.image', pbimage=imagename+'.flux', outfile=imagename+'.image.pbcor', overwrite=True)
exportfits(imagename=imagename+'.image.pbcor', fitsimage=imagename+'.image.pbcor.fits', overwrite=True)
exportfits(imagename=imagename+'.flux',        fitsimage=imagename+'.flux.fits',        overwrite=True)


print "# running clean on a line in spw3"
input = 'calibrated.ms.av60s.split'
output = 'calibrated.ms.av60s.split.line.spw3'
#---------------------------------------------------
# LINE IMAGING (MOSAIC MODE)
os.system('rm -rf ' + output + '*')
clean(vis = input,
       imagename = output,
       field = '3~151', # SgrB2
       spw = '3',
       imagermode = 'mosaic',
       mode = 'velocity',
       nchan = 100,
       start = '-100km/s',
       width = '2.0km/s',
       restfreq = '104.015GHz',
       outframe = 'LSRK',
       interactive = T,
       niter = 200,
       imsize = [1296,1296],
       cell = '0.45arcsec',
       weighting = 'briggs',
       phasecenter = 3,
       robust = 0.5,
       threshold = '34.1mJy',
       pbcor = F,
       usescratch= T)

# RMS 0.0391 Jy/beam
# Beam    3.57 arcsec x    1.75 arcsec pa= 73.59 deg

print "#creating Fits files"
imagename = output
impbcor(imagename=imagename+'.image', pbimage=imagename+'.flux', outfile=imagename+'.image.pbcor', overwrite=True)
exportfits(imagename=imagename+'.image.pbcor', fitsimage=imagename+'.image.pbcor.fits', overwrite=True)
exportfits(imagename=imagename+'.flux',        fitsimage=imagename+'.flux.fits',        overwrite=True)


# Spectral line imaging
print "# running clean on a line in spw1"
clean(vis = 'calibrated.ms.av60s.split',
       imagename = 'calibrated.ms.av60s.split.line.spw1',
       field = '3~151', # SgrB2
       spw = '1',
       imagermode = 'mosaic',
       mode = 'velocity',
       nchan = 200,
       start = '-200km/s',
       width = '2km/s',
       restfreq = '89.200013GHz',
       outframe = 'LSRK',
       interactive = T,
       niter = 200,
       imsize = [1296,1296],
       cell = '0.45arcsec',
       weighting = 'briggs',
       phasecenter = 3,
       robust = 0.5,
       threshold = '35mJy',
       pbcor = F,
       usescratch= T)

# RMS 0.03569205 Jy/beam
# Beam    4.17 arcsec x    1.99 arcsec pa= 74.54 deg

myimagebase = 'calibrated.ms.av60s.split.line.spw1'
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits')
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits')

# Spectral line imaging
clean(vis = 'calibrated.ms.av60s.split',
       imagename = 'calibrated.ms.av60s.split.line.spw2',
       field = '3~151', # SgrB2
       spw = '2',
       imagermode = 'mosaic',
       mode = 'velocity',
       nchan = 100,
       start = '-100km/s',
       width = '2km/s',
       restfreq = '101.37GHz',
       outframe = 'LSRK',
       interactive = T,
       niter = 200,
       imsize = [1296,1296],
       cell = '0.45arcsec',
       weighting = 'briggs',
       phasecenter = 3,
       robust = 0.5,
       threshold = '35mJy',
       pbcor = F,
       usescratch= T)

# RMS 0.03676751 Jy/beam
# Beam    3.67 arcsec x    1.75 arcsec pa= 74.32 deg
  
myimagebase = 'calibrated.ms.av60s.split.line.spw2'
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits')
exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits')
