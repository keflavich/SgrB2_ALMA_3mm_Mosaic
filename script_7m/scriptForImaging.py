import re

if re.search('^4.2.1', casadef.casa_version) == None:
 sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.2.1')


print "# Split out science data."

os.system('rm -rf SgrB2_a_03_7M.lowres.cal*')
split(vis = 'calibrated.ms', outputvis = 'SgrB2_a_03_7M.lowres.cal', field = '4~56', keepflags = T, width = 10, datacolumn = 'data')
os.system('rm -rf SgrB2_a_03_7M.cal*')
split(vis = 'calibrated.ms', outputvis = 'SgrB2_a_03_7M.cal', field = '4~56', keepflags = T, datacolumn = 'data')

print "# Running clean."

# You have not specified a source Id, I will assume you want to clean the science target(s).

# Creating cube for spw 0 with low velocity resolution, to check for lines/line free regions
os.system('rm -rf SgrB2_a_03_7M.lowres.spw0.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.lowres.cal',
  imagename = 'SgrB2_a_03_7M.lowres.spw0',
  field = '0~52', # SgrB2
  spw = '0,4,8,12',
  mode = 'channel', width = 1, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 1000, threshold = '300mJy',
  robust = 0.5, usescratch = True)

# Creating cube for spw 1 with low velocity resolution, to check for lines/line free regions
os.system('rm -rf SgrB2_a_03_7M.lowres.spw1.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.lowres.cal',
  imagename = 'SgrB2_a_03_7M.lowres.spw1',
  field = '0~52', # SgrB2
  spw = '1,5,9,13',
  mode = 'channel', width = 1, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 1000, threshold = '300mJy',
  robust = 0.5, usescratch = True)

# Creating cube for spw 2 with low velocity resolution, to check for lines/line free regions
os.system('rm -rf SgrB2_a_03_7M.lowres.spw2.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.lowres.cal',
  imagename = 'SgrB2_a_03_7M.lowres.spw2',
  field = '0~52', # SgrB2
  spw = '2,6,10,14',
  mode = 'channel', width = 1, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 1000, threshold = '300mJy',
  robust = 0.5, usescratch = True)

# Creating cube for spw 2 with low velocity resolution, to check for lines/line free regions
os.system('rm -rf SgrB2_a_03_7M.lowres.spw3.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.lowres.cal',
  imagename = 'SgrB2_a_03_7M.lowres.spw3',
  field = '0~52', # SgrB2
  spw = '3,7,11,15',
  mode = 'channel', width = 1, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 1000, threshold = '300mJy',
  robust = 0.5, usescratch = True)


# Continuum subtraction on low spectral resolution data. Not easy to spot the continuum, but...
os.system('rm -rf SgrB2_a_03_7M.lowres.cal.cont*')
uvcontsub(vis = 'SgrB2_a_03_7M.lowres.cal',
          field = '',
          fitspw = '0:65~80;175~255;375~410;430~460;638~652;730~740,1:0~30;110~130;165~175;185~195;210~235;550~580;650~670;715~750,2:480~500;615~640,3:15~80;415~450;550~600',
          combine = 'spw',
          solint= 'int',
          fitorder = 1,
          want_cont = True)

# Continuum subtraction on full spectral resolution data. Not easy to spot the continuum, but...
os.system('rm -rf SgrB2_a_03_7M.cal.cont*')
uvcontsub(vis = 'SgrB2_a_03_7M.cal',
          field = '',
          fitspw = '0:600~1200;1600~1800;2040~2200;2300~2400;2450~2550;2650~2800;2900~3000;3150~3900;4050~4450;4550~4750;4870~5300;5500~6400;7000~7500;7650~8079,1:300~3000;4100~4600;5850~6300;6800~7050;7850~8050,2:550~850;1800~3800;3950~4450;4500~5050;5200~6700;7000~8000,3:200~1150;1400~4500;5500~7200,4:600~1200;1600~1800;2000~2500;2700~3000;3150~3400;3500~3850;4050~4400;4600~5300;5700~6400;7300~7500;7650~8079,5:250~2800;3900~4010;4100~4600;5500~6300;6500~7300,6:400~850;1400~1500;1600~1900;2100~3800;4000~5050;5300~6000;6150~6800;7000~8079,7:200~1100;1300~1600;1800~2800;3000~4550;4700~5350;5500~6400;6600~7250;7500~8000,8:350~450;550~1200;1600~1800;1950~2600;2750~2800;3150~3900;4000~4400;4550~5350;5500~6650;7050~7500;7650~8079,9:250~3000;4100~4600;5700~6300;6450~7250;7850~8079,10:450~800;1600~3700;3950~4400;4500~5050;5250~6800;6900~8079,11:200~1100;1750~4600;4750~5300;5500~7300;7500~7900,12:750~1200;1900~2550;3100~3350;3550~3900;4000~4410;4550~4750;4850~5350;5500~6350;7000~7450;7650~7850,13:250~3000;4100~4700;5500~6350;6450~7350;7800~8000,14:450~800;1400~3750;4000~5050;5150~5500;5600~6800;7100~7450;7570~8079,15:0~4500;4700~7300;7650~7950',
          combine = 'spw',
          solint= 'int',
          fitorder = 0,
          want_cont = True)

start = '-100km/s'
nchan = 250


# Spw 0, HNC (K=3)
os.system('rm -rf SgrB2_a_03_7M.HNC.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.HNC',
  field = '0~52', # SgrB2
  spw = '0,4,8,12',
  mode = 'velocity', restfreq='90663.564MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 20000, threshold = '300mJy',
  robust = 0.5, usescratch = True)

# Spw 0, HC3N
os.system('rm -rf SgrB2_a_03_7M.HC3N.dirty*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.HC3N.dirty',
  field = '0~52', # SgrB2
  spw = '0,4,8,12',
  mode = 'velocity', restfreq='90979.02MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 0, threshold = '300mJy',
  robust = 0.5, usescratch = True,
  #minpb = 0.8, # attempt to mitigate 'inf' artifacts
     )

maskname = 'SgrB2_a_03_7M.HC3N.mask'
os.system('rm -rf {0}'.format(maskname))
immath(imagename = 'SgrB2_a_03_7M.HC3N.dirty.image',
       outfile = maskname,
       expr = 'iif(IM0 > 1.0, 1.0, 0.0)')
exportfits("SgrB2_a_03_7M.HC3N.dirty.image",
           "SgrB2_a_03_7M.HC3N.dirty.image.fits", dropdeg=True,
           overwrite=True)
exportfits(maskname,
           maskname+".fits", dropdeg=True,
           overwrite=True)

# Problems exist starting around 2000-3000
for niter in (2000,20000):
    os.system('rm -rf SgrB2_a_03_7M.HC3N.clean.niter{0}*'.format(niter))
    clean(vis = 'SgrB2_a_03_7M.cal.contsub',
      imagename = 'SgrB2_a_03_7M.HC3N.clean.niter{0}'.format(niter),
      field = '0~52', # SgrB2
      spw = '0,4,8,12',
      mode = 'velocity', restfreq='90979.02MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
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



# Spw 0, H41a
os.system('rm -rf SgrB2_a_03_7M.H41a.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.H41a',
  field = '0~52', # SgrB2
  spw = '0,4,8,12',
  mode = 'velocity', restfreq='92034.43MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 20000, threshold = '300mJy',
  robust = 0.5, usescratch = True)


# Spw 0, CH3CN (K=3)
os.system('rm -rf SgrB2_a_03_7M.CH3CN_5-4_3.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.CH3CN_5-4_3',
  field = '0~52', # SgrB2
  spw = '0,4,8,12',
  mode = 'velocity', restfreq='91971.465MHz', start = '-250km/s', width = '1.0km/s', nchan=500, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 2000, threshold = '450mJy',
  robust = 0.5, usescratch = True)

# Spw 1, HCO+
os.system('rm -rf SgrB2_a_03_7M.HCOp.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.HCOp',
  field = '0~52', # SgrB2
  spw = '1,5,9,13',
  mode = 'velocity', restfreq='89188.526MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 2000, threshold = '450mJy',
  robust = 0.5, usescratch = True)

# Spw 1, HCN
os.system('rm -rf SgrB2_a_03_7M.HCN.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.HCN',
  field = '0~52', # SgrB2
  spw = '1,5,9,13',
  mode = 'velocity', restfreq='88.631847GHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.5arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 2000, threshold = '450mJy',
  robust = 0.5, usescratch = True)

# Spw 2, H2CS # WAS PREVIOUSLY .H2CS. - hasn't been run since renaming
os.system('rm -rf SgrB2_a_03_7M.H2CS313212.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.H2CS313212',
  field = '0~52', # SgrB2
  spw = '2,6,10,14',
  mode = 'velocity', restfreq='101477.885MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.2arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 2000, threshold = '450mJy',
  robust = 0.5, usescratch = True)

# Spw 2, H2CO
os.system('rm -rf SgrB2_a_03_7M.H2CO615-616.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.H2CO615-616',
  field = '0~52', # SgrB2
  spw = '2,6,10,14',
  mode = 'velocity', restfreq='101332.993MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.2arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 2000, threshold = '450mJy',
  robust = 0.5, usescratch = True)

# Spw 3, H2CS 303-202
os.system('rm -rf SgrB2_a_03_7M.H2CS303-202.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.H2CS303-202',
  field = '0~52', # SgrB2
  spw = '3,7,11,15',
  mode = 'velocity', restfreq='103040.548MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.2arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 1500, threshold = '650mJy',
  robust = 0.5, usescratch = True)

#103.03993
# Spw 3, H2CS 322-221
os.system('rm -rf SgrB2_a_03_7M.H2CS322-221.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.H2CS322-221',
  field = '0~52', # SgrB2
  spw = '3,7,11,15',
  mode = 'velocity', restfreq='103039.93MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.2arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'uniform',
  niter = 1500, threshold = '650mJy',
  usescratch = True)


# Spw 3, H2CS 321-220
os.system('rm -rf SgrB2_a_03_7M.H2CS321-220.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.H2CS321-220',
  field = '0~52', # SgrB2
  spw = '3,7,11,15',
  mode = 'velocity', restfreq='103051.867MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.2arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'briggs',
  niter = 1500, threshold = '650mJy',
  robust = 0.5, usescratch = True)

# Spw 3, CF+ (continuum subtracted data seem dominated by artifacts, try on non
# continuum subtracted data, but CF+ seems largely undetected, except possibly
# on SgrB2(N))
os.system('rm -rf SgrB2_a_03_7M.CFp.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal',
  imagename = 'SgrB2_a_03_7M.CFp',
  field = '0~52', # SgrB2
  spw = '3,7,11,15',
  mode = 'velocity', restfreq='102587.476MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = [216, 216],
  cell = '2.2arcsec',
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
  mask = ['box[[83pix,95pix],[135pix,156pix]]','box[[62pix,143pix],[71pix,149pix]]','box[[142pix,69pix],[149pix,78pix]]'],
  weighting = 'briggs',
  niter = 1000, threshold = '250mJy',
  robust = 0.5, usescratch = True)

myimages = ['SgrB2_a_03_7M.HNC', 'SgrB2_a_03_7M.CH3CN_5-4_3',
            'SgrB2_a_03_7M.HCOp', 'SgrB2_a_03_7M.HCN', 'SgrB2_a_03_7M.H2CS313212',
            'SgrB2_a_03_7M.H2CO615-616', 'SgrB2_a_03_7M.H2CS303-202',
            'SgrB2_a_03_7M.H2CS321-220', 'SgrB2_a_03_7M.CFp',
            'SgrB2_a_03_7M.HC3N.dirty',
            'SgrB2_a_03_7M.HC3N.clean',
            'SgrB2_a_03_7M.H41a',
            'SgrB2_a_03_7M.lowres.spw0',
            'SgrB2_a_03_7M.lowres.spw1',
            'SgrB2_a_03_7M.lowres.spw2',
            'SgrB2_a_03_7M.lowres.spw3',
           ]
for myimagebase in myimages:
  impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
  exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
  exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', dropdeg=True, overwrite=True)



