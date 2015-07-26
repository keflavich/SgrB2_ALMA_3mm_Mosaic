vis = '/Volumes/passport/alma/sgrb2_b3/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.lowres.cal.cont'
#imagename = '/Volumes/passport/alma/sgrb2_b3/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.continuum'
#clean(vis = vis,
#      imagename = imagename,
#      field = '0~52', # SgrB2
#      #spw = '3,7,11,15',
#      mode = 'mfs', width = 1, outframe = 'lsrk',
#      imagermode = 'mosaic',
#      interactive = F,
#      imsize = [216, 216],
#      cell = '2.5arcsec',
#      phasecenter = 'J2000 17h47m19.4 -28d23m29',
#    #  mask = 'box[[83pix,95pix],[135pix,156pix]]',
#      weighting = 'briggs',
#      niter = 10000, threshold = '30mJy',
#      robust = 0.5, usescratch = True)

imagename = '/Volumes/passport/alma/sgrb2_b3/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.uniform.continuum'
clean(vis = vis,
      imagename = imagename,
      field = '0~52', # SgrB2
      #spw = '3,7,11,15',
      mode = 'mfs', width = 1, outframe = 'lsrk',
      imagermode = 'mosaic',
      interactive = F,
      imsize = [512, 512],
      cell = '1.1arcsec',
      phasecenter = 'J2000 17h47m19.4 -28d23m29',
    #  mask = 'box[[83pix,95pix],[135pix,156pix]]',
      weighting = 'uniform',
      niter = 50000, threshold = '10mJy',
      usescratch = True)

start = '-100km/s'
nchan = 250

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


#103.03993
# Spw 3, H2CS 322-221
os.system('rm -rf SgrB2_a_03_7M.H2CS322-221.uniform.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.cal.contsub',
  imagename = 'SgrB2_a_03_7M.H2CS322-221.uniform',
  field = '0~52', # SgrB2
  spw = '3,7,11,15',
  mode = 'velocity', restfreq='103039.93MHz', start = start, width = '1.0km/s', nchan=nchan, outframe = 'lsrk',
  imagermode = 'mosaic',
  interactive = F,
  imsize = imsize,
  cell = cell,
  phasecenter = 'J2000 17h47m19.4 -28d23m29',
#  mask = 'box[[83pix,95pix],[135pix,156pix]]',
  weighting = 'uniform',
  niter = 1500, threshold = '650mJy',
  usescratch = True)
