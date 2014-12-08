"""
December 8: try masking on bright stuff
"""

# Creating cube for spw 0 with low velocity resolution, to check for lines/line free regions
os.system('rm -rf SgrB2_a_03_7M.lowres.spw0.*')
#default(clean)
clean(vis = 'SgrB2_a_03_7M.lowres.cal',
  imagename = 'SgrB2_a_03_7M.lowres.spw0.dirty',
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
  niter = 0, threshold = '300mJy',
  robust = 0.5, usescratch = True)

ia.open('SgrB2_a_03_7M.lowres.spw0.dirty')
ia.calcmask(mask='SgrB2_a_03_7M.lowres.spw0.dirty > 1.0', name='mask_gt1jy')
ia.done()

clean(vis = 'SgrB2_a_03_7M.lowres.cal',
  imagename = 'SgrB2_a_03_7M.lowres.spw0.clean',
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
  niter = 10000, threshold = '300mJy',
  robust = 0.5, usescratch = True)
