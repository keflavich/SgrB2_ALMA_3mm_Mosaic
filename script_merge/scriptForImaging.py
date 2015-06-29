
# scp /Volumes/passport/alma/sgrb2_b3/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/calibrated.ms cleese:/scratch/aginsbur/sgrb2/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.calibrated.ms
inputvis = ['../member.uid___A001_X121_X4ba/calibrated/SgrB2_a_03_TC.calibrated.ms',
            '../member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.calibrated.ms',
            ]

for line, restfreq in (('HCOp','89.18853GHz'),('HCN','88633.9360MHz')):
    output = 'SgrB2_b3_7M_12M.{0}'.format(line)
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
          start = '-100km/s',
          nchan = 200,
          restfreq = restfreq,
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
