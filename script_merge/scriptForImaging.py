"""
Joint imaging of 7m and 12m data
"""

# scp /Volumes/passport/alma/sgrb2_b3/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/calibrated.ms cleese:/scratch/aginsbur/sgrb2/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.calibrated.ms
inputvis = ['../member.uid___A001_X121_X4ba/calibrated/SgrB2_a_03_TC.calibrated.ms',
            '../member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.calibrated.ms',
            ]

for line, restfreq in (
                       ('HNC','90.663574GHz'),
                       ('H41a','92034.43MHz'),
                       ('CH3CN','91971.465MHz'),
                       ('HC3N','90979.02MHz'),
                       ('HCOp','89.18853GHz'),
                       ('HCN','88.631847GHz'),
                       ('H2CS303-202','103040.548MHz'),
                       ('H2CO615-616','101332.993MHz'),
                       ('H2CS322-221','103039.93MHz'),
                       ('H2CS321-220','103051.867MHz'),
                       ('H2CS313212','101477.885MHz'),
                       ('CFp','102587.476MHz'),
                      ):
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
          robust = 0.5,
          phasecenter = 'J2000 17h47m19.4 -28d23m29',
          threshold = '34.1mJy',
          pbcor = F,
          usescratch= T)
    myimagebase = output
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True)
    exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', overwrite=True)
