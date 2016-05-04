"""
Joint imaging of 7m and 12m data
"""
import os

calpath = '../calibrated/'

vistemplate = 'SgrB2_a_03_{0}.calibrated.ms'
vis_7m='SgrB2_a_03_7M.calibrated.ms'
vis_tc='SgrB2_a_03_TC2.calibrated.ms'
vis_te='SgrB2_a_03_TE.calibrated.ms'

velocity_range = -165,130
velocity_res = 2.0
nchans = int((velocity_range[1]-velocity_range[0])/velocity_res)

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
                       ('H15NC','88.86569GHz'),
                      ):

    outms_template = "{line}_SgrB2_a_03_{array}.cvel.ms"
    concatvis = "{line}_SgrB2_a_03_cvel_merge.ms".format(line=line)

    # cvel the data first
    for array in ('7M','TC2','TE'):
        outputvis = outms_template.format(line=line, array=array)
        if not os.path.exists(outputvis):
            # maybe cvel2?
            print("cveling line {0} array {1}".format(line, array))
            cvel2(vis=os.path.join(calpath, vistemplate.format(array)),
                  outputvis=outputvis,
                  datacolumn='data',
                  field='SgrB2',
                  mode='velocity',
                  nchan=nchans,
                  start='{0}km/s'.format(velocity_range[0]),
                  width='{0}km/s'.format(velocity_res),
                  interpolation='linear',
                  phasecenter='',
                  restfreq=restfreq,
                  outframe='LSRK',
                  veltype='radio',
                  hanning=False,)
            

    if not os.path.exists(concatvis):
        concat(vis=[outms_template.format(line=line, array=array) for array in ('7M','TC2','TE')],
               concatvis=concatvis)


    output = 'SgrB2_b3_7M_12M.{0}'.format(line)
    #---------------------------------------------------
    # LINE IMAGING (MOSAIC MODE)
    os.system('rm -rf ' + output + '*/')
    print("Imaging {0}".format(output))
    tclean(vis=concatvis,
           imagename=output,
           field='SgrB2', # SgrB2
           spw='',
           gridder='mosaic',
           specmode='cube',
           start='{0}km/s'.format(velocity_range[0]),
           width='{0}km/s'.format(velocity_res), interpolation='linear',
           nchan=nchans,
           restfreq=restfreq,
           veltype='radio',
           outframe='LSRK',
           interactive=F,
           niter=200,
           imsize=[4096,4096],
           cell='0.125arcsec',
           weighting='briggs',
           robust=0.5,
           phasecenter='J2000 17:47:19.242 -28.23.33.22',
           threshold='15mJy',
           savemodel='none',
          )
    myimagebase = output
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True, dropdeg=True)
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', overwrite=True, dropdeg=True)
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', overwrite=True, dropdeg=True)
