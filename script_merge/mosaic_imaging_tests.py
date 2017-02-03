vis_7m='SgrB2_a_03_7M.calibrated.ms'
vis_tc='SgrB2_a_03_TC2.calibrated.ms'
vis_te='SgrB2_a_03_TE.calibrated.ms'

velocity_range = -165,130
velocity_range = 45,47
velocity_res = 2.0
nchans = int((velocity_range[1]-velocity_range[0])/velocity_res)
restfreq = '90979.02MHz'

threshold = {'12m': {-2: '40mJy', 0: '15mJy', 2: '15mJy',},
             '7m': {-2: '200mJy', 0: '200mJy', 2: '200mJy'},
             '7m12m': {-2: '40mJy', 0: '15mJy', 2: '15mJy',},
            }

robust = 0
for arrays, suffix in ((['TC2','TE'], '12m'),
                       (['7M'], '7m'),
                       (['TC2','TE','7M'], '7m12m'),
                      ):

    outms_template = "{line}_SgrB2_a_03_{array}.cvel.ms"
    vis = [outms_template.format(array=arr, line='HC3N')
           for arr in arrays]

    for field in range(150):
        output = 'singlechannelfield{3}_SgrB2_b3_{0}_robust{2}.{1}.fits'.format(suffix, 'HC3N',
                                                          robust, field)
        os.system('rm -rf ' + output + '*/')
        print("Imaging {0}".format(output))
        tclean(vis=vis,
               imagename=output,
               field='{0}'.format(field),
               spw='',
               gridder='mosaic',
               specmode='cube',
               start='{0}km/s'.format(velocity_range[0]),
               width='{0}km/s'.format(velocity_res),
               interpolation='linear',
               nchan=nchans,
               restfreq=restfreq,
               veltype='radio',
               outframe='LSRK',
               interactive=F,
               niter=0,
               imsize=[1024,1024],
               cell='0.5arcsec',
               weighting='briggs',
               robust=robust,
               phasecenter='J2000 17:47:19.242 -28.23.33.22',
               threshold=threshold[suffix],
               savemodel='none',
              )
        myimagebase = output
        impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
        exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', overwrite=True, dropdeg=True)

        for fnsuffix in ('pb', 'weight', 'sumwt', 'psf', 'model', 'mask',
                       'image', 'residual'):
            os.system('rm -rf {0}.{1}'.format(output, fnsuffix))
        

        print("Done with {0}".format(output))

