"""
Image the whole of Sgr B2 with coarse spatial resolution merging all data
"""
import numpy as np

#ptg_numbers = [3, 53, 54, 55, 56, 57, 65, 66, 67, 68, 69, 76, 77, 78, 79, 87,
#               88, 89, 90, 91, 98, 99, 100, 101, 102, 110, 111, 112, 113]
rootpath = '../../'

vis_7m = os.path.join(rootpath, 'calibrated/SgrB2_a_03_7M.calibrated.ms')
vis_TC = os.path.join(rootpath, 'calibrated/SgrB2_a_03_TC2.calibrated.ms')
vis_TE = os.path.join(rootpath, 'calibrated/SgrB2_a_03_TE.calibrated.ms')


#phasecenter='J2000 17:47:20.048 -28.23.35.72'
phasecenter='J2000 17:47:19.434 -28.23.29.22'
cell='0.15arcsec' # cell size for imaging (0.5" resoln)
imsize = [3200,3200] # size of image in pixels.

# imaging control
# ----------------

weighting = 'briggs'
robust=0.0
# use a high threshold (100 is pretty high!) to make convergence faster
threshold = '100.0mJy'

spws = {'TE':
        {0: '0,4,8,12',
         1: '1,5,9,13',
         2: '2,6,10,14',
         3: '3,7,11,15',
        },
        'TC':
        {0: '0,4,8',
         1: '1,5,9',
         2: '2,6,10',
         3: '3,7,11',
        },
        '7m':
        {0: '0,4,8,12',
         1: '1,5,9,13',
         2: '2,6,10,14',
         3: '3,7,11,15',
        },
       }


nchans_total = {0: 7680, 1: 7680, 2: 7680, 3: 7680}
frange = {0: [90357.27912, 92220.858],
          1: [88552.69612, 90417.088],
          2: [100440.526,102304.10488],
          3: [102301.456,104164.21788],
         }
fstep = {0:250., # kHz
         1:250., # kHz
         2:250., # kHz
         3:250., # kHz
        }
nchans_total = {ii: int(np.abs(np.diff(frange[ii])/fstep[ii]*1000.)[0])
                for ii in frange}

ncubes_per_window = 20

assert 'spwlist' in locals(), "Specify 'spwlist' as an iterable."

for spwnum in spwlist:
    spwnum = int(spwnum)

    concatvis = os.path.join(rootpath,
                             'SgrB2_TC_TE_7m.spw{0}.merge'.format(spwnum))
    if not os.path.exists(concatvis):
        print "# running cvel on all lines in spw{0}".format(spwnum)
        cvelvises = []

        spw = spws['TC'][spwnum]
        for ss in spw.split(","):
            ss = int(ss)

            cvelvis = os.path.join(rootpath, 'SgrB2_a_03_TC.spw{0}.cvel'.format(ss))
            cvelvises.append(cvelvis)
            if not os.path.exists(cvelvis):
                print("cvel'ing TC spw {0} to {1}".format(spw, ss, ))
                cvel(vis=vis_TC,
                     outputvis=cvelvis,
                     passall=False, field='', spw=str(ss), selectdata=True,
                     timerange='', array='', antenna='', scan='', mode='frequency',
                     nchan=nchans_total[spwnum],
                     start='{0}MHz'.format(frange[spwnum][0]),
                     width='{0}kHz'.format(fstep[spwnum]), interpolation='linear',
                     phasecenter='', restfreq='', outframe='LSRK', veltype='radio',
                     hanning=False,)

        spw = spws['TE'][spwnum]
        for ss in spw.split(","):
            cvelvis = os.path.join(rootpath, 'SgrB2_a_03_TE.spw{0}.cvel'.format(ss))
            cvelvises.append(cvelvis)
            if not os.path.exists(cvelvis):
                print("cvel'ing TE spw {0} to {1}".format(spw, ss, ))
                cvel(vis=vis_TE,
                     outputvis=cvelvis,
                     passall=False, field='', spw=str(ss), selectdata=True,
                     timerange='', array='', antenna='', scan='', mode='frequency',
                     nchan=nchans_total[spwnum],
                     start='{0}MHz'.format(frange[spwnum][0]),
                     width='{0}kHz'.format(fstep[spwnum]), interpolation='linear',
                     phasecenter='', restfreq='', outframe='LSRK', veltype='radio',
                     hanning=False,)

        spw = spws['7m'][spwnum]
        for ss in spw.split(","):
            cvelvis = os.path.join(rootpath, 'SgrB2_a_03_7m.spw{0}.cvel'.format(ss))
            cvelvises.append(cvelvis)
            if not os.path.exists(cvelvis):
                print("cvel'ing 7m spw {0} to {1}".format(spw, ss, ))
                cvel(vis=vis_7m,
                     outputvis=cvelvis,
                     passall=False, field='', spw=str(ss), selectdata=True,
                     timerange='', array='', antenna='', scan='', mode='frequency',
                     nchan=nchans_total[spwnum],
                     start='{0}MHz'.format(frange[spwnum][0]),
                     width='{0}kHz'.format(fstep[spwnum]), interpolation='linear',
                     phasecenter='', restfreq='', outframe='LSRK', veltype='radio',
                     hanning=False,)

        assert concat(vis=cvelvises, concatvis=concatvis,)
    else:
        print "Already cvel'd spw {0} to {1}".format(spwnum, concatvis)

    print "# running clean on all lines in spw{0}".format(spwnum)
    nchans_total_thiscube = nchans_total[spwnum]
    nchans_per_cube = int(nchans_total_thiscube/ncubes_per_window)
    for ii in range(ncubes_per_window):
        # include a 1-pixel buffer
        start = nchans_per_cube*ii -1
        if start <= 0:
            start = 0
        end = nchans_per_cube*(ii+1) +1
        if end > nchans_total_thiscube:
            end = nchans_total_thiscube
        output = 'piece_of_full_SgrB2_TETC7m_r0_cube.spw{0}.channels{1}to{2}'.format(spwnum, start, end)

        # Channel-based gridding has major bugs when dealing with CVEL'd data
        # It is therefore necessary to compute the frequency gridding by hand
        startfreq = "{0}GHz".format(frange[spwnum][0]/1e3 + start * fstep[spwnum]/1e6)
        width = "{0}kHz".format(fstep[spwnum])


        # LINE IMAGING (MOSAIC MODE)
        if (not (os.path.exists(output+".image.fits") or
                 os.path.exists(output+".image.pbcor.fits"))):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            # by Kumar Golap & Steve Myers' suggestion, use the three
            # independent MSes rather than the contat'd MS: concat sometimes
            # puts the data in a stupid order 
            # Also per Kumar's suggestion, each tclean can be run in parallel
            # as long as savemodel='none' because this will open them in
            # read-only mode, and since Lustre is a parallel filesystem, it
            # should be happy with multiple reads of different sections
            tclean(vis = cvelvises,
                   imagename = output,
                   field = '',
                   spw = '', # there should be only one
                   gridder = 'mosaic',
                   specmode = 'cube',
                   width = width,
                   start = startfreq,
                   nchan = nchans_per_cube + 2, # 1 channel at either end for buffer
                   veltype = 'radio',
                   outframe = 'LSRK',
                   deconvolver='clark',
                   interactive = F,
                   niter = 1000000, # force the clean to go to threshold
                   imsize = imsize,
                   cell = cell,
                   weighting = weighting,
                   phasecenter = phasecenter,
                   robust = robust,
                   threshold = threshold,
                   savemodel='none')

            myimagebase = output
            exportfits(myimagebase+'.image', myimagebase+'.image.fits',
                       dropdeg=True, overwrite=True)
            impbcor(imagename=myimagebase+'.image',pbimage=myimagebase+'.pb',
                    outfile=myimagebase+'.image.pbcor', overwrite=True)
            exportfits(myimagebase+'.image.pbcor',
                       myimagebase+'.image.pbcor.fits', dropdeg=True,
                       overwrite=True)

            for suffix in ('psf', 'weight', 'sumwt', 'pb', 'model', 'residual',
                           'mask', 'image'):
                os.system('rm -rf {0}.{1}'.format(myimagebase, suffix))
        print("Completed {0}".format(output))
