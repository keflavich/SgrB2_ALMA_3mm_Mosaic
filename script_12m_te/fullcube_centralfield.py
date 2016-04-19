import numpy as np

# TODO: check that pointing #'s are the same for each (should be, because I did not split out Sgr B2)
ptg_numbers = [3, 53, 54, 55, 56, 57, 65, 66, 67, 68, 69, 76, 77, 78, 79, 87,
               88, 89, 90, 91, 98, 99, 100, 101, 102, 110, 111, 112, 113]
rootpath = '../'

TC_Center = os.path.join(rootpath,'center/SgrB2_a_03_TC_center.ms')
if not os.path.exists(TC_Center):
    assert split(vis=os.path.join(rootpath,
                                  'calibrated/SgrB2_a_03_TC2.calibrated.ms'),
                 outputvis=TC_Center,
                 field=",".join(map(str,ptg_numbers)),
                 datacolumn='data',
                )
TE_Center = os.path.join(rootpath,'center/SgrB2_a_03_TE_center.ms')
if not os.path.exists(TE_Center):
    assert split(vis=os.path.join(rootpath, 'calibrated/SgrB2_a_03_TE.calibrated.ms'),
                 outputvis=TE_Center,
                 field=",".join(map(str,ptg_numbers)),
                 datacolumn='data',
                )



phasecenter='J2000 17:47:20.048 -28.23.35.72'
cell='0.15arcsec' # cell size for imaging (0.5" resoln)
imsize = [1024,1024] # size of image in pixels.

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean.

weighting = 'briggs'
robust=-0.5
threshold = '50.0mJy'

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
        }
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


for spwnum in '3210':
    spwnum = int(spwnum)

    concatvis = os.path.join(rootpath,
                             'center/SgrB2_TC_TE.spw{0}.merge'.format(spwnum))
    if not os.path.exists(concatvis):
        print "# running cvel on all lines in spw{0}".format(spwnum)
        cvelvises = []

        spw = spws['TC'][spwnum]
        for ss in spw.split(","):
            ss = int(ss)

            cvelvis = os.path.join(rootpath, 'center/SgrB2_a_03_TC_center.spw{0}.cvel'.format(ss))
            cvelvises.append(cvelvis)
            if not os.path.exists(cvelvis):
                print("cvel'ing TC spw {0} to {1}".format(spw, ss, ))
                cvel(vis=TC_Center,
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
            cvelvis = os.path.join(rootpath, 'center/SgrB2_a_03_TE_center.spw{0}.cvel'.format(ss))
            cvelvises.append(cvelvis)
            if not os.path.exists(cvelvis):
                print("cvel'ing TE spw {0} to {1}".format(spw, ss, ))
                cvel(vis=TE_Center,
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
        start = nchans_per_cube*ii
        end = nchans_per_cube*(ii+1)
        if end > nchans_total_thiscube:
            end = nchans_total_thiscube
        output = 'piece_of_full_SgrB2_12m_r-0.5_cube.spw{0}.channels{1}to{2}'.format(spwnum, start, end)

        # Channel-based gridding has major bugs when dealing with CVEL'd data
        # It is therefore necessary to compute the frequency gridding by hand
        startfreq = "{0}GHz".format(frange[spwnum][0]/1e3 + start * fstep[spwnum]/1e6)
        width = "{0}kHz".format(fstep[spwnum])


        # LINE IMAGING (MOSAIC MODE)
        if not os.path.exists(output+".image"):
            print "Imaging {0}".format(output)
            os.system('rm -rf ' + output + '*')
            tclean(vis = concatvis,
                   imagename = output,
                   field = '',
                   spw = '', # there should be only one
                   gridder = 'mosaic',
                   specmode = 'cube',
                   width = width,
                   start = startfreq,
                   nchan = nchans_per_cube,
                   veltype = 'radio',
                   outframe = 'LSRK',
                   deconvolver='clark',
                   interactive = F,
                   niter = 10000,
                   imsize = imsize,
                   cell = cell,
                   weighting = weighting,
                   phasecenter = phasecenter,
                   robust = robust,
                   threshold = threshold,
                   savemodel='none')

          
        myimagebase = output
        # I've given up on primary beam correction, at least for now
        exportfits(imagename=myimagebase+'.image',
                   fitsimage=myimagebase+'.image.fits',
                   overwrite=True,
                   dropdeg=True)

