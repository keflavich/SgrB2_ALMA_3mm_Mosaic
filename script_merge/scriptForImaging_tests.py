"""
Joint imaging of 7m and 12m data
Test parameters to figure out how best to reduce the continuum
"""
import os

# scp /Volumes/passport/alma/sgrb2_b3/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/calibrated.ms cleese:/scratch/aginsbur/sgrb2/2013.1.00269.S/science_goal.uid___A001_X121_X4b6/group.uid___A001_X121_X4b7/member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.calibrated.ms
inputvis = ['../member.uid___A001_X121_X4ba/calibrated/SgrB2_a_03_TC.calibrated.ms',
            '../member.uid___A001_X121_X4bc/calibrated/SgrB2_a_03_7M.calibrated.ms',
            ]
concatvis = 'SgrB2_a_03_merge_7m_12m.ms'
if not os.path.exists(concatvis):
    concat(vis=inputvis, concatvis=concatvis)
    plotms(vis=concatvis,yaxis='wt',xaxis='uvdist',spw='0~2:200',
           coloraxis='spw',plotfile='combine_WT.png')

# parameter exploration
niters = [0, 100, 1000, 10000, 100000]
thresholds = [30, 10, 1]

for niter in niters:
    for threshold in thresholds:

        # spectral window 3 from the concatenated data....
        outfilename = 'SgrB2_a_03_7M_12M_concat.continuum.spw3.niter{0}.thresh{1}mjy'.format(niter, threshold)
        os.system('rm -rf {0}.*'.format(outfilename))
        clean(vis = concatvis,
              imagename = outfilename,
              field = 'SgrB2',
              spw = '3,7,11,15,19,23,27,31',
              mode = 'mfs', outframe = 'lsrk',
              psfmode = 'clark',
              imagermode = 'mosaic',
              gain = 0.05,
              interactive = F,
              imsize = [1296, 1296],
              cell = '0.45 arcsec',
              phasecenter = 'J2000 17h47m19.4 -28d23m29',
              weighting = 'briggs',
              negcomponent=1,
              niter = niter, threshold = '{0}mJy'.format(threshold),
              robust = 0.5, usescratch = True)

        exportfits(imagename=outfilename+'.image', fitsimage=outfilename+'.image.fits', overwrite=True)
