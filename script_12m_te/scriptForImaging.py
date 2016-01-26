# You should run this with CASA 4.4

mysteps = [0,1,2,3,4,5]
thesteps = []
step_title = {0: 'Image continuum',
              1: 'Image line cube SPW 0',
              2: 'Image line cube SPW 1',
              3: 'Image line cube SPW 2',
              4: 'Image line cube SPW 3',
              5: 'Export FITS images'}

try:
  print 'List of steps to be executed ...', mysteps
  thesteps = mysteps
except:
  print 'global variable mysteps not set.'
if (thesteps==[]):
  thesteps = range(0,len(step_title))
  print 'Executing all steps: ', thesteps

rootname = 'SgrB2_a_03_TE'

thevis = 'calibrated.ms'

print thevis

v0 = '91.98709GHz' # CH3CN v=0 5( 0)- 4( 0)
v1 = '89.18853GHz' # HCO+ v=0  	1 - 0
v2 = '101.47789GHz' # H2CS 3( 1, 3)- 2( 1, 2)
v3 = '103.04055GHz' # H2CS 3( 0, 3)- 2( 0, 2)

myimages = set([])

myimsize = [4096, 4096]
phasecenter = 'J2000 17:47:19.438000 -28.23.29.78000'

# Image continuum of the target source and continuum subtraction
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  contvis='SgrB2_TE_contsplit.ms'
  split(vis=thevis, outputvis=contvis,
        field='SgrB2', width=384, datacolumn='data')

  os.system('rm -rf '+rootname+'.cont_robust0.5_tclean.*')
  myimagebase = rootname+'.cont_robust0.5_tclean'
  tclean(vis=contvis,
         imagename=myimagebase,
         #field="3,4,15,14,26,140,151,129,141", # corner and center fields for test image
         field='SgrB2',
         gridder='mosaic',
         phasecenter=phasecenter,
         spw="",
         #spw="0",
         specmode="mfs",
         niter=100000,
         threshold="0.15mJy",
         deconvolver="clark",
         interactive=False,
         imsize=myimsize,
         cell="0.125arcsec",
         outframe='LSRK',
         weighting="briggs",
         robust = 0.5, 
         savemodel='none')
  impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
  exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
  exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
  exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
  exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image


  os.system('rm -rf '+rootname+'.cont_robust0.5.*')
  myimagebase = rootname+'.cont_robust0.5'
  clean(vis=contvis,
        imagename=myimagebase,
        #field="3,4,15,14,26,140,151,129,141", # corner and center fields for test image
        field='3~151',
        imagermode='mosaic',
        phasecenter='3',
        spw="",
        #spw="0",
        mode="mfs",
        niter=10000,
        threshold="0.15mJy",
        psfmode="clark",
        interactive=False,
        mask = [],
        imsize=myimsize,
        cell="0.125arcsec",
        outframe='LSRK',
        weighting="briggs",
        robust = 0.5, 
        usescratch=True)
  impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
  exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
  exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', dropdeg=True, overwrite=True) # export the PB image
  exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image

  os.system('rm -rf '+rootname+'.cont_uniform.*')
  myimagebase = rootname+'.cont_uniform'
  clean(vis=contvis,
        imagename=myimagebase,
        #field="3,4,15,14,26,140,151,129,141", # corner and center fields for test image
        field='3~151',
        imagermode='mosaic',
        phasecenter='3',
        spw="",
        #spw="0",
        mode="mfs",
        niter=10000,
        threshold="0.15mJy",
        psfmode="clark",
        interactive=False,
        mask = [],
        imsize=myimsize,
        cell="0.125arcsec",
        outframe='LSRK',
        weighting="uniform",
        usescratch=True)

  impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
  exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
  exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', dropdeg=True, overwrite=True) # export the PB image
  exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image

# Image target SPW 0
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf '+rootname+'.SPW0*')
  clean(vis=thevis,
        imagename=''+rootname+'.SPW0',
        field="3~151",
        imagermode='mosaic',
        phasecenter='3',
        spw="0,4,8,12",
        mode="velocity",
        start='-40km/s',
        nchan=50,
        niter=1000,
        interpolation="linear",
        threshold="1mJy",
        interactive=False,
        mask=[],
        imsize=myimsize,
        cell=['0.125arcsec'],
        restfreq=v0,
        outframe='LSRK',
        weighting="briggs",
        robust=0.5,
        usescratch=True)

  myimages.add(rootname+'.SPW0')

# Image target SPW 1
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf '+rootname+'.SPW1*')
  clean(vis=thevis,
        imagename=rootname+'.SPW1',
        field="3~151",
        imagermode='mosaic',
        phasecenter='3',
        spw="1,5,9,13",
        mode="velocity",
        start='-40km/s',
        nchan=50,
        niter=1000,
        interpolation="linear",
        threshold="1mJy",
        interactive=False,
        mask=[],
        imsize=myimsize,
        cell=['0.125arcsec'],
        restfreq=v1,
        outframe='LSRK',
        weighting="briggs",
        robust=0.5,
        usescratch=True)
  myimages.add(rootname+'.SPW1')

# Image target SPW 2
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf '+rootname+'.SPW2*')
  clean(vis=thevis,
        imagename=rootname+'.SPW2',
        field="3~151",
        imagermode='mosaic',
        phasecenter='3',
        spw="2,6,10,14",
        mode="velocity",
        start='-40km/s',
        nchan=50,
        niter=1000,
        interpolation="linear",
        threshold="1mJy",
        interactive=False,
        mask=[],
        imsize=myimsize,
        cell=['0.125arcsec'],
        restfreq=v2,
        outframe='LSRK',
        weighting="briggs",
        robust=0.5,
        usescratch=True)
  myimages.add(rootname+'.SPW2')

# Image target SPW 3
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf '+rootname+'.SPW3*')
  clean(vis=thevis,
        imagename=rootname+'.SPW3',
        field="3~151",
        imagermode='mosaic',
        phasecenter='3',
        spw='3,7,11,15',
        mode="velocity",
        start='-40km/s',
        nchan=50,
        niter=1000,
        interpolation="linear",
        threshold="1mJy",
        interactive=False,
        mask=[],
        imsize=myimsize,
        cell=['0.125arcsec'],
        restfreq=v3,
        outframe='LSRK',
        weighting="briggs",
        robust=0.5,
        usescratch=True)

  myimages.add(rootname+'.SPW3')


# export fits
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  for myimagebase in myimages:
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
