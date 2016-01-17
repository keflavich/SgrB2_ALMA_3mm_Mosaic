# You should run this with CASA 4.4 

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

myimsize = [2560, 2560]

# Image continuum of the target source and continuum subtraction
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf '+rootname+'.cont.*')
  clean(vis=thevis,
        imagename=rootname+'.cont',
        #field="3,4,15,14,26,140,151,129,141", # corner and center fields for test image
        field='3~151',
        imagermode='mosaic',
        phasecenter='3',
        spw="",
        #spw="0",
        mode="mfs",
        niter=100,
        threshold="0.15mJy",
        psfmode="clark",
        interactive=True,
        mask = [],
        imsize=myimsize,
        cell="0.2arcsec",
        outframe='LSRK',
        weighting="briggs",
        robust = 0.5, 
        usescratch=True)


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
        interpolation="linear",
        threshold="1mJy",
        interactive=True,
        mask=[],
        imsize=myimsize,
        cell=['0.2arcsec'],
        restfreq=v0,
        outframe='LSRK',
        weighting="briggs",
        robust=0.5,
        usescratch=True)

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
        interpolation="linear",
        threshold="1mJy",
        interactive=True,
        mask=[],
        imsize=myimsize,
        cell=['0.2arcsec'],
        restfreq=v1,
        outframe='LSRK',
        weighting="briggs",
        robust=0.5,
        usescratch=True)

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
        interpolation="linear",
        threshold="1mJy",
        interactive=True,
        mask=[],
        imsize=myimsize,
        cell=['0.2arcsec'],
        restfreq=v2,
        outframe='LSRK',
        weighting="briggs",
        robust=0.5,
        usescratch=True)

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
        interpolation="linear",
        threshold="1mJy",
        interactive=True,
        mask=[],
        imsize=myimsize,
        cell=['0.2arcsec'],
        restfreq=v3,
        outframe='LSRK',
        weighting="briggs",
        robust=0.5,
        usescratch=True)



# export fits
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  myimages.add(rootname+'.cont')
  myimages.add(rootname+'.SPW0')
  #myimages.add(rootname+'.SPW1')
  #myimages.add(rootname+'.SPW2')
  #myimages.add(rootname+'.SPW3')
  
  for myimagebase in myimages:
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.flux', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits') # export the corrected image
    exportfits(imagename=myimagebase+'.flux', fitsimage=myimagebase+'.flux.fits') # export the PB image
