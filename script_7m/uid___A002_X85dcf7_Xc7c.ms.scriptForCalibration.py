# ALMA Data Reduction Script

# Calibration

thesteps = []
step_title = {0: 'Import of the ASDM',
              1: 'Fix of SYSCAL table times',
              2: 'listobs',
              3: 'A priori flagging',
              4: 'Generation and time averaging of the WVR cal table',
              5: 'Generation of the Tsys cal table',
              6: 'Generation of the antenna position cal table',
              7: 'Application of the WVR, Tsys and antpos cal tables',
              8: 'Split out science SPWs and time average',
              9: 'Listobs, clear pointing table, and save original flags',
              10: 'Initial flagging',
              11: 'Putting a model for the flux calibrator(s)',
              12: 'Save flags before bandpass cal',
              13: 'Bandpass calibration',
              14: 'Save flags before gain cal',
              15: 'Gain calibration',
              16: 'Save flags before applycal',
              17: 'Application of the bandpass and gain cal tables',
              18: 'Split out corrected column'}

if 'applyonly' not in globals(): applyonly = False
try:
  print 'List of steps to be executed ...', mysteps
  thesteps = mysteps
except:
  print 'global variable mysteps not set.'
if (thesteps==[]):
  thesteps = range(0,len(step_title))
  print 'Executing all steps: ', thesteps

# The Python variable 'mysteps' will control which steps
# are executed when you start the script using
#   execfile('scriptForCalibration.py')
# e.g. setting
#   mysteps = [2,3,4]# before starting the script will make the script execute
# only steps 2, 3, and 4
# Setting mysteps = [] will make it execute all steps.

import re

import os

if applyonly != True: es = aU.stuffForScienceDataReduction() 


#if re.search('^4.2.1', casadef.casa_version) == None:
# sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.2.1')


# CALIBRATE_AMPLI: J1733-130
# CALIBRATE_ATMOSPHERE: J1700-2610,J1733-130,SgrB2
# CALIBRATE_BANDPASS: J1700-2610
# CALIBRATE_FLUX: J1733-130
# CALIBRATE_FOCUS: 
# CALIBRATE_PHASE: J1744-3116
# CALIBRATE_POINTING: J1700-2610
# OBSERVE_TARGET: SgrB2

# Using reference antenna = CM03

# Import of the ASDM
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  if os.path.exists('uid___A002_X85dcf7_Xc7c.ms') == False:
    importasdm('uid___A002_X85dcf7_Xc7c', asis='Antenna Station Receiver Source CalAtmosphere CalWVR')
  if applyonly != True: es.fixForCSV2555('uid___A002_X85dcf7_Xc7c.ms')

# Fix of SYSCAL table times
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  from recipes.almahelpers import fixsyscaltimes
  fixsyscaltimes(vis = 'uid___A002_X85dcf7_Xc7c.ms')

print "# A priori calibration"

# listobs
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.listobs')
  listobs(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    listfile = 'uid___A002_X85dcf7_Xc7c.ms.listobs')
  
  

# A priori flagging
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  flagdata(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    mode = 'manual',
    spw = '0~23',
    autocorr = T,
    flagbackup = F)
  
  flagdata(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    mode = 'manual',
    intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*',
    flagbackup = F)
  
  flagcmd(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'plot',
    plotfile = 'uid___A002_X85dcf7_Xc7c.ms.flagcmd.png')
  
  flagcmd(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'apply')
  

# Generation and time averaging of the WVR cal table
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  

# Generation of the Tsys cal table
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.tsys') 
  gencal(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.tsys',
    caltype = 'tsys')
  
  if applyonly != True: aU.plotbandpass(caltable='uid___A002_X85dcf7_Xc7c.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='5~123',showfdm=True, 
    field='', figfile='uid___A002_X85dcf7_Xc7c.ms.tsys.plots.overlayTime/uid___A002_X85dcf7_Xc7c.ms.tsys') 
  
  
  if applyonly != True: es.checkCalTable('uid___A002_X85dcf7_Xc7c.ms.tsys', msName='uid___A002_X85dcf7_Xc7c.ms', interactive=False) 
  

# Generation of the antenna position cal table
mystep = 6
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Position for antenna CM01 is derived from baseline run made on 2014-06-02 04:15:36.
  
  # Position for antenna CM12 is derived from baseline run made on 2014-06-02 04:15:36.
  
  # Position for antenna CM11 is derived from baseline run made on 2014-06-02 04:15:36.
  
  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.antpos') 
  gencal(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.antpos',
    caltype = 'antpos',
    antenna = 'CM01,CM11,CM12',
    parameter = [0,0,0,0,0,0,0,0,0])
  #  parameter = [0.000254319605787,-0.000394619449271,-6.89422103618e-05,0.000187085476452,-0.00047325931429,-0.00032447903775,0.000360958191864,-0.000491564156966,-1.53252508184e-05])
  

# Application of the WVR, Tsys and antpos cal tables
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  
  from recipes.almahelpers import tsysspwmap
  tsysmap = tsysspwmap(vis = 'uid___A002_X85dcf7_Xc7c.ms', tsystable = 'uid___A002_X85dcf7_Xc7c.ms.tsys', tsysChanTol = 1)
  
  
  
  applycal(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    field = '0',
    spw = '16,18,20,22',
    gaintable = ['uid___A002_X85dcf7_Xc7c.ms.tsys', 'uid___A002_X85dcf7_Xc7c.ms.antpos'],
    gainfield = ['0', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    field = '1',
    spw = '16,18,20,22',
    gaintable = ['uid___A002_X85dcf7_Xc7c.ms.tsys', 'uid___A002_X85dcf7_Xc7c.ms.antpos'],
    gainfield = ['1', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1744-3116 didn't have any Tsys measurement, so I used the one made on SgrB2. This is probably Ok.
  
  applycal(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    field = '2',
    spw = '16,18,20,22',
    gaintable = ['uid___A002_X85dcf7_Xc7c.ms.tsys', 'uid___A002_X85dcf7_Xc7c.ms.antpos'],
    gainfield = ['3', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    field = '3~55',
    spw = '16,18,20,22',
    gaintable = ['uid___A002_X85dcf7_Xc7c.ms.tsys', 'uid___A002_X85dcf7_Xc7c.ms.antpos'],
    gainfield = ['3', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[]],
    calwt = T,
    flagbackup = F)
  
  if applyonly != True: es.getCalWeightStats('uid___A002_X85dcf7_Xc7c.ms') 
  

# Split out science SPWs and time average
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split') 
  split(vis = 'uid___A002_X85dcf7_Xc7c.ms',
    outputvis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    datacolumn = 'corrected',
    spw = '16,18,20,22',
    keepflags = T)
  
  

print "# Calibration"

# Listobs, clear pointing table, and save original flags
mystep = 9
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.listobs')
  listobs(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    listfile = 'uid___A002_X85dcf7_Xc7c.ms.split.listobs')
  
  tb.open('uid___A002_X85dcf7_Xc7c.ms.split/POINTING', nomodify = False)
  a = tb.rownumbers()
  tb.removerows(a)
  tb.close()
  
  if not os.path.exists('uid___A002_X85dcf7_Xc7c.ms.split.flagversions/Original.flags'):
    flagmanager(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
      mode = 'save',
      versionname = 'Original')
  
  

# Initial flagging
mystep = 10
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Flagging shadowed data
  
  flagdata(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    mode = 'shadow',
    flagbackup = F)
  
  # Note TST: Flagging edge channels (needs more tweaking...)
  
  flagdata(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    mode = 'manual',
    spw = '0:0~79;8080~8159,1:0~79;8080~8159,2:0~79;8080~8159,3:0~79;8080~8159',
    flagbackup = F)
  
  

# Putting a model for the flux calibrator(s)
mystep = 11
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

### Note TST: lines below use flux values as interpolated from measurements in the calibrator database

  setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    field = '1', # source name = J1733-130
    spw = '0', # center frequency of spw = 91.2832960712GHz
    standard = 'manual',
    fluxdensity = [2.21283, 0, 0, 0]) # frequency of measurement = 91.2812420154GHz
  
  setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    field = '1', # source name = J1733-130
    spw = '1', # center frequency of spw = 89.4795398714GHz
    standard = 'manual',
    fluxdensity = [2.24225, 0, 0, 0]) # frequency of measurement = 91.2812420154GHz
  
  setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    field = '1', # source name = J1733-130
    spw = '2', # center frequency of spw = 101.366296071GHz
    standard = 'manual',
    fluxdensity = [2.06462, 0, 0, 0]) # frequency of measurement = 91.2812420154GHz
  
  setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    field = '1', # source name = J1733-130
    spw = '3', # center frequency of spw = 103.226396556GHz
    standard = 'manual',
    fluxdensity = [2.039925, 0, 0, 0]) # frequency of measurement = 91.2812420154GHz
  
### Note TST: lines below use flux values as stored in the asdm; these may not be accurate, better use the ones above
###           interpolated from measurements found in the calibrator database

#  setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
#    field = '1', # source name = J1733-130
#    spw = '0', # center frequency of spw = 91.2832960712GHz
#    standard = 'manual',
#    fluxdensity = [2.37679537166, 0, 0, 0]) # frequency of measurement = 91.2812420154GHz
  
#  setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
#    field = '1', # source name = J1733-130
#    spw = '1', # center frequency of spw = 89.4795398714GHz
#    standard = 'manual',
#    fluxdensity = [2.37679537166, 0, 0, 0]) # frequency of measurement = 91.2812420154GHz
  
#  setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
#    field = '1', # source name = J1733-130
#    spw = '2', # center frequency of spw = 101.366296071GHz
#    standard = 'manual',
#    fluxdensity = [2.37679537166, 0, 0, 0]) # frequency of measurement = 91.2812420154GHz
  
#  setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
#    field = '1', # source name = J1733-130
#    spw = '3', # center frequency of spw = 103.226396556GHz
#    standard = 'manual',
#    fluxdensity = [2.37679537166, 0, 0, 0]) # frequency of measurement = 91.2812420154GHz
  
  

# Save flags before bandpass cal
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    mode = 'save',
    versionname = 'BeforeBandpassCalibration')
  
  

# Bandpass calibration
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.ap_pre_bandpass') 
  
  gaincal(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.ap_pre_bandpass',
    field = '0', # J1700-2610
    spw = '0:3264~4896,1:3264~4896,2:3264~4896,3:3264~4896',
    solint = 'int',
    refant = 'CM03',
    calmode = 'p')
  
  if applyonly != True: es.checkCalTable('uid___A002_X85dcf7_Xc7c.ms.split.ap_pre_bandpass', msName='uid___A002_X85dcf7_Xc7c.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.bandpass') 
  bandpass(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.bandpass',
    field = '0', # J1700-2610
    solint = 'inf',
    combine = 'scan',
    refant = 'CM03',
    solnorm = T,
    bandtype = 'B',
    gaintable = 'uid___A002_X85dcf7_Xc7c.ms.split.ap_pre_bandpass')
  
  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.bandpass_smooth40ch') 
  
  bandpass(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.bandpass_smooth40ch',
    field = '0', # J1700-2610
    solint = 'inf,40ch',
    combine = 'scan',
    refant = 'CM03',
    solnorm = T,
    bandtype = 'B',
    gaintable = 'uid___A002_X85dcf7_Xc7c.ms.split.ap_pre_bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_X85dcf7_Xc7c.ms.split.bandpass_smooth40ch', msName='uid___A002_X85dcf7_Xc7c.ms.split', interactive=False) 
  
  if applyonly != True: es.checkCalTable('uid___A002_X85dcf7_Xc7c.ms.split.bandpass', msName='uid___A002_X85dcf7_Xc7c.ms.split', interactive=False) 
  

# Save flags before gain cal
mystep = 14
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    mode = 'save',
    versionname = 'BeforeGainCalibration')
  
  

# Gain calibration
mystep = 15
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.phase_int') 
  gaincal(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.phase_int',
    field = '0~3', # J1700-2610,J1733-130,J1744-3116,SgrB2
    solint = 'int',
    refant = 'CM03',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X85dcf7_Xc7c.ms.split.bandpass_smooth40ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_X85dcf7_Xc7c.ms.split.phase_int', msName='uid___A002_X85dcf7_Xc7c.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.ampli_inf') 
  gaincal(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.ampli_inf',
    field = '0~3', # J1700-2610,J1733-130,J1744-3116,SgrB2
    solint = 'inf',
    refant = 'CM03',
    gaintype = 'T',
    calmode = 'a',
    gaintable = ['uid___A002_X85dcf7_Xc7c.ms.split.bandpass_smooth40ch', 'uid___A002_X85dcf7_Xc7c.ms.split.phase_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_X85dcf7_Xc7c.ms.split.ampli_inf', msName='uid___A002_X85dcf7_Xc7c.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.flux_inf') 
  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.fluxscale') 
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_X85dcf7_Xc7c.ms.split.fluxscale')
  
  fluxscale(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.ampli_inf',
    fluxtable = 'uid___A002_X85dcf7_Xc7c.ms.split.flux_inf',
    reference = '1') # J1733-130
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: es.fluxscale2(caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.ampli_inf', removeOutliers=True, msName='uid___A002_X85dcf7_Xc7c.ms', writeToFile=True, preavg=10000)
  
  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.phase_inf') 
  gaincal(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.phase_inf',
    field = '0~3', # J1700-2610,J1733-130,J1744-3116,SgrB2
    solint = 'inf',
    refant = 'CM03',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X85dcf7_Xc7c.ms.split.bandpass_smooth40ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_X85dcf7_Xc7c.ms.split.phase_inf', msName='uid___A002_X85dcf7_Xc7c.ms.split', interactive=False) 
  

# Save flags before applycal
mystep = 16
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    mode = 'save',
    versionname = 'BeforeApplycal')
  
  

# Application of the bandpass and gain cal tables
mystep = 17
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i in ['0', '1']: # J1700-2610,J1733-130
    applycal(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
      field = i,
      gaintable = ['uid___A002_X85dcf7_Xc7c.ms.split.bandpass_smooth40ch', 'uid___A002_X85dcf7_Xc7c.ms.split.phase_int', 'uid___A002_X85dcf7_Xc7c.ms.split.flux_inf'],
      gainfield = ['', i, i],
      interp = 'linear,linear',
      calwt = F,
      flagbackup = F)
  
  applycal(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    field = '2,4~55', # SgrB2
    gaintable = ['uid___A002_X85dcf7_Xc7c.ms.split.bandpass_smooth40ch', 'uid___A002_X85dcf7_Xc7c.ms.split.phase_inf', 'uid___A002_X85dcf7_Xc7c.ms.split.flux_inf'],
    gainfield = ['', '2', '2'], # J1744-3116
    interp = 'linear,linear',
    calwt = F,
    flagbackup = F)
  

# Split out corrected column
mystep = 18
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.cal') 
  split(vis = 'uid___A002_X85dcf7_Xc7c.ms.split',
    outputvis = 'uid___A002_X85dcf7_Xc7c.ms.split.cal',
    datacolumn = 'corrected',
    keepflags = T)
  
  

