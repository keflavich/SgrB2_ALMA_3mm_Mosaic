# ALMA Data Reduction Script

# Calibration

thesteps = []
step_title = {0: 'Import of the ASDM',
              1: 'Fix of SYSCAL table times',
              2: 'Running fixplanets on fields with 0,0 coordinates',
              3: 'listobs',
              4: 'A priori flagging',
              5: 'Generation and time averaging of the WVR cal table',
              6: 'Generation of the Tsys cal table',
              7: 'Generation of the antenna position cal table',
              8: 'Application of the WVR, Tsys and antpos cal tables',
              9: 'Split out science SPWs and time average',
              10: 'Listobs, clear pointing table, and save original flags',
              11: 'Initial flagging',
              12: 'Putting a model for the flux calibrator(s)',
              13: 'Save flags before bandpass cal',
              14: 'Bandpass calibration',
              15: 'Save flags before gain cal',
              16: 'Gain calibration',
              17: 'Save flags before applycal',
              18: 'Application of the bandpass and gain cal tables',
              19: 'Split out corrected column',
              20: 'Save flags after applycal'}

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


#if re.search('^4.3.1', casadef.casa_version) == None:
# sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.3.1')


# CALIBRATE_AMPLI: Venus
# CALIBRATE_ATMOSPHERE: J1700-2610,SgrB2,Venus
# CALIBRATE_BANDPASS: J1700-2610
# CALIBRATE_FLUX: Venus
# CALIBRATE_FOCUS: 
# CALIBRATE_PHASE: J1744-3116
# CALIBRATE_POINTING: J1700-2610
# OBSERVE_TARGET: SgrB2

# Using reference antenna = DV18

# Import of the ASDM
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  if os.path.exists('uid___A002_X95e355_X1f13.ms') == False:
    importasdm('uid___A002_X95e355_X1f13', asis='Antenna Station Receiver Source CalAtmosphere CalWVR', bdfflags=True, lazy=False)
  if applyonly != True: es.fixForCSV2555('uid___A002_X95e355_X1f13.ms')

# Fix of SYSCAL table times
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  from recipes.almahelpers import fixsyscaltimes
  fixsyscaltimes(vis = 'uid___A002_X95e355_X1f13.ms')

print "# A priori calibration"

# Running fixplanets on fields with 0,0 coordinates
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  fixplanets(vis = 'uid___A002_X95e355_X1f13.ms',
    field = '1', # Venus
    fixuvw = T)
  

# listobs
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X95e355_X1f13.ms.listobs')
  listobs(vis = 'uid___A002_X95e355_X1f13.ms',
    listfile = 'uid___A002_X95e355_X1f13.ms.listobs')
  
  

# A priori flagging
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  flagdata(vis = 'uid___A002_X95e355_X1f13.ms',
    mode = 'manual',
    spw = '1~24',
    autocorr = T,
    flagbackup = F)
  
  flagdata(vis = 'uid___A002_X95e355_X1f13.ms',
    mode = 'manual',
    intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*',
    flagbackup = F)
  
  flagcmd(vis = 'uid___A002_X95e355_X1f13.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'plot',
    plotfile = 'uid___A002_X95e355_X1f13.ms.flagcmd.png')
  
  flagcmd(vis = 'uid___A002_X95e355_X1f13.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'apply')
  

# Generation and time averaging of the WVR cal table
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X95e355_X1f13.ms.wvr') 
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.wvrgcal') 
  
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_X95e355_X1f13.ms.wvrgcal')
  
  wvrgcal(vis = 'uid___A002_X95e355_X1f13.ms',
    caltable = 'uid___A002_X95e355_X1f13.ms.wvr',
    spw = [17, 19, 21, 23],
    smooth = '6.048s',
    toffset = 0,
    tie = ['SgrB2,J1744-3116'],
    statsource = 'SgrB2')
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: aU.plotWVRSolutions(caltable='uid___A002_X95e355_X1f13.ms.wvr', spw='17', antenna='DV18',
    yrange=[-199,199],subplot=22, interactive=False,
    figfile='uid___A002_X95e355_X1f13.ms.wvr.plots/uid___A002_X95e355_X1f13.ms.wvr') 
  
  #Note: If you see wraps in these plots, try changing yrange or unwrap=True 
  #Note: If all plots look strange, it may be a bad WVR on the reference antenna.
  #      To check, you can set antenna='' to show all baselines.
  

# Generation of the Tsys cal table
mystep = 6
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X95e355_X1f13.ms.tsys') 
  gencal(vis = 'uid___A002_X95e355_X1f13.ms',
    caltable = 'uid___A002_X95e355_X1f13.ms.tsys',
    caltype = 'tsys')
  
  # Flagging edge channels
  
  flagdata(vis = 'uid___A002_X95e355_X1f13.ms.tsys',
    mode = 'manual',
    spw = '9:0~3;124~127,11:0~3;124~127,13:0~3;124~127,15:0~3;124~127',
    flagbackup = F)
  
  if applyonly != True: aU.plotbandpass(caltable='uid___A002_X95e355_X1f13.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='92.1875%',showfdm=True, showBasebandNumber=True, 
    field='', figfile='uid___A002_X95e355_X1f13.ms.tsys.plots.overlayTime/uid___A002_X95e355_X1f13.ms.tsys') 
  
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.tsys', msName='uid___A002_X95e355_X1f13.ms', interactive=False) 
  

# Generation of the antenna position cal table
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Position for antenna DA63 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DA48 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DA61 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Note: no baseline run found for antenna DA45.
  
  # Note: no baseline run found for antenna DV13.
  
  # Note: no baseline run found for antenna DA41.
  
  # Position for antenna DA43 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DV15 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DA62 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DV11 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DV10 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DV22 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DV04 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DA52 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DA50 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Note: no baseline run found for antenna DV24.
  
  # Note: no baseline run found for antenna DA54.
  
  # Note: no baseline run found for antenna DA55.
  
  # Note: no baseline run found for antenna DV07.
  
  # Position for antenna DV17 is derived from baseline run made on 2014-12-02 07:19:25.
  
  # Position for antenna DV16 is derived from baseline run made on 2014-12-02 07:19:25.
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.antpos') 
  gencal(vis = 'uid___A002_X95e355_X1f13.ms',
    caltable = 'uid___A002_X95e355_X1f13.ms.antpos',
    caltype = 'antpos',
    antenna = 'DA52,DV22,DA63,DA48,DA61,DV11,DV10,DV04,DV15,DA50,DA43,DV16,DV17,DA62',
    parameter = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  #  parameter = [3.80642939475e-05,0.000208067074874,-0.000100619093019,-0.000438846177427,0.000608006242819,0.000170992518114,-0.000578157130418,0.000314672396318,7.21083597271e-05,8.36793333292e-07,7.72997736931e-08,-1.10128894448e-06,-0.000476244700753,0.000646645340667,0.000222433784303,-0.000503372919646,0.000499006970657,0.000187546652871,-0.000411181905591,0.000281534048615,0.000194200978332,-0.000597593147508,0.000539165930576,0.000205271842583,4.27011400461e-07,-5.01982867718e-07,-5.47152012587e-07,7.03148543835e-07,1.39698386192e-08,1.25262886286e-07,3.88827174902e-07,-9.77888703346e-08,-6.7101791501e-07,1.17253512144e-06,9.49949026108e-08,-6.10016286373e-07,-2.06753611565e-07,-5.96046447754e-08,-1.28522515297e-07,5.33182173967e-07,-8.19563865662e-07,-7.23171979189e-07])
  

# Application of the WVR, Tsys and antpos cal tables
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  
  from recipes.almahelpers import tsysspwmap
  tsysmap = tsysspwmap(vis = 'uid___A002_X95e355_X1f13.ms', tsystable = 'uid___A002_X95e355_X1f13.ms.tsys', tsysChanTol = 1)
  
  
  
  applycal(vis = 'uid___A002_X95e355_X1f13.ms',
    field = '0',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X95e355_X1f13.ms.tsys', 'uid___A002_X95e355_X1f13.ms.wvr', 'uid___A002_X95e355_X1f13.ms.antpos'],
    gainfield = ['0', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_X95e355_X1f13.ms',
    field = '1',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X95e355_X1f13.ms.tsys', 'uid___A002_X95e355_X1f13.ms.wvr', 'uid___A002_X95e355_X1f13.ms.antpos'],
    gainfield = ['1', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1744-3116 didn't have any Tsys measurement, so I used the one made on SgrB2. This is probably Ok.
  
  applycal(vis = 'uid___A002_X95e355_X1f13.ms',
    field = '2',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X95e355_X1f13.ms.tsys', 'uid___A002_X95e355_X1f13.ms.wvr', 'uid___A002_X95e355_X1f13.ms.antpos'],
    gainfield = ['3', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_X95e355_X1f13.ms',
    field = '3~151',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_X95e355_X1f13.ms.tsys', 'uid___A002_X95e355_X1f13.ms.wvr', 'uid___A002_X95e355_X1f13.ms.antpos'],
    gainfield = ['3', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  if applyonly != True: es.getCalWeightStats('uid___A002_X95e355_X1f13.ms') 
  

# Split out science SPWs and time average
mystep = 9
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split') 
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.flagversions') 
  split(vis = 'uid___A002_X95e355_X1f13.ms',
    outputvis = 'uid___A002_X95e355_X1f13.ms.split',
    datacolumn = 'corrected',
    spw = '17,19,21,23',
    keepflags = T)
  
  

print "# Calibration"

# Listobs, clear pointing table, and save original flags
mystep = 10
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.listobs')
  listobs(vis = 'uid___A002_X95e355_X1f13.ms.split',
    listfile = 'uid___A002_X95e355_X1f13.ms.split.listobs')
  
  tb.open('uid___A002_X95e355_X1f13.ms.split/POINTING', nomodify = False)
  a = tb.rownumbers()
  tb.removerows(a)
  tb.close()
  
  if not os.path.exists('uid___A002_X95e355_X1f13.ms.split.flagversions/Original.flags'):
    flagmanager(vis = 'uid___A002_X95e355_X1f13.ms.split',
      mode = 'save',
      versionname = 'Original')
  
  

# Initial flagging
mystep = 11
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Flagging shadowed data
  
  flagdata(vis = 'uid___A002_X95e355_X1f13.ms.split',
    mode = 'shadow',
    flagbackup = F)
  
  

# Putting a model for the flux calibrator(s)
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  setjy(vis = 'uid___A002_X95e355_X1f13.ms.split',
    field = '1', # Venus
    spw = '0,1,2,3',
    standard = 'Butler-JPL-Horizons 2012')
  
  if applyonly != True:
    os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.setjy.field*.png') 
    for i in ['1']:
      plotms(vis = 'uid___A002_X95e355_X1f13.ms.split',
        xaxis = 'uvdist',
        yaxis = 'amp',
        ydatacolumn = 'model',
        field = str(i),
        spw = '0,1,2,3',
        avgchannel = '9999',
        coloraxis = 'spw',
        plotfile = 'uid___A002_X95e355_X1f13.ms.split.setjy.field'+i+'.png')
  

# Save flags before bandpass cal
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X95e355_X1f13.ms.split',
    mode = 'save',
    versionname = 'BeforeBandpassCalibration')
  
  

# Bandpass calibration
mystep = 14
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.ap_pre_bandpass') 
  
  gaincal(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.ap_pre_bandpass',
    field = '0', # J1700-2610
    spw = '0:3072~4608,1:3072~4608,2:3072~4608,3:3072~4608',
    scan = '1,2,4',
    solint = 'int',
    refant = 'DV18',
    calmode = 'p')
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.split.ap_pre_bandpass', msName='uid___A002_X95e355_X1f13.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.bandpass') 
  bandpass(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.bandpass',
    field = '0', # J1700-2610
    scan = '1,2,4',
    solint = 'inf',
    combine = 'scan',
    refant = 'DV18',
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_X95e355_X1f13.ms.split.ap_pre_bandpass')
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch') 
  
  bandpass(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch',
    field = '0', # J1700-2610
    scan = '1,2,4',
    solint = 'inf,20ch',
    combine = 'scan',
    refant = 'DV18',
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_X95e355_X1f13.ms.split.ap_pre_bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch', msName='uid___A002_X95e355_X1f13.ms.split', interactive=False) 
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.split.bandpass', msName='uid___A002_X95e355_X1f13.ms.split', interactive=False) 
  

# Save flags before gain cal
mystep = 15
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X95e355_X1f13.ms.split',
    mode = 'save',
    versionname = 'BeforeGainCalibration')
  
  

# Gain calibration
mystep = 16
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Note: the Solar system object used for flux calibration is highly resolved on some baselines.
  # Note: we will first determine the flux of the phase calibrator(s) on a subset of antennas.
  
  delmod('uid___A002_X95e355_X1f13.ms.split',field='2')
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.phase_short_int') 
  gaincal(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.phase_short_int',
    field = '1', # Venus
    selectdata = T,
    antenna = 'DV04,DV18&',
    solint = 'int',
    refant = 'DV18',
    minblperant = 1,
    minsnr = 2.0,
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch')
  
  gaincal(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.phase_short_int',
    field = '0,2', # J1700-2610,J1744-3116
    selectdata = T,
    solint = 'int',
    refant = 'DV18',
    minblperant = 1,
    minsnr = 2.0,
    gaintype = 'G',
    calmode = 'p',
    append = T,
    gaintable = 'uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.split.phase_short_int', msName='uid___A002_X95e355_X1f13.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.ampli_short_inf') 
  gaincal(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.ampli_short_inf',
    field = '0,1,2', # J1700-2610,Venus,J1744-3116
    selectdata = T,
    solint = 'inf',
    refant = 'DV18',
    minblperant = 1,
    minsnr = 2.0,
    gaintype = 'T',
    calmode = 'a',
    gaintable = ['uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch', 'uid___A002_X95e355_X1f13.ms.split.phase_short_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.split.ampli_short_inf', msName='uid___A002_X95e355_X1f13.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.flux_short_inf') 
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.fluxscale') 
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_X95e355_X1f13.ms.split.fluxscale')
  
  fluxscaleDict = fluxscale(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.ampli_short_inf',
    fluxtable = 'uid___A002_X95e355_X1f13.ms.split.flux_short_inf',
    reference = '1') # Venus
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: es.fluxscale2(caltable = 'uid___A002_X95e355_X1f13.ms.split.ampli_short_inf', removeOutliers=True, msName='uid___A002_X95e355_X1f13.ms', writeToFile=True, preavg=10000)
  
  f = open('uid___A002_X95e355_X1f13.ms.split.fluxscale')
  fc = f.readlines()
  f.close()
  
  for phaseCalName in ['J1744-3116']:
    for i in range(len(fc)):
      if fc[i].find('Flux density for '+phaseCalName) != -1 and re.search('in SpW=[0-9]+(?: \(.*?\))? is: [0-9]+\.[0-9]+', fc[i], re.DOTALL|re.IGNORECASE) != None:
        line = (re.search('in SpW=[0-9]+(?: \(.*?\))? is: [0-9]+\.[0-9]+', fc[i], re.DOTALL|re.IGNORECASE)).group(0)
        spwId = (line.split('='))[1].split()[0]
        flux = float((line.split(':'))[1].split()[0])
        setjy(vis = 'uid___A002_X95e355_X1f13.ms.split',
          field = phaseCalName.replace(';','*;').split(';')[0],
          spw = spwId,
          standard = 'manual',
          fluxdensity = [flux,0,0,0])
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.phase_int') 
  gaincal(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.phase_int',
    field = '0,1,2', # J1700-2610,Venus,J1744-3116
    solint = 'int',
    refant = 'DV18',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.split.phase_int', msName='uid___A002_X95e355_X1f13.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.flux_inf') 
  gaincal(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.flux_inf',
    field = '0,1,2', # J1700-2610,Venus,J1744-3116
    solint = 'inf',
    refant = 'DV18',
    gaintype = 'T',
    calmode = 'a',
    gaintable = ['uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch', 'uid___A002_X95e355_X1f13.ms.split.phase_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.split.flux_inf', msName='uid___A002_X95e355_X1f13.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.phase_inf') 
  gaincal(vis = 'uid___A002_X95e355_X1f13.ms.split',
    caltable = 'uid___A002_X95e355_X1f13.ms.split.phase_inf',
    field = '0,1,2', # J1700-2610,Venus,J1744-3116
    solint = 'inf',
    refant = 'DV18',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_X95e355_X1f13.ms.split.phase_inf', msName='uid___A002_X95e355_X1f13.ms.split', interactive=False) 
  

# Save flags before applycal
mystep = 17
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X95e355_X1f13.ms.split',
    mode = 'save',
    versionname = 'BeforeApplycal')
  
  

# Application of the bandpass and gain cal tables
mystep = 18
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i in ['0', '1']: # J1700-2610,Venus
    applycal(vis = 'uid___A002_X95e355_X1f13.ms.split',
      field = str(i),
      gaintable = ['uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch', 'uid___A002_X95e355_X1f13.ms.split.phase_int', 'uid___A002_X95e355_X1f13.ms.split.flux_inf'],
      gainfield = ['', i, i],
      interp = 'linear,linear',
      calwt = T,
      flagbackup = F)
  
  applycal(vis = 'uid___A002_X95e355_X1f13.ms.split',
    field = '2,3~151', # SgrB2
    gaintable = ['uid___A002_X95e355_X1f13.ms.split.bandpass_smooth20ch', 'uid___A002_X95e355_X1f13.ms.split.phase_inf', 'uid___A002_X95e355_X1f13.ms.split.flux_inf'],
    gainfield = ['', '2', '2'], # J1744-3116
    interp = 'linear,linear',
    calwt = T,
    flagbackup = F)
  

# Split out corrected column
mystep = 19
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.cal') 
  os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.cal.flagversions') 
  split(vis = 'uid___A002_X95e355_X1f13.ms.split',
    outputvis = 'uid___A002_X95e355_X1f13.ms.split.cal',
    datacolumn = 'corrected',
    keepflags = T)
  
  

# Save flags after applycal
mystep = 20
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
    mode = 'save',
    versionname = 'AfterApplycal')
  
  

