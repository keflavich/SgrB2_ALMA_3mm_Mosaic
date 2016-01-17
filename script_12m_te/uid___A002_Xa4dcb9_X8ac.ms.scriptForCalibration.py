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
              9: 'Listobs, clear pointing table and save original flags',
              10: 'Initial flagging',
              11: 'Putting a model for the flux calibrator(s)',
              12: 'Save flags before bandpass cal',
              13: 'Bandpass calibration',
              14: 'Save flags before gain cal',
              15: 'Gain calibration',
              16: 'Save flags before applycal',
              17: 'Application of the bandpass and gain cal tables',
              18: 'Split out corrected column',
              19: 'Save flags after applycal'}

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


#if re.search('^4.4.0', casadef.casa_version) == None:
# sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.4.0')


# CALIBRATE_AMPLI: Titan
# CALIBRATE_ATMOSPHERE: J1733-1304,SgrB2,Titan
# CALIBRATE_BANDPASS: J1733-1304
# CALIBRATE_FLUX: Titan
# CALIBRATE_FOCUS: 
# CALIBRATE_PHASE: J1752-2956
# CALIBRATE_POINTING: J1733-1304
# OBSERVE_CHECK: 
# OBSERVE_TARGET: SgrB2

# Using reference antenna = DA55 chnage to DV04

# Import of the ASDM
mystep = 0
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  if os.path.exists('uid___A002_Xa4dcb9_X8ac.ms') == False:
    importasdm('uid___A002_Xa4dcb9_X8ac', asis='Antenna Station Receiver Source CalAtmosphere CalWVR CorrelatorMode SBSummary', bdfflags=True, lazy=True, process_caldevice=False)
  if applyonly != True: es.fixForCSV2555('uid___A002_Xa4dcb9_X8ac.ms')

# Fix of SYSCAL table times
mystep = 1
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  from recipes.almahelpers import fixsyscaltimes
  fixsyscaltimes(vis = 'uid___A002_Xa4dcb9_X8ac.ms')

print "# A priori calibration"

# listobs
mystep = 2
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.listobs')
  listobs(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    listfile = 'uid___A002_Xa4dcb9_X8ac.ms.listobs')
  
  

# A priori flagging
mystep = 3
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    mode = 'manual',
    spw = '1~24',
    autocorr = T,
    flagbackup = F)
  
  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    mode = 'manual',
    intent = '*POINTING*,*SIDEBAND_RATIO*,*ATMOSPHERE*',
    flagbackup = F)
  
  flagcmd(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'plot',
    plotfile = 'uid___A002_Xa4dcb9_X8ac.ms.flagcmd.png')
  
  flagcmd(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    inpmode = 'table',
    useapplied = True,
    action = 'apply')

# Generation and time averaging of the WVR cal table
mystep = 4
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.wvr') 
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.wvrgcal') 
  
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_Xa4dcb9_X8ac.ms.wvrgcal')
  
  wvrgcal(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.wvr',
    spw = [17, 19, 21, 23],
    smooth = '6.048s',
    toffset = 0,
    tie = ['SgrB2,J1752-2956'],
    statsource = 'SgrB2')
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: aU.plotWVRSolutions(caltable='uid___A002_Xa4dcb9_X8ac.ms.wvr', spw='17', antenna='DV04',
    yrange=[-199,199],subplot=22, interactive=False,
    figfile='uid___A002_Xa4dcb9_X8ac.ms.wvr.plots/uid___A002_Xa4dcb9_X8ac.ms.wvr') 
  
  #Note: If you see wraps in these plots, try changing yrange or unwrap=True 
  #Note: If all plots look strange, it may be a bad WVR on the reference antenna.
  #      To check, you can set antenna='' to show all baselines.
  

# Generation of the Tsys cal table
mystep = 5
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.tsys') 
  gencal(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.tsys',
    caltype = 'tsys')
  
  # Flagging edge channels
  
  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.tsys',
    mode = 'manual',
    spw = '9:0~3;124~127,11:0~3;124~127,13:0~3;124~127,15:0~3;124~127',
    flagbackup = F)


#   Flag part of data with bad Tsys data

  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.tsys',  # high Tsys values for antenna DA41 in this scan
    mode = 'manual',
    antenna='DA41',
    scan='18',
    spw='9,11',       
    flagbackup=F)

  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.tsys',  # large oscilation in the end of the freq. range
    mode = 'manual',
    antenna='DV03',
    spw='11:0~30', #'11:90.0~90.5GHz',       
    flagbackup=F)
 
  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.tsys',  # large oscilation in the end of the freq. range
    mode = 'manual',
    antenna='DA51',
    spw='11:0~14', #'11:90.25~90.5GHz',       
    flagbackup=F)

  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.tsys',  # large oscilation in the end of the freq. range
    mode = 'manual',
    antenna='DA43,DA52,DA58,DV22,DV24',
    spw='13:0~11', #'13:100.3~100.55GHz',       
    flagbackup=F)

  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.tsys',  # large oscilation in the end of the freq. range
    mode = 'manual',
    antenna='DA54,DV16',
    spw='13:0~14', #'13:100.3~100.6GHz',       
    flagbackup=F)


  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.tsys',  # large spike present in almost all antennas
    mode = 'manual',
    spw='9:38~41', #'9:91.63~91.69GHz',
    antenna='DA42,DA43,DA44,DA45,DA46,DA47,DA54,DA55,DV01,DV05,DV16,DV25',       
    flagbackup=F)

  
  if applyonly != True: aU.plotbandpass(caltable='uid___A002_Xa4dcb9_X8ac.ms.tsys', overlay='time', 
    xaxis='freq', yaxis='amp', subplot=22, buildpdf=False, interactive=False,
    showatm=True,pwv='auto',chanrange='92.1875%',showfdm=True, showBasebandNumber=True, showimage=False, 
    field='', figfile='uid___A002_Xa4dcb9_X8ac.ms.tsys.plots.overlayTime/uid___A002_Xa4dcb9_X8ac.ms.tsys') 
  
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.tsys', msName='uid___A002_Xa4dcb9_X8ac.ms', interactive=False) 
  

# Generation of the antenna position cal table
mystep = 6
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Position for antenna DV19 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Position for antenna DV18 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Position for antenna DA64 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Note: no baseline run found for antenna DA62.
  
  # Position for antenna DA45 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Note: no baseline run found for antenna DA44.
  
  # Note: no baseline run found for antenna DA47.
  
  # Position for antenna DV14 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Position for antenna DA63 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Note: no baseline run found for antenna DV15.
  
  # Note: no baseline run found for antenna DV11.
  
  # Position for antenna DV10 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Note: no baseline run found for antenna DV22.
  
  # Note: no baseline run found for antenna DV04.
  
  # Note: no baseline run found for antenna DA52.
  
  # Position for antenna DA53 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Note: no baseline run found for antenna DA51.
  
  # Note: no baseline run found for antenna DV24.
  
  # Note: no baseline run found for antenna DA57.
  
  # Note: no baseline run found for antenna DA55.
  
  # Note: no baseline run found for antenna DV07.
  
  # Position for antenna DA58 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Note: no baseline run found for antenna DV03.
  
  # Note: no baseline run found for antenna DV01.
  
  # Position for antenna DV17 is derived from baseline run made on 2015-06-18 07:40:18.
  
  # Note: no baseline run found for antenna DV05.
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.antpos') 
  gencal(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.antpos',
    caltype = 'antpos',
    antenna = 'DA45,DA53,DA58,DA63,DA64,DV10,DV14,DV17,DV18,DV19',
    parameter = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  #  parameter = [-2.38629e-04,6.51773e-04,2.24732e-04,-2.70048e-05,-1.74361e-04,-2.25566e-04,3.44785e-04,-4.44180e-04,-5.83638e-04,7.35842e-05,-1.11237e-04,-3.48316e-04,2.61424e-04,-8.20123e-04,-4.17703e-04,1.93743e-05,2.94065e-05,-3.65355e-04,1.60241e-04,-9.06593e-05,-3.67151e-04,3.14085e-04,-7.30734e-05,-5.97357e-04,-1.65702e-04,4.36814e-04,-1.03487e-04,3.04381e-04,-4.77407e-04,-4.33876e-04])
  
  
  # antenna x_offset y_offset z_offset total_offset baseline_date
  # DA64     2.61424e-04   -8.20123e-04   -4.17703e-04    9.56776e-04      2015-06-18 07:40:18
  # DA58     3.44785e-04   -4.44180e-04   -5.83638e-04    8.10436e-04      2015-06-18 07:40:18
  # DA45    -2.38629e-04    6.51773e-04    2.24732e-04    7.29559e-04      2015-06-18 07:40:18
  # DV19     3.04381e-04   -4.77407e-04   -4.33876e-04    7.13312e-04      2015-06-18 07:40:18
  # DV17     3.14085e-04   -7.30734e-05   -5.97357e-04    6.78841e-04      2015-06-18 07:40:18
  # DV18    -1.65702e-04    4.36814e-04   -1.03487e-04    4.78511e-04      2015-06-18 07:40:18
  # DV14     1.60241e-04   -9.06593e-05   -3.67151e-04    4.10726e-04      2015-06-18 07:40:18
  # DA63     7.35842e-05   -1.11237e-04   -3.48316e-04    3.72977e-04      2015-06-18 07:40:18
  # DV10     1.93743e-05    2.94065e-05   -3.65355e-04    3.67048e-04      2015-06-18 07:40:18
  # DA53    -2.70048e-05   -1.74361e-04   -2.25566e-04    2.86376e-04      2015-06-18 07:40:18
  

# Application of the WVR, Tsys and antpos cal tables
mystep = 7
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  
  from recipes.almahelpers import tsysspwmap
  tsysmap = tsysspwmap(vis = 'uid___A002_Xa4dcb9_X8ac.ms', tsystable = 'uid___A002_Xa4dcb9_X8ac.ms.tsys', tsysChanTol = 1)
  
    
  applycal(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    field = '0',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_Xa4dcb9_X8ac.ms.tsys', 'uid___A002_Xa4dcb9_X8ac.ms.wvr', 'uid___A002_Xa4dcb9_X8ac.ms.antpos'],
    gainfield = ['0', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  applycal(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    field = '1',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_Xa4dcb9_X8ac.ms.tsys', 'uid___A002_Xa4dcb9_X8ac.ms.wvr', 'uid___A002_Xa4dcb9_X8ac.ms.antpos'],
    gainfield = ['1', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  # Note: J1752-2956 didn't have any Tsys measurement, so I used the one made on SgrB2. This is probably Ok.
  
  applycal(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    field = '2',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_Xa4dcb9_X8ac.ms.tsys', 'uid___A002_Xa4dcb9_X8ac.ms.wvr', 'uid___A002_Xa4dcb9_X8ac.ms.antpos'],
    gainfield = ['3', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  
  
  applycal(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    field = '3~151',
    spw = '17,19,21,23',
    gaintable = ['uid___A002_Xa4dcb9_X8ac.ms.tsys', 'uid___A002_Xa4dcb9_X8ac.ms.wvr', 'uid___A002_Xa4dcb9_X8ac.ms.antpos'],
    gainfield = ['3', '', ''],
    interp = 'linear,linear',
    spwmap = [tsysmap,[],[]],
    calwt = T,
    flagbackup = F)
  
  if applyonly != True: es.getCalWeightStats('uid___A002_Xa4dcb9_X8ac.ms') 
  

# Split out science SPWs and time average
mystep = 8
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split') 
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.flagversions') 
  split(vis = 'uid___A002_Xa4dcb9_X8ac.ms',
    outputvis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    datacolumn = 'corrected',
    spw = '17,19,21,23',
    keepflags = T)
  
  

print "# Calibration"

# Listobs, clear pointing table and save original flags
mystep = 9
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.listobs')
  listobs(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    listfile = 'uid___A002_Xa4dcb9_X8ac.ms.split.listobs')
  
  tb.open('uid___A002_Xa4dcb9_X8ac.ms.split/POINTING', nomodify = False)
  a = tb.rownumbers()
  tb.removerows(a)
  tb.close()
  
  if not os.path.exists('uid___A002_Xa4dcb9_X8ac.ms.split.flagversions/Original.flags'):
    flagmanager(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
      mode = 'save',
      versionname = 'Original')
  
  

# Initial flagging
mystep = 10
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Flagging shadowed data
  
  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    mode = 'shadow',
    flagbackup = F)
  
  # Flagging atmospheric line(s)
  
  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    mode = 'manual',
    spw = '1:6623~7679',
    field = '1', #Titan
    flagbackup = F)

  # Flagging spike due to WVR local oscillator leakage, for all sources
  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    mode = 'manual',
    spw = '0:2174~2388',
    field='',
    flagbackup=F)
    
# Flagging features in the bandpass calibrator
#  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
#    mode = 'manual',
#    spw='1:4969~5067;7279~7336',
#    field='0',
#    flagbackup=F)

# high values amplitude scan of phase calibrator

  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split', 
    mode = 'manual',
    field = '2', #phase cal
    timerange= '03:15:50~03:15:55',       
    flagbackup = F)

# phase outliers in phase cal

  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split', 
    mode = 'manual',
    field = '2', #phase cal
    antenna='DV22,DV23',
    timerange= '03:25:10~03:25:18',       
    flagbackup = F)


# Putting a model for the flux calibrator(s)
mystep = 11
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  setjy(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    field = '1', # Titan
    spw = '0,1,2,3',
    standard = 'Butler-JPL-Horizons 2012')
  
  if applyonly != True:
    os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.setjy.field*.png') 
    for i in ['1']:
      plotms(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
        xaxis = 'uvdist',
        yaxis = 'amp',
        ydatacolumn = 'model',
        field = str(i),
        spw = '0,1,2,3',
        avgchannel = '9999',
        coloraxis = 'spw',
        plotfile = 'uid___A002_Xa4dcb9_X8ac.ms.split.setjy.field'+i+'.png')
  

# Save flags before bandpass cal
mystep = 12
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    mode = 'save',
    versionname = 'BeforeBandpassCalibration')
  
  

# Bandpass calibration
mystep = 13
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.ap_pre_bandpass') 
  
  gaincal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.ap_pre_bandpass',
    field = '0', # J1733-1304
    spw = '0:3072~4608,1:3072~4608,2:3072~4608,3:3072~4608',
    scan = '1,2,4',
    solint = 'int',
    refant = 'DV04',
    calmode = 'p')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.split.ap_pre_bandpass', msName='uid___A002_Xa4dcb9_X8ac.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.bandpass') 
  bandpass(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.bandpass',
    field = '0', # J1733-1304
    scan = '1,2,4',
    solint = 'inf',
    combine = 'scan',
    refant = 'DV04',
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_Xa4dcb9_X8ac.ms.split.ap_pre_bandpass')
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch') 
  
  bandpass(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch',
    field = '0', # J1733-1304
    scan = '1,2,4',
    solint = 'inf,20ch',
    combine = 'scan',
    refant = 'DV04',
    solnorm = True,
    bandtype = 'B',
    gaintable = 'uid___A002_Xa4dcb9_X8ac.ms.split.ap_pre_bandpass')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch', msName='uid___A002_Xa4dcb9_X8ac.ms.split', interactive=False) 
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.split.bandpass', msName='uid___A002_Xa4dcb9_X8ac.ms.split', interactive=False) 
  

# Save flags before gain cal
mystep = 14
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    mode = 'save',
    versionname = 'BeforeGainCalibration')
  
  
# Titan still present

# Gain calibration
mystep = 15
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  # Note: the Solar system object used for flux calibration is highly resolved on some baselines.
  # Note: we will first determine the flux of the phase calibrator(s) on a subset of antennas.
  
  delmod('uid___A002_Xa4dcb9_X8ac.ms.split',field='2')
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.phase_short_int') 
  gaincal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.phase_short_int',
    field = '1', # Titan
    selectdata = T,
    antenna = 'DA41,DA42,DA43,DA46,DA49,DA50,DA51,DA53,DA54,DA55,DA58,DA59,DA60,DA61,DA63,DV01,DV04,DV05,DV13,DV14,DV15,DV16,DV18,DV19,DV20,DV21,DV23,DV24,DV25&',
    solint = 'int',
    refant = 'DV04',
    gaintype = 'G',
    calmode = 'p',
    uvrange= '0~950m', #Titan resolved, see amp vs uvdistance plot 
    gaintable = 'uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch')
  
  gaincal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.phase_short_int',
    field = '0,2', # J1733-1304,J1752-2956
    selectdata = T,
    solint = 'int',
    refant = 'DV04',
    gaintype = 'G',
    calmode = 'p',
    append = T,
    gaintable = 'uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.split.phase_short_int', msName='uid___A002_Xa4dcb9_X8ac.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.ampli_short_inf') 
  gaincal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.ampli_short_inf',
    field = '0,1,2', # J1733-1304,Titan,J1752-2956
    selectdata = T,
    solint = 'inf',
    refant = 'DV04',
    gaintype = 'G', # was T
    calmode = 'ap', #was a
    gaintable = ['uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch', 'uid___A002_Xa4dcb9_X8ac.ms.split.phase_short_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.split.ampli_short_inf', msName='uid___A002_Xa4dcb9_X8ac.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.flux_short_inf') 
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.fluxscale') 
  mylogfile = casalog.logfile()
  casalog.setlogfile('uid___A002_Xa4dcb9_X8ac.ms.split.fluxscale')
  
  fluxscaleDict = fluxscale(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.ampli_short_inf',
    fluxtable = 'uid___A002_Xa4dcb9_X8ac.ms.split.flux_short_inf',
    reference = '1') # Titan
  
  casalog.setlogfile(mylogfile)
  
  if applyonly != True: es.fluxscale2(caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.ampli_short_inf', removeOutliers=True, msName='uid___A002_Xa4dcb9_X8ac.ms', writeToFile=True, preavg=10000)
  
  f = open('uid___A002_Xa4dcb9_X8ac.ms.split.fluxscale')
  fc = f.readlines()
  f.close()
  
  for phaseCalName in ['J1752-2956']:
    for i in range(len(fc)):
      if fc[i].find('Flux density for '+phaseCalName) != -1 and re.search('in SpW=[0-9]+(?: \(.*?\))? is: [0-9]+\.[0-9]+', fc[i], re.DOTALL|re.IGNORECASE) != None:
        line = (re.search('in SpW=[0-9]+(?: \(.*?\))? is: [0-9]+\.[0-9]+', fc[i], re.DOTALL|re.IGNORECASE)).group(0)
        spwId = (line.split('='))[1].split()[0]
        flux = float((line.split(':'))[1].split()[0])
        setjy(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
          field = phaseCalName.replace(';','*;').split(';')[0],
          spw = spwId,
          standard = 'manual',
          fluxdensity = [flux,0,0,0])
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.phase_int') 
  gaincal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.phase_int',
    field = '0,1,2', # J1733-1304,Titan,J1752-2956
    solint = 'int',
    refant = 'DV04',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.split.phase_int', msName='uid___A002_Xa4dcb9_X8ac.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.flux_inf') 
  gaincal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.flux_inf',
    field = '0,1,2', # J1733-1304,Titan,J1752-2956
    solint = 'inf',
    refant = 'DV04',
    gaintype = 'G', #was T
    calmode = 'ap', #was a
    gaintable = ['uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch', 'uid___A002_Xa4dcb9_X8ac.ms.split.phase_int'])
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.split.flux_inf', msName='uid___A002_Xa4dcb9_X8ac.ms.split', interactive=False) 
  
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.phase_inf') 
  gaincal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    caltable = 'uid___A002_Xa4dcb9_X8ac.ms.split.phase_inf',
    field = '0,1,2', # J1733-1304,Titan,J1752-2956
    solint = 'inf',
    refant = 'DV04',
    gaintype = 'G',
    calmode = 'p',
    gaintable = 'uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch')
  
  if applyonly != True: es.checkCalTable('uid___A002_Xa4dcb9_X8ac.ms.split.phase_inf', msName='uid___A002_Xa4dcb9_X8ac.ms.split', interactive=False) 
  

# Save flags before applycal, 4/12 titan dissapearance, continue from here
mystep = 16
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    mode = 'save',
    versionname = 'BeforeApplycal')
  
  

# Application of the bandpass and gain cal tables
mystep = 17
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  for i in ['0', '1']: # J1733-1304,Titan
    applycal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
      field = str(i),
      gaintable = ['uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch', 'uid___A002_Xa4dcb9_X8ac.ms.split.phase_int', 'uid___A002_Xa4dcb9_X8ac.ms.split.flux_inf'],
      gainfield = ['', i, i],
      interp = 'linear,linear',
      calwt = T,
      flagbackup = F)
  
  applycal(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    field = '2,3~151', # SgrB2
    gaintable = ['uid___A002_Xa4dcb9_X8ac.ms.split.bandpass_smooth20ch', 'uid___A002_Xa4dcb9_X8ac.ms.split.phase_inf', 'uid___A002_Xa4dcb9_X8ac.ms.split.flux_inf'],
    gainfield = ['', '2', '2'], # J1752-2956
    interp = 'linear,linear',
    calwt = T,
    flagbackup = F)
  

# Split out corrected column
mystep = 18
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.cal') 
  os.system('rm -rf uid___A002_Xa4dcb9_X8ac.ms.split.cal.flagversions') 
  split(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split',
    outputvis = 'uid___A002_Xa4dcb9_X8ac.ms.split.cal',
    datacolumn = 'corrected',
    keepflags = T)

#  Flag the residual atmospheric lines in the target data
  flagdata(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split.cal', mode = 'manual', spw = '2:5193~5718', field = '3~151', flagbackup = F)

# Save flags after applycal
mystep = 19
if(mystep in thesteps):
  casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
  print 'Step ', mystep, step_title[mystep]

  
  flagmanager(vis = 'uid___A002_Xa4dcb9_X8ac.ms.split.cal',
    mode = 'save',
    versionname = 'AfterApplycal')
  
  
#es.generateQA2Report(ms1='uid___A002_Xa4dcb9_X8ac.ms',refAnt='DV04',target='5',phase_cal='2')

