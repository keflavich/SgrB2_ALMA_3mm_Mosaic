import re

if re.search('^4.2.1', casadef.casa_version) == None:
 sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.2.1')


print "# Flux calibration of the data."

#                                       "J1744-3116"     0      91.28      0.57291      0.64488       2014-07-01T06:27:00       uid___A002_X85b7b2_Xb3.ms.split.cal
#                                       "J1744-3116"     0      91.28      0.66541      0.64488       2014-07-02T06:14:17     uid___A002_X85c183_X1434.ms.split.cal
#                                       "J1744-3116"     0      91.28      0.63996      0.64488       2014-07-03T04:06:55      uid___A002_X85dcf7_Xc7c.ms.split.cal
#                                       "J1744-3116"     0      91.28      0.63921      0.64488       2014-07-03T05:27:41      uid___A002_X85dcf7_Xefe.ms.split.cal
#                                       "J1744-3116"     1      89.48      0.57260      0.64662       2014-07-01T06:27:00       uid___A002_X85b7b2_Xb3.ms.split.cal
#                                       "J1744-3116"     1      89.48      0.66318      0.64662       2014-07-02T06:14:17     uid___A002_X85c183_X1434.ms.split.cal
#                                       "J1744-3116"     1      89.48      0.64136      0.64662       2014-07-03T04:06:55      uid___A002_X85dcf7_Xc7c.ms.split.cal
#                                       "J1744-3116"     1      89.48      0.64546      0.64662       2014-07-03T05:27:41      uid___A002_X85dcf7_Xefe.ms.split.cal
#                                       "J1744-3116"     2     101.37      0.55682      0.61645       2014-07-01T06:27:00       uid___A002_X85b7b2_Xb3.ms.split.cal
#                                       "J1744-3116"     2     101.37      0.65152      0.61645       2014-07-02T06:14:17     uid___A002_X85c183_X1434.ms.split.cal
#                                       "J1744-3116"     2     101.37      0.61134      0.61645       2014-07-03T04:06:55      uid___A002_X85dcf7_Xc7c.ms.split.cal
#                                       "J1744-3116"     2     101.37      0.61632      0.61645       2014-07-03T05:27:41      uid___A002_X85dcf7_Xefe.ms.split.cal
#                                       "J1744-3116"     3     103.23      0.55175      0.62013       2014-07-01T06:27:00       uid___A002_X85b7b2_Xb3.ms.split.cal
#                                       "J1744-3116"     3     103.23      0.64895      0.62013       2014-07-02T06:14:17     uid___A002_X85c183_X1434.ms.split.cal
#                                       "J1744-3116"     3     103.23      0.60971      0.62013       2014-07-03T04:06:55      uid___A002_X85dcf7_Xc7c.ms.split.cal
#                                       "J1744-3116"     3     103.23      0.60886      0.62013       2014-07-03T05:27:41      uid___A002_X85dcf7_Xefe.ms.split.cal

setjy(vis = 'uid___A002_X85b7b2_Xb3.ms.split.cal',
  field = 'J1744-3116',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.64488, 0, 0, 0])

setjy(vis = 'uid___A002_X85c183_X1434.ms.split.cal',
  field = 'J1744-3116',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.64488, 0, 0, 0])

setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split.cal',
  field = 'J1744-3116',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.64488, 0, 0, 0])

setjy(vis = 'uid___A002_X85dcf7_Xefe.ms.split.cal',
  field = 'J1744-3116',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.64488, 0, 0, 0])

setjy(vis = 'uid___A002_X85b7b2_Xb3.ms.split.cal',
  field = 'J1744-3116',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.64662, 0, 0, 0])

setjy(vis = 'uid___A002_X85c183_X1434.ms.split.cal',
  field = 'J1744-3116',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.64662, 0, 0, 0])

setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split.cal',
  field = 'J1744-3116',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.64662, 0, 0, 0])

setjy(vis = 'uid___A002_X85dcf7_Xefe.ms.split.cal',
  field = 'J1744-3116',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.64662, 0, 0, 0])

setjy(vis = 'uid___A002_X85b7b2_Xb3.ms.split.cal',
  field = 'J1744-3116',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.61645, 0, 0, 0])

setjy(vis = 'uid___A002_X85c183_X1434.ms.split.cal',
  field = 'J1744-3116',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.61645, 0, 0, 0])

setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split.cal',
  field = 'J1744-3116',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.61645, 0, 0, 0])

setjy(vis = 'uid___A002_X85dcf7_Xefe.ms.split.cal',
  field = 'J1744-3116',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.61645, 0, 0, 0])

setjy(vis = 'uid___A002_X85b7b2_Xb3.ms.split.cal',
  field = 'J1744-3116',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.62013, 0, 0, 0])

setjy(vis = 'uid___A002_X85c183_X1434.ms.split.cal',
  field = 'J1744-3116',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.62013, 0, 0, 0])

setjy(vis = 'uid___A002_X85dcf7_Xc7c.ms.split.cal',
  field = 'J1744-3116',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.62013, 0, 0, 0])

setjy(vis = 'uid___A002_X85dcf7_Xefe.ms.split.cal',
  field = 'J1744-3116',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.62013, 0, 0, 0])

os.system('rm -rf uid___A002_X85b7b2_Xb3.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X85b7b2_Xb3.ms.split.cal',
  caltable = 'uid___A002_X85b7b2_Xb3.ms.split.cal.ampli_inf',
  field = 'J1744-3116',
  solint = 'inf',
  combine = 'scan',
  refant = '',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X85b7b2_Xb3.ms.split.cal',
  field = '3,5~56', # J1744-3116,SgrB2
  gaintable = 'uid___A002_X85b7b2_Xb3.ms.split.cal.ampli_inf',
  gainfield = '3', # J1744-3116
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_X85c183_X1434.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X85c183_X1434.ms.split.cal',
  caltable = 'uid___A002_X85c183_X1434.ms.split.cal.ampli_inf',
  field = 'J1744-3116',
  solint = 'inf',
  combine = 'scan',
  refant = '',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X85c183_X1434.ms.split.cal',
  field = '3,5~56', # J1744-3116,SgrB2
  gaintable = 'uid___A002_X85c183_X1434.ms.split.cal.ampli_inf',
  gainfield = '3', # J1744-3116
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_X85dcf7_Xc7c.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X85dcf7_Xc7c.ms.split.cal',
  caltable = 'uid___A002_X85dcf7_Xc7c.ms.split.cal.ampli_inf',
  field = 'J1744-3116',
  solint = 'inf',
  combine = 'scan',
  refant = '',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X85dcf7_Xc7c.ms.split.cal',
  field = '2,4~55', # J1744-3116,SgrB2
  gaintable = 'uid___A002_X85dcf7_Xc7c.ms.split.cal.ampli_inf',
  gainfield = '2', # J1744-3116
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_X85dcf7_Xefe.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X85dcf7_Xefe.ms.split.cal',
  caltable = 'uid___A002_X85dcf7_Xefe.ms.split.cal.ampli_inf',
  field = 'J1744-3116',
  solint = 'inf',
  combine = 'scan',
  refant = '',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X85dcf7_Xefe.ms.split.cal',
  field = '2,4~55', # J1744-3116,SgrB2
  gaintable = 'uid___A002_X85dcf7_Xefe.ms.split.cal.ampli_inf',
  gainfield = '2', # J1744-3116
  calwt = F,
  flagbackup = F)

print "# Concatenating the data."

concat(vis = ['uid___A002_X85b7b2_Xb3.ms.split.cal', 'uid___A002_X85c183_X1434.ms.split.cal', 'uid___A002_X85dcf7_Xc7c.ms.split.cal', 'uid___A002_X85dcf7_Xefe.ms.split.cal'],
  concatvis = 'calibrated.ms')


