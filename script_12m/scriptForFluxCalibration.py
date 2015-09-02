import re

if re.search('^4.3.1', casadef.casa_version) == None:
 sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.3.1')


print "# Flux calibration of the data."

#                                       "J1744-3116"     0      91.28            -     5.55588       2014-12-06T14:36:04     uid___A002_X95e355_X1f13.ms.split.cal
#                                       "J1744-3116"     0      91.28      5.55588      5.55588       2014-12-06T16:04:17     uid___A002_X95e355_X220a.ms.split.cal
#                                       "J1744-3116"     0      91.28            -     5.55588       2015-04-01T07:28:11      uid___A002_X9cffbd_Xefe.ms.split.cal
#                                       "J1744-3116"     0      91.28            -     5.55588       2015-04-02T08:41:09      uid___A002_X9d13e3_Xd4f.ms.split.cal
#                                       "J1744-3116"     1      89.48            -     5.46806       2014-12-06T14:36:04     uid___A002_X95e355_X1f13.ms.split.cal
#                                       "J1744-3116"     1      89.48      5.46806      5.46806       2014-12-06T16:04:17     uid___A002_X95e355_X220a.ms.split.cal
#                                       "J1744-3116"     1      89.48            -     5.46806       2015-04-01T07:28:11      uid___A002_X9cffbd_Xefe.ms.split.cal
#                                       "J1744-3116"     1      89.48            -     5.46806       2015-04-02T08:41:09      uid___A002_X9d13e3_Xd4f.ms.split.cal
#                                       "J1744-3116"     2     101.37      5.45036      5.45036       2014-12-06T14:36:04     uid___A002_X95e355_X1f13.ms.split.cal
#                                       "J1744-3116"     2     101.37            -     5.45036       2014-12-06T16:04:17     uid___A002_X95e355_X220a.ms.split.cal
#                                       "J1744-3116"     2     101.37            -     5.45036       2015-04-01T07:28:11      uid___A002_X9cffbd_Xefe.ms.split.cal
#                                       "J1744-3116"     2     101.37            -     5.45036       2015-04-02T08:41:09      uid___A002_X9d13e3_Xd4f.ms.split.cal
#                                       "J1744-3116"     3     103.23      0.56742      0.56742       2014-12-06T14:36:04     uid___A002_X95e355_X1f13.ms.split.cal
#                                       "J1744-3116"     3     103.23            -     0.56742       2014-12-06T16:04:17     uid___A002_X95e355_X220a.ms.split.cal
#                                       "J1744-3116"     3     103.23            -     0.56742       2015-04-01T07:28:11      uid___A002_X9cffbd_Xefe.ms.split.cal
#                                       "J1744-3116"     3     103.23            -     0.56742       2015-04-02T08:41:09      uid___A002_X9d13e3_Xd4f.ms.split.cal
#                                       "J1752-2956"     0      91.28            -     0.05556       2014-12-06T14:36:04     uid___A002_X95e355_X1f13.ms.split.cal
#                                       "J1752-2956"     0      91.28            -     0.05556       2014-12-06T16:04:17     uid___A002_X95e355_X220a.ms.split.cal
#                                       "J1752-2956"     0      91.28      0.05556      0.05556       2015-04-01T07:28:11      uid___A002_X9cffbd_Xefe.ms.split.cal
#                                       "J1752-2956"     0      91.28      0.06456      0.05556       2015-04-02T08:41:09      uid___A002_X9d13e3_Xd4f.ms.split.cal
#                                       "J1752-2956"     1      89.48            -     0.05595       2014-12-06T14:36:04     uid___A002_X95e355_X1f13.ms.split.cal
#                                       "J1752-2956"     1      89.48            -     0.05595       2014-12-06T16:04:17     uid___A002_X95e355_X220a.ms.split.cal
#                                       "J1752-2956"     1      89.48      0.05571      0.05595       2015-04-01T07:28:11      uid___A002_X9cffbd_Xefe.ms.split.cal
#                                       "J1752-2956"     1      89.48      0.06434      0.05595       2015-04-02T08:41:09      uid___A002_X9d13e3_Xd4f.ms.split.cal
#                                       "J1752-2956"     2     101.37            -     0.05754       2014-12-06T14:36:04     uid___A002_X95e355_X1f13.ms.split.cal
#                                       "J1752-2956"     2     101.37            -     0.05754       2014-12-06T16:04:17     uid___A002_X95e355_X220a.ms.split.cal
#                                       "J1752-2956"     2     101.37      0.05078      0.05754       2015-04-01T07:28:11      uid___A002_X9cffbd_Xefe.ms.split.cal
#                                       "J1752-2956"     2     101.37      0.05988      0.05754       2015-04-02T08:41:09      uid___A002_X9d13e3_Xd4f.ms.split.cal
#                                       "J1752-2956"     3     103.23            -     0.05090       2014-12-06T14:36:04     uid___A002_X95e355_X1f13.ms.split.cal
#                                       "J1752-2956"     3     103.23            -     0.05090       2014-12-06T16:04:17     uid___A002_X95e355_X220a.ms.split.cal
#                                       "J1752-2956"     3     103.23      0.05049      0.05090       2015-04-01T07:28:11      uid___A002_X9cffbd_Xefe.ms.split.cal
#                                       "J1752-2956"     3     103.23      0.06001      0.05090       2015-04-02T08:41:09      uid___A002_X9d13e3_Xd4f.ms.split.cal

setjy(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = 'J1744-3116',
  spw = '0',
  standard = 'manual',
  fluxdensity = [5.55588, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = 'J1744-3116',
  spw = '0',
  standard = 'manual',
  fluxdensity = [5.55588, 0, 0, 0])

setjy(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = 'J1744-3116',
  spw = '0',
  standard = 'manual',
  fluxdensity = [5.55588, 0, 0, 0])

setjy(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = 'J1744-3116',
  spw = '0',
  standard = 'manual',
  fluxdensity = [5.55588, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = 'J1744-3116',
  spw = '1',
  standard = 'manual',
  fluxdensity = [5.46806, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = 'J1744-3116',
  spw = '1',
  standard = 'manual',
  fluxdensity = [5.46806, 0, 0, 0])

setjy(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = 'J1744-3116',
  spw = '1',
  standard = 'manual',
  fluxdensity = [5.46806, 0, 0, 0])

setjy(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = 'J1744-3116',
  spw = '1',
  standard = 'manual',
  fluxdensity = [5.46806, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = 'J1744-3116',
  spw = '2',
  standard = 'manual',
  fluxdensity = [5.45036, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = 'J1744-3116',
  spw = '2',
  standard = 'manual',
  fluxdensity = [5.45036, 0, 0, 0])

setjy(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = 'J1744-3116',
  spw = '2',
  standard = 'manual',
  fluxdensity = [5.45036, 0, 0, 0])

setjy(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = 'J1744-3116',
  spw = '2',
  standard = 'manual',
  fluxdensity = [5.45036, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = 'J1744-3116',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.56742, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = 'J1744-3116',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.56742, 0, 0, 0])

setjy(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = 'J1744-3116',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.56742, 0, 0, 0])

setjy(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = 'J1744-3116',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.56742, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = 'J1752-2956',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.05556, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = 'J1752-2956',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.05556, 0, 0, 0])

setjy(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = 'J1752-2956',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.05556, 0, 0, 0])

setjy(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = 'J1752-2956',
  spw = '0',
  standard = 'manual',
  fluxdensity = [0.05556, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = 'J1752-2956',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.05595, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = 'J1752-2956',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.05595, 0, 0, 0])

setjy(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = 'J1752-2956',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.05595, 0, 0, 0])

setjy(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = 'J1752-2956',
  spw = '1',
  standard = 'manual',
  fluxdensity = [0.05595, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = 'J1752-2956',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.05754, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = 'J1752-2956',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.05754, 0, 0, 0])

setjy(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = 'J1752-2956',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.05754, 0, 0, 0])

setjy(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = 'J1752-2956',
  spw = '2',
  standard = 'manual',
  fluxdensity = [0.05754, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = 'J1752-2956',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.05090, 0, 0, 0])

setjy(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = 'J1752-2956',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.05090, 0, 0, 0])

setjy(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = 'J1752-2956',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.05090, 0, 0, 0])

setjy(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = 'J1752-2956',
  spw = '3',
  standard = 'manual',
  fluxdensity = [0.05090, 0, 0, 0])

os.system('rm -rf uid___A002_X95e355_X1f13.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  caltable = 'uid___A002_X95e355_X1f13.ms.split.cal.ampli_inf',
  field = 'J1744-3116,J1752-2956',
  solint = 'inf',
  combine = 'scan',
  refant = 'DV18',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X95e355_X1f13.ms.split.cal',
  field = '2,3~151', # J1744-3116,SgrB2
  gaintable = 'uid___A002_X95e355_X1f13.ms.split.cal.ampli_inf',
  gainfield = '2', # J1744-3116
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_X95e355_X220a.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  caltable = 'uid___A002_X95e355_X220a.ms.split.cal.ampli_inf',
  field = 'J1744-3116,J1752-2956',
  solint = 'inf',
  combine = 'scan',
  refant = 'DV18',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X95e355_X220a.ms.split.cal',
  field = '2,3~151', # J1744-3116,SgrB2
  gaintable = 'uid___A002_X95e355_X220a.ms.split.cal.ampli_inf',
  gainfield = '2', # J1744-3116
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_X9cffbd_Xefe.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  caltable = 'uid___A002_X9cffbd_Xefe.ms.split.cal.ampli_inf',
  field = 'J1744-3116,J1752-2956',
  solint = 'inf',
  combine = 'scan',
  refant = 'DA64',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X9cffbd_Xefe.ms.split.cal',
  field = '2,3~151', # J1752-2956,SgrB2
  gaintable = 'uid___A002_X9cffbd_Xefe.ms.split.cal.ampli_inf',
  gainfield = '2', # J1752-2956
  calwt = F,
  flagbackup = F)

os.system('rm -rf uid___A002_X9d13e3_Xd4f.ms.split.cal.ampli_inf') 
gaincal(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  caltable = 'uid___A002_X9d13e3_Xd4f.ms.split.cal.ampli_inf',
  field = 'J1744-3116,J1752-2956',
  solint = 'inf',
  combine = 'scan',
  refant = 'DV18',
  gaintype = 'T',
  calmode = 'a')

applycal(vis = 'uid___A002_X9d13e3_Xd4f.ms.split.cal',
  field = '2,3~151', # J1752-2956,SgrB2
  gaintable = 'uid___A002_X9d13e3_Xd4f.ms.split.cal.ampli_inf',
  gainfield = '2', # J1752-2956
  calwt = F,
  flagbackup = F)

print "# Concatenating the data."

concat(vis = ['uid___A002_X95e355_X1f13.ms.split.cal', 'uid___A002_X95e355_X220a.ms.split.cal', 'uid___A002_X9cffbd_Xefe.ms.split.cal', 'uid___A002_X9d13e3_Xd4f.ms.split.cal'],
  concatvis = 'calibrated.ms')


