import re

if re.search('^4.2.2', casadef.casa_version) == None: 
    sys.exit('ERROR: PLEASE USE THE SAME VERSION OF CASA THAT YOU USED FOR GENERATING THE SCRIPT: 4.2.2')

concat(vis=['uid___A002_X9d13e3_Xd4f.ms.split.cal','uid___A002_X9cffbd_Xefe.ms.split.cal','uid___A002_X95e355_X1f13.ms.split.cal'],
       concatvis='SgrB2_a_03_TC.calibrated.ms', copypointing = False)
