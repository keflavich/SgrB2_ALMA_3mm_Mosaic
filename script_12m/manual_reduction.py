"""
This is my guess at the manual reduction process from September 1. I did this
all successfully once in the past (June or July?), but I do not presently have
access to logs from that (and I can't remember what machine I did it on...)
EDIT:
    based on log casapy-20150617-094231.log

1. Some machines have this error and need their libxmls fixed before importasdm
   can be run:
https://help.almascience.org/index.php?/Knowledgebase/Article/View/266/5/what-should-i-do-if-scriptforpi-aborts-at-the-importasdm-step-with-error
2. Just run scriptForPI... apparently it worked.
3. Nope, that's wrong!  I definitely modified the scriptForCalibration's somehow.  I guess I produced them all from the QA2 generator?
   AHA! generateReducScript!

ipython-20150616-163341.log:es.es.generateReducScript('name_of_asdm')
ipython-20150616-163341.log:es.generateReducScript('name_of_asdm')
ipython-20150616-163341.log:es.generateReducScript('uid___A002_X9d13e3_Xd4f')
ipython-20150616-171839.log:es.generateReducScript('uid___A002_X95e355_X220a')
ipython-20150616-171843.log:es.generateReducScript('uid___A002_X95e355_X1f13')
ipython-20150616-171843.log:es.generateReducScript(['uid___A002_X95e355_X1f13.ms.split.cal', 'uid___A002_X95e355_X220a.ms.split.cal', 'uid___A002_X9cffbd_Xefe.ms.split.cal', 'uid___A002_X9d13e3_Xd4f.ms.split.cal'], step='fluxcal')
corresponds to 171818
ipython-20150616-171843.log:es.generateReducScript('calibrated.ms',step='imaging')
ipython-20150616-171848.log:es.generateReducScript('uid___A002_X9cffbd_Xefe')

es.generateReducScript('uid___A002_X95e355_X1f13')
mysteps=[]
execfile('uid___A002_X95e355_X1f13.ms.scriptForCalibration.py')
es.generateReducScript(['uid___A002_X95e355_X1f13.ms.split.cal',
                        'uid___A002_X95e355_X220a.ms.split.cal',
                        'uid___A002_X9cffbd_Xefe.ms.split.cal',
                        'uid___A002_X9d13e3_Xd4f.ms.split.cal'], step='fluxcal')
execfile("scriptForFluxCalibration.py") # make sure this is the es-generated one!!
es.generateReducScript('calibrated.ms',step='imaging')
execfile("scriptForImaging.py")
"""
# Start this in the Script directory
import os
import shutil
import glob
cwd = os.path.split(os.getcwd())[1]
if cwd != 'script':
    raise ValueError("Start this script in the script/ directory")

# AnalysisUtils required (probably because I didn't disable the QA2 stuff...)
import sys
sys.path.append(os.path.join(os.getcwd(), 'analysis_scripts'))
sys.path.append("/users/thunter/AIV/science/analysis_scripts/")
import analysisUtils as au
import analysisUtils as aU
es = au.stuffForScienceDataReduction()

# yep, start it in script, then immediately chdir
os.chdir('../calibrated')

sourcefiles = [
               '../raw/uid___A002_X95e355_X1f13.asdm.sdm',
               '../raw/uid___A002_X9cffbd_Xefe.asdm.sdm',
               '../raw/uid___A002_X95e355_X220a.asdm.sdm',
               '../raw/uid___A002_X9d13e3_Xd4f.asdm.sdm',
]



uidnames = [os.path.basename(sf[:-9]) for sf in sourcefiles]

for source in sourcefiles:
    for target in (',', os.path.basename(source[:-9])):
        try:
            os.symlink(source, target)
        except OSError:
            continue

for uidname in uidnames:
    es.generateReducScript(uidname)

for fn in glob.glob("uid*scriptForCalibration.py"):
    shutil.copy(fn, '../script/')

calnames = [unm+".ms.split.cal" for unm in uidnames]

# Normally shouldn't need to do this?  But, if there are '-'s in the flux
# column of allFluxes.txt (which gets used by es.generateReducScript), it
# is needed.
es.generateFluxFile(calnames)

es.generateReducScript(calnames,
                       step='fluxcal')

shutil.copy('scriptForFluxCalibration.py', '../script/')


for uidname in uidnames:
    importasdm(uidname, uidname+".ms",
               asis='Antenna Station Receiver Source CalAtmosphere CalWVR',
               bdfflags=True, lazy=False)

for uidname in uidnames:
    execfile("../script/{0}.ms.scriptForCalibration.py".format(uidname))

execfile("../script/scriptForFluxCalibration.py") # make sure this is the es-generated one!!
es.generateReducScript('calibrated.ms',step='imaging')
#execfile("../script/scriptForImaging.py")
