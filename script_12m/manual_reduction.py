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
es.generateReducScript('uid___A002_X95e355_X1f13')
mysteps=[]
execfile('uid___A002_X95e355_X1f13.ms.scriptForCalibration.py')
es.generateReducScript(['uid___A002_X95e355_X1f13.ms.split.cal',
                        'uid___A002_X95e355_X220a.ms.split.cal',
                        'uid___A002_X9cffbd_Xefe.ms.split.cal',
                        'uid___A002_X9d13e3_Xd4f.ms.split.cal'], step='fluxcal')
execfile("scriptForFluxCalibration.py")
es.generateReducScript('calibrated.ms',step='imaging')
execfile("scriptForImaging.py")
"""
# Start this in the Script directory
import os
cwd = os.path.split(os.path.getcwd())[1]
if cwd != 'script':
    raise ValueError("Start this script in the script/ directory")

# AnalysisUtils required (probably because I didn't disable the QA2 stuff...)
import sys
sys.path.append(os.path.join(os.getcwd(), 'analysis_scripts'))
import analysisUtils as aU

# yep, start it in script, then immediately chdir
os.chdir('../calibrated')

os.symlink('../raw/uid___A002_X95e355_X1f13.asdm.sdm/', '.')
os.symlink('../raw/uid___A002_X9cffbd_Xefe.asdm.sdm/', '.')
os.symlink('../raw/uid___A002_X95e355_X220a.asdm.sdm/', '.')
os.symlink('../raw/uid___A002_X9d13e3_Xd4f.asdm.sdm/', '.')

importasdm('uid___A002_X95e355_X1f13.asdm.sdm', 'uid___A002_X95e355_X1f13.ms', asis='Antenna Station Receiver Source CalAtmosphere CalWVR', bdfflags=True, lazy=False)
importasdm('uid___A002_X9cffbd_Xefe.asdm.sdm',  'uid___A002_X9cffbd_Xefe.ms', asis='Antenna Station Receiver Source CalAtmosphere CalWVR', bdfflags=True, lazy=False)
importasdm('uid___A002_X95e355_X220a.asdm.sdm', 'uid___A002_X95e355_X220a.ms', asis='Antenna Station Receiver Source CalAtmosphere CalWVR', bdfflags=True, lazy=False)
importasdm('uid___A002_X9d13e3_Xd4f.asdm.sdm',  'uid___A002_X9d13e3_Xd4f.ms', asis='Antenna Station Receiver Source CalAtmosphere CalWVR', bdfflags=True, lazy=False)

execfile('../script/uid___A002_X95e355_X1f13.ms.scriptForCalibration.py')
execfile('../script/uid___A002_X9cffbd_Xefe.ms.scriptForCalibration.py')
execfile('../script/uid___A002_X95e355_X220a.ms.scriptForCalibration.py')
execfile('../script/uid___A002_X9d13e3_Xd4f.ms.scriptForCalibration.py')
