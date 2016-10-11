"""
Based on TC manual_reduction: do the same for 7m

"""
# Start this in the Script directory
import os
cwd = os.path.split(os.getcwd())[1]
if cwd != 'script':
    raise ValueError("Start this script in the script/ directory")

# AnalysisUtils required (probably because I didn't disable the QA2 stuff...)
import sys
sys.path.append(os.path.join(os.getcwd(), 'analysis_scripts'))
sys.path.append("/users/thunter/AIV/science/analysis_scripts/")
import analysisUtils as au
es = au.stuffForScienceDataReduction()

# yep, start it in script, then immediately chdir
os.chdir('../calibrated')

sourcefiles = [
    '../raw/uid___A002_X85b7b2_Xb3.asdm.sdm',
    '../raw/uid___A002_X85c183_X1434.asdm.sdm',
    '../raw/uid___A002_X85dcf7_Xc7c.asdm.sdm',
    '../raw/uid___A002_X85dcf7_Xefe.asdm.sdm',
    ]

uidnames = [os.path.basename(sf[:-9]) for sf in sourcefiles]

for source in sourcefiles:
    for target in (',', source[:-9]):
        try:
            os.symlink(source, target)
        except OSError:
            continue

for uidname in uidnames:
    es.generateReducScript(uidname)

calnames = [unm+".ms.split.cal" for unm in uidnames]

es.generateReducScript(calnames,
                       step='fluxcal')
es.generateReducScript('calibrated.ms',step='imaging')


for uidname in uidnames:
    importasdm(uidname, uidname+".ms",
               asis='Antenna Station Receiver Source CalAtmosphere CalWVR',
               bdfflags=True, lazy=False)

for uidname in uidnames:
    execfile("../script/{0}.ms.scriptForCalibration.py".format(uidname))

execfile("../script/scriptForFluxCalibration.py") # make sure this is the es-generated one!!
#execfile("../script/scriptForImaging.py")
