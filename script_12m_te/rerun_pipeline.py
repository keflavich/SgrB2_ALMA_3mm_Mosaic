# run from script directory
# (or just cut out the chdir and run from raw...)
# This script is intended to re-run the pipeline so that you can use the latest
# version of the CASA pipeline on your data.
#
# I ran into a problem where the 4.2 pipeline produced files with a .tar.gz
# suffix but the new pipeline wants a .tgz suffix.  Not sure if this fixes it,
# but I figured it can't hurt to try...
# (it's also not clear whether the pipeline script itself needs to be
# regenerated for the newest version of CASA, but I can't find docs on how to
# do that)
#
# Based on directions in the ALMA pipeline quickstart guide:
# https://almascience.nrao.edu/documents-and-tools/cycle-2/alma-pipeline-quickstart-guide/at_download/file
import os
import shutil
import glob

os.chdir('../raw')

shutil.copy('../script/casa_pipescript.py', '.')
shutil.copy('../calibration/flux.csv', '.')
for fn in glob.glob("../calibration/*flagtemplate.txt"):
    shutil.copy(fn, '.')

for fn in glob.glob("*.asdm.sdm"):
    os.symlink(fn, fn[:-9])

execfile('casa_pipescript.py')
