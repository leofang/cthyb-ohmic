# /bin/bash
#
# To access the pyalps module provided in VisTrails, one has to put symbolic links of
# its dylib files in a certain directory (/usr/local/lib/alps_dylib/ in this case) and
# specify the environment variable DYLD_LIBRARY_PATH so that python can find them at
# runtime.
#
# The notebook IPython+pyalps.ipynb provides an example to import pyalps in IPython 
# and then use it to perform the analysis of the tutorial hybridization-01, so put both
# files in the hybridization-01 directory.
#
# Note that although IPython is also provided in VisTrails, it is so outdated that can
# not be used (no notebokk interface). Thus, using user-installed IPython is necessary.
#
# Last updated: Nov. 30, 2015


export DYLD_LIBRARY_PATH=/usr/local/lib/alps_dylib/
ipython notebook


################## working area ##################
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/usr/local/lib/alps_dylib/:/Applications/VisTrails/VisTrails.app/Contents/Resources/lib/
#export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH:/Applications/VisTrails/VisTrails.app/Contents/Resources/lib/
##export DYLD_FRAMEWORK_PATH=$DYLD_FRAMEWORK_PATH:/usr/local/lib/alps_dylib/:/Applications/VisTrails/VisTrails.app/Contents/Frameworks
#export DYLD_FALLBACK_FRAMEWORK_PATH=$DYLD_FALLBACK_FRAMEWORK_PATH:/Applications/VisTrails/VisTrails.app/Contents/Frameworks
#export PYTHONPATH=$PYTHONPATH:/Applications/VisTrails/VisTrails.app/Contents/Resources/:/Applications/VisTrails/VisTrails.app/Contents/Resources/lib/:/Applications/VisTrails/VisTrails.app/Contents/Resources/lib/python2.7


