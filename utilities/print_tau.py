import pyalps
from pyalps.hdf5 import archive
import sys
import numpy as np

#print sys.argv[1]
#ar=archive(sys.argv[1])
TIME = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 200, 240, 300]) * 180
Mu = 0.5
#print TIME

for time in TIME:
#    print time
   ar=archive("output_semicircular_MAX_TIME_"+str(time)+"/hyb.param_BETA100.000_Mu_"+str(Mu)+"/hyb.param_BETA100.000_Mu_"+str(Mu)+".out.h5")
   print ar["/simulation/results/density_0/tau/value"]
