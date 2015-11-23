 #############################################################################/
 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations
 #
 # ALPS Libraries
 #
 # Copyright (C) 2012 by Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 #
 #
 # This software is part of the ALPS Applications, published under the ALPS
 # Application License; you can use, redistribute it and/or modify it under
 # the terms of the license, either version 1 or (at your option) any later
 # version.
 # 
 # You should have received a copy of the ALPS Application License along with
 # the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 # available from http://alps.comp-phys.org/.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 # FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 # SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 # FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 # ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 # DEALINGS IN THE SOFTWARE.
 #
 #############################################################################/

 # This tutorial is a minimal example illustrating the use of the python 
 # interface to hybridization expansion solver.
 #
 # The hybridization function is chosen such that the impurity model is 
 # equivalent to a correlated site couple to a single bath site with coupling 
 # V=1 and energy epsilon=0 (Delta(tau)=-V**2/2=const.) . It can therefore be 
 # compared to the exact result obtained by exact diagonalization (see 
 # subdirectory ED).
 #
 # Run this script as:
 # alpspython tutorial1.py
 #
 # For multi-core machines the mpi can be used by, for example, setting MPI=2 
 # in "runApplication".
 #
 # Modified by Leo Fang (leofang@phy.duke.edu)
 #
 # Reason: The script was incompatible with the ALPS binary package v2.2.0b3 
 # for Mac OS X, so it is rewritten to fit the style of other ALPS Python 
 # tutorial scripts. In particular it relies heavily on the pyalps modules.
 #
 # Last updated: May, 11, 2015

# Import the pyalps and plotting libraries
import pyalps
import matplotlib.pyplot as plt
import pyalps.plot
from pyalps.hdf5 import archive # hdf5 interface
import numpy as np

# Specify solver parameters as a dictionary
parms={
'SWEEPS'              : 5000,
'MAX_TIME'            : 60,
'THERMALIZATION'      : 100,
'SEED'                : 7,
'N_MEAS'              : 50,
'N_HISTOGRAM_ORDERS'  : 50,
'N_ORBITALS'          : 2,
'U'                   : 4.0,
'MU'                  : 2.0,
'DELTA'               : "delta-01.dat",
'N_TAU'               : 1000,
'BETA'                : 45,
'TEXT_OUTPUT'         : 1,
'DEBUGGER'            : 50000, # Leo Fang: this number is solely for debug purpose
}

# Write the parameters to a HDF5 input file
input_file = pyalps.writeInputH5Files('test_input',[parms])

# Write a simple (constant) hybridization function to file "delta-01.dat".
# This corresponds to the impurity site with on-site interaction U
# coupled to a single noninteracting bath site with hybridization
# V=1 at energy epsilon=0: Delta(tau)=-V^2/2=const.
f=open("delta-01.dat","w")
for i in range(parms["N_TAU"]+1):
	f.write("%i %f %f\n"%(i,-0.5,-0.5))
f.close()

## Feed the HDF5 file into the CT-HYB application to solve the impurity model;
## the parameter "MPI=2" can be modified or removed
##pyalps.runApplication('hybridization',input_file,MPI=2)
#
## Load the output HDF5 file
#result_files = pyalps.getResultFiles(prefix='parms-hyb-01')
#
## Extract the measurement results. We first look at the histograms.
#data = pyalps.loadMeasurements(result_files,['order_histogram_0','order_histogram_1','order_histogram_total'])
#data = pyalps.flatten(data)
#
## Make one plot with all data
#for dataset in data:
#    dataset.props['label'] = str(dataset.props['observable'])
#
## Generate histograms
#plt.figure()
#plt.xlabel('Order')
#plt.title('Order histogram')
#pyalps.plot.plot(data)
#plt.legend()
#
## We next look at the imaginary time Green functions. 
## Load the exact diagonalization result
#ED_file = open("./exact_diagonalization/Gt.dat","r")
#lines = ED_file.readlines()
#ED_file.close()
#del lines[0] # Remove the header
#tau=[]
#up=[]
#down=[]
#for line in lines:
#	tau.append(float(line.split()[0]))
#	up.append(float(line.split()[1]))
#	down.append(float(line.split()[2]))
#
#
## Load the imaginary time Green functions calculated by CT-HYB 
#data = pyalps.loadMeasurements(result_files,['g_0','g_1'])
#data = pyalps.flatten(data)
#
## Generate labels and rescale the imaginary time
#for dataset in data:
#    dataset.props['label'] = str(dataset.props['observable'])
#    # Rescale the imaginary time according to BETA
#    dataset.x = dataset.x*dataset.props['BETA']/float(dataset.props['N_TAU'])
#
#
## Plot Green functions of CT-HYB vs ED 
#plt.figure()
#plt.xlabel(r'$\tau$')
#plt.ylabel(r'$G_i(\tau)$')
#plt.title('Imaginary-time Green functions: CT-HYB vs ED ')
## CT-HYB results
#pyalps.plot.plot(data)
## ED results
#a=plt.plot(tau,up,'r--') # red dashed curve
#b=plt.plot(tau,down,'g:') # green dotted curve
#plt.setp(a,label='G_up_ED')
#plt.setp(b,label='G_down_ED')
#
#
## The end points of the CT-HYB Green functions require post processing (one 
## can zoom in the above figure and see that they deviate from the ED results; 
## see the documentation). Because loadMeasurements returns MCVectorData objects
## that cannot be further manipulated, we read the data directly from the output 
## HDF5 file: 
# 
#ar=archive('parms-hyb-01.task1.out.h5','r')
#g_0_mean=ar['simulation/results/g_0/mean/value']
#g_0_error=ar['simulation/results/g_0/mean/error']
#g_1_mean=ar['simulation/results/g_1/mean/value']
#g_1_error=ar['simulation/results/g_1/mean/error']
#for i in [-1,0]:
#    g_0_mean[i]=g_0_mean[i]*2.0
#    g_1_mean[i]=g_1_mean[i]*2.0
#tau_CTHYB=np.array(range(parms['N_TAU']+1))*float(parms['BETA'])/float(parms['N_TAU'])
#
## Plot corrected Green functions of CT-HYB vs ED 
#plt.figure()
#plt.xlabel(r'$\tau$')
#plt.ylabel(r'$G_i(\tau)$')
#plt.title('Imaginary-time Green functions: CT-HYB vs ED (corrected)')
## CT-HYB results
#G=[]
#G.append(plt.errorbar(tau_CTHYB, np.array(g_0_mean), g_0_error, label="G_%i"%i))
#G.append(plt.errorbar(tau_CTHYB, np.array(g_1_mean), g_1_error, label="G_%i"%i))
## ED results
#a=plt.plot(tau,up,'r--') # red dashed curve
#b=plt.plot(tau,down,'g:') # green dotted curve
#plt.setp(a,label='G_up_ED')
#plt.setp(b,label='G_down_ED')
#
#
### Show all figures
#plt.legend()
#plt.show()



################################# Working area #################################
#import subprocess
# export the parameters to the text file "parms-hyb-01"
#input_file = pyalps.writeParameterFile('parms-hyb-01', parms)
# convert the generated text file into an HDF5 file
#subprocess.call('p2h5 parms-hyb-01.h5', stdin=open('parms-hyb-01'), shell=True)
#subprocess.call("p2h5 parms-hyb-01.h5 < parms-hyb-01", shell=True)
#result_files = pyalps.getResultFiles(pattern='*.out.h5')
#result_files = pyalps.getResultFiles(prefix='parms-hyb-01.out')
#print pyalps.loadObservableList(result_files)
#print data2[0].x[3]
#print data2[0].props['BETA']
#print type(data).__name__
