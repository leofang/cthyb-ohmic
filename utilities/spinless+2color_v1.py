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
 #
 # This tutorial shows how to repeatedly solve a single impurity Anderson model
 # using the hybridization expansion solver and how to extract the calculated
 # observables and their error via the python interface.
 #
 # The script can be used to reproduce the results of Figs. 1 and 2 of the ALPS
 # CT-HYB paper.
 #
 # The results show the decrease of the effective local moment of the impurity
 # with decreasing temperature due to Kondo screening. For simplicity, a 
 # hybridization function corresponding to a semielliptical density of states 
 # is used.
 #
 # The impurity model is solved by calling the hybridization executable. For 
 # details on how to use the solver see the documentation.
 #
 # Run this script as:
 # alpspython tutorial2.py
 #
 # Modified by Leo Fang (leofang@phy.duke.edu)
 #
 # Reason: The script was incompatible with the ALPS binary package v2.2.0b3 
 # for Mac OS X, so some of which is rewritten to fit the style of other ALPS 
 # Python tutorial scripts.
 #
 # Last updated: May, 11, 2015

import pyalps
from pyalps.hdf5 import archive   # hdf5 interface
import pyalps.cthyb
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np
from numpy import exp, log, sqrt, pi, arctan # some math
import subprocess
import os


# generate a sequence of temperatures between a and b which are equidistant on a logarithmic scale
#N_T  = 5    # number of temperature points
#Tmin = 0.01 # minimum temperature
#Tmax = 100.0 # maximum temperature
#Tdiv = exp(log(Tmax/Tmin)/N_T)
#T=Tmax
#Tvalues=[]
#for i in range(N_T+1):
#  Tvalues.append(T)
#  T/=Tdiv

#Tvalues=[0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
Tvalues=[0.001]

N_Mu     = 20
Mu_min   = -5.
Mu_max   = 5.
Mu_div   = (Mu_max-Mu_min)/N_Mu
MuValues = []
for i in range(N_Mu+1):
  Mu = Mu_min + Mu_div*i
  MuValues.append(Mu)

#print Tvalues;
#print MuValues;

V = [1., 1.];   # coupling strength to colors 0 and 1
N_TAU = 1000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
N_ORBITALS = 1
#runtime = 1     # solver runtime (in seconds)
N_HISTOGRAM_ORDERS = 200 # number of histogram orders to be measured


#values=[[] for Mu in MuValues]
#errors=[[] for Mu in MuValues]
values=[[0 for Mu in MuValues] for T in Tvalues] 
errors=[[0 for Mu in MuValues] for T in Tvalues] 
histogram_values=[]
histogram_errors=[]
x=0 # counter for the histogram


#for i in [0, 1]:
#     # Create the hybridization function according to the given parameters
#     print "creating initial hybridization..." 
#     delta=[]
#     for j in range(N_TAU+1):
#       g0tau = -0.5*V[i]**2
#       delta.append(g0tau) # delta=t**2 g
#   
#     # Write hybridization function to hdf5 archive (solver input)
#     ar=archive('DELTA'+str(i)+".h5",'w')
#     for m in range(N_ORBITALS):
#         ar['/Delta_%i'%m]=delta
#     del ar


#f=open("condor_submit_file","w")
#f.write("Universe = parallel\n") # run MPI jobs
#f.write("Executable = ../cthyb_ohmic\n")
#f.write("Initialdir = /home/yf30/cthyb-ohmic-color/test_spinless+2color\n")
#f.write("output = condor-output/$(Cluster).$(Process).out\n")
#f.write("error = condor-output/$(Cluster).$(Process).err\n")
#f.write("Log = condor-output/$(Cluster).$(Process).log\n")
#f.write("notification=Error\n")
#f.write("notify_user = leofang@phy.duke.edu\n")
#f.write("+Department = \"Physics\"\n")
#f.write("should_transfer_files=NO\n") # use Physics shared files system
##f.write("when_to_transfer_output = ON_EXIT\n")
#f.write("requirements = ( OpSys == \"LINUX\" && Arch ==\"X86_64\" && FileSystemDomain != \"\" )\n")
#f.write("machine_count = 1") # number of cores for a MPI job
#f.write("request_memory = 50")
##f.write("periodic_release = (NumGlobusSubmits < 5) && ((CurrentTime - EnteredCurrentStatus) > (60*60))\n")
##f.write("periodic_hold =  (JobStatus==2) && ((CurrentTime - EnteredCurrentStatus) > (%i*60*60))" % MAX_WALL_TIME)
##f.write("+ProjectName=\"duke\"\n")
#f.write("\n")


for Mu_counter, Mu in enumerate(MuValues):
  parameters=[]
  for t in Tvalues:
      # prepare the input parameters; they can be used inside the script and are passed to the solver
      parameters.append(
             {
               # solver parameters
               'SWEEPS'                     : 1000000000,                       
               'DEBUGGER'                   : 50000, #for debug purpose
               'THERMALIZATION'             : 1000,  
               'SEED'                       : 43,                               
               'N_MEAS'                     : 10,      
               'N_ORBITALS'                 : 1,  
               'BASENAME'                   : "hyb.param_Mu%.1f_BETA%.3f"%(Mu,1/t), # base name of the h5 output file
               'MAX_TIME'                   : 120,                          
               'VERBOSE'                    : 0,                            
               'TEXT_OUTPUT'                : 1,                                
               'N_ENV' 			    : 2,
               'SPINFLIP'                   : 0,
               'DELTA0'                     : "DELTA0.h5",                    
               'DELTA1'			    : "DELTA1.h5",
               'DELTA_IN_HDF5'              : 1,                               
               # physical parameters
               'U'                          : 4,                               
               'MU'                         : Mu,                            
               'BETA'                       : 1/t,                   
               # measurements
               'MEASURE_nnw'                : 0,                               
               #'MEASURE_time'               : 1,
               # measurement parameters
               'N_HISTOGRAM_ORDERS'         : N_HISTOGRAM_ORDERS,           
               'N_TAU'                      : N_TAU,      
               'N_MATSUBARA'                : int(N_TAU/(2*pi)), # number of Matsubara frequencies
               #'N_W'                        : 1, # number of bosonic Matsubara frequencies for the local susceptibility
               # additional parameters (used outside the solver only)
               't'                          : t, 
             }
          )


  for parms in parameters:
#        # Write the parameters into the output HDF5 file
#        ar=archive(parms['BASENAME']+'.out.h5','a')
#        ar['/parameters']=parms
#        del ar
#     
#        # Write the parameters to a HDF5 input file
#        input_file = pyalps.writeInputH5Files(parms['BASENAME']+'.in.h5',[parms])
#
#        # Write the parameters to the condor submit file
#        f.write("Arguments = " + parms['BASENAME']+ ".in.h5" + "\n")
#        f.write("Queue\n\n")
  
        # Solve the impurity model in parallel
        pyalps.runApplication('../cthyb_ohmic', parms['BASENAME']+'.in.h5')

        # rename the summary file
        for filename in os.listdir("."):
            if filename.startswith("simulation"):
                  os.rename(filename, "simulation_" + parms['BASENAME'] + ".dat") 
  
#        # extract the local occupation
#        ar=archive(parms['BASENAME']+'.out.h5','w')
#        n_mean=ar['simulation/results/density_0/mean/value']
#        n_error=ar['simulation/results/density_0/mean/error']
#  
#        del ar
#  
#        #T=1.0/parms['BETA']
#        #t_counter = Tvalues.index(1.0/parms['BETA'])
#        #print t_counter
#        t_counter = Tvalues.index(parms['t'])
#        values[t_counter][Mu_counter] = n_mean
#        errors[t_counter][Mu_counter] = n_error
#  
#        # The for loop stops here

#f.close()

#
#
#
## Generate the plot of susceptibility vs temperature
#plt.figure()
#plt.xlabel(r'$\mu$')
#plt.ylabel(r'$n$')
#plt.title(r'cthyb_ohmic: density vs $\mu$')
#a=[]
#for i in range(len(Tvalues)):
#    a.append(plt.errorbar(np.array(MuValues), np.array(values[i]), errors[i], label="T=%.2f"%Tvalues[i]))
#plt.xlim([Mu_min, Mu_max])
##plt.ylim([0.4, 0.6])
#plt.legend(loc='lower right')
#
#theory = (0.5+arctan(MuValues)/pi)
#a.append(plt.plot(np.array(MuValues), theory, 'k-', linewidth=2.0, label="theory (T=0)"))
#
#plt.show()


################################# Working area #################################
#      input_file = pyalps.writeParameterFile(parms['BASENAME']+'.txt', parms)
#      subprocess.call('p2h5 '+parms['BASENAME']+'.in.h5', stdin=open(parms['BASENAME']+'.txt'), shell=True)
#a=plt.errorbar(np.array(Tvalues), np.array(values[0]), errors[0], label="U=%.1f"%MuValues[0])
#ddplt.setp(a,'r')
#b=plt.errorbar(np.array(Tvalues), np.array(values[1]), errors[1], label="U=%.1f"%MuValues[1])
#plt.setp(b,'g')
