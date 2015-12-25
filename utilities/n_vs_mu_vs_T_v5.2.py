import pyalps
from pyalps.hdf5 import archive   # hdf5 interface
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyalps.plot
import numpy as np
from numpy import exp, log, sqrt, pi, arctan # some math
import subprocess
import os, sys, shutil
from hyb_config import * # read parameters



if goal == 1:
    print "Executing " + sys.argv[0] + " to generate input files..."
    if DOS == 1:
       print "density of states: delta function (flat hybridization)"
    elif DOS == 2:
       print "density of states: semicircular"
    else:
       sys.exit("Invalid DOS! Abort!")
elif goal == 2:
    print "Executing " + sys.argv[0] + " to analyze output files..."
#elif goal == 0:
#    print "Executing " + sys.argv[0] + " for test purpose..."
#    print Tvalues
#    sys.exit()
else:
    sys.exit("Invalid goal! Abort!")


initial_dir = os.getcwd()


if 'MuValues' not in locals():
    MuValues = []
    Mu_div   = (Mu_max-Mu_min)/N_Mu
    for i in range(N_Mu+1):
      Mu = Mu_min + Mu_div*i
      MuValues.append(Mu)


values=[[[0 for Mu in MuValues] for orb in range(N_ORBITALS)] for T in Tvalues] 
errors=[[[0 for Mu in MuValues] for orb in range(N_ORBITALS)] for T in Tvalues] 

#histogram_values=[]
#histogram_errors=[]
#x=0 # counter for the histogram

# create condor log folder and prepare condor submit file
if goal == 1:
    if os.path.exists(output_dir): # if the directory exists, remove and then create
       shutil.rmtree(output_dir)
       os.makedirs(output_dir)
    else:     			   # if it doesn't exist, simply create
       os.makedirs(output_dir)
    f=open("condor_submit_file","w")
    
    #f.write("Universe = parallel\n") # run MPI jobs ---- seems not working!
    #f.write("machine_count = 1\n") # number of cores for a MPI job
    f.write("Universe = vanilla\n")
    f.write("Executable = " + executable + "\n")
    f.write("notification=Error\n")
    f.write("notify_user = leofang@phy.duke.edu\n")
    f.write("+Department = \"Physics\"\n")
    f.write("should_transfer_files=NO\n") # use Physics shared files system
    #f.write("should_transfer_files=YES\n") # test if output text files are transferred
    #f.write("when_to_transfer_output = ON_EXIT\n")
    requirement = "requirements = (( OpSys == \"LINUX\" && Arch ==\"X86_64\" && FileSystemDomain != \"\" ) && ("
    # do not use atl machines because they do not have OpenMPI installed
    for i in [11, 12, 13, 14, 15, 16]:
       requirement += "(TARGET.Machine != \"atl0%i.phy.duke.edu\")" %i
       requirement = (requirement+")" if i==16 else requirement+" && ")
    # only use nano machines 
    if OnlyUseNanoMachines is True:
       requirement += " && ("
       for i in [6, 7, 8, 9]:
           requirement += "(TARGET.Machine == \"nano0%i.internal.phy.duke.edu\")" %i
           requirement = (requirement+")" if i==9 else requirement+" || ")
    f.write(requirement+")\n")
    f.write("request_memory = 50\n")
    f.write("getenv=True\n") # important for finding shared libraries!
    #f.write("periodic_release = (NumGlobusSubmits < 5) && ((CurrentTime - EnteredCurrentStatus) > (60*60))\n")
    #f.write("periodic_hold =  (JobStatus==2) && ((CurrentTime - EnteredCurrentStatus) > (%i*60*60))" % MAX_WALL_TIME)
    #f.write("+ProjectName=\"duke\"\n")
    f.write("\n")


## create a list of input hdf5 filenames
#numberCores = 4
#g=open("mpirun_test.sh","w")
#g.write("#!/bin/bash\n")


for T in Tvalues:
    if goal == 1: # create hybridization function
        # flat density of states
        if DOS == 1:
            delta=[]
            for i in range(N_TAU+1):
               g0tau = -0.5
               delta.append(g0tau)
    
        # semicircular density of states
        if DOS == 2:  
            g=[]
            I=complex(0., 1.)
            for n in range(N_MATSUBARA):
               wn=(2*n+1)*pi*T
               g.append(2.0/(I*wn + I*sqrt(4*W**2+wn**2))) # use GF with semielliptical DOS
            delta=[]
            for i in range(N_TAU/2+1): # generate half of the array G(0) to G(beta/2)
               tau=i/T/N_TAU
               g0tau=0.0;
               for n in range(N_MATSUBARA):
                  iwn=complex(0.0, (2*n+1)*pi*T)
                  #g0tau+=((g[n])*exp(-iw*tau)).real # Fourier transform without tail subtracted
                  g0tau+=((g[n]-1.0/iwn)*exp(-iwn*tau)).real # Fourier transform with tail subtracted
               g0tau *= 2.0*T
               g0tau += -1.0/2.0 # add back contribution of the tail
               delta.append(g0tau) 
            temp=delta[::-1] # reverse the array to reduce the calculation time
            del temp[0]      # remove the middile point G(beta/2) to avoid double counting
            delta=delta+temp # full array
            if len(delta) != N_TAU+1:
        	    sys.exit("Wrong array size! Abort!")
    
        for i in range(N_ENV if N_ENV in locals() else 1):
           ar=archive('DELTA' + str(i) + ".h5", 'w')
           for m in range(N_ORBITALS):
              ar['/Delta_%i'%m]=np.array(delta)*V[i]**2 # need to conver delta to a numpy array!
           del ar
       
        print "done"


    for Mu_counter, Mu in enumerate(MuValues):
          # prepare the input parameters; they can be used inside the script and are passed to the solver
          parms = {
                   # solver parameters
                   'SWEEPS'             : SWEEPS,                       
                   'DEBUGGER'           : DEBUGGER, #for debug purpose
                   'THERMALIZATION'     : THERMALIZATION,  
                   'SEED'               : SEED,                               
                   'N_MEAS'             : N_MEAS,      
                   'N_ORBITALS'         : N_ORBITALS,  
                   'BASENAME'           : "hyb.param_BETAt%.3f_Mu_%.2f_U_%.3f"%(W/T, Mu, U), # base name of the h5 output file
                   'MAX_TIME'           : MAX_TIME,                         
                   'VERBOSE'            : 0,                            
                   'VERY_VERBOSE'       : VERY_VERBOSE,
                   'TEXT_OUTPUT'        : 1,                                
                   'N_ENV' 		: N_ENV, # number of colors
                   'SPINFLIP'           : SPINFLIP,
    #               'DELTA'              : "delta-00.dat", # for N_ENV=1
    #               'DELTA0'             : "delta-00.dat",                    
    #               'DELTA1'		    : "delta-01.dat",
    #               'DELTA_IN_HDF5'      : 0,                               
                   'DELTA'              : "DELTA0.h5", # for N_ENV=1
                   'DELTA0'             : "DELTA0.h5",                    
                   'DELTA1'		    : "DELTA1.h5",
                   'DELTA_IN_HDF5'      : 1,                               
                   # physical parameters
                   'U'                  : U,                               
                   'MU'                 : Mu,                            
                   'BETA'               : 1/T, # inverse temperature 
                   'T'			: T,   # temperature              
                   # measurements
                   'MEASURE_nnw'        : MEASURE_nnw,                               
		   'N_W'                : 1,   # for static susceptibility chi(0) 
                   'MEASURE_time'       : MEASURE_time,
                   # measurement parameters
                   'N_HISTOGRAM_ORDERS' : N_HISTOGRAM_ORDERS,           
                   'N_TAU'              : N_TAU,      
                   'N_MATSUBARA'        : int(N_TAU/(2*pi)), # number of Matsubara frequencies
                   #'N_W'                   : 1, # number of bosonic Matsubara frequencies for the local susceptibility
                   # additional parameters (used outside the solver only)
                   't'                  : W,   # hopping strength/bandwidth
                 }
              
          # Write the parameters to a HDF5 input & output files as well as the condor submit file
          # and put the hybridization functions in the same directory (in order to avoid file transfer)
          output_path = output_dir + "/" + parms['BASENAME'] + "/"
    
          if goal == 1:
              print output_path
              os.makedirs(output_path)
        
              # input & output files
              input_file = pyalps.writeInputH5Files(output_path + parms['BASENAME']+'.in.h5', [parms])
              #shutil.copyfile(output_path + parms['BASENAME']+'.in.h5', output_path + parms['BASENAME']+'.out.h5')
         
              # copy hybridization function to the path
              for i in range(N_ENV if N_ENV in locals() else 1):
                 shutil.copyfile('DELTA' + str(i) + ".h5", output_path + 'DELTA' + str(i) + ".h5")
        
              # condor submit file
              f.write("Initialdir = " + initial_dir + "/" + output_path + "\n")
              f.write("output = " + "$(Cluster).$(Process).out\n")
              f.write("error = "  + "$(Cluster).$(Process).err\n")
              f.write("Log = "    + "$(Cluster).$(Process).log\n")
        #      f.write("transfer_input_files = " + parms['BASENAME'] + ".in.h5," + parms['BASENAME'] + ".out.h5\n")
        #      f.write("transfer_input_files = " + parms['BASENAME'] + ".in.h5," + parms['BASENAME'] + ".out.h5\n")
              f.write("Arguments = " + parms['BASENAME']+ ".in.h5" + "\n")
              f.write("Queue\n\n")
    
    #      # Write mpirun commands
    #      g.write("mpirun -np " + str(numberCores) + " ../cthyb_ohmic " + parms['BASENAME']+ ".in.h5\n")
    #      g.write("mv simulation.dat simulation_" + parms['BASENAME'] + ".dat\n")
    #  
    #      # Solve the impurity model in parallel
    #      pyalps.runApplication('../cthyb_ohmic',parms['BASENAME']+'.in.h5')
      
          if goal == 2:   
              ar=archive(output_path + parms['BASENAME']+'.out.h5')

              t_counter = Tvalues.tolist().index(parms['T'])
              for orb in range(N_ORBITALS):
                 n_mean=ar['simulation/results/density_%i/mean/value'%orb]
                 n_error=ar['simulation/results/density_%i/mean/error'%orb]
                 values[t_counter][orb][Mu_counter] = n_mean
                 errors[t_counter][orb][Mu_counter] = n_error
    
              del ar
    
          # The triple for loop stops here

if goal == 1:
    f.close()
    for i in range(N_ENV if N_ENV in locals() else 1):
       os.remove('DELTA' + str(i) + '.h5')

#g.close()
#os.chmod("mpirun_test.sh", 0700)


if goal == 2:
       print "Start plotting..."
       plt.figure()
       plt.xlabel(r'$\mu/t$')
       plt.ylabel(r'$n$')
       plt.title(r'cthyb_ohmic: density vs $\mu$')
       a=[]
       for i in range(len(Tvalues)):
          for orb in range(N_ORBITALS):
             a.append(plt.errorbar(np.array(MuValues)/parms['t'], np.array(values[i][orb]), errors[i][orb], \
                      label=r"orbital %i, $\beta t$=%.3f"%(orb, parms['t']/Tvalues[i])))
       plt.xlim(np.array([Mu_min, Mu_max])/parms['t'])
 
       theory=[]
       for mu in MuValues:
          # Note that in this expression the bandwidth is set to be 1.
          theory.append((pi*(V[0]**2-1)+V[0]**2*arctan(mu/sqrt(4-mu**2))+(V[0]**2-2)*arctan((-V[0]**2+2)*mu/V[0]**2/sqrt(4-mu**2)))/(2*pi*(V[0]**2-1)))
       a.append(plt.plot(np.array(MuValues)/parms['t'], np.array(theory), 'k-', linewidth=2.0, label="theory (T=0)"))
       
       plt.legend(loc='lower right', prop={'size':10})
       plt.savefig('n_vs_mu_vs_T_' + output_dir + '.pdf')
       #plt.show()
   


################################# Working area #################################
#      input_file = pyalps.writeParameterFile(parms['BASENAME']+'.txt', parms)
#      subprocess.call('p2h5 '+parms['BASENAME']+'.in.h5', stdin=open(parms['BASENAME']+'.txt'), shell=True)
#a=plt.errorbar(np.array(Tvalues), np.array(values[0]), errors[0], label="U=%.1f"%MuValues[0])
#ddplt.setp(a,'r')
#b=plt.errorbar(np.array(Tvalues), np.array(values[1]), errors[1], label="U=%.1f"%MuValues[1])
#plt.setp(b,'g')
