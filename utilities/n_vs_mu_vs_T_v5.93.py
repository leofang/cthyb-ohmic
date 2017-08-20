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
from os.path import expanduser #get home path
from hyb_config_v1 import * # read parameters



if goal == 1:
    print "Executing " + sys.argv[0] + " to generate input files..."
    if DOS == 1:
       print "density of states: delta function (flat hybridization)"
    elif DOS == 2:
       print "density of states: semicircular"
    elif DOS == 3:
       print "density of states: flat"
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

# always compute symmetry pt
if U/2. not in MuValues:
  MuValues.append(U/2.)
  MuValues.sort()


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
    f.write("Executable = " + expanduser("~") + "/" + executable + "\n")
    f.write("notification=Error\n")
    f.write("notify_user = leofang@phy.duke.edu\n")
    if (RunOnOSG if "RunOnOSG" in locals() else False): # run on Duke Ci-Connect
       f.write("+ProjectName=\"duke-CMT\"\n")
       f.write("should_transfer_files=YES\n") # output files are transferred
       f.write("when_to_transfer_output = ON_EXIT\n")
       # Periodically retry the jobs every 60 seconds, up to a maximum of 5 retries.
       f.write("periodic_release =  (NumJobStarts < 5) && ((CurrentTime - EnteredCurrentStatus) > 60)\n")
       # Send the job to Held state on failure. 
       f.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
       # Reconnect jobs faster if disconnected (default is 2400 sec)
       f.write("JobLeaseDuration = 90\n")
    else:  # run on physics condor
       f.write("+Department = \"Physics\"\n")
       f.write("should_transfer_files=NO\n") # use Physics shared files system
       requirement = "requirements = (( OpSys == \"LINUX\" && Arch ==\"X86_64\" && FileSystemDomain != \"\" ) && ("
       # do not use atl or phy-compute machines because they do not have OpenMPI installed
       for i in [11, 12, 13, 14, 15, 16]:
          requirement += "(TARGET.Machine != \"atl0%i.phy.duke.edu\")" %i
          requirement = (requirement+")" if i==16 else requirement+" && ")
       requirement += " && ("
       for i in [1, 4]:
          requirement += "(TARGET.Machine != \"phy-compute-0%i.phy.duke.edu\")" %i
          requirement = (requirement+")" if i==4 else requirement+" && ")
       # only use nano machines 
       if OnlyUseNanoMachines is True:
          requirement += " && ("
          for i in [6, 7, 8, 9]:
              requirement += "(TARGET.Machine == \"nano0%i.internal.phy.duke.edu\")" %i
              requirement = (requirement+")" if i==9 else requirement+" || ")
       # my code can only run on SL6 machines since the condor was back online on Sep 21, 2016 
       requirement += "&& OPSYSMAJORVER == 6" 
       f.write(requirement+")\n")
       #f.write("request_memory = 50\n")
       f.write("getenv=True\n") # important for finding shared libraries!
    #f.write("periodic_release = (NumGlobusSubmits < 5) && ((CurrentTime - EnteredCurrentStatus) > (60*60))\n")
    #f.write("periodic_hold =  (JobStatus==2) && ((CurrentTime - EnteredCurrentStatus) > (%i*60*60))" % MAX_WALL_TIME)
    f.write("\n")

## create a list of input hdf5 filenames
#numberCores = 4
#g=open("mpirun_test.sh","w")
#g.write("#!/bin/bash\n")


for T in Tvalues:
    if goal == 1: # create hybridization function
        # flat hybridization (\delta[\tau]=const)
        if DOS == 1:
            delta=[]
            for i in range(N_TAU+1):
               g0tau = -0.5
               delta.append(g0tau)
    
        # semicircular or flat density of states
        if DOS == 2 or DOS == 3: 
            g=np.zeros(N_MATSUBARA, dtype=complex) # create an array filled with 0.+0.I
            I=complex(0., 1.)
            for n in range(N_MATSUBARA):
               wn=(2*n+1)*pi*T
               if DOS == 2:
                  g[n]=2.0/(I*wn + I*sqrt(4*W**2+wn**2)) # use GF with semielliptical DOS
               if DOS == 3:
                  g[n]=log((2*W+I*wn)/(-2*W+I*wn))/(4*W) # use GF with flat DOS


            # NOTE: the beta factor is cancelled out in the exponent
            m = np.array(range(N_TAU/2+1))     # index the imaginary time (from tau=0 to tau=beta/2)
            n = np.array(range(N_MATSUBARA))   # index the Matsubara frequency
            iwn = 1j*(2.*n+1.)*pi*T            # compute iwn; note "1j" is the imaginary number i in python
            F = exp(-1j*pi*np.outer(m, 2.*n+1.)/N_TAU)  # the Fourier transform matrix
            delta = np.dot(F, g-1./iwn).real    # convert g(iwn) to g(tau), with tail subtracted
            delta *= 2.*T # the factor of 2 appears because only the real part is taken
            delta += -0.5 # add back contribution of the tail
#            delta=[]
#            for i in range(N_TAU/2+1): # generate half of the array G(0) to G(beta/2)
#               tau=i/T/N_TAU
#               g0tau=0.0;
#               for n in range(N_MATSUBARA):
#                  iwn=complex(0.0, (2*n+1)*pi*T)
#                  #g0tau+=((g[n])*exp(-iw*tau)).real # Fourier transform without tail subtracted
#                  g0tau+=((g[n]-1.0/iwn)*exp(-iwn*tau)).real # Fourier transform with tail subtracted
#               g0tau *= 2.0*T
#               g0tau += -1.0/2.0 # add back contribution of the tail
#               delta.append(g0tau) 
            temp=delta[::-1] # reverse the array to reduce the calculation time
            temp=np.delete(temp, 0)      # remove the middile point G(beta/2) to avoid double counting
            delta=np.append(delta, temp) # full array
            if len(delta) != N_TAU+1:
        	    sys.exit("Wrong array size! Abort!")
    
        for i in range(N_ENV if "N_ENV" in locals() else 1):
           ar=archive('DELTA' + str(i) + ".h5", 'w')
           for m in range(N_ORBITALS):
              ar['/Delta_%i'%m]=delta*V[i]**2
           del ar
       
        print "hybridization generated..."


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
                   'BASENAME'           : "hyb.param_BETAt%.3f_Mu_%.3f_U_%.3f"%(1./T, Mu, U), # base name of the h5 output file
                   'MAX_TIME'           : MAX_TIME,                         
                   'VERBOSE'            : 0,                            
                   'VERY_VERBOSE'       : VERY_VERBOSE,
                   'TEXT_OUTPUT'        : 1,                                
                   'N_ENV' 		: N_ENV, # number of colors
                   'SPINFLIP'           : SPINFLIP,
                   'Dissipation'        : Dissipation, # trun dissipation on or off 
                   'N_W'                : N_W,   # =1 for static susceptibility chi(0) 
                   'COLORFLIP'          : COLORFLIP,  # color flip update
                   'V0'                 : V[0],       # hybridization strength of lead 0
                   'V1'                 : V[1],       # hybridization strength of lead 1
                   'WORM'               : WORM,       # for worm update
                   'ETA'                : ETA,        # for worm update
    #               'DELTA'              : "delta-00.dat", # for N_ENV=1
    #               'DELTA0'             : "delta-00.dat",                    
    #               'DELTA1'		    : "delta-01.dat",
    #               'DELTA_IN_HDF5'      : 0,                               
                   'DELTA'              : "DELTA0.h5", # for N_ENV=1
                   'DELTA0'             : "DELTA0.h5", # color 1 for N_ENV=2                   
                   'DELTA1'		: "DELTA1.h5", # color 2 for N_ENV=2
                   'DELTA_IN_HDF5'      : 1,                               
                   # physical parameters
                   'U'                  : U,                               
                   'MU'                 : Mu,                            
                   'BETA'               : 1./T, # inverse temperature 
                   'T'			: T,   # temperature              
                   'r'                  : r,   # dissipation strength
                   'C0'                 : C0,  # dissipation capacitance
                   # measurements
                   'MEASURE_nnw'        : MEASURE_nnw,                               
                   'MEASURE_time'       : MEASURE_time,
                   'MEASURE_time_worm'  : MEASURE_time_worm,   #for worm update
                   'MEASURE_conductance': MEASURE_conductance,
                   # measurement parameters
                   'N_HISTOGRAM_ORDERS' : N_HISTOGRAM_ORDERS,           
                   'N_TAU'              : N_TAU,      
                   'N_MATSUBARA'        : N_MATSUBARA, # number of Matsubara frequencies
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
              for i in range(N_ENV if "N_ENV" in locals() else 1):
                 shutil.copyfile('DELTA' + str(i) + ".h5", output_path + 'DELTA' + str(i) + ".h5")
        
              # condor submit file
              f.write("Initialdir = " + initial_dir + "/" + output_path + "\n")
              f.write("output = " + "$(Cluster).$(Process).out\n")
              f.write("error = "  + "$(Cluster).$(Process).err\n")
              f.write("Log = "    + "$(Cluster).$(Process).log\n")
              if (RunOnOSG if "RunOnOSG" in locals() else False): # run on OSG nodes
                  f.write("transfer_input_files = " + parms['BASENAME'] + ".in.h5")
                  for i in range(N_ENV if "N_ENV" in locals() else 1):
                      f.write(', DELTA' + str(i) + ".h5")
                  f.write('\n')
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
    for i in range(N_ENV if "N_ENV" in locals() else 1):
       os.remove('DELTA' + str(i) + '.h5')

#g.close()
#os.chmod("mpirun_test.sh", 0700)


if goal == 2:
       print "Start plotting..."
       #print values[0][0]

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
       plt.ylim(np.array([0., 1.]))

       if ZeroTresult is True: 
         theory=[]
         if N_ENV==1:
             V = V[0]
         if N_ENV==2:
             V = sqrt(V[0]**2 + V[1]**2) # rescale the width
             #print V
         # Note that in this expression the bandwidth is set to be 1.
         for mu in MuValues:
            if DOS == 2:
               # my approach
               theory.append( (pi*((V**2)-1)+(V**2)*arctan(mu/sqrt(4-mu**2))+((V**2)-2)*arctan((-(V**2)+2)*mu/(V**2)/sqrt(4-mu**2)))/(2*pi*((V**2)-1)) )
               #print "mu=%.3f, n=%.10f"%(mu, (pi*((V**2)-1)+(V**2)*arctan(mu/sqrt(4-mu**2))+((V**2)-2)*arctan((-(V**2)+2)*mu/(V**2)/sqrt(4-mu**2)))/(2*pi*((V**2)-1)) )
               ## Gu's approach
               #theory.append( 0.5-arctan(mu*(V**2-2)/(sqrt(4-mu**2)*V**2))/pi  )
            if DOS == 3:
               theory.append(0.5+arctan(4*W*mu/(pi*(V**2)))/pi)
         a.append(plt.plot(np.array(MuValues)/parms['t'], np.array(theory), 'k-', linewidth=1.0, label="theory (T=0)"))
       
       plt.legend(loc='lower right', prop={'size':10})
       plt.savefig('n_vs_mu_vs_T_' + output_dir + '.pdf')
       #plt.show()
   

       print "interpolating width as a function of temperature..."
       Gamma_vs_T = "{"
       for t_counter, T in enumerate(Tvalues):
          #for Mu_counter, Mu in enumerate(MuValues):
          #t_counter = Tvalues.tolist().index(T)
          for orb in range(N_ORBITALS):
              x=0 # switch
              for Mu_counter, Mu in enumerate(MuValues):
                 if values[t_counter][orb][Mu_counter]>0.2 and x==0:
                     Mu_min = (Mu-Mu_div) + (0.2-values[t_counter][orb][Mu_counter-1])*Mu_div/(values[t_counter][orb][Mu_counter]-values[t_counter][orb][Mu_counter-1])
                     x+=1
                 #values[t_counter][orb][Mu_counter]
              x=0 # switch
              for Mu_counter, Mu in enumerate(MuValues):
                 if values[t_counter][orb][Mu_counter]>0.8 and x==0:
                     Mu_max = (Mu-Mu_div) + (0.8-values[t_counter][orb][Mu_counter-1])*Mu_div/(values[t_counter][orb][Mu_counter]-values[t_counter][orb][Mu_counter-1])
                     x+=1
              Gamma_vs_T += "{%.6f, %.3f}"%(T, Mu_max-Mu_min) 
          if t_counter != len(Tvalues)-1:
              Gamma_vs_T += ", "
       Gamma_vs_T += "}" 
       print Gamma_vs_T

################################# Working area #################################
#      input_file = pyalps.writeParameterFile(parms['BASENAME']+'.txt', parms)
#      subprocess.call('p2h5 '+parms['BASENAME']+'.in.h5', stdin=open(parms['BASENAME']+'.txt'), shell=True)
#a=plt.errorbar(np.array(Tvalues), np.array(values[0]), errors[0], label="U=%.1f"%MuValues[0])
#ddplt.setp(a,'r')
#b=plt.errorbar(np.array(Tvalues), np.array(values[1]), errors[1], label="U=%.1f"%MuValues[1])
#plt.setp(b,'g')
