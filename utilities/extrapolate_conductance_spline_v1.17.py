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
from scipy.interpolate import InterpolatedUnivariateSpline # extrapolation
from scipy.interpolate import BarycentricInterpolator      # extrapolation
from hyb_config_v1 import * # read parameters



print "Executing " + sys.argv[0] + " to analyze output files..."
interpolation_order = int(sys.argv[1])
print "The interpolation order is " + str(interpolation_order) + "."
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


#values=[[] for u in Uvalues]
#errors=[[] for u in Uvalues]

for T in Tvalues:
    sigma_file = open("giwn0_T_%.4f_%s_spline_order_%i_with_min_max_spread.dat"%(T, output_dir, interpolation_order), "w")
    for Mu_counter, Mu in enumerate(MuValues):
#    for Ucounter, U in enumerate(Uvalues):
          # prepare the input parameters; they can be used inside the script and are passed to the solver
          parms = {
                   # solver parameters
                   'SWEEPS'             : SWEEPS,
                   'DEBUGGER'           : DEBUGGER, #for debug purpose
                   'THERMALIZATION'     : THERMALIZATION,
                   'SEED'               : SEED,
                   'N_MEAS'             : N_MEAS,
                   'N_ORBITALS'         : N_ORBITALS,  
                   'BASENAME'           : "hyb.param_BETAt%.3f_Mu_%.2f_U_%.3f"%(1/T, Mu, U), # base name of the h5 output file
                   'MAX_TIME'           : MAX_TIME,
                   'VERBOSE'            : 0,
                   'VERY_VERBOSE'       : VERY_VERBOSE,
                   'TEXT_OUTPUT'        : 1,
                   'N_ENV'              : N_ENV, # number of colors
                   'SPINFLIP'           : SPINFLIP,
                   'Dissipation'        : Dissipation, # trun dissipation on or off 
                   'N_W'                : N_W,   # =1 for static susceptibility chi(0) 
                   'DELTA'              : "DELTA0.h5", # for N_ENV=1
                   'DELTA0'             : "DELTA0.h5", # color 1 for N_ENV=2                   
                   'DELTA1'             : "DELTA1.h5", # color 2 for N_ENV=2
                   'DELTA_IN_HDF5'      : 1,
                   # physical parameters
                   'U'                  : U, 
                   'MU'                 : Mu,
                   'BETA'               : 1/T, # inverse temperature 
                   'T'                  : T,   # temperature              
                   'r'                  : r,   # dissipation strength
                   'C0'                 : C0,  # dissipation capacitance
                   # measurements       
                   'MEASURE_nnw'        : MEASURE_nnw, 
                   'MEASURE_time'       : MEASURE_time,
                   'MEASURE_conductance': MEASURE_conductance,
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
    
          print "extracting result from " + output_path + "..."

          ED_file = open(output_path+"/giwn.dat", "r")
          lines = ED_file.readlines()
          ED_file.close()
          iwn=[]
          values=[[] for i in range(N_ENV if "N_ENV" in locals() else 1)]

          for line in lines:
               iwn.append(float(line.split()[0]))
               for i in range(N_ENV if "N_ENV" in locals() else 1):
                   values[i].append(float(line.split()[i+1]))
          iwn=np.array(iwn)
#          giwn=np.zeros(len(values[0]))
#          for i in range(N_ENV if "N_ENV" in locals() else 1):
#               giwn = giwn + np.array(values[i])
          sigma_file.write("%.8f"%(Mu))

          iwn = iwn[0:interpolation_order]
          # positions to inter/extrapolate
          x = np.linspace(0, max(iwn), 3*len(iwn))
          for i in range(N_ENV if "N_ENV" in locals() else 1):
             giwn = np.array(values[i][0:interpolation_order])

#             # do Barycentric interpolation
#             s = BarycentricInterpolator(iwn, giwn)
#             y = s(x)
#             sigma_file.write("   %.8f"%(y[0]))
        
             # do spline inter/extrapolation
             a=[]
             for j in range(5): # 1<=k<=5
                 s = InterpolatedUnivariateSpline(iwn, giwn, k=j+1) # spline order k: 1 linear, 2 quadratic, 3 cubic ...
                 y = s(x)
                 a.append( y[0] )
             a=np.array(a)
             sigma_file.write("   %.8f   %.8f   %.8f"%(a.mean(), a.min()-a.mean(), a.max()-a.mean() ))
#             sigma_file.write("   %.8f   %.8f"%(a.mean(), a.std(ddof=1)/sqrt(len(a))) )

          sigma_file.write("\n")
    sigma_file.write("\n")
    sigma_file.close()

          #a=plt.plot(tau, values)
          #plt.setp(a, label=r"N_ENV=1, V=%.3f"%(V*sqrt(2)) )


          # below is for chi(0)  
#          ar=archive(output_path + parms['BASENAME']+'.out.h5')
#
#          nn_0_0=ar['simulation/results/nnw_re_0_0/mean/value']
#          nn_1_1=ar['simulation/results/nnw_re_1_1/mean/value']
#          nn_1_0=ar['simulation/results/nnw_re_1_0/mean/value']
#          dnn_0_0=ar['simulation/results/nnw_re_0_0/mean/error']
#          dnn_1_1=ar['simulation/results/nnw_re_1_1/mean/error']
#          dnn_1_0=ar['simulation/results/nnw_re_1_0/mean/error']
#
#          nn  = nn_0_0 + nn_1_1 - 2*nn_1_0
#          dnn = sqrt(dnn_0_0**2 + dnn_1_1**2 + ((2*dnn_1_0)**2) )
#
#          T = parms['T']
#          values[Ucounter].append(T*nn)
#          errors[Ucounter].append(T*dnn)
#
#          del ar
    
          # The triple for loop stops here

#print "Start plotting..."
#plt.figure()
#plt.xlabel(r'$T/t$')
#plt.ylabel(r'$T\chi$')
#plt.title(r'cthyb_ohmic: $\chi$ vs T')
#plt.xscale('log')
#a=[]
#for i in range(len(Uvalues)):
#   a.append(plt.errorbar(np.array(Tvalues)/parms['t'], np.array(values[i]), errors[i], \
#               label="U=%.3f"%(Uvalues[i]/parms['t'])))
#plt.xlim(np.array([Tmin, Tmax])/parms['t'])
#plt.ylim([0.0, 1.0])
#plt.legend(loc='lower right', prop={'size':10})
#
#plt.savefig('chi_vs_T_vs_U_' + output_dir + '.pdf') 


################################# Working area #################################
#      input_file = pyalps.writeParameterFile(parms['BASENAME']+'.txt', parms)
#      subprocess.call('p2h5 '+parms['BASENAME']+'.in.h5', stdin=open(parms['BASENAME']+'.txt'), shell=True)
#a=plt.errorbar(np.array(Tvalues), np.array(values[0]), errors[0], label="U=%.1f"%MuValues[0])
#ddplt.setp(a,'r')
#b=plt.errorbar(np.array(Tvalues), np.array(values[1]), errors[1], label="U=%.1f"%MuValues[1])
#plt.setp(b,'g')
