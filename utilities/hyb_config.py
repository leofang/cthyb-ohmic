import numpy as np
from numpy import sqrt, exp, log

############## specify the job to be done ############
goal = 1 # 1: generate input files; 2: analyze the output
DOS  = 2 # 1: constant hybridization; 2: semicircular DOS


################# physical parameters ################
W        = 1.0 # bandwidth
V        = np.array([0.255, 0.255])*W # coupling strength to colors 0 and 1
N_ORBITALS = 2

# on-site interaction (for chi_vs_T_vs_U.py)
Uvalues  = [0., 2.5]

# on-site interaction (for others)
#U = 0

# temperature (for k_vs_T_vs_mu.py)
#Tvalues = []
#for i in range(11):
#    Tvalues.append(sqrt(0.0001+i*0.00029))

# temperaturea (for n_vs_mu_vs_T.py)
#Tvalues  = W/np.array([0.01, 0.1, 1.0, 10.0, 100.0])[::-1] # the values here are dimensionless (=beta*t)
#Tvalues  = W/np.array([100.0]) # the values here are dimensionless (=beta*t)

# temperature (for chi_vs_T_vs_U.py)
N_T  = 20    # number of temperature points
Tmin = 0.006 # minimum temperature
Tmax = 100.0 # maximum temperature
Tdiv = exp(log(Tmax/Tmin)/N_T)
T=Tmax
Tvalues=[]
for i in range(N_T+1):
  Tvalues.append(T)
  T/=Tdiv

# chemical potential (for n_vs_mu_vs_T.py)
#N_Mu     = 25
#Mu_min   = U/2-1.0*W
#Mu_max   = U/2+1.0*W

# chemical potential (for test)
#MuValues = np.array([0.0, 0.5])*W


############## cthyb solver parameters ###############
THERMALIZATION = 2000
SWEEPS = 1000000000
MAX_TIME = 600
SEED = 88
N_MEAS = 250
N_TAU = 3000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
N_MATSUBARA = 480
SPINFLIP = 1
N_HISTOGRAM_ORDERS = 150 # number of histogram orders to be measured
DEBUGGER = 10000
VERY_VERBOSE = 0
N_ENV = 1
MEASURE_time = 0
MEASURE_nnw = 1


################### condor parameters ################
OnlyUseNanoMachines = False
executable = "/home/yf30/cthyb-ohmic-color/cthyb_ohmic"
#executable = "/home/yf30/alps/bin/hybridization"
output_dir = "output_chi2"


