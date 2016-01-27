import numpy as np
from numpy import sqrt, exp, log, pi

############## specify the job to be done ############
goal = 2            # 1: generate input files; 2: analyze the output
DOS  = 2            # 1: constant hybridization; 2: semicircular DOS; 3: flat DOS
ZeroTresult = False  # plot the T=0 theory curve or not (for n_vs_mu_vs_T.py) 

################# physical parameters ################
W        = 1.0 # bandwidth
V        = np.array([0.18, 0.18])*W # coupling strength to colors 0 and 1
N_ORBITALS = 2
N_ENV = 2

## on-site interaction (for chi_vs_T_vs_U.py)
Uvalues  = [0., 0.105, 0.21, 0.42, 0.84, 2.5]

# on-site interaction (for others)
#U = 0.105

# temperature (for k_vs_T_vs_mu.py)
#Tvalues = []
#for i in range(11):
#    Tvalues.append(sqrt(0.0001+i*0.00029))

# temperaturea (for n_vs_mu_vs_T.py)
#Tvalues  = W/np.array([0.01, 0.1, 1.0, 10.0, 100.0])[::-1] # the values here are dimensionless (=beta*t)
#Tvalues  = W/np.array([100.0]) # the values here are dimensionless (=beta*t)

# temperature (for chi_vs_T_vs_U.py)
N_T  = 25    # number of temperature points
Tmin = 0.002 # minimum temperature
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

# chemical potential (for k_vs_T_vs_mu.py)
#MuValues = np.array([0.0, 0.5])*W

# chemical potential (for test purpose)
#MuValues = [0.12]


############## cthyb solver parameters ###############
THERMALIZATION = 2000
SWEEPS = 1000000000
MAX_TIME = 600
SEED = 88
N_MEAS = 250
N_TAU = 3000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
N_MATSUBARA = int(N_TAU/(2*pi))
SPINFLIP = 0
N_HISTOGRAM_ORDERS = 150 # number of histogram orders to be measured
DEBUGGER = 10000
VERY_VERBOSE = 0
MEASURE_time = 0
MEASURE_nnw = 1


################### condor parameters ################
OnlyUseNanoMachines = False
executable = "/home/yf30/cthyb-ohmic-color/cthyb_ohmic"
#output_dir = "cthyb_Ekin_U_%.3f"%U
output_dir = "cthyb_semicircular_spinful_chi_V_0.18_lead_%i"%N_ENV
#output_dir = "cthyb_semicircular_test2"
#output_dir = "cthyb_semicircular_spinful_n_V_0.18_U_%.3f_lead_%i"%(U, N_ENV)

