import numpy as np
from numpy import sqrt, exp, log, pi

############## specify the job to be done ############
goal = 1            # 1: generate input files; 2: analyze the output
DOS  = 2            # 1: constant hybridization; 2: semicircular DOS; 3: flat DOS
ZeroTresult = False  # plot the T=0 theory curve or not (for n_vs_mu_vs_T.py) 

################# physical parameters ################
W        = 1.0 # bandwidth
V        = np.array([0.18, 0.18])*sqrt(2.)*W # coupling strength to colors 0 and 1
N_ORBITALS = 2
N_ENV = 1

## on-site interaction (for chi_vs_T_vs_U.py)
#Uvalues  = [0., 0.105, 0.21, 0.42, 0.84, 2.5]

# on-site interaction (for others)
U = 0.

# temperature (for k_vs_T_vs_mu.py)
#Tvalues = []
#for i in range(11):
#    Tvalues.append(sqrt(0.0001+i*0.00029))

# temperaturea (for n_vs_mu_vs_T.py)
#Tvalues  = W/np.array([0.01, 0.1, 1.0, 10.0, 100.0])[::-1] # the values here are dimensionless (=beta*t)
Tvalues  = W/np.array([100., 45.0]) # the values here are dimensionless (=beta*t)

# temperature (for chi_vs_T_vs_U.py)
#N_T  = 25    # number of temperature points
#Tmin = 0.002 # minimum temperature
#Tmax = 100.0 # maximum temperature
#Tdiv = exp(log(Tmax/Tmin)/N_T)
#T=Tmax
#Tvalues=[]
#for i in range(N_T+1):
#  Tvalues.append(T)
#  T/=Tdiv

## chemical potential (for n_vs_mu_vs_T.py)
#N_Mu     = 25
#Mu_min   = U/2-0.2*W
#Mu_max   = U/2+0.2*W

# chemical potential (for k_vs_T_vs_mu.py)
#MuValues = np.array([0.0, 0.5])*W

# chemical potential (for test purpose)
MuValues = [-0.1, 0., 0.1]


############## cthyb solver parameters ###############
THERMALIZATION = 10000
SWEEPS = 10000000000000.
MAX_TIME = 60
SEED = 88
N_MEAS = 400
N_TAU = 3000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
N_MATSUBARA = int(N_TAU/(2*pi))
SPINFLIP = 0
N_HISTOGRAM_ORDERS = 150 # number of histogram orders to be measured
DEBUGGER = 10000
VERY_VERBOSE = 0
MEASURE_time = 1
MEASURE_nnw = 0
MEASURE_conductance = 1
Dissipation = 0
r=0
C0=0
N_W = 10

################### condor parameters ################
OnlyUseNanoMachines = False
executable = "/home/yf30/cthyb-ohmic-ohmic/cthyb_ohmic"
#output_dir = "cthyb_Ekin_U_%.3f"%U
#output_dir = "cthyb_semicircular_spinful_chi_V_0.180_lead_%i"%N_ENV
#output_dir = "cthyb_semicircular_test2"
#output_dir = "cthyb_semicircular_spinless_n_V_%.3f_U_%.3f_lead_%i_long"%(V[0], U, N_ENV)
output_dir = "20160805_TEST_lead_%i_orb_%i_V_%.3f_MAX_TIME_%i_mins"%(N_ENV, N_ORBITALS, V[0], MAX_TIME/60.)
