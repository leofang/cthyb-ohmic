import numpy as np
from numpy import sqrt, exp, log, pi

############## specify the job to be done ############
goal = 1            # 1: generate input files; 2: analyze the output
DOS  = 2            # 1: constant hybridization; 2: semicircular DOS; 3: flat DOS
ZeroTresult = False  # plot the T=0 theory curve or not (for n_vs_mu_vs_T.py) 

################# physical parameters ################
W        = 1.0 # bandwidth

## on-site interaction
U        = 16    # for others 
Uvalues  = [U]
#Uvalues  = [16., 32., 50.]   # for chi_vs_T_vs_U.py

v        = sqrt(U/8.)     # symmetric coupling strength; U/\Gamma=8; semicircular
#v        = sqrt(U/(2.*pi))     # symmetric coupling strength; U/\Gamma=8; flat
#v        = sqrt(0.02)
V        = np.array([v, v])*W # coupling strength to colors 0 and 1
N_ORBITALS = 2
N_ENV = 2

# temperature (for k_vs_T_vs_mu.py)
#Tvalues = []
#for i in range(11):
#    Tvalues.append(sqrt(0.0001+i*0.00029))

# temperaturea (for n_vs_mu_vs_T.py)
#Tvalues  = W/np.array([0.01, 1.0, 100.0])[::-1] # the values here are dimensionless (=beta*t)
#Tvalues  = W/np.array([100.0]) # the values here are dimensionless (=beta*t)

# temperature (for chi_vs_T_vs_U.py)
N_T  = 25    # number of temperature points
Tmin = 0.005 # minimum temperature
Tmax = 1.0   # maximum temperature
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
THERMALIZATION = 10000
SWEEPS = 10000000000000.
MAX_TIME = 3600 # 60mins
SEED = 88
N_MEAS = 400
N_TAU = 5000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
N_MATSUBARA = int(N_TAU/(2*pi))
SPINFLIP = 0
N_HISTOGRAM_ORDERS = 150 # number of histogram orders to be measured
DEBUGGER = 0
VERY_VERBOSE = 0
MEASURE_time = 0
MEASURE_nnw = 0
MEASURE_conductance = 1
N_W=11
r=0
Dissipation=0
C0=0


################### condor parameters ################
RunOnOSG = True   # true if submitting jobs on Ci-Connect or OSGConnect
OnlyUseNanoMachines = False
executable = "cthyb-ohmic-conductance/cthyb_ohmic"  # relative to home directory
#output_dir = "cthyb_Ekin_U_%.3f"%U
#output_dir = "cthyb_semicircular_spinful_chi_V_0.180_lead_%i"%N_ENV
#output_dir = "cthyb_semicircular_test2"
#output_dir = "cthyb_semicircular_spinful_n_V_%.3f_U_%.3f_lead_%i"%(V[0], U, N_ENV)
output_dir = "20160810_semicircular_sigma_V_%.3f_U_%.3f_MAX_TIME_%i_hrs"%(sqrt(V[0]**2+V[1]**2), U, MAX_TIME/3600.)
