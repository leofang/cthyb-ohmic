import numpy as np
from numpy import sqrt, exp, log, pi

############## specify the job to be done ############
goal = 1            # 1: generate input files; 2: analyze the output
DOS  = 3            # 1: constant hybridization; 2: semicircular DOS; 3: flat DOS
ZeroTresult = False  # plot the T=0 theory curve or not (for n_vs_mu_vs_T.py) 

################# physical parameters ################
W        = 1.0 # bandwidth

## on-site interaction
U        = 0    # for others 
#Uvalues  = [U]
#Uvalues  = [16., 32., 50.]   # for chi_vs_T_vs_U.py

#v        = sqrt(U/16.)     # symmetric coupling strength; U/\Gamma=8; semicircular; two leads
#v        = sqrt(U/(2.*pi))     # symmetric coupling strength; U/\Gamma=8; flat
V0        = sqrt(0.25)
V1        = sqrt(0.25)
V        = np.array([V0, V1])*W # coupling strength to colors 0 and 1
N_ORBITALS = 1
N_ENV = 2

# temperature (for k_vs_T_vs_mu.py)
#Tvalues = []
#for i in range(11):
#    Tvalues.append(sqrt(0.0001+i*0.00029))

# temperaturea (for n_vs_mu_vs_T.py)
Tvalues  = W/np.array([100.])[::-1] # the values here are dimensionless (=beta*t)
#Tvalues  = W/np.array([20., 40., 60., 80., 100.0, 120., 140., 160., 180., 190., 200.])[::-1] # the values here are dimensionless (=beta*t)
#Tvalues  = W/np.array([160., 180., 190., 200.]) # the values here are dimensionless (=beta*t)

# temperature (for chi_vs_T_vs_U.py)
#N_T  = 25    # number of temperature points
#Tmin = 0.005 # minimum temperature
#Tmax = 1.0   # maximum temperature
#Tdiv = exp(log(Tmax/Tmin)/N_T)
#T=Tmax
#Tvalues=[]
#for i in range(N_T+1):
#  Tvalues.append(T)
#  T/=Tdiv

# chemical potential (for n_vs_mu_vs_T.py)
N_Mu     = 1
if U > sqrt(V[0]**2+V[1]**2):
   Mu_min = -0.5*U
   Mu_max = 1.5*U
else:
   Mu_min = U/2-1.2*W
   Mu_max = U/2+1.2*W

# chemical potential (for k_vs_T_vs_mu.py)
#MuValues = np.array([0.0, 0.5])*W

# chemical potential (for test purpose)
#MuValues = [0.12]


############## cthyb solver parameters ###############
THERMALIZATION = 10000
SWEEPS = 10000000000000.
MAX_TIME = 60 # 20mins
SEED = 88
N_MEAS = 400
N_TAU = 3000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
N_MATSUBARA = int(N_TAU/(2*pi))
SPINFLIP = 0
N_HISTOGRAM_ORDERS = 150 # number of histogram orders to be measured
DEBUGGER = 50000
VERY_VERBOSE = 0
MEASURE_time = 0
MEASURE_nnw = 0

MEASURE_conductance = 1
N_W=10

r=0
Dissipation=0
C0=0

COLORFLIP=1

WORM = 1
ETA = 0.2
MEASURE_time_worm = 0

################### condor parameters ################
RunOnOSG = True
OnlyUseNanoMachines = False
executable = "cthyb-ohmic-conductance/cthyb_ohmic"
#output_dir = "cthyb_Ekin_U_%.3f"%U
#output_dir = "cthyb_semicircular_spinful_chi_V_0.180_lead_%i"%N_ENV
#output_dir = "cthyb_semicircular_test2"
#output_dir = "cthyb_semicircular_spinful_n_V_%.3f_U_%.3f_lead_%i"%(V[0], U, N_ENV)
#output_dir = "20170127_semicircular_sigma_%i_lead_%i_orb_V_%.3f_U_%.3f_MAX_TIME_%.1f_hrs"%(N_ENV, N_ORBITALS, (V[0] if N_ENV is 1 else sqrt(V[0]**2+V[1]**2)), U, MAX_TIME/3600.)
output_dir = "20170212_worm"
