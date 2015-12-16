# specify the job to be done
goal = 1 # 1: generate input files; 2: analyze the output
DOS  = 1 # 1: constant hybridization; 2: semicircular DOS

# physical parameters
#Tvalues=[0.01, 0.1, 1.0, 10.0, 100.0]
Tvalues  = [0.01, 0.1]
MuValues = [0.5]
N_Mu     = 1
Mu_min   = 0.5
Mu_max   = 0.5
V        = [1., 1.]   # coupling strength to colors 0 and 1

# cthyb solver parameters
THERMALIZATION = 1000
SWEEPS = 1000000000
MAX_TIME = 1200
SEED = 88
N_MEAS = 100
N_TAU = 10000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
N_ORBITALS = 1
N_HISTOGRAM_ORDERS = 100 # number of histogram orders to be measured
DEBUGGER = 10000
VERY_VERBOSE = 0
N_ENV = 1

# condor parameters
executable = "../cthyb_ohmic"
#executable = "/home/yf30/alps/bin/hybridization"
output_dir = "output_flat"


