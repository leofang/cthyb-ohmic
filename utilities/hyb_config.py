#Tvalues=[0.01, 0.1, 1.0, 10.0, 100.0]
Tvalues=[0.01, 0.02]
N_Mu     = 20
Mu_min   = -5.
Mu_max   = 5.
V = [1., 1.];   # coupling strength to colors 0 and 1

MAX_TIME = 600
N_TAU = 5000    # number of tau-points; must be large enough for the lowest temperature (set to at least 5*BETA*U)
N_ORBITALS = 1
#runtime = 1     # solver runtime (in seconds)
N_HISTOGRAM_ORDERS = 250 # number of histogram orders to be measured

executable = "/home/yf30/alps/bin/hybridization"
output_dir = "output_flat"

