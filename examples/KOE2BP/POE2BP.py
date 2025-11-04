# POE2BP.py, https://github.com/bhanson10/gbees/tree/main/examples/POE2BP
# Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

import sys
sys.path.append('../../')
import gbeespy as gbees # type: ignore

DIM_f = 2 # State dimension

# This function defines the dynamics model - required
def POE2BP(x, t, dx, coef):
    f0 = 0
    f1 = coef[0]**2.0 / x[0]**3.0; 
    return [f0, f1]

#==================================== Read in initial discrete measurement ==================================#
print("Reading in initial discrete measurement...\n")

P_DIR = "./results/python"      # Saved PDFs path
M_DIR = "./measurements"     # Measurement path
M_FILE = "measurement0.txt"  # Measurement file
M = gbees.Meas_create(DIM_f, M_DIR, M_FILE)
#============================================================================================================#

#=========================================== Read in user inputs ============================================#
print("Reading in user inputs...\n")

dx = [None] * DIM_f                            # Grid width, default is half of the std. dev. from the initial measurement 
dx[0] = M.cov[0][0]**0.5/8;                    # dL
dx[1] = M.cov[1][1]**0.5/5;                    # dl
G = gbees.Grid_create(DIM_f, 1E-9, M.mean, dx) # Inputs: (dimension, probability threshold, center, grid width)    

coef = [19.910350621818949]                    # POE2BP trajectory attributes (mu)
T = gbees.Traj_create(len(coef), coef)         # Inputs: (# of coefficients, coefficients)

NUM_DIST = 5                                   # Number of distributions recorded per measurement
NUM_MEAS = 1                                   # Number of measurements
DEL_STEP = 20                                  # Number of steps per deletion procedure
OUTPUT_FREQ = 20                               # Number of steps per output to terminal
CAPACITY = int(2**13);                         # Size of hash table (power of 2 for optimal hashing)
OUTPUT = False                                 # Write info to terminal
RECORD = True                                  # Write PDFs to .txt file
MEASURE = False                                # Take discrete measurement updates
BOUNDS = False                                 # Add inadmissible regions to grid
COLLISIONS = False;                            # Track collisions
TV = False;                                    # Time-invariant dynamics 
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(POE2BP, None, None, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, CAPACITY, DIM_f, OUTPUT, RECORD, MEASURE, BOUNDS, COLLISIONS, TV)