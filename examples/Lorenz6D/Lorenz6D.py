# Lorenz6D.py, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz6D
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

import sys
sys.path.append('../../')
import gbeespy as gbees # type: ignore

DIM_f = 6 # State dimension

# This function defines the dynamics model - required
def Lorenz6D(x, t, dx, coef):
    f0 = (x[1] - x[4]) * x[5] - x[0] + coef[0]
    f1 = (x[2] - x[5]) * x[0] - x[1] + coef[0]
    f2 = (x[3] - x[0]) * x[1] - x[2] + coef[0]
    f3 = (x[4] - x[1]) * x[2] - x[3] + coef[0]
    f4 = (x[5] - x[2]) * x[3] - x[4] + coef[0]
    f5 = (x[0] - x[3]) * x[4] - x[5] + coef[0]
    return [f0, f1, f2, f3, f4, f5]

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
for i in range(DIM_f):
    dx[i] = (M.cov[i][i]**(0.5))/2
G = gbees.Grid_create(DIM_f, 0.0, 8E-9, M.mean, dx) # Inputs: (dimension, initial time, probability threshold, center, grid width)    

coef = [4.0]                                   # Lorenz6D trajectory attributes (F)
T = gbees.Traj_create(len(coef), coef)         # Inputs: (# of coefficients, coefficients)

NUM_DIST = 2                                   # Number of distributions recorded per measurement
NUM_MEAS = 1                                   # Number of measurements
DEL_STEP = 20                                  # Number of steps per deletion procedure
OUTPUT_FREQ = 20                               # Number of steps per output to terminal
CAPACITY = int(2**26);                         # Size of hash table (power of 2 for optimal hashing)
OUTPUT = True                                  # Write info to terminal
RECORD = True                                  # Write PDFs to .txt file
MEASURE = False                                # Take discrete measurement updates
BOUNDS = False                                 # Add inadmissible regions to grid
COLLISIONS = False;                            # Track collisions
TV = False;                                    # Time-invariant dynamics 
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(Lorenz6D, None, None, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, CAPACITY, DIM_f, OUTPUT, RECORD, MEASURE, BOUNDS, COLLISIONS, TV)