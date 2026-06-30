# duffing.py, https://github.com/bhanson10/gbees/tree/main/examples/duffing
# Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

import math
import sys
sys.path.append('../../')
import gbeespy as gbees # type: ignore

DIM_f = 2 # State dimension

# This function defines the dynamics model - required
def duffing(x, t, coef):

    f0 = x[1];
    f1 = x[0] - x[0] * x[0] * x[0] - coef[0] * x[1] + coef[1] * math.cos(coef[2] * t); 
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

factor = [2.0, 2.0]
G = gbees.Grid_create(DIM_f, 0.0, 1E-9, M, factor, False) # Inputs: (dimension, initial time, probability threshold, measure, grid width factor, rotate grid)    

coef = [0.35, 0.3, 1]                                     # duffing trajectory attributes (mu)
T = gbees.Traj_create(len(coef), coef)                    # Inputs: (# of coefficients, coefficients)

NUM_DIST = 11                                             # Number of distributions recorded per measurement
NUM_MEAS = 1                                              # Number of measurements
DEL_STEP = 20                                             # Number of steps per deletion procedure
OUTPUT_FREQ = 20                                          # Number of steps per output to terminal
CAPACITY = int(2**13);                                    # Size of hash table (power of 2 for optimal hashing)
OUTPUT = True                                             # Write info to terminal
RECORD = True                                             # Write PDFs to .txt file
MEASURE = False                                           # Take discrete measurement updates
BOUNDS = False                                            # Add inadmissible regions to grid
COLLISIONS = False;                                       # Track collisions
TV = True;                                                # Time-invariant dynamics 
BINARY = False;                                           # Binary output file
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(duffing, None, None, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, CAPACITY, DIM_f, OUTPUT, RECORD, MEASURE, BOUNDS, COLLISIONS, TV, BINARY)