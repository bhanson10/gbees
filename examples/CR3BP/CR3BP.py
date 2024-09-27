# CR3BP.py, https://github.com/bhanson10/gbees/tree/main/examples/CR3BP
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

import sys
sys.path.append('../../')
import gbeespy as gbees  # type: ignore

DIM_f = 6 # State dimension
DIM_h = 6 # Measurement dimension

# This function defines the dynamics model - required
def CR3BP(x, dx, coef):
    r1 = ((x[0] + coef[0])**2 + (x[1])**2 + (x[2])**2)**1.5
    r2 = ((x[0] - 1 + coef[0])**2 + (x[1])**2 + (x[2])**2)**1.5
    f1 = x[3]
    f2 = x[4]
    f3 = x[5]
    f4 = 2*x[4]+x[0]-(coef[0]*(x[0]-1+coef[0])/r2)-((1-coef[0])*(x[0]+coef[0])/r1)
    f5 = -2*x[3]+x[1]-(coef[0]*x[1]/r2)-((1-coef[0])*x[1]/r1)
    f6 = -(coef[0]*x[2]/r2)-((1-coef[0])*x[2]/r1)
    return [f1, f2, f3, f4, f5, f6]

# This function defines the initial grid boundaries - optional
def CR3BP_J(x, coef):
    r1 = ((x[0]+coef[0])**2+(x[1])**2+(x[2])**2)**0.5
    r2 = ((x[0]-1+coef[0])**2+(x[1])**2+(x[1])**2)**0.5
    J = (x[0])**2.0 + (x[1])**2.0 + (2*(1-coef[0])/r1) + (2*coef[0]/r2) + coef[0]*(1 - coef[0]) - ((x[3])**2.0 + (x[4])**2.0 + (x[5])**2.0)
    return J

#==================================== Read in initial discrete measurement ==================================#
print("Reading in initial discrete measurement...\n")

P_DIR = "./results/python"      # Saved PDFs path
M_DIR = "./measurements"     # Measurement path
M_FILE = "measurement0.txt"; # Measurement file
M = gbees.Meas_create(DIM_f, M_DIR, M_FILE)
#============================================================================================================#

#=========================================== Read in user inputs ============================================#
print("Reading in user inputs...\n")

dx = [None] * DIM_f                             # Grid width, default is half of the std. dev. from the initial measurement 
for i in range(DIM_f):
    dx[i] = (M.cov[i][i]**(0.5))/2
G = gbees.Grid_create(DIM_f, 1E-7, M.mean, dx); # Inputs: (dimension, probability threshold, center, grid width)    
 
coef = [2.528017528540000E-5]                   # PCR3BP trajectory attributes (mu)
T = gbees.Traj_create(len(coef), coef);         # Inputs: (# of coefficients, coefficients)

NUM_DIST = 8;                                   # Number of distributions recorded per measurement
NUM_MEAS = 4;                                   # Number of measurements
DEL_STEP = 20;                                  # Number of steps per deletion procedure
OUTPUT_FREQ = 20;                               # Number of steps per output to terminal
OUTPUT = True;                                 # Write info to terminal
RECORD = True;                                  # Write PDFs to .txt file
MEASURE = False;                                # Take discrete measurement updates
BOUNDS = True;                                  # Add inadmissible regions to grid
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(CR3BP, None, CR3BP_J, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM_h, OUTPUT, RECORD, MEASURE, BOUNDS)