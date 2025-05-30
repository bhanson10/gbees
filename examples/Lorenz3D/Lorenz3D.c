// Lorenz3D.c, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz3D
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#include "../../gbees.h" 
#include "Lorenz3D.h"

// This function defines the dynamics model - required
void Lorenz3D(double* f, double* x, double t, double* coef){
    f[0] = coef[0] * (x[1] - x[0]);
    f[1] = -x[1] - x[0] * x[2];
    f[2] = -coef[1] * x[2] + x[0] * x[1] - coef[1] * coef[2];
}

// This function defines the measurement model - required if MEASURE == true
void z(double* h, double* x, double t, double* coef){
    h[0] = x[2];
}

int main(void){
    //=================================== Read in initial discrete measurement =================================//
    printf("Reading in initial discrete measurement...\n\n");

    char* P_DIR = "./results/c";     // Saved PDFs path
    char* M_DIR = "./measurements";    // Measurement path
    char* M_FILE = "measurement0.txt"; // Measurement file
    Meas M = Meas_create(DIM_f, M_DIR, M_FILE);
    //==========================================================================================================//

    //========================================== Read in user inputs ===========================================//
    printf("Reading in user inputs...\n\n");

    double factor[DIM_f] = {1.0, 1.0, 1.0}; 
    Grid G = Grid_create(DIM_f, 0.0, 1E-7, M, factor, false); // Inputs: (dimension, initial time, probability threshold, measurement, grid width factor, rotate grid)       

    double coef[] = {4.0, 1.0, 48.0};                         // Lorenz3D trajectory attributes (sigma, beta, r)
    Traj T = Traj_create(3, coef);                            // Inputs: (# of coefficients, coefficients)

    int NUM_DIST = 4;                                         // Number of distributions recorded per measurement
    int NUM_MEAS = 2;                                         // Number of measurements
    int DEL_STEP = 20;                                        // Number of steps per deletion procedure
    int OUTPUT_FREQ = 20;                                     // Number of steps per output to terminal
    int CAPACITY = (int)pow(2,12);                            // Size of hash table (power of 2 for optimal hashing)
    bool OUTPUT = true;                                       // Write info to terminal
    bool RECORD = true;                                       // Write PDFs to .txt file
    bool MEASURE = true;                                      // Take discrete measurement updates
    bool BOUNDS = false;                                      // Add inadmissible regions to grid
    bool COLLISIONS = false;                                  // Track collisions
    bool TV = false;                                          // Time-invariant dynamics 
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    run_gbees(Lorenz3D, z, NULL, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, CAPACITY, DIM_h, OUTPUT, RECORD, MEASURE, BOUNDS, COLLISIONS, TV);

    return 0;
}
