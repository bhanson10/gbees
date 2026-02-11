// Lorenz6D.c, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz6D
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#include "../../gbees.h" 
#include "Lorenz6D.h"

// This function defines the dynamics model - required
void Lorenz6D(double* f, double* x, double t, double* coef){
    f[0] = (x[1] - x[4]) * x[5] - x[0] + coef[0];
    f[1] = (x[2] - x[5]) * x[0] - x[1] + coef[0];
    f[2] = (x[3] - x[0]) * x[1] - x[2] + coef[0];
    f[3] = (x[4] - x[1]) * x[2] - x[3] + coef[0];
    f[4] = (x[5] - x[2]) * x[3] - x[4] + coef[0];
    f[5] = (x[0] - x[3]) * x[4] - x[5] + coef[0];
}

int main(void){
    //=================================== Read in initial discrete measurement =================================//
    printf("Reading in initial discrete measurement...\n\n");

    char* P_DIR = "./results/gbees/c"; // Saved PDFs path
    char* M_DIR = "./measurements";    // Measurement path
    char* M_FILE = "measurement0.txt"; // Measurement file
    Meas M = Meas_create(DIM_f, M_DIR, M_FILE);
    //==========================================================================================================//

    //========================================== Read in user inputs ===========================================//
    printf("Reading in user inputs...\n\n");

    double factor[DIM_f] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};              
    Grid G = Grid_create(DIM_f, 0.0, 1E-8, M, factor, false); // Inputs: (dimension, initial time, probability threshold, measure, grid width factor, rotate grid)       
  
    double coef[] = {4.0};                                    // Lorenz6D trajectory attributes (F)
    Traj T = Traj_create(1, coef);                            // Inputs: (# of coefficients, coefficients)
  
    int NUM_DIST = 2;                                         // Number of distributions recorded per measurement
    int NUM_MEAS = 1;                                         // Number of measurements
    int DEL_STEP = 20;                                        // Number of steps per deletion procedure
    int OUTPUT_FREQ = 20;                                     // Number of steps per output to terminal
    int CAPACITY = (int)pow(2,26);                            // Size of hash table (power of 2 for optimal hashing)
    bool OUTPUT = true;                                       // Write info to terminal
    bool RECORD = true;                                       // Write PDFs to .txt file
    bool MEASURE = false;                                     // Take discrete measurement updates
    bool BOUNDS = false;                                      // Add inadmissible regions to grid
    bool COLLISIONS = false;                                  // Track collisions
    bool TV = false;                                          // Time-invariant dynamics 
    bool BINARY = false;                                      // Binary output file
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    run_gbees(Lorenz6D, NULL, NULL, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, CAPACITY, DIM_f, OUTPUT, RECORD, MEASURE, BOUNDS, COLLISIONS, TV, BINARY);

    return 0;
}
