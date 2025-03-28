// POE2BP.c, https://github.com/bhanson10/gbees/tree/main/examples/POE2BP
// Copyright 2025 by Benjamin Hanson, published under BSD 3-Clause License.

#include "../../gbees.h" 
#include "POE2BP.h"

// This function defines the dynamics model - required
void POE2BP(double* f, double* x, double t, double* coef){
    f[0] = 0;
    f[1] = pow(coef[0], 2.0) / pow(x[0], 3.0); 
}

int main(void){
    //=================================== Read in initial discrete measurement =================================//
    printf("Reading in initial discrete measurement...\n\n");

    char* P_DIR = "./results/c";       // Saved PDFs path
    char* M_DIR = "./measurements";    // Measurement path
    char* M_FILE = "measurement0.txt"; // Measurement file
    Meas M = Meas_create(DIM_f, M_DIR, M_FILE);
    //==========================================================================================================//

    //========================================== Read in user inputs ===========================================//
    printf("Reading in user inputs...\n\n");

    double factor[DIM_f] = {4.0, 2.5};
    Grid G = Grid_create(DIM_f, 0.0, 1E-9, M, factor); // Inputs: (dimension, initial time, probability threshold, measure, grid width factor)       

    double coef[] = {19.910350621818949};          // POE2BP trajectory attributes (mu)
    Traj T = Traj_create(1, coef);                 // Inputs: (# of coefficients, coefficients)

    int NUM_DIST = 5;                              // Number of distributions recorded per measurement
    int NUM_MEAS = 1;                              // Number of measurements
    int DEL_STEP = 20;                             // Number of steps per deletion procedure
    int OUTPUT_FREQ = 20;                          // Number of steps per output to terminal
    int CAPACITY = (int)pow(2,13);                 // Size of hash table (power of 2 for optimal hashing)
    bool OUTPUT = false;                           // Write info to terminal
    bool RECORD = true;                            // Write PDFs to .txt file
    bool MEASURE = false;                          // Take discrete measurement updates
    bool BOUNDS = false;                           // Add inadmissible regions to grid
    bool COLLISIONS = false;                       // Track collisions
    bool TV = false;                               // Time-invariant dynamics 
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    run_gbees(POE2BP, NULL, NULL, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, CAPACITY, DIM_f, OUTPUT, RECORD, MEASURE, BOUNDS, COLLISIONS, TV);

    return 0;
}
