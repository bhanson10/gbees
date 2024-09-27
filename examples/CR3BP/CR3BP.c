// CR3BP.c, https://github.com/bhanson10/gbees/tree/main/examples/CR3BP
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#include "../../gbees.h" 
#include "CR3BP.h"

// This function defines the dynamics model - required
void CR3BP(double* f, double* x, double* dx, double* coef){
    double r1 = pow(pow(x[0]+coef[0],2) + pow(x[1],2) + pow(x[2],2), 1.5);
    double r2 = pow(pow(x[0]-1+coef[0],2) + pow(x[1],2) + pow(x[2],2), 1.5);
    f[0] = x[3];
    f[1] = x[4];
    f[2] = x[5];
    f[3] = 2*x[4]+x[0]-(coef[0]*(x[0]-1+coef[0])/r2)-((1-coef[0])*(x[0]+coef[0])/r1);
    f[4] = -2*x[3]+x[1]-(coef[0]*x[1]/r2)-((1-coef[0])*x[1]/r1);
    f[5] = -(coef[0]*x[2]/r2)-((1-coef[0])*x[2]/r1);
}

// This function defines the boundaries - optional
double CR3BP_J(double* x, double* coef){
    double r1 = pow(pow(x[0]+coef[0],2)+pow(x[1],2)+pow(x[2],2), 0.5);
    double r2 = pow(pow(x[0]-1+coef[0],2)+pow(x[1],2)+pow(x[2],2), 0.5);
    double J = pow(x[0], 2.0) + pow(x[1], 2.0) + (2*(1-coef[0])/r1) + (2*coef[0]/r2) + coef[0]*(1 - coef[0]) - (pow(x[3], 2.0) + pow(x[4], 2.0) + pow(x[5], 2.0));
    return J;
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

    double del[DIM_f];                              // Grid width, default is half of the std. dev. from the initial measurement 
    for(int i = 0; i < DIM_f; i ++){
        del[i] = pow(M.cov[i][i],0.5)/2;
    }
    Grid G = Grid_create(DIM_f, 1E-7, M.mean, del); // Inputs: (dimension, probability threshold, center, grid width)       

    double coef[] = {2.528017528540000E-5};         // CR3BP trajectory attributes (mu)
    Traj T = Traj_create(1, coef);                  // Inputs: (# of coefficients, coefficients)

    int NUM_DIST = 17;                              // Number of distributions recorded per measurement
    int NUM_MEAS = 1;                               // Number of measurements
    int DEL_STEP = 20;                              // Number of steps per deletion procedure
    int OUTPUT_FREQ = 20;                           // Number of steps per output to terminal
    bool OUTPUT = true;                            // Write info to terminal after certain amount of steps
    bool RECORD = true;                             // Write PDFs to .txt file
    bool MEASURE = false;                           // Take discrete measurement updates
    bool BOUNDS = true;                             // Add inadmissible regions to grid
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    run_gbees(CR3BP, NULL, CR3BP_J, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM_h, OUTPUT, RECORD, MEASURE, BOUNDS);

    return 0;
}
