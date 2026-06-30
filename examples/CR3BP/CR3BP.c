// CR3BP.c, https://github.com/bhanson10/gbees/tree/main/examples/CR3BP
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#include "../../gbees.h" 
#include "CR3BP.h"

// This function defines the dynamics model - required
void CR3BP(double* f, double* x, double t, double* coef){
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

double* init_v(int N){
    double* v = (double*)malloc(N * sizeof(double));
    if (v == NULL) {
        perror("init_v: row malloc failed");
        return NULL;
    }

    return v; 
}

double** cast_matrix(int N, int M, double matrix[N][M]){

    double** result = init_m(N, M); 
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            result[i][j] = matrix[i][j];
        }
    }

    return result; 
}

double* mvtimes(int N, int M, double** matrix_a, const double* vector_b){

    double* result = init_v(N); 

    for (int i = 0; i < N; i++) {
        result[i] = 0.0;
        for (int j = 0; j < M; j++) {
            result[i] += matrix_a[i][j] * vector_b[j];
        }
    }
    return result;
}

double* synodic_to_inertial(double* rv_S, double t, double mu){

  double* rv_I = init_v(6);
  double* r_S = init_v(3);
  double* v_S = init_v(3);
  double Q_IS[3][3] = {{cos(t), -sin(t), 0}, 
                       {sin(t), cos(t),  0}, 
                       {0,      0,       1}};
  double** Q_IS_m = cast_matrix(3, 3, Q_IS);
  double Qdot_IS[3][3] = {{-sin(t), -cos(t), 0}, 
                          {cos(t), -sin(t),  0}, 
                          {0,      0,        0}};
  double** Qdot_IS_m = cast_matrix(3, 3, Qdot_IS);

  r_S[0] = rv_S[0] + mu; r_S[1] = rv_S[1]; r_S[2] = rv_S[2]; 
  v_S[0] = rv_S[3]; v_S[1] = rv_S[4]; v_S[2] = rv_S[5];

  double* r_I = mvtimes(3, 3, Q_IS_m, r_S); 
  double* v_I_a = mvtimes(3, 3, Q_IS_m, v_S); 
  double* v_I_b = mvtimes(3, 3, Qdot_IS_m, r_S); 

  rv_I[0] = r_I[0];
  rv_I[1] = r_I[1]; 
  rv_I[2] = r_I[2]; 
  rv_I[3] = (v_I_a[0] + v_I_b[0]); 
  rv_I[4] = (v_I_a[1] + v_I_b[1]); 
  rv_I[5] = (v_I_a[2] + v_I_b[2]); 
  
  for(int i = 0; i < 3; i++){
    free(Q_IS_m[i]); 
    free(Qdot_IS_m[i]); 
  }
  free(r_S); 
  free(v_S); 
  free(r_I); 
  free(v_I_a); 
  free(v_I_b); 

  return rv_I; 
}

void R_RR(double* h, double* x, double t, double* coef){

    double* x_I = synodic_to_inertial(x, t, coef[0]);
    double x_g[6] = {-7.0348559842e-03, -5.4507116728e-03, 1.3714177099e-02, 2.9202056747e-01, -1.8122240010e-01, 7.7768543989e-02};
    double* rho = init_v(3); 
    double* rho_dot = init_v(3); 
    for(int i = 0; i < 3; i++){
        rho[i] = x_I[i] - x_g[i]; 
        rho_dot[i] = x_I[i + 3] - x_g[i + 3]; 
    }

    h[0] = sqrt(pow(rho[0], 2.0) + pow(rho[1], 2.0) + pow(rho[2], 2.0)); // range
    h[1] = dot(rho, rho_dot, 3) / h[0]; 

    return;
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

    double factor[DIM_f] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    Grid G = Grid_create(DIM_f, 0.0, 1E-9, M, factor, false); // Inputs: (dimension, initial time, probability threshold, measurement, grid width factor, rotate grid)       

    double coef[] = {1.215058446919780E-2};                   // CR3BP trajectory attributes (mu)
    Traj T = Traj_create(1, coef);                            // Inputs: (# of coefficients, coefficients)

    int NUM_DIST = 2;                                         // Number of distributions recorded per measurement
    int NUM_MEAS = 1;                                         // Number of measurements
    int DEL_STEP = 20;                                        // Number of steps per deletion procedure
    int OUTPUT_FREQ = 20;                                     // Number of steps per output to terminal
    int CAPACITY = (int)pow(2,27);                            // Size of hash table (power of 2 for optimal hashing)
    bool OUTPUT = true;                                       // Write info to terminal after certain amount of steps
    bool RECORD = true;                                       // Write PDFs to .txt file
    bool MEASURE = false;                                     // Take discrete measurement updates
    bool BOUNDS = true;                                       // Add inadmissible regions to grid
    bool COLLISIONS = false;                                  // Track collisions
    bool TV = false;                                          // Time-invariant dynamics 
    bool BINARY = true;                                       // Binary output file
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    run_gbees(CR3BP, NULL, CR3BP_J, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, CAPACITY, DIM_h, OUTPUT, RECORD, MEASURE, BOUNDS, COLLISIONS, TV, BINARY);
    return 0;
}
