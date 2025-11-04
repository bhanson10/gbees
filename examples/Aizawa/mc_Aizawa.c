#include <stdio.h>
#include <stdlib.h>
#include "frame_conversion.h"

int dim = 3; 

void Lorenz3D(double t, double y[], double dydt[], int dim, double coef[]){
    dydt[0] = (y[2] - coef[1]) * y[0] - coef[3] * y[1]; 
    dydt[1] = coef[3] * y[0] + (y[2] - coef[1]) * y[1]; 
    dydt[2] = coef[2] + coef[0] * y[2] - (pow(y[2], 3.0) / 3.0) - pow(y[0], 2.0) + coef[5] * y[2] * pow(y[0], 3.0); 
    return; 
}

// 4th-order Runge-Kutta method
void RK4(void (*f)(double, double[], double[], int, double[]), double t0, double tf, int Nmc, int dim, double coef[], double mc[Nmc][dim + 1], double h, int NF, int NM, char* subfolder_name){
    double t = t0;
    int count = 1; 
    int rcount = 1; 
    size_t FILE_NAME_SIZE;
    char* FILE_NAME;
    FILE* record_file;
    double k1[dim], k2[dim], k3[dim], k4[dim], y_temp[dim];

    while (t < tf) {
        for(int j = 0; j < Nmc; j++){
            for (int i = 0; i < dim; i++) y_temp[i] = mc[j][i + 1];

            f(t, y_temp, k1, dim, coef);  // k1 = f(t, y)
            for (int i = 0; i < dim; i++) y_temp[i] = mc[j][i + 1] + h * k1[i] / 2;

            f(t + h / 2, y_temp, k2, dim, coef);  // k2 = f(t + h/2, y + h*k1/2)
            for (int i = 0; i < dim; i++) y_temp[i] = mc[j][i + 1] + h * k2[i] / 2;
            
            f(t + h / 2, y_temp, k3, dim, coef);  // k3 = f(t + h/2, y + h*k2/2)
            for (int i = 0; i < dim; i++) y_temp[i] = mc[j][i + 1] + h * k3[i];
            
            f(t + h, y_temp, k4, dim, coef);  // k4 = f(t + h, y + h*k3)
            for (int i = 0; i < dim; i++) mc[j][i + 1] += h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;

            if((fmod(count - 1, NF) == 0)&&(count != 1)){

                if(j == 0){
                    FILE_NAME_SIZE = snprintf(NULL, 0, "%s%s%s%d%s", "./results/mc/", subfolder_name, "/mc_", rcount, ".txt") + 1;
                    FILE_NAME = malloc(FILE_NAME_SIZE);
                    if (FILE_NAME == NULL) {
                        perror("malloc failed");
                        return;
                    }

                    // Create the combined string
                    snprintf(FILE_NAME, FILE_NAME_SIZE, "%s%s%s%d%s", "./results/mc/", subfolder_name, "/mc_", rcount, ".txt");
                    record_file = fopen(FILE_NAME, "w");
                    fprintf(record_file, "%.10f\n", t - t0);

                    rcount++;
                }

                fprintf(record_file, "%.10e %.10e %.10e %.10e\n", mc[j][0], mc[j][1], mc[j][2], mc[j][3]);
            }
        }
        
        t += h; count++; 
    }

    if(rcount != NM){
        FILE_NAME_SIZE = snprintf(NULL, 0, "%s%s%s%d%s", "./results/mc/", subfolder_name, "/mc_", rcount, ".txt") + 1;
        FILE_NAME = malloc(FILE_NAME_SIZE);
        if (FILE_NAME == NULL){
            perror("malloc failed");
            return;
        }

        // Create the combined string
        snprintf(FILE_NAME, FILE_NAME_SIZE, "%s%s%s%d%s", "./results/mc/", subfolder_name, "/mc_", rcount, ".txt");
        record_file = fopen(FILE_NAME, "w");
        fprintf(record_file, "%.10f\n", t - t0);

        for(int j = 0; j < Nmc; j++){
            fprintf(record_file, "%.10e %.10e %.10e %.10e\n", mc[j][0], mc[j][1], mc[j][2], mc[j][3]);
        }
    }
}

int main(int argc, char* argv[]){


    double mu[3] = {0.1, 0.0, 0.0};
    double coef[6] = {0.95, 0.7, 0.6, 3.5, 0.25, 0.1};   
    int Nmc = 1E5;
    double Sigma[3][3] = {{0.0025, 0.0, 0.0},
                          {0.0, 0.0025, 0.0},
                          {0.0, 0.0, 0.0025}};
    double L[3][3] = {0};  // Cholesky factor
    double states[Nmc][3]; // Output samples
    if (cholesky_decomposition(3, Sigma, L) != 0) {
        fprintf(stderr, "Covariance matrix is not positive definite\n");
        return -1;
    }
    mvnrnd(3, Nmc, mu, L, states);
    double mc[Nmc][4];  // Output samples + probability
    double sum = 0; 
    for(int i = 0; i < Nmc; i++){
        mc[i][0] = gauss_probability(3, states[i], mu, Sigma); sum += mc[i][0];
        for(int j = 0; j < 3; j++){
            mc[i][j + 1] = states[i][j]; 
        }
    }
    for(int i = 0; i < Nmc; i++){
        mc[i][0] /= sum;
    }

    size_t FILE_NAME_SIZE;
    char* FILE_NAME;
    FILE* record_file;

    char* subfolder_name = "C1E5";

    FILE_NAME_SIZE = snprintf(NULL, 0, "%s%s%s", "./results/mc/", subfolder_name, "/mc_0.txt") + 1;
    FILE_NAME = malloc(FILE_NAME_SIZE);
    if (FILE_NAME == NULL) {
        perror("malloc failed");
        return -1;
    }
    // Create the combined string
    snprintf(FILE_NAME, FILE_NAME_SIZE,  "%s%s%s", "./results/mc/", subfolder_name, "/mc_0.txt");
    record_file = fopen(FILE_NAME, "w");
    fprintf(record_file, "0\n");
    for(int i = 0; i < Nmc; i++){
        fprintf(record_file, "%.10e %.10e %.10e %.10e\n", mc[i][0], mc[i][1], mc[i][2], mc[i][3]);
    }

    // Final integration time (relative to ephem->jd_ref)
    double tstart = 0;
    double tend = tstart + 1;  
    int NM = 5; 
    double rt = (tend - tstart) / (NM - 1); 
    int NF = 1000; 
    double dt = rt / NF; 

    RK4(Lorenz3D, tstart, tend, Nmc, dim, coef, mc, dt, NF, NM, subfolder_name);
}
