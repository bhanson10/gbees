// PCR3BP.h, https://github.com/bhanson10/gbees/tree/main/examples/PCR3BP
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef PCR3BP_H
#define PCR3BP_H

#define DIM_f 4 // State dimension
#define DIM_h 3 // Measurement dimension

// This function defines the dynamics model - required
void PCR3BP(double* f, double* x, double t, double* dx, double* coef);

// This function defines the measurement model - required if MEASURE == true
void rtrr(double* h, double* x, double t, double* dx, double* coef);

// This function defines the initial grid boundaries - optional
double PCR3BP_J(double* x, double* coef);

#endif
