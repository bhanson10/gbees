// CR3BP.h, https://github.com/bhanson10/gbees/tree/main/examples/CR3BP
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef CR3BP_H
#define CR3BP_H

#define DIM_f 6 // State dimension
#define DIM_h 6 // Measurement dimension

// This function defines the dynamics model - required
void CR3BP(double* f, double* x, double t, double* dx, double* coef);

// This function defines the boundaries - optional
double CR3BP_J(double* x, double* coef);

#endif
