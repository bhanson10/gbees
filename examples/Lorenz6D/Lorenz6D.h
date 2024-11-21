// Lorenz3D.h, https://github.com/bhanson10/gbees/tree/main/examples/Lorenz3D
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef LORENZ6D_H
#define LORENZ6D_H

#define DIM_f 6 // State dimension

// This function defines the dynamics model - required
void Lorenz6D(double* f, double* x, double t, double* dx, double* coef);


#endif
