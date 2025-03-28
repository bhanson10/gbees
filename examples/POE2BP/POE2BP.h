// POE2BP.h, https://github.com/bhanson10/gbees/tree/main/examples/POE2BP
// Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

#ifndef POE2BP_H
#define POE2BP_H

#define DIM_f 2 // State dimension

// This function defines the dynamics model - required
void POE2BP(double* f, double* x, double t, double* coef);


#endif
