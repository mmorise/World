//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef WORLD_D4C_H_
#define WORLD_D4C_H_

#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// Struct for D4C
//-----------------------------------------------------------------------------
typedef struct {
  double threshold;
} D4COption;

//-----------------------------------------------------------------------------
// D4C() calculates the aperiodicity estimated by D4C.
//
// Input:
//   x                  : Input signal
//   x_length           : Length of x
//   fs                 : Sampling frequency
//   temporal_positions : Time axis
//   f0                 : F0 contour
//   f0_length          : Length of F0 contour
//   fft_size           : Number of samples of the aperiodicity in one frame.
//                      : It is given by the equation fft_size / 2 + 1.
// Output:
//   aperiodicity  : Aperiodicity estimated by D4C.
//-----------------------------------------------------------------------------
void D4C(const double *x, int x_length, int fs,
    const double *temporal_positions, const double *f0, int f0_length,
    int fft_size, const D4COption *option, double **aperiodicity);

//-----------------------------------------------------------------------------
// InitializeD4COption allocates the memory to the struct and sets the
// default parameters.
//
// Output:
//   option   : Struct for the optional parameter.
//-----------------------------------------------------------------------------
void InitializeD4COption(D4COption *option);

WORLD_END_C_DECLS

#endif  // WORLD_D4C_H_
