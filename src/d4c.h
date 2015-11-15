//-----------------------------------------------------------------------------
// Copyright 2012-2015 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_D4C_H_
#define WORLD_D4C_H_

//-----------------------------------------------------------------------------
// Struct for D4C
//-----------------------------------------------------------------------------
typedef struct {
  double dummy;  // This is the future update.
} D4COption;

//-----------------------------------------------------------------------------
// D4C() calculates the aperiodicity estimated by D4C.
// Input:
//   x            : Input signal
//   x_length     : Length of x
//   fs           : Sampling frequency
//   time_axis    : Time axis
//   f0           : F0 contour
//   f0_length    : Length of F0 contour
//   fft_size     : The number of samples of the aperiodicity in one frame.
//                : It is given by the equation fft_size / 2 + 1.
// Output:
//   aperiodicity  : Aperiodicity estimated by D4C.
//-----------------------------------------------------------------------------
void D4C(double *x, int x_length, int fs, double *time_axis, double *f0,
  int f0_length, int fft_size, D4COption *option, double **aperiodicity);

//-----------------------------------------------------------------------------
// InitializeD4COption allocates the memory to the struct and sets the
// default parameters.
// Output:
//   option   : Struct for the optional parameter.
//-----------------------------------------------------------------------------
void InitializeD4COption(D4COption *option);

#endif  // WORLD_D4C_H_
