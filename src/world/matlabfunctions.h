//-----------------------------------------------------------------------------
// Copyright 2012-2016 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_MATLABFUNCTIONS_H_
#define WORLD_MATLABFUNCTIONS_H_

#include "world/common.h"
#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// fftshift() swaps the left and right halves of input vector.
// http://www.mathworks.com/help/matlab/ref/fftshift.html
// Input:
//   x              : Input vector
//   x_length       : Length of x
// Output:
//   y              : Swapped vector x
// Caution:
//   Lengths of index and edges must be the same.
//-----------------------------------------------------------------------------
void fftshift(const double *x, int x_length, double *y);

//-----------------------------------------------------------------------------
// histc() counts the number of values in vector x that fall between the
// elements in the edges vector (which must contain monotonically
// nondecreasing values). n is a length(edges) vector containing these counts.
// No elements of x can be complex.
// http://www.mathworks.co.jp/help/techdoc/ref/histc.html
// Input:
//   x              : Input vector
//   x_length       : Length of x
//   edges          : Input matrix (1-dimension)
//   edges_length   : Length of edges
// Output:
//   index          : Result counted in vector x
// Caution:
//   Lengths of index and edges must be the same.
//-----------------------------------------------------------------------------
void histc(const double *x, int x_length, const double *edges,
  int edges_length, int *index);

//-----------------------------------------------------------------------------
// interp1() interpolates to find yi, the values of the underlying function Y
// at the points in the vector or array xi. x must be a vector.
// http://www.mathworks.co.jp/help/techdoc/ref/interp1.html
// Input:
//   x          : Input vector (Time axis)
//   y          : Values at x[n]
//   x_length   : Length of x (Length of y must be the same)
//   xi         : Required vector
//   xi_length  : Length of xi (Length of yi must be the same)
// Output:
//   yi         : Interpolated vector
//-----------------------------------------------------------------------------
void interp1(const double *x, const double *y, int x_length, const double *xi,
  int xi_length, double *yi);

//-----------------------------------------------------------------------------
// decimate() carries out down sampling by both IIR and FIR filters.
// Filter coeffiencts are based on FilterForDecimate().
// Input:
//   x          : Input signal
//   x_length   : Length of x
//   r          : Coefficient used for down sampling
//                (fs after down sampling is fs/r)
// Output:
//   y          : Output signal
//-----------------------------------------------------------------------------
void decimate(const double *x, int x_length, int r, double *y);

//-----------------------------------------------------------------------------
// matlab_round() calculates rounding.
// Input:
//   x    : Input value
// Output:
//   y    : Rounded value
//-----------------------------------------------------------------------------
int matlab_round(double x);

//-----------------------------------------------------------------------------
// diff() calculates differences and approximate derivatives
// http://www.mathworks.co.jp/help/techdoc/ref/diff.html
// Input:
//   x          : Input signal
//   x_length   : Length of x
// Output:
//   y          : Output signal
//-----------------------------------------------------------------------------
void diff(const double *x, int x_length, double *y);

//-----------------------------------------------------------------------------
// interp1Q() is the special case of interp1().
// We can use this function, provided that All periods of x-axis is the same.
// Input:
//   x          : Origin of the x-axis
//   shift      : Period of the x-axis
//   y          : Values at x[n]
//   x_length   : Length of x (Length of y must be the same)
//   xi         : Required vector
//   xi_length  : Length of xi (Length of yi must be the same)
// Output:
//   yi         : Interpolated vector
// Caution:
//   Length of xi and yi must be the same.
//-----------------------------------------------------------------------------
void interp1Q(double x, double shift, const double *y, int x_length,
  const double *xi, int xi_length, double *yi);

//-----------------------------------------------------------------------------
// randn() generates pseudorandom numbers based on xorshift method.
// Output:
//   A generated pseudorandom number
//-----------------------------------------------------------------------------
double randn(void);

//-----------------------------------------------------------------------------
// fast_fftfilt() carries out the convolution on the frequency domain.
// Input:
//   x                : Input signal
//   x_length         : Length of x
//   h                : Impulse response
//   h_length         : Length of h
//   fft_size         : Length of FFT
//   forward_real_fft : Struct to speed up the forward FFT
//   inverse_real_fft : Struct to speed up the inverse FFT
// Output:
//   y                : Calculated result.
//-----------------------------------------------------------------------------
void fast_fftfilt(const double *x, int x_length, const double *h, int h_length,
  int fft_size, const ForwardRealFFT *forward_real_fft,
  const InverseRealFFT *inverse_real_fft, double *y);

//-----------------------------------------------------------------------------
// matlab_std() calculates the standard deviation of the input vector.
// Input:
//   x          : Input vector
//   x_length   : Length of x
// Output:
//   Calculated standard deviation
//-----------------------------------------------------------------------------
double matlab_std(const double *x, int x_length);

WORLD_END_C_DECLS

#endif  // WORLD_MATLABFUNCTIONS_H_
