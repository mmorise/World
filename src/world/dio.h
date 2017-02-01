//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//-----------------------------------------------------------------------------
#ifndef WORLD_DIO_H_
#define WORLD_DIO_H_

#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// Struct for DIO
//-----------------------------------------------------------------------------
typedef struct {
  double f0_floor;
  double f0_ceil;
  double channels_in_octave;
  double frame_period;  // msec
  int speed;  // (1, 2, ..., 12)
  double allowed_range;  // Threshold used for fixing the F0 contour.
} DioOption;

//-----------------------------------------------------------------------------
// DIO
//
// Input:
//   x                    : Input signal
//   x_length             : Length of x
//   fs                   : Sampling frequency
//   option               : Struct to order the parameter for DIO
//
// Output:
//   temporal_positions   : Temporal positions.
//   f0                   : F0 contour.
//-----------------------------------------------------------------------------
void Dio(const double *x, int x_length, int fs, const DioOption *option,
  double *temporal_positions, double *f0);

//-----------------------------------------------------------------------------
// InitializeDioOption allocates the memory to the struct and sets the
// default parameters.
//
// Output:
//   option   : Struct for the optional parameter.
//-----------------------------------------------------------------------------
void InitializeDioOption(DioOption *option);

//-----------------------------------------------------------------------------
// GetSamplesForDIO() calculates the number of samples required for Dio().
//
// Input:
//   fs             : Sampling frequency [Hz]
//   x_length       : Length of the input signal [Sample].
//   frame_period   : Frame shift [msec]
//
// Output:
//   The number of samples required to store the results of Dio()
//-----------------------------------------------------------------------------
int GetSamplesForDIO(int fs, int x_length, double frame_period);

WORLD_END_C_DECLS

#endif  // WORLD_DIO_H_
