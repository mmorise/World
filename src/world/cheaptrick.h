//-----------------------------------------------------------------------------
// Copyright 2012-2016 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_CHEAPTRICK_H_
#define WORLD_CHEAPTRICK_H_

#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// Struct for CheapTrick
//-----------------------------------------------------------------------------
typedef struct {
  // This is defined as the struct for future update.
  double q1;
  double f0_floor;
} CheapTrickOption;

//-----------------------------------------------------------------------------
// CheapTrick() calculates the spectrogram that consists of spectral envelopes
// estimated by CheapTrick.
// Input:
//   x            : Input signal
//   x_length     : Length of x
//   fs           : Sampling frequency
//   time_axis    : Time axis
//   f0           : F0 contour
//   f0_length    : Length of F0 contour
//   option       : Struct to order the parameter for CheapTrick
// Output:
//   spectrogram  : Spectrogram estimated by CheapTrick.
//-----------------------------------------------------------------------------
void CheapTrick(const double *x, int x_length, int fs, const double *time_axis,
  const double *f0, int f0_length, const CheapTrickOption *option,
  double **spectrogram);

//-----------------------------------------------------------------------------
// InitializeCheapTrickOption allocates the memory to the struct and sets the
// default parameters.
// Output:
//   option   : Struct for the optional parameter.
//-----------------------------------------------------------------------------
void InitializeCheapTrickOption(CheapTrickOption *option);

//-----------------------------------------------------------------------------
// GetFFTSizeForCheapTrick() calculates the FFT size based on the sampling
// frequency and the lower limit of f0 (It is defined in world.h).
// Input:
//   fs      : Sampling frequency
// Output:
//   FFT size
//-----------------------------------------------------------------------------
int GetFFTSizeForCheapTrick(int fs, const CheapTrickOption *option);

WORLD_END_C_DECLS

#endif  // WORLD_CHEAPTRICK_H_
