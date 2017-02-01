//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//-----------------------------------------------------------------------------
#ifndef WORLD_STONEMASK_H_
#define WORLD_STONEMASK_H_

#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// StoneMask() refines the estimated F0 by Dio()
//
// Input:
//   x                      : Input signal
//   x_length               : Length of the input signal
//   fs                     : Sampling frequency
//   time_axis              : Temporal information
//   f0                     : f0 contour
//   f0_length              : Length of f0
//
// Output:
//   refined_f0             : Refined F0
//-----------------------------------------------------------------------------
void StoneMask(const double *x, int x_length, int fs,
    const double *temporal_positions, const double *f0, int f0_length,
    double *refined_f0);

WORLD_END_C_DECLS

#endif  // WORLD_STONEMASK_H_
