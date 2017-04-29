//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/04/29
//-----------------------------------------------------------------------------
#ifndef WORLD_CODEC_H_
#define WORLD_CODEC_H_

#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// GetNumberOfAperiodicities provides the number of dimensions for aperiodicity
// coding. It is determined by only fs.
//
// Input:
//   fs       : Sampling frequency
//
// Output:
//   Number of aperiodicities
//-----------------------------------------------------------------------------
int GetNumberOfAperiodicities(int fs);

//-----------------------------------------------------------------------------
// CodeAperiodicity codes the aperiodicity. The number of dimensions is
// determined by fs.
//
// Input:
//   aperiodicity       : Aperiodicity before coding
//   f0_length          : Length of F0 contour
//   fs                 : Sampling frequency
//   fft_size           : FFT size
//
// Output:
//   coded_aperiodicity : Coded aperiodicity
//-----------------------------------------------------------------------------
void CodeAperiodicity(const double * const *aperiodicity, int f0_length,
  int fs, int fft_size, double **coded_aperiodicity);

//-----------------------------------------------------------------------------
// DecodeAperiodicity decodes the coded aperiodicity.
//
// Input:
//   coded_aperiodicity : Coded aperiodicity
//   f0_length          : Length of F0 contour
//   fs                 : Sampling frequency
//   fft_size           : FFT size
//
// Output:
//   aperiodicity       : Decoded aperiodicity
//-----------------------------------------------------------------------------
void DecodeAperiodicity(const double * const *coded_aperiodicity,
  int f0_length, int fs, int fft_size, double **aperiodicity);

//-----------------------------------------------------------------------------
// CodeSpectralEnvelope codes the spectral envelope.
//
// Input:
//   aperiodicity         : Aperiodicity before coding
//   f0_length            : Length of F0 contour
//   fs                   : Sampling frequency
//   fft_size             : FFT size
//   number_of_dimensions : Parameter for compression
//
// Output:
//   coded_spectral_envelope
//-----------------------------------------------------------------------------
void CodeSpectralEnvelope(const double * const *spectrogram, int f0_length,
  int fs, int fft_size, int number_of_dimensions,
  double **coded_spectral_envelope);

//-----------------------------------------------------------------------------
// DecodeSpectralEnvelope decodes the coded spectral envelope.
//
// Input:
//   coded_aperiodicity   : Coded aperiodicity
//   f0_length            : Length of F0 contour
//   fs                   : Sampling frequency
//   fft_size             : FFT size
//   number_of_dimensions : Parameter for compression
//
// Output:
//   spectrogram
//-----------------------------------------------------------------------------
void DecodeSpectralEnvelope(const double * const *coded_spectral_envelope,
  int f0_length, int fs, int fft_size, int number_of_dimensions,
  double **spectrogram);

WORLD_END_C_DECLS

#endif  // WORLD_CODEC_H_
