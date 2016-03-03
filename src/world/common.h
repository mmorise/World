//-----------------------------------------------------------------------------
// Copyright 2012-2016 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//-----------------------------------------------------------------------------
#ifndef WORLD_COMMON_H_
#define WORLD_COMMON_H_

#include "world/fft.h"
#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// Structs on FFT
//-----------------------------------------------------------------------------
// Forward FFT in the real sequence
typedef struct {
  int fft_size;
  double *waveform;
  fft_complex *spectrum;
  fft_plan forward_fft;
} ForwardRealFFT;

// Inverse FFT in the real sequence
typedef struct {
  int fft_size;
  double *waveform;
  fft_complex *spectrum;
  fft_plan inverse_fft;
} InverseRealFFT;

// Minimum phase analysis from logarithmic power spectrum
typedef struct {
  int fft_size;
  double *log_spectrum;
  fft_complex *minimum_phase_spectrum;
  fft_complex *cepstrum;
  fft_plan inverse_fft;
  fft_plan forward_fft;
} MinimumPhaseAnalysis;

//-----------------------------------------------------------------------------
// GetSuitableFFTSize() calculates the suitable FFT size.
// The size is defined as the minimum length whose length is longer than
// the input sample.
// Input:
//   sample     : Length of the input signal
// Output:
//   Suitable FFT size
//-----------------------------------------------------------------------------
int GetSuitableFFTSize(int sample);

//-----------------------------------------------------------------------------
// These four functions are simple max() and min() function
// for "int" and "double" type.
//-----------------------------------------------------------------------------
inline int MyMaxInt(int x, int y) {
  return x > y ? x : y;
}

inline double MyMaxDouble(double x, double y) {
  return x > y ? x : y;
}

inline int MyMinInt(int x, int y) {
  return x < y ? x : y;
}

inline double MyMinDouble(double x, double y) {
  return x < y ? x : y;
}

//-----------------------------------------------------------------------------
// These functions are used in at least two different .cpp files
// DCCorrection is used in CheapTrick() and D4C().
void DCCorrection(const double *input, double current_f0, int fs, int fft_size,
    double *output);

// LinearSmoothing is used in CheapTrick() and D4C().
void LinearSmoothing(const double *input, double width, int fs, int fft_size,
    double *output);

// NuttallWindow is used in Dio() and D4C().
void NuttallWindow(int y_length, double *y);

//-----------------------------------------------------------------------------
// These functions are used to speed up the processing.
// Forward FFT
void InitializeForwardRealFFT(int fft_size, ForwardRealFFT *forward_real_fft);
void DestroyForwardRealFFT(ForwardRealFFT *forward_real_fft);

// Inverse FFT
void InitializeInverseRealFFT(int fft_size, InverseRealFFT *inverse_real_fft);
void DestroyInverseRealFFT(InverseRealFFT *inverse_real_fft);

// Minimum phase analysis (This analysis uses FFT)
void InitializeMinimumPhaseAnalysis(int fft_size,
  MinimumPhaseAnalysis *minimum_phase);
void GetMinimumPhaseSpectrum(const MinimumPhaseAnalysis *minimum_phase);
void DestroyMinimumPhaseAnalysis(MinimumPhaseAnalysis *minimum_phase);

WORLD_END_C_DECLS

#endif  // WORLD_COMMON_H_
