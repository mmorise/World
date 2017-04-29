//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/04/29
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

// Inverse FFT in the complex sequence
typedef struct {
  int fft_size;
  fft_complex *input;
  fft_complex *output;
  fft_plan inverse_fft;
} InverseComplexFFT;

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
//
// Input:
//   sample : Length of the input signal
//
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

//-----------------------------------------------------------------------------
// DCCorrection interpolates the power under f0 Hz
// and is used in CheapTrick() and D4C().
//-----------------------------------------------------------------------------
void DCCorrection(const double *input, double current_f0, int fs, int fft_size,
    double *output);

//-----------------------------------------------------------------------------
// LinearSmoothing() carries out the spectral smoothing by rectangular window
// whose length is width Hz and is used in CheapTrick() and D4C().
//-----------------------------------------------------------------------------
void LinearSmoothing(const double *input, double width, int fs, int fft_size,
    double *output);

//-----------------------------------------------------------------------------
// NuttallWindow() calculates the coefficients of Nuttall window whose length
// is y_length and is used in Dio(), Harvest() and D4C().
//-----------------------------------------------------------------------------
void NuttallWindow(int y_length, double *y);

//-----------------------------------------------------------------------------
// GetSafeAperiodicity() limit the range of aperiodicity from 0.001 to
// 0.999999999999 (1 - world::kMySafeGuardMinimum).
//-----------------------------------------------------------------------------
inline double GetSafeAperiodicity(double x) {
  return MyMaxDouble(0.001, MyMinDouble(0.999999999999, x));
}

//-----------------------------------------------------------------------------
// These functions are used to speed up the processing.
// Forward FFT
void InitializeForwardRealFFT(int fft_size, ForwardRealFFT *forward_real_fft);
void DestroyForwardRealFFT(ForwardRealFFT *forward_real_fft);

// Inverse FFT
void InitializeInverseRealFFT(int fft_size, InverseRealFFT *inverse_real_fft);
void DestroyInverseRealFFT(InverseRealFFT *inverse_real_fft);

// Inverse FFT (Complex)
void InitializeInverseComplexFFT(int fft_size,
  InverseComplexFFT *inverse_complex_fft);
void DestroyInverseComplexFFT(InverseComplexFFT *inverse_complex_fft);

// Minimum phase analysis (This analysis uses FFT)
void InitializeMinimumPhaseAnalysis(int fft_size,
  MinimumPhaseAnalysis *minimum_phase);
void GetMinimumPhaseSpectrum(const MinimumPhaseAnalysis *minimum_phase);
void DestroyMinimumPhaseAnalysis(MinimumPhaseAnalysis *minimum_phase);

WORLD_END_C_DECLS

#endif  // WORLD_COMMON_H_
