//-----------------------------------------------------------------------------
// Copyright 2012-2016 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// Spectral envelope estimation on the basis of the idea of CheapTrick.
//-----------------------------------------------------------------------------
#include "world/cheaptrick.h"

#include <math.h>

#include "world/common.h"
#include "world/constantnumbers.h"
#include "world/matlabfunctions.h"

namespace {

//-----------------------------------------------------------------------------
// SmoothingWithRecovery() carries out the spectral smoothing and spectral
// recovery on the Cepstrum domain.
//-----------------------------------------------------------------------------
static void SmoothingWithRecovery(double current_f0, int fs, int fft_size, double q1,
    const ForwardRealFFT *forward_real_fft,
    const InverseRealFFT *inverse_real_fft, double *spectral_envelope) {
// We can control q1 as the parameter. 2015/9/22 by M. Morise
// const double q1 = -0.09;  // Please see the reference in CheapTrick.
  double *smoothing_lifter = new double[fft_size];
  double *compensation_lifter = new double[fft_size];

  smoothing_lifter[0] = 1;
  compensation_lifter[0] = (1.0 - 2.0 * q1) + 2.0 * q1;
  double quefrency;
  for (int i = 1; i <= forward_real_fft->fft_size / 2; ++i) {
    quefrency = static_cast<double>(i) / fs;
    smoothing_lifter[i] = sin(world::kPi * current_f0 * quefrency) /
      (world::kPi * current_f0 * quefrency);
    compensation_lifter[i] = (1.0 - 2.0 * q1) + 2.0 * q1 *
      cos(2.0 * world::kPi * quefrency * current_f0);
  }

  for (int i = 0; i <= fft_size / 2; ++i)
    forward_real_fft->waveform[i] = log(forward_real_fft->waveform[i]);
  for (int i = 1; i < fft_size / 2; ++i)
    forward_real_fft->waveform[fft_size - i] = forward_real_fft->waveform[i];
  fft_execute(forward_real_fft->forward_fft);

  for (int i = 0; i <= fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] = forward_real_fft->spectrum[i][0] *
      smoothing_lifter[i] * compensation_lifter[i] / fft_size;
    inverse_real_fft->spectrum[i][1] = 0.0;
  }
  fft_execute(inverse_real_fft->inverse_fft);

  for (int i = 0; i <= fft_size / 2; ++i)
    spectral_envelope[i] = exp(inverse_real_fft->waveform[i]);

  delete[] smoothing_lifter;
  delete[] compensation_lifter;
}

//-----------------------------------------------------------------------------
// GetPowerSpectrum() calculates the power_spectrum with DC correction.
// DC stands for Direct Current. In this case, the component from 0 to F0 Hz
// is corrected.
//-----------------------------------------------------------------------------
static void GetPowerSpectrum(int fs, double current_f0, int fft_size,
    const ForwardRealFFT *forward_real_fft) {
  int half_window_length = matlab_round(1.5 * fs / current_f0);

  // FFT
  for (int i = half_window_length * 2 + 1; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  // Calculation of the power spectrum.
  double *power_spectrum = forward_real_fft->waveform;
  for (int i = 0; i <= fft_size / 2; ++i)
    power_spectrum[i] =
      forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
      forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];

  // DC correction
  DCCorrection(power_spectrum, current_f0, fs, fft_size, power_spectrum);
}

//-----------------------------------------------------------------------------
// SetParametersForGetWindowedWaveform()
//-----------------------------------------------------------------------------
static void SetParametersForGetWindowedWaveform(int half_window_length, int x_length,
    double temporal_position, int fs, double current_f0, int *base_index,
    int *index, double *window) {
  for (int i = -half_window_length; i <= half_window_length; ++i)
    base_index[i + half_window_length] = i;
  for (int i = 0; i <= half_window_length * 2; ++i)
    index[i] = MyMinInt(x_length - 1, MyMaxInt(0,
        matlab_round(temporal_position * fs + base_index[i])));

  // Designing of the window function
  double average = 0.0;
  double position;
  double bias = temporal_position * fs - matlab_round(temporal_position * fs);
  for (int i = 0; i <= half_window_length * 2; ++i) {
    position = (static_cast<double>(base_index[i]) / 1.5 + bias) / fs;
    window[i] = 0.5 * cos(world::kPi * position * current_f0) + 0.5;
    average += window[i] * window[i];
  }
  average = sqrt(average);
  for (int i = 0; i <= half_window_length * 2; ++i) window[i] /= average;
}

//-----------------------------------------------------------------------------
// GetWindowedWaveform() windows the waveform by F0-adaptive window
//-----------------------------------------------------------------------------
static void GetWindowedWaveform(const double *x, int x_length, int fs,
    double current_f0, double temporal_position,
    const ForwardRealFFT *forward_real_fft) {
  int half_window_length = matlab_round(1.5 * fs / current_f0);

  int *base_index = new int[half_window_length * 2 + 1];
  int *index = new int[half_window_length * 2 + 1];
  double *window  = new double[half_window_length * 2 + 1];

  SetParametersForGetWindowedWaveform(half_window_length, x_length,
      temporal_position, fs, current_f0, base_index, index, window);

  // F0-adaptive windowing
  double *waveform = forward_real_fft->waveform;
  for (int i = 0; i <= half_window_length * 2; ++i)
    waveform[i] = x[index[i]] * window[i] + randn() * 0.000000000000001;
  double tmp_weight1 = 0;
  double tmp_weight2 = 0;
  for (int i = 0; i <= half_window_length * 2; ++i) {
    tmp_weight1 += waveform[i];
    tmp_weight2 += window[i];
  }
  double weighting_coefficient = tmp_weight1 / tmp_weight2;
  for (int i = 0; i <= half_window_length * 2; ++i)
    waveform[i] -= window[i] * weighting_coefficient;

  delete[] base_index;
  delete[] index;
  delete[] window;
}

//-----------------------------------------------------------------------------
// CheapTrickGeneralBody() calculates a spectral envelope at a temporal
// position. This function is only used in CheapTrick().
// Caution:
//   forward_fft is allocated in advance to speed up the processing.
//-----------------------------------------------------------------------------
static void CheapTrickGeneralBody(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double temporal_position, double q1,
    const ForwardRealFFT *forward_real_fft,
    const InverseRealFFT *inverse_real_fft, double *spectral_envelope) {
  // F0-adaptive windowing
  GetWindowedWaveform(x, x_length, fs, current_f0, temporal_position,
    forward_real_fft);

  // Calculate power spectrum with DC correction
  // Note: The calculated power spectrum is stored in an array for waveform.
  // In this imprementation, power spectrum is transformed by FFT (NOT IFFT).
  // However, the same result is obtained.
  // This is tricky but important for simple implementation.
  GetPowerSpectrum(fs, current_f0, fft_size, forward_real_fft);

  // Smoothing of the power (linear axis)
  // forward_real_fft.waveform is the power spectrum.
  LinearSmoothing(forward_real_fft->waveform, current_f0 * 2.0 / 3.0,
      fs, fft_size, forward_real_fft->waveform);

  // Smoothing (log axis) and spectral recovery on the cepstrum domain.
  SmoothingWithRecovery(current_f0, fs, fft_size, q1, forward_real_fft,
      inverse_real_fft, spectral_envelope);
}

}  // namespace

int GetFFTSizeForCheapTrick(int fs, const CheapTrickOption *option) {
  return static_cast<int>(pow(2.0, 1.0 +
      static_cast<int>(log(3.0 * fs / option->f0_floor + 1) / world::kLog2)));
}

void CheapTrick(const double *x, int x_length, int fs, const double *time_axis,
    const double *f0, int f0_length, const CheapTrickOption *option,
    double **spectrogram) {
  int fft_size = GetFFTSizeForCheapTrick(fs, option);
  double *spectral_envelope = new double[fft_size];

  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);
  InverseRealFFT inverse_real_fft = {0};
  InitializeInverseRealFFT(fft_size, &inverse_real_fft);

  double current_f0;
  for (int i = 0; i < f0_length; ++i) {
    current_f0 = f0[i] <= option->f0_floor ? world::kDefaultF0 : f0[i];
    CheapTrickGeneralBody(x, x_length, fs, current_f0, fft_size,
        time_axis[i], option->q1, &forward_real_fft, &inverse_real_fft,
        spectral_envelope);
    for (int j = 0; j <= fft_size / 2; ++j)
      spectrogram[i][j] = spectral_envelope[j];
  }

  DestroyForwardRealFFT(&forward_real_fft);
  DestroyInverseRealFFT(&inverse_real_fft);
  delete[] spectral_envelope;
}

void InitializeCheapTrickOption(CheapTrickOption *option) {
  // q1 is the parameter used for the spectral recovery.
  // Since The parameter is optimized, you don't need to change the parameter.
  option->q1 = -0.09;
  // f0_floor and fs is used to determine fft_size;
  // We strongly recommend not to change this value unless you have enough
  // knowledge of the signal processing in CheapTrick.
  option->f0_floor = world::kFloorF0;
}
