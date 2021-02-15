//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//
// Band-aperiodicity estimation on the basis of the idea of D4C.
//-----------------------------------------------------------------------------
#include "world/d4c.h"

#include <math.h>
#include <algorithm>  // for std::sort()

#include "world/common.h"
#include "world/constantnumbers.h"
#include "world/matlabfunctions.h"

namespace {
//-----------------------------------------------------------------------------
// SetParametersForGetWindowedWaveform()
//-----------------------------------------------------------------------------
static void SetParametersForGetWindowedWaveform(int half_window_length,
    int x_length, double current_position, int fs, double current_f0,
    int window_type, double window_length_ratio, int *base_index,
    int *safe_index, double *window) {
  for (int i = -half_window_length; i <= half_window_length; ++i)
    base_index[i + half_window_length] = i;
  int origin = matlab_round(current_position * fs + 0.001);
  for (int i = 0; i <= half_window_length * 2; ++i)
    safe_index[i] =
      MyMinInt(x_length - 1, MyMaxInt(0, origin + base_index[i]));

  // Designing of the window function
  double position;
  if (window_type == world::kHanning) {  // Hanning window
    for (int i = 0; i <= half_window_length * 2; ++i) {
      position = (2.0 * base_index[i] / window_length_ratio) / fs;
      window[i] = 0.5 * cos(world::kPi * position * current_f0) + 0.5;
    }
  } else {  // Blackman window
    for (int i = 0; i <= half_window_length * 2; ++i) {
      position = (2.0 * base_index[i] / window_length_ratio) / fs;
      window[i] = 0.42 + 0.5 * cos(world::kPi * position * current_f0) +
        0.08 * cos(world::kPi * position * current_f0 * 2);
    }
  }
}

//-----------------------------------------------------------------------------
// GetWindowedWaveform() windows the waveform by F0-adaptive window
// In the variable window_type, 1: hanning, 2: blackman
//-----------------------------------------------------------------------------
static void GetWindowedWaveform(const double *x, int x_length, int fs,
    double current_f0, double current_position, int window_type,
    double window_length_ratio, double *waveform) {
  int half_window_length =
    matlab_round(window_length_ratio * fs / current_f0 / 2.0);

  int *base_index = new int[half_window_length * 2 + 1];
  int *safe_index = new int[half_window_length * 2 + 1];
  double *window  = new double[half_window_length * 2 + 1];

  SetParametersForGetWindowedWaveform(half_window_length, x_length,
      current_position, fs, current_f0, window_type, window_length_ratio,
      base_index, safe_index, window);

  // F0-adaptive windowing
  for (int i = 0; i <= half_window_length * 2; ++i)
    waveform[i] =
      x[safe_index[i]] * window[i] + randn() * world::kMySafeGuardMinimum;

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
  delete[] safe_index;
  delete[] window;
}

//-----------------------------------------------------------------------------
// GetCentroid() calculates the energy centroid (see the book, time-frequency
// analysis written by L. Cohen).
//-----------------------------------------------------------------------------
static void GetCentroid(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double current_position,
    const ForwardRealFFT *forward_real_fft, double *centroid) {
  for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;
  GetWindowedWaveform(x, x_length, fs, current_f0,
      current_position, world::kBlackman, 4.0, forward_real_fft->waveform);
  double power = 0.0;
  for (int i = 0; i <= matlab_round(2.0 * fs / current_f0) * 2; ++i)
    power += forward_real_fft->waveform[i] * forward_real_fft->waveform[i];
  for (int i = 0; i <= matlab_round(2.0 * fs / current_f0) * 2; ++i)
    forward_real_fft->waveform[i] /= sqrt(power);

  fft_execute(forward_real_fft->forward_fft);
  double *tmp_real = new double[fft_size / 2 + 1];
  double *tmp_imag = new double[fft_size / 2 + 1];
  for (int i = 0; i <= fft_size / 2; ++i) {
    tmp_real[i] = forward_real_fft->spectrum[i][0];
    tmp_imag[i] = forward_real_fft->spectrum[i][1];
  }

  for (int i = 0; i < fft_size; ++i)
    forward_real_fft->waveform[i] *= i + 1.0;
  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i)
    centroid[i] = forward_real_fft->spectrum[i][0] * tmp_real[i] +
      tmp_imag[i] * forward_real_fft->spectrum[i][1];

  delete[] tmp_real;
  delete[] tmp_imag;
}

//-----------------------------------------------------------------------------
// GetStaticCentroid() calculates the temporally static energy centroid.
// Basic idea was proposed by H. Kawahara.
//-----------------------------------------------------------------------------
static void GetStaticCentroid(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double current_position,
    const ForwardRealFFT *forward_real_fft, double *static_centroid) {
  double *centroid1 = new double[fft_size / 2 + 1];
  double *centroid2 = new double[fft_size / 2 + 1];

  GetCentroid(x, x_length, fs, current_f0, fft_size,
      current_position - 0.25 / current_f0, forward_real_fft, centroid1);
  GetCentroid(x, x_length, fs, current_f0, fft_size,
      current_position + 0.25 / current_f0, forward_real_fft, centroid2);

  for (int i = 0; i <= fft_size / 2; ++i)
    static_centroid[i] = centroid1[i] + centroid2[i];

  DCCorrection(static_centroid, current_f0, fs, fft_size, static_centroid);
  delete[] centroid1;
  delete[] centroid2;
}

//-----------------------------------------------------------------------------
// GetSmoothedPowerSpectrum() calculates the smoothed power spectrum.
// The parameters used for smoothing are optimized in davance.
//-----------------------------------------------------------------------------
static void GetSmoothedPowerSpectrum(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double current_position,
    const ForwardRealFFT *forward_real_fft, double *smoothed_power_spectrum) {
  for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;
  GetWindowedWaveform(x, x_length, fs, current_f0,
      current_position, world::kHanning, 4.0, forward_real_fft->waveform);

  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i)
    smoothed_power_spectrum[i] =
      forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
      forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];
  DCCorrection(smoothed_power_spectrum, current_f0, fs, fft_size,
      smoothed_power_spectrum);
  LinearSmoothing(smoothed_power_spectrum, current_f0, fs, fft_size,
      smoothed_power_spectrum);
}

//-----------------------------------------------------------------------------
// GetStaticGroupDelay() calculates the temporally static group delay.
// This is the fundamental parameter in D4C.
//-----------------------------------------------------------------------------
static void GetStaticGroupDelay(const double *static_centroid,
    const double *smoothed_power_spectrum, int fs, double f0,
    int fft_size, double *static_group_delay) {
  for (int i = 0; i <= fft_size / 2; ++i)
    static_group_delay[i] = static_centroid[i] / smoothed_power_spectrum[i];
  LinearSmoothing(static_group_delay, f0 / 2.0, fs, fft_size,
      static_group_delay);

  double *smoothed_group_delay = new double[fft_size / 2 + 1];
  LinearSmoothing(static_group_delay, f0, fs, fft_size,
      smoothed_group_delay);

  for (int i = 0; i <= fft_size / 2; ++i)
    static_group_delay[i] -= smoothed_group_delay[i];

  delete[] smoothed_group_delay;
}

//-----------------------------------------------------------------------------
// GetCoarseAperiodicity() calculates the aperiodicity in multiples of 3 kHz.
// The upper limit is given based on the sampling frequency.
//-----------------------------------------------------------------------------
static void GetCoarseAperiodicity(const double *static_group_delay, int fs,
    int fft_size, int number_of_aperiodicities, const double *window,
    int window_length, const ForwardRealFFT *forward_real_fft,
    double *coarse_aperiodicity) {
  int boundary =
    matlab_round(fft_size * 8.0 / window_length);
  int half_window_length = window_length / 2;

  for (int i = 0; i < fft_size; ++i) forward_real_fft->waveform[i] = 0.0;

  double *power_spectrum = new double[fft_size / 2 + 1];
  int center;
  for (int i = 0; i < number_of_aperiodicities; ++i) {
    center =
      static_cast<int>(world::kFrequencyInterval * (i + 1) * fft_size / fs);
    for (int j = 0; j <= half_window_length * 2; ++j)
      forward_real_fft->waveform[j] =
        static_group_delay[center - half_window_length + j] * window[j];
    fft_execute(forward_real_fft->forward_fft);
    for (int j = 0 ; j <= fft_size / 2; ++j)
      power_spectrum[j] =
        forward_real_fft->spectrum[j][0] * forward_real_fft->spectrum[j][0] +
        forward_real_fft->spectrum[j][1] * forward_real_fft->spectrum[j][1];
    std::sort(power_spectrum, power_spectrum + fft_size / 2 + 1);
    for (int j = 1 ; j <= fft_size / 2; ++j)
      power_spectrum[j] += power_spectrum[j - 1];
    coarse_aperiodicity[i] =
      10 * log10(power_spectrum[fft_size / 2 - boundary - 1] /
                 power_spectrum[fft_size / 2]);
  }
  delete[] power_spectrum;
}

static double D4CLoveTrainSub(const double *x, int fs, int x_length,
    double current_f0, double current_position, int f0_length, int fft_size,
    int boundary0, int boundary1, int boundary2,
    ForwardRealFFT *forward_real_fft) {
  double *power_spectrum = new double[fft_size];

  int window_length = matlab_round(1.5 * fs / current_f0) * 2 + 1;
  GetWindowedWaveform(x, x_length, fs, current_f0, current_position,
    world::kBlackman, 3.0, forward_real_fft->waveform);

  for (int i = window_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  for (int i = 0; i <= boundary0; ++i) power_spectrum[i] = 0.0;
  for (int i = boundary0 + 1; i < fft_size / 2 + 1; ++i)
    power_spectrum[i] =
    forward_real_fft->spectrum[i][0] * forward_real_fft->spectrum[i][0] +
    forward_real_fft->spectrum[i][1] * forward_real_fft->spectrum[i][1];
  for (int i = boundary0; i <= boundary2; ++i)
    power_spectrum[i] += +power_spectrum[i - 1];

  double aperiodicity0 = power_spectrum[boundary1] / power_spectrum[boundary2];
  delete[] power_spectrum;
  return aperiodicity0;
}

//-----------------------------------------------------------------------------
// D4CLoveTrain() determines the aperiodicity with VUV detection.
// If a frame was determined as the unvoiced section, aperiodicity is set to
// very high value as the safeguard.
// If it was voiced section, the aperiodicity of 0 Hz is set to -60 dB.
//-----------------------------------------------------------------------------
static void D4CLoveTrain(const double *x, int fs, int x_length,
    const double *f0, int f0_length, const double *temporal_positions,
    double *aperiodicity0) {
  double lowest_f0 = 40.0;
  int fft_size = static_cast<int>(pow(2.0, 1.0 +
    static_cast<int>(log(3.0 * fs / lowest_f0 + 1) / world::kLog2)));
  ForwardRealFFT forward_real_fft = { 0 };
  InitializeForwardRealFFT(fft_size, &forward_real_fft);

  // Cumulative powers at 100, 4000, 7900 Hz are used for VUV identification.
  int boundary0 = static_cast<int>(ceil(100.0 * fft_size / fs));
  int boundary1 = static_cast<int>(ceil(4000.0 * fft_size / fs));
  int boundary2 = static_cast<int>(ceil(7900.0 * fft_size / fs));
  for (int i = 0; i < f0_length; ++i) {
    if (f0[i] == 0.0) {
      aperiodicity0[i] = 0.0;
      continue;
    }
    aperiodicity0[i] = D4CLoveTrainSub(x, fs, x_length,
      MyMaxDouble(f0[i], lowest_f0), temporal_positions[i], f0_length,
      fft_size, boundary0, boundary1, boundary2, &forward_real_fft);
  }

  DestroyForwardRealFFT(&forward_real_fft);
}

//-----------------------------------------------------------------------------
// D4CGeneralBody() calculates a spectral envelope at a temporal
// position. This function is only used in D4C().
// Caution:
//   forward_fft is allocated in advance to speed up the processing.
//-----------------------------------------------------------------------------
static void D4CGeneralBody(const double *x, int x_length, int fs,
    double current_f0, int fft_size, double current_position,
    int number_of_aperiodicities, const double *window, int window_length,
    const ForwardRealFFT *forward_real_fft, double *coarse_aperiodicity) {
  double *static_centroid = new double[fft_size / 2 + 1];
  double *smoothed_power_spectrum = new double[fft_size / 2 + 1];
  double *static_group_delay = new double[fft_size / 2 + 1];
  GetStaticCentroid(x, x_length, fs, current_f0, fft_size, current_position,
      forward_real_fft, static_centroid);
  GetSmoothedPowerSpectrum(x, x_length, fs, current_f0, fft_size,
      current_position, forward_real_fft, smoothed_power_spectrum);
  GetStaticGroupDelay(static_centroid, smoothed_power_spectrum,
      fs, current_f0, fft_size, static_group_delay);

  GetCoarseAperiodicity(static_group_delay, fs, fft_size,
      number_of_aperiodicities, window, window_length, forward_real_fft,
      coarse_aperiodicity);

  // Revision of the result based on the F0
  for (int i = 0; i < number_of_aperiodicities; ++i)
    coarse_aperiodicity[i] = MyMinDouble(0.0,
      coarse_aperiodicity[i] + (current_f0 - 100) / 50.0);

  delete[] static_centroid;
  delete[] smoothed_power_spectrum;
  delete[] static_group_delay;
}

static void InitializeAperiodicity(int f0_length, int fft_size,
    double **aperiodicity) {
  for (int i = 0; i < f0_length; ++i)
    for (int j = 0; j < fft_size / 2 + 1; ++j)
      aperiodicity[i][j] = 1.0 - world::kMySafeGuardMinimum;
}

static void GetAperiodicity(const double *coarse_frequency_axis,
    const double *coarse_aperiodicity, int number_of_aperiodicities,
    const double *frequency_axis, int fft_size, double *aperiodicity) {
  interp1(coarse_frequency_axis, coarse_aperiodicity,
    number_of_aperiodicities + 2, frequency_axis, fft_size / 2 + 1,
    aperiodicity);
  for (int i = 0; i <= fft_size / 2; ++i)
    aperiodicity[i] = pow(10.0, aperiodicity[i] / 20.0);
}

}  // namespace

void D4C(const double *x, int x_length, int fs,
    const double *temporal_positions, const double *f0, int f0_length,
    int fft_size, const D4COption *option, double **aperiodicity) {
  randn_reseed();

  InitializeAperiodicity(f0_length, fft_size, aperiodicity);

  int fft_size_d4c = static_cast<int>(pow(2.0, 1.0 +
    static_cast<int>(log(4.0 * fs / world::kFloorF0D4C + 1) /
      world::kLog2)));

  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size_d4c, &forward_real_fft);

  int number_of_aperiodicities =
    static_cast<int>(MyMinDouble(world::kUpperLimit, fs / 2.0 -
      world::kFrequencyInterval) / world::kFrequencyInterval);
  // Since the window function is common in D4CGeneralBody(),
  // it is designed here to speed up.
  int window_length =
    static_cast<int>(world::kFrequencyInterval * fft_size_d4c / fs) * 2 + 1;
  double *window =  new double[window_length];
  NuttallWindow(window_length, window);

  // D4C Love Train (Aperiodicity of 0 Hz is given by the different algorithm)
  double *aperiodicity0 = new double[f0_length];
  D4CLoveTrain(x, fs, x_length, f0, f0_length, temporal_positions,
      aperiodicity0);

  double *coarse_aperiodicity = new double[number_of_aperiodicities + 2];
  coarse_aperiodicity[0] = -60.0;
  coarse_aperiodicity[number_of_aperiodicities + 1] =
    -world::kMySafeGuardMinimum;
  double *coarse_frequency_axis = new double[number_of_aperiodicities + 2];
  for (int i = 0; i <= number_of_aperiodicities; ++i)
    coarse_frequency_axis[i] = i * world::kFrequencyInterval;
  coarse_frequency_axis[number_of_aperiodicities + 1] = fs / 2.0;

  double *frequency_axis = new double[fft_size / 2 + 1];
  for (int i = 0; i <= fft_size / 2; ++i)
    frequency_axis[i] = static_cast<double>(i) * fs / fft_size;

  for (int i = 0; i < f0_length; ++i) {
    if (f0[i] == 0 || aperiodicity0[i] <= option->threshold) continue;
    D4CGeneralBody(x, x_length, fs, MyMaxDouble(world::kFloorF0D4C, f0[i]),
        fft_size_d4c, temporal_positions[i], number_of_aperiodicities, window,
        window_length, &forward_real_fft, &coarse_aperiodicity[1]);

    // Linear interpolation to convert the coarse aperiodicity into its
    // spectral representation.
    GetAperiodicity(coarse_frequency_axis, coarse_aperiodicity,
        number_of_aperiodicities, frequency_axis, fft_size, aperiodicity[i]);
  }

  DestroyForwardRealFFT(&forward_real_fft);
  delete[] aperiodicity0;
  delete[] coarse_frequency_axis;
  delete[] coarse_aperiodicity;
  delete[] window;
  delete[] frequency_axis;
}

void InitializeD4COption(D4COption *option) {
  option->threshold = world::kThreshold;
}
