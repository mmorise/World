//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//
// F0 estimation based on instantaneous frequency.
// This method is carried out by using the output of Dio().
//-----------------------------------------------------------------------------
#include "world/stonemask.h"

#include <math.h>

#include "world/common.h"
#include "world/constantnumbers.h"
#include "world/fft.h"
#include "world/matlabfunctions.h"

namespace {
//-----------------------------------------------------------------------------
// GetBaseIndex() calculates the temporal positions for windowing.
// Since the result includes negative value and the value that exceeds the
// length of the input signal, it must be modified appropriately.
//-----------------------------------------------------------------------------
static void GetBaseIndex(double current_position, const double *base_time,
    int base_time_length, int fs, int *index_raw) {
  for (int i = 0; i < base_time_length; ++i)
    index_raw[i] = matlab_round((current_position + base_time[i]) * fs);
}

//-----------------------------------------------------------------------------
// GetMainWindow() generates the window function.
//-----------------------------------------------------------------------------
static void GetMainWindow(double current_position, const int *index_raw,
    int base_time_length, int fs, double window_length_in_time,
    double *main_window) {
  double tmp = 0.0;
  for (int i = 0; i < base_time_length; ++i) {
    tmp = (index_raw[i] - 1.0) / fs - current_position;
    main_window[i] = 0.42 +
      0.5 * cos(2.0 * world::kPi * tmp / window_length_in_time) +
      0.08 * cos(4.0 * world::kPi * tmp / window_length_in_time);
  }
}

//-----------------------------------------------------------------------------
// GetDiffWindow() generates the differentiated window.
// Diff means differential.
//-----------------------------------------------------------------------------
static void GetDiffWindow(const double *main_window, int base_time_length,
    double *diff_window) {
  diff_window[0] = -main_window[1] / 2.0;
  for (int i = 1; i < base_time_length - 1; ++i)
    diff_window[i] = -(main_window[i + 1] - main_window[i - 1]) / 2.0;
  diff_window[base_time_length - 1] = main_window[base_time_length - 2] / 2.0;
}

//-----------------------------------------------------------------------------
// GetSpectra() calculates two spectra of the waveform windowed by windows
// (main window and diff window).
//-----------------------------------------------------------------------------
static void GetSpectra(const double *x, int x_length, int fft_size,
    const int *index_raw, const double *main_window, const double *diff_window,
    int base_time_length, const ForwardRealFFT *forward_real_fft,
    fft_complex *main_spectrum, fft_complex *diff_spectrum) {
  int *index = new int[base_time_length];

  for (int i = 0; i < base_time_length; ++i)
    index[i] = MyMaxInt(0, MyMinInt(x_length - 1, index_raw[i] - 1));
  for (int i = 0; i < base_time_length; ++i)
    forward_real_fft->waveform[i] = x[index[i]] * main_window[i];
  for (int i = base_time_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;

  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    main_spectrum[i][0] = forward_real_fft->spectrum[i][0];
    main_spectrum[i][1] = forward_real_fft->spectrum[i][1];
  }

  for (int i = 0; i < base_time_length; ++i)
    forward_real_fft->waveform[i] = x[index[i]] * diff_window[i];
  for (int i = base_time_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    diff_spectrum[i][0] = forward_real_fft->spectrum[i][0];
    diff_spectrum[i][1] = forward_real_fft->spectrum[i][1];
  }

  delete[] index;
}

//-----------------------------------------------------------------------------
// FixF0() fixed the F0 by instantaneous frequency.
//-----------------------------------------------------------------------------
static double FixF0(const double *power_spectrum, const double *numerator_i,
    int fft_size, int fs, double initial_f0, int number_of_harmonics) {
  double *amplitude_list = new double[number_of_harmonics];
  double *instantaneous_frequency_list = new double[number_of_harmonics];
  int index;
  for (int i = 0; i < number_of_harmonics; ++i) {
    index = MyMinInt(matlab_round(initial_f0 * fft_size / fs * (i + 1)), 
            fft_size / 2);
    instantaneous_frequency_list[i] = power_spectrum[index] == 0.0 ? 0.0 :
      static_cast<double>(index) * fs / fft_size +
      numerator_i[index] / power_spectrum[index] * fs / 2.0 / world::kPi;
    amplitude_list[i] = sqrt(power_spectrum[index]);
  }
  double denominator = 0.0;
  double numerator = 0.0;
  for (int i = 0; i < number_of_harmonics; ++i) {
    numerator += amplitude_list[i] * instantaneous_frequency_list[i];
    denominator += amplitude_list[i] * (i + 1);
  }
  delete[] amplitude_list;
  delete[] instantaneous_frequency_list;
  return numerator / (denominator + world::kMySafeGuardMinimum);
}

//-----------------------------------------------------------------------------
// GetTentativeF0() calculates the F0 based on the instantaneous frequency.
//-----------------------------------------------------------------------------
static double GetTentativeF0(const double *power_spectrum,
    const double *numerator_i, int fft_size, int fs, double initial_f0) {
  double tentative_f0 =
    FixF0(power_spectrum, numerator_i, fft_size, fs, initial_f0, 2);

  // If the fixed value is too large, the result will be rejected.
  if (tentative_f0 <= 0.0 || tentative_f0 > initial_f0 * 2) return 0.0;

  return FixF0(power_spectrum, numerator_i, fft_size, fs, tentative_f0, 6);
}

//-----------------------------------------------------------------------------
// GetMeanF0() calculates the instantaneous frequency.
//-----------------------------------------------------------------------------
static double GetMeanF0(const double *x, int x_length, int fs,
    double current_position, double initial_f0, int fft_size,
    double window_length_in_time, const double *base_time,
    int base_time_length) {
  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);
  fft_complex *main_spectrum = new fft_complex[fft_size];
  fft_complex *diff_spectrum = new fft_complex[fft_size];

  int *index_raw = new int[base_time_length];
  double *main_window = new double[base_time_length];
  double *diff_window = new double[base_time_length];

  GetBaseIndex(current_position, base_time, base_time_length, fs, index_raw);
  GetMainWindow(current_position, index_raw, base_time_length, fs,
      window_length_in_time, main_window);
  GetDiffWindow(main_window, base_time_length, diff_window);
  GetSpectra(x, x_length, fft_size, index_raw, main_window, diff_window,
      base_time_length, &forward_real_fft, main_spectrum, diff_spectrum);

  double *power_spectrum = new double[fft_size / 2 + 1];
  double *numerator_i = new double[fft_size / 2 + 1];
  for (int j = 0; j <= fft_size / 2; ++j) {
    numerator_i[j] = main_spectrum[j][0] * diff_spectrum[j][1] -
      main_spectrum[j][1] * diff_spectrum[j][0];
    power_spectrum[j] = main_spectrum[j][0] * main_spectrum[j][0] +
      main_spectrum[j][1] * main_spectrum[j][1];
  }

  double tentative_f0 = GetTentativeF0(power_spectrum, numerator_i,
    fft_size, fs, initial_f0);

  delete[] diff_spectrum;
  delete[] diff_window;
  delete[] main_window;
  delete[] index_raw;
  delete[] numerator_i;
  delete[] power_spectrum;
  delete[] main_spectrum;
  DestroyForwardRealFFT(&forward_real_fft);

  return tentative_f0;
}

//-----------------------------------------------------------------------------
// GetRefinedF0() fixes the F0 estimated by Dio(). This function uses
// instantaneous frequency.
//-----------------------------------------------------------------------------
static double GetRefinedF0(const double *x, int x_length, int fs,
    double current_potision, double initial_f0) {
  if (initial_f0 <= world::kFloorF0StoneMask || initial_f0 > fs / 12.0)
    return 0.0;

  int half_window_length = static_cast<int>(1.5 * fs / initial_f0 + 1.0);
  double window_length_in_time = (2.0 * half_window_length + 1.0) / fs;
  double *base_time = new double[half_window_length * 2 + 1];
  for (int i = 0; i < half_window_length * 2 + 1; i++)
    base_time[i] = static_cast<double>(-half_window_length + i) / fs;
  int fft_size = static_cast<int>(pow(2.0, 2.0 +
    static_cast<int>(log(half_window_length * 2.0 + 1.0) / world::kLog2)));

  double mean_f0 = GetMeanF0(x, x_length, fs, current_potision,
      initial_f0, fft_size, window_length_in_time, base_time,
      half_window_length * 2 + 1);

  // If amount of correction is overlarge (20 %), initial F0 is employed.
  if (fabs(mean_f0 - initial_f0) > initial_f0 * 0.2) mean_f0 = initial_f0;

  delete[] base_time;

  return mean_f0;
}

}  // namespace

void StoneMask(const double *x, int x_length, int fs,
    const double *temporal_positions, const double *f0, int f0_length,
    double *refined_f0) {
  for (int i = 0; i < f0_length; i++)
    refined_f0[i] =
      GetRefinedF0(x, x_length, fs, temporal_positions[i], f0[i]);
}
