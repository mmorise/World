//-----------------------------------------------------------------------------
// Copyright 2012-2015 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// F0 estimation based on instantaneous frequency.
// This method is carried out by using the output of Dio().
//-----------------------------------------------------------------------------
#include "./stonemask.h"

#include <math.h>

#include "./common.h"
#include "./constantnumbers.h"
#include "./fft.h"
#include "./matlabfunctions.h"

namespace {

//-----------------------------------------------------------------------------
// GetIndexRaw() calculates the temporal positions for windowing.
// Since the result includes negative value and the value that exceeds the
// length of the input signal, it must be modified appropriately.
//-----------------------------------------------------------------------------
void GetIndexRaw(double current_time, double *base_time, int base_time_length,
    int fs, int *index_raw) {
  for (int i = 0; i < base_time_length; ++i)
    index_raw[i] = matlab_round((current_time + base_time[i]) * fs);
}

//-----------------------------------------------------------------------------
// GetMainWindow() generates the window function.
//-----------------------------------------------------------------------------
void GetMainWindow(double current_time, int *index_raw,
    int base_time_length, int fs, double window_length_in_time,
    double *main_window) {
  double tmp = 0.0;
  for (int i = 0; i < base_time_length; ++i) {
    tmp = static_cast<double>(index_raw[i] - 1.0) / fs - current_time;
    main_window[i] = 0.42 +
      0.5 * cos(2.0 * world::kPi * tmp / window_length_in_time) +
      0.08 * cos(4.0 * world::kPi * tmp / window_length_in_time);
  }
}

//-----------------------------------------------------------------------------
// GetDiffWindow() generates the differentiated window.
// Diff means differential.
//-----------------------------------------------------------------------------
void GetDiffWindow(double *main_window, int base_time_length,
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
void GetSpectra(double *x, int x_length, int fft_size, int *index_raw,
    double *main_window, double *diff_window, int base_time_length,
    ForwardRealFFT *forward_real_fft, fft_complex *main_spectrum,
    fft_complex *diff_spectrum) {
  int *index = new int[base_time_length];

  for (int i = 0; i < base_time_length; ++i)
    index[i] = MyMax(0, MyMin(x_length - 1, index_raw[i] - 1));
  for (int i = 0; i < base_time_length; ++i)
    forward_real_fft->waveform[i] = x[index[i]] * main_window[i];
  for (int i = base_time_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;

  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    main_spectrum[i][0] = forward_real_fft->spectrum[i][0];
    main_spectrum[i][1] = -forward_real_fft->spectrum[i][1];
  }

  for (int i = 0; i < base_time_length; ++i)
    forward_real_fft->waveform[i] = x[index[i]] * diff_window[i];
  for (int i = base_time_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    diff_spectrum[i][0] = forward_real_fft->spectrum[i][0];
    diff_spectrum[i][1] = -forward_real_fft->spectrum[i][1];
  }

  delete[] index;
}

//-----------------------------------------------------------------------------
// GetTentativeF0() calculates the F0 based on the instantaneous frequency.
// Calculated value is tentative because it is fixed as needed.
//-----------------------------------------------------------------------------
double GetTentativeF0(double *power_spectrum, double *numerator_i,
    int fft_size, int fs, double f0_initial) {
  double power_list[6];
  double fixp_list[6];
  int index;
  for (int i = 0; i < 2; ++i) {
    index = matlab_round(f0_initial * fft_size / fs * (i + 1));

    fixp_list[i] = static_cast<double>(index) * fs / fft_size +
      numerator_i[index] / power_spectrum[index] * fs / 2.0 / world::kPi;
    power_list[i] = power_spectrum[index];
  }

  double tmp1 = 0.0;
  double tmp2 = 0.0;
  for (int i = 0; i < 2; ++i) {
    tmp1 += power_list[i] * fixp_list[i] / (i + 1);
    tmp2 += power_list[i];
  }

  // block division by zero which will cause segmentatin fault.
  if (tmp2 == 0) return 0;
  f0_initial = tmp1 / tmp2;

  for (int i = 0; i < 6; ++i) {
    index = matlab_round(f0_initial * fft_size / fs * (i + 1));
    fixp_list[i] = static_cast<double>(index) * fs / fft_size +
      numerator_i[index] / power_spectrum[index] * fs / 2.0 / world::kPi;
    power_list[i] = sqrt(power_spectrum[index]);
  }
  // A minor bug was fixed (2015/12/02)
  tmp1 = 0.0;
  tmp2 = 0.0;
  for (int i = 0; i < 6; ++i) {
    tmp1 += power_list[i] * fixp_list[i];
    tmp2 += power_list[i] * (i + 1);
  }

  return tmp1 / tmp2;
}

//-----------------------------------------------------------------------------
// GetMeanF0() calculates the instantaneous frequency.
//-----------------------------------------------------------------------------
double GetMeanF0(double *x, int x_length, int fs, double current_time,
    double f0_initial, int fft_size, double window_length_in_time,
    double *base_time, int base_time_length) {
  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);
  fft_complex *main_spectrum = new fft_complex[fft_size];
  fft_complex *diff_spectrum = new fft_complex[fft_size];

  int *index_raw = new int[base_time_length];
  double *main_window = new double[base_time_length];
  double *diff_window = new double[base_time_length];

  GetIndexRaw(current_time, base_time, base_time_length, fs, index_raw);
  GetMainWindow(current_time, index_raw, base_time_length, fs,
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
      fft_size, fs, f0_initial);

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
double GetRefinedF0(double *x, int x_length, int fs, double current_time,
    double current_f0) {
  // A safeguard was added (2015/12/02).
  if (current_f0 <= 0.0 || current_f0 > fs / 12.0)
    return 0.0;

  double f0_initial = current_f0;  // bug fix 2015/11/29
  int half_window_length = static_cast<int>(3.0 * static_cast<double>(fs)
    / f0_initial / 2.0 + 1.0);
  double window_length_in_time = (2.0 *
    static_cast<double>(half_window_length) + 1) /
    static_cast<double>(fs);
  double *base_time = new double[half_window_length * 2 + 1];
  for (int i = 0; i < half_window_length * 2 + 1; i++) {
    base_time[i] = static_cast<double>(-half_window_length + i) / fs;
  }
  int fft_size = static_cast<int>(pow(2.0, 2.0 +
    static_cast<int>(log(half_window_length * 2.0 + 1.0) / world::kLog2)));

  double mean_f0 = GetMeanF0(x, x_length, fs, current_time,
      f0_initial, fft_size, window_length_in_time, base_time,
      half_window_length * 2 + 1);

  // If amount of correction is overlarge (20 %), initial F0 is employed.
  if (fabs(mean_f0 - f0_initial) / f0_initial > 0.2) mean_f0 = f0_initial;

  delete[] base_time;

  return mean_f0;
}

}  // namespace

void StoneMask(double *x, int x_length, int fs, double *time_axis, double *f0,
    int f0_length, double *refined_f0) {
  for (int i = 0; i < f0_length; i++)
    refined_f0[i] = GetRefinedF0(x, x_length, fs, time_axis[i], f0[i]);
}
