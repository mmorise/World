//-----------------------------------------------------------------------------
// Copyright 2012-2016 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// F0 estimation based on Harvest.
//-----------------------------------------------------------------------------
#include "world/harvest.h"

#include <math.h>
#include <cstdlib>

#include "world/common.h"
#include "world/constantnumbers.h"
#include "world/fft.h"
#include "world/matlabfunctions.h"

//-----------------------------------------------------------------------------
// struct for RawEventByHarvest()
// "negative" means "zero-crossing point going from positive to negative"
// "positive" means "zero-crossing point going from negative to positive"
//-----------------------------------------------------------------------------
typedef struct {
  double *negative_interval_locations;
  double *negative_intervals;
  int number_of_negatives;
  double *positive_interval_locations;
  double *positive_intervals;
  int number_of_positives;
  double *peak_interval_locations;
  double *peak_intervals;
  int number_of_peaks;
  double *dip_interval_locations;
  double *dip_intervals;
  int number_of_dips;
} ZeroCrossings;

namespace {
//-----------------------------------------------------------------------------
// Since the waveform of beginning and ending after decimate include noise,
// the input waveform is extended. This is the processing for the
// compatibility with MATLAB version.
//-----------------------------------------------------------------------------
static void GetWaveformAndSpectrumSub(const double *x, int x_length,
    int y_length, double actual_fs, int decimation_ratio, double *y) {
  if (decimation_ratio == 1) {
    for (int i = 0; i < x_length; ++i) y[i] = x[i];
    return;
  }

  int lag =
    static_cast<int>(ceil(140.0 / decimation_ratio) * decimation_ratio);
  int new_x_length = x_length + lag * 2;
  double *new_y = new double[new_x_length];
  for (int i = 0; i < new_x_length; ++i) new_y[i] = 0.0;
  double *new_x = new double[new_x_length];
  for (int i = 0; i < lag; ++i) new_x[i] = x[0];
  for (int i = lag; i < lag + x_length; ++i) new_x[i] = x[i - lag];
  for (int i = lag + x_length; i < new_x_length; ++i) {
    new_x[i] = x[x_length - 1];
  }

  decimate(new_x, new_x_length, decimation_ratio, new_y);
  for (int i = 0; i < y_length; ++i) y[i] = new_y[lag / decimation_ratio + i];

  delete[] new_x;
  delete[] new_y;
}

//-----------------------------------------------------------------------------
// GetDownsampledSignal() calculates the spectrum for estimation.
// This function carries out downsampling to speed up the estimation process
// and calculates the spectrum of the downsampled signal.
//-----------------------------------------------------------------------------
static void GetWaveformAndSpectrum(const double *x, int x_length,
  int y_length, double actual_fs, int fft_size, int decimation_ratio,
  double *y, fft_complex *y_spectrum) {
  // Initialization
  for (int i = 0; i < fft_size; ++i) y[i] = 0.0;

  // Processing for the compatibility with MATLAB version
  GetWaveformAndSpectrumSub(x, x_length, y_length, actual_fs,
      decimation_ratio, y);

  // Removal of the DC component (y = y - mean value of y)
  double mean_y = 0.0;
  for (int i = 0; i < y_length; ++i) mean_y += y[i];
  mean_y /= y_length;
  for (int i = 0; i < y_length; ++i) y[i] -= mean_y;
  for (int i = y_length; i < fft_size; ++i) y[i] = 0.0;

  fft_plan forwardFFT =
    fft_plan_dft_r2c_1d(fft_size, y, y_spectrum, FFT_ESTIMATE);
  fft_execute(forwardFFT);

  fft_destroy_plan(forwardFFT);
}

//-----------------------------------------------------------------------------
// GetFilteredSignal() calculates the signal that is the convolution of the
// input signal and low-pass filter.
// This function is only used in RawEventByDio()
//-----------------------------------------------------------------------------
static void GetFilteredSignal(double boundary_f0, int fft_size, double fs,
    const fft_complex *y_spectrum, int y_length, double *filtered_signal) {
  int filter_length_half = matlab_round(fs / boundary_f0 * 2.0);
  double *low_pass_filter = new double[fft_size];
  NuttallWindow(filter_length_half * 2 + 1, low_pass_filter);
  for (int i = -filter_length_half; i <= filter_length_half; ++i)
    low_pass_filter[i + filter_length_half] *=
      cos(2 * world::kPi * boundary_f0 * i / fs);
  for (int i = filter_length_half * 2 + 1; i < fft_size; ++i)
    low_pass_filter[i] = 0.0;

  fft_complex *low_pass_filter_spectrum = new fft_complex[fft_size];
  fft_plan forwardFFT = fft_plan_dft_r2c_1d(fft_size, low_pass_filter,
      low_pass_filter_spectrum, FFT_ESTIMATE);
  fft_execute(forwardFFT);

  // Convolution
  double tmp = y_spectrum[0][0] * low_pass_filter_spectrum[0][0] -
    y_spectrum[0][1] * low_pass_filter_spectrum[0][1];
  low_pass_filter_spectrum[0][1] =
    y_spectrum[0][0] * low_pass_filter_spectrum[0][1] +
    y_spectrum[0][1] * low_pass_filter_spectrum[0][0];
  low_pass_filter_spectrum[0][0] = tmp;
  for (int i = 1; i <= fft_size / 2; ++i) {
    tmp = y_spectrum[i][0] * low_pass_filter_spectrum[i][0] -
      y_spectrum[i][1] * low_pass_filter_spectrum[i][1];
    low_pass_filter_spectrum[i][1] =
      y_spectrum[i][0] * low_pass_filter_spectrum[i][1] +
      y_spectrum[i][1] * low_pass_filter_spectrum[i][0];
    low_pass_filter_spectrum[i][0] = tmp;
    low_pass_filter_spectrum[fft_size - i - 1][0] =
      low_pass_filter_spectrum[i][0];
    low_pass_filter_spectrum[fft_size - i - 1][1] =
      low_pass_filter_spectrum[i][1];
  }

  fft_plan inverseFFT = fft_plan_dft_c2r_1d(fft_size,
      low_pass_filter_spectrum, filtered_signal, FFT_ESTIMATE);
  fft_execute(inverseFFT);

  // Compensation of the delay.
  int index_bias = filter_length_half + 1;
  for (int i = 0; i < y_length; ++i)
    filtered_signal[i] = filtered_signal[i + index_bias];

  fft_destroy_plan(inverseFFT);
  fft_destroy_plan(forwardFFT);
  delete[] low_pass_filter_spectrum;
  delete[] low_pass_filter;
}

//-----------------------------------------------------------------------------
// CheckEvent() returns 1, provided that the input value is over 1.
// This function is for RawEventByDio().
//-----------------------------------------------------------------------------
static inline int CheckEvent(int x) {
  return x > 0 ? 1 : 0;
}

//-----------------------------------------------------------------------------
// ZeroCrossingEngine() calculates the zero crossing points from positive to
// negative.
//-----------------------------------------------------------------------------
static int ZeroCrossingEngine(const double *filtered_signal, int y_length,
    double fs, double *interval_locations, double *intervals) {
  int *negative_going_points = new int[y_length];

  for (int i = 0; i < y_length - 1; ++i)
    negative_going_points[i] =
      0.0 < filtered_signal[i] && filtered_signal[i + 1] <= 0.0 ? i + 1 : 0;
  negative_going_points[y_length - 1] = 0;

  int *edges = new int[y_length];
  int count = 0;
  for (int i = 0; i < y_length; ++i)
    if (negative_going_points[i] > 0)
      edges[count++] = negative_going_points[i];

  if (count < 2) {
    delete[] edges;
    delete[] negative_going_points;
    return 0;
  }

  double *fine_edges = new double[count];
  for (int i = 0; i < count; ++i)
    fine_edges[i] =
      edges[i] - filtered_signal[edges[i] - 1] /
      (filtered_signal[edges[i]] - filtered_signal[edges[i] - 1]);

  for (int i = 0; i < count - 1; ++i) {
    intervals[i] = fs / (fine_edges[i + 1] - fine_edges[i]);
    interval_locations[i] = (fine_edges[i] + fine_edges[i + 1]) / 2.0 / fs;
  }

  delete[] fine_edges;
  delete[] edges;
  delete[] negative_going_points;
  return count - 1;
}

//-----------------------------------------------------------------------------
// GetFourZeroCrossingIntervals() calculates four zero-crossing intervals.
// (1) Zero-crossing going from negative to positive.
// (2) Zero-crossing going from positive to negative.
// (3) Peak, and (4) dip. (3) and (4) are calculated from the zero-crossings of
// the differential of waveform.
//-----------------------------------------------------------------------------
static void GetFourZeroCrossingIntervals(double *filtered_signal, int y_length,
    double actual_fs, ZeroCrossings *zero_crossings) {
  const int kMaximumNumber = y_length;
  zero_crossings->negative_interval_locations = new double[kMaximumNumber];
  zero_crossings->positive_interval_locations = new double[kMaximumNumber];
  zero_crossings->peak_interval_locations = new double[kMaximumNumber];
  zero_crossings->dip_interval_locations = new double[kMaximumNumber];
  zero_crossings->negative_intervals = new double[kMaximumNumber];
  zero_crossings->positive_intervals = new double[kMaximumNumber];
  zero_crossings->peak_intervals = new double[kMaximumNumber];
  zero_crossings->dip_intervals = new double[kMaximumNumber];

  zero_crossings->number_of_negatives = ZeroCrossingEngine(filtered_signal,
      y_length, actual_fs, zero_crossings->negative_interval_locations,
      zero_crossings->negative_intervals);

  for (int i = 0; i < y_length; ++i) filtered_signal[i] = -filtered_signal[i];
  zero_crossings->number_of_positives = ZeroCrossingEngine(filtered_signal,
      y_length, actual_fs, zero_crossings->positive_interval_locations,
      zero_crossings->positive_intervals);

  for (int i = 0; i < y_length - 1; ++i) filtered_signal[i] =
    filtered_signal[i] - filtered_signal[i + 1];
  zero_crossings->number_of_peaks = ZeroCrossingEngine(filtered_signal,
      y_length - 1, actual_fs, zero_crossings->peak_interval_locations,
      zero_crossings->peak_intervals);

  for (int i = 0; i < y_length - 1; ++i)
    filtered_signal[i] = -filtered_signal[i];
  zero_crossings->number_of_dips = ZeroCrossingEngine(filtered_signal,
      y_length - 1, actual_fs, zero_crossings->dip_interval_locations,
      zero_crossings->dip_intervals);
}

//-----------------------------------------------------------------------------
// GetF0CandidatesSub() calculates the f0 candidates and deviations.
// This is the sub-function of GetF0Candidates() and assumes the calculation.
//-----------------------------------------------------------------------------
static void GetF0CandidatesSub(double **const interpolated_f0_set,
    int time_axis_length, double f0_floor, double f0_ceil, double boundary_f0,
    double *f0_candidates) {
  for (int i = 0; i < time_axis_length; ++i) {
    f0_candidates[i] = (interpolated_f0_set[0][i] +
      interpolated_f0_set[1][i] + interpolated_f0_set[2][i] +
      interpolated_f0_set[3][i]) / 4.0;

    // 1.1 and 0.9 are optimized values.
    if (f0_candidates[i] > boundary_f0 * 1.1 ||
        f0_candidates[i] < boundary_f0 * 0.9 ||
        f0_candidates[i] > f0_ceil || f0_candidates[i] < f0_floor) {
      f0_candidates[i] = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
// GetF0Candidates() calculates the F0 candidates based on the zero-crossings.
// Calculation of F0 candidates is carried out in GetF0CandidatesSub().
//-----------------------------------------------------------------------------
static void GetF0Candidates(const ZeroCrossings *zero_crossings,
    double boundary_f0, double f0_floor, double f0_ceil,
    const double *time_axis, int time_axis_length, double *f0_candidates) {
  if (0 == CheckEvent(zero_crossings->number_of_negatives - 2) *
      CheckEvent(zero_crossings->number_of_positives - 2) *
      CheckEvent(zero_crossings->number_of_peaks - 2) *
      CheckEvent(zero_crossings->number_of_dips - 2)) {
    for (int i = 0; i < time_axis_length; ++i) {
      f0_candidates[i] = 0.0;
    }
    return;
  }

  double *interpolated_f0_set[4];
  for (int i = 0; i < 4; ++i)
    interpolated_f0_set[i] = new double[time_axis_length];

  interp1(zero_crossings->negative_interval_locations,
      zero_crossings->negative_intervals,
      zero_crossings->number_of_negatives,
      time_axis, time_axis_length, interpolated_f0_set[0]);
  interp1(zero_crossings->positive_interval_locations,
      zero_crossings->positive_intervals,
      zero_crossings->number_of_positives,
      time_axis, time_axis_length, interpolated_f0_set[1]);
  interp1(zero_crossings->peak_interval_locations,
      zero_crossings->peak_intervals, zero_crossings->number_of_peaks,
      time_axis, time_axis_length, interpolated_f0_set[2]);
  interp1(zero_crossings->dip_interval_locations,
      zero_crossings->dip_intervals, zero_crossings->number_of_dips,
      time_axis, time_axis_length, interpolated_f0_set[3]);

  GetF0CandidatesSub(interpolated_f0_set, time_axis_length, f0_floor,
      f0_ceil, boundary_f0, f0_candidates);
  for (int i = 0; i < 4; ++i) delete[] interpolated_f0_set[i];
}

//-----------------------------------------------------------------------------
// DestroyZeroCrossings() frees the memory of array in the struct
//-----------------------------------------------------------------------------
static void DestroyZeroCrossings(ZeroCrossings *zero_crossings) {
  delete[] zero_crossings->negative_interval_locations;
  delete[] zero_crossings->positive_interval_locations;
  delete[] zero_crossings->peak_interval_locations;
  delete[] zero_crossings->dip_interval_locations;
  delete[] zero_crossings->negative_intervals;
  delete[] zero_crossings->positive_intervals;
  delete[] zero_crossings->peak_intervals;
  delete[] zero_crossings->dip_intervals;
}

//-----------------------------------------------------------------------------
// CalculateRawEvent() calculates the zero-crossings.
//-----------------------------------------------------------------------------
static void CalculateRawEvent(double boundary_f0, double fs,
  const fft_complex *y_spectrum, int y_length, int fft_size, double f0_floor,
  double f0_ceil, const double *time_axis, int time_axis_length,
  double *f0_candidates) {
  double *filtered_signal = new double[fft_size];
  GetFilteredSignal(boundary_f0, fft_size, fs, y_spectrum,
    y_length, filtered_signal);

  ZeroCrossings zero_crossings = { 0 };
  GetFourZeroCrossingIntervals(filtered_signal, y_length, fs,
    &zero_crossings);

  GetF0Candidates(&zero_crossings, boundary_f0, f0_floor, f0_ceil,
    time_axis, time_axis_length, f0_candidates);

  DestroyZeroCrossings(&zero_crossings);
  delete[] filtered_signal;
}

//-----------------------------------------------------------------------------
// GetF0CandidateAndStabilityMap() calculates all f0 candidates and
// their stabilities.
//-----------------------------------------------------------------------------
static void GetF0CandidateMap(double *boundary_f0_list, int number_of_bands,
  double actual_fs, int y_length, double *time_axis, int f0_length,
  fft_complex *y_spectrum, int fft_size, double f0_floor, double f0_ceil,
  double **f0_candidate_map) {
  // Calculation of the acoustics events (zero-crossing)
  for (int i = 0; i < number_of_bands; ++i) {
    CalculateRawEvent(boundary_f0_list[i], actual_fs, y_spectrum, y_length,
      fft_size, f0_floor, f0_ceil, time_axis, f0_length, f0_candidate_map[i]);
  }
}

//-----------------------------------------------------------------------------
// DetectF0CandidatesSub1() calculates VUV areas.
//-----------------------------------------------------------------------------
static int DetectF0CandidatesSub1(const int *vuv, int number_of_channels,
  int *st, int *ed) {
  int number_of_voiced_sections = 0;
  int tmp;
  for (int i = 1; i < number_of_channels; ++i) {
    tmp = vuv[i] - vuv[i - 1];
    if (tmp == 1) st[number_of_voiced_sections] = i;
    if (tmp == -1) ed[number_of_voiced_sections++] = i;
  }

  return number_of_voiced_sections;
}

//-----------------------------------------------------------------------------
// DetectF0CandidatesSub2() calculates F0 candidates in a frame
//-----------------------------------------------------------------------------
static int DetectF0CandidatesSub2(const int *vuv,
    double **const raw_f0_candidates, int index, int number_of_voiced_sections,
    const int *st, const int *ed, int max_candidates, double *f0_list) {
  int number_of_candidates = 0;
  double tmp_f0;
  for (int i = 0; i < number_of_voiced_sections; ++i) {
    if (ed[i] - st[i] < 10) continue;

    tmp_f0 = 0.0;
    for (int j = st[i]; j < ed[i]; ++j) {
      tmp_f0 += raw_f0_candidates[j][index];
    }
    tmp_f0 /= (ed[i] - st[i]);
    f0_list[number_of_candidates++] = tmp_f0;
  }

  for (int i = number_of_candidates; i < max_candidates; ++i) {
    f0_list[i] = 0.0;
  }
  return number_of_candidates;
}

//-----------------------------------------------------------------------------
// DetectF0Candidates() detectes F0 candidates from multi-channel candidates.
//-----------------------------------------------------------------------------
static int DetectF0Candidates(double **const raw_f0_candidates,
    int number_of_channels, int f0_length, int max_candidates,
    double **f0_candidates) {
  int number_of_candidates = 0;

  int *vuv = new int[number_of_channels];
  int *st = new int[number_of_channels];
  int *ed = new int[number_of_channels];
  int number_of_voiced_sections;
  for (int i = 0; i < f0_length; ++i) {
    for (int j = 0; j < number_of_channels; ++j) {
      vuv[j] = raw_f0_candidates[j][i] > 0 ? 1 : 0;
    }
    vuv[0] = vuv[number_of_channels - 1] = 0;
    number_of_voiced_sections = DetectF0CandidatesSub1(vuv, number_of_channels,
      st, ed);
    number_of_candidates =
      MyMaxInt(number_of_candidates, DetectF0CandidatesSub2(vuv,
          raw_f0_candidates, i, number_of_voiced_sections, st, ed,
          max_candidates, f0_candidates[i]));
  }

  delete[] vuv;
  delete[] st;
  delete[] ed;
  return number_of_candidates;
}

//-----------------------------------------------------------------------------
// OverlapF0Candidates() spreads the candidates to anteroposterior frames.
//-----------------------------------------------------------------------------
static void OverlapF0Candidates(int f0_length, int number_of_candidates,
    double **f0_candidates) {
  const int n = 3;
  for (int i = 1; i <= n; ++i) {
    for (int j = 0; j < number_of_candidates; ++j) {
      for (int k = i; k < f0_length; ++k) {
        f0_candidates[k][j + (number_of_candidates * i)] =
          f0_candidates[k - i][j];
      }
      for (int k = 0; k < f0_length - i; ++k) {
        f0_candidates[k][j + (number_of_candidates * (i + n))] =
          f0_candidates[k + i][j];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// GetIndexRaw() calculates the temporal positions for windowing.
// Since the result includes negative value and the value that exceeds the
// length of the input signal, it must be modified appropriately.
//-----------------------------------------------------------------------------
static void GetIndexRaw(double current_time, const double *base_time,
    int base_time_length, double fs, int *index_raw) {
  // First-aid treatment
  int basic_index = matlab_round((current_time + base_time[0]) * fs + 0.001);

  for (int i = 0; i < base_time_length; ++i) index_raw[i] = basic_index + i;
}

//-----------------------------------------------------------------------------
// GetMainWindow() generates the window function.
//-----------------------------------------------------------------------------
static void GetMainWindow(double current_time, const int *index_raw,
    int base_time_length, double fs, double window_length_in_time,
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
// FixF0() fixed the F0 by instantaneous frequency.
//-----------------------------------------------------------------------------
static void FixF0(const double *power_spectrum, const double *numerator_i,
  int fft_size, double fs, double current_f0, int number_of_harmonics,
  double *refined_f0, double *score) {
  double *amp_list = new double[number_of_harmonics];
  double *fixp_list = new double[number_of_harmonics];

  int index;
  for (int i = 0; i < number_of_harmonics; ++i) {
    index = matlab_round(current_f0 * fft_size / fs * (i + 1));
    fixp_list[i] = power_spectrum[index] == 0.0 ? 0.0 :
      static_cast<double>(index) * fs / fft_size +
      numerator_i[index] / power_spectrum[index] * fs / 2.0 / world::kPi;
    amp_list[i] = sqrt(power_spectrum[index]);
  }
  double denominator = 0.0;
  double numerator = 0.0;
  *score = 0.0;
  for (int i = 0; i < number_of_harmonics; ++i) {
    numerator += amp_list[i] * fixp_list[i];
    denominator += amp_list[i] * (i + 1.0);
    *score += fabs((fixp_list[i] / (i + 1.0) - current_f0) / current_f0);
  }

  *refined_f0 = numerator / (denominator + world::kMySafeGuardMinimum);
  *score = 1.0 / (*score / number_of_harmonics + world::kMySafeGuardMinimum);

  delete[] amp_list;
  delete[] fixp_list;
}

//-----------------------------------------------------------------------------
// GetMeanF0() calculates the instantaneous frequency.
//-----------------------------------------------------------------------------
static void GetMeanF0(const double *x, int x_length, double fs,
    double current_time, double current_f0, int fft_size,
    double window_length_in_time, const double *base_time,
    int base_time_length, double *refined_f0, double *score) {
  ForwardRealFFT forward_real_fft = { 0 };
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

  int max_trim = MyMinInt(static_cast<int>(fs / 2.0 / current_f0), 6);
  FixF0(power_spectrum, numerator_i, fft_size, fs, current_f0, max_trim,
      refined_f0, score);

  delete[] diff_spectrum;
  delete[] diff_window;
  delete[] main_window;
  delete[] index_raw;
  delete[] numerator_i;
  delete[] power_spectrum;
  delete[] main_spectrum;
  DestroyForwardRealFFT(&forward_real_fft);
}

//-----------------------------------------------------------------------------
// GetRefinedF0() calculates F0 and its score based on instantaneous frequency.
//-----------------------------------------------------------------------------
static void GetRefinedF0(const double *x, int x_length, double fs,
  double current_time, double current_f0, double f0_floor, double f0_ceil,
  double *refined_f0, double *score) {
  if (current_f0 <= 0.0) {
    *refined_f0 = 0.0;
    *score = 0.0;
    return;
  }

  int half_window_length = static_cast<int>(1.5 * fs / current_f0 + 1.0);
  double window_length_in_time =
    (2.0 * static_cast<double>(half_window_length) + 1) / fs;
  double *base_time = new double[half_window_length * 2 + 1];
  for (int i = 0; i < half_window_length * 2 + 1; i++) {
    base_time[i] = static_cast<double>(-half_window_length + i) / fs;
  }
  int fft_size = static_cast<int>(pow(2.0, 2.0 +
    static_cast<int>(log(half_window_length * 2.0 + 1.0) / world::kLog2)));

  GetMeanF0(x, x_length, fs, current_time,
    current_f0, fft_size, window_length_in_time, base_time,
    half_window_length * 2 + 1, refined_f0, score);

  if (*refined_f0 < f0_floor || *refined_f0 > f0_ceil || *score < 2.5) {
    *refined_f0 = 0.0;
    *score = 0.0;
  }

  delete[] base_time;
}

//-----------------------------------------------------------------------------
// RefineF0() modifies the F0 by instantaneous frequency.
//-----------------------------------------------------------------------------
static void RefineF0Candidates(const double *x, int x_length, double fs,
  const double *time_axis, int f0_length, int max_candidates,
  double f0_floor, double f0_ceil,
  double **refined_f0_candidates, double **f0_candidates_score) {
  for (int i = 0; i < f0_length; i++) {
    for (int j = 0; j < max_candidates; ++j)
      GetRefinedF0(x, x_length, fs, time_axis[i], refined_f0_candidates[i][j],
        f0_floor, f0_ceil,
        &refined_f0_candidates[i][j], &f0_candidates_score[i][j]);
  }
}

//-----------------------------------------------------------------------------
// SelectBestF0() obtains the best one based on the candidates.
//-----------------------------------------------------------------------------
static double SelectBestF0(double reference_f0, const double *f0_candidates,
  int number_of_candidates, double allowed_range, double *best_error) {
  double best_f0 = 0.0;
  *best_error = allowed_range;

  double tmp;
  for (int i = 0; i < number_of_candidates; ++i) {
    tmp = fabs(reference_f0 - f0_candidates[i]) / reference_f0;
    if (tmp > *best_error) continue;
    best_f0 = f0_candidates[i];
    *best_error = tmp;
  }

  return best_f0;
}

//-----------------------------------------------------------------------------
// Subfunction of RemoveUnreliableCandidates().
//-----------------------------------------------------------------------------
static void RemoveUnreliableCandidatesSub(int i, int j,
    double **const tmp_f0_candidates, int number_of_candidates,
    double **f0_candidates, double **f0_candidates_score) {
  double reference_f0 = f0_candidates[i][j];
  double error1, error2, min_error;
  double threshold = 0.05;
  if (reference_f0 == 0) return;
  SelectBestF0(reference_f0, tmp_f0_candidates[i + 1],
    number_of_candidates, 1.0, &error1);
  SelectBestF0(reference_f0, tmp_f0_candidates[i - 1],
    number_of_candidates, 1.0, &error2);
  min_error = MyMinDouble(error1, error2);
  if (min_error <= threshold) return;
  f0_candidates[i][j] = 0;
  f0_candidates_score[i][j] = 0;
}

//-----------------------------------------------------------------------------
// RemoveUnreliableCandidates().
//-----------------------------------------------------------------------------
static void RemoveUnreliableCandidates(int f0_length, int number_of_candidates,
    double **f0_candidates, double **f0_candidates_score) {
  double **tmp_f0_candidates = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i) {
    tmp_f0_candidates[i] = new double[number_of_candidates];
  }
  for (int i = 1; i < f0_length - 1; ++i) {
    for (int j = 0; j < number_of_candidates; ++j) {
      tmp_f0_candidates[i][j] = f0_candidates[i][j];
    }
  }
  for (int i = 1; i < f0_length - 1; ++i) {
    for (int j = 0; j < number_of_candidates; ++j) {
      RemoveUnreliableCandidatesSub(i, j, tmp_f0_candidates,
          number_of_candidates, f0_candidates, f0_candidates_score);
    }
  }
  for (int i = 0; i < f0_length; ++i) delete[] tmp_f0_candidates[i];

  delete[] tmp_f0_candidates;
}

//-----------------------------------------------------------------------------
// SearchF0Base() gets the F0 with the highest score.
//-----------------------------------------------------------------------------
static void SearchF0Base(double **const f0_candidates,
    double **const f0_candidates_score, int f0_length,
    int number_of_candidates, double *base_f0_contour) {
  double tmp_best_score;
  for (int i = 0; i < f0_length; ++i) {
    base_f0_contour[i] = tmp_best_score = 0.0;
    for (int j = 0; j < number_of_candidates; ++j) {
      if (f0_candidates_score[i][j] > tmp_best_score) {
        base_f0_contour[i] = f0_candidates[i][j];
        tmp_best_score = f0_candidates_score[i][j];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Step 1: Rapid change of F0 contour is replaced by 0.
//-----------------------------------------------------------------------------
static void FixStep1(const double *f0_base, int f0_length,
    double allowed_range, double *f0_step1) {
  f0_step1[0] = f0_step1[1] = 0.0;
  double reference_f0;
  for (int i = 2; i < f0_length; ++i) {
    if (f0_base[i] == 0.0) continue;
    reference_f0 = f0_base[i - 1] * 2 - f0_base[i - 2];
    f0_step1[i] =
      fabs((f0_base[i] - reference_f0) / reference_f0) > allowed_range &&
      fabs((f0_base[i] - f0_base[i - 1])) / f0_base[i - 1] > allowed_range ?
      0.0 : f0_base[i];
  }
}

//-----------------------------------------------------------------------------
// GetBoundaryList() detects boundaries between voiced and unvoiced sections.
//-----------------------------------------------------------------------------
static int GetBoundaryList(const double *f0, int f0_length,
    int *boundary_list) {
  int number_of_boundaries = 0;
  int *vuv = new int[f0_length];
  for (int i = 0; i < f0_length; ++i)
    vuv[i] = f0[i] > 0 ? 1 : 0;
  vuv[0] = vuv[f0_length - 1] = 0;

  for (int i = 1; i < f0_length; ++i) {
    if (vuv[i] - vuv[i - 1] != 0) {
      boundary_list[number_of_boundaries] = i - number_of_boundaries % 2;
      number_of_boundaries++;
    }
  }

  delete[] vuv;
  return number_of_boundaries;
}

//-----------------------------------------------------------------------------
// Step 2: Voiced sections with a short period are removed.
//-----------------------------------------------------------------------------
static void FixStep2(const double *f0_step1, int f0_length,
    int voice_range_minimum, double *f0_step2) {
  for (int i = 0; i < f0_length; ++i) f0_step2[i] = f0_step1[i];
  int *boundary_list = new int[f0_length];
  int number_of_boundaries =
    GetBoundaryList(f0_step1, f0_length, boundary_list);

  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    if (boundary_list[i * 2 + 1] - boundary_list[i * 2] >= voice_range_minimum)
      continue;
    for (int j = boundary_list[i * 2]; j <= boundary_list[(i * 2) + 1]; ++j) {
      f0_step2[j] = 0.0;
    }
  }
  delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// GetMultiChannelF0() separates each voiced section into independent channel.
//-----------------------------------------------------------------------------
static void GetMultiChannelF0(const double *f0, int f0_length,
    const int *boundary_list, int number_of_boundaries,
    double **multi_channel_f0) {
  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    for (int j = 0; j < boundary_list[i * 2]; ++j) {
      multi_channel_f0[i][j] = 0.0;
    }
    for (int j = boundary_list[i * 2]; j <= boundary_list[i * 2 + 1]; ++j) {
      multi_channel_f0[i][j] = f0[j];
    }
    for (int j = boundary_list[i * 2 + 1] + 1; j < f0_length; ++j) {
      multi_channel_f0[i][j] = 0.0;
    }
  }
}

//-----------------------------------------------------------------------------
// ExtendF0() : The Hand erasing the Space.
// The subfunction of Extend().
//-----------------------------------------------------------------------------
static int ExtendF0(const double *f0, int f0_length, int origin,
    int last_point, int shift, double **const f0_candidates,
    int number_of_candidates, double allowed_range, double *extended_f0) {
  const int threshold = 4;
  double tmp_f0 = extended_f0[origin];
  int shifted_origin = origin;

  int distance = abs(last_point - origin);
  int *index_list = new int[distance + 1];
  for (int i = 0; i <= distance; ++i) index_list[i] = origin + shift * i;

  int count = 0;
  double dammy;
  for (int i = 0; i <= distance; ++i) {
    extended_f0[index_list[i] + shift] =
      SelectBestF0(tmp_f0, f0_candidates[index_list[i] + shift],
      number_of_candidates, allowed_range, &dammy);
    if (extended_f0[index_list[i] + shift] == 0.0) {
      count++;
    } else {
      tmp_f0 = extended_f0[index_list[i] + shift];
      count = 0;
      shifted_origin = index_list[i] + shift;
    }
    if (count == threshold) break;
  }

  delete[] index_list;
  return shifted_origin;
}

//-----------------------------------------------------------------------------
// Swap the f0 contour and boundary.
// It is used in ExtendSub() and MergeF0();
//-----------------------------------------------------------------------------
static void Swap(int index1, int index2, double **f0, int *boundary) {
  double *tmp_pointer;
  int tmp_index;
  tmp_pointer = f0[index1];
  f0[index1] = f0[index2];
  f0[index2] = tmp_pointer;
  tmp_index = boundary[index1 * 2];
  boundary[index1 * 2] = boundary[index2 * 2];
  boundary[index2 * 2] = tmp_index;
  tmp_index = boundary[index1 * 2 + 1];
  boundary[index1 * 2 + 1] = boundary[index2 * 2 + 1];
  boundary[index2 * 2 + 1] = tmp_index;
}

//-----------------------------------------------------------------------------
// Subfunction of ExtendF0().
//-----------------------------------------------------------------------------
static int ExtendSub(double **const extended_f0, const int *boundary_list,
    int number_of_sections, double **selected_extended_f0,
    int *selected_boundary_list) {
  const double threshold = 2200.0;
  int count = 0;
  double mean_f0 = 0.0;
  int st, ed;
  for (int i = 0; i < number_of_sections; ++i) {
    st = boundary_list[i * 2];
    ed = boundary_list[i * 2 + 1];
    for (int j = st; j < ed; ++j) mean_f0 += extended_f0[i][j];
    mean_f0 /= static_cast<double>(ed - st);
    if (threshold / mean_f0 < ed - st) {
      Swap(count++, i, selected_extended_f0, selected_boundary_list);
    }
  }
  return count;
}

//-----------------------------------------------------------------------------
// Extend() : The Hand erasing the Space.
//-----------------------------------------------------------------------------
static int Extend(double **const multi_channel_f0, int number_of_sections,
    int f0_length, const int *boundary_list, double **const f0_candidates,
    int number_of_candidates, double allowed_range, double **extended_f0,
    int *shifted_boundary_list) {
  const int threshold = 100;
  for (int i = 0; i < number_of_sections; ++i) {
    shifted_boundary_list[i * 2 + 1] = ExtendF0(multi_channel_f0[i],
      f0_length, boundary_list[i * 2 + 1],
      MyMinInt(f0_length - 2, boundary_list[i * 2 + 1] + threshold), 1,
      f0_candidates, number_of_candidates, allowed_range, extended_f0[i]);
    shifted_boundary_list[i * 2] = ExtendF0(multi_channel_f0[i], f0_length,
      boundary_list[i * 2], MyMaxInt(1, boundary_list[i * 2] - threshold), -1,
      f0_candidates, number_of_candidates, allowed_range, extended_f0[i]);
  }

  return ExtendSub(multi_channel_f0, shifted_boundary_list,
    number_of_sections, extended_f0, shifted_boundary_list);
}

//-----------------------------------------------------------------------------
// Indices are sorted.
//-----------------------------------------------------------------------------
static void MakeSortedOrder(const int *boundary_list, int number_of_sections,
    int *order) {
  for (int i = 0; i < number_of_sections; ++i) order[i] = i;
  int tmp;
  for (int i = 1; i < number_of_sections; ++i) {
    for (int j = i - 1; j >= 0; --j) {
      if (boundary_list[order[j] * 2] > boundary_list[order[i] * 2]) {
        tmp = order[i];
        order[i] = order[j];
        order[j] = tmp;
      } else {
        break;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Serach the highest score with the candidate F0.
//-----------------------------------------------------------------------------
static double SearchScore(double f0, const double *f0_candidates,
    const double  *f0_candidates_score, int number_of_candidates) {
  double score = 0.0;
  for (int i = 0; i < number_of_candidates; ++i) {
    if (f0 == f0_candidates[i] && score < f0_candidates_score[i]) {
      score = f0_candidates_score[i];
    }
  }
  return score;
}

//-----------------------------------------------------------------------------
// Subfunction of MergeF0()
//-----------------------------------------------------------------------------
static int MergeF0Sub(const double *f0_1, int f0_length, int st1, int ed1,
    const double *f0_2, int st2, int ed2, double **const f0_candidates,
    double **const f0_candidates_score, int number_of_candidates,
    double *merged_f0) {
  if (st1 <= st2 && ed1 >= ed2)
    return ed1;

  double score1 = 0.0;
  double score2 = 0.0;
  for (int i = st2; i <= ed1; ++i) {
    score1 += SearchScore(f0_1[i], f0_candidates[i], f0_candidates_score[i],
      number_of_candidates);
    score2 += SearchScore(f0_2[i], f0_candidates[i], f0_candidates_score[i],
      number_of_candidates);
  }
  if (score1 > score2) {
    for (int i = ed1; i <= ed2; ++i) merged_f0[i] = f0_2[i];
  } else {
    for (int i = st2; i <= ed2; ++i) merged_f0[i] = f0_2[i];
  }

  return ed2;
}

//-----------------------------------------------------------------------------
// Overlapped F0 contours are merged by the likability score.
//-----------------------------------------------------------------------------
static void MergeF0(double **const multi_channel_f0, int *boundary_list,
    int number_of_channels, int f0_length, double **const f0_candidates,
    double **const f0_candidates_score, int number_of_candidates,
    double *merged_f0) {
  int *order = new int[number_of_channels];
  MakeSortedOrder(boundary_list, number_of_channels, order);

  for (int i = 0; i < f0_length; ++i)
    merged_f0[i] = multi_channel_f0[0][i];

  for (int i = 1; i < number_of_channels; ++i) {
    if (boundary_list[order[i] * 2] - boundary_list[1] > 0) {
      for (int j = boundary_list[order[i] * 2];
        j <= boundary_list[order[i] * 2 + 1]; ++j) {
        merged_f0[j] = multi_channel_f0[order[i]][j];
      }
      boundary_list[0] = boundary_list[order[i] * 2];
      boundary_list[1] = boundary_list[order[i] * 2 + 1];
    } else {
      boundary_list[1] =
        MergeF0Sub(merged_f0, f0_length, boundary_list[0], boundary_list[1],
        multi_channel_f0[order[i]], boundary_list[order[i] * 2],
        boundary_list[order[i] * 2 + 1], f0_candidates, f0_candidates_score,
        number_of_candidates, merged_f0);
    }
  }

  delete[] order;
}

//-----------------------------------------------------------------------------
// Step 3: Voiced sections are extended based on the continuity of F0 contour
//-----------------------------------------------------------------------------
static void FixStep3(const double *f0_step2, int f0_length,
    int number_of_candidates, double **const f0_candidates,
    double allowed_range, double **const f0_candidates_score,
    double *f0_step3) {
  for (int i = 0; i < f0_length; ++i) f0_step3[i] = f0_step2[i];
  int *boundary_list = new int[f0_length];
  int number_of_boundaries =
    GetBoundaryList(f0_step2, f0_length, boundary_list);

  double **multi_channel_f0 = new double *[number_of_boundaries / 2];
  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    multi_channel_f0[i] = new double[f0_length];
  }
  GetMultiChannelF0(f0_step2, f0_length, boundary_list, number_of_boundaries,
      multi_channel_f0);

  int number_of_channels =
    Extend(multi_channel_f0, number_of_boundaries / 2, f0_length,
    boundary_list, f0_candidates, number_of_candidates, allowed_range,
    multi_channel_f0, boundary_list);

  MergeF0(multi_channel_f0, boundary_list, number_of_channels, f0_length,
      f0_candidates, f0_candidates_score, number_of_candidates, f0_step3);

  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    delete[] multi_channel_f0[i];
  }
  delete[] multi_channel_f0;
  delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// Step 4: F0s in short unvoiced section are faked
//-----------------------------------------------------------------------------
static void FixStep4(const double *f0_step3, int f0_length, int threshold,
    double *f0_step4) {
  for (int i = 0; i < f0_length; ++i) f0_step4[i] = f0_step3[i];
  int *boundary_list = new int[f0_length];
  int number_of_boundaries =
    GetBoundaryList(f0_step3, f0_length, boundary_list);

  int distance;
  double tmp0, tmp1, coefficient;
  int count;
  for (int i = 0; i < number_of_boundaries / 2 - 1; ++i) {
    distance = boundary_list[(i + 1) * 2] - boundary_list[i * 2 + 1] - 1;
    if (distance >= threshold) continue;
    tmp0 = f0_step3[boundary_list[i * 2 + 1]] + 1;
    tmp1 = f0_step3[boundary_list[(i + 1) * 2]] - 1;
    coefficient = static_cast<double>(tmp1 - tmp0) / (distance + 1);
    count = 1;
    for (int j = boundary_list[i * 2 + 1] + 1;
        j <= boundary_list[(i + 1) * 2] - 1; ++j) {
      f0_step4[j] = tmp0 + coefficient * count;
      count++;
    }
  }
  delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// FixF0Contour() obtains the likely F0 contour.
//-----------------------------------------------------------------------------
static void FixF0Contour(double **const f0_candidates,
    double **const f0_candidates_score, int f0_length,
    int number_of_candidates, double *best_f0_contour) {
  double *tmp_f0_contour1 = new double[f0_length];
  double *tmp_f0_contour2 = new double[f0_length];

  SearchF0Base(f0_candidates, f0_candidates_score, f0_length,
      number_of_candidates, tmp_f0_contour1);
  FixStep1(tmp_f0_contour1, f0_length, 0.008, tmp_f0_contour2);
  FixStep2(tmp_f0_contour2, f0_length, 6, tmp_f0_contour1);
  FixStep3(tmp_f0_contour1, f0_length, number_of_candidates, f0_candidates,
      0.18, f0_candidates_score, tmp_f0_contour2);
  FixStep4(tmp_f0_contour2, f0_length, 9, best_f0_contour);

  delete[] tmp_f0_contour1;
  delete[] tmp_f0_contour2;
}

//-----------------------------------------------------------------------------
// This function uses zero-lag Butterworth filter.
//-----------------------------------------------------------------------------
static void FilteringF0(const double *a, const double *b, double *x,
    int x_length, int st, int ed, double *y) {
  double w[2] = { 0.0, 0.0 };
  double wt;
  double *tmp_x = new double[x_length];

  for (int i = 0; i < st; ++i) x[i] = x[st];
  for (int i = ed + 1; i < x_length; ++i) x[i] = x[ed];

  for (int i = 0; i < x_length; ++i) {
    wt = x[i] + a[0] * w[0] + a[1] * w[1];
    tmp_x[x_length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
    w[1] = w[0];
    w[0] = wt;
  }

  w[0] = w[1] = 0.0;
  for (int i = 0; i < x_length; ++i) {
    wt = tmp_x[i] + a[0] * w[0] + a[1] * w[1];
    y[x_length - i - 1] = b[0] * wt + b[1] * w[0] + b[0] * w[1];
    w[1] = w[0];
    w[0] = wt;
  }

  delete[] tmp_x;
}

//-----------------------------------------------------------------------------
// SmoothF0Contour() uses the zero-lag Butterworth filter for smoothing.
//-----------------------------------------------------------------------------
static void SmoothF0Contour(const double *f0, int f0_length,
    double *smoothed_f0) {
  const double b[2] =
    { 0.0078202080334971724, 0.015640416066994345 };
  const double a[2] =
    { 1.7347257688092754, -0.76600660094326412 };
  int lag = 300;
  int new_f0_length = f0_length + lag * 2;
  double *f0_contour = new double[new_f0_length];
  for (int i = 0; i < lag; ++i) f0_contour[i] = 0.0;
  for (int i = lag; i < lag + f0_length; ++i) f0_contour[i] = f0[i - lag];
  for (int i = lag + f0_length; i < new_f0_length; ++i) f0_contour[i] = 0.0;

  int *boundary_list = new int[new_f0_length];
  int number_of_boundaries =
    GetBoundaryList(f0_contour, new_f0_length, boundary_list);
  double **multi_channel_f0 = new double *[number_of_boundaries / 2];
  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    multi_channel_f0[i] = new double[new_f0_length];
  }
  GetMultiChannelF0(f0_contour, new_f0_length, boundary_list,
      number_of_boundaries, multi_channel_f0);

  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    FilteringF0(a, b, multi_channel_f0[i], new_f0_length,
      boundary_list[i * 2], boundary_list[i * 2 + 1], f0_contour);
    for (int j = boundary_list[i * 2]; j <= boundary_list[i * 2 + 1]; ++j) {
      smoothed_f0[j - lag] = f0_contour[j];
    }
  }

  for (int i = 0; i < number_of_boundaries / 2; ++i) {
    delete[] multi_channel_f0[i];
  }
  delete[] f0_contour;
  delete[] boundary_list;
}

//-----------------------------------------------------------------------------
// HarvestGeneralBody() estimates the F0 contour based on Harvest.
//-----------------------------------------------------------------------------
static void HarvestGeneralBody(const double *x, int x_length, int fs,
    int frame_period, double f0_floor, double f0_ceil,
    double channels_in_octave, int speed, double *time_axis, double *f0) {
  double f0_floor_adjusted = f0_floor * 0.9;
  double f0_ceil_adjusted = f0_ceil * 1.1;
  int number_of_channels =
    1 + static_cast<int>(log(f0_ceil_adjusted / f0_floor_adjusted) /
    world::kLog2 * channels_in_octave);
  double *boundary_f0_list = new double[number_of_channels];
  for (int i = 0; i < number_of_channels; ++i) {
    boundary_f0_list[i] =
      f0_floor_adjusted * pow(2.0, (i + 1) / channels_in_octave);
  }

  // normalization
  int decimation_ratio = MyMaxInt(MyMinInt(speed, 12), 1);
  int y_length = (1 + static_cast<int>(x_length / decimation_ratio));
  double actual_fs = static_cast<double>(fs) / decimation_ratio;
  int fft_size = GetSuitableFFTSize(y_length +
      (4 * static_cast<int>(1.0 + actual_fs / boundary_f0_list[0] / 2.0)));

  // Calculation of the spectrum used for the f0 estimation
  double *y = new double[fft_size];
  fft_complex *y_spectrum = new fft_complex[fft_size];
  GetWaveformAndSpectrum(x, x_length, y_length, actual_fs, fft_size,
      decimation_ratio, y, y_spectrum);

  int f0_length = GetSamplesForHarvest(fs, x_length, frame_period);
  for (int i = 0; i < f0_length; ++i) {
    time_axis[i] = i * frame_period / 1000.0;
  }
  double **raw_f0_candidates = new double *[number_of_channels];
  for (int i = 0; i < number_of_channels; ++i) {
    raw_f0_candidates[i] = new double[f0_length];
  }
  GetF0CandidateMap(boundary_f0_list, number_of_channels,
      actual_fs, y_length, time_axis, f0_length, y_spectrum,
      fft_size, f0_floor, f0_ceil, raw_f0_candidates);

  int max_candidates = matlab_round(number_of_channels / 10) * 7;
  double **f0_candidates = new double *[f0_length];
  double **f0_candidates_score = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i) {
    f0_candidates[i] = new double[max_candidates];
    f0_candidates_score[i] = new double[max_candidates];
  }

  int number_of_candidates = DetectF0Candidates(raw_f0_candidates,
      number_of_channels, f0_length, max_candidates, f0_candidates);
  OverlapF0Candidates(f0_length, number_of_candidates, f0_candidates);
  number_of_candidates = number_of_candidates * 7;
  RefineF0Candidates(y, y_length, actual_fs, time_axis, f0_length,
      number_of_candidates, f0_floor, f0_ceil, f0_candidates,
      f0_candidates_score);
  RemoveUnreliableCandidates(f0_length, number_of_candidates,
      f0_candidates, f0_candidates_score);

  double *best_f0_contour = new double[f0_length];
  FixF0Contour(f0_candidates, f0_candidates_score, f0_length,
      number_of_candidates, best_f0_contour);
  SmoothF0Contour(best_f0_contour, f0_length, f0);

  delete[] y;
  delete[] best_f0_contour;
  delete[] y_spectrum;
  for (int i = 0; i < f0_length; ++i) {
    delete[] f0_candidates[i];
    delete[] f0_candidates_score[i];
  }
  delete[] f0_candidates;
  delete[] f0_candidates_score;
  for (int i = 0; i < number_of_channels; ++i) {
    delete[] raw_f0_candidates[i];
  }
  delete[] raw_f0_candidates;
  delete[] boundary_f0_list;
}

}  // namespace

int GetSamplesForHarvest(int fs, int x_length, double frame_period) {
  return static_cast<int>(x_length / static_cast<double>(fs) /
    (frame_period / 1000.0)) + 1;
}

void Harvest(const double *x, int x_length, int fs,
    const HarvestOption *option, double *time_axis, double *f0) {
  // Several parameters will be controllable for debug.
  double target_fs = 8000.0;
  int dimension_ratio = matlab_round(fs / target_fs);
  double channels_in_octave = 40;

  if (option->frame_period == 1.0) {
    HarvestGeneralBody(x, x_length, fs, 1, option->f0_floor,
      option->f0_ceil, channels_in_octave, dimension_ratio, time_axis, f0);
    return;
  }

  int basic_frame_period = 1;
  int basic_f0_length =
    GetSamplesForHarvest(fs, x_length, basic_frame_period);
  double *basic_f0 = new double[basic_f0_length];
  double *basic_time_axis = new double[basic_f0_length];
  HarvestGeneralBody(x, x_length, fs, basic_frame_period, option->f0_floor,
      option->f0_ceil, channels_in_octave, dimension_ratio, basic_time_axis,
      basic_f0);

  int f0_length = GetSamplesForHarvest(fs, x_length, option->frame_period);
  for (int i = 0; i < f0_length; ++i) {
    time_axis[i] = i * option->frame_period / 1000.0;
    f0[i] = basic_f0[MyMinInt(basic_f0_length - 1,
      matlab_round(time_axis[i] * 1000.0))];
  }

  delete[] basic_f0;
  delete[] basic_time_axis;
}

void InitializeHarvestOption(HarvestOption *option) {
  // You can change default parameters.
  option->f0_ceil = world::kCeilF0;
  option->f0_floor = world::kFloorF0;
  option->frame_period = 5;
}
