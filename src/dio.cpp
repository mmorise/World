//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/03/04
//
// F0 estimation based on DIO (Distributed Inline-filter Operation).
//-----------------------------------------------------------------------------
#include "world/dio.h"

#include <math.h>

#include "world/common.h"
#include "world/constantnumbers.h"
#include "world/matlabfunctions.h"

//-----------------------------------------------------------------------------
// struct for GetFourZeroCrossingIntervals()
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
// DesignLowCutFilter() calculates the coefficients the filter.
//-----------------------------------------------------------------------------
static void DesignLowCutFilter(int N, int fft_size, double *low_cut_filter) {
  for (int i = 1; i <= N; ++i)
    low_cut_filter[i - 1] = 0.5 - 0.5 * cos(i * 2.0 * world::kPi / (N + 1));
  for (int i = N; i < fft_size; ++i) low_cut_filter[i] = 0.0;
  double sum_of_amplitude = 0.0;
  for (int i = 0; i < N; ++i) sum_of_amplitude += low_cut_filter[i];
  for (int i = 0; i < N; ++i)
    low_cut_filter[i] = -low_cut_filter[i] / sum_of_amplitude;
  for (int i = 0; i < (N - 1) / 2; ++i)
    low_cut_filter[fft_size - (N - 1) / 2 + i] = low_cut_filter[i];
  for (int i = 0; i < N; ++i)
    low_cut_filter[i] = low_cut_filter[i + (N - 1) / 2];
  low_cut_filter[0] += 1.0;
}

//-----------------------------------------------------------------------------
// GetSpectrumForEstimation() calculates the spectrum for estimation.
// This function carries out downsampling to speed up the estimation process
// and calculates the spectrum of the downsampled signal.
//-----------------------------------------------------------------------------
static void GetSpectrumForEstimation(const double *x, int x_length,
    int y_length, double actual_fs, int fft_size, int decimation_ratio,
    fft_complex *y_spectrum) {
  double *y = new double[fft_size];

  // Initialization
  for (int i = 0; i < fft_size; ++i) y[i] = 0.0;

  // Downsampling
  if (decimation_ratio != 1)
    decimate(x, x_length, decimation_ratio, y);
  else
    for (int i = 0; i < x_length; ++i) y[i] = x[i];

  // Removal of the DC component (y = y - mean value of y)
  double mean_y = 0.0;
  for (int i = 0; i < y_length; ++i) mean_y += y[i];
  mean_y /= y_length;
  for (int i = 0; i < y_length; ++i) y[i] -= mean_y;
  for (int i = y_length; i < fft_size; ++i) y[i] = 0.0;

  fft_plan forwardFFT =
    fft_plan_dft_r2c_1d(fft_size, y, y_spectrum, FFT_ESTIMATE);
  fft_execute(forwardFFT);

  // Low cut filtering (from 0.1.4)
  int cutoff_in_sample = matlab_round(actual_fs / 50.0);  // Cutoff is 50.0 Hz
  DesignLowCutFilter(cutoff_in_sample * 2 + 1, fft_size, y);

  fft_complex *filter_spectrum = new fft_complex[fft_size];
  forwardFFT.c_out = filter_spectrum;
  fft_execute(forwardFFT);

  double tmp = 0;
  for (int i = 0; i <= fft_size / 2; ++i) {
    // Complex number multiplications.
    tmp = y_spectrum[i][0] * filter_spectrum[i][0] -
      y_spectrum[i][1] * filter_spectrum[i][1];
    y_spectrum[i][1] = y_spectrum[i][0] * filter_spectrum[i][1] +
      y_spectrum[i][1] * filter_spectrum[i][0];
    y_spectrum[i][0] = tmp;
  }

  fft_destroy_plan(forwardFFT);
  delete[] y;
  delete[] filter_spectrum;
}

//-----------------------------------------------------------------------------
// GetBestF0Contour() calculates the best f0 contour based on scores of
// all candidates. The F0 with highest score is selected.
//-----------------------------------------------------------------------------
static void GetBestF0Contour(int f0_length,
    const double * const * f0_candidates, const double * const * f0_scores,
    int number_of_bands, double *best_f0_contour) {
  double tmp;
  for (int i = 0; i < f0_length; ++i) {
    tmp = f0_scores[0][i];
    best_f0_contour[i] = f0_candidates[0][i];
    for (int j = 1; j < number_of_bands; ++j) {
      if (tmp > f0_scores[j][i]) {
        tmp = f0_scores[j][i];
        best_f0_contour[i] = f0_candidates[j][i];
      }
    }
  }
}

//-----------------------------------------------------------------------------
// FixStep1() is the 1st step of the postprocessing.
// This function eliminates the unnatural change of f0 based on allowed_range.
//-----------------------------------------------------------------------------
static void FixStep1(const double *best_f0_contour, int f0_length,
    int voice_range_minimum, double allowed_range, double *f0_step1) {
  double *f0_base = new double[f0_length];
  // Initialization
  for (int i = 0; i < voice_range_minimum; ++i) f0_base[i] = 0.0;
  for (int i = voice_range_minimum; i < f0_length - voice_range_minimum; ++i)
    f0_base[i] = best_f0_contour[i];
  for (int i = f0_length - voice_range_minimum; i < f0_length; ++i)
    f0_base[i] = 0.0;

  // Processing to prevent the jumping of f0
  for (int i = 0; i < voice_range_minimum; ++i) f0_step1[i] = 0.0;
  for (int i = voice_range_minimum; i < f0_length; ++i)
    f0_step1[i] = fabs((f0_base[i] - f0_base[i - 1]) /
    (world::kMySafeGuardMinimum + f0_base[i])) <
    allowed_range ? f0_base[i] : 0.0;

  delete[] f0_base;
}

//-----------------------------------------------------------------------------
// FixStep2() is the 2nd step of the postprocessing.
// This function eliminates the suspected f0 in the anlaut and auslaut.
//-----------------------------------------------------------------------------
static void FixStep2(const double *f0_step1, int f0_length,
    int voice_range_minimum, double *f0_step2) {
  for (int i = 0; i < f0_length; ++i) f0_step2[i] = f0_step1[i];

  int center = (voice_range_minimum - 1) / 2;
  for (int i = center; i < f0_length - center; ++i) {
    for (int j = -center; j <= center; ++j) {
      if (f0_step1[i + j] == 0) {
        f0_step2[i] = 0.0;
        break;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// GetNumberOfVoicedSections() counts the number of voiced sections.
//-----------------------------------------------------------------------------
static void GetNumberOfVoicedSections(const double *f0, int f0_length,
    int *positive_index, int *negative_index, int *positive_count,
    int *negative_count) {
  *positive_count = *negative_count = 0;
  for (int i = 1; i < f0_length; ++i)
    if (f0[i] == 0 && f0[i - 1] != 0)
      negative_index[(*negative_count)++] = i - 1;
    else
      if (f0[i - 1] == 0 && f0[i] != 0)
        positive_index[(*positive_count)++] = i;
}

//-----------------------------------------------------------------------------
// SelectOneF0() corrects the f0[current_index] based on
// f0[current_index + sign].
//-----------------------------------------------------------------------------
static double SelectBestF0(double current_f0, double past_f0,
    const double * const * f0_candidates, int number_of_candidates,
    int target_index, double allowed_range) {
  double reference_f0 = (current_f0 * 3.0 - past_f0) / 2.0;

  double minimum_error = fabs(reference_f0 - f0_candidates[0][target_index]);
  double best_f0 = f0_candidates[0][target_index];

  double current_error;
  for (int i = 1; i < number_of_candidates; ++i) {
    current_error = fabs(reference_f0 - f0_candidates[i][target_index]);
    if (current_error < minimum_error) {
      minimum_error = current_error;
      best_f0 = f0_candidates[i][target_index];
    }
  }
  if (fabs(1.0 - best_f0 / reference_f0) > allowed_range)
    return 0.0;
  return best_f0;
}

//-----------------------------------------------------------------------------
// FixStep3() is the 3rd step of the postprocessing.
// This function corrects the f0 candidates from backward to forward.
//-----------------------------------------------------------------------------
static void FixStep3(const double *f0_step2, int f0_length,
    const double * const * f0_candidates, int number_of_candidates,
    double allowed_range, const int *negative_index, int negative_count,
    double *f0_step3) {
  for (int i = 0; i < f0_length; i++) f0_step3[i] = f0_step2[i];

  int limit;
  for (int i = 0; i < negative_count; ++i) {
    limit = i == negative_count - 1 ? f0_length - 1 : negative_index[i + 1];
    for (int j = negative_index[i]; j < limit; ++j) {
      f0_step3[j + 1] =
        SelectBestF0(f0_step3[j], f0_step3[j - 1], f0_candidates,
            number_of_candidates, j + 1, allowed_range);
      if (f0_step3[j + 1] == 0) break;
    }
  }
}

//-----------------------------------------------------------------------------
// FixStep4() is the 4th step of the postprocessing.
// This function corrects the f0 candidates from forward to backward.
//-----------------------------------------------------------------------------
static void FixStep4(const double *f0_step3, int f0_length,
    const double * const * f0_candidates, int number_of_candidates,
    double allowed_range, const int *positive_index, int positive_count,
    double *f0_step4) {
  for (int i = 0; i < f0_length; ++i) f0_step4[i] = f0_step3[i];

  int limit;
  for (int i = positive_count - 1; i >= 0; --i) {
    limit = i == 0 ? 1 : positive_index[i - 1];
    for (int j = positive_index[i]; j > limit; --j) {
      f0_step4[j - 1] =
        SelectBestF0(f0_step4[j], f0_step4[j + 1], f0_candidates,
            number_of_candidates, j - 1, allowed_range);
      if (f0_step4[j - 1] == 0) break;
    }
  }
}

//-----------------------------------------------------------------------------
// FixF0Contour() calculates the definitive f0 contour based on all f0
// candidates. There are four steps.
//-----------------------------------------------------------------------------
static void FixF0Contour(double frame_period, int number_of_candidates,
    int fs, const double * const * f0_candidates,
    const double *best_f0_contour, int f0_length, double f0_floor,
    double allowed_range, double *fixed_f0_contour) {
  int voice_range_minimum =
    static_cast<int>(0.5 + 1000.0 / frame_period / f0_floor) * 2 + 1;

  if (f0_length <= voice_range_minimum) return;

  double *f0_tmp1 = new double[f0_length];
  double *f0_tmp2 = new double[f0_length];

  FixStep1(best_f0_contour, f0_length, voice_range_minimum,
      allowed_range, f0_tmp1);
  FixStep2(f0_tmp1, f0_length, voice_range_minimum, f0_tmp2);

  int positive_count, negative_count;
  int *positive_index = new int[f0_length];
  int *negative_index = new int[f0_length];
  GetNumberOfVoicedSections(f0_tmp2, f0_length, positive_index,
      negative_index, &positive_count, &negative_count);
  FixStep3(f0_tmp2, f0_length, f0_candidates, number_of_candidates,
      allowed_range, negative_index, negative_count, f0_tmp1);
  FixStep4(f0_tmp1, f0_length, f0_candidates, number_of_candidates,
      allowed_range, positive_index, positive_count, fixed_f0_contour);

  delete[] f0_tmp1;
  delete[] f0_tmp2;
  delete[] positive_index;
  delete[] negative_index;
}

//-----------------------------------------------------------------------------
// GetFilteredSignal() calculates the signal that is the convolution of the
// input signal and low-pass filter.
// This function is only used in RawEventByDio()
//-----------------------------------------------------------------------------
static void GetFilteredSignal(int half_average_length, int fft_size,
    const fft_complex *y_spectrum, int y_length, double *filtered_signal) {
  double *low_pass_filter = new double[fft_size];
  // Nuttall window is used as a low-pass filter.
  // Cutoff frequency depends on the window length.
  NuttallWindow(half_average_length * 4, low_pass_filter);
  for (int i = half_average_length * 4; i < fft_size; ++i)
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
  int index_bias = half_average_length * 2;
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
// negative. Thanks to Custom.Maid http://custom-made.seesaa.net/ (2012/8/19)
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
  // x_length / 4 (old version) is fixed at 2013/07/14
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
// GetF0CandidateContourSub() calculates the f0 candidates and deviations.
// This is the sub-function of GetF0Candidates() and assumes the calculation.
//-----------------------------------------------------------------------------
static void GetF0CandidateContourSub(
    const double * const * interpolated_f0_set, int f0_length, double f0_floor,
    double f0_ceil, double boundary_f0, double *f0_candidate,
    double *f0_score) {
  for (int i = 0; i < f0_length; ++i) {
    f0_candidate[i] = (interpolated_f0_set[0][i] +
      interpolated_f0_set[1][i] + interpolated_f0_set[2][i] +
      interpolated_f0_set[3][i]) / 4.0;

    f0_score[i] = sqrt(((interpolated_f0_set[0][i] - f0_candidate[i]) *
      (interpolated_f0_set[0][i] - f0_candidate[i]) +
      (interpolated_f0_set[1][i] - f0_candidate[i]) *
      (interpolated_f0_set[1][i] - f0_candidate[i]) +
      (interpolated_f0_set[2][i] - f0_candidate[i]) *
      (interpolated_f0_set[2][i] - f0_candidate[i]) +
      (interpolated_f0_set[3][i] - f0_candidate[i]) *
      (interpolated_f0_set[3][i] - f0_candidate[i])) / 3.0);

    if (f0_candidate[i] > boundary_f0 || f0_candidate[i] < boundary_f0 / 2.0 ||
        f0_candidate[i] > f0_ceil || f0_candidate[i] < f0_floor) {
      f0_candidate[i] = 0.0;
      f0_score[i] = world::kMaximumValue;
    }
  }
}

//-----------------------------------------------------------------------------
// GetF0CandidateContour() calculates the F0 candidates based on the
// zero-crossings.
//-----------------------------------------------------------------------------
static void GetF0CandidateContour(const ZeroCrossings *zero_crossings,
    double boundary_f0, double f0_floor, double f0_ceil,
    const double *temporal_positions, int f0_length,
    double *f0_candidate, double *f0_score) {
  if (0 == CheckEvent(zero_crossings->number_of_negatives - 2) *
      CheckEvent(zero_crossings->number_of_positives - 2) *
      CheckEvent(zero_crossings->number_of_peaks - 2) *
      CheckEvent(zero_crossings->number_of_dips - 2)) {
    for (int i = 0; i < f0_length; ++i) {
      f0_score[i] = world::kMaximumValue;
      f0_candidate[i] = 0.0;
    }
    return;
  }

  double *interpolated_f0_set[4];
  for (int i = 0; i < 4; ++i)
    interpolated_f0_set[i] = new double[f0_length];

  interp1(zero_crossings->negative_interval_locations,
      zero_crossings->negative_intervals,
      zero_crossings->number_of_negatives,
      temporal_positions, f0_length, interpolated_f0_set[0]);
  interp1(zero_crossings->positive_interval_locations,
      zero_crossings->positive_intervals,
      zero_crossings->number_of_positives,
      temporal_positions, f0_length, interpolated_f0_set[1]);
  interp1(zero_crossings->peak_interval_locations,
      zero_crossings->peak_intervals, zero_crossings->number_of_peaks,
      temporal_positions, f0_length, interpolated_f0_set[2]);
  interp1(zero_crossings->dip_interval_locations,
      zero_crossings->dip_intervals, zero_crossings->number_of_dips,
      temporal_positions, f0_length, interpolated_f0_set[3]);

  GetF0CandidateContourSub(interpolated_f0_set, f0_length, f0_floor,
      f0_ceil, boundary_f0, f0_candidate, f0_score);
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
// GetF0CandidateFromRawEvent() calculates F0 candidate contour in 1-ch signal
//-----------------------------------------------------------------------------
static void GetF0CandidateFromRawEvent(double boundary_f0, double fs,
    const fft_complex *y_spectrum, int y_length, int fft_size, double f0_floor,
    double f0_ceil, const double *temporal_positions, int f0_length,
    double *f0_score, double *f0_candidate) {
  double *filtered_signal = new double[fft_size];
  GetFilteredSignal(matlab_round(fs / boundary_f0 / 2.0), fft_size, y_spectrum,
      y_length, filtered_signal);

  ZeroCrossings zero_crossings = {0};
  GetFourZeroCrossingIntervals(filtered_signal, y_length, fs,
      &zero_crossings);

  GetF0CandidateContour(&zero_crossings, boundary_f0, f0_floor, f0_ceil,
      temporal_positions, f0_length, f0_candidate, f0_score);

  DestroyZeroCrossings(&zero_crossings);
  delete[] filtered_signal;
}

//-----------------------------------------------------------------------------
// GetF0CandidatesAndScores() calculates all f0 candidates and their scores.
//-----------------------------------------------------------------------------
static void GetF0CandidatesAndScores(const double *boundary_f0_list,
    int number_of_bands, double actual_fs, int y_length,
    const double *temporal_positions, int f0_length,
    const fft_complex *y_spectrum, int fft_size, double f0_floor,
    double f0_ceil, double **raw_f0_candidates, double **raw_f0_scores) {
  double *f0_candidate = new double[f0_length];
  double *f0_score = new double[f0_length];

  // Calculation of the acoustics events (zero-crossing)
  for (int i = 0; i < number_of_bands; ++i) {
    GetF0CandidateFromRawEvent(boundary_f0_list[i], actual_fs, y_spectrum,
        y_length, fft_size, f0_floor, f0_ceil, temporal_positions, f0_length,
        f0_score, f0_candidate);
    for (int j = 0; j < f0_length; ++j) {
      // A way to avoid zero division
      raw_f0_scores[i][j] = f0_score[j] /
        (f0_candidate[j] + world::kMySafeGuardMinimum);
      raw_f0_candidates[i][j] = f0_candidate[j];
    }
  }

  delete[] f0_candidate;
  delete[] f0_score;
}

//-----------------------------------------------------------------------------
// DioGeneralBody() estimates the F0 based on Distributed Inline-filter
// Operation.
//-----------------------------------------------------------------------------
static void DioGeneralBody(const double *x, int x_length, int fs,
    double frame_period, double f0_floor, double f0_ceil,
    double channels_in_octave, int speed, double allowed_range,
    double *temporal_positions, double *f0) {
  int number_of_bands = 1 + static_cast<int>(log(f0_ceil / f0_floor) /
    world::kLog2 * channels_in_octave);
  double *boundary_f0_list = new double[number_of_bands];
  for (int i = 0; i < number_of_bands; ++i)
    boundary_f0_list[i] = f0_floor * pow(2.0, (i + 1) / channels_in_octave);

  // normalization
  int decimation_ratio = MyMaxInt(MyMinInt(speed, 12), 1);
  int y_length = (1 + static_cast<int>(x_length / decimation_ratio));
  double actual_fs = static_cast<double>(fs) / decimation_ratio;
  int fft_size = GetSuitableFFTSize(y_length +
      (4 * static_cast<int>(1.0 + actual_fs / boundary_f0_list[0] / 2.0)));

  // Calculation of the spectrum used for the f0 estimation
  fft_complex *y_spectrum = new fft_complex[fft_size];
  GetSpectrumForEstimation(x, x_length, y_length, actual_fs, fft_size,
      decimation_ratio, y_spectrum);

  double **f0_candidates = new double *[number_of_bands];
  double **f0_scores = new double *[number_of_bands];
  int f0_length = GetSamplesForDIO(fs, x_length, frame_period);
  for (int i = 0; i < number_of_bands; ++i) {
    f0_candidates[i] = new double[f0_length];
    f0_scores[i] = new double[f0_length];
  }

  for (int i = 0; i < f0_length; ++i)
    temporal_positions[i] = i * frame_period / 1000.0;

  GetF0CandidatesAndScores(boundary_f0_list, number_of_bands,
      actual_fs, y_length, temporal_positions, f0_length, y_spectrum,
      fft_size, f0_floor, f0_ceil, f0_candidates, f0_scores);

  // Selection of the best value based on fundamental-ness.
  // This function is related with SortCandidates() in MATLAB.
  double *best_f0_contour = new double[f0_length];
  GetBestF0Contour(f0_length, f0_candidates, f0_scores,
      number_of_bands, best_f0_contour);

  // Postprocessing to find the best f0-contour.
  FixF0Contour(frame_period, number_of_bands, fs, f0_candidates,
      best_f0_contour, f0_length, f0_floor, allowed_range, f0);

  delete[] best_f0_contour;
  delete[] y_spectrum;
  for (int i = 0; i < number_of_bands; ++i) {
    delete[] f0_scores[i];
    delete[] f0_candidates[i];
  }
  delete[] f0_scores;
  delete[] f0_candidates;
  delete[] boundary_f0_list;
}

}  // namespace

int GetSamplesForDIO(int fs, int x_length, double frame_period) {
  return static_cast<int>(1000.0 * x_length / fs / frame_period) + 1;
}

void Dio(const double *x, int x_length, int fs, const DioOption *option,
    double *temporal_positions, double *f0) {
  DioGeneralBody(x, x_length, fs, option->frame_period, option->f0_floor,
      option->f0_ceil, option->channels_in_octave, option->speed,
      option->allowed_range, temporal_positions, f0);
}

void InitializeDioOption(DioOption *option) {
  // You can change default parameters.
  option->channels_in_octave = 2.0;
  option->f0_ceil = world::kCeilF0;
  option->f0_floor = world::kFloorF0;
  option->frame_period = 5;

  // You can use the value from 1 to 12.
  // Default value 11 is for the fs of 44.1 kHz.
  // The lower value you use, the better performance you can obtain.
  option->speed = 1;

  // You can give a positive real number as the threshold.
  // The most strict value is 0, and there is no upper limit.
  // On the other hand, I think that the value from 0.02 to 0.2 is reasonable.
  option->allowed_range = 0.1;
}
