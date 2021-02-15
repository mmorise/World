//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//
// Voice synthesis based on f0, spectrogram and aperiodicity.
// forward_real_fft, inverse_real_fft and minimum_phase are used to speed up.
//-----------------------------------------------------------------------------
#include "world/synthesis.h"

#include <math.h>

#include "world/common.h"
#include "world/constantnumbers.h"
#include "world/matlabfunctions.h"

namespace {

static void GetNoiseSpectrum(int noise_size, int fft_size,
    const ForwardRealFFT *forward_real_fft) {
  double average = 0.0;
  for (int i = 0; i < noise_size; ++i) {
    forward_real_fft->waveform[i] = randn();
    average += forward_real_fft->waveform[i];
  }

  average /= noise_size;
  for (int i = 0; i < noise_size; ++i)
    forward_real_fft->waveform[i] -= average;
  for (int i = noise_size; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);
}

//-----------------------------------------------------------------------------
// GetAperiodicResponse() calculates an aperiodic response.
//-----------------------------------------------------------------------------
static void GetAperiodicResponse(int noise_size, int fft_size,
    const double *spectrum, const double *aperiodic_ratio, double current_vuv,
    const ForwardRealFFT *forward_real_fft,
    const InverseRealFFT *inverse_real_fft,
    const MinimumPhaseAnalysis *minimum_phase, double *aperiodic_response) {
  GetNoiseSpectrum(noise_size, fft_size, forward_real_fft);

  if (current_vuv != 0.0)
    for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
      minimum_phase->log_spectrum[i] =
        log(spectrum[i] * aperiodic_ratio[i]) / 2.0;
  else
    for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
      minimum_phase->log_spectrum[i] = log(spectrum[i]) / 2.0;
  GetMinimumPhaseSpectrum(minimum_phase);

  for (int i = 0; i <= fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] =
      minimum_phase->minimum_phase_spectrum[i][0] *
      forward_real_fft->spectrum[i][0] -
      minimum_phase->minimum_phase_spectrum[i][1] *
      forward_real_fft->spectrum[i][1];
    inverse_real_fft->spectrum[i][1] =
      minimum_phase->minimum_phase_spectrum[i][0] *
      forward_real_fft->spectrum[i][1] +
      minimum_phase->minimum_phase_spectrum[i][1] *
      forward_real_fft->spectrum[i][0];
  }
  fft_execute(inverse_real_fft->inverse_fft);
  fftshift(inverse_real_fft->waveform, fft_size, aperiodic_response);
}

//-----------------------------------------------------------------------------
// RemoveDCComponent()
//-----------------------------------------------------------------------------
static void RemoveDCComponent(const double *periodic_response, int fft_size,
    const double *dc_remover, double *new_periodic_response) {
  double dc_component = 0.0;
  for (int i = fft_size / 2; i < fft_size; ++i)
    dc_component += periodic_response[i];
  for (int i = 0; i < fft_size / 2; ++i)
    new_periodic_response[i] = -dc_component * dc_remover[i];
  for (int i = fft_size / 2; i < fft_size; ++i)
    new_periodic_response[i] -= dc_component * dc_remover[i];
}

//-----------------------------------------------------------------------------
// GetSpectrumWithFractionalTimeShift() calculates a periodic spectrum with
// the fractional time shift under 1/fs.
//-----------------------------------------------------------------------------
static void GetSpectrumWithFractionalTimeShift(int fft_size,
    double coefficient, const InverseRealFFT *inverse_real_fft) {
  double re, im, re2, im2;
  for (int i = 0; i <= fft_size / 2; ++i) {
    re = inverse_real_fft->spectrum[i][0];
    im = inverse_real_fft->spectrum[i][1];
    re2 = cos(coefficient * i);
    im2 = sqrt(1.0 - re2 * re2);  // sin(pshift)

    inverse_real_fft->spectrum[i][0] = re * re2 + im * im2;
    inverse_real_fft->spectrum[i][1] = im * re2 - re * im2;
  }
}

//-----------------------------------------------------------------------------
// GetPeriodicResponse() calculates a periodic response.
//-----------------------------------------------------------------------------
static void GetPeriodicResponse(int fft_size, const double *spectrum,
    const double *aperiodic_ratio, double current_vuv,
    const InverseRealFFT *inverse_real_fft,
    const MinimumPhaseAnalysis *minimum_phase, const double *dc_remover,
    double fractional_time_shift, int fs, double *periodic_response) {
  if (current_vuv <= 0.5 || aperiodic_ratio[0] > 0.999) {
    for (int i = 0; i < fft_size; ++i) periodic_response[i] = 0.0;
    return;
  }

  for (int i = 0; i <= minimum_phase->fft_size / 2; ++i)
    minimum_phase->log_spectrum[i] =
      log(spectrum[i] * (1.0 - aperiodic_ratio[i]) +
      world::kMySafeGuardMinimum) / 2.0;
  GetMinimumPhaseSpectrum(minimum_phase);

  for (int i = 0; i <= fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] =
      minimum_phase->minimum_phase_spectrum[i][0];
    inverse_real_fft->spectrum[i][1] =
      minimum_phase->minimum_phase_spectrum[i][1];
  }

  // apply fractional time delay of fractional_time_shift seconds
  // using linear phase shift
  double coefficient =
    2.0 * world::kPi * fractional_time_shift * fs / fft_size;
  GetSpectrumWithFractionalTimeShift(fft_size, coefficient, inverse_real_fft);

  fft_execute(inverse_real_fft->inverse_fft);
  fftshift(inverse_real_fft->waveform, fft_size, periodic_response);
  RemoveDCComponent(periodic_response, fft_size, dc_remover,
      periodic_response);
}

static void GetSpectralEnvelope(double current_time, double frame_period,
    int f0_length, const double * const *spectrogram, int fft_size,
    double *spectral_envelope) {
  int current_frame_floor = MyMinInt(f0_length - 1,
    static_cast<int>(floor(current_time / frame_period)));
  int current_frame_ceil = MyMinInt(f0_length - 1,
    static_cast<int>(ceil(current_time / frame_period)));
  double interpolation = current_time / frame_period - current_frame_floor;

  if (current_frame_floor == current_frame_ceil)
    for (int i = 0; i <= fft_size / 2; ++i)
      spectral_envelope[i] = fabs(spectrogram[current_frame_floor][i]);
  else
    for (int i = 0; i <= fft_size / 2; ++i)
      spectral_envelope[i] =
        (1.0 - interpolation) * fabs(spectrogram[current_frame_floor][i]) +
        interpolation * fabs(spectrogram[current_frame_ceil][i]);
}

static void GetAperiodicRatio(double current_time, double frame_period,
    int f0_length, const double * const *aperiodicity, int fft_size,
    double *aperiodic_spectrum) {
  int current_frame_floor = MyMinInt(f0_length - 1,
    static_cast<int>(floor(current_time / frame_period)));
  int current_frame_ceil = MyMinInt(f0_length - 1,
    static_cast<int>(ceil(current_time / frame_period)));
  double interpolation = current_time / frame_period - current_frame_floor;

  if (current_frame_floor == current_frame_ceil)
    for (int i = 0; i <= fft_size / 2; ++i)
      aperiodic_spectrum[i] =
        pow(GetSafeAperiodicity(aperiodicity[current_frame_floor][i]), 2.0);
  else
    for (int i = 0; i <= fft_size / 2; ++i)
      aperiodic_spectrum[i] = pow((1.0 - interpolation) *
          GetSafeAperiodicity(aperiodicity[current_frame_floor][i]) +
          interpolation *
          GetSafeAperiodicity(aperiodicity[current_frame_ceil][i]), 2.0);
}

//-----------------------------------------------------------------------------
// GetOneFrameSegment() calculates a periodic and aperiodic response at a time.
//-----------------------------------------------------------------------------
static void GetOneFrameSegment(double current_vuv, int noise_size,
    const double * const *spectrogram, int fft_size,
    const double * const *aperiodicity, int f0_length, double frame_period,
    double current_time, double fractional_time_shift, int fs,
    const ForwardRealFFT *forward_real_fft,
    const InverseRealFFT *inverse_real_fft,
    const MinimumPhaseAnalysis *minimum_phase, const double *dc_remover,
    double *response) {
  double *aperiodic_response = new double[fft_size];
  double *periodic_response = new double[fft_size];

  double *spectral_envelope = new double[fft_size];
  double *aperiodic_ratio = new double[fft_size];
  GetSpectralEnvelope(current_time, frame_period, f0_length, spectrogram,
      fft_size, spectral_envelope);
  GetAperiodicRatio(current_time, frame_period, f0_length, aperiodicity,
      fft_size, aperiodic_ratio);

  // Synthesis of the periodic response
  GetPeriodicResponse(fft_size, spectral_envelope, aperiodic_ratio,
      current_vuv, inverse_real_fft, minimum_phase, dc_remover,
      fractional_time_shift, fs, periodic_response);

  // Synthesis of the aperiodic response
  GetAperiodicResponse(noise_size, fft_size, spectral_envelope,
      aperiodic_ratio, current_vuv, forward_real_fft,
      inverse_real_fft, minimum_phase, aperiodic_response);

  double sqrt_noise_size = sqrt(static_cast<double>(noise_size));
  for (int i = 0; i < fft_size; ++i)
    response[i] =
      (periodic_response[i] * sqrt_noise_size + aperiodic_response[i]) /
      fft_size;

  delete[] spectral_envelope;
  delete[] aperiodic_ratio;
  delete[] periodic_response;
  delete[] aperiodic_response;
}

static void GetTemporalParametersForTimeBase(const double *f0, int f0_length,
    int fs, int y_length, double frame_period, double lowest_f0,
    double *time_axis, double *coarse_time_axis, double *coarse_f0,
    double *coarse_vuv) {
  for (int i = 0; i < y_length; ++i)
    time_axis[i] = i / static_cast<double>(fs);
  // the array 'coarse_time_axis' is supposed to have 'f0_length + 1' positions
  for (int i = 0; i < f0_length; ++i) {
    coarse_time_axis[i] = i * frame_period;
    coarse_f0[i] = f0[i] < lowest_f0 ? 0.0 : f0[i];
    coarse_vuv[i] = coarse_f0[i] == 0.0 ? 0.0 : 1.0;
  }
  coarse_time_axis[f0_length] = f0_length * frame_period;
  coarse_f0[f0_length] = coarse_f0[f0_length - 1] * 2 -
    coarse_f0[f0_length - 2];
  coarse_vuv[f0_length] = coarse_vuv[f0_length - 1] * 2 -
    coarse_vuv[f0_length - 2];
}

static int GetPulseLocationsForTimeBase(const double *interpolated_f0,
    const double *time_axis, int y_length, int fs, double *pulse_locations,
    int *pulse_locations_index, double *pulse_locations_time_shift) {
  double *total_phase = new double[y_length];
  double *wrap_phase = new double[y_length];
  double *wrap_phase_abs = new double[y_length - 1];
  total_phase[0] = 2.0 * world::kPi * interpolated_f0[0] / fs;
  wrap_phase[0] = fmod(total_phase[0], 2.0 * world::kPi);
  for (int i = 1; i < y_length; ++i) {
    total_phase[i] = total_phase[i - 1] +
      2.0 * world::kPi * interpolated_f0[i] / fs;
    wrap_phase[i] = fmod(total_phase[i], 2.0 * world::kPi);
    wrap_phase_abs[i - 1] = fabs(wrap_phase[i] - wrap_phase[i - 1]);
  }

  int number_of_pulses = 0;
  for (int i = 0; i < y_length - 1; ++i) {
    if (wrap_phase_abs[i] > world::kPi) {
      pulse_locations[number_of_pulses] = time_axis[i];
      pulse_locations_index[number_of_pulses] = i;

      // calculate the time shift in seconds between exact fractional pulse
      // position and the integer pulse position (sample i)
      // as we don't have access to the exact pulse position, we infer it
      // from the point between sample i and sample i + 1 where the
      // accummulated phase cross a multiple of 2pi
      // this point is found by solving y1 + x * (y2 - y1) = 0 for x, where y1
      // and y2 are the phases corresponding to sample i and i + 1, offset so
      // they cross zero; x >= 0
      double y1 = wrap_phase[i] - 2.0 * world::kPi;
      double y2 = wrap_phase[i + 1];
      double x = -y1 / (y2 - y1);
      pulse_locations_time_shift[number_of_pulses] = x / fs;

      ++number_of_pulses;
    }
  }

  delete[] wrap_phase_abs;
  delete[] wrap_phase;
  delete[] total_phase;

  return number_of_pulses;
}

static int GetTimeBase(const double *f0, int f0_length, int fs,
    double frame_period, int y_length, double lowest_f0,
    double *pulse_locations, int *pulse_locations_index,
    double *pulse_locations_time_shift, double *interpolated_vuv) {
  double *time_axis = new double[y_length];
  double *coarse_time_axis = new double[f0_length + 1];
  double *coarse_f0 = new double[f0_length + 1];
  double *coarse_vuv = new double[f0_length + 1];
  GetTemporalParametersForTimeBase(f0, f0_length, fs, y_length, frame_period,
      lowest_f0, time_axis, coarse_time_axis, coarse_f0, coarse_vuv);
  double *interpolated_f0 = new double[y_length];
  interp1(coarse_time_axis, coarse_f0, f0_length + 1,
      time_axis, y_length, interpolated_f0);
  interp1(coarse_time_axis, coarse_vuv, f0_length + 1,
      time_axis, y_length, interpolated_vuv);

  for (int i = 0; i < y_length; ++i) {
    interpolated_vuv[i] = interpolated_vuv[i] > 0.5 ? 1.0 : 0.0;
    interpolated_f0[i] =
      interpolated_vuv[i] == 0.0 ? world::kDefaultF0 : interpolated_f0[i];
  }

  int number_of_pulses = GetPulseLocationsForTimeBase(interpolated_f0,
      time_axis, y_length, fs, pulse_locations, pulse_locations_index,
      pulse_locations_time_shift);

  delete[] coarse_vuv;
  delete[] coarse_f0;
  delete[] coarse_time_axis;
  delete[] time_axis;
  delete[] interpolated_f0;

  return number_of_pulses;
}

static void GetDCRemover(int fft_size, double *dc_remover) {
  double dc_component = 0.0;
  for (int i = 0; i < fft_size / 2; ++i) {
    dc_remover[i] = 0.5 -
      0.5 * cos(2.0 * world::kPi * (i + 1.0) / (1.0 + fft_size));
    dc_remover[fft_size - i - 1] = dc_remover[i];
    dc_component += dc_remover[i] * 2.0;
  }
  for (int i = 0; i < fft_size / 2; ++i) {
    dc_remover[i] /= dc_component;
    dc_remover[fft_size - i - 1] = dc_remover[i];
  }
}

}  // namespace

void Synthesis(const double *f0, int f0_length,
    const double * const *spectrogram, const double * const *aperiodicity,
    int fft_size, double frame_period, int fs, int y_length, double *y) {
  randn_reseed();

  double *impulse_response = new double[fft_size];

  for (int i = 0; i < y_length; ++i) y[i] = 0.0;

  MinimumPhaseAnalysis minimum_phase = {0};
  InitializeMinimumPhaseAnalysis(fft_size, &minimum_phase);
  InverseRealFFT inverse_real_fft = {0};
  InitializeInverseRealFFT(fft_size, &inverse_real_fft);
  ForwardRealFFT forward_real_fft = {0};
  InitializeForwardRealFFT(fft_size, &forward_real_fft);

  double *pulse_locations = new double[y_length];
  int *pulse_locations_index = new int[y_length];
  double *pulse_locations_time_shift = new double[y_length];
  double *interpolated_vuv = new double[y_length];
  int number_of_pulses = GetTimeBase(f0, f0_length, fs, frame_period / 1000.0,
      y_length, fs / fft_size + 1.0, pulse_locations, pulse_locations_index,
      pulse_locations_time_shift, interpolated_vuv);

  double *dc_remover = new double[fft_size];
  GetDCRemover(fft_size, dc_remover);

  frame_period /= 1000.0;
  int noise_size;
  int index, offset, lower_limit, upper_limit;
  for (int i = 0; i < number_of_pulses; ++i) {
    noise_size = pulse_locations_index[MyMinInt(number_of_pulses - 1, i + 1)] -
      pulse_locations_index[i];

    GetOneFrameSegment(interpolated_vuv[pulse_locations_index[i]], noise_size,
        spectrogram, fft_size, aperiodicity, f0_length, frame_period,
        pulse_locations[i], pulse_locations_time_shift[i], fs,
        &forward_real_fft, &inverse_real_fft, &minimum_phase, dc_remover,
        impulse_response);
    offset = pulse_locations_index[i] - fft_size / 2 + 1;
    lower_limit = MyMaxInt(0, -offset);
    upper_limit = MyMinInt(fft_size, y_length - offset);
    for (int j = lower_limit; j < upper_limit; ++j) {
      index = j + offset;
      y[index] += impulse_response[j];
    }
  }

  delete[] dc_remover;
  delete[] pulse_locations;
  delete[] pulse_locations_index;
  delete[] pulse_locations_time_shift;
  delete[] interpolated_vuv;

  DestroyMinimumPhaseAnalysis(&minimum_phase);
  DestroyInverseRealFFT(&inverse_real_fft);
  DestroyForwardRealFFT(&forward_real_fft);

  delete[] impulse_response;
}
