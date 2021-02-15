//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//
// Voice synthesis based on f0, spectrogram and aperiodicity.
// This is an implementation for real-time applications.
// Note: Several functions are same as those of synthesis.cpp.
//
// Caution: This is an implementation as a prototype.
//          Specifications may change. There may be a bug.
//
// Caution: DC removal was re-implemented. However, this implementation is
//          different from the implementation in synthesis.cpp. The sound
//          quality is almost all the same in both implementations.
//-----------------------------------------------------------------------------
#include "world/synthesisrealtime.h"

#include <math.h>
#include <stdlib.h>

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
        log(spectrum[i] * aperiodic_ratio[i] +
        world::kMySafeGuardMinimum) / 2.0;
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

static void ClearRingBuffer(int start, int end, WorldSynthesizer *synth) {
  int pointer;
  for (int i = start; i < end; ++i) {
    pointer = i % synth->number_of_pointers;
    synth->number_of_pulses[pointer] = 0;
    if (synth->pulse_locations[pointer] != NULL) {
      delete[] synth->pulse_locations[pointer];
      synth->pulse_locations[pointer] = NULL;
    }
    if (synth->interpolated_vuv[pointer] != NULL) {
      delete[] synth->interpolated_vuv[pointer];
      synth->interpolated_vuv[pointer] = NULL;
    }
    if (synth->pulse_locations_index[pointer] != NULL) {
      delete[] synth->pulse_locations_index[pointer];
      synth->pulse_locations_index[pointer] = NULL;
    }
  }
}

static int SeekSynthesizer(double current_location, WorldSynthesizer *synth) {
  int frame_number = static_cast<int>(current_location / synth->frame_period);

  int tmp_pointer = synth->current_pointer2;
  int tmp;
  for (int i = 0; i < synth->head_pointer - synth->current_pointer2; ++i) {
    tmp = (tmp_pointer + i) % synth->number_of_pointers;
    if (synth->f0_origin[tmp] <= frame_number &&
        frame_number < synth->f0_origin[tmp] + synth->f0_length[tmp]) {
      tmp_pointer += i;
      break;
    }
  }
  ClearRingBuffer(synth->current_pointer2, tmp_pointer, synth);
  synth->current_pointer2 = tmp_pointer;
  return 1;
}

static void SearchPointer(int frame,  WorldSynthesizer *synth, int flag,
    double **front, double **next) {
  int pointer = synth->current_pointer2 % synth->number_of_pointers;
  int index = -1;
  for (int i = 0; i < synth->f0_length[pointer]; ++i)
    if (synth->f0_origin[pointer] + i == frame) {
      index = i;
      break;
    }

  double ***tmp_pointer =
    flag == 0 ? synth->spectrogram : synth->aperiodicity;

  *front = tmp_pointer[pointer][index];
  *next = index == synth->f0_length[pointer] - 1 ?
    tmp_pointer[(synth->current_pointer2 + 1) %
    synth->number_of_pointers][0] : tmp_pointer[pointer][index + 1];
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
    new_periodic_response[i] = 0.0;
  for (int i = fft_size / 2; i < fft_size; ++i)
    new_periodic_response[i] -= dc_component * dc_remover[i - fft_size / 2];
}

//-----------------------------------------------------------------------------
// GetPeriodicResponse() calculates an aperiodic response.
//-----------------------------------------------------------------------------
static void GetPeriodicResponse(int fft_size, const double *spectrum,
    const double *aperiodic_ratio, double current_vuv,
    const InverseRealFFT *inverse_real_fft,
    const MinimumPhaseAnalysis *minimum_phase,
    const double *dc_remover, double *periodic_response) {
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

  fft_execute(inverse_real_fft->inverse_fft);
  fftshift(inverse_real_fft->waveform, fft_size, periodic_response);
  RemoveDCComponent(periodic_response, fft_size, dc_remover,
      periodic_response);
}

static void GetSpectralEnvelope(double current_location,
    WorldSynthesizer *synth, double *spectral_envelope) {
  int current_frame_floor =
    static_cast<int>(current_location / synth->frame_period);

  int current_frame_ceil =
    static_cast<int>(ceil(current_location / synth->frame_period));
  double interpolation =
    current_location / synth->frame_period - current_frame_floor;
  double *front = NULL;
  double *next = NULL;
  SearchPointer(current_frame_floor, synth, 0, &front, &next);

  if (current_frame_floor == current_frame_ceil)
    for (int i = 0; i <= synth->fft_size / 2; ++i)
      spectral_envelope[i] = fabs(front[i]);
  else
    for (int i = 0; i <= synth->fft_size / 2; ++i)
      spectral_envelope[i] =
      (1.0 - interpolation) * fabs(front[i]) + interpolation * fabs(next[i]);
}

static void GetAperiodicRatio(double current_location,
    WorldSynthesizer *synth, double *aperiodic_spectrum) {
  int current_frame_floor =
    static_cast<int>(current_location / synth->frame_period);

  int current_frame_ceil =
    static_cast<int>(ceil(current_location / synth->frame_period));
  double interpolation =
    current_location / synth->frame_period - current_frame_floor;

  double *front = NULL;
  double *next = NULL;
  SearchPointer(current_frame_floor, synth, 1, &front, &next);

  if (current_frame_floor == current_frame_ceil)
    for (int i = 0; i <= synth->fft_size / 2; ++i)
      aperiodic_spectrum[i] = pow(GetSafeAperiodicity(front[i]), 2.0);
  else
    for (int i = 0; i <= synth->fft_size / 2; ++i)
      aperiodic_spectrum[i] =
        pow((1.0 - interpolation) * GetSafeAperiodicity(front[i]) +
        interpolation * GetSafeAperiodicity(next[i]), 2.0);
}

static double GetCurrentVUV(int current_location, WorldSynthesizer *synth) {
  double current_vuv = 0.0;
  int pointer = synth->current_pointer % synth->number_of_pointers;

  int start_sample = MyMaxInt(0,
    static_cast<int>(ceil((synth->f0_origin[pointer] - 1) *
    synth->frame_period * synth->fs)));

  current_vuv =
    synth->interpolated_vuv[pointer][current_location - start_sample + 1];
  return current_vuv;
}

//-----------------------------------------------------------------------------
// GetOneFrameSegment() calculates a periodic and aperiodic response at a time.
//-----------------------------------------------------------------------------
static void GetOneFrameSegment(int noise_size, int current_location,
    WorldSynthesizer *synth) {
  double *aperiodic_response = new double[synth->fft_size];
  double *periodic_response = new double[synth->fft_size];
  double *spectral_envelope = new double[synth->fft_size];
  double *aperiodic_ratio = new double[synth->fft_size];

  double tmp_location = static_cast<double>(current_location) / synth->fs;
  SeekSynthesizer(tmp_location, synth);
  GetSpectralEnvelope(tmp_location, synth, spectral_envelope);
  GetAperiodicRatio(tmp_location, synth, aperiodic_ratio);

  double current_vuv = GetCurrentVUV(current_location, synth);

  // Synthesis of the periodic response
  GetPeriodicResponse(synth->fft_size, spectral_envelope, aperiodic_ratio,
      current_vuv, &synth->inverse_real_fft, &synth->minimum_phase,
      synth->dc_remover, periodic_response);

  // Synthesis of the aperiodic response
  GetAperiodicResponse(noise_size, synth->fft_size, spectral_envelope,
      aperiodic_ratio, current_vuv, &synth->forward_real_fft,
      &synth->inverse_real_fft, &synth->minimum_phase, aperiodic_response);

  double sqrt_noise_size = sqrt(static_cast<double>(noise_size));
  for (int i = 0; i < synth->fft_size; ++i)
    synth->impulse_response[i] =
    (periodic_response[i] * sqrt_noise_size + aperiodic_response[i]) /
      synth->fft_size;

  delete[] spectral_envelope;
  delete[] aperiodic_ratio;
  delete[] periodic_response;
  delete[] aperiodic_response;
}

static void GetTemporalParametersForTimeBase(const double *f0, int f0_length,
    WorldSynthesizer *synth, double *coarse_time_axis, double *coarse_f0,
    double *coarse_vuv) {
  int cumulative_frame = MyMaxInt(0, synth->cumulative_frame - f0_length);
  coarse_f0[0] = synth->handoff_f0;
  coarse_time_axis[0] = cumulative_frame * synth->frame_period;
  coarse_vuv[0] = synth->handoff_f0 == 0 ? 0.0 : 1.0;
  for (int i = 0; i < f0_length; ++i) {
    coarse_time_axis[i + synth->handoff] =
      (i + cumulative_frame + synth->handoff) * synth->frame_period;
    coarse_f0[i + synth->handoff] = f0[i];
    coarse_vuv[i + synth->handoff] = f0[i] == 0.0 ? 0.0 : 1.0;
  }
}

static void GetPulseLocationsForTimeBase(const double *interpolated_f0,
    const double *time_axis, int number_of_samples, double origin,
    WorldSynthesizer *synth) {
  double *total_phase = new double[number_of_samples + synth->handoff];
  total_phase[0] = synth->handoff == 1 ? synth->handoff_phase :
    2.0 * world::kPi * interpolated_f0[0] / synth->fs;

  total_phase[1] = total_phase[0] + 2.0 * world::kPi *
    interpolated_f0[0] / synth->fs;
  for (int i = 1 + synth->handoff; i < number_of_samples + synth->handoff; ++i)
    total_phase[i] = total_phase[i - 1] +
      2.0 * world::kPi * interpolated_f0[i - synth->handoff] / synth->fs;
  synth->handoff_phase = total_phase[number_of_samples - 1 + synth->handoff];

  double *wrap_phase = new double[number_of_samples + synth->handoff];
  for (int i = 0; i < number_of_samples + synth->handoff; ++i)
    wrap_phase[i] = fmod(total_phase[i], 2.0 * world::kPi);

  double *wrap_phase_abs = new double[number_of_samples + synth->handoff];
  for (int i = 0; i < number_of_samples - 1 + synth->handoff; ++i)
    wrap_phase_abs[i] = fabs(wrap_phase[i + 1] - wrap_phase[i]);

  int pointer = synth->head_pointer % synth->number_of_pointers;
  int number_of_pulses = 0;
  for (int i = 0; i < number_of_samples - 1 + synth->handoff; ++i)
    if (wrap_phase_abs[i] > world::kPi) {
      synth->pulse_locations[pointer][number_of_pulses] =
        time_axis[i] - static_cast<double>(synth->handoff) / synth->fs;
      synth->pulse_locations_index[pointer][number_of_pulses] =
        matlab_round(synth->pulse_locations[pointer][number_of_pulses] *
            synth->fs);
      ++number_of_pulses;
    }
  synth->number_of_pulses[pointer] = number_of_pulses;

  if (number_of_pulses != 0)
    synth->last_location =
      synth->pulse_locations_index[pointer][number_of_pulses - 1];
  delete[] wrap_phase_abs;
  delete[] wrap_phase;
  delete[] total_phase;
}

static void GetTimeBase(const double *f0, int f0_length, int start_sample,
    int number_of_samples, WorldSynthesizer *synth) {
  double *coarse_time_axis = new double[f0_length + synth->handoff];
  double *coarse_f0 = new double[f0_length + synth->handoff];
  double *coarse_vuv = new double[f0_length + synth->handoff];

  GetTemporalParametersForTimeBase(f0, f0_length, synth,
      coarse_time_axis, coarse_f0, coarse_vuv);

  double *interpolated_f0 = new double[number_of_samples];
  double *time_axis = new double[number_of_samples];

  for (int i = 0; i < number_of_samples; ++i)
    time_axis[i] = (i + start_sample) / static_cast<double>(synth->fs);

  int pointer = synth->head_pointer % synth->number_of_pointers;
  interp1(coarse_time_axis, coarse_f0, f0_length + synth->handoff, time_axis,
    number_of_samples, interpolated_f0);
  interp1(coarse_time_axis, coarse_vuv, f0_length + synth->handoff, time_axis,
    number_of_samples, synth->interpolated_vuv[pointer]);
  for (int i = 0; i < number_of_samples; ++i) {
    synth->interpolated_vuv[pointer][i] =
      synth->interpolated_vuv[pointer][i] > 0.5 ? 1.0 : 0.0;
    interpolated_f0[i] =
      synth->interpolated_vuv[pointer][i] == 0.0 ?
      world::kDefaultF0 : interpolated_f0[i];
  }

  GetPulseLocationsForTimeBase(interpolated_f0, time_axis, number_of_samples,
      coarse_time_axis[0], synth);

  synth->handoff_f0 = interpolated_f0[number_of_samples - 1];
  delete[] time_axis;
  delete[] interpolated_f0;
  delete[] coarse_time_axis;
  delete[] coarse_f0;
  delete[] coarse_vuv;
}

static int GetNextPulseLocationIndex(WorldSynthesizer *synth) {
  int pointer = synth->current_pointer % synth->number_of_pointers;
  if (synth->i < synth->number_of_pulses[pointer] - 1)
    return synth->pulse_locations_index[pointer][synth->i + 1];
  else if (synth->current_pointer == synth->head_pointer - 1)
    return 0;

  for (int i = 1; i < synth->number_of_pointers; ++i) {
    pointer = (i + synth->current_pointer) % synth->number_of_pointers;
    if (synth->number_of_pulses[pointer] != 0)
      return synth->pulse_locations_index[pointer][0];
  }
  return 0;
}

static int UpdateSynthesizer(int current_location, WorldSynthesizer *synth) {
  int pointer = synth->current_pointer % synth->number_of_pointers;
  if (synth->i < synth->number_of_pulses[pointer] - 1) {
    synth->i++;
    return 1;
  } else {
    if (synth->current_pointer == synth->head_pointer - 1) return 0;
  }

  for (int i = 1; i < synth->number_of_pointers; ++i) {
    pointer = (i + synth->current_pointer) % synth->number_of_pointers;
    if (synth->number_of_pulses[pointer] != 0) {
      synth->i = 0;
      synth->current_pointer += i;
      return 1;
    }
  }
  return 0;
}

static int CheckSynthesizer(WorldSynthesizer *synth) {
  if (synth->synthesized_sample + synth->buffer_size >= synth->last_location)
    return 0;

  int pointer = synth->current_pointer % synth->number_of_pointers;
  while (synth->number_of_pulses[pointer] == 0) {
    if (synth->current_pointer == synth->head_pointer) break;
    synth->current_pointer++;
    pointer = synth->current_pointer % synth->number_of_pointers;
  }
  return 1;
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

void InitializeSynthesizer(int fs, double frame_period, int fft_size,
    int buffer_size, int number_of_pointers, WorldSynthesizer *synth) {
  // Set basic parameters
  synth->fs = fs;
  synth->frame_period = frame_period / 1000.0;
  synth->buffer_size = buffer_size;
  synth->number_of_pointers = number_of_pointers;
  synth->fft_size = fft_size;

  // Memory allocations
  synth->f0_length = new int[number_of_pointers];
  synth->spectrogram = new double**[number_of_pointers];
  synth->aperiodicity = new double**[number_of_pointers];
  synth->interpolated_vuv = new double*[number_of_pointers];
  synth->pulse_locations = new double*[number_of_pointers];
  synth->pulse_locations_index = new int*[number_of_pointers];
  synth->number_of_pulses = new int[number_of_pointers];
  synth->f0_origin = new int[number_of_pointers];
  for (int i = 0; i < synth->number_of_pointers; ++i) {
    synth->interpolated_vuv[i] = NULL;
    synth->pulse_locations[i] = NULL;
    synth->pulse_locations_index[i] = NULL;
  }

  synth->buffer = new double[buffer_size * 2 + fft_size];
  synth->impulse_response = new double[synth->fft_size];
  synth->dc_remover = new double[synth->fft_size / 2];

  // Initilize internal parameters
  RefreshSynthesizer(synth);

  InitializeMinimumPhaseAnalysis(fft_size, &synth->minimum_phase);
  InitializeInverseRealFFT(fft_size, &synth->inverse_real_fft);
  InitializeForwardRealFFT(fft_size, &synth->forward_real_fft);
}

int AddParameters(double *f0, int f0_length, double **spectrogram,
    double **aperiodicity, WorldSynthesizer *synth) {
  if (synth->head_pointer - synth->current_pointer2 ==
    synth->number_of_pointers)
    return 0;  // Since the queue is full, we cannot add the parameters.
  int pointer = synth->head_pointer % synth->number_of_pointers;
  synth->f0_length[pointer] = f0_length;
  synth->f0_origin[pointer] = synth->cumulative_frame + 1;
  synth->cumulative_frame += f0_length;

  synth->spectrogram[pointer] = spectrogram;
  synth->aperiodicity[pointer] = aperiodicity;
  if (synth->cumulative_frame < 1) {
    synth->handoff_f0 = f0[f0_length - 1];
    synth->number_of_pulses[pointer] = 0;
    synth->head_pointer++;
    synth->handoff = 1;
    return 1;
  }

  int start_sample =
    MyMaxInt(0, static_cast<int>(ceil((synth->cumulative_frame - f0_length) *
      synth->frame_period * synth->fs)));
  int end_sample =
    static_cast<int>(ceil((synth->cumulative_frame) *
      synth->frame_period * synth->fs));
  int number_of_samples = end_sample - start_sample;

  // Memory allocation
  synth->interpolated_vuv[pointer] = new double[number_of_samples + 1];
  synth->pulse_locations[pointer] = new double[number_of_samples];
  synth->pulse_locations_index[pointer] = new int[number_of_samples];

  GetTimeBase(f0, f0_length, start_sample, number_of_samples, synth);

  synth->handoff_f0 = f0[f0_length - 1];
  synth->head_pointer++;
  synth->handoff = 1;
  return 1;
}

void RefreshSynthesizer(WorldSynthesizer *synth) {
  ClearRingBuffer(0, synth->number_of_pointers, synth);
  synth->handoff_phase = 0;
  synth->handoff_f0 = 0;
  synth->cumulative_frame = -1;
  synth->last_location = 0;

  synth->current_pointer = 0;
  synth->current_pointer2 = 0;
  synth->head_pointer = 0;
  synth->handoff = 0;

  synth->i = 0;
  synth->current_frame = 0;

  synth->synthesized_sample = 0;

  for (int i = 0; i < synth->buffer_size * 2 + synth->fft_size; ++i)
    synth->buffer[i] = 0;
  GetDCRemover(synth->fft_size / 2, synth->dc_remover);
}

void DestroySynthesizer(WorldSynthesizer *synth) {
  RefreshSynthesizer(synth);
  delete[] synth->f0_length;

  delete[] synth->spectrogram;
  delete[] synth->aperiodicity;

  delete[] synth->buffer;
  delete[] synth->impulse_response;
  delete[] synth->dc_remover;

  delete[] synth->interpolated_vuv;
  delete[] synth->pulse_locations;
  delete[] synth->pulse_locations_index;
  delete[] synth->number_of_pulses;
  delete[] synth->f0_origin;

  DestroyMinimumPhaseAnalysis(&synth->minimum_phase);
  DestroyInverseRealFFT(&synth->inverse_real_fft);
  DestroyForwardRealFFT(&synth->forward_real_fft);
}

int IsLocked(WorldSynthesizer *synth) {
  int judge = 0;
  if (synth->head_pointer - synth->current_pointer2 ==
      synth->number_of_pointers)
    judge++;
  if (synth->synthesized_sample + synth->buffer_size >= synth->last_location)
    judge++;

  return judge == 2 ? 1 : 0;
}

int Synthesis2(WorldSynthesizer *synth) {
  if (CheckSynthesizer(synth) == 0)
    return 0;
  for (int i = 0; i < synth->buffer_size + synth->fft_size; ++i)
    synth->buffer[i] = synth->buffer[i + synth->buffer_size];

  int pointer = synth->current_pointer % synth->number_of_pointers;
  int noise_size, offset, tmp, index;
  int current_location = synth->pulse_locations_index[pointer][synth->i];
  while (current_location < synth->synthesized_sample + synth->buffer_size) {
    tmp = GetNextPulseLocationIndex(synth);
    noise_size = tmp - current_location;

    GetOneFrameSegment(noise_size, current_location, synth);
    offset =
      current_location - synth->synthesized_sample - synth->fft_size / 2 + 1;
    for (int i = MyMaxInt(0, -offset); i < synth->fft_size; ++i) {
      index = i + offset;
      synth->buffer[index] += synth->impulse_response[i];
    }
    current_location = tmp;
    UpdateSynthesizer(current_location, synth);
  }
  synth->synthesized_sample += synth->buffer_size;
  SeekSynthesizer(synth->synthesized_sample, synth);
  return 1;
}
