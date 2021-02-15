//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//-----------------------------------------------------------------------------
#ifndef WORLD_SYNTHESISREALTIME_H_
#define WORLD_SYNTHESISREALTIME_H_

#include "world/common.h"
#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

//-----------------------------------------------------------------------------
// A struct for real-time synthesis.
// I use no class for compatibility in C language.
// Please make a class for encapsulating it as needed.
// This synthesizer uses a ring buffer.
//-----------------------------------------------------------------------------
typedef struct {
  // Basic parameters
  int fs;
  double frame_period;
  int buffer_size;
  int number_of_pointers;
  int fft_size;

  // Sound buffer for output. The length is buffer_size [sample].
  double *buffer;
  int current_pointer;
  int i;

  // For DC removal
  double *dc_remover;

  //---------------------------------------------------------------------------
  // Followings are internal parameters.
  // You should not modify them if you are not expert.

  // Speech parameters in each pointer.
  int *f0_length;
  int *f0_origin;
  double ***spectrogram;
  double ***aperiodicity;


  // Note:
  // This is an extremely rough implementation.
  // I should optimize this implementation.
  int current_pointer2;
  int head_pointer;
  int synthesized_sample;

  // Internal parameters.
  int handoff;
  double handoff_phase;
  double handoff_f0;
  int last_location;

  int cumulative_frame;
  int current_frame;

  double **interpolated_vuv;
  double **pulse_locations;
  int **pulse_locations_index;
  int *number_of_pulses;

  double *impulse_response;

  // FFT
  MinimumPhaseAnalysis minimum_phase;
  InverseRealFFT inverse_real_fft;
  ForwardRealFFT forward_real_fft;
} WorldSynthesizer;

//-----------------------------------------------------------------------------
// InitializeSynthesizer() initializes the synthesizer based on basic
// parameters.
//
// Input:
//   fs                   : Sampling frequency
//   frame_period         : Frame period (ms)
//   fft_size             : FFT size
//   buffer_size          : Buffer size (sample)
//   number_of_pointers   : The number of elements in the ring buffer
//
// Output:
//   synth                : Initialized synthesizer
//-----------------------------------------------------------------------------
void InitializeSynthesizer(int fs, double frame_period, int fft_size,
  int buffer_size, int number_of_pointers, WorldSynthesizer *synth);

//-----------------------------------------------------------------------------
// AddParameters() attempts to add speech parameters.
// You can add several frames at the same time.
//
// Input:
//   f0                   : F0 contour with length of f0_length
//   f0_length            : This is associated with the number of frames
//   spectrogram          : Spectrogram
//   aperiodicity         : Aperiodicity
//
// Output:
//   synth                : Synthesizer
//
// Return value:
//   1: True, 0: False.
//-----------------------------------------------------------------------------
int AddParameters(double *f0, int f0_length, double **spectrogram,
  double **aperiodicity, WorldSynthesizer *synth);

//-----------------------------------------------------------------------------
// RefreshSynthesizer() sets the parameters to default.
//-----------------------------------------------------------------------------
void RefreshSynthesizer(WorldSynthesizer *synth);

//-----------------------------------------------------------------------------
// DestroySynthesizer() release the memory.
//-----------------------------------------------------------------------------
void DestroySynthesizer(WorldSynthesizer *synth);

//-----------------------------------------------------------------------------
// IsLocked() checks whether the synthesizer is locked or not.
// "Lock" is defined as the situation that the ring buffer cannot add
// parameters and cannot synthesize the waveform.
// It will be caused when the duration calculated by the number of added frames
// is below 1 / F0 + buffer_size / fs.
// If this function returns True, please refresh the synthesizer.
//
// Input:
//   Synth            : Synthesizer (pointer)
//
// Output:
//   1: True, 0: False.
//-----------------------------------------------------------------------------
int IsLocked(WorldSynthesizer *synth);

//-----------------------------------------------------------------------------
// Synthesis2() generates speech with length of synth->buffer_size sample.
// The parameters are automatially updated, and memory is also released.
//
// Input:
//   Synth            : Synthesizer (pointer)
//
// Output:
//   1: True, 0: False.
//-----------------------------------------------------------------------------
int Synthesis2(WorldSynthesizer *synth);

WORLD_END_C_DECLS

#endif  // WORLD_SYNTHESISREALTIME_H_
