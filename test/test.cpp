//-----------------------------------------------------------------------------
// Copyright 2012-2016 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// Test program for WORLD 0.1.2 (2012/08/19)
// Test program for WORLD 0.1.3 (2013/07/26)
// Test program for WORLD 0.1.4 (2014/04/29)
// Test program for WORLD 0.1.4_3 (2015/03/07)
// Test program for WORLD 0.2.0 (2015/05/29)
// Test program for WORLD 0.2.0_1 (2015/05/31)
// Test program for WORLD 0.2.0_2 (2015/06/06)
// Test program for WORLD 0.2.0_3 (2015/07/28)
// Test program for WORLD 0.2.0_4 (2015/11/15)
// Test program for WORLD in GitHub (2015/11/16-)
// Latest update: 2016/03/04

// test.exe input.wav outout.wav f0 spec
// input.wav  : Input file
// output.wav : Output file
// f0         : F0 scaling (a positive number)
// spec       : Formant scaling (a positive number)
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#if (defined (__WIN32__) || defined (_WIN32)) && !defined (__MINGW32__)
#include <conio.h>
#include <windows.h>
#pragma comment(lib, "winmm.lib")
#pragma warning(disable : 4996)
#endif
#if (defined (__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
#include <stdint.h>
#include <sys/time.h>
#endif

// For .wav input/output functions.
#include "audioio.h"

// WORLD core functions.
// Note: win.sln uses an option in Additional Include Directories.
// To compile the program, the option "-I $(SolutionDir)..\src" was set.
#include "world/d4c.h"
#include "world/dio.h"
#include "world/matlabfunctions.h"
#include "world/cheaptrick.h"
#include "world/stonemask.h"
#include "world/synthesis.h"

#if (defined (__linux__) || defined(__CYGWIN__) || defined(__APPLE__))
// Linux porting section: implement timeGetTime() by gettimeofday(),
#ifndef DWORD
#define DWORD uint32_t
#endif
DWORD timeGetTime() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  DWORD ret = static_cast<DWORD>(tv.tv_usec / 1000 + tv.tv_sec * 1000);
  return ret;
}
#endif

//-----------------------------------------------------------------------------
// struct for WORLD
// This struct is an option.
// Users are NOT forced to use this struct.
//-----------------------------------------------------------------------------
typedef struct {
  double frame_period;
  int fs;

  double *f0;
  double *time_axis;
  int f0_length;

  double **spectrogram;
  double **aperiodicity;
  int fft_size;
} WorldParameters;

namespace {

void DisplayInformation(int fs, int nbit, int x_length) {
  printf("File information\n");
  printf("Sampling : %d Hz %d Bit\n", fs, nbit);
  printf("Length %d [sample]\n", x_length);
  printf("Length %f [sec]\n", static_cast<double>(x_length) / fs);
}

void F0Estimation(double *x, int x_length, WorldParameters *world_parameters) {
  DioOption option = {0};
  InitializeDioOption(&option);

  // Modification of the option
  // When you You must set the same value.
  // If a different value is used, you may suffer a fatal error because of a
  // illegal memory access.
  option.frame_period = world_parameters->frame_period;

  // Valuable option.speed represents the ratio for downsampling.
  // The signal is downsampled to fs / speed Hz.
  // If you want to obtain the accurate result, speed should be set to 1.
  option.speed = 1;

  // You should not set option.f0_floor to under world::kFloorF0.
  // If you want to analyze such low F0 speech, please change world::kFloorF0.
  // Processing speed may sacrify, provided that the FFT length changes.
  option.f0_floor = 71.0;

  // You can give a positive real number as the threshold.
  // Most strict value is 0, but almost all results are counted as unvoiced.
  // The value from 0.02 to 0.2 would be reasonable.
  option.allowed_range = 0.1;

  // Parameters setting and memory allocation.
  world_parameters->f0_length = GetSamplesForDIO(world_parameters->fs,
    x_length, world_parameters->frame_period);
  world_parameters->f0 = new double[world_parameters->f0_length];
  world_parameters->time_axis = new double[world_parameters->f0_length];
  double *refined_f0 = new double[world_parameters->f0_length];

  printf("\nAnalysis\n");
  DWORD elapsed_time = timeGetTime();
  Dio(x, x_length, world_parameters->fs, &option, world_parameters->time_axis,
      world_parameters->f0);
  printf("DIO: %d [msec]\n", timeGetTime() - elapsed_time);

  // StoneMask is carried out to improve the estimation performance.
  elapsed_time = timeGetTime();
  StoneMask(x, x_length, world_parameters->fs, world_parameters->time_axis,
      world_parameters->f0, world_parameters->f0_length, refined_f0);
  printf("StoneMask: %d [msec]\n", timeGetTime() - elapsed_time);

  for (int i = 0; i < world_parameters->f0_length; ++i)
    world_parameters->f0[i] = refined_f0[i];

  delete[] refined_f0;
  return;
}

void SpectralEnvelopeEstimation(double *x, int x_length,
    WorldParameters *world_parameters) {
  CheapTrickOption option = {0};
  InitializeCheapTrickOption(&option);

  // This value may be better one for HMM speech synthesis.
  // Default value is -0.09.
  option.q1 = -0.15;

  // Important notice (2016/02/02)
  // You can control a parameter used for the lowest F0 in speech.
  // You must not set the f0_floor to 0.
  // It will cause a fatal error because fft_size indicates the infinity.
  // You must not change the f0_floor after memory allocation.
  // You should check the fft_size before excucing the analysis/synthesis.
  // The default value (71.0) is strongly recommended.
  // On the other hand, setting the lowest F0 of speech is a good choice
  // to reduce the fft_size.
  option.f0_floor = 71.0;

  // Parameters setting and memory allocation.
  world_parameters->fft_size =
    GetFFTSizeForCheapTrick(world_parameters->fs, &option);
  world_parameters->spectrogram = new double *[world_parameters->f0_length];
  for (int i = 0; i < world_parameters->f0_length; ++i) {
    world_parameters->spectrogram[i] =
      new double[world_parameters->fft_size / 2 + 1];
  }

  DWORD elapsed_time = timeGetTime();
  CheapTrick(x, x_length, world_parameters->fs, world_parameters->time_axis,
      world_parameters->f0, world_parameters->f0_length, &option,
      world_parameters->spectrogram);
  printf("CheapTrick: %d [msec]\n", timeGetTime() - elapsed_time);
}

void AperiodicityEstimation(double *x, int x_length,
    WorldParameters *world_parameters) {
  D4COption option = {0};
  InitializeD4COption(&option);

  // Parameters setting and memory allocation.
  world_parameters->aperiodicity = new double *[world_parameters->f0_length];
  for (int i = 0; i < world_parameters->f0_length; ++i) {
    world_parameters->aperiodicity[i] =
      new double[world_parameters->fft_size / 2 + 1];
  }

  DWORD elapsed_time = timeGetTime();
  // option is not implemented in this version. This is for future update.
  // We can use "NULL" as the argument.
  D4C(x, x_length, world_parameters->fs, world_parameters->time_axis,
      world_parameters->f0, world_parameters->f0_length,
      world_parameters->fft_size, &option, world_parameters->aperiodicity);
  printf("D4C: %d [msec]\n", timeGetTime() - elapsed_time);
}

void ParameterModification(int argc, char *argv[], int fs, int f0_length,
    int fft_size, double *f0, double **spectrogram) {
  // F0 scaling
  if (argc >= 4) {
    double shift = atof(argv[3]);
    for (int i = 0; i < f0_length; ++i) f0[i] *= shift;
  }
  if (argc < 5) return;

  // Spectral stretching
  double ratio = atof(argv[4]);
  double *freq_axis1 = new double[fft_size];
  double *freq_axis2 = new double[fft_size];
  double *spectrum1 = new double[fft_size];
  double *spectrum2 = new double[fft_size];

  for (int i = 0; i <= fft_size / 2; ++i) {
    freq_axis1[i] = ratio * i / fft_size * fs;
    freq_axis2[i] = static_cast<double>(i) / fft_size * fs;
  }
  for (int i = 0; i < f0_length; ++i) {
    for (int j = 0; j <= fft_size / 2; ++j)
      spectrum1[j] = log(spectrogram[i][j]);
    interp1(freq_axis1, spectrum1, fft_size / 2 + 1, freq_axis2,
      fft_size / 2 + 1, spectrum2);
    for (int j = 0; j <= fft_size / 2; ++j)
      spectrogram[i][j] = exp(spectrum2[j]);
    if (ratio >= 1.0) continue;
    for (int j = static_cast<int>(fft_size / 2.0 * ratio);
        j <= fft_size / 2; ++j)
      spectrogram[i][j] =
      spectrogram[i][static_cast<int>(fft_size / 2.0 * ratio) - 1];
  }
  delete[] spectrum1;
  delete[] spectrum2;
  delete[] freq_axis1;
  delete[] freq_axis2;
}

void WaveformSynthesis(WorldParameters *world_parameters, int fs,
    int y_length, double *y) {
  DWORD elapsed_time;
  // Synthesis by the aperiodicity
  printf("\nSynthesis\n");
  elapsed_time = timeGetTime();
  Synthesis(world_parameters->f0, world_parameters->f0_length,
      world_parameters->spectrogram, world_parameters->aperiodicity,
      world_parameters->fft_size, world_parameters->frame_period, fs,
      y_length, y);
  printf("WORLD: %d [msec]\n", timeGetTime() - elapsed_time);
}

void DestroyMemory(WorldParameters *world_parameters) {
  delete[] world_parameters->time_axis;
  delete[] world_parameters->f0;
  for (int i = 0; i < world_parameters->f0_length; ++i) {
    delete[] world_parameters->spectrogram[i];
    delete[] world_parameters->aperiodicity[i];
  }
  delete[] world_parameters->spectrogram;
  delete[] world_parameters->aperiodicity;
}

}  // namespace

//-----------------------------------------------------------------------------
// Test program.
// test.exe input.wav outout.wav f0 spec flag
// input.wav  : argv[1] Input file
// output.wav : argv[2] Output file
// f0         : argv[3] F0 scaling (a positive number)
// spec       : argv[4] Formant shift (a positive number)
//-----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
  if (argc != 2 && argc != 3 && argc != 4 && argc != 5) {
    printf("error\n");
    return -2;
  }

  // 2016/01/28: Important modification.
  // Memory allocation is carried out in advanse.
  // This is for compatibility with C language.
  int x_length = GetAudioLength(argv[1]);
  if (x_length <= 0) {
    if (x_length == 0)
      printf("error: File not found.\n");
    else
      printf("error: The file is not .wav format.\n");
    return -1;
  }
  double *x = new double[x_length];
  // wavread() must be called after GetAudioLength().
  int fs, nbit;
  wavread(argv[1], &fs, &nbit, x);
  DisplayInformation(fs, nbit, x_length);

  //---------------------------------------------------------------------------
  // Analysis part
  //---------------------------------------------------------------------------
  // 2016/02/02
  // A new struct is introduced to implement safe program.
  WorldParameters world_parameters = { 0 };
  // You must set fs and frame_period before analysis/synthesis.
  world_parameters.fs = fs;

  // 5.0 ms is the default value.
  // Generally, the inverse of the lowest F0 of speech is the best.
  // However, the more elapsed time is required.
  world_parameters.frame_period = 5.0;

  // F0 estimation
  F0Estimation(x, x_length, &world_parameters);

  // Spectral envelope estimation
  SpectralEnvelopeEstimation(x, x_length, &world_parameters);

  // Aperiodicity estimation by D4C
  AperiodicityEstimation(x, x_length, &world_parameters);

  // Note that F0 must not be changed until all parameters are estimated.
  ParameterModification(argc, argv, fs, world_parameters.f0_length,
    world_parameters.fft_size, world_parameters.f0,
    world_parameters.spectrogram);

  //---------------------------------------------------------------------------
  // Synthesis part
  //---------------------------------------------------------------------------
  // The length of the output waveform
  int y_length = static_cast<int>((world_parameters.f0_length - 1) *
    world_parameters.frame_period / 1000.0 * fs) + 1;
  double *y = new double[y_length];
  // Synthesis
  WaveformSynthesis(&world_parameters, fs, y_length, y);

  // Output
  wavwrite(y, y_length, fs, 16, argv[2]);

  delete[] y;
  delete[] x;
  DestroyMemory(&world_parameters);

  printf("complete.\n");
  return 0;
}
