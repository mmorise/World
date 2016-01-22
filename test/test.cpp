// Copyright 2012-2015 Masanori Morise. All Rights Reserved.
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

// test.exe input.wav outout.wav f0 spec
// input.wav  : Input file
// output.wav : Output file
// f0         : F0 scaling (a positive number)
// spec       : Formant scaling (a positive number)

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

#include "./../src/d4c.h"
#include "./../src/dio.h"
#include "./../src/matlabfunctions.h"
#include "./../src/cheaptrick.h"
#include "./../src/stonemask.h"
#include "./../src/synthesis.h"

// Frame shift [msec]
#define FRAMEPERIOD 5.0

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

namespace {
bool CheckLoadedFile(double *x, int fs, int nbit, int x_length) {
  if (x == NULL) {
    printf("error: File not found.\n");
    return false;
  }

  printf("File information\n");
  printf("Sampling : %d Hz %d Bit\n", fs, nbit);
  printf("Length %d [sample]\n", x_length);
  printf("Length %f [sec]\n", static_cast<double>(x_length) / fs);
  return true;
}

void F0Estimation(double *x, int x_length, int fs, int f0_length, double *f0,
    double *time_axis) {
  double *refined_f0 = new double[f0_length];

  DioOption option;
  InitializeDioOption(&option);  // Initialize the option
  // Modification of the option
  option.frame_period = FRAMEPERIOD;
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

  printf("\nAnalysis\n");
  DWORD elapsed_time = timeGetTime();
  Dio(x, x_length, fs, option, time_axis, f0);
  printf("DIO: %d [msec]\n", timeGetTime() - elapsed_time);

  // StoneMask is carried out to improve the estimation performance.
  elapsed_time = timeGetTime();
  StoneMask(x, x_length, fs, time_axis, f0, f0_length, refined_f0);
  printf("StoneMask: %d [msec]\n", timeGetTime() - elapsed_time);

  for (int i = 0; i < f0_length; ++i) f0[i] = refined_f0[i];

  delete[] refined_f0;
  return;
}

void SpectralEnvelopeEstimation(double *x, int x_length, int fs,
  double *time_axis, double *f0, int f0_length, double **spectrogram) {
  CheapTrickOption option;
  InitializeCheapTrickOption(&option);  // Initialize the option
  option.q1 = -0.15; // This value may be better one for HMM speech synthesis.

  DWORD elapsed_time = timeGetTime();
  CheapTrick(x, x_length, fs, time_axis, f0, f0_length, &option, spectrogram);
  printf("CheapTrick: %d [msec]\n", timeGetTime() - elapsed_time);
}

void AperiodicityEstimation(double *x, int x_length, int fs, double *time_axis,
    double *f0, int f0_length, int fft_size, double **aperiodicity) {
  D4COption option;
  InitializeD4COption(&option);  // Initialize the option

  DWORD elapsed_time = timeGetTime();
  // option is not implemented in this version. This is for future update.
  // We can use "NULL" as the argument.
  D4C(x, x_length, fs, time_axis, f0, f0_length, fft_size, &option,
      aperiodicity);
  printf("D4C: %d [msec]\n", timeGetTime() - elapsed_time);
}

void ParameterModification(int argc, char *argv[], int fs, double *f0,
    int f0_length, double **spectrogram) {
  int fft_size = GetFFTSizeForCheapTrick(fs);
  // F0 scaling
  if (argc >= 4) {
    double shift = atof(argv[3]);
    for (int i = 0; i < f0_length; ++i) f0[i] *= shift;
  }
  // Spectral stretching
  if (argc >= 5) {
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
      if (ratio < 1.0) {
        for (int j = static_cast<int>(fft_size / 2.0 * ratio);
            j <= fft_size / 2; ++j)
          spectrogram[i][j] =
          spectrogram[i][static_cast<int>(fft_size / 2.0 * ratio) - 1];
      }
    }
    delete[] spectrum1;
    delete[] spectrum2;
    delete[] freq_axis1;
    delete[] freq_axis2;
  }
}

void WaveformSynthesis(double *f0, int f0_length, double **spectrogram,
    double **aperiodicity, int fft_size, double frame_period, int fs,
    int y_length, double *y) {
  DWORD elapsed_time;
  // Synthesis by the aperiodicity
  printf("\nSynthesis\n");
  elapsed_time = timeGetTime();
  Synthesis(f0, f0_length, spectrogram, aperiodicity,
      fft_size, FRAMEPERIOD, fs, y_length, y);
  printf("WORLD: %d [msec]\n", timeGetTime() - elapsed_time);
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
  int fs, nbit, x_length;
  double *x = wavread(argv[1], &fs, &nbit, &x_length);

  if (CheckLoadedFile(x, fs, nbit, x_length) == false) {
    printf("error: File not found.\n");
    return -1;
  }

  // Allocate memories
  // The number of samples for F0
  int f0_length = GetSamplesForDIO(fs, x_length, FRAMEPERIOD);
  double *f0 = new double[f0_length];
  double *time_axis = new double[f0_length];

  // FFT size for CheapTrick
  int fft_size = GetFFTSizeForCheapTrick(fs);
  double **spectrogram = new double *[f0_length];
  double **aperiodicity = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i) {
    spectrogram[i] = new double[fft_size / 2 + 1];
    aperiodicity[i] = new double[fft_size / 2 + 1];
  }

  // F0 estimation
  F0Estimation(x, x_length, fs, f0_length, f0, time_axis);

  // Spectral envelope estimation
  SpectralEnvelopeEstimation(x, x_length, fs, time_axis, f0, f0_length,
      spectrogram);

  // Aperiodicity estimation by D4C
  AperiodicityEstimation(x, x_length, fs, time_axis, f0, f0_length,
      fft_size, aperiodicity);

  // Note that F0 must not be changed until all parameters are estimated.
  ParameterModification(argc, argv, fs, f0, f0_length, spectrogram);

  // The length of the output waveform
  int y_length =
    static_cast<int>((f0_length - 1) * FRAMEPERIOD / 1000.0 * fs) + 1;
  double *y = new double[y_length];
  // Synthesis
  WaveformSynthesis(f0, f0_length, spectrogram, aperiodicity, fft_size,
      FRAMEPERIOD, fs, y_length, y);

  // Output
  wavwrite(y, y_length, fs, 16, argv[2]);

  printf("complete.\n");

  delete[] x;
  delete[] time_axis;
  delete[] f0;
  delete[] y;
  for (int i = 0; i < f0_length; ++i) {
    delete[] spectrogram[i];
    delete[] aperiodicity[i];
  }
  delete[] spectrogram;
  delete[] aperiodicity;

  return 0;
}
