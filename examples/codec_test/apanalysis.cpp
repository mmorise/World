//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/05/01
//
// Summary:
// This example estimates the aperiodicity from an audio file
// and then saves the result to a file.
//
// How to use:
// % apanalysis -h
//
// Related works: f0analysis.cpp, spanalysis.cpp, readandsynthesis.cpp
//-----------------------------------------------------------------------------
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "../../tools/audioio.h"
#include "../../tools/parameterio.h"
#include "world/cheaptrick.h"  // used for determining the default FFT size.
#include "world/codec.h"
#include "world/d4c.h"

namespace {

//-----------------------------------------------------------------------------
// Display how to use this program
//-----------------------------------------------------------------------------
void usage(char *argv) {
  printf("\n");
  printf(" %s - aperiodicity estimation by D4C\n", argv);
  printf("\n");
  printf("  usage:\n");
  printf("   %s input.wav input.f0 [options]\n", argv);
  printf("  options:\n");
  printf("   -f f    : FFT size (samples)            [variable]\n");
  printf("           : Default depends on fs (44100 -> 2048, 16000 -> 1024)\n");
  printf("   -t t    : threshhold used in D4C Lovetrain  [0.85]\n");
  printf("   -c      : compression                       [false]\n");
  printf("   -o name : filename used for output          [output.ap]\n");
  printf("\n");
}

//-----------------------------------------------------------------------------
// Set parameters from command line options
//-----------------------------------------------------------------------------
int SetOption(int argc, char **argv, int *fft_size, double *threshold,
    int *compression_flag, char *filename) {
  while (--argc) {
    if (strcmp(argv[argc], "-f") == 0) *fft_size = atoi(argv[argc + 1]);
    if (strcmp(argv[argc], "-t") == 0) *threshold = atof(argv[argc + 1]);
    if (strcmp(argv[argc], "-o") == 0)
      snprintf(filename, 200, argv[argc + 1]);
    if (strcmp(argv[argc], "-c") == 0) *compression_flag = 1;
    if (strcmp(argv[argc], "-h") == 0) {
      usage(argv[0]);
      return 0;
    }
  }
  return 1;
}

//-----------------------------------------------------------------------------
// Write coded aperiodicity
//-----------------------------------------------------------------------------
void WriteCodedAperiodicity(const char *filename,
    const double * const *aperiodicity, int fs, int f0_length,
    double frame_period, int fft_size) {
  int number_of_aperiodicities = GetNumberOfAperiodicities(fs);
  double **coded_aperiodicity = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i)
    coded_aperiodicity[i] = new double[number_of_aperiodicities];

  CodeAperiodicity(aperiodicity, f0_length, fs, fft_size,
    coded_aperiodicity);
  WriteAperiodicity(filename, fs, f0_length, frame_period,
    fft_size, number_of_aperiodicities, coded_aperiodicity);

  for (int i = 0; i < f0_length; ++i) delete[] coded_aperiodicity[i];
  delete[] coded_aperiodicity;
}

}  // namespace

//-----------------------------------------------------------------------------
// This example estimates the aperiodicity from an audio file
// and then saves the result to a file.
//-----------------------------------------------------------------------------
int main(int argc, char **argv) {
  // Command check
  if (argc < 2) return 0;
  if (0 == strcmp(argv[1], "-h")) {
    usage(argv[0]);
    return -1;
  }

  // Read F0 information
  int f0_length = static_cast<int>(GetHeaderInformation(argv[2], "NOF "));
  double frame_period = GetHeaderInformation(argv[2], "FP  ");
  double *f0 = new double[f0_length];
  double *temporal_positions = new double[f0_length];
  ReadF0(argv[2], temporal_positions, f0);

  // Read an audio file
  int x_length = GetAudioLength(argv[1]);
  if (x_length <= 0) {
    if (x_length == 0) {
      printf("error: File not found.\n");
    } else {
      printf("error: File is not .wav format.\n");
    }
    return -1;
  }
  double *x = new double[x_length];
  int fs, nbit;
  wavread(argv[1], &fs, &nbit, x);

  // Default parameters
  D4COption option = { 0 };
  InitializeD4COption(&option);
  option.threshold = 0.85;
  char filename[200] = "output.ap";
  CheapTrickOption c_option = { 0 };
  InitializeCheapTrickOption(fs, &c_option);
  int fft_size = c_option.fft_size;
  int compression_flag = 0;

  // Options from command line
  if (SetOption(argc, argv, &fft_size, &option.threshold, &compression_flag,
    filename) == 0) return 0;

  // Aperiodicity analysis
  double **aperiodicity = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i)
    aperiodicity[i] = new double[fft_size / 2 + 1];
  D4C(x, x_length, fs, temporal_positions, f0, f0_length, fft_size,
      &option, aperiodicity);

  // File output
  if (compression_flag == 0) {
    // If you want to write the raw aperiodicity.
    WriteAperiodicity(filename, fs, f0_length, frame_period,
        fft_size, 0, aperiodicity);
  } else {
    // If you want to write the coded aperiodicity.
    WriteCodedAperiodicity(filename, aperiodicity, fs, f0_length, frame_period,
        fft_size);
  }

  // Memory deallocation
  for (int i = 0; i < f0_length; ++i) delete[] aperiodicity[i];
  delete[] aperiodicity;
  delete[] f0;
  delete[] temporal_positions;
  delete[] x;
  return 0;
}
