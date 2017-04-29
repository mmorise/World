//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/04/29
//
// Summary:
// This example reads three parameters and generates a waveform from them.
//
// How to use:
// % readandsynthesis -h
//
// Related works: f0analysis.cpp, spanalysis.cpp, apanalysis.cpp
//-----------------------------------------------------------------------------
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "../../tools/audioio.h"
#include "../../tools/parameterio.h"
#include "world/codec.h"
#include "world/synthesis.h"

namespace {
void usage(char *argv) {
  printf("\n");
  printf(" %s - synthesis from three parameters\n", argv);
  printf("\n");
  printf("  usage:\n");
  printf("   %s input.f0 input.sp input.ap [option]\n", argv);
  printf("  option:\n");
  printf("   -o name : filename used for output          [output.wav]\n");
  printf("\n");
}

int SetOption(int argc, char **argv, char *filename) {
  while (--argc) {
    if (strcmp(argv[argc], "-o") == 0)
      snprintf(filename, 200, argv[argc + 1]);
    if (strcmp(argv[argc], "-h") == 0) {
      usage(argv[0]);
      return 0;
    }
  }
  return 1;
}

void ReadCodedAperiodicity(const char *filename, int f0_length, int fs,
    int fft_size, int number_of_dimensions, double **aperiodicity) {
  double **coded_aperiodicity = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i)
    coded_aperiodicity[i] = new double[number_of_dimensions];
  ReadAperiodicity(filename, coded_aperiodicity);
  DecodeAperiodicity(coded_aperiodicity, f0_length, fs, fft_size,
      aperiodicity);
  for (int i = 0; i < f0_length; ++i) delete[] coded_aperiodicity[i];
  delete[] coded_aperiodicity;
}

void ReadCodedSpectralEnvelope(const char *filename, int f0_length, int fs,
    int fft_size, int number_of_dimensions, double **spectrogram) {
  double **coded_spectral_envelope = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i)
    coded_spectral_envelope[i] = new double[number_of_dimensions];
  ReadSpectralEnvelope(filename, coded_spectral_envelope);
  DecodeSpectralEnvelope(coded_spectral_envelope, f0_length, fs,
      fft_size, number_of_dimensions, spectrogram);
  for (int i = 0; i < f0_length; ++i) delete[] coded_spectral_envelope[i];
  delete[] coded_spectral_envelope;
}

}  // namespace

//-----------------------------------------------------------------------------
// This example reads three parameters and generates a waveform.
//-----------------------------------------------------------------------------
int main(int argc, char **argv) {
  // Command check
  if (argc < 2) return 0;
  if (0 == strcmp(argv[1], "-h")) {
    usage(argv[0]);
    return -1;
  }

  // Get parameters required for memory allocation
  int f0_length = static_cast<int>(GetHeaderInformation(argv[1], "NOF "));
  int fft_size = static_cast<int>(GetHeaderInformation(argv[2], "FFT "));
  int fs = static_cast<int>(GetHeaderInformation(argv[2], "FS  "));
  double frame_period = GetHeaderInformation(argv[2], "FP  ");
  char filename[200] = "output.wav";

  // Option from command line
  if (SetOption(argc, argv, filename) == 0) return 0;

  // Memory allocation
  double *f0 = new double[f0_length];
  double *temporal_positions = new double[f0_length];
  double **spectrogram = new double *[f0_length];
  double **aperiodicity = new double *[f0_length];
  for (int i = 0; i < f0_length; ++i) {
    spectrogram[i] = new double[fft_size / 2 + 1];
    aperiodicity[i] = new double[fft_size / 2 + 1];
  }

  // Read F0
  ReadF0(argv[1], temporal_positions, f0);

  // Read spectral envelope
  int number_of_dimensions =
    static_cast<int>(GetHeaderInformation(argv[2], "NOD "));
  if (number_of_dimensions == 0) {
    // If the spectral envelope is not coded.
    ReadSpectralEnvelope(argv[2], spectrogram);
  } else {
    // If the spectral envelope is coded.
    ReadCodedSpectralEnvelope(argv[2], f0_length, fs, fft_size,
      number_of_dimensions, spectrogram);
  }

  // Read aperiodicity
  number_of_dimensions =
    static_cast<int>(GetHeaderInformation(argv[3], "NOD "));
  if (number_of_dimensions == 0) {
    // If the aperiodicity is not coded.
    ReadAperiodicity(argv[3], aperiodicity);
  } else {
    // If the aperiodicity is coded.
    ReadCodedAperiodicity(argv[3], f0_length, fs, fft_size,
        number_of_dimensions, aperiodicity);
  }

  // Synthesis
  int y_length = static_cast<int>(f0_length * frame_period / 1000.0 * fs);
  double *y = new double[y_length];
  Synthesis(f0, f0_length, spectrogram, aperiodicity, fft_size, frame_period,
      fs, y_length, y);

  // File output
  wavwrite(y, y_length, fs, 16, filename);

  // Memory deallocation
  delete[] y;
  for (int i = 0; i < f0_length; ++i) {
    delete[] spectrogram[i];
    delete[] aperiodicity[i];
  }
  delete[] spectrogram;
  delete[] aperiodicity;
  delete[] f0;
  delete[] temporal_positions;
  return 0;
}
