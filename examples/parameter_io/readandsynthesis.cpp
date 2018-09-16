//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/03/12
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
      snprintf(filename, strlen(argv[argc + 1]) + 1, "%s", argv[argc + 1]);
    if (strcmp(argv[argc], "-h") == 0) {
      usage(argv[0]);
      return 0;
    }
  }
  return 1;
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

  // Read three parameters
  ReadF0(argv[1], temporal_positions, f0);
  ReadSpectralEnvelope(argv[2], spectrogram);
  ReadAperiodicity(argv[3], aperiodicity);

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
