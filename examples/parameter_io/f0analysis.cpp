//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/03/12
//
// Summary:
// This example estimates F0 from an audio file and saves it to a file.
//
// How to use:
// % f0analysis -h
//
// Related works: spanalysis.cpp, apanalysis.cpp, readandsynthesis.cpp
//-----------------------------------------------------------------------------
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "../../tools/audioio.h"
#include "../../tools/parameterio.h"
#include "world/harvest.h"
#include "world/constantnumbers.h"

namespace {

//-----------------------------------------------------------------------------
// Display how to use this program
//-----------------------------------------------------------------------------
void usage(char *argv) {
  printf("\n");
  printf(" %s - F0 estimation by Harvest\n", argv);
  printf("\n");
  printf("  usage:\n");
  printf("   %s input.wav [options]\n", argv);
  printf("  options:\n");
  printf("   -f f    : floor of frequency range (Hz) [40]\n");
  printf("   -c c    : ceil of frequency range (Hz)  [800]\n");
  printf("   -s s    : shift length (ms)             [5]\n");
  printf("   -o name : filename used for output      [output.f0]\n");
  printf("   -t      : text file is given            [binary]\n");
  printf("\n");
}

//-----------------------------------------------------------------------------
// Set parameters from command line options
//-----------------------------------------------------------------------------
int SetOption(int argc, char **argv, double *f0_floor, double *f0_ceil,
    double *frame_period, char *filename, int *text_flag) {
  while (--argc) {
    if (strcmp(argv[argc], "-f") == 0) *f0_floor = atof(argv[argc + 1]);
    if (strcmp(argv[argc], "-c") == 0) *f0_ceil = atof(argv[argc + 1]);
    if (strcmp(argv[argc], "-s") == 0) *frame_period = atof(argv[argc + 1]);
    if (strcmp(argv[argc], "-o") == 0)
      snprintf(filename, strlen(argv[argc + 1]) + 1, "%s", argv[argc + 1]);
    if (strcmp(argv[argc], "-t") == 0) *text_flag = 1;
    if (strcmp(argv[argc], "-h") == 0) {
      usage(argv[0]);
      return 0;
    }
  }
  return 1;
}

}  // namespace

//-----------------------------------------------------------------------------
// This example estimates F0 from an audio file and saves it to a file.
//-----------------------------------------------------------------------------
int main(int argc, char **argv) {
  // Command check
  if (argc < 2) return 0;

  // Default parameters
  HarvestOption option = { 0 };
  InitializeHarvestOption(&option);
  option.frame_period = 5.0;
  option.f0_floor = world::kFloorF0;
  option.f0_ceil = world::kCeilF0;
  char filename[200] = "output.f0";
  int text_flag = 0;

  // Options from command line
  if (SetOption(argc, argv, &option.f0_floor, &option.f0_ceil,
    &option.frame_period, filename, &text_flag) == 0) return -1;

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

  // F0 analysis
  int number_of_frames =
    GetSamplesForHarvest(fs, x_length, option.frame_period);
  double *f0 = new double[number_of_frames];
  double *temporal_positions = new double[number_of_frames];
  Harvest(x, x_length, fs, &option, temporal_positions, f0);

  // File output
  WriteF0(filename, number_of_frames, option.frame_period, temporal_positions,
      f0, text_flag);

  // Memory deallocation
  delete[] f0;
  delete[] temporal_positions;
  delete[] x;
  return 0;
}
