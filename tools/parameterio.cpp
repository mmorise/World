//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/03/12
//
// Save/Load functions for three speech parameters.
//-----------------------------------------------------------------------------
#include "./parameterio.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

namespace {
static void WriteOneParameter(FILE *fp, const char *text,
    double parameter, int size) {
  fwrite(text, 1, 4, fp);
  if (size == 4) {
    int parameter_int = static_cast<int>(parameter);
    fwrite(&parameter_int, 4, 1, fp);
  } else {
    fwrite(&parameter, 8, 1, fp);
  }
}

static void LoadParameters(FILE *fp, int *number_of_frames, int *fft_size,
    int *number_of_dimensions) {
  char data_check[12];
  fread(&data_check, 1, 4, fp);  // NOF
  fread(number_of_frames, 4, 1, fp);

  fread(&data_check, 1, 12, fp);  // FP (skipped)

  fread(&data_check, 1, 4, fp);  // FFT
  fread(fft_size, 4, 1, fp);

  fread(&data_check, 1, 4, fp);  // NOD
  fread(number_of_dimensions, 4, 1, fp);
  *number_of_dimensions =
    *number_of_dimensions == 0 ? *fft_size / 2 + 1 : *number_of_dimensions;

  fread(&data_check, 1, 8, fp);  // FS (skipped)
}

static int CheckHeader(FILE *fp, const char *text) {
  char data_check[5];
  fread(data_check, 1, 4, fp);  // "F0  "
  data_check[4] = '\0';
  if (0 != strcmp(data_check, text)) {
    printf("Header error.\n");
    fclose(fp);
    return 0;
  }
  return 1;
}

}  // namespace

void WriteF0(const char *filename, int f0_length, double frame_period,
    const double *temporal_positions, const double *f0, int text_flag) {
  if (text_flag == 1) {
    FILE *fp = fopen(filename, "w");
    if (NULL == fp) {
      printf("File cannot be opened.\n");
      return;
    }
    for (int i = 0; i < f0_length; ++i)
      fprintf(fp, "%.5f %.5f\r\n", temporal_positions[i], f0[i]);
    fclose(fp);
  } else {
    FILE *fp = fopen(filename, "wb");
    if (NULL == fp) {
      printf("File cannot be opened.\n");
      return;
    }

    // Header
    fwrite("F0  ", 1, 4, fp);

    // Parameters
    WriteOneParameter(fp, "NOF ", f0_length, 4);
    WriteOneParameter(fp, "FP  ", frame_period, 8);

    // Data
    fwrite(f0, 8, f0_length, fp);
    fclose(fp);
  }
}

int ReadF0(const char *filename, double *temporal_positions, double *f0) {
  FILE *fp = fopen(filename, "rb");
  if (NULL == fp) {
    printf("File cannot be opened.\n");
    return 0;
  }

  // Header
  if (CheckHeader(fp, "F0  ") == 0) return 0;

  // Parameters
  char data_check[5];
  fread(data_check, 1, 4, fp);  // "NOF "
  int number_of_frames;
  fread(&number_of_frames, 4, 1, fp);

  fread(data_check, 1, 4, fp);  // "FP  "
  double frame_period;
  fread(&frame_period, 8, 1, fp);

  // Data
  fread(f0, 8, number_of_frames, fp);

  fclose(fp);
  for (int i = 0; i < number_of_frames; ++i)
    temporal_positions[i] = i / 1000.0 * frame_period;
  return 1;
}

double GetHeaderInformation(const char *filename, const char *parameter) {
  FILE *fp = fopen(filename, "rb");
  if (NULL == fp) {
    printf("File cannot be opened.\n");
    return 0;
  }

  char data_check[5];
  data_check[4] = '\0';
  for (int i = 0; i < 13; ++i) {
    fread(data_check, 1, 4, fp);
    if (0 != strcmp(data_check, parameter)) continue;
    if (0 == strcmp(parameter, "FP  ")) {
      double answer;
      fread(&answer, 8, 1, fp);
      fclose(fp);
      return answer;
    } else {
      int answer;
      fread(&answer, 4, 1, fp);
      fclose(fp);
      return static_cast<double>(answer);
    }
  }
  return 0;
}

void WriteSpectralEnvelope(const char *filename, int fs, int f0_length,
    double frame_period, int fft_size, int number_of_dimensions,
    const double * const *spectrogram) {
  FILE *fp = fopen(filename, "wb");
  if (NULL == fp) {
    printf("File cannot be opened.\n");
    return;
  }

  // Header
  fwrite("SPEC", 1, 4, fp);

  // Parameters
  WriteOneParameter(fp, "NOF ", f0_length, 4);
  WriteOneParameter(fp, "FP  ", frame_period, 8);
  WriteOneParameter(fp, "FFT ", fft_size, 4);
  WriteOneParameter(fp, "NOD ", number_of_dimensions, 4);
  WriteOneParameter(fp, "FS  ", fs, 4);

  number_of_dimensions =
    number_of_dimensions == 0 ? fft_size / 2 + 1 : number_of_dimensions;

  // Data
  for (int i = 0; i < f0_length; ++i)
    fwrite(spectrogram[i], 8, number_of_dimensions, fp);
  fclose(fp);
}

int ReadSpectralEnvelope(const char *filename, double **spectrogram) {
  FILE *fp = fopen(filename, "rb");
  if (NULL == fp) {
    printf("File cannot be opened.\n");
    return 0;
  }

  // Header
  if (CheckHeader(fp, "SPEC") == 0) return 0;

  // Parameters
  int number_of_frames, fft_size, number_of_dimensions;
  LoadParameters(fp, &number_of_frames, &fft_size, &number_of_dimensions);

  // Data
  for (int i = 0; i < number_of_frames; ++i)
    fread(spectrogram[i], 8, number_of_dimensions, fp);

  fclose(fp);
  return 1;
}

void WriteAperiodicity(const char *filename, int fs, int f0_length,
    double frame_period, int fft_size, int number_of_dimensions,
    const double * const *aperiodicity) {
  FILE *fp = fopen(filename, "wb");
  if (NULL == fp) {
    printf("File cannot be opened.\n");
    return;
  }

  // Header
  fwrite("AP  ", 1, 4, fp);

  // Parameters
  WriteOneParameter(fp, "NOF ", f0_length, 4);
  WriteOneParameter(fp, "FP  ", frame_period, 8);
  WriteOneParameter(fp, "FFT ", fft_size, 4);
  WriteOneParameter(fp, "NOD ", number_of_dimensions, 4);
  WriteOneParameter(fp, "FS  ", fs, 4);
  number_of_dimensions =
    number_of_dimensions == 0 ? fft_size / 2 + 1 : number_of_dimensions;

  // Data
  for (int i = 0; i < f0_length; ++i)
    fwrite(aperiodicity[i], 8, number_of_dimensions, fp);
  fclose(fp);
}

int ReadAperiodicity(const char *filename, double **aperiodicity) {
  FILE *fp = fopen(filename, "rb");
  if (NULL == fp) {
    printf("File cannot be opened.\n");
    return 0;
  }

  // Header
  if (CheckHeader(fp, "AP  ") == 0) return 0;

  // Parameters
  int number_of_frames, fft_size, number_of_dimensions;
  LoadParameters(fp, &number_of_frames, &fft_size, &number_of_dimensions);

  // Data
  for (int i = 0; i < number_of_frames; ++i)
    fread(aperiodicity[i], 8, number_of_dimensions, fp);

  fclose(fp);
  return 1;
}
