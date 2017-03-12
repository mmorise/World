//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/03/12
//
// .wav input/output functions were modified for compatibility with C language.
// Since these functions (wavread() and wavwrite()) are roughly implemented,
// we recommend more suitable functions provided by other organizations.
// This file is independent of WORLD project and for the test.cpp.
//-----------------------------------------------------------------------------
#include "./audioio.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#if (defined (__WIN32__) || defined (_WIN32)) && !defined (__MINGW32__)
#pragma warning(disable : 4996)
#endif

namespace {

static inline int MyMaxInt(int x, int y) {
  return x > y ? x : y;
}

static inline int MyMinInt(int x, int y) {
  return x < y ? x : y;
}

//-----------------------------------------------------------------------------
// CheckHeader() checks the .wav header. This function can only support the
// monaural wave file. This function is only used in waveread().
//-----------------------------------------------------------------------------
static int CheckHeader(FILE *fp) {
  char data_check[5];
  fread(data_check, 1, 4, fp);  // "RIFF"
  data_check[4] = '\0';
  if (0 != strcmp(data_check, "RIFF")) {
    printf("RIFF error.\n");
    return 0;
  }
  fseek(fp, 4, SEEK_CUR);
  fread(data_check, 1, 4, fp);  // "WAVE"
  if (0 != strcmp(data_check, "WAVE")) {
    printf("WAVE error.\n");
    return 0;
  }
  fread(data_check, 1, 4, fp);  // "fmt "
  if (0 != strcmp(data_check, "fmt ")) {
    printf("fmt error.\n");
    return 0;
  }
  fread(data_check, 1, 4, fp);  // 1 0 0 0
  if (!(16 == data_check[0] && 0 == data_check[1] &&
      0 == data_check[2] && 0 == data_check[3])) {
    printf("fmt (2) error.\n");
    return 0;
  }
  fread(data_check, 1, 2, fp);  // 1 0
  if (!(1 == data_check[0] && 0 == data_check[1])) {
    printf("Format ID error.\n");
    return 0;
  }
  fread(data_check, 1, 2, fp);  // 1 0
  if (!(1 == data_check[0] && 0 == data_check[1])) {
    printf("This function cannot support stereo file\n");
    return 0;
  }
  return 1;
}

//-----------------------------------------------------------------------------
// GetParameters() extracts fp, nbit, wav_length from the .wav file
// This function is only used in wavread().
//-----------------------------------------------------------------------------
static int GetParameters(FILE *fp, int *fs, int *nbit, int *wav_length) {
  char data_check[5] = {0};
  data_check[4] = '\0';
  unsigned char for_int_number[4];
  fread(for_int_number, 1, 4, fp);
  *fs = 0;
  for (int i = 3; i >= 0; --i) *fs = *fs * 256 + for_int_number[i];
  // Quantization
  fseek(fp, 6, SEEK_CUR);
  fread(for_int_number, 1, 2, fp);
  *nbit = for_int_number[0];

  // Skip until "data" is found. 2011/03/28
  while (0 != fread(data_check, 1, 1, fp)) {
    if (data_check[0] == 'd') {
      fread(&data_check[1], 1, 3, fp);
      if (0 != strcmp(data_check, "data"))
        fseek(fp, -3, SEEK_CUR);
      else
        break;
    }
  }
  if (0 != strcmp(data_check, "data")) {
    printf("data error.\n");
    return 0;
  }

  fread(for_int_number, 1, 4, fp);  // "data"
  *wav_length = 0;
  for (int i = 3; i >= 0; --i)
    *wav_length = *wav_length * 256 + for_int_number[i];
  *wav_length /= (*nbit / 8);
  return 1;
}

}  // namespace

void wavwrite(const double *x, int x_length, int fs, int nbit,
    const char *filename) {
  FILE *fp = fopen(filename, "wb");
  if (NULL == fp) {
    printf("File cannot be opened.\n");
    return;
  }

  char text[4] = {'R', 'I', 'F', 'F'};
  uint32_t long_number = 36 + x_length * 2;
  fwrite(text, 1, 4, fp);
  fwrite(&long_number, 4, 1, fp);

  text[0] = 'W';
  text[1] = 'A';
  text[2] = 'V';
  text[3] = 'E';
  fwrite(text, 1, 4, fp);
  text[0] = 'f';
  text[1] = 'm';
  text[2] = 't';
  text[3] = ' ';
  fwrite(text, 1, 4, fp);

  long_number = 16;
  fwrite(&long_number, 4, 1, fp);
  int16_t short_number = 1;
  fwrite(&short_number, 2, 1, fp);
  short_number = 1;
  fwrite(&short_number, 2, 1, fp);
  long_number = fs;
  fwrite(&long_number, 4, 1, fp);
  long_number = fs * 2;
  fwrite(&long_number, 4, 1, fp);
  short_number = 2;
  fwrite(&short_number, 2, 1, fp);
  short_number = 16;
  fwrite(&short_number, 2, 1, fp);

  text[0] = 'd';
  text[1] = 'a';
  text[2] = 't';
  text[3] = 'a';
  fwrite(text, 1, 4, fp);
  long_number = x_length * 2;
  fwrite(&long_number, 4, 1, fp);

  int16_t tmp_signal;
  for (int i = 0; i < x_length; ++i) {
    tmp_signal = static_cast<int16_t>(MyMaxInt(-32768,
        MyMinInt(32767, static_cast<int>(x[i] * 32767))));
    fwrite(&tmp_signal, 2, 1, fp);
  }

  fclose(fp);
}

int GetAudioLength(const char *filename) {
  FILE *fp = fopen(filename, "rb");
  if (NULL == fp) {
    return 0;
  }

  if (0 == CheckHeader(fp)) {
    fclose(fp);
    return -1;
  }

  char data_check[5] = { 0 };
  data_check[4] = '\0';
  unsigned char for_int_number[4];

  // Quantization
  fseek(fp, 10, SEEK_CUR);
  fread(for_int_number, 1, 2, fp);
  int nbit = for_int_number[0];

  while (0 != fread(data_check, 1, 1, fp)) {
    if ('d' == data_check[0]) {
      fread(&data_check[1], 1, 3, fp);
      if (0 != strcmp(data_check, "data"))
        fseek(fp, -3, SEEK_CUR);
      else
        break;
    }
  }
  if (0 != strcmp(data_check, "data")) {
    fclose(fp);
    return -1;
  }

  fread(for_int_number, 1, 4, fp);  // "data"
  fclose(fp);

  int wav_length = 0;
  for (int i = 3; i >= 0; --i)
    wav_length = wav_length * 256 + for_int_number[i];
  wav_length /= (nbit / 8);

  return wav_length;
}

void wavread(const char* filename, int *fs, int *nbit, double *x) {
  FILE *fp = fopen(filename, "rb");
  if (NULL == fp) {
    printf("File not found.\n");
    return;
  }

  if (0 == CheckHeader(fp)) {
    fclose(fp);
    return;
  }

  int x_length;
  if (0 == GetParameters(fp, fs, nbit, &x_length)) {
    fclose(fp);
    return;
  }

  int quantization_byte = *nbit / 8;
  double zero_line = pow(2.0, *nbit - 1);
  double tmp, sign_bias;
  unsigned char for_int_number[4];
  for (int i = 0; i < x_length; ++i) {
    sign_bias = tmp = 0.0;
    fread(for_int_number, 1, quantization_byte, fp);  // "data"
    if (for_int_number[quantization_byte-1] >= 128) {
      sign_bias = pow(2.0, *nbit - 1);
      for_int_number[quantization_byte - 1] =
        for_int_number[quantization_byte - 1] & 0x7F;
    }
    for (int j = quantization_byte - 1; j >= 0; --j)
      tmp = tmp * 256.0 + for_int_number[j];
    x[i] = (tmp - sign_bias) / zero_line;
  }
  fclose(fp);
}
