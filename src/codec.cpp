//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/03/22
//
// coder/decoder functions for the spectral envelope and aperiodicity.
// Note: Codec for spectral envelope is not implemented yet.
//-----------------------------------------------------------------------------
#include "world/codec.h"

#include <math.h>

#include "world/matlabfunctions.h"
#include "world/constantnumbers.h"

namespace {
//-----------------------------------------------------------------------------
// Aperiodicity is initialized by the value 1.0 - world::kMySafeGuardMinimum.
// This value means the frame/frequency index is aperiodic.
//-----------------------------------------------------------------------------
static void InitializeAperiodicity(int f0_length, int fft_size,
    double **aperiodicity) {
  for (int i = 0; i < f0_length; ++i)
    for (int j = 0; j < fft_size / 2 + 1; ++j)
      aperiodicity[i][j] = 1.0 - world::kMySafeGuardMinimum;
}

//-----------------------------------------------------------------------------
// This function identifies whether this frame is voiced or unvoiced.
//-----------------------------------------------------------------------------
static int CheckVUV(const double *coarse_aperiodicity,
    int number_of_aperiodicities, double *tmp_aperiodicity) {
  double tmp = 0.0;
  for (int i = 0; i < number_of_aperiodicities; ++i) {
    tmp += coarse_aperiodicity[i];
    tmp_aperiodicity[i + 1] = coarse_aperiodicity[i];
  }
  tmp /= number_of_aperiodicities;

  return tmp > -0.5 ? 1 : 0;  // -0.5 is not optimized, but okay.
}

//-----------------------------------------------------------------------------
// Aperiodicity is obtained from the coded aperiodicity.
//-----------------------------------------------------------------------------
static void GetAperiodicity(const double *coarse_frequency_axis,
    const double *coarse_aperiodicity, int number_of_aperiodicities,
    const double *frequency_axis, int fft_size, double *aperiodicity) {
  interp1(coarse_frequency_axis, coarse_aperiodicity,
      number_of_aperiodicities + 2, frequency_axis, fft_size / 2 + 1,
      aperiodicity);
  for (int i = 0; i <= fft_size / 2; ++i)
    aperiodicity[i] = pow(10.0, aperiodicity[i] / 20.0);
}

}  // namespace

int GetNumberOfAperiodicities(int fs) {
  return  static_cast<int>(MyMinDouble(world::kUpperLimit, fs / 2.0 -
    world::kFrequencyInterval) / world::kFrequencyInterval);
}

void CodeAperiodicity(const double * const *aperiodicity, int f0_length,
    int fs, int fft_size, double **coded_aperiodicity) {
  int number_of_aperiodicities = GetNumberOfAperiodicities(fs);
  double *coarse_frequency_axis = new double[number_of_aperiodicities];
  for (int i = 0; i < number_of_aperiodicities; ++i)
    coarse_frequency_axis[i] = world::kFrequencyInterval * (i + 1.0);

  double *log_aperiodicity = new double[fft_size / 2 + 1];

  for (int i = 0; i < f0_length; ++i) {
    for (int j = 0; j < fft_size / 2 + 1; ++j)
      log_aperiodicity[j] = 20 * log10(aperiodicity[i][j]);
    interp1Q(0, static_cast<double>(fs) / fft_size, log_aperiodicity,
        fft_size / 2 + 1, coarse_frequency_axis, number_of_aperiodicities,
        coded_aperiodicity[i]);
  }

  delete[] log_aperiodicity;
}

void DecodeAperiodicity(const double * const *coded_aperiodicity,
    int f0_length, int fs, int fft_size, double **aperiodicity) {
  InitializeAperiodicity(f0_length, fft_size, aperiodicity);
  int number_of_aperiodicities = GetNumberOfAperiodicities(fs);
  double *frequency_axis = new double[fft_size / 2 + 1];
  for (int i = 0; i <= fft_size / 2; ++i)
    frequency_axis[i] = static_cast<double>(fs) / fft_size * i;

  double *coarse_frequency_axis = new double[number_of_aperiodicities + 2];
  for (int i = 0; i <= number_of_aperiodicities; ++i)
    coarse_frequency_axis[i] = i * world::kFrequencyInterval;
  coarse_frequency_axis[number_of_aperiodicities + 1] = fs / 2.0;

  double *coarse_aperiodicity = new double[number_of_aperiodicities + 2];
  coarse_aperiodicity[0] = -60.0;
  coarse_aperiodicity[number_of_aperiodicities + 1] =
    -world::kMySafeGuardMinimum;

  for (int i = 0; i < f0_length; ++i) {
    if (CheckVUV(coded_aperiodicity[i], number_of_aperiodicities,
      coarse_aperiodicity) == 1) continue;
    GetAperiodicity(coarse_frequency_axis, coarse_aperiodicity,
        number_of_aperiodicities, frequency_axis, fft_size, aperiodicity[i]);
  }

  delete[] coarse_aperiodicity;
  delete[] coarse_frequency_axis;
  delete[] frequency_axis;
}
