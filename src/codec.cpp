//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/05/09
//
// Coder/decoder functions for the spectral envelope and aperiodicity.
//-----------------------------------------------------------------------------
#include "world/codec.h"

#include <math.h>

#include "world/constantnumbers.h"
#include "world/fft.h"
#include "world/matlabfunctions.h"

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

//-----------------------------------------------------------------------------
// Frequency is converted into its mel representation.
//-----------------------------------------------------------------------------
static inline double FrequencyToMel(double frequency) {
  return world::kM0 * log(frequency / world::kF0 + 1.0);
}

//-----------------------------------------------------------------------------
// Mel is converted into frequency.
//-----------------------------------------------------------------------------
static inline double MelToFrequency(double mel) {
  return world::kF0 * (exp(mel / world::kM0) - 1.0);
}

//-----------------------------------------------------------------------------
// DCT for spectral envelope coding
//-----------------------------------------------------------------------------
static void DCTForCodec(const double *mel_spectrum, int max_dimension,
    const fft_complex *weight, const ForwardRealFFT *forward_real_fft,
    int number_of_dimensions, double *mel_cepstrum) {
  int bias = max_dimension / 2;
  for (int i = 0; i < max_dimension / 2; ++i) {
    forward_real_fft->waveform[i] = mel_spectrum[i * 2];
    forward_real_fft->waveform[i + bias] =
      mel_spectrum[max_dimension - (i * 2) - 1];
  }
  fft_execute(forward_real_fft->forward_fft);

  double normalization = sqrt(forward_real_fft->fft_size);
  for (int i = 0; i < number_of_dimensions; ++i)
    mel_cepstrum[i] = (forward_real_fft->spectrum[i][0] * weight[i][0] -
      forward_real_fft->spectrum[i][1] * weight[i][1]) / normalization;
}

//-----------------------------------------------------------------------------
// IDCT for spectral envelope decoding
//-----------------------------------------------------------------------------
static void IDCTForCodec(const double *mel_cepstrum, int max_dimension,
    const fft_complex *weight, const InverseComplexFFT *inverse_complex_fft,
    int number_of_dimensions, double *mel_spectrum) {
  double normalization = sqrt(inverse_complex_fft->fft_size);
  for (int i = 0; i < number_of_dimensions; ++i) {
    inverse_complex_fft->input[i][0] =
      mel_cepstrum[i] * weight[i][0] * normalization;
    inverse_complex_fft->input[i][1] =
      -mel_cepstrum[i] * weight[i][1] * normalization;
  }
  for (int i = number_of_dimensions; i < max_dimension; ++i) {
    inverse_complex_fft->input[i][0] = 0.0;
    inverse_complex_fft->input[i][1] = 0.0;
  }

  fft_execute(inverse_complex_fft->inverse_fft);

  for (int i = 0; i < max_dimension / 2; ++i) {
    mel_spectrum[i * 2] = inverse_complex_fft->output[i][0];
    mel_spectrum[(i * 2) + 1] =
      inverse_complex_fft->output[max_dimension - i - 1][0];
  }
}

//-----------------------------------------------------------------------------
// Spectral envelope in a frame is coded
//-----------------------------------------------------------------------------
static void CodeOneFrame(const double *log_spectral_envelope,
    const double *frequency_axis, int fft_size, const double *mel_axis,
    const fft_complex *weight, int max_dimension, int number_of_dimensions,
    const ForwardRealFFT *forward_real_fft, double *coded_spectral_envelope) {
  double *mel_spectrum = new double[max_dimension];
  interp1(frequency_axis, log_spectral_envelope, fft_size / 2 + 1,
      mel_axis, max_dimension, mel_spectrum);

  // DCT
  DCTForCodec(mel_spectrum, max_dimension, weight, forward_real_fft,
      number_of_dimensions, coded_spectral_envelope);

  delete[] mel_spectrum;
}

//-----------------------------------------------------------------------------
// Coded spectral envelope in a frame is decoded
//-----------------------------------------------------------------------------
static void DecodeOneFrame(const double *coded_spectral_envelope,
    const double *frequency_axis, int fft_size, const double *mel_axis,
    const fft_complex *weight, int max_dimension, int number_of_dimensions,
    const InverseComplexFFT *inverse_complex_fft, double *spectral_envelope) {
  double *mel_spectrum = new double[max_dimension + 2];

  // IDCT
  IDCTForCodec(coded_spectral_envelope, max_dimension, weight,
      inverse_complex_fft, number_of_dimensions, &mel_spectrum[1]);
  mel_spectrum[0] = mel_spectrum[1];
  mel_spectrum[max_dimension + 1] = mel_spectrum[max_dimension];

  interp1(mel_axis, mel_spectrum, max_dimension + 2, frequency_axis,
      fft_size / 2 + 1, spectral_envelope);

  for (int i = 0; i < fft_size / 2 + 1; ++i)
    spectral_envelope[i] = exp(spectral_envelope[i] / max_dimension);

  delete[] mel_spectrum;
}

//-----------------------------------------------------------------------------
// GetParameters() generates the required parameters.
//-----------------------------------------------------------------------------
static void GetParametersForCoding(double floor_frequency,
    double ceil_frequency, int fs, int fft_size, double *mel_axis,
    double *frequency_axis, fft_complex *weight) {
  int max_dimension = fft_size / 2;
  double floor_mel = FrequencyToMel(floor_frequency);
  double ceil_mel = FrequencyToMel(ceil_frequency);

  // Generate the mel axis and the weighting vector for DCT.
  for (int i = 0; i < max_dimension; ++i) {
    mel_axis[i] = (ceil_mel - floor_mel) * i / max_dimension + floor_mel;
    weight[i][0] = 2.0 * cos(i * world::kPi / fft_size) / sqrt(fft_size);
    weight[i][1] = 2.0 * sin(i * world::kPi / fft_size) / sqrt(fft_size);
  }
  weight[0][0] /= sqrt(2.0);

  // Generate the frequency axis on mel scale
  for (int i = 0; i < max_dimension; ++i)
    frequency_axis[i] = FrequencyToMel(static_cast<double>(i) * fs / fft_size);
}

//-----------------------------------------------------------------------------
// GetParameters() generates the required parameters.
//-----------------------------------------------------------------------------
static void GetParametersForDecoding(double floor_frequency,
    double ceil_frequency, int fs, int fft_size, int number_of_dimensions,
    double *mel_axis, double *frequency_axis, fft_complex *weight) {
  int max_dimension = fft_size / 2;
  double floor_mel = FrequencyToMel(floor_frequency);
  double ceil_mel = FrequencyToMel(ceil_frequency);

  // Generate the weighting vector for IDCT.
  for (int i = 0; i < number_of_dimensions; ++i) {
    weight[i][0] = cos(i * world::kPi / fft_size) * sqrt(fft_size);
    weight[i][1] = sin(i * world::kPi / fft_size) * sqrt(fft_size);
  }
  weight[0][0] /= sqrt(2.0);
  // Generate the mel axis for IDCT.
  for (int i = 0; i < max_dimension; ++i)
    mel_axis[i + 1] =
      MelToFrequency((ceil_mel - floor_mel) * i / max_dimension + floor_mel);
  mel_axis[0] = 0;
  mel_axis[max_dimension + 1] = fs / 2.0;

  // Generate the frequency axis
  for (int i = 0; i < fft_size / 2 + 1; ++i)
    frequency_axis[i] = static_cast<double>(i) * fs / fft_size;
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

  delete[] coarse_frequency_axis;
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

void CodeSpectralEnvelope(const double * const *spectrogram, int f0_length,
    int fs, int fft_size, int number_of_dimensions,
    double **coded_spectral_envelope) {
  double *mel_axis = new double[fft_size / 2];
  double *frequency_axis = new double[fft_size / 2 + 1];
  double *tmp_spectrum = new double[fft_size / 2 + 1];
  fft_complex *weight = new fft_complex[fft_size / 2];

  // Generation of the required parameters
  GetParametersForCoding(world::kFloorFrequency,
      MyMinDouble(fs / 2.0, world::kCeilFrequency), fs, fft_size,
      mel_axis, frequency_axis, weight);

  ForwardRealFFT forward_real_fft = { 0 };
  InitializeForwardRealFFT(fft_size / 2, &forward_real_fft);

  for (int i = 0; i < f0_length; ++i) {
    for (int j = 0; j < fft_size / 2 + 1; ++j)
      tmp_spectrum[j] = log(spectrogram[i][j]);
    CodeOneFrame(tmp_spectrum, frequency_axis, fft_size, mel_axis, weight,
        fft_size / 2, number_of_dimensions, &forward_real_fft,
        coded_spectral_envelope[i]);
  }

  DestroyForwardRealFFT(&forward_real_fft);
  delete[] weight;
  delete[] tmp_spectrum;
  delete[] frequency_axis;
  delete[] mel_axis;
}

void DecodeSpectralEnvelope(const double * const *coded_spectral_envelope,
    int f0_length, int fs, int fft_size, int number_of_dimensions,
    double **spectrogram) {
  double *mel_axis = new double[fft_size / 2 + 2];
  double *frequency_axis = new double[fft_size / 2 + 1];
  double *tmp_spectrum = new double[fft_size / 2 + 1];
  fft_complex *weight = new fft_complex[fft_size / 2];

  // Generation of the required parameters
  GetParametersForDecoding(world::kFloorFrequency,
      MyMinDouble(fs / 2.0, world::kCeilFrequency),
      fs, fft_size, number_of_dimensions, mel_axis, frequency_axis, weight);

  InverseComplexFFT inverse_complex_fft = { 0 };
  InitializeInverseComplexFFT(fft_size / 2, &inverse_complex_fft);

  for (int i = 0; i < f0_length; ++i) {
    DecodeOneFrame(coded_spectral_envelope[i], frequency_axis, fft_size,
        mel_axis, weight, fft_size / 2, number_of_dimensions,
        &inverse_complex_fft, spectrogram[i]);
  }

  DestroyInverseComplexFFT(&inverse_complex_fft);
  delete[] weight;
  delete[] tmp_spectrum;
  delete[] frequency_axis;
  delete[] mel_axis;
}
