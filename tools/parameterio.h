//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/03/12
//-----------------------------------------------------------------------------
#ifndef WORLD_PARAMETERIO_H_
#define WORLD_PARAMETERIO_H_

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
// WriteF0() writes the F0 contour.
//
// Input:
//   filename           : Filename used for file output
//   f0_length          : Length of F0 contour
//   frame_period       : Frame shift used for analysis
//   temporal_positions : Time axis
//   f0                 : F0 contour
//   text_flag          : The file is written as text (NOT binary).
//-----------------------------------------------------------------------------
void WriteF0(const char *filename, int f0_length, double frame_period,
  const double *temporal_positions, const double *f0, int text_flag);

//-----------------------------------------------------------------------------
// ReadF0() reads the F0 contour from a file..
//
// Input:
//   filename           : Filename used for file output
//
// Output:
//   temporal_positions : Time axis
//   f0                 : F0 contour
//
// Note: Memory must be allocated before calling this function.
//-----------------------------------------------------------------------------
int ReadF0(const char *filename, double *temporal_positions, double *f0);

//-----------------------------------------------------------------------------
// GetHeaderInformation() reads a parameter from a file.
//
// Input:
//   filename             : Filename used for file output
//   parameter            : Required parameter
//                        : "NOF ": number of samples (int)
//                        : "FP  ": frame shift (double)
//                        : "FFT ": FFT size (int)
//                        : "NOD ": number of dimensions (int)
//                        : "FS  ": sampling frequency
//
// Note: F0 does not include "FFT ", "NOD ", and "FS  ".
//       These are ignored in cases where they are used for F0 file.
//-----------------------------------------------------------------------------
double GetHeaderInformation(const char *filename, const char *parameter);

//-----------------------------------------------------------------------------
// WriteSpectralEnvelope() writes the spectral envelope.
//
// Input:
//   filename             : Filename used for file output
//   fs                   : Sampling frequency
//   f0_length            : Length of F0 contour
//   frame_period         : Frame shift used for analysis
//   fft_size             : FFT size used for analysis
//   number_of_dimensions : Number of dimensions per frame
//   spectrogram          : Spectral envelope estimated by CheapTrick
//-----------------------------------------------------------------------------
void WriteSpectralEnvelope(const char *filename, int fs, int f0_length,
  double frame_period, int fft_size, int number_of_dimensions,
  const double * const *spectrogram);

//-----------------------------------------------------------------------------
// ReadSpectralEnvelope() reads the spectral envelope from a file.
//
// Input:
//   filename           : Filename used for file output
//
// Output:
//   spectrogram        : Spectral envelope estimated by CheapTrick
//
// Note: Memory must be allocated before calling this function.
//-----------------------------------------------------------------------------
int ReadSpectralEnvelope(const char *filename, double **spectrogram);

//-----------------------------------------------------------------------------
// WriteAperiodicity() writes the aperiodicity.
//
// Input:
//   filename             : Filename used for file output
//   fs                   : Sampling frequency
//   f0_length            : Length of F0 contour
//   frame_period         : Frame shift used for analysis
//   fft_size             : FFT size used for analysis
//   number_of_dimensions : Number of dimensions per frame
//   spectrogram          : Spectral envelope estimated by CheapTrick
//-----------------------------------------------------------------------------
void WriteAperiodicity(const char *filename, int fs, int f0_length,
  double frame_period, int fft_size, int number_of_dimensions,
  const double * const *aperiodicity);

//-----------------------------------------------------------------------------
// ReadAperiodicity() reads the aperiodicity from a file.
//
// Input:
//   filename           : Filename used for file output
//
// Output:
//   spectrogram        : Spectral envelope estimated by CheapTrick
//
// Note: Memory must be allocated before calling this function.
//-----------------------------------------------------------------------------
int ReadAperiodicity(const char *filename, double **aperiodicity);

#ifdef __cplusplus
}
#endif

#endif  // WORLD_PARAMETERIO_H_
