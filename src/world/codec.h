//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/03/22
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GetNumberOfAperiodicities provides the number of dimensions for aperiodicity
// coding. It is determined by only fs.
//
// Input:
//   fs       : Sampling frequency
//
// Output:
//   Number of aperiodicities
//-----------------------------------------------------------------------------
int GetNumberOfAperiodicities(int fs);

//-----------------------------------------------------------------------------
// CodeAperiodicity codes the aperiodicity. The number of dimensions is
// determined by fs.
//
// Input:
//   aperiodicity       : Aperiodicity before coding
//   f0_length          : Length of F0 contour
//   fs                 : Sampling frequency
//   fft_size           : FFT size
//
// Output:
//   coded_aperiodicity : Coded aperiodicity
//-----------------------------------------------------------------------------
void CodeAperiodicity(const double * const *aperiodicity, int f0_length,
  int fs, int fft_size, double **coded_aperiodicity);

//-----------------------------------------------------------------------------
// DecodeAperiodicity decodes the codedaperiodicity.
//
// Input:
//   coded_aperiodicity : Coded aperiodicity
//   f0_length          : Length of F0 contour
//   fs                 : Sampling frequency
//   fft_size           : FFT size
//
// Output:
//   aperiodicity       : Decoded aperiodicity
//-----------------------------------------------------------------------------
void DecodeAperiodicity(const double * const *coded_aperiodicity,
  int f0_length, int fs, int fft_size, double **aperiodicity);
