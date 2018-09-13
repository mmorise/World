//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//
// Matlab functions implemented for WORLD
// Since these functions are implemented as the same function of Matlab,
// the source code does not follow the style guide (Names of variables
// and functions).
// Please see the reference of Matlab to show the usage of functions.
// Caution:
//   The functions wavread() and wavwrite() were removed to the /src.
//   they were moved to the test/audioio.cpp. (2016/01/28)
//-----------------------------------------------------------------------------
#include "world/matlabfunctions.h"

#include <math.h>
#include <stdint.h>

#include "world/constantnumbers.h"

namespace {
//-----------------------------------------------------------------------------
// FilterForDecimate() calculates the coefficients of low-pass filter and
// carries out the filtering. This function is only used for decimate().
//-----------------------------------------------------------------------------
static void FilterForDecimate(const double *x, int x_length, int r, double *y) {
  double a[3], b[2];  // filter Coefficients
  switch (r) {
    case 11:  // fs : 44100 (default)
      a[0] = 2.450743295230728;
      a[1] = -2.06794904601978;
      a[2] = 0.59574774438332101;
      b[0] = 0.0026822508007163792;
      b[1] = 0.0080467524021491377;
      break;
    case 12:  // fs : 48000
      a[0] = 2.4981398605924205;
      a[1] = -2.1368928194784025;
      a[2] = 0.62187513816221485;
      b[0] = 0.0021097275904709001;
      b[1] = 0.0063291827714127002;
      break;
    case 10:
      a[0] = 2.3936475118069387;
      a[1] = -1.9873904075111861;
      a[2] = 0.5658879979027055;
      b[0] = 0.0034818622251927556;
      b[1] = 0.010445586675578267;
      break;
    case 9:
      a[0] = 2.3236003491759578;
      a[1] = -1.8921545617463598;
      a[2] = 0.53148928133729068;
      b[0] = 0.0046331164041389372;
      b[1] = 0.013899349212416812;
      break;
    case 8:  // fs : 32000
      a[0] = 2.2357462340187593;
      a[1] = -1.7780899984041358;
      a[2] = 0.49152555365968692;
      b[0] = 0.0063522763407111993;
      b[1] = 0.019056829022133598;
      break;
    case 7:
      a[0] = 2.1225239019534703;
      a[1] = -1.6395144861046302;
      a[2] = 0.44469707800587366;
      b[0] = 0.0090366882681608418;
      b[1] = 0.027110064804482525;
      break;
    case 6:  // fs : 24000 and 22050
      a[0] = 1.9715352749512141;
      a[1] = -1.4686795689225347;
      a[2] = 0.3893908434965701;
      b[0] = 0.013469181309343825;
      b[1] = 0.040407543928031475;
      break;
    case 5:
      a[0] = 1.7610939654280557;
      a[1] = -1.2554914843859768;
      a[2] = 0.3237186507788215;
      b[0] = 0.021334858522387423;
      b[1] = 0.06400457556716227;
      break;
    case 4:  // fs : 16000
      a[0] = 1.4499664446880227;
      a[1] = -0.98943497080950582;
      a[2] = 0.24578252340690215;
      b[0] = 0.036710750339322612;
      b[1] = 0.11013225101796784;
      break;
    case 3:
      a[0] = 0.95039378983237421;
      a[1] = -0.67429146741526791;
      a[2] = 0.15412211621346475;
      b[0] = 0.071221945171178636;
      b[1] = 0.21366583551353591;
      break;
    case 2:  // fs : 8000
      a[0] = 0.041156734567757189;
      a[1] = -0.42599112459189636;
      a[2] = 0.041037215479961225;
      b[0] = 0.16797464681802227;
      b[1] = 0.50392394045406674;
      break;
    default:
      a[0] = 0.0;
      a[1] = 0.0;
      a[2] = 0.0;
      b[0] = 0.0;
      b[1] = 0.0;
  }

  // Filtering on time domain.
  double w[3] = {0.0, 0.0, 0.0};
  double wt;
  for (int i = 0; i < x_length; ++i) {
    wt = x[i] + a[0] * w[0] + a[1] * w[1] + a[2] * w[2];
    y[i] = b[0] * wt + b[1] * w[0] + b[1] * w[1] + b[0] * w[2];
    w[2] = w[1];
    w[1] = w[0];
    w[0] = wt;
  }
}

}  // namespace

void fftshift(const double *x, int x_length, double *y) {
  for (int i = 0; i < x_length / 2; ++i) {
    y[i] = x[i + x_length / 2];
    y[i + x_length / 2] = x[i];
  }
}

void histc(const double *x, int x_length, const double *edges,
    int edges_length, int *index) {
  int count = 1;

  int i = 0;
  for (; i < edges_length; ++i) {
    index[i] = 1;
    if (edges[i] >= x[0]) break;
  }
  for (; i < edges_length; ++i) {
    if (edges[i] < x[count]) {
      index[i] = count;
    } else {
      index[i--] = count++;
    }
    if (count == x_length) break;
  }
  count--;
  for (i++; i < edges_length; ++i) index[i] = count;
}

void interp1(const double *x, const double *y, int x_length, const double *xi,
    int xi_length, double *yi) {
  double *h = new double[x_length - 1];
  double *p = new double[xi_length];
  double *s = new double[xi_length];
  int *k = new int[xi_length];

  for (int i = 0; i < x_length - 1; ++i) h[i] = x[i + 1] - x[i];
  for (int i = 0; i < xi_length; ++i) {
    p[i] = i;
    k[i] = 0;
  }

  histc(x, x_length, xi, xi_length, k);

  for (int i = 0; i < xi_length; ++i)
    s[i] = (xi[i] - x[k[i] - 1]) / h[k[i] - 1];

  for (int i = 0; i < xi_length; ++i)
    yi[i] = y[k[i] - 1] + s[i] * (y[k[i]] - y[k[i] - 1]);

  delete[] k;
  delete[] s;
  delete[] p;
  delete[] h;
}

void decimate(const double *x, int x_length, int r, double *y) {
  const int kNFact = 9;
  double *tmp1 = new double[x_length + kNFact * 2];
  double *tmp2 = new double[x_length + kNFact * 2];

  for (int i = 0; i < kNFact; ++i) tmp1[i] = 2 * x[0] - x[kNFact - i];
  for (int i = kNFact; i < kNFact + x_length; ++i) tmp1[i] = x[i - kNFact];
  for (int i = kNFact + x_length; i < 2 * kNFact + x_length; ++i)
    tmp1[i] = 2 * x[x_length - 1] - x[x_length - 2 - (i - (kNFact + x_length))];

  FilterForDecimate(tmp1, 2 * kNFact + x_length, r, tmp2);
  for (int i = 0; i < 2 * kNFact + x_length; ++i)
    tmp1[i] = tmp2[2 * kNFact + x_length - i - 1];
  FilterForDecimate(tmp1, 2 * kNFact + x_length, r, tmp2);
  for (int i = 0; i < 2 * kNFact + x_length; ++i)
    tmp1[i] = tmp2[2 * kNFact + x_length - i - 1];

  int nout = (x_length - 1) / r + 1;
  int nbeg = r - r * nout + x_length;

  int count = 0;
  for (int i = nbeg; i < x_length + kNFact; i += r)
    y[count++] = tmp1[i + kNFact - 1];

  delete[] tmp1;
  delete[] tmp2;
}

int matlab_round(double x) {
  return x > 0 ? static_cast<int>(x + 0.5) : static_cast<int>(x - 0.5);
}

void diff(const double *x, int x_length, double *y) {
  for (int i = 0; i < x_length - 1; ++i) y[i] = x[i + 1] - x[i];
}

void interp1Q(double x, double shift, const double *y, int x_length,
    const double *xi, int xi_length, double *yi) {
  double *xi_fraction = new double[xi_length];
  double *delta_y = new double[x_length];
  int *xi_base = new int[xi_length];

  double delta_x = shift;
  for (int i = 0; i < xi_length; ++i) {
    xi_base[i] = static_cast<int>((xi[i] - x) / delta_x);
    xi_fraction[i] = (xi[i] - x) / delta_x - xi_base[i];
  }
  diff(y, x_length, delta_y);
  delta_y[x_length - 1] = 0.0;

  for (int i = 0; i < xi_length; ++i)
    yi[i] = y[xi_base[i]] + delta_y[xi_base[i]] * xi_fraction[i];

  // Bug was fixed at 2013/07/14 by M. Morise
  delete[] xi_fraction;
  delete[] xi_base;
  delete[] delta_y;
}

// You must not use these variables.
// Note:
// I have no idea to implement the randn() and randn_reseed() without the
// global variables. If you have a good idea, please give me the information.
static uint32_t g_randn_x = 123456789;
static uint32_t g_randn_y = 362436069;
static uint32_t g_randn_z = 521288629;
static uint32_t g_randn_w = 88675123;

void randn_reseed() {
    g_randn_x = 123456789;
    g_randn_y = 362436069;
    g_randn_z = 521288629;
    g_randn_w = 88675123;
}

double randn(void) {
  uint32_t t;
  t = g_randn_x ^ (g_randn_x << 11);
  g_randn_x = g_randn_y;
  g_randn_y = g_randn_z;
  g_randn_z = g_randn_w;
  g_randn_w = (g_randn_w ^ (g_randn_w >> 19)) ^ (t ^ (t >> 8));

  uint32_t tmp = g_randn_w >> 4;
  for (int i = 0; i < 11; ++i) {
    t = g_randn_x ^ (g_randn_x << 11);
    g_randn_x = g_randn_y;
    g_randn_y = g_randn_z;
    g_randn_z = g_randn_w;
    g_randn_w = (g_randn_w ^ (g_randn_w >> 19)) ^ (t ^ (t >> 8));
    tmp += g_randn_w >> 4;
  }
  return tmp / 268435456.0 - 6.0;
}

void fast_fftfilt(const double *x, int x_length, const double *h, int h_length,
    int fft_size, const ForwardRealFFT *forward_real_fft,
    const InverseRealFFT *inverse_real_fft, double *y) {
  fft_complex *x_spectrum = new fft_complex[fft_size];

  for (int i = 0; i < x_length; ++i)
    forward_real_fft->waveform[i] = x[i] / fft_size;
  for (int i = x_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);
  for (int i = 0; i <= fft_size / 2; ++i) {
    x_spectrum[i][0] = forward_real_fft->spectrum[i][0];
    x_spectrum[i][1] = forward_real_fft->spectrum[i][1];
  }

  for (int i = 0; i < h_length; ++i)
    forward_real_fft->waveform[i] = h[i] / fft_size;
  for (int i = h_length; i < fft_size; ++i)
    forward_real_fft->waveform[i] = 0.0;
  fft_execute(forward_real_fft->forward_fft);

  for (int i = 0; i <= fft_size / 2; ++i) {
    inverse_real_fft->spectrum[i][0] =
      x_spectrum[i][0] * forward_real_fft->spectrum[i][0] -
      x_spectrum[i][1] * forward_real_fft->spectrum[i][1];
    inverse_real_fft->spectrum[i][1] =
      x_spectrum[i][0] * forward_real_fft->spectrum[i][1] +
      x_spectrum[i][1] * forward_real_fft->spectrum[i][0];
  }
  fft_execute(inverse_real_fft->inverse_fft);

  for (int i = 0; i < fft_size; ++i)
    y[i] = inverse_real_fft->waveform[i];

  delete[] x_spectrum;
}

double matlab_std(const double *x, int x_length) {
  double average = 0.0;
  for (int i = 0; i < x_length; ++i) average += x[i];
  average /= x_length;

  double s = 0.0;
  for (int i = 0; i < x_length; ++i) s += pow(x[i] - average, 2.0);
  s /= (x_length - 1);

  return sqrt(s);
}
