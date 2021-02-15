//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] meiji.ac.jp (Masanori Morise)
// Last update: 2021/02/15
//
// This header file only defines constant numbers used for several function.
//-----------------------------------------------------------------------------
#ifndef WORLD_CONSTANT_NUMBERS_H_
#define WORLD_CONSTANT_NUMBERS_H_

namespace world {
  // for Dio()
  const double kCutOff = 50.0;

  // for StoneMask()
  const double kFloorF0StoneMask = 40.0;

  const double kPi = 3.1415926535897932384;
  const double kMySafeGuardMinimum = 0.000000000001;
  const double kEps = 0.00000000000000022204460492503131;
  const double kFloorF0 = 71.0;
  const double kCeilF0 = 800.0;
  const double kDefaultF0 = 500.0;
  const double kLog2 = 0.69314718055994529;
  // Maximum standard deviation not to be selected as a best f0.
  const double kMaximumValue = 100000.0;

  // Note to me (fs: 48000)
  // 71 Hz is the limit to maintain the FFT size at 2048.
  // If we use 70 Hz as FLOOR_F0, the FFT size of 4096 is required.

  // for D4C()
  const int kHanning = 1;
  const int kBlackman = 2;
  const double kFrequencyInterval = 3000.0;
  const double kUpperLimit = 15000.0;
  const double kThreshold = 0.85;
  const double kFloorF0D4C = 47.0;

  // for Codec (Mel scale)
  // S. Stevens & J. Volkmann,
  // The Relation of Pitch to Frequency: A Revised Scale,
  // American Journal of Psychology, vol. 53, no. 3, pp. 329-353, 1940.
  const double kM0 = 1127.01048;
  const double kF0 = 700.0;
  const double kFloorFrequency = 40.0;
  const double kCeilFrequency = 20000.0;

}  // namespace world

#endif  // WORLD_CONSTANT_NUMBERS_H_
