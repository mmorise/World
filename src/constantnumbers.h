//-----------------------------------------------------------------------------
// Copyright 2012-2016 Masanori Morise. All Rights Reserved.
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
//
// This header file only defines constant numbers used for several function.
//-----------------------------------------------------------------------------
#ifndef WORLD_CONSTANT_NUMBERS_H_
#define WORLD_CONSTANT_NUMBERS_H_

namespace world {
  const double kPi = 3.1415926535897932384;
  const double kMySafeGuardMinimum = 0.000000000001;
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
}  // namespace world

#endif
