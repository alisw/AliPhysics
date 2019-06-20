//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file GPUReconstructionTimeframe.h
/// \author David Rohr

#ifndef GPURECONSTRUCTIONTIMEFRAME_H
#define GPURECONSTRUCTIONTIMEFRAME_H

#include "GPUChainTracking.h"
#include "GPUDataTypes.h"
#include <vector>
#include <random>
#include <tuple>

namespace o2
{
namespace tpc
{
struct ClusterNative;
}
} // namespace o2

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct ClusterNativeAccessFullTPC;

class GPUReconstructionTimeframe
{
 public:
  GPUReconstructionTimeframe(GPUChainTracking* rec, int (*read)(int), int nEvents);
  int LoadCreateTimeFrame(int iEvent);
  int LoadMergedEvents(int iEvent);
  int ReadEventShifted(int i, float shiftZ, float minZ = -1e6, float maxZ = -1e6, bool silent = false);
  void MergeShiftedEvents();

  static constexpr int ORBIT_RATE = 11245;
  static constexpr int DRIFT_TIME = 93000;
  static constexpr int TPCZ = 250;
  static constexpr int TIME_ORBIT = 1000000000 / ORBIT_RATE;

 private:
  constexpr static unsigned int NSLICES = GPUReconstruction::NSLICES;

  void SetDisplayInformation(int iCol);

  GPUChainTracking* mChain;
  int (*mReadEvent)(int);
  int mNEventsInDirectory;

  std::uniform_real_distribution<double> mDisUniReal;
  std::uniform_int_distribution<unsigned long long int> mDisUniInt;
  std::mt19937_64 mRndGen1;
  std::mt19937_64 mRndGen2;

  int mTrainDist = 0;
  float mCollisionProbability = 0.;
  int mMaxBunchesFull;
  int mMaxBunches;

  int mNTotalCollisions = 0;

  long long int mEventStride;
  int mSimBunchNoRepeatEvent;
  std::vector<char> mEventUsed;
  std::vector<std::tuple<GPUTrackingInOutPointers, GPUChainTracking::InOutMemory, o2::tpc::ClusterNativeAccessFullTPC>> mShiftedEvents;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
