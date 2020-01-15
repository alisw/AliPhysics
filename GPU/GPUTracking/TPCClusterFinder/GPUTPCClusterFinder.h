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

/// \file GPUTPCClusterFinder.h
/// \author David Rohr

#ifndef O2_GPU_GPUTPCCLUSTERFINDER_H
#define O2_GPU_GPUTPCCLUSTERFINDER_H

#include "GPUDef.h"
#include "GPUProcessor.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

namespace deprecated
{
class PackedDigit;
class ClusterNative;
} // namespace deprecated

class GPUTPCClusterFinder : public GPUProcessor
{
 public:
  struct Memory {
    size_t nDigits = 0;
    size_t nPeaks = 0;
    size_t nClusters = 0;
  };

#ifndef GPUCA_GPUCODE
  void InitializeProcessor();
  void RegisterMemoryAllocation();
  void SetMaxData(const GPUTrackingInOutPointers& io);

  void* SetPointersInput(void* mem);
  void* SetPointersOutput(void* mem);
  void* SetPointersScratch(void* mem);
  void* SetPointersMemory(void* mem);

  size_t getNSteps(size_t items) const;
  void SetNMaxDigits(size_t n);
#endif

  deprecated::PackedDigit* mPdigits = nullptr;
  deprecated::PackedDigit* mPpeaks = nullptr;
  deprecated::PackedDigit* mPfilteredPeaks = nullptr;
  unsigned char* mPisPeak = nullptr;
  ushort* mPchargeMap = nullptr;
  unsigned char* mPpeakMap = nullptr;
  uint* mPclusterInRow = nullptr;
  deprecated::ClusterNative* mPclusterByRow = nullptr;
  int* mPbuf = nullptr;
  Memory* mPmemory = nullptr;

  int mISlice = 0;
  constexpr static int mScanWorkGroupSize = GPUCA_THREAD_COUNT_SCAN;
  size_t mNMaxClusterPerRow = 0;
  size_t mNMaxDigits = 0;
  size_t mNMaxPeaks = 0;
  size_t mNMaxClusters = 0;
  size_t mBufSize = 0;
  size_t mNBufs = 0;

  unsigned short mMemoryId = 0;

#ifndef GPUCA_GPUCODE
  void DumpDigits(std::ostream& out);
  void DumpChargeMap(std::ostream& out);
  void DumpPeaks(std::ostream& out);
  void DumpPeaksCompacted(std::ostream& out);
  void DumpSuppressedPeaks(std::ostream& out);
  void DumpSuppressedPeaksCompacted(std::ostream& out);
  void DumpCountedPeaks(std::ostream& out);
  void DumpClusters(std::ostream& out);
#endif
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
