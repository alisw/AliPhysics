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

/// \file GPUTPCGlobalTracking.h
/// \author David Rohr

#ifndef GPUTPCGLOBALTRACKING_H
#define GPUTPCGLOBALTRACKING_H

#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
MEM_CLASS_PRE()
class GPUTPCTracker;

#if !defined(__OPENCL__) || defined(__OPENCLCPP__)
class GPUTPCGlobalTracking : public GPUKernelTemplate
{
 public:
  struct GPUSharedMemory {
    CA_SHARED_STORAGE(MEM_LG(GPUTPCRow) mRows[GPUCA_ROW_COUNT]);
  };

  typedef GPUconstantref() MEM_GLOBAL(GPUTPCTracker) processorType;
  GPUhdi() CONSTEXPRRET static GPUDataTypes::RecoStep GetRecoStep() { return GPUCA_RECO_STEP::TPCSliceTracking; }
  GPUhdi() static processorType* Processor(MEM_TYPE(GPUConstantMem) & processors)
  {
    return processors.tpcTrackers;
  }
  template <int iKernel = GPUKernelTemplate::defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & smem, processorType& tracker);

  GPUd() static int GlobalTrackingSliceOrder(int iSlice);

 private:
  GPUd() static int PerformGlobalTrackingRun(GPUTPCTracker& tracker, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() smem, const GPUTPCTracker& sliceSource, int iTrack, int rowIndex, float angle, int direction);
  GPUd() static void PerformGlobalTracking(int nBlocks, int nThreads, int iBlock, int iThread, const GPUTPCTracker& tracker, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() smem, GPUTPCTracker& sliceTarget, bool right);
};
#endif

class GPUTPCGlobalTrackingCopyNumbers : public GPUKernelTemplate
{
 public:
  typedef GPUconstantref() MEM_GLOBAL(GPUTPCTracker) processorType;
  GPUhdi() CONSTEXPRRET static GPUDataTypes::RecoStep GetRecoStep() { return GPUCA_RECO_STEP::TPCSliceTracking; }
  MEM_TEMPLATE()
  GPUhdi() static processorType* Processor(MEM_TYPE(GPUConstantMem) & processors)
  {
    return processors.tpcTrackers;
  }
  template <int iKernel = GPUKernelTemplate::defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & smem, processorType& tracker, int n);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCTRACKLETCONSTRUCTOR_H
