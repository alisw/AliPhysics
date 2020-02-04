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

/// \file GPUTPCClusterFinderKernels.h
/// \author David Rohr

#ifndef O2_GPU_GPUTPCCLUSTERFINDERKERNEL_H
#define O2_GPU_GPUTPCCLUSTERFINDERKERNEL_H

#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"
#include "GPUTPCSharedMemoryData.h"

#include "clusterFinderDefs.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCClusterFinderKernels : public GPUKernelTemplate
{
 public:
  class GPUTPCSharedMemory : public GPUKernelTemplate::GPUTPCSharedMemoryScan64<int, GPUCA_THREAD_COUNT_SCAN>
  {
   public:
    union {
      GPUTPCSharedMemoryData::search_t search;
      GPUTPCSharedMemoryData::noise_t noise;
      GPUTPCSharedMemoryData::count_t count;
      GPUTPCSharedMemoryData::build_t build;
      GPUTPCSharedMemoryData::zs_t zs;
    };
  };

  enum K : int {
    fillChargeMap = 0,
    resetMaps = 1,
    findPeaks = 2,
    noiseSuppression = 3,
    updatePeaks = 4,
    countPeaks = 5,
    computeClusters = 6,
    nativeScanUpStart = 7,
    nativeScanUp = 8,
    nativeScanTop = 9,
    nativeScanDown = 10,
    compactDigit = 11,
    decodeZS = 12
  };

#ifdef HAVE_O2HEADERS
  typedef GPUTPCClusterFinder processorType;
  GPUhdi() static processorType* Processor(GPUConstantMem& processors)
  {
    return processors.tpcClusterer;
  }
#endif

  GPUhdi() CONSTEXPR static GPUDataTypes::RecoStep GetRecoStep()
  {
    return GPUDataTypes::RecoStep::TPCClusterFinding;
  }
  template <int iKernel = 0, typename... Args>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUTPCSharedMemory& smem, processorType& clusterer, Args... args);

 private:
  GPUd() static int compactionElems(processorType& clusterer, int stage);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
