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

/// \file Deconvolution.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_DECONVOLUTION_H
#define O2_GPU_DECONVOLUTION_H

#include "clusterFinderDefs.h"
#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"
#include "GPUTPCClusterFinder.h"
#include "Array2D.h"
#include "PackedCharge.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class GPUTPCCFDeconvolution : public GPUKernelTemplate
{

 public:
  struct GPUSharedMemory : public GPUKernelTemplate::GPUSharedMemoryScan64<short, GPUCA_THREAD_COUNT_CLUSTERER> {
    ChargePos posBcast1[SCRATCH_PAD_WORK_GROUP_SIZE];
    uchar aboveThresholdBcast[SCRATCH_PAD_WORK_GROUP_SIZE];
    uchar buf[SCRATCH_PAD_WORK_GROUP_SIZE * SCRATCH_PAD_COUNT_N];
  };

  enum K : int {
    countPeaks,
  };

  static GPUd() void countPeaksImpl(int, int, int, int, GPUSharedMemory&, const Array2D<uchar>&, Array2D<PackedCharge>&, const ChargePos*, const uint);

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

  template <int iKernel = defaultKernel, typename... Args>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer, Args... args);

 private:
  static GPUd() char countPeaksAroundDigit(const ChargePos&, const Array2D<uchar>&);
  static GPUd() char countPeaksScratchpadInner(ushort, const uchar*, uchar*);
  static GPUd() char countPeaksScratchpadOuter(ushort, ushort, uchar, const uchar*);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
