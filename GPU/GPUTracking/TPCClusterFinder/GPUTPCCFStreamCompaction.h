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

/// \file StreamCompaction.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_STREAM_COMPACTION_H
#define O2_GPU_STREAM_COMPACTION_H

#include "clusterFinderDefs.h"
#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"
#include "GPUTPCClusterFinder.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class GPUTPCCFStreamCompaction : public GPUKernelTemplate
{

 public:
  enum K : int {
    scanStart = 0,
    scanUp = 1,
    scanTop = 2,
    scanDown = 3,
    compactDigits = 4,
  };

  struct GPUSharedMemory : public GPUKernelTemplate::GPUSharedMemoryScan64<int, GPUCA_THREAD_COUNT_SCAN> {
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

  template <int iKernel = GPUKernelTemplate::defaultKernel, typename... Args>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer, Args... args);

 private:
  static GPUd() void nativeScanUpStartImpl(int, int, int, int, GPUSharedMemory&,
                                           const uchar*, int*, int*,
                                           int);

  static GPUd() void nativeScanUpImpl(int, int, int, int, GPUSharedMemory&,
                                      int*, int*, int);

  static GPUd() void nativeScanTopImpl(int, int, int, int, GPUSharedMemory&,
                                       int*, int);

  static GPUd() void nativeScanDownImpl(int, int, int, int, GPUSharedMemory&,
                                        int*, const int*, unsigned int, int);

  static GPUd() void compactImpl(int, int, int, int, GPUSharedMemory&,
                                 const ChargePos*, ChargePos*,
                                 const uchar*, int*, const int*,
                                 int, tpccf::SizeT);
  static GPUd() int compactionElems(processorType& clusterer, int stage);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
