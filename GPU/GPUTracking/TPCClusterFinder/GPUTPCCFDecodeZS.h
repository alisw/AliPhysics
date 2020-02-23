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

/// \file DecodeZS.h
/// \author David Rohr

#ifndef O2_GPU_DECODE_ZS_H
#define O2_GPU_DECODE_ZS_H

#include "clusterFinderDefs.h"
#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"
#include "DataFormatsTPC/ZeroSuppression.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class GPUTPCClusterFinder;

class GPUTPCCFDecodeZS : public GPUKernelTemplate
{
 public:
  struct GPUSharedMemory {
    CA_SHARED_STORAGE(unsigned int ZSPage[o2::tpc::TPCZSHDR::TPC_ZS_PAGE_SIZE / sizeof(unsigned int)]);
    unsigned int RowClusterOffset[o2::tpc::TPCZSHDR::TPC_MAX_ZS_ROW_IN_ENDPOINT];
    unsigned int nRowsRegion;
    unsigned int regionStartRow;
    unsigned int nThreadsPerRow;
    unsigned int rowStride;
    unsigned int decodeBits;
    float decodeBitsFactor;
  };

  enum K : int {
    decodeZS,
  };

  static GPUd() void decode(GPUTPCClusterFinder& clusterer, GPUSharedMemory& s, int nBlocks, int nThreads, int iBlock, int iThread);

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
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
