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

/// \file GPUGeneralKernels.h
/// \author David Rohr

#ifndef GPUGENERALKERNELS_H
#define GPUGENERALKERNELS_H

#include "GPUDef.h"
#include "GPUDataTypes.h"

#ifdef GPUCA_GPUCODE
#ifdef __CUDACC__
#include <cub/cub.cuh>
#define GPUCA_CUB cub
#elif defined(__HIPCC__)
#include <hipcub/hipcub.hpp>
#define GPUCA_CUB hipcub
#endif
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
MEM_CLASS_PRE()
struct GPUConstantMem;

class GPUKernelTemplate
{
 public:
  enum K { defaultKernel = 0 };

  MEM_CLASS_PRE()
  struct GPUSharedMemory {
  };

  template <class T, int I>
  struct GPUSharedMemoryScan64 {
    // Provides the shared memory resources for CUB collectives
#if (defined(__CUDACC__) || defined(__HIPCC__)) && defined(GPUCA_GPUCODE)
    typedef GPUCA_CUB::BlockScan<T, I> BlockScan;
    union {
      typename BlockScan::TempStorage cubTmpMem;
      int tmpBroadcast;
    };
#endif
  };

  typedef GPUconstantref() MEM_CONSTANT(GPUConstantMem) processorType;
  GPUhdi() CONSTEXPRRET static GPUDataTypes::RecoStep GetRecoStep() { return GPUCA_RECO_STEP::NoRecoStep; }
  MEM_TEMPLATE()
  GPUhdi() static processorType* Processor(MEM_TYPE(GPUConstantMem) & processors)
  {
    return &processors;
  }
#ifdef GPUCA_NOCOMPAT
  template <int iKernel, typename... Args>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & smem, processorType& processors, Args... args)
  {
  }
#else
  template <int iKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & smem, processorType& processors)
  {
  }
#endif
};

// Clean memory, ptr multiple of 16, size will be extended to multiple of 16
class GPUMemClean16 : public GPUKernelTemplate
{
 public:
  GPUhdi() CONSTEXPRRET static GPUDataTypes::RecoStep GetRecoStep() { return GPUCA_RECO_STEP::NoRecoStep; }
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & smem, processorType& processors, GPUglobalref() void* ptr, unsigned long size);
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
