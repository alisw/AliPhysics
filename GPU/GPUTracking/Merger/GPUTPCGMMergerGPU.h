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

/// \file GPUTPCGMMergerGPU.h
/// \author David Rohr

#ifndef GPUTPCGMMERGERGPUCA_H
#define GPUTPCGMMERGERGPUCA_H

#include "GPUGeneralKernels.h"
#include "GPUConstantMem.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCGMMergerGeneral : public GPUKernelTemplate
{
 public:
  GPUhdi() CONSTEXPR static GPUDataTypes::RecoStep GetRecoStep() { return GPUDataTypes::RecoStep::TPCMerging; }
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  typedef GPUTPCGMMerger processorType;
  GPUhdi() static processorType* Processor(GPUConstantMem& processors)
  {
    return &processors.tpcMerger;
  }
#endif
};

class GPUTPCGMMergerTrackFit : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, int mode);
#endif
};

class GPUTPCGMMergerFollowLoopers : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerSliceRefit : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, int iSlice);
#endif
};

class GPUTPCGMMergerUnpackGlobal : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, int iSlice);
#endif
};

class GPUTPCGMMergerUnpackSaveNumber : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, int id);
#endif
};

class GPUTPCGMMergerUnpackResetIds : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, int id);
#endif
};

class GPUTPCGMMergerResolve : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, char useOrigTrackParam, char mergeAll);
#endif
};

class GPUTPCGMMergerClearLinks : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, char nOutput);
#endif
};

class GPUTPCGMMergerMergeWithinPrepare : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerMergeSlicesPrepare : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, int border0, int border1, char useOrigTrackParam);
#endif
};

class GPUTPCGMMergerMergeBorders : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger, int iSlice, char withinSlice, char mergeMode);
#endif
};

class GPUTPCGMMergerMergeCE : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerLinkGlobalTracks : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerCollect : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerPrepareClusters : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerSortTracks : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerSortTracksQPt : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerSortTracksPrepare : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

class GPUTPCGMMergerFinalize : public GPUTPCGMMergerGeneral
{
 public:
#if !defined(GPUCA_ALIROOT_LIB) || !defined(GPUCA_GPUCODE)
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& merger);
#endif
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
