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

/// \file GPUTPCCFMCLabelFlattener.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_GPUTPCCF_MCLABEL_FLATTENER_H
#define O2_GPU_GPUTPCCF_MCLABEL_FLATTENER_H

#include "GPUGeneralKernels.h"

#include "clusterFinderDefs.h"
#include "GPUTPCClusterFinder.h"
#include "GPUConstantMem.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

struct GPUTPCLinearLabels;

class GPUTPCCFMCLabelFlattener : public GPUKernelTemplate
{

 public:
  struct GPUSharedMemory {
  };

  enum K : int {
    setRowOffsets,
    flatten,
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

  template <int iKernel = defaultKernel, typename... Args>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer, Args... args);

  static void setGlobalOffsetsAndAllocate(GPUTPCClusterFinder&, GPUTPCLinearLabels&);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
