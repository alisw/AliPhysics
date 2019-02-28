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

/// \file GPUTRDTrackerGPU.h
/// \author David Rohr

#ifndef GPUTRDTRACKERGPUCA_H
#define GPUTRDTRACKERGPUCA_H

#include "GPUGeneralKernels.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class GPUTRDTrackerGPU : public GPUKernelTemplate
{
 public:
  GPUhdi() static GPUDataTypes::RecoStep GetRecoStep() { return GPUCA_RECO_STEP::TRDTracking; }
#if defined(GPUCA_BUILD_TRD) || !defined(GPUCA_GPUCODE)
  template <int iKernel = 0>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUTPCSharedMemory& smem, processorType& processors);
#endif
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTRDTRACKERGPUCA_H
