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

/// \file GPUReconstructionOCL1.h
/// \author David Rohr

#ifndef GPURECONSTRUCTIONOCL1_H
#define GPURECONSTRUCTIONOCL1_H

#include "GPUReconstructionOCL.h"

#ifdef _WIN32
extern "C" __declspec(dllexport) GPUCA_NAMESPACE::gpu::GPUReconstruction* GPUReconstruction_Create_OCL(const GPUCA_NAMESPACE::gpu::GPUSettingsProcessing& cfg);
#else
extern "C" GPUCA_NAMESPACE::gpu::GPUReconstruction* GPUReconstruction_Create_OCL(const GPUCA_NAMESPACE::gpu::GPUSettingsProcessing& cfg);
#endif

namespace GPUCA_NAMESPACE::gpu
{
struct GPUReconstructionOCL1Internals;

class GPUReconstructionOCL1Backend : public GPUReconstructionOCL
{
 public:
  ~GPUReconstructionOCL1Backend() override = default;

 protected:
  GPUReconstructionOCL1Backend(const GPUSettingsProcessing& cfg);

  template <class T, int I = 0, typename... Args>
  int runKernelBackend(krnlSetup& _xyz, const Args&... args);
  template <class S, class T, int I, bool MULTI>
  S& getKernelObject();

  RecoStepField AvailableRecoSteps() override { return (RecoStep::TPCSliceTracking); }
  bool ContextForAllPlatforms() override { return true; }
  bool CheckPlatform(unsigned int i) override;
  int GetOCLPrograms() override;
};

using GPUReconstructionOCL1 = GPUReconstructionKernels<GPUReconstructionOCL1Backend>;
} // namespace GPUCA_NAMESPACE::gpu

#endif
