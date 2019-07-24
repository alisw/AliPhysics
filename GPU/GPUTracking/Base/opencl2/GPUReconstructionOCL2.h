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

/// \file GPUReconstructionOCL2.h
/// \author David Rohr

#ifndef GPURECONSTRUCTIONOCL2_H
#define GPURECONSTRUCTIONOCL2_H

#include "GPUReconstructionOCL.h"

#ifdef _WIN32
extern "C" __declspec(dllexport) GPUCA_NAMESPACE::gpu::GPUReconstruction* GPUReconstruction_Create_OCL2(const GPUCA_NAMESPACE::gpu::GPUSettingsProcessing& cfg);
#else
extern "C" GPUCA_NAMESPACE::gpu::GPUReconstruction* GPUReconstruction_Create_OCL2(const GPUCA_NAMESPACE::gpu::GPUSettingsProcessing& cfg);
#endif

namespace GPUCA_NAMESPACE::gpu
{
struct GPUReconstructionOCL2Internals;

class GPUReconstructionOCL2Backend : public GPUReconstructionOCL
{
 public:
  ~GPUReconstructionOCL2Backend() override = default;

 protected:
  GPUReconstructionOCL2Backend(const GPUSettingsProcessing& cfg);

  template <class T, int I = 0, typename... Args>
  int runKernelBackend(const krnlExec& x, const krnlRunRange& y, const krnlEvent& z, const Args&... args);
  template <class S, class T, int I = 0>
  S& getKernelObject(int num);

  int GetOCLPrograms() override;
  bool CheckPlatform(unsigned int i) override;
};

using GPUReconstructionOCL2 = GPUReconstructionKernels<GPUReconstructionOCL2Backend>;
} // namespace GPUCA_NAMESPACE::gpu

#endif
