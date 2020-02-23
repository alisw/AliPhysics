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

/// \file GPUITSFitterKernels.h
/// \author David Rohr, Maximiliano Puccio

#ifndef GPUITSFITTERKERNELS_H
#define GPUITSFITTERKERNELS_H

#include "GPUGeneralKernels.h"
namespace o2
{
namespace its
{
struct TrackingFrameInfo;
}
} // namespace o2

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCGMPropagator;
class GPUITSFitter;
class GPUITSTrack;

class GPUITSFitterKernel : public GPUKernelTemplate
{
 public:
  GPUhdi() CONSTEXPR static GPUDataTypes::RecoStep GetRecoStep() { return GPUDataTypes::RecoStep::ITSTracking; }
  template <int iKernel = defaultKernel>
  GPUd() static void Thread(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& processors);

 protected:
  GPUd() static bool fitTrack(GPUITSFitter& Fitter, GPUTPCGMPropagator& prop, GPUITSTrack& track, int start, int end, int step);
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
