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

/// \file GPUConstantMem.h
/// \author David Rohr

#ifndef GPUCONSTANTMEM_H
#define GPUCONSTANTMEM_H

#include "GPUTPCTracker.h"
#include "GPUParam.h"
#include "GPUDataTypes.h"

// Dummies for stuff not supported in legacy code (ROOT 5 / OPENCL1.2)
#if defined(GPUCA_NOCOMPAT_ALLCINT) && (!defined(GPUCA_GPULIBRARY) || !defined(GPUCA_ALIROOT_LIB))
#include "GPUTPCGMMerger.h"
#include "GPUTRDTracker.h"
#else
namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCGMMerger
{
};
class GPUTRDTracker
{
  void SetMaxData(const GPUTrackingInOutPointers& io) {}
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE
#endif

// Dummies for stuff not suppored in legacy code, or for what requires O2 headers while not available
#if defined(GPUCA_NOCOMPAT_ALLCINT) && (!defined(GPUCA_GPULIBRARY) || !defined(GPUCA_ALIROOT_LIB)) && defined(HAVE_O2HEADERS)
#include "GPUTPCConvert.h"
#include "GPUTPCCompression.h"
#include "GPUITSFitter.h"
#include "GPUTPCClusterFinder.h"
#else
#include "GPUO2FakeClasses.h"
#endif

#ifdef GPUCA_KERNEL_DEBUGGER_OUTPUT
#include "GPUKernelDebugOutput.h"
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
MEM_CLASS_PRE()
struct GPUConstantMem {
  MEM_CONSTANT(GPUParam)
  param;
  MEM_GLOBAL(GPUTPCTracker)
  tpcTrackers[GPUCA_NSLICES];
  GPUTPCConvert tpcConverter;
  GPUTPCCompression tpcCompressor;
  GPUTPCGMMerger tpcMerger;
  GPUTRDTracker trdTracker;
  GPUTPCClusterFinder tpcClusterer[GPUCA_NSLICES];
  GPUITSFitter itsFitter;
  GPUTrackingInOutPointers ioPtrs;
  GPUCalibObjectsConst calibObjects;
#ifdef GPUCA_KERNEL_DEBUGGER_OUTPUT
  GPUKernelDebugOutput debugOutput;
#endif
};

// Must be placed here, to avoid circular header dependency
GPUdi() GPUconstantref() const MEM_CONSTANT(GPUParam) & GPUProcessor::Param() const { return mConstantMem->param; }

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
