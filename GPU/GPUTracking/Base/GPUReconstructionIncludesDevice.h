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

/// \file GPUReconstructionIncludesDevice.h
/// \author David Rohr

#ifndef GPURECONSTRUCTIONINCLUDESDEVICE_H
#define GPURECONSTRUCTIONINCLUDESDEVICE_H

#include "GPUDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
}
} // namespace GPUCA_NAMESPACE
using namespace GPUCA_NAMESPACE::gpu;

#include "GPUTPCTrackParam.cxx"
#include "GPUTPCTrack.cxx"
#include "GPUTPCHitArea.cxx"
#include "GPUTPCGrid.cxx"
#include "GPUTPCRow.cxx"
#include "GPUParam.cxx"
#include "GPUTPCTracker.cxx"

#include "GPUGeneralKernels.cxx"

#include "GPUTPCTrackletSelector.cxx"
#include "GPUTPCNeighboursFinder.cxx"
#include "GPUTPCNeighboursCleaner.cxx"
#include "GPUTPCStartHitsFinder.cxx"
#include "GPUTPCStartHitsSorter.cxx"
#include "GPUTPCTrackletConstructor.cxx"

#ifdef GPUCA_BUILD_MERGER
#include "GPUTPCGMMergerGPU.cxx"
#include "GPUTPCGMMerger.h"
#include "GPUTPCGMTrackParam.cxx"
#include "GPUTPCGMPhysicalTrackModel.cxx"
#include "GPUTPCGMPropagator.cxx"
#ifdef HAVE_O2HEADERS
#include "MatLayerCylSet.cxx"
#include "MatLayerCyl.cxx"
#include "Ray.cxx"
#endif
#endif

#ifdef GPUCA_BUILD_DEDX
#include "GPUdEdx.cxx"
#endif

#ifdef GPUCA_BUILD_TPCCONVERT
#include "GPUTPCConvertKernel.cxx"
#endif

#ifdef GPUCA_BUILD_TPCCOMPRESSION
#include "GPUTPCCompressionKernels.cxx"
#include "GPUTPCCompressionTrackModel.cxx"
#endif

#ifdef GPUCA_BUILD_TRD
#include "GPUTRDTrackerGPU.cxx"
#include "GPUTRDTrack.cxx"
#include "GPUTRDTracker.cxx"
#include "GPUTRDTrackletWord.cxx"
#include "TRDGeometryBase.cxx"
#endif

#ifdef GPUCA_BUILD_ITS
#include "GPUITSFitterKernels.cxx"
#if !defined(GPUCA_O2_LIB) && defined(__CUDACC__)
#include "TrackerTraitsNV.cu"
#include "Context.cu"
#include "Stream.cu"
#include "DeviceStoreNV.cu"
#include "Utils.cu"
#endif
#endif

#endif
