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

#if !defined(GPUCA_OPENCL1) && !defined(GPUCA_ALIROOT_LIB)
// Files for TPC Merger
#include "GPUTPCGMMergerGPU.cxx"
#include "GPUTPCGMMerger.h"
#include "GPUTPCGMTrackParam.cxx"
#include "GPUTPCGMPhysicalTrackModel.cxx"
#include "GPUTPCGMPropagator.cxx"

#if defined(HAVE_O2HEADERS)
// Files for propagation with material
#include "MatLayerCylSet.cxx"
#include "MatLayerCyl.cxx"
#include "Ray.cxx"

// Files for GPU dEdx
#include "GPUdEdx.cxx"

// Files for TPC Transformation
#include "GPUTPCConvertKernel.cxx"

// Files for TPC Compression
#include "GPUTPCCompressionKernels.cxx"
#include "GPUTPCCompressionTrackModel.cxx"

// Files for TPC Cluster Finder
#include "ClusterAccumulator.cxx"
#include "GPUTPCCFStreamCompaction.cxx"
#include "GPUTPCCFChargeMapFiller.cxx"
#include "GPUTPCCFPeakFinder.cxx"
#include "GPUTPCCFNoiseSuppression.cxx"
#include "GPUTPCCFClusterizer.cxx"
#include "GPUTPCCFDeconvolution.cxx"
#include "GPUTPCCFMCLabelFlattener.cxx"
#include "GPUTPCCFDecodeZS.cxx"

// Files for TRD Tracking
#include "GPUTRDTrackerGPU.cxx"
#include "GPUTRDTrack.cxx"
#include "GPUTRDTracker.cxx"
#include "GPUTRDTrackletWord.cxx"
#include "TRDGeometryBase.cxx"

// Files for ITS Track Fit
#include "GPUITSFitterKernels.cxx"

#if !defined(GPUCA_O2_LIB) && defined(__HIPCC__) && !defined(GPUCA_NO_ITS_TRAITS)
#include "VertexerTraitsHIP.hip.cxx"
#include "ContextHIP.hip.cxx"
#include "DeviceStoreVertexerHIP.hip.cxx"
#include "ClusterLinesHIP.hip.cxx"
#include "UtilsHIP.hip.cxx"
#elif !defined(GPUCA_O2_LIB) && defined(__CUDACC__) && !defined(GPUCA_NO_ITS_TRAITS)
#include "TrackerTraitsNV.cu"
#include "VertexerTraitsGPU.cu"
#include "Context.cu"
#include "Stream.cu"
#include "DeviceStoreNV.cu"
#include "DeviceStoreVertexerGPU.cu"
#include "ClusterLinesGPU.cu"
#include "Utils.cu"
#endif // !defined(GPUCA_O2_LIB) && defined(__CUDACC__) && !defined(GPUCA_NO_ITS_TRAITS)

#endif // HAVE_O2HEADERS
#endif // (!defined(__OPENCL__) || defined(__OPENCLCPP__)) && !defined(GPUCA_ALIROOT_LIB)

#endif // GPURECONSTRUCTIONINCLUDESDEVICE_H
