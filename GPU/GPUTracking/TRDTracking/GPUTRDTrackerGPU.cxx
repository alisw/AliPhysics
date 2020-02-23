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

/// \file GPUTRDTrackerGPU.cxx
/// \author David Rohr

#include "GPUTRDTrackerGPU.h"
#include "GPUTRDGeometry.h"
#include "GPUConstantMem.h"
#if defined(WITH_OPENMP) && !defined(GPUCA_GPUCODE)
#include "GPUReconstruction.h"
#endif

using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUTRDTrackerGPU::Thread<0>(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& smem, processorType& processors)
{
#if defined(WITH_OPENMP) && !defined(GPUCA_GPUCODE)
#pragma omp parallel for num_threads(processors.trdTracker.GetRec().GetDeviceProcessingSettings().nThreads)
#endif
  for (int i = get_global_id(0); i < processors.trdTracker.NTracks(); i += get_global_size(0)) {
    processors.trdTracker.DoTrackingThread(i, get_global_id(0));
  }
}
