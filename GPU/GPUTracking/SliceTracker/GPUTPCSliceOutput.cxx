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

/// \file GPUTPCSliceOutput.cxx
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#include "GPUOutputControl.h"
#include "GPUTPCSliceOutput.h"
#include "GPUCommonMath.h"

using namespace GPUCA_NAMESPACE::gpu;

unsigned int GPUTPCSliceOutput::EstimateSize(unsigned int nOfTracks, unsigned int nOfTrackClusters)
{
  // calculate the amount of memory [bytes] needed for the event
  return sizeof(GPUTPCSliceOutput) + sizeof(GPUTPCTrack) * nOfTracks + sizeof(GPUTPCSliceOutCluster) * nOfTrackClusters;
}

#ifndef GPUCA_GPUCODE
void GPUTPCSliceOutput::Allocate(GPUTPCSliceOutput*& ptrOutput, int nTracks, int nTrackHits, GPUOutputControl* outputControl, void*& internalMemory)
{
  // Allocate All memory needed for slice output
  const size_t memsize = EstimateSize(nTracks, nTrackHits);

  if (outputControl && outputControl->OutputType != GPUOutputControl::AllocateInternal) {
    if (outputControl->OutputMaxSize - ((char*)outputControl->OutputPtr - (char*)outputControl->OutputBase) < memsize) {
      outputControl->EndOfSpace = 1;
      ptrOutput = nullptr;
      return;
    }
    ptrOutput = reinterpret_cast<GPUTPCSliceOutput*>(outputControl->OutputPtr);
    outputControl->OutputPtr = (char*)outputControl->OutputPtr + memsize;
  } else {
    if (internalMemory) {
      free(internalMemory);
    }
    internalMemory = malloc(memsize);
    ptrOutput = reinterpret_cast<GPUTPCSliceOutput*>(internalMemory);
  }
  ptrOutput->SetMemorySize(memsize);
}
#endif
