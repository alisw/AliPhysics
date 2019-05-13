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

/// \file GPUITSFitter.cxx
/// \author David Rohr, Maximiliano Puccio

#include "GPUITSFitter.h"

#include "ITStracking/Road.h"
#include "ITStracking/Cluster.h"
#include "GPUITSTrack.h"
#include "GPUReconstruction.h"

using namespace GPUCA_NAMESPACE::gpu;

#ifndef GPUCA_GPUCODE
void GPUITSFitter::InitializeProcessor()
{
}

void* GPUITSFitter::SetPointersInput(void* mem)
{
  computePointerWithAlignment(mem, mRoads, mNumberOfRoads);
  for (int i = 0; i < 7; i++) {
    computePointerWithAlignment(mem, mTF[i], mNTF[i]);
  }
  return mem;
}

void* GPUITSFitter::SetPointersTracks(void* mem)
{
  computePointerWithAlignment(mem, mTracks, mNMaxTracks);
  return mem;
}

void* GPUITSFitter::SetPointersMemory(void* mem)
{
  computePointerWithAlignment(mem, mMemory, 1);
  return mem;
}

void GPUITSFitter::RegisterMemoryAllocation()
{
  AllocateAndInitializeLate();
  mMemoryResInput = mRec->RegisterMemoryAllocation(this, &GPUITSFitter::SetPointersInput, GPUMemoryResource::MEMORY_INPUT, "ITSInput");
  mMemoryResTracks = mRec->RegisterMemoryAllocation(this, &GPUITSFitter::SetPointersTracks, GPUMemoryResource::MEMORY_OUTPUT, "ITSTracks");
  mMemoryResMemory = mRec->RegisterMemoryAllocation(this, &GPUITSFitter::SetPointersMemory, GPUMemoryResource::MEMORY_PERMANENT, "ITSMemory");
}

void GPUITSFitter::SetMaxData() { mNMaxTracks = mNumberOfRoads; }
#endif

void GPUITSFitter::clearMemory()
{
  new (mMemory) Memory;
}
