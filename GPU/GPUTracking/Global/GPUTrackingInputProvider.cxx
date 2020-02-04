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

/// \file GPUTrackingInputProvider.cxx
/// \author David Rohr

#include "GPUTrackingInputProvider.h"
#include "GPUDataTypes.h"
#include "GPUReconstruction.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace o2::tpc;

void GPUTrackingInputProvider::InitializeProcessor() {}
void* GPUTrackingInputProvider::SetPointersInputZS(void* mem)
{
  if (holdsTPCZS && (mRec->GetRecoStepsGPU() & GPUDataTypes::RecoStep::TPCClusterFinding)) {
    computePointerWithAlignment(mem, mPzsMeta, 1);
    computePointerWithAlignment(mem, mPzsSizes, GPUTrackingInOutZS::NSLICES * GPUTrackingInOutZS::NENDPOINTS);
    computePointerWithAlignment(mem, mPzsPtrs, GPUTrackingInOutZS::NSLICES * GPUTrackingInOutZS::NENDPOINTS);
  }
  return mem;
}

void* GPUTrackingInputProvider::SetPointersInputGPUOnly(void* mem)
{
  return mem;
}

void GPUTrackingInputProvider::RegisterMemoryAllocation()
{
  mResourceZS = mRec->RegisterMemoryAllocation(this, &GPUTrackingInputProvider::SetPointersInputZS, GPUMemoryResource::MEMORY_INPUT, "InputZS");
  mRec->RegisterMemoryAllocation(this, &GPUTrackingInputProvider::SetPointersInputGPUOnly, GPUMemoryResource::MEMORY_INPUT | GPUMemoryResource::MEMORY_GPU, "InputGPU"); // TODO: move more here (example)
}

void GPUTrackingInputProvider::SetMaxData(const GPUTrackingInOutPointers& io)
{
  holdsTPCZS = io.tpcZS;
}
