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
  if (mHoldTPCZS) {
    computePointerWithAlignment(mem, mPzsMeta);
    computePointerWithAlignment(mem, mPzsSizes, GPUTrackingInOutZS::NSLICES * GPUTrackingInOutZS::NENDPOINTS);
    computePointerWithAlignment(mem, mPzsPtrs, GPUTrackingInOutZS::NSLICES * GPUTrackingInOutZS::NENDPOINTS);
  }
  return mem;
}

void* GPUTrackingInputProvider::SetPointersInputClusterNativeAccess(void* mem)
{
  if (mHoldTPCClusterNative) {
    computePointerWithAlignment(mem, mPclusterNativeAccess);
  }
  return mem;
}

void* GPUTrackingInputProvider::SetPointersInputClusterNativeBuffer(void* mem)
{
  if (mHoldTPCClusterNative) {
    computePointerWithAlignment(mem, mPclusterNativeBuffer, mNClusterNative);
  }
  return mem;
}

void* GPUTrackingInputProvider::SetPointersInputClusterNativeOutput(void* mem)
{
  if (mHoldTPCClusterNativeOutput) {
    computePointerWithAlignment(mem, mPclusterNativeOutput, mNClusterNative);
  }
  return mem;
}

void GPUTrackingInputProvider::RegisterMemoryAllocation()
{
  mResourceZS = mRec->RegisterMemoryAllocation(this, &GPUTrackingInputProvider::SetPointersInputZS, GPUMemoryResource::MEMORY_INPUT, "InputZS");
  mResourceClusterNativeAccess = mRec->RegisterMemoryAllocation(this, &GPUTrackingInputProvider::SetPointersInputClusterNativeAccess, GPUMemoryResource::MEMORY_INPUT, "ClusterNativeAccess");
  mResourceClusterNativeBuffer = mRec->RegisterMemoryAllocation(this, &GPUTrackingInputProvider::SetPointersInputClusterNativeBuffer, GPUMemoryResource::MEMORY_INPUT_FLAG | GPUMemoryResource::MEMORY_GPU | GPUMemoryResource::MEMORY_EXTERNAL | GPUMemoryResource::MEMORY_CUSTOM, "ClusterNativeBuffer");
  mResourceClusterNativeOutput = mRec->RegisterMemoryAllocation(this, &GPUTrackingInputProvider::SetPointersInputClusterNativeOutput, GPUMemoryResource::MEMORY_OUTPUT_FLAG | GPUMemoryResource::MEMORY_HOST | GPUMemoryResource::MEMORY_CUSTOM, "ClusterNativeOutput");
}

void GPUTrackingInputProvider::SetMaxData(const GPUTrackingInOutPointers& io)
{
  mHoldTPCZS = io.tpcZS && (mRec->GetRecoStepsGPU() & GPUDataTypes::RecoStep::TPCClusterFinding);
  mHoldTPCClusterNative = (io.tpcZS || io.tpcPackedDigits || io.clustersNative) && mRec->IsGPU();
  mHoldTPCClusterNativeOutput = (io.tpcZS || io.tpcPackedDigits);
}
