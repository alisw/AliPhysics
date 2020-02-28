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

/// \file GPUTPCConvert.cxx
/// \author David Rohr

#include "GPUTPCConvert.h"
#include "TPCFastTransform.h"
#include "GPUTPCClusterData.h"
#include "GPUReconstruction.h"
#include "GPUO2DataTypes.h"

using namespace GPUCA_NAMESPACE::gpu;

void GPUTPCConvert::InitializeProcessor() {}

void* GPUTPCConvert::SetPointersInput(void* mem)
{
  computePointerWithAlignment(mem, mInputClusters, mNClustersTotal);
  return mem;
}

void* GPUTPCConvert::SetPointersOutput(void* mem)
{
  computePointerWithAlignment(mem, mClusters, mNClustersTotal);
  return mem;
}

void* GPUTPCConvert::SetPointersMemory(void* mem)
{
  computePointerWithAlignment(mem, mMemory, 1);
  if (mRec->GetRecoStepsGPU() & GPUDataTypes::RecoStep::TPCConversion) {
    computePointerWithAlignment(mem, mClustersNativeBuffer, 1);
  }
  return mem;
}

void GPUTPCConvert::RegisterMemoryAllocation()
{
  mMemoryResMemory = mRec->RegisterMemoryAllocation(this, &GPUTPCConvert::SetPointersMemory, GPUMemoryResource::MEMORY_INPUT | GPUMemoryResource::MEMORY_PERMANENT, "TPCConvertMemory");
  mMemoryResInput = mRec->RegisterMemoryAllocation(this, &GPUTPCConvert::SetPointersInput, GPUMemoryResource::MEMORY_INPUT | GPUMemoryResource::MEMORY_EXTERNAL | GPUMemoryResource::MEMORY_CUSTOM_TRANSFER, "TPCConvertInput");
  mMemoryResOutput = mRec->RegisterMemoryAllocation(this, &GPUTPCConvert::SetPointersOutput, GPUMemoryResource::MEMORY_OUTPUT, "TPCConvertOutput");
}

void GPUTPCConvert::SetMaxData(const GPUTrackingInOutPointers& io)
{
  if (io.tpcPackedDigits || io.tpcZS) {
    // TODO: Don't do anything for now, set from the outside, should be fixed
  } else if (mClustersNative) {
    mNClustersTotal = mClustersNative->nClustersTotal;
  } else {
    mNClustersTotal = 0;
  }
}
