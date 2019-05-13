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

/// \file GPUProcessor.cxx
/// \author David Rohr

#include "GPUProcessor.h"
#include "GPUReconstruction.h"
#include "GPUReconstructionDeviceBase.h"

using namespace GPUCA_NAMESPACE::gpu;

GPUProcessor::GPUProcessor() : mRec(nullptr), mGPUProcessorType(PROCESSOR_TYPE_CPU), mDeviceProcessor(nullptr), mCAParam(nullptr), mAllocateAndInitializeLate(false) {}

GPUProcessor::~GPUProcessor()
{
  if (mRec && mRec->GetDeviceProcessingSettings().memoryAllocationStrategy == GPUMemoryResource::ALLOCATION_INDIVIDUAL) {
    Clear();
  }
}

void GPUProcessor::InitGPUProcessor(GPUReconstruction* rec, GPUProcessor::ProcessorType type, GPUProcessor* slaveProcessor)
{
  mRec = rec;
  mGPUProcessorType = type;
  if (slaveProcessor) {
    slaveProcessor->mDeviceProcessor = this;
  }
  mCAParam = type == PROCESSOR_TYPE_DEVICE ? (reinterpret_cast<GPUReconstructionDeviceBase*>(rec))->DeviceParam() : &rec->GetParam();
}

void GPUProcessor::Clear() { mRec->FreeRegisteredMemory(this, true); }
