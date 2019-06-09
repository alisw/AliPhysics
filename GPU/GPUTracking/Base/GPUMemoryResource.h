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

/// \file GPUMemoryResource.h
/// \author David Rohr

#ifndef GPUMEMORYRESOURCE_H
#define GPUMEMORYRESOURCE_H

#include "GPUCommonDef.h"
#include "GPUProcessor.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUMemoryResource
{
  friend class GPUReconstruction;
  friend class GPUReconstructionCPU;

 public:
  enum MemoryType {
    MEMORY_HOST = 1,
    MEMORY_GPU = 2,
    MEMORY_INPUT_FLAG = 4,
    MEMORY_INPUT = 7,
    MEMORY_OUTPUT_FLAG = 8,
    MEMORY_OUTPUT = 11,
    MEMORY_INOUT = 15,
    MEMORY_SCRATCH = 16,
    MEMORY_SCRATCH_HOST = 17,
    MEMORY_EXTERNAL = 32,
    MEMORY_PERMANENT = 64,
    MEMORY_CUSTOM = 128,
    MEMORY_CUSTOM_TRANSFER = 256
  };
  enum AllocationType { ALLOCATION_AUTO = 0,
                        ALLOCATION_INDIVIDUAL = 1,
                        ALLOCATION_GLOBAL = 2 };

#ifndef GPUCA_GPUCODE
  GPUMemoryResource(GPUProcessor* proc, void* (GPUProcessor::*setPtr)(void*), MemoryType type, const char* name = "") : mProcessor(proc), mPtr(nullptr), mPtrDevice(nullptr), mSetPointers(setPtr), mType(type), mSize(0), mName(name)
  {
  }
  GPUMemoryResource(const GPUMemoryResource&) CON_DEFAULT;
#endif

#ifndef __OPENCL__
  void* SetPointers(void* ptr)
  {
    return (mProcessor->*mSetPointers)(ptr);
  }
  void* SetDevicePointers(void* ptr) { return (mProcessor->mDeviceProcessor->*mSetPointers)(ptr); }
  void* Ptr() { return mPtr; }
  void* PtrDevice() { return mPtrDevice; }
  size_t Size() const { return mSize; }
  const char* Name() const { return mName; }
  MemoryType Type() const { return mType; }
#endif

 private:
  GPUProcessor* mProcessor;
  void* mPtr;
  void* mPtrDevice;
  void* (GPUProcessor::*mSetPointers)(void*);
  MemoryType mType;
  size_t mSize;
  const char* mName;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
