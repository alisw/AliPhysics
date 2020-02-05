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

/// \file GPUTPCClusterFinder.cxx
/// \author David Rohr

#include "GPUTPCClusterFinder.h"
#include "GPUReconstruction.h"

#include "DataFormatsTPC/ZeroSuppression.h"
#include "Digit.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace o2::tpc;

void GPUTPCClusterFinder::InitializeProcessor() {}

void* GPUTPCClusterFinder::SetPointersMemory(void* mem)
{
  computePointerWithAlignment(mem, mPmemory, 1);
  return mem;
}

void* GPUTPCClusterFinder::SetPointersInput(void* mem)
{
  if (mNMaxPages && (mRec->GetRecoStepsGPU() & GPUDataTypes::RecoStep::TPCClusterFinding)) {
    computePointerWithAlignment(mem, mPzs, mNMaxPages * TPCZSHDR::TPC_ZS_PAGE_SIZE);
  }
  if (mNMaxPages || (mRec->GetRecoStepsGPU() & GPUDataTypes::RecoStep::TPCClusterFinding)) {
    computePointerWithAlignment(mem, mPdigits, mNMaxDigits);
  }
  return mem;
}

void* GPUTPCClusterFinder::SetPointersOutput(void* mem)
{
  computePointerWithAlignment(mem, mPclusterInRow, GPUCA_ROW_COUNT);
  computePointerWithAlignment(mem, mPclusterByRow, GPUCA_ROW_COUNT * mNMaxClusterPerRow);
  return mem;
}

void* GPUTPCClusterFinder::SetPointersScratch(void* mem)
{
  computePointerWithAlignment(mem, mPpeaks, mNMaxPeaks);
  computePointerWithAlignment(mem, mPfilteredPeaks, mNMaxPeaks);
  computePointerWithAlignment(mem, mPisPeak, mNMaxDigits);
  computePointerWithAlignment(mem, mPchargeMap, TPC_NUM_OF_PADS * TPC_MAX_TIME_PADDED);
  computePointerWithAlignment(mem, mPpeakMap, TPC_NUM_OF_PADS * TPC_MAX_TIME_PADDED);
  computePointerWithAlignment(mem, mPbuf, mBufSize * mNBufs);
  return mem;
}

void GPUTPCClusterFinder::RegisterMemoryAllocation()
{
  mRec->RegisterMemoryAllocation(this, &GPUTPCClusterFinder::SetPointersInput, GPUMemoryResource::MEMORY_INPUT | GPUMemoryResource::MEMORY_GPU, "TPCClustererInput");
  mRec->RegisterMemoryAllocation(this, &GPUTPCClusterFinder::SetPointersScratch, GPUMemoryResource::MEMORY_SCRATCH, "TPCClustererScratch", GPUMemoryReuse{GPUMemoryReuse::REUSE_1TO1, GPUMemoryReuse::ClustererScratch, (unsigned short)(mISlice % mRec->GetDeviceProcessingSettings().nTPCClustererLanes)});
  mMemoryId = mRec->RegisterMemoryAllocation(this, &GPUTPCClusterFinder::SetPointersMemory, GPUMemoryResource::MEMORY_PERMANENT, "TPCClustererMemory");
  mRec->RegisterMemoryAllocation(this, &GPUTPCClusterFinder::SetPointersOutput, GPUMemoryResource::MEMORY_OUTPUT, "TPCClustererOutput");
}

void GPUTPCClusterFinder::SetMaxData(const GPUTrackingInOutPointers& io)
{
  mNMaxPeaks = 0.5f * mNMaxDigits;
  mNMaxClusters = 0.2f * mNMaxPeaks;
  mNMaxClusterPerRow = 0.01f * mNMaxDigits;
  mBufSize = nextMultipleOf<std::max<int>(GPUCA_MEMALIGN_SMALL, mScanWorkGroupSize)>(mNMaxDigits);
  mNBufs = getNSteps(mBufSize);
}

void GPUTPCClusterFinder::SetNMaxDigits(size_t nDigits, size_t nPages)
{
  mNMaxDigits = nextMultipleOf<std::max<int>(GPUCA_MEMALIGN_SMALL, mScanWorkGroupSize)>(nDigits);
  mNMaxPages = nPages;
}

size_t GPUTPCClusterFinder::getNSteps(size_t items) const
{
  if (items == 0) {
    return 0;
  }
  size_t c = 1;
  size_t capacity = mScanWorkGroupSize;
  while (items > capacity) {
    capacity *= mScanWorkGroupSize;
    c++;
  }
  return c;
}
