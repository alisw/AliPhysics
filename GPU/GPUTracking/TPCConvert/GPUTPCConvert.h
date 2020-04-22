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

/// \file GPUTPCConvert.h
/// \author David Rohr

#ifndef GPUTPCCONVERT_H
#define GPUTPCCONVERT_H

#include "GPUDef.h"
#include "GPUProcessor.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUTPCClusterData;

class GPUTPCConvert : public GPUProcessor
{
  friend class GPUTPCConvertKernel;
  friend class GPUChainTracking;

 public:
#ifndef GPUCA_GPUCODE
  void InitializeProcessor();
  void RegisterMemoryAllocation();
  void SetMaxData(const GPUTrackingInOutPointers& io);

  void* SetPointersOutput(void* mem);
  void* SetPointersMemory(void* mem);
#endif

  constexpr static unsigned int NSLICES = GPUCA_NSLICES;

  struct Memory {
    GPUTPCClusterData* clusters[NSLICES];
  };

 protected:
  Memory* mMemory = nullptr;
  GPUTPCClusterData* mClusters = nullptr;
  unsigned int mNClustersTotal = 0;

  short mMemoryResOutput = -1;
  short mMemoryResMemory = -1;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
