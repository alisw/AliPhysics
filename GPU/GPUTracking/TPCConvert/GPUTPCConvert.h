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

namespace o2
{
namespace tpc
{
struct ClusterNativeAccessFullTPC;
struct ClusterNative;
} // namespace tpc
} // namespace o2

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUTPCClusterData;
struct ClusterNativeAccessExt;
class TPCFastTransform;

class GPUTPCConvert : public GPUProcessor
{
  friend class GPUTPCConvertKernel;
  friend class GPUChainTracking;

 public:
#ifndef GPUCA_GPUCODE
  void InitializeProcessor();
  void RegisterMemoryAllocation();
  void SetMaxData();

  void* SetPointersInput(void* mem);
  void* SetPointersOutput(void* mem);
  void* SetPointersMemory(void* mem);

  void set(ClusterNativeAccessExt* clustersNative, const TPCFastTransform* transform)
  {
    mClustersNative = clustersNative;
    mTransform = transform;
  }
#endif
  GPUd() const ClusterNativeAccessExt* getClustersNative()
  {
    return mClustersNative;
  }

  constexpr static unsigned int NSLICES = GPUCA_NSLICES;

  struct Memory {
    GPUTPCClusterData* clusters[NSLICES];
  };

 protected:
  ClusterNativeAccessExt* mClustersNative = nullptr;
  ClusterNativeAccessExt* mClustersNativeBuffer;

  const TPCFastTransform* mTransform = nullptr;
  Memory* mMemory = nullptr;
  o2::tpc::ClusterNative* mInputClusters;
  GPUTPCClusterData* mClusters = nullptr;
  unsigned int mNClustersTotal = 0;

  short mMemoryResInput = -1;
  short mMemoryResOutput = -1;
  short mMemoryResMemory = -1;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
