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

/// \file GPUTPCCFGather.cxx
/// \author David Rohr

#include "GPUTPCCFGather.h"
using namespace GPUCA_NAMESPACE::gpu;
using namespace GPUCA_NAMESPACE::gpu::tpccf;

template <>
GPUdii() void GPUTPCCFGather::Thread<0>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer, o2::tpc::ClusterNative* ptr)
{
  for (int i = 0; i < iBlock; i++) {
    ptr += clusterer.mPclusterInRow[i];
  }
  for (unsigned int i = iThread; i < clusterer.mPclusterInRow[iBlock]; i += nThreads) {
    ptr[i] = clusterer.mPclusterByRow[iBlock * clusterer.mNMaxClusterPerRow + i];
  }
}
