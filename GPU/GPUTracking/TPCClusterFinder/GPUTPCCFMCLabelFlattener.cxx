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

/// \file GPUTPCCFMCLabelFlattener.cxx
/// \author Felix Weiglhofer

#include "GPUTPCCFMCLabelFlattener.h"

#if !defined(GPUCA_GPUCODE)
#include "GPUHostDataTypes.h"
#endif

using namespace GPUCA_NAMESPACE::gpu;

#if !defined(GPUCA_GPUCODE)
void GPUTPCCFMCLabelFlattener::setGlobalOffsetsAndAllocate(
  GPUTPCClusterFinder& cls,
  GPUTPCLinearLabels& labels)
{
  uint headerOffset = labels.header.size();
  uint dataOffset = labels.data.size();

  for (Row row = 0; row < GPUCA_ROW_COUNT; row++) {
    cls.mPlabelHeaderOffset[row] = headerOffset;
    headerOffset += cls.mPclusterInRow[row];

    cls.mPlabelDataOffset[row] = dataOffset;
    auto& lastInterim = cls.mPlabelsByRow[cls.mNMaxClusterPerRow * row + cls.mPclusterInRow[row] - 1];
    uint labelsInRow = lastInterim.offset + lastInterim.labels.size();
    dataOffset += labelsInRow;
  }

  labels.header.resize(headerOffset);
  labels.data.resize(dataOffset);
}
#endif

template <>
GPUd() void GPUTPCCFMCLabelFlattener::Thread<GPUTPCCFMCLabelFlattener::setRowOffsets>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory&, processorType& clusterer)
{
#if !defined(GPUCA_GPUCODE)
  Row row = get_global_id(0);

  uint clusterInRow = clusterer.mPclusterInRow[row];
  uint offset = 0;

  for (uint i = 0; i < clusterInRow; i++) {
    auto& interim = clusterer.mPlabelsByRow[row * clusterer.mNMaxClusterPerRow + i];
    interim.offset = offset;
    offset += interim.labels.size();
  }
#endif
}

template <>
GPUd() void GPUTPCCFMCLabelFlattener::Thread<GPUTPCCFMCLabelFlattener::flatten>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory&, processorType& clusterer, uint row, GPUTPCLinearLabels* out)
{
#if !defined(GPUCA_GPUCODE)
  uint idx = get_global_id(0);

  GPUTPCClusterMCInterim& interim = clusterer.mPlabelsByRow[row * clusterer.mNMaxClusterPerRow + idx];

  uint headerOffset = clusterer.mPlabelHeaderOffset[row] + idx;
  uint dataOffset = clusterer.mPlabelDataOffset[row] + interim.offset;

  out->header[headerOffset] = dataOffset;
  std::copy(interim.labels.begin(), interim.labels.end(), &out->data[dataOffset]);
#endif
}
