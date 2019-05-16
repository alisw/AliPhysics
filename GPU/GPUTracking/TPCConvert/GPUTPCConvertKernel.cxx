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

/// \file GPUTPCConvertKernel.cxx
/// \author David Rohr

#include "GPUTPCConvertKernel.h"
#include "GPUConstantMem.h"
#include "ClusterNativeAccessExt.h"
#include "TPCFastTransform.h"
#include "GPUTPCClusterData.h"

using namespace GPUCA_NAMESPACE::gpu;

#if defined(GPUCA_BUILD_TPCCONVERT) && !defined(GPUCA_ALIROOT_LIB) && defined(HAVE_O2HEADERS)
template <>
GPUd() void GPUTPCConvertKernel::Thread<0>(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUTPCSharedMemory& smem, processorType& processors)
{
  const int iSlice = iBlock / GPUCA_ROW_COUNT;
  const int iRow = iBlock % GPUCA_ROW_COUNT;
  GPUTPCConvert& convert = processors.tpcConverter;
  const ClusterNativeAccessExt* native = convert.mClustersNative;
  GPUTPCClusterData* clusters = convert.mMemory->clusters[iSlice];
  const TPCFastTransform* transform = convert.mTransform;
  const int continuousMaxTimeBin = processors.param.continuousMaxTimeBin;
  const int idOffset = native->clusterOffset[iSlice][iRow];
  const int indexOffset = native->clusterOffset[iSlice][iRow] - native->clusterOffset[iSlice][0];

  for (unsigned int k = get_local_id(0); k < native->nClusters[iSlice][iRow]; k += get_local_size(0)) {
    const auto& cin = native->clusters[iSlice][iRow][k];
    float x = 0, y = 0, z = 0;
    if (continuousMaxTimeBin == 0) {
      transform->Transform(iSlice, iRow, cin.getPad(), cin.getTime(), x, y, z);
    } else {
      transform->TransformInTimeFrame(iSlice, iRow, cin.getPad(), cin.getTime(), x, y, z, continuousMaxTimeBin);
    }
    auto& cout = clusters[indexOffset + k];
    cout.x = x;
    cout.y = y;
    cout.z = z;
    cout.row = iRow;
    cout.amp = cin.qMax;
    cout.flags = cin.getFlags();
    cout.id = idOffset + k;
  }
}
#endif
