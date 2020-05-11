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
#include "TPCFastTransform.h"
#include "GPUTPCClusterData.h"
#include "GPUO2DataTypes.h"
#include "GPUTPCConvertImpl.h"

using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUTPCConvertKernel::Thread<0>(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() GPUSharedMemory& GPUrestrict() smem, processorType& GPUrestrict() processors)
{
  const int iSlice = iBlock / GPUCA_ROW_COUNT;
  const int iRow = iBlock % GPUCA_ROW_COUNT;
  GPUTPCConvert& GPUrestrict() convert = processors.tpcConverter;
  const o2::tpc::ClusterNativeAccess* GPUrestrict() native = processors.ioPtrs.clustersNative;
  GPUTPCClusterData* GPUrestrict() clusters = convert.mMemory->clusters[iSlice];
  const int idOffset = native->clusterOffset[iSlice][iRow];
  const int indexOffset = native->clusterOffset[iSlice][iRow] - native->clusterOffset[iSlice][0];

  for (unsigned int k = get_local_id(0); k < native->nClusters[iSlice][iRow]; k += get_local_size(0)) {
    const auto& GPUrestrict() clin = native->clusters[iSlice][iRow][k];
    float x, y, z;
    GPUTPCConvertImpl::convert(processors, iSlice, iRow, clin.getPad(), clin.getTime(), x, y, z);
    auto& GPUrestrict() clout = clusters[indexOffset + k];
    clout.x = x;
    clout.y = y;
    clout.z = z;
    clout.row = iRow;
    clout.amp = clin.qTot;
    clout.flags = clin.getFlags();
    clout.id = idOffset + k;
#ifdef GPUCA_TPC_RAW_PROPAGATE_PAD_ROW_TIME
    clout.pad = clin.getPad();
    clout.time = clin.getTime();
#endif
  }
}
