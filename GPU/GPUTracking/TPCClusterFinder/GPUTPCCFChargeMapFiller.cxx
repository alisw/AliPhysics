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

/// \file ChargeMapFiller.cxx
/// \author Felix Weiglhofer

#include "GPUTPCCFChargeMapFiller.h"
#include "ChargePos.h"
#include "Array2D.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace GPUCA_NAMESPACE::gpu::deprecated;

template <>
GPUdii() void GPUTPCCFChargeMapFiller::Thread<GPUTPCCFChargeMapFiller::fillChargeMap>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer)
{
  Array2D<PackedCharge> chargeMap(reinterpret_cast<PackedCharge*>(clusterer.mPchargeMap));
  Array2D<uint> indexMap(clusterer.mPindexMap);
  fillChargeMapImpl(get_num_groups(0), get_local_size(0), get_group_id(0), get_local_id(0), clusterer.mPdigits, chargeMap, indexMap, clusterer.mPmemory->counters.nDigits);
}

GPUd() void GPUTPCCFChargeMapFiller::fillChargeMapImpl(int nBlocks, int nThreads, int iBlock, int iThread,
                                                       const Digit* digits,
                                                       Array2D<PackedCharge>& chargeMap,
                                                       Array2D<uint>& indexMap,
                                                       size_t maxDigit)
{
  size_t idx = get_global_id(0);
  if (idx >= maxDigit) {
    return;
  }
  Digit myDigit = digits[idx];

  ChargePos pos(myDigit);
  chargeMap[pos] = PackedCharge(myDigit.charge);
  CPU_ONLY(indexMap.safeWrite(pos, idx));
}

template <>
GPUdii() void GPUTPCCFChargeMapFiller::Thread<GPUTPCCFChargeMapFiller::resetMaps>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer)
{
  Array2D<PackedCharge> chargeMap(reinterpret_cast<PackedCharge*>(clusterer.mPchargeMap));
  Array2D<uchar> isPeakMap(clusterer.mPpeakMap);
  resetMapsImpl(get_num_groups(0), get_local_size(0), get_group_id(0), get_local_id(0), clusterer.mPdigits, chargeMap, isPeakMap);
}

GPUd() void GPUTPCCFChargeMapFiller::resetMapsImpl(int nBlocks, int nThreads, int iBlock, int iThread,
                                                   const Digit* digits,
                                                   Array2D<PackedCharge>& chargeMap,
                                                   Array2D<uchar>& isPeakMap)
{
  size_t idx = get_global_id(0);
  Digit myDigit = digits[idx];

  ChargePos pos(myDigit);

  chargeMap[pos] = PackedCharge(0);
  isPeakMap[pos] = 0;
}
