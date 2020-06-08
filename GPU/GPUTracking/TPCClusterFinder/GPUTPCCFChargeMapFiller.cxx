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
#include "DataFormatsTPC/Digit.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace GPUCA_NAMESPACE::gpu::tpccf;

template <>
GPUdii() void GPUTPCCFChargeMapFiller::Thread<GPUTPCCFChargeMapFiller::fillIndexMap>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer)
{
  Array2D<uint> indexMap(clusterer.mPindexMap);
  fillIndexMapImpl(get_num_groups(0), get_local_size(0), get_group_id(0), get_local_id(0), clusterer.mPmemory->fragment, clusterer.mPdigits, indexMap, clusterer.mPmemory->counters.nPositions);
}

GPUd() void GPUTPCCFChargeMapFiller::fillIndexMapImpl(int nBlocks, int nThreads, int iBlock, int iThread,
                                                      const CfFragment& fragment,
                                                      const tpc::Digit* digits,
                                                      Array2D<uint>& indexMap,
                                                      size_t maxDigit)
{
  size_t idx = get_global_id(0);
  if (idx >= maxDigit) {
    return;
  }
  CPU_ONLY(idx += fragment.digitsStart);
  CPU_ONLY(tpc::Digit digit = digits[idx]);
  CPU_ONLY(ChargePos pos(digit.getRow(), digit.getPad(), fragment.toLocal(digit.getTimeStamp())));
  CPU_ONLY(indexMap.safeWrite(pos, idx));
}

template <>
GPUdii() void GPUTPCCFChargeMapFiller::Thread<GPUTPCCFChargeMapFiller::fillFromDigits>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer)
{
  Array2D<PackedCharge> chargeMap(reinterpret_cast<PackedCharge*>(clusterer.mPchargeMap));
  fillFromDigitsImpl(get_num_groups(0), get_local_size(0), get_group_id(0), get_local_id(0), clusterer.mPmemory->fragment, clusterer.mPmemory->counters.nPositions, clusterer.mPdigits, clusterer.mPpositions, chargeMap);
}

GPUd() void GPUTPCCFChargeMapFiller::fillFromDigitsImpl(int nBlocks, int nThreads, int iBlock, int iThread, const CfFragment& fragment, size_t digitNum,
                                                        const tpc::Digit* digits,
                                                        ChargePos* positions,
                                                        Array2D<PackedCharge>& chargeMap)
{
  size_t idx = get_global_id(0);
  if (idx >= digitNum) {
    return;
  }
  tpc::Digit digit = digits[fragment.digitsStart + idx];

  ChargePos pos(digit.getRow(), digit.getPad(), fragment.toLocal(digit.getTimeStamp()));
  positions[idx] = pos;
  chargeMap[pos] = PackedCharge(digit.getChargeFloat());
}

template <>
GPUdii() void GPUTPCCFChargeMapFiller::Thread<GPUTPCCFChargeMapFiller::findFragmentStart>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer)
{
  // TODO: use binary search
  if (get_global_id(0) != 0) {
    return;
  }

  size_t nDigits = clusterer.mPmemory->counters.nDigits;
  const tpc::Digit* digits = clusterer.mPdigits;
  size_t st = 0;
  for (; st < nDigits && digits[st].getTimeStamp() < clusterer.mPmemory->fragment.first(); st++) {
  }

  size_t end = st;
  for (; end < nDigits && digits[end].getTimeStamp() < clusterer.mPmemory->fragment.last(); end++) {
  }

  clusterer.mPmemory->fragment.digitsStart = st;
  clusterer.mPmemory->counters.nPositions = end - st;
}
