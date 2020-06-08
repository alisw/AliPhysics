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

/// \file GPUTPCCFDeconvolution.cxx
/// \author Felix Weiglhofer

#include "GPUTPCCFDeconvolution.h"
#include "CfConsts.h"
#include "CfUtils.h"
#include "ChargePos.h"
#include "GPUDefMacros.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace GPUCA_NAMESPACE::gpu::tpccf;

template <>
GPUdii() void GPUTPCCFDeconvolution::Thread<0>(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem, processorType& clusterer)
{
  Array2D<PackedCharge> chargeMap(reinterpret_cast<PackedCharge*>(clusterer.mPchargeMap));
  Array2D<uchar> isPeakMap(clusterer.mPpeakMap);
  GPUTPCCFDeconvolution::countPeaksImpl(get_num_groups(0), get_local_size(0), get_group_id(0), get_local_id(0), smem, isPeakMap, chargeMap, clusterer.mPpositions, clusterer.mPmemory->counters.nPositions);
}

GPUdii() void GPUTPCCFDeconvolution::countPeaksImpl(int nBlocks, int nThreads, int iBlock, int iThread, GPUSharedMemory& smem,
                                                    const Array2D<uchar>& peakMap,
                                                    Array2D<PackedCharge>& chargeMap,
                                                    const ChargePos* positions,
                                                    const uint digitnum)
{
  SizeT idx = get_global_id(0);

  bool iamDummy = (idx >= digitnum);
  idx = iamDummy ? digitnum - 1 : idx;

  ChargePos pos = positions[idx];

  bool iamPeak = CfUtils::isPeak(peakMap[pos]);

  char peakCount = (iamPeak) ? 1 : 0;

  ushort ll = get_local_id(0);
  ushort partId = ll;

  ushort in3x3 = 0;
  partId = CfUtils::partition(smem, ll, iamPeak, SCRATCH_PAD_WORK_GROUP_SIZE, &in3x3);

  if (partId < in3x3) {
    smem.posBcast1[partId] = pos;
  }
  GPUbarrier();

  CfUtils::blockLoad(
    peakMap,
    in3x3,
    SCRATCH_PAD_WORK_GROUP_SIZE,
    ll,
    0,
    8,
    CfConsts::InnerNeighbors,
    smem.posBcast1,
    smem.buf);

  uchar aboveThreshold = 0;
  if (partId < in3x3) {
    peakCount = countPeaksScratchpadInner(partId, smem.buf, &aboveThreshold);
  }

  ushort in5x5 = 0;
  partId = CfUtils::partition(smem, partId, peakCount > 0 && !iamPeak, in3x3, &in5x5);

  if (partId < in5x5) {
    smem.posBcast1[partId] = pos;
    smem.aboveThresholdBcast[partId] = aboveThreshold;
  }
  GPUbarrier();

  CfUtils::condBlockLoad<uchar, true>(
    peakMap,
    in5x5,
    SCRATCH_PAD_WORK_GROUP_SIZE,
    ll,
    0,
    16,
    CfConsts::OuterNeighbors,
    smem.posBcast1,
    smem.aboveThresholdBcast,
    smem.buf);

  if (partId < in5x5) {
    peakCount = countPeaksScratchpadOuter(partId, aboveThreshold, smem.buf);
    peakCount *= -1;
  }

  if (iamDummy) {
    return;
  }

  bool has3x3 = (peakCount > 0);
  peakCount = CAMath::Abs(int(peakCount));
  bool split = (peakCount > 1);

  peakCount = (peakCount == 0) ? 1 : peakCount;

  PackedCharge charge = chargeMap[pos];
  PackedCharge p(charge.unpack() / peakCount, has3x3, split);

  chargeMap[pos] = p;
}

GPUd() char GPUTPCCFDeconvolution::countPeaksScratchpadInner(
  ushort ll,
  const uchar* isPeak,
  uchar* aboveThreshold)
{
  char peaks = 0;
  GPUCA_UNROLL(U(), U())
  for (uchar i = 0; i < 8; i++) {
    uchar p = isPeak[ll * 8 + i];
    peaks += CfUtils::isPeak(p);
    *aboveThreshold |= uchar(CfUtils::isAboveThreshold(p)) << i;
  }

  return peaks;
}

GPUd() char GPUTPCCFDeconvolution::countPeaksScratchpadOuter(
  ushort ll,
  uchar aboveThreshold,
  const uchar* isPeak)
{
  char peaks = 0;
  GPUCA_UNROLL(U(), U())
  for (uchar i = 0; i < 16; i++) {
    uchar p = isPeak[ll * 16 + i];
    peaks += CfUtils::isPeak(p);
  }

  return peaks;
}
