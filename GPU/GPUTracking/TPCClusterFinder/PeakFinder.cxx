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

/// \file PeakFinder.cxx
/// \author Felix Weiglhofer

#include "PeakFinder.h"

#include "Array2D.h"
#include "CfUtils.h"
#include "Digit.h"
#include "PackedCharge.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace GPUCA_NAMESPACE::gpu::deprecated;

GPUd() bool PeakFinder::isPeakScratchPad(
  GPUTPCClusterFinderKernels::GPUTPCSharedMemory& smem,
  const Digit* digit,
  ushort N,
  GPUglobalref() const PackedCharge* chargeMap,
  GPUsharedref() ChargePos* posBcast,
  GPUsharedref() PackedCharge* buf)
{
  ushort ll = get_local_id(0);

  const Timestamp time = digit->time;
  const Row row = digit->row;
  const Pad pad = digit->pad;

  const GlobalPad gpad = Array2D::tpcGlobalPadIdx(row, pad);
  ChargePos pos = {gpad, time};

  bool belowThreshold = (digit->charge <= QMAX_CUTOFF);

  ushort lookForPeaks;
  ushort partId = CfUtils::partition(
    smem,
    ll,
    belowThreshold,
    SCRATCH_PAD_WORK_GROUP_SIZE,
    &lookForPeaks);

  if (partId < lookForPeaks) {
    posBcast[partId] = pos;
  }
  GPUbarrier();

  CfUtils::fillScratchPad_PackedCharge(
    chargeMap,
    lookForPeaks,
    SCRATCH_PAD_WORK_GROUP_SIZE,
    ll,
    0,
    N,
    CfConsts::InnerNeighbors,
    posBcast,
    buf);

  if (belowThreshold) {
    return false;
  }

  bool peak = true;
  for (ushort i = 0; i < N; i++) {
    Charge q = buf[N * partId + i].unpack();
    peak &= (digit->charge > q) || (CfConsts::InnerTestEq[i] && digit->charge == q);
  }

  return peak;
}

GPUd() bool PeakFinder::isPeak(
  const Digit* digit,
  GPUglobalref() const PackedCharge* chargeMap)
{
  if (digit->charge <= QMAX_CUTOFF) {
    return false;
  }

  const Charge myCharge = digit->charge;
  const Timestamp time = digit->time;
  const Row row = digit->row;
  const Pad pad = digit->pad;

  const GlobalPad gpad = Array2D::tpcGlobalPadIdx(row, pad);

  bool peak = true;

#define CMP_NEIGHBOR(dp, dt, cmpOp)                           \
  do {                                                        \
    PackedCharge p = CHARGE(chargeMap, gpad + dp, time + dt); \
    const Charge otherCharge = p.unpack();                    \
    peak &= (otherCharge cmpOp myCharge);                     \
  } while (false)

#define CMP_LT CMP_NEIGHBOR(-1, -1, <=)
#define CMP_T CMP_NEIGHBOR(-1, 0, <=)
#define CMP_RT CMP_NEIGHBOR(-1, 1, <=)

#define CMP_L CMP_NEIGHBOR(0, -1, <=)
#define CMP_R CMP_NEIGHBOR(0, 1, <)

#define CMP_LB CMP_NEIGHBOR(1, -1, <)
#define CMP_B CMP_NEIGHBOR(1, 0, <)
#define CMP_RB CMP_NEIGHBOR(1, 1, <)

#if defined(CHARGEMAP_TILING_LAYOUT)
  CMP_LT;
  CMP_T;
  CMP_RT;
  CMP_R;
  CMP_RB;
  CMP_B;
  CMP_LB;
  CMP_L;
#else
  CMP_LT;
  CMP_T;
  CMP_RT;
  CMP_L;
  CMP_R;
  CMP_LB;
  CMP_B;
  CMP_RB;
#endif

#undef CMP_LT
#undef CMP_T
#undef CMP_RT
#undef CMP_L
#undef CMP_R
#undef CMP_LB
#undef CMP_B
#undef CMP_RB
#undef CMP_NEIGHBOR

  return peak;
}

GPUd() void PeakFinder::findPeaksImpl(int nBlocks, int nThreads, int iBlock, int iThread, GPUTPCClusterFinderKernels::GPUTPCSharedMemory& smem,
                                      GPUglobalref() const PackedCharge* chargeMap,
                                      GPUglobalref() const Digit* digits,
                                      uint digitnum,
                                      GPUglobalref() uchar* isPeakPredicate,
                                      GPUglobalref() uchar* peakMap)
{
  size_t idx = get_global_id(0);

  // For certain configurations dummy work items are added, so the total
  // number of work items is dividable by 64.
  // These dummy items also compute the last digit but discard the result.
  Digit myDigit = digits[CAMath::Min(idx, (size_t)(digitnum - 1))];

  uchar peak;
#if defined(BUILD_CLUSTER_SCRATCH_PAD)
  peak = isPeakScratchPad(smem, &myDigit, SCRATCH_PAD_SEARCH_N, chargeMap, smem.search.posBcast, smem.search.buf);
#else
  peak = isPeak(&myDigit, chargeMap);
#endif

  // Exit early if dummy. See comment above.
  bool iamDummy = (idx >= digitnum);
  if (iamDummy) {
    return;
  }

  isPeakPredicate[idx] = peak;

  const GlobalPad gpad = Array2D::tpcGlobalPadIdx(myDigit.row, myDigit.pad);

  IS_PEAK(peakMap, gpad, myDigit.time) =
    ((myDigit.charge > CHARGE_THRESHOLD) << 1) | peak;
}
