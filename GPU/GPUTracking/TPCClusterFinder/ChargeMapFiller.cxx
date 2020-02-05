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

#include "ChargeMapFiller.h"
#include "ChargePos.h"
#include "Array2D.h"

using namespace GPUCA_NAMESPACE::gpu;
using namespace GPUCA_NAMESPACE::gpu::deprecated;

GPUd() void ChargeMapFiller::fillChargeMapImpl(int nBlocks, int nThreads, int iBlock, int iThread, GPUTPCClusterFinderKernels::GPUTPCSharedMemory& smem,
                                               const Digit* digits,
                                               Array2D<PackedCharge>& chargeMap,
                                               size_t maxDigit)
{
  size_t idx = get_global_id(0);
  if (idx >= maxDigit) {
    return;
  }
  Digit myDigit = digits[idx];

  chargeMap[ChargePos(myDigit)] = PackedCharge(myDigit.charge);
}

GPUd() void ChargeMapFiller::resetMapsImpl(int nBlocks, int nThreads, int iBlock, int iThread, GPUTPCClusterFinderKernels::GPUTPCSharedMemory& smem,
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
