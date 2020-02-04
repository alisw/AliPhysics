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

/// \file PackedCharge.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_PACKED_CHARGE_H
#define O2_GPU_PACKED_CHARGE_H

#include "clusterFinderDefs.h"
#include "GPUCommonMath.h"

// Ugly workaround because cuda doesn't like member constants
#define O2_GPU_PACKED_CHARGE_ADC_BITS 10
#define O2_GPU_PACKED_CHARGE_DECIMAL_BITS 4
#define O2_GPU_PACKED_CHARGE_CHARGE_BITS (O2_GPU_PACKED_CHARGE_ADC_BITS + O2_GPU_PACKED_CHARGE_DECIMAL_BITS)

namespace GPUCA_NAMESPACE
{
namespace gpu
{

namespace PackedChargeDefs
{
using BasicType = unsigned short;

static_assert(sizeof(BasicType) == 2);

GPUconstexpr() BasicType Has3x3PeakMask = 1 << O2_GPU_PACKED_CHARGE_CHARGE_BITS;
GPUconstexpr() BasicType IsSplitMask = 1 << (O2_GPU_PACKED_CHARGE_CHARGE_BITS + 1);
GPUconstexpr() BasicType ChargeMask = (1 << O2_GPU_PACKED_CHARGE_CHARGE_BITS) - 1;
GPUconstexpr() BasicType MaxVal = (1 << O2_GPU_PACKED_CHARGE_CHARGE_BITS) - 1;
GPUconstexpr() Charge Shift = 1 << O2_GPU_PACKED_CHARGE_DECIMAL_BITS;

#ifdef GPUCA_CPUCODE
static_assert(Has3x3PeakMask == 0x4000);
static_assert(IsSplitMask == 0x8000);
static_assert(ChargeMask == 0x3FFF);
static_assert(MaxVal == 0x3FFF);
static_assert(Shift == 16.f);
#endif

} // namespace PackedChargeDefs

class PackedCharge
{
 public:
  using BasicType = PackedChargeDefs::BasicType;

  PackedCharge() = default;
  GPUdi() explicit PackedCharge(Charge q) : PackedCharge(q, false, false) {}
  GPUdi() PackedCharge(Charge q, bool peak3x3, bool wasSplit)
  {
    using namespace PackedChargeDefs;

    mVal = q * Shift;
    mVal = CAMath::Min(MaxVal, mVal); // ensure only lower 14 bits are set
    mVal |= (peak3x3) ? Has3x3PeakMask : BasicType(0);
    mVal |= (wasSplit) ? IsSplitMask : BasicType(0);
  }

  GPUdi() Charge unpack() const { return Charge(mVal & PackedChargeDefs::ChargeMask) / PackedChargeDefs::Shift; }
  GPUdi() bool has3x3Peak() const { return mVal & PackedChargeDefs::Has3x3PeakMask; }
  GPUdi() bool isSplit() const { return mVal & PackedChargeDefs::IsSplitMask; }

 private:
  BasicType mVal;
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
