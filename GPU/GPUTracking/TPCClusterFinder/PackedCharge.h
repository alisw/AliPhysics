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

namespace GPUCA_NAMESPACE
{
namespace gpu
{

class PackedCharge
{
 public:
  using BasicType = unsigned short;
  static_assert(sizeof(BasicType) == 2);

  enum Constants : BasicType {
    Null = 0,
    ADCBits = 10,
    DecimalBits = 4,
    ChargeBits = ADCBits + DecimalBits,
    ChargeMask = (1 << ChargeBits) - 1,
    MaxVal = ChargeMask,
    Has3x3PeakMask = 1 << ChargeBits,
    IsSplitMask = 1 << (ChargeBits + 1),
  };

  GPUdDefault() PackedCharge() CON_DEFAULT;
  GPUdi() explicit PackedCharge(Charge q) : PackedCharge(q, false, false) {}
  GPUdi() PackedCharge(Charge q, bool peak3x3, bool wasSplit)
  {
    mVal = q * Charge(1 << DecimalBits);
    mVal = CAMath::Min<BasicType>(MaxVal, mVal); // ensure only lower 14 bits are set
    mVal |= (peak3x3) ? Has3x3PeakMask : Null;
    mVal |= (wasSplit) ? IsSplitMask : Null;
  }

  GPUdi() Charge unpack() const { return Charge(mVal & ChargeMask) / Charge(1 << DecimalBits); }
  GPUdi() bool has3x3Peak() const { return mVal & Has3x3PeakMask; }
  GPUdi() bool isSplit() const { return mVal & IsSplitMask; }

 private:
  BasicType mVal;
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
