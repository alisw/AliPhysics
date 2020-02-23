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

/// \file ChargePos.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_CHARGE_POS_H
#define O2_GPU_CHARGE_POS_H

#include "clusterFinderDefs.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{

struct ChargePos {
  GlobalPad gpad;
  Timestamp time;

  GPUdDefault() ChargePos() CON_DEFAULT;

  GPUdi() explicit ChargePos(const deprecated::Digit& d)
  {
    gpad = tpcGlobalPadIdx(d.row, d.pad);
    time = d.time + PADDING_TIME;
  }

  GPUdi() ChargePos(const GlobalPad& p, const Timestamp& t) : gpad(p), time(t) {}

  GPUdi() ChargePos delta(const Delta2& d) const
  {
    return {GlobalPad(gpad + d.x), Timestamp(time + d.y)};
  }

 private:
  // Maps the position of a pad given as row and index in that row to a unique
  // index between 0 and TPC_NUM_OF_PADS.
  static GPUdi() GlobalPad tpcGlobalPadIdx(Row row, Pad pad)
  {
    return TPC_PADS_PER_ROW_PADDED * row + pad + PADDING_PAD;
  }
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
