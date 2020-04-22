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

#define INVALID_TIME_BIN (-PADDING_TIME - 1)

struct ChargePos {
  tpccf::GlobalPad gpad;
  tpccf::TPCFragmentTime timePadded;

  GPUdDefault() ChargePos() CON_DEFAULT;

  GPUdi() ChargePos(tpccf::Row row, tpccf::Pad pad, tpccf::TPCFragmentTime t)
  {
    gpad = tpcGlobalPadIdx(row, pad);
    timePadded = t + PADDING_TIME;
  }

  GPUdi() ChargePos(const tpccf::GlobalPad& p, const tpccf::TPCFragmentTime& t) : gpad(p), timePadded(t) {}

  GPUdi() ChargePos delta(const tpccf::Delta2& d) const
  {
    return {tpccf::GlobalPad(gpad + d.x), tpccf::TPCFragmentTime(timePadded + d.y)};
  }

  GPUdi() bool valid() const { return timePadded >= 0; }

  GPUdi() tpccf::Row row() const { return gpad / TPC_PADS_PER_ROW_PADDED; }
  GPUdi() tpccf::Pad pad() const { return gpad % TPC_PADS_PER_ROW_PADDED - PADDING_PAD; }
  GPUdi() tpccf::TPCFragmentTime time() const { return timePadded - PADDING_TIME; }

  GPUdi() bool isPadding() const
  {
    tpccf::Pad pad = gpad % TPC_PADS_PER_ROW_PADDED;
    return timePadded < PADDING_TIME || timePadded >= TPC_MAX_FRAGMENT_LEN || pad < PADDING_PAD;
  }

 private:
  // Maps the position of a pad given as row and index in that row to a unique
  // index between 0 and TPC_NUM_OF_PADS.
  static GPUdi() tpccf::GlobalPad tpcGlobalPadIdx(tpccf::Row row, tpccf::Pad pad)
  {
    return TPC_PADS_PER_ROW_PADDED * row + pad + PADDING_PAD;
  }
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
