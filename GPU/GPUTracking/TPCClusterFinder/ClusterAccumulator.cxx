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

/// \file ClusterAccumulator.cxx
/// \author Felix Weiglhofer

#include "ClusterAccumulator.h"
#include "CfUtils.h"
#include "ClusterNative.h"

#if !defined(__OPENCL__)
#include <cmath>
using namespace std;
#endif

using namespace GPUCA_NAMESPACE::gpu;

GPUd() void ClusterAccumulator::toNative(const deprecated::Digit& d, deprecated::ClusterNative& cn) const
{
  bool isEdgeCluster = CfUtils::isAtEdge(&d);
  bool wasSplitInTime = mSplitInTime >= MIN_SPLIT_NUM;
  bool wasSplitInPad = mSplitInPad >= MIN_SPLIT_NUM;

  uchar flags = 0;
  flags |= (isEdgeCluster) ? deprecated::CN_FLAG_IS_EDGE_CLUSTER : 0;
  flags |= (wasSplitInTime) ? deprecated::CN_FLAG_SPLIT_IN_TIME : 0;
  flags |= (wasSplitInPad) ? deprecated::CN_FLAG_SPLIT_IN_PAD : 0;

  cn.qmax = d.charge;
  cn.qtot = mQtot;
  deprecated::cnSetTimeFlags(&cn, mTimeMean, flags);
  deprecated::cnSetPad(&cn, mPadMean);
  deprecated::cnSetSigmaTime(&cn, mTimeSigma);
  deprecated::cnSetSigmaPad(&cn, mPadSigma);
}

GPUd() void ClusterAccumulator::update(Charge splitCharge, Delta dp, Delta dt)
{
  mQtot += splitCharge;
  mPadMean += splitCharge * dp;
  mTimeMean += splitCharge * dt;
  mPadSigma += splitCharge * dp * dp;
  mTimeSigma += splitCharge * dt * dt;
}

GPUd() Charge ClusterAccumulator::updateInner(PackedCharge charge, Delta dp, Delta dt)
{
  Charge q = charge.unpack();

  update(q, dp, dt);

  bool split = charge.isSplit();
  mSplitInTime += (dt != 0 && split);
  mSplitInPad += (dp != 0 && split);

  return q;
}

GPUd() Charge ClusterAccumulator::updateOuter(PackedCharge charge, Delta dp, Delta dt)
{
  Charge q = charge.unpack();

  bool split = charge.isSplit();
  bool has3x3 = charge.has3x3Peak();

  update((has3x3) ? 0.f : q, dp, dt);

  mSplitInTime += (dt != 0 && split && !has3x3);
  mSplitInPad += (dp != 0 && split && !has3x3);

  return q;
}

GPUd() void ClusterAccumulator::finalize(const deprecated::Digit& myDigit)
{
  mQtot += myDigit.charge;
  if (mQtot == 0) {
    return; // TODO: Why does this happen?
  }

  mPadMean /= mQtot;
  mTimeMean /= mQtot;
  mPadSigma /= mQtot;
  mTimeSigma /= mQtot;

  mPadSigma = sqrt(mPadSigma - mPadMean * mPadMean);
  mTimeSigma = sqrt(mTimeSigma - mTimeMean * mTimeMean);

  mPadMean += myDigit.pad;
  mTimeMean += myDigit.time;

#if defined(CORRECT_EDGE_CLUSTERS)
  if (CfUtils::isAtEdge(myDigit)) {
    float s = (myDigit->pad < 2) ? 1.f : -1.f;
    bool c = s * (mPadMean - myDigit->pad) > 0.f;
    mPadMean = (c) ? myDigit->pad : mPadMean;
  }
#endif
}
