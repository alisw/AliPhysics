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

/// \file ClusterAccumulator.h
/// \author Felix Weiglhofer

#ifndef O2_GPU_CLUSTER_ACCUMULATOR_H
#define O2_GPU_CLUSTER_ACCUMULATOR_H

#include "clusterFinderDefs.h"
#include "PackedCharge.h"

namespace GPUCA_NAMESPACE
{

namespace tpc
{
struct ClusterNative;
}

namespace gpu
{

class ClusterAccumulator
{

 public:
  GPUd() Charge updateInner(PackedCharge, Delta2);
  GPUd() Charge updateOuter(PackedCharge, Delta2);

  GPUd() void finalize(const deprecated::Digit&);
  GPUd() void toNative(const deprecated::Digit&, tpc::ClusterNative&) const;

 private:
  float mQtot = 0;
  float mPadMean = 0;
  float mPadSigma = 0;
  float mTimeMean = 0;
  float mTimeSigma = 0;
  uchar mSplitInTime = 0;
  uchar mSplitInPad = 0;

  GPUd() void update(Charge, Delta2);
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
