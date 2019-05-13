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

/// \file GPUITSTrack.h
/// \author David Rohr, Maximiliano Puccio

#ifndef GPUITSTRACK_H
#define GPUITSTRACK_H

#include "GPUTPCGMTrackParam.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUITSTrack : public GPUTPCGMTrackParam
{
 public:
  GPUTPCGMTrackParam::GPUTPCOuterParam mOuterParam;
  float mAlpha;
  int mClusters[7];
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
