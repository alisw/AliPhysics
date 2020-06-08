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

/// \file GPUTPCGMMergerTypes.h
/// \author David Rohr

#ifndef GPUTPCGMMERGERTYPES_H
#define GPUTPCGMMERGERTYPES_H

#include "GPUTPCDef.h"
#include "GPUGeneralKernels.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
namespace GPUTPCGMMergerTypes
{

enum attachTypes { attachAttached = 0x40000000,
                   attachGood = 0x20000000,
                   attachGoodLeg = 0x10000000,
                   attachTube = 0x08000000,
                   attachHighIncl = 0x04000000,
                   attachTrackMask = 0x03FFFFFF,
                   attachFlagMask = 0xFC000000 };

struct InterpolationErrorHit {
  float posY;
  float errorY;
  float posZ;
  float errorZ;
};

struct InterpolationErrors {
  InterpolationErrorHit hit[GPUCA_MERGER_MAX_TRACK_CLUSTERS];
};

struct GPUResolveSharedMemory : public GPUKernelTemplate::GPUSharedMemoryScan64<short, GPUCA_GET_THREAD_COUNT(GPUCA_LB_GPUTPCGMMergerResolve_step3)> {
  int iTrack1[GPUCA_GET_THREAD_COUNT(GPUCA_LB_GPUTPCGMMergerResolve_step3)];
  int iTrack2[GPUCA_GET_THREAD_COUNT(GPUCA_LB_GPUTPCGMMergerResolve_step3)];
};

struct GPUTPCGMBorderRange {
  int fId;
  float fMin, fMax;
};

} // namespace GPUTPCGMMergerTypes
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
