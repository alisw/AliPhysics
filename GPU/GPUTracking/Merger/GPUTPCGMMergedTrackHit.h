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

/// \file GPUTPCGMMergedTrackHit.h
/// \author David Rohr

#ifndef GPUTPCGMMERGEDTRACKHIT_H
#define GPUTPCGMMERGEDTRACKHIT_H

#include "GPUCommonDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUTPCGMMergedTrackHit {
  unsigned int num;
  unsigned char slice, row, leg, state;
#ifdef GPUCA_ALIROOT_LIB
  float x, y, z;
  unsigned short amp;
#endif

  // NOTE: the lower states must match those from ClusterNative!
  enum hitState { flagSplitPad = 0x1,
                  flagSplitTime = 0x2,
                  flagSplit = 0x3,
                  flagEdge = 0x4,
                  flagSingle = 0x8,
                  flagShared = 0x10,
                  clustererAndSharedFlags = 0x1F,
                  flagRejectDistance = 0x20,
                  flagRejectErr = 0x40,
                  flagReject = 0x60,
                  flagNotFit = 0x80 };
};

struct GPUTPCGMMergedTrackHitXYZ {
  float x, y, z;
  unsigned short amp;
#ifdef GPUCA_TPC_RAW_PROPAGATE_PAD_ROW_TIME
  float pad;
  float time;
#endif
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
