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

/// \file GPUTPCHit.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCHIT_H
#define GPUTPCHIT_H

#include "GPUTPCDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
/**
 * @class GPUTPCHit
 *
 * The GPUTPCHit class is the internal representation
 * of the TPC clusters for the GPUTPCTracker algorithm.
 *
 */
class GPUTPCHit
{
 public:
  GPUhd() float Y() const { return mY; }
  GPUhd() float Z() const { return mZ; }

  GPUhd() void SetY(float v) { mY = v; }
  GPUhd() void SetZ(float v) { mZ = v; }

 protected:
  float mY, mZ; // Y and Z position of the TPC cluster

 private:
  friend class GPUTPCNeighboursFinder;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCHIT_H
