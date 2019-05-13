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

/// \file GPUTPCHitId.h
/// \author Matthias Kretz, Sergey Gorbunov, David Rohr

#ifndef GPUTPCHITID_H
#define GPUTPCHITID_H

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCHitId
{
 public:
  GPUhd() void Set(int row, int hit) { mId = (hit << 8) | row; }
  GPUhd() int RowIndex() const { return mId & 0xff; }
  GPUhd() int HitIndex() const { return mId >> 8; }

 private:
  int mId;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCHITID_H
