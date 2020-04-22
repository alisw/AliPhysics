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

/// \file GPUTPCRow.h
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#ifndef GPUTPCROW_H
#define GPUTPCROW_H

#include "GPUTPCDef.h"
#include "GPUTPCGrid.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
/**
 * @class GPUTPCRow
 *
 * The GPUTPCRow class is a hit and cells container for one TPC row.
 * It is the internal class of the GPUTPCTracker algorithm.
 *
 */
MEM_CLASS_PRE()
class GPUTPCRow
{
  MEM_CLASS_PRE2()
  friend class GPUTPCSliceData;

 public:
#if !defined(GPUCA_GPUCODE)
  GPUTPCRow();
#endif //! GPUCA_GPUCODE

  GPUhd() int NHits() const
  {
    return mNHits;
  }
  GPUhd() float X() const { return mX; }
  GPUhd() float MaxY() const { return mMaxY; }
  GPUhd() MakeType(const MEM_LG(GPUTPCGrid) &) Grid() const { return mGrid; }

  GPUhd() float Hy0() const { return mHy0; }
  GPUhd() float Hz0() const { return mHz0; }
  GPUhd() float HstepY() const { return mHstepY; }
  GPUhd() float HstepZ() const { return mHstepZ; }
  GPUhd() float HstepYi() const { return mHstepYi; }
  GPUhd() float HstepZi() const { return mHstepZi; }
  GPUhd() int HitNumberOffset() const { return mHitNumberOffset; }
  GPUhd() unsigned int FirstHitInBinOffset() const { return mFirstHitInBinOffset; }

 private:
  friend class GPUTPCNeighboursFinder;
  friend class GPUTPCStartHitsFinder;

  int mNHits;       // number of hits
  float mX;         // X coordinate of the row
  float mMaxY;      // maximal Y coordinate of the row
  MEM_LG(GPUTPCGrid)
  mGrid; // grid of hits

  // hit packing:
  float mHy0;     // offset
  float mHz0;     // offset
  float mHstepY;  // step size
  float mHstepZ;  // step size
  float mHstepYi; // inverse step size
  float mHstepZi; // inverse step size

  int mHitNumberOffset; // index of the first hit in the hit array, used as
  // offset in GPUTPCSliceData::LinkUp/DownData/HitDataY/...
  unsigned int mFirstHitInBinOffset; // offset in Tracker::mRowData to find the FirstHitInBin
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCROW_H
