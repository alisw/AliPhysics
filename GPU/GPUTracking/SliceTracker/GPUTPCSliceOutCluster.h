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

/// \file GPUTPCSliceOutCluster.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCSLICEOUTCLUSTER_H
#define GPUTPCSLICEOUTCLUSTER_H

#include "GPUTPCDef.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
/**
 * @class GPUTPCSliceOutCluster
 * GPUTPCSliceOutCluster class contains clusters which are assigned to slice tracks.
 * It is used to send the data from TPC slice trackers to the GlobalMerger
 */
class GPUTPCSliceOutCluster
{
 public:
  GPUhd() void Set(unsigned int id, unsigned char row, unsigned char flags, unsigned short amp, float x, float y, float z)
  {
    mRow = row;
    mFlags = flags;
    mAmp = amp;
    mId = id;
    mX = x;
    mY = y;
    mZ = z;
  }

  GPUhd() float GetX() const { return mX; }
  GPUhd() float GetY() const { return mY; }
  GPUhd() float GetZ() const { return mZ; }
  GPUhd() unsigned int GetId() const { return mId; }
  GPUhd() unsigned char GetRow() const { return mRow; }
  GPUhd() unsigned char GetFlags() const { return mFlags; }
  GPUhd() unsigned short GetAmp() const { return mAmp; }

 private:
  unsigned int mId;     // Id ( slice, patch, cluster )
  unsigned char mRow;   // row
  unsigned char mFlags; // flags
  unsigned short mAmp;  // amplitude
  float mX;             // coordinates
  float mY;             // coordinates
  float mZ;             // coordinates

#ifdef GPUCA_TPC_RAW_PROPAGATE_PAD_ROW_TIME
 public:
  float mPad;
  float mTime;
#endif
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
