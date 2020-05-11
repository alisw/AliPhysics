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

/// \file GPUTPCTrack.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCTRACK_H
#define GPUTPCTRACK_H

#include "GPUTPCBaseTrackParam.h"
#include "GPUTPCDef.h"
#include "GPUTPCSliceOutCluster.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
/**
 * @class GPUTPCTrack
 *
 * The class describes the [partially] reconstructed TPC track [candidate].
 * The class is dedicated for internal use by the GPUTPCTracker algorithm.
 * The track parameters at both ends are stored separately in the GPUTPCEndPoint class
 */
MEM_CLASS_PRE()
class GPUTPCTrack
{
 public:
#if !defined(GPUCA_GPUCODE)
  GPUTPCTrack() : mFirstHitID(0), mNHits(0), mLocalTrackId(-1), mParam()
  {
  }
  ~GPUTPCTrack() CON_DEFAULT;
#endif //! GPUCA_GPUCODE

  GPUhd() int NHits() const { return mNHits; }
  GPUhd() int LocalTrackId() const { return mLocalTrackId; }
  GPUhd() int FirstHitID() const { return mFirstHitID; }
  GPUhd() MakeType(const MEM_LG(GPUTPCBaseTrackParam) &) Param() const { return mParam; }

  GPUhd() void SetNHits(int v) { mNHits = v; }
  GPUhd() void SetLocalTrackId(int v) { mLocalTrackId = v; }
  GPUhd() void SetFirstHitID(int v) { mFirstHitID = v; }

  MEM_TEMPLATE()
  GPUhd() void SetParam(const MEM_TYPE(GPUTPCBaseTrackParam) & v) { mParam = v; }

  // Only if used as replacement for SliceOutTrack
  GPUhd() static int GetSize(int nClust) { return sizeof(GPUTPCTrack) + nClust * sizeof(GPUTPCSliceOutCluster); }
  GPUhd() const GPUTPCTrack* GetNextTrack() const { return (const GPUTPCTrack*)(((char*)this) + GetSize(mNHits)); }
  GPUhd() GPUTPCTrack* NextTrack() { return (GPUTPCTrack*)(((char*)this) + GetSize(mNHits)); }
  GPUhd() void SetOutTrackCluster(int i, const GPUTPCSliceOutCluster& v) { ((GPUTPCSliceOutCluster*)((char*)this + sizeof(*this)))[i] = v; }
  GPUhd() const GPUTPCSliceOutCluster* OutTrackClusters() const { return (const GPUTPCSliceOutCluster*)((char*)this + sizeof(*this)); }
  GPUhd() const GPUTPCSliceOutCluster& OutTrackCluster(int i) const { return OutTrackClusters()[i]; }

 private:
  int mFirstHitID;   // index of the first track cell in the track->cell pointer array
  int mNHits;        // number of track cells
  int mLocalTrackId; // Id of local track this global track belongs to, index of this track itself if it is a local track
  MEM_LG(GPUTPCBaseTrackParam)
  mParam; // track parameters

 private:
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif // GPUTPCTRACK_H
