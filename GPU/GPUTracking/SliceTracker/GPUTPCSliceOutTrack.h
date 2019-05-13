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

/// \file GPUTPCSliceOutTrack.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCSLICEOUTTRACK_H
#define GPUTPCSLICEOUTTRACK_H

#include "GPUTPCBaseTrackParam.h"
#include "GPUTPCSliceOutCluster.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
/**
 * @class GPUTPCSliceOutTrack
 * GPUTPCSliceOutTrack class is used to store TPC tracks,
 * which are reconstructed by the TPCCATracker slice tracker.
 *
 * The class contains:
 * - fitted track parameters at its first row, the covariance matrix, chi^2, NDF (number of degrees of freedom )
 * - n of clusters assigned to the track
 * - clusters in corresponding cluster arrays
 *
 * The class is used to transport the data between GPUTPCTracker{Component} and GPUTPCGBMerger{Component}
 *
 */
class GPUTPCSliceOutTrack
{
 public:
  GPUhd() int NClusters() const { return mNClusters; }
  GPUhd() const GPUTPCBaseTrackParam& Param() const { return mParam; }
  GPUhd() const GPUTPCSliceOutCluster& Cluster(int i) const { return mClusters[i]; }
  GPUhd() const GPUTPCSliceOutCluster* Clusters() const { return mClusters; }

  GPUhd() void SetNClusters(int v) { mNClusters = v; }
  GPUhd() void SetParam(const GPUTPCBaseTrackParam& v) { mParam = v; }
  GPUhd() void SetCluster(int i, const GPUTPCSliceOutCluster& v) { mClusters[i] = v; }

  GPUhd() static int GetSize(int nClust) { return sizeof(GPUTPCSliceOutTrack) + nClust * sizeof(GPUTPCSliceOutCluster); }

  GPUhd() int LocalTrackId() const { return mLocalTrackId; }
  GPUhd() void SetLocalTrackId(int v) { mLocalTrackId = v; }

  GPUhd() GPUTPCSliceOutTrack* NextTrack()
  {
    return (GPUTPCSliceOutTrack*)(((char*)this) + GetSize(mNClusters));
  }

  GPUhd() const GPUTPCSliceOutTrack* GetNextTrack() const { return (GPUTPCSliceOutTrack*)(((char*)this) + GetSize(mNClusters)); }

 private:
  GPUTPCBaseTrackParam mParam; //* fitted track parameters at its innermost cluster
  int mNClusters;              //* number of track clusters
  int mLocalTrackId;           // See AliHLTPCCATrack.h
#ifdef __OPENCL__
  GPUTPCSliceOutCluster mClusters[1]; //* track clusters
#else
  GPUTPCSliceOutCluster mClusters[0]; //* track clusters
#endif
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
