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

/// \file GPUTPCGMSliceTrack.h
/// \author Sergey Gorbunov, David Rohr

#ifndef GPUTPCGMSLICETRACK_H
#define GPUTPCGMSLICETRACK_H

#include "GPUTPCSliceOutTrack.h"
#include "GPUTPCGMTrackParam.h"
#include "GPUCommonMath.h"
#include "GPUO2DataTypes.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
/**
 * @class GPUTPCGMSliceTrack
 *
 * The class describes TPC slice tracks used in GPUTPCGMMerger
 */
class GPUTPCGMMerger;
class GPUTPCGMSliceTrack
{
 public:
  float Alpha() const { return mAlpha; }
  char Slice() const { return (char)mSlice; }
  char CSide() const { return mSlice >= 18; }
  int NClusters() const { return mNClusters; }
  int PrevNeighbour() const { return mNeighbour[0]; }
  int NextNeighbour() const { return mNeighbour[1]; }
  int Neighbour(int i) const { return mNeighbour[i]; }
  int PrevSegmentNeighbour() const { return mSegmentNeighbour[0]; }
  int NextSegmentNeighbour() const { return mSegmentNeighbour[1]; }
  int SegmentNeighbour(int i) const { return mSegmentNeighbour[i]; }
  const GPUTPCSliceOutTrack* OrigTrack() const { return mOrigTrack; }
  float X() const { return mX; }
  float Y() const { return mY; }
  float Z() const { return mZ; }
  float SinPhi() const { return mSinPhi; }
  float CosPhi() const { return mCosPhi; }
  float SecPhi() const { return mSecPhi; }
  float DzDs() const { return mDzDs; }
  float QPt() const { return mQPt; }
  float TZOffset() const { return mTZOffset; }
  float Leg() const { return mLeg; }

  int LocalTrackId() const { return mLocalTrackId; }
  void SetLocalTrackId(int v) { mLocalTrackId = v; }
  int GlobalTrackId(int n) const { return mGlobalTrackIds[n]; }
  void SetGlobalTrackId(int n, int v) { mGlobalTrackIds[n] = v; }

  float MaxClusterZ() { return CAMath::Max(mOrigTrack->Clusters()->GetZ(), (mOrigTrack->Clusters() + mOrigTrack->NClusters() - 1)->GetZ()); }
  float MinClusterZ() { return CAMath::Min(mOrigTrack->Clusters()->GetZ(), (mOrigTrack->Clusters() + mOrigTrack->NClusters() - 1)->GetZ()); }
  GPUd() float MaxClusterT(const o2::tpc::ClusterNative* cls) { return CAMath::Max(cls[mOrigTrack->Clusters()->GetId()].getTime(), cls[(mOrigTrack->Clusters() + mOrigTrack->NClusters() - 1)->GetId()].getTime()); }
  GPUd() float MinClusterT(const o2::tpc::ClusterNative* cls) { return CAMath::Min(cls[mOrigTrack->Clusters()->GetId()].getTime(), cls[(mOrigTrack->Clusters() + mOrigTrack->NClusters() - 1)->GetId()].getTime()); }

  void Set(const GPUTPCGMTrackParam& trk, const GPUTPCSliceOutTrack* sliceTr, float alpha, int slice);
  void Set(const GPUTPCGMMerger* merger, const GPUTPCSliceOutTrack* sliceTr, float alpha, int slice);

  void SetGlobalSectorTrackCov()
  {
    mC0 = 1;
    mC2 = 1;
    mC3 = 0;
    mC5 = 1;
    mC7 = 0;
    mC9 = 1;
    mC10 = 0;
    mC12 = 0;
    mC14 = 10;
  }

  void SetNClusters(int v) { mNClusters = v; }
  void SetPrevNeighbour(int v) { mNeighbour[0] = v; }
  void SetNextNeighbour(int v) { mNeighbour[1] = v; }
  void SetNeighbor(int v, int i) { mNeighbour[i] = v; }
  void SetPrevSegmentNeighbour(int v) { mSegmentNeighbour[0] = v; }
  void SetNextSegmentNeighbour(int v) { mSegmentNeighbour[1] = v; }
  void SetLeg(unsigned char v) { mLeg = v; }

  void CopyParamFrom(const GPUTPCGMSliceTrack& t)
  {
    mX = t.mX;
    mY = t.mY;
    mZ = t.mZ;
    mSinPhi = t.mSinPhi;
    mDzDs = t.mDzDs;
    mQPt = t.mQPt;
    mCosPhi = t.mCosPhi, mSecPhi = t.mSecPhi;
    mAlpha = t.mAlpha;
  }

  bool FilterErrors(const GPUTPCGMMerger* merger, int iSlice, float maxSinPhi = GPUCA_MAX_SIN_PHI, float sinPhiMargin = 0.f);
  bool TransportToX(GPUTPCGMMerger* merger, float x, float Bz, GPUTPCGMBorderTrack& b, float maxSinPhi, bool doCov = true) const;
  bool TransportToXAlpha(GPUTPCGMMerger* merger, float x, float sinAlpha, float cosAlpha, float Bz, GPUTPCGMBorderTrack& b, float maxSinPhi) const;
  void CopyBaseTrackCov();

 private:
  const GPUTPCSliceOutTrack* mOrigTrack;                    // pointer to original slice track
  float mX, mY, mZ, mSinPhi, mDzDs, mQPt, mCosPhi, mSecPhi; // parameters
  float mTZOffset;                                          // Z offset with early transform, T offset otherwise
  float mC0, mC2, mC3, mC5, mC7, mC9, mC10, mC12, mC14;     // covariances
  float mAlpha;                                             // alpha angle
  int mSlice;                                               // slice of this track segment
  int mNClusters;                                           // N clusters
  int mNeighbour[2];                                        //
  int mSegmentNeighbour[2];                                 //
  int mLocalTrackId;                                        // Corrected local track id in terms of GMSliceTracks array
  int mGlobalTrackIds[2];                                   // IDs of associated global tracks
  unsigned char mLeg;                                       // Leg of this track segment
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
