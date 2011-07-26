// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTPCTrackGeometry.cxx
/// @author Matthias Richter
/// @date   2011-05-20
/// @brief  Desciption of a track by a sequence of track points
///

#include "AliHLTTPCTrackGeometry.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCSpacePointContainer.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTComponent.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "TMath.h"
#include <memory>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCTrackGeometry)

AliHLTTPCTrackGeometry::AliHLTTPCTrackGeometry()
  : AliHLTTrackGeometry()
{
  /// standard constructor
}

AliHLTTPCTrackGeometry::AliHLTTPCTrackGeometry(const AliHLTTPCTrackGeometry& src)
  : AliHLTTrackGeometry(src)
{
  /// copy constructor
}

AliHLTTPCTrackGeometry& AliHLTTPCTrackGeometry::operator=(const AliHLTTPCTrackGeometry& src)
{
  /// assignment operator
  AliHLTTrackGeometry::operator=(src);
  return *this;
}

AliHLTTPCTrackGeometry::~AliHLTTPCTrackGeometry()
{
  /// destructor
}

float AliHLTTPCTrackGeometry::GetPlaneAlpha(AliHLTUInt32_t planeId) const
{
  /// alpha of the plane
  UInt_t slice=AliHLTTPCSpacePointData::GetSlice(planeId);
  float alpha=( slice + 0.5 ) * TMath::Pi() / 9.0;
  if (alpha>TMath::TwoPi()) alpha-=TMath::TwoPi();
  return alpha;
}

float AliHLTTPCTrackGeometry::GetPlaneR(AliHLTUInt32_t planeId) const
{
  /// radial distance from global {0,0,0}
  UInt_t partition=AliHLTTPCSpacePointData::GetPatch(planeId);
  UInt_t number=AliHLTTPCSpacePointData::GetNumber(planeId);
  Int_t row=AliHLTTPCTransform::GetFirstRow(partition)+number;
  return AliHLTTPCTransform::Row2X(row);
}

float AliHLTTPCTrackGeometry::GetPlaneTheta(AliHLTUInt32_t /*planeId*/) const
{
  /// theta of the plane
  return 0.0;
}

bool AliHLTTPCTrackGeometry::CheckBounds(AliHLTUInt32_t planeId, float u, float /*v*/) const
{
  /// check bounds in u and v coordinate
  float r=GetPlaneR(planeId);
  if (r<AliHLTTPCTransform::GetFirstRow(0)) return false;

  // TODO: check if the pad width needs to be considered here
  return TMath::Abs(TMath::ASin(u/r))<=TMath::Pi()/18;
}

int AliHLTTPCTrackGeometry::CalculateTrackPoints(const AliHLTExternalTrackParam& track)
{
  /// calculate the track points, expects the global magnetic field to be initialized
  AliHLTGlobalBarrelTrack bt(track);
  return CalculateTrackPoints(bt);
}

int AliHLTTPCTrackGeometry::CalculateTrackPoints(AliHLTGlobalBarrelTrack& track)
{
  /// calculate the track points, expects the global magnetic field to be initialized
  int iResult=0;
  int firstpadrow=0;
  for (;
       firstpadrow<AliHLTTPCTransform::GetNRows() && 
	 AliHLTTPCTransform::Row2X(firstpadrow)+AliHLTTPCTransform::GetPadLength(firstpadrow)<track.GetX();
       firstpadrow++);
  if (firstpadrow>=AliHLTTPCTransform::GetNRows()) return 0;
  iResult=CalculateTrackPoints(track, firstpadrow, 1);
  if (iResult>=0 && firstpadrow>0)
    iResult=CalculateTrackPoints(track, firstpadrow-1, -1);
  return iResult;
}

int AliHLTTPCTrackGeometry::CalculateTrackPoints(AliHLTGlobalBarrelTrack& track, int firstpadrow, int step)
{
  /// calculate the track points, expects the global magnetic field to be initialized
  float offsetAlpha=0.0;
  for (int padrow=firstpadrow; padrow>=0 && padrow<AliHLTTPCTransform::GetNRows(); padrow+=step) {
    float x=AliHLTTPCTransform::Row2X(padrow);
    float y=0.0;
    float z=0.0;

    int maxshift=9;
    int shift=0;
    int result=0;
    do {
      // start calculation of crossing points with padrow planes in the slice of the first point
      // plane alpha corresponds to alpha of the track, switch to neighboring slice if the result
      // is out of bounds
      if ((result=track.CalculateCrossingPoint(x, track.GetAlpha()-offsetAlpha, y, z))<1) break;
      float pointAlpha=TMath::ATan(y/x);
      if (TMath::Abs(pointAlpha)>TMath::Pi()/18) {
	offsetAlpha+=(pointAlpha>0?-1:1)*TMath::Pi()/9;
      	result=0;
      }
    } while (result==0 && shift++<maxshift);
    if (result<1) continue;
    float planealpha=track.GetAlpha()-offsetAlpha;
    if (planealpha<0) planealpha+=TMath::TwoPi();
    int slice=int(9*planealpha/TMath::Pi());
    //if (z<0) slice+=18;
    int partition=AliHLTTPCTransform::GetPatch(padrow);
    int row=padrow-AliHLTTPCTransform::GetFirstRow(partition);
    UInt_t id=AliHLTTPCSpacePointData::GetID(slice, partition, row);
    if (TMath::Abs(planealpha-GetPlaneAlpha(id))>0.0001) {
      HLTError("alpha missmatch for plane %08x (slice %d): alpha from id %f (%.0f), expected %f (%.0f)", id, slice, GetPlaneAlpha(id), 180*GetPlaneAlpha(id)/TMath::Pi(), planealpha, 180*planealpha/TMath::Pi());
    }
    AddTrackPoint(AliHLTTrackPoint(id, y, z));
  }
  return 0;
}

int AliHLTTPCTrackGeometry::FindMatchingTrackPoint(AliHLTUInt32_t spacepointId, float spacepoint[3], AliHLTUInt32_t& planeId)
{
  /// find the track point which can be associated to a spacepoint with coordinates and id
  UInt_t slice=AliHLTTPCSpacePointData::GetSlice(spacepointId);
  UInt_t partition=AliHLTTPCSpacePointData::GetPatch(spacepointId);
  int row=AliHLTTPCTransform::GetPadRow(spacepoint[0]);
  if (row<AliHLTTPCTransform::GetFirstRow(partition) || row>AliHLTTPCTransform::GetLastRow(partition)) {
    HLTError("row number %d calculated from x value %f is outside slice %d partition %d", row, spacepoint[0], slice, partition);
    return -EINVAL;
  }
  row-=AliHLTTPCTransform::GetFirstRow(partition);
  UInt_t id=AliHLTTPCSpacePointData::GetID(slice, partition, row);
  const AliHLTTrackPoint* point=GetTrackPoint(id);
  if (!point && slice<18) {
    id=AliHLTTPCSpacePointData::GetID(slice+18, partition, row);
    point=GetTrackPoint(id);
  } else if (!point && slice>=18) {
    id=AliHLTTPCSpacePointData::GetID(slice-18, partition, row);
    point=GetTrackPoint(id);
  }
  if (point) {
    planeId=id;
    if (point->HaveAssociatedSpacePoint()) return 0; // already occupied
    float maxdy=5;
    //float maxdz=5;
    if (TMath::Abs(point->GetU()-spacepoint[1])>maxdy) return -ENOENT;
    //if (TMath::Abs(point->GetV()-spacepoint[2])>maxdz) return -ENOENT;
    return 1;
  }
  return -ENOENT;
}

AliHLTSpacePointContainer* AliHLTTPCTrackGeometry::ConvertToSpacePoints() const
{
  /// create a collection of all points
  std::auto_ptr<AliHLTTPCSpacePointContainer> spacepoints(new AliHLTTPCSpacePointContainer);
  if (!spacepoints.get()) return NULL;

  const vector<AliHLTTrackPoint>& trackPoints=GetTrackPoints();
  unsigned i=0;
  while (i<trackPoints.size()) {
    // allocate buffer for all points, even though the buffer might not be filled
    // completely because of a partition change
    int nofPoints=trackPoints.size()-i;
    int blocksize=sizeof(AliHLTTPCClusterData)+nofPoints*sizeof(AliHLTTPCSpacePointData);
    AliHLTUInt8_t* pBuffer=spacepoints->Alloc(blocksize);
    if (!pBuffer) return NULL;
    AliHLTTPCClusterData* pClusterData=reinterpret_cast<AliHLTTPCClusterData*>(pBuffer);
    pClusterData->fSpacePointCnt=0;
    AliHLTTPCSpacePointData* pClusters=pClusterData->fSpacePoints;
    int currentSlice=-1;
    int currentPartition=-1;
    for (; i<trackPoints.size(); i++) {
      AliHLTUInt32_t planeId=trackPoints[i].GetId();
      int slice=AliHLTTPCSpacePointData::GetSlice(planeId);
      int partition=AliHLTTPCSpacePointData::GetPatch(planeId);
      int number=AliHLTTPCSpacePointData::GetNumber(planeId);
      if ((currentSlice>=0 && currentSlice!=slice) || (currentPartition>=0 && currentPartition!=partition)) {
	// change of partition or slice, need to go to next block
	// 2011-07-26 currently all spacepoints go into one block, if separated
	// blocks per partition are needed one has to leave the inner loop here
	// and set the data block specification below
	// Caution: not tested, only the last block seems to make it through
	//break;
      }
      currentSlice=slice;
      currentPartition=partition;
      pClusters[pClusterData->fSpacePointCnt].fX=GetPlaneR(planeId);
      pClusters[pClusterData->fSpacePointCnt].fY=trackPoints[i].GetU();
      pClusters[pClusterData->fSpacePointCnt].fZ=trackPoints[i].GetV();
      pClusters[pClusterData->fSpacePointCnt].fID=planeId;
      pClusters[pClusterData->fSpacePointCnt].fPadRow=AliHLTTPCTransform::GetFirstRow(partition)+number;
      pClusters[pClusterData->fSpacePointCnt].fSigmaY2=0.;
      pClusters[pClusterData->fSpacePointCnt].fSigmaZ2=0.;
      pClusters[pClusterData->fSpacePointCnt].fCharge=0;
      pClusters[pClusterData->fSpacePointCnt].fQMax=0;
      pClusters[pClusterData->fSpacePointCnt].fUsed=0;
      pClusters[pClusterData->fSpacePointCnt].fTrackN=0;
      pClusterData->fSpacePointCnt++;
    }
    AliHLTComponentBlockData bd;
    AliHLTComponent::FillBlockData(bd);
    bd.fPtr=pBuffer;
    bd.fSize=blocksize;
    AliHLTComponent::SetDataType(bd.fDataType, "CLUSTERS", "TPC ");
    bd.fSpecification=//AliHLTTPCDefinitions::EncodeDataSpecification(currentSlice, currentSlice, currentPartition, currentPartition);
    spacepoints->AddInputBlock(&bd);
  }

  return spacepoints.release();
}
