// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timur Pocheptsov <Timur.Pocheptsov@cern.ch>           *
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

/// @file  AliHLTGlobalTrackResidualsComponent.cxx 
/// @author Timur Pocheptsov
/// @date   
/// @brief A histogramming component for plotting the Y and Z track residual
///        

#if __GNUC__>= 3
using namespace std;
#endif

#include <algorithm>

#include <TMath.h>

#include "AliHLTGlobalTrackResidualsComponent.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTDataTypes.h"

ClassImp(AliHLTGlobalTrackResidualsComponent)

//_______________________________________________________________________________________________

AliHLTGlobalTrackResidualsComponent::AliHLTGlobalTrackResidualsComponent()
                                : AliHLTProcessor(),
                                  fResY("y_residuals", "y residuals", kNBins, -1., 1.),
                                  fResZ("z_residuals", "z residuals", kNBins, -1., 1.),
                                  fSortedX()
{
  //Ctor.
  fResY.SetMarkerStyle(8);
  fResY.SetMarkerSize(0.4);
  fResY.SetXTitle("Y [cm]");
  fResY.SetDirectory(0);
  
  fResZ.SetMarkerStyle(8);
  fResZ.SetMarkerSize(0.4);
  fResZ.SetXTitle("Z [cm]");
  fResZ.SetDirectory(0);
  
  CleanClusters();
}

//_______________________________________________________________________________________________
const char * AliHLTGlobalTrackResidualsComponent::GetComponentID()
{
  //Component's name.
  return "GlobalTrackResiduals";
}

//_______________________________________________________________________________________________
void AliHLTGlobalTrackResidualsComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  //Possible input data types
  list.clear();
  list.push_back(AliHLTTPCDefinitions::fgkClustersDataType|kAliHLTDataOriginTPC);
  list.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
}

//_______________________________________________________________________________________________
AliHLTComponentDataType AliHLTGlobalTrackResidualsComponent::GetOutputDataType()
{
  //Output's data type(s).
  return kAliHLTMultipleDataType;
}

//_______________________________________________________________________________________________
int AliHLTGlobalTrackResidualsComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  //Output's data type(s).
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeTNtuple|kAliHLTDataOriginTPC);
  tgtList.push_back(kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
  return tgtList.size();
}

//_______________________________________________________________________________________________
void AliHLTGlobalTrackResidualsComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //Approximate output histograms sizes.
  constBase  = sizeof(TH1F) * 2; //Size of histogram objects.
  constBase += sizeof(Float_t) * kNBins * 2; //Size of memory, allocated by 1D histograms.
  //Forget about strings (name and title), just multiply by 2.
  constBase *= 2;
  inputMultiplier = 1;
}

//_______________________________________________________________________________________________
AliHLTComponent* AliHLTGlobalTrackResidualsComponent::Spawn()
{
  //Create the component.
  return new AliHLTGlobalTrackResidualsComponent;
}

//_______________________________________________________________________________________________
int AliHLTGlobalTrackResidualsComponent::DoInit(int /*argc*/, const char** /*argv*/)
{
  //(Re)Initialize component.
  ResetHistograms();
  return 0;
}

//_______________________________________________________________________________________________
int AliHLTGlobalTrackResidualsComponent::DoDeinit()
{
  //DoNothing will be better name.
  return 0;
}

//_______________________________________________________________________________________________
int AliHLTGlobalTrackResidualsComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  //Process global barrel tracks and clusters, calculate residuals.
  if (GetFirstInputBlock(kAliHLTDataTypeSOR) || GetFirstInputBlock(kAliHLTDataTypeEOR))
    return 0;

  //Read input data, find residuals, fill histgrams.
  ProcessBlocks();

  //Do output now.
  PushBack(&fResY, kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC, 0);
  PushBack(&fResZ, kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC, 0);

  return 0;
}

//_______________________________________________________________________________________________
void AliHLTGlobalTrackResidualsComponent::ProcessBlocks()
{
  //1. Read cluster blocks.
  ReadClusterBlocks();
  //2. Loop over merged tracks, calculate residuals, fill histogramms.
  std::vector<AliHLTGlobalBarrelTrack> ts;

  Int_t totalTracks = 0;
  const AliHLTComponentBlockData * i = GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);

  for (; i; i = GetNextInputBlock()) {
    if (i->fDataType != (kAliHLTDataTypeTrack | kAliHLTDataOriginTPC))
      continue;

    ts.clear();
    AliHLTGlobalBarrelTrack::ConvertTrackDataArray((AliHLTTracksData*)i->fPtr, i->fSize, ts);

    totalTracks += Int_t(ts.size());

    std::vector<AliHLTGlobalBarrelTrack>::size_type j = 0, e = ts.size();
    for (; j != e; ++j) {
      fSortedX.clear();
      SortHitsX(ts[j]);
      FillResiduals(ts[j]);
    }

    HLTDebug("TrackResiduals found %d tracks", totalTracks);
  }
}

//_______________________________________________________________________________________________
void AliHLTGlobalTrackResidualsComponent::ReadClusterBlocks()
{
  //Loop over blocks, find cluster blocks, extract space points.
  CleanClusters();

  Int_t totalSpacePoints = 0;
  const AliHLTComponentBlockData* iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType);

  for (; iter; iter = GetNextInputBlock()) {
    if (iter->fDataType != AliHLTTPCDefinitions::fgkClustersDataType)
      continue;

    const AliHLTUInt8_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter);
    const AliHLTUInt8_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);

    const AliHLTTPCClusterData * clusterData = (AliHLTTPCClusterData *)iter->fPtr;
    const Int_t nSpacepoint = Int_t(clusterData->fSpacePointCnt);
    totalSpacePoints += nSpacepoint;

    //This part is from AliHLTTPCTrackHistoComponent. Logic is not clear -
    //is it possible that I can have two blocks with same minSlice and minPartition???
    //and one of them with 0 spacepoint?
    if (nSpacepoint) {
      HLTDebug("TrackResiduals component found %d spacepoints in slice %d partition %d", nSpacepoint, minSlice, minPartition);
      fClustersArray[minSlice][minPartition] = (AliHLTTPCSpacePointData*)clusterData->fSpacePoints;
      fNSpacePoints[minSlice][minPartition]  = nSpacepoint;
    }
  }

  HLTDebug("TrackResiduals found %d spacepoints", totalSpacePoints);
}

namespace {

void Rotate(Float_t* xy, Float_t alpha);
Bool_t CmpX(const std::pair<Float_t, UInt_t>& rhs, const std::pair<Float_t, UInt_t>& lhs);

}

//_______________________________________________________________________________________________
void AliHLTGlobalTrackResidualsComponent::SortHitsX(const AliHLTGlobalBarrelTrack& gt)
{
  //Extract hits' Xs for the track gt, sort them.
  fSortedX.clear();

  const UInt_t * hitnum = gt.GetPoints();
  Int_t prevSlice = -1;
  Float_t rotAngle = 0.f;

  for (UInt_t i = 0; i < gt.GetNumberOfPoints(); ++i) {
    const UInt_t idTrack    = hitnum[i];
    const UInt_t pos        = idTrack & 0x3fffff;
    const Int_t  sliceTrack = (idTrack >> 25) & 0x7f;
    const UInt_t patchTrack = (idTrack >> 22) & 0x7;

    if (!fClustersArray[sliceTrack][patchTrack])
      continue;

    //The following conditional is from the original code.
    if (sliceTrack > 36 || patchTrack > 5) {
      HLTError("Corrupted TPC cluster Id: slice %d, patch %d, cluster %d", sliceTrack, patchTrack, idTrack);
      continue;
    }

    if (fNSpacePoints[sliceTrack][patchTrack] <= pos) {
      HLTError("Space point array out of boundaries!");
      continue;
    }

    if (sliceTrack != prevSlice) {
      if (prevSlice != -1)
        prevSlice < sliceTrack ? rotAngle += 0.349066 : rotAngle -= 0.349066;
      prevSlice = sliceTrack;
    }

    Float_t clusterXY[] = {fClustersArray[sliceTrack][patchTrack][pos].fX,
                           fClustersArray[sliceTrack][patchTrack][pos].fY};

    Rotate(clusterXY, rotAngle);

    fSortedX.push_back(std::pair<Float_t, UInt_t>(clusterXY[0], i));
  }

  std::sort(fSortedX.begin(), fSortedX.end(), CmpX);
}

//_______________________________________________________________________________________________
void AliHLTGlobalTrackResidualsComponent::FillResiduals(const AliHLTGlobalBarrelTrack& gt)
{
  //Find residuals using clusters and helix approximation.
  const UInt_t * hitnum = gt.GetPoints();
  AliExternalTrackParam track(gt);
  Int_t prevSlice = -1;
  Float_t rotAngle = 0.f;

  std::vector<std::pair<Float_t, UInt_t> >::size_type i = 0;
  for (; i < fSortedX.size(); ++i) {
    const UInt_t idTrack = hitnum[fSortedX[i].second];
    const UInt_t pos = idTrack & 0x3fffff;
    const Int_t sliceTrack = (idTrack >> 25) & 0x7f;
    const UInt_t patchTrack = (idTrack >> 22) & 0x7;

    if(!fClustersArray[sliceTrack][patchTrack])
      continue;

    //The following conditionals are from the original code.
    if (sliceTrack > 36 || patchTrack > 5) {
      HLTError("Corrupted TPC cluster Id: slice %d, patch %d, cluster %d", sliceTrack, patchTrack, idTrack);
      continue;
    }

    if (fNSpacePoints[sliceTrack][patchTrack] <= pos) {
      HLTError("Space point array out of boundaries!");
      continue;
    }

    if (sliceTrack != prevSlice) {
      if (prevSlice != -1)
        prevSlice < sliceTrack ? rotAngle += 0.349066 : rotAngle -= 0.349066;
      prevSlice = sliceTrack;
    }

    if (track.PropagateTo(fSortedX[i].first, GetBz())) {
      Float_t clusterXYZ[] = {fClustersArray[sliceTrack][patchTrack][pos].fX,
                              fClustersArray[sliceTrack][patchTrack][pos].fY,
                              fClustersArray[sliceTrack][patchTrack][pos].fZ};
      Rotate(clusterXYZ, rotAngle);

      fResY.Fill(clusterXYZ[1] - track.GetY());
      fResZ.Fill(clusterXYZ[2] - track.GetZ());
    } else
      break;
  }
}

//_______________________________________________________________________________________________
void AliHLTGlobalTrackResidualsComponent::CleanClusters()
{
  //Set pointers and counters to zero.
  for (int i = 0; i < 36; ++i) {
    for (int j = 0; j < 6; ++j) {
      fClustersArray[i][j] = 0;
      fNSpacePoints[i][j]  = 0;
    }
  }
}

//_______________________________________________________________________________________________
void AliHLTGlobalTrackResidualsComponent::ResetHistograms()
{
  //Set default values.
  fResY.Reset();
  fResZ.Reset();

  fResY.SetBins(kNBins, -1., 1.);
  fResZ.SetBins(kNBins, -1., 1.);
}

namespace
{

//_______________________________________________________________________________________________
void Rotate(Float_t* xy, Float_t alpha)
{
  //From hits's local to track's _local_.
  const Float_t cosA = TMath::Cos(alpha);
  const Float_t sinA = TMath::Sin(alpha);
  const Float_t xPrim = xy[0] * cosA - xy[1] * sinA;
  const Float_t yPrim = xy[0] * sinA + xy[1] * cosA;
  xy[0] = xPrim;
  xy[1] = yPrim;
}

//_______________________________________________________________________________________________
Bool_t CmpX(const std::pair<Float_t, UInt_t>& rhs, const std::pair<Float_t, UInt_t>& lhs)
{
  //Sort "hits" along x.
  return rhs.first < lhs.first;
}


}
