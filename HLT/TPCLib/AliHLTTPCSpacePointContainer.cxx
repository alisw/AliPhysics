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

/// @file   AliHLTTPCSpacePointContainer.h
/// @author Matthias Richter
/// @date   2011-04-29
/// @brief  Helper class for handling of HLT TPC space point data blocks
///

#include "AliHLTTPCSpacePointContainer.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTComponent.h"
#include "AliHLTTemplates.h"
#include "TMath.h"
#include <memory>
#include <algorithm>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCSpacePointContainer)

AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointContainer()
  : AliHLTSpacePointContainer()
  , fClusters()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointContainer(const AliHLTTPCSpacePointContainer& c)
  : AliHLTSpacePointContainer(c)
  , fClusters(c.fClusters.begin(), c.fClusters.end())
{
  /// copy constructor
}

AliHLTTPCSpacePointContainer& AliHLTTPCSpacePointContainer::operator=(const AliHLTTPCSpacePointContainer& c)
{
  /// assignment operator
  if (&c==this) return *this;
  AliHLTSpacePointContainer::operator=(c);
  fClusters=c.fClusters;

  return *this;
}

AliHLTTPCSpacePointContainer::~AliHLTTPCSpacePointContainer()
{
  // destructor
}

int AliHLTTPCSpacePointContainer::AddInputBlock(const AliHLTComponentBlockData* pDesc)
{
  // add input block to the collection
  if (!pDesc) return -EINVAL;
  int count=0;
  if (pDesc->fDataType!=AliHLTTPCDefinitions::fgkClustersDataType) {
    HLTWarning("ignoring data block of type %s", AliHLTComponent::DataType2Text(pDesc->fDataType).c_str());
    return 0;
  }
  if (!pDesc->fPtr || pDesc->fSize<sizeof(AliHLTTPCClusterData)) return -ENODATA;

  // consistency check of the input block
  const AliHLTTPCClusterData* pClusterData=reinterpret_cast<AliHLTTPCClusterData*>(pDesc->fPtr);
  if (pClusterData->fSpacePointCnt*sizeof(AliHLTTPCSpacePointData)+sizeof(AliHLTTPCClusterData)!=pDesc->fSize) {
    HLTError("data block of type %s corrupted: number of entries %d is not consistent with block size %d",
	     AliHLTComponent::DataType2Text(pDesc->fDataType).c_str(), pClusterData->fSpacePointCnt, pDesc->fSize);
    return -EBADMSG;
  }

  AliHLTUInt8_t minslice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
  AliHLTUInt8_t maxslice = AliHLTTPCDefinitions::GetMaxSliceNr( pDesc->fSpecification );
  AliHLTUInt8_t minpart  = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );
  AliHLTUInt8_t maxpart  = AliHLTTPCDefinitions::GetMaxPatchNr( pDesc->fSpecification );
  bool bIsSinglePartition=(pDesc->fSpecification==kAliHLTVoidDataSpec?false:minslice==maxslice && minpart==maxpart);

  for (UInt_t i=0; i<pClusterData->fSpacePointCnt; i++) {
    AliHLTUInt32_t clusterID=~(AliHLTUInt32_t)0;
    if (bIsSinglePartition) {
      // cluster ID has to match ID encoded from slice, partition and index
      clusterID=AliHLTTPCSpacePointData::GetID(minslice, minpart, i);
      if (clusterID!=pClusterData->fSpacePoints[i].fID) {
	HLTWarning("cluster ID 0x%08x does not match slice %d partition %d index %d",
		   pClusterData->fSpacePoints[i].fID, minslice, minpart, i);
      }
    } else {
      // check the cluster ID for correct bounds
      clusterID=pClusterData->fSpacePoints[i].fID;
      UInt_t clusterSlice =AliHLTTPCSpacePointData::GetSlice(clusterID);
      UInt_t clusterPart  =AliHLTTPCSpacePointData::GetPatch(clusterID);
      UInt_t clusterNumber=AliHLTTPCSpacePointData::GetNumber(clusterID);
      if (pDesc->fSpecification!=kAliHLTVoidDataSpec) {
	if (clusterSlice<minslice || clusterSlice>maxslice ||
	    clusterPart<minpart || clusterPart>maxpart) {
	  HLTWarning("cluster ID 0x%08d out of range indicated by data specification 0x%08x", clusterID, pDesc->fSpecification);
	} else if (clusterNumber!=i) {
	  HLTWarning("cluster ID 0x%08d does not match index %d in block 0x%08x", clusterID, i, pDesc->fSpecification);
	}
      }
    }
    if (fClusters.find(clusterID)==fClusters.end()) {
      // new cluster
      fClusters[clusterID]=AliHLTTPCSpacePointProperties(&pClusterData->fSpacePoints[i]);
      count++;
    } else {
      HLTError("cluster with ID 0x%08x already existing, skipping cluster %d of data block 0x%08x",
	       clusterID, i, pDesc->fSpecification);
    }
  }

  return count;
}

int AliHLTTPCSpacePointContainer::GetClusterIDs(vector<AliHLTUInt32_t>& tgt) const
{
  // get array of cluster IDs
  tgt.clear();
  transform(fClusters.begin(), fClusters.end(), back_inserter(tgt), HLT::AliGetKey());
  return tgt.size();
}

float AliHLTTPCSpacePointContainer::GetX(AliHLTUInt32_t clusterID) const
{
  // get X coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->fX;
}

float AliHLTTPCSpacePointContainer::GetXWidth(AliHLTUInt32_t clusterID) const
{
  // get error for X coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return 0.0; // fixed in padrow number
}

float AliHLTTPCSpacePointContainer::GetY(AliHLTUInt32_t clusterID) const
{
  // get Y coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->fY;
}

float AliHLTTPCSpacePointContainer::GetYWidth(AliHLTUInt32_t clusterID) const
{
  // get error for Y coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->fSigmaY2;
}

float AliHLTTPCSpacePointContainer::GetZ(AliHLTUInt32_t clusterID) const
{
  // get Z coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->fZ;
}

float AliHLTTPCSpacePointContainer::GetZWidth(AliHLTUInt32_t clusterID) const
{
  // get error for Z coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->fSigmaZ2;
}

float AliHLTTPCSpacePointContainer::GetCharge(AliHLTUInt32_t clusterID) const
{
  // get charge
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->fCharge;
}

float AliHLTTPCSpacePointContainer::GetPhi(AliHLTUInt32_t clusterID) const
{
  // get charge
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  int slice=AliHLTTPCSpacePointData::GetSlice(fClusters.find(clusterID)->second.Data()->fID);
  return ( slice + 0.5 ) * TMath::Pi() / 9.0;
}

void AliHLTTPCSpacePointContainer::Clear(Option_t * /*option*/)
{
  // clear the object and reset pointer references
}

void AliHLTTPCSpacePointContainer::Print(ostream& out, Option_t */*option*/) const
{
  // print to stream
  out << "AliHLTTPCSpacePointContainer::Print" << endl;
  out << "n clusters: " << fClusters.size() << endl;
  for (std::map<AliHLTUInt32_t, AliHLTTPCSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    out << " " << cl->first << cl->second << endl;
  }
}

AliHLTSpacePointContainer* AliHLTTPCSpacePointContainer::SelectByTrack(int trackId, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a specific track
  std::auto_ptr<AliHLTTPCSpacePointContainer> c(new AliHLTTPCSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties, int>(&AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::GetTrackId,trackId));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCSpacePointContainer::SelectByMC(int mcId, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a specific MC track
  std::auto_ptr<AliHLTTPCSpacePointContainer> c(new AliHLTTPCSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties, int>(&AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::GetMCId,mcId));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCSpacePointContainer::UsedClusters(bool /*bAlloc*/) const
{
  /// create a collection of all used clusters
  std::auto_ptr<AliHLTTPCSpacePointContainer> c(new AliHLTTPCSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties, bool>(&AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::IsUsed,true));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCSpacePointContainer::UnusedClusters(bool /*bAlloc*/) const
{
  /// create a collection of all unused clusters
  std::auto_ptr<AliHLTTPCSpacePointContainer> c(new AliHLTTPCSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties, bool>(&AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::IsUsed,false));
  return c.release();
}

int AliHLTTPCSpacePointContainer::MarkUsed(const AliHLTUInt32_t* clusterIDs, int arraySize)
{
  /// mark the clusters with specified IDs as used
  if (!clusterIDs) return -EINVAL;
  int iCount=0;
  for (int i=0; i<arraySize; i++) {
    if (fClusters.find(clusterIDs[i])==fClusters.end()) continue;
    fClusters[clusterIDs[i]].MarkUsed();
    iCount++;
  }
  return iCount;
}

int AliHLTTPCSpacePointContainer::SetTrackID(int trackID, const AliHLTUInt32_t* clusterIDs, int arraySize)
{
  /// set track id for specified clusters
  if (!clusterIDs) return -EINVAL;
  int iCount=0;
  for (int i=0; i<arraySize; i++) {
    if (fClusters.find(clusterIDs[i])==fClusters.end()) continue;
    fClusters[clusterIDs[i]].SetTrackId(trackID);
    iCount++;
  }
  return iCount;
}

int AliHLTTPCSpacePointContainer::SetMCID(int mcID, const AliHLTUInt32_t* clusterIDs, int arraySize)
{
  /// set mc id for specified clusters
  if (!clusterIDs) return -EINVAL;
  int iCount=0;
  for (int i=0; i<arraySize; i++) {
    if (fClusters.find(clusterIDs[i])==fClusters.end()) continue;
    fClusters[clusterIDs[i]].SetMCId(mcID);
    iCount++;
  }
  return iCount;
}

AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::AliHLTTPCSpacePointProperties()
  : fCluster(NULL)
  , fUsed(false)
  , fTrackId(-1)
  , fMCId(-1)
{
  // constructor
}

AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::AliHLTTPCSpacePointProperties(const AliHLTTPCSpacePointData* pCluster)
  : fCluster(pCluster)
  , fUsed(false)
  , fTrackId(-1)
  , fMCId(-1)
{
  // constructor
}

AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::AliHLTTPCSpacePointProperties(const AliHLTTPCSpacePointProperties& src)
  : fCluster(src.fCluster)
  , fUsed(src.fUsed)
  , fTrackId(src.fTrackId)
  , fMCId(src.fMCId)
{
  // copy constructor
}

AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties& AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::operator=(const AliHLTTPCSpacePointProperties& src)
{
  // assignment operator
  if (&src==this) return *this;
  fCluster=src.fCluster;
  fUsed=src.fUsed;
  fTrackId=src.fTrackId;
  fMCId=src.fMCId;

  return *this;
}

void AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties::Print(ostream& out, Option_t */*option*/) const
{
  // print info
  if (!Data()) {
    out << "no data";
    return;
  }
  const AliHLTTPCSpacePointData* data=Data();
  out << " " << data->fX << " " << data->fY << " " << data->fZ << " " << (UInt_t)data->fPadRow
      << " " << fTrackId << " " << fMCId << " " << fUsed;
}

ostream& operator<<(ostream &out, const AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties& p)
{
  p.Print(out);
  return out;
}
