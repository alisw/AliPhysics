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
#include "AliHLTTPCGeometry.h"
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
  , fSelections()
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
  , fSelections()
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

  Clear();
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
    // check the cluster ID for correct bounds
    AliHLTUInt32_t clusterID = pClusterData->fSpacePoints[i].GetRawID();
    UInt_t clusterSlice  = AliHLTTPCSpacePointData::GetSlice(clusterID);
    UInt_t clusterPart   = AliHLTTPCSpacePointData::GetPatch(clusterID);
    UInt_t clusterNumber = AliHLTTPCSpacePointData::GetNumber(clusterID);
    if (pDesc->fSpecification!=kAliHLTVoidDataSpec) {
      if (clusterSlice<minslice || clusterSlice>maxslice ||
	  clusterPart<minpart || clusterPart>maxpart) {
	HLTWarning("cluster ID 0x%08d out of range indicated by data specification 0x%08x", clusterID, pDesc->fSpecification);
      } 
    }
    {
      // consistency check for x and row number
      // UInt_t clusterSlice =AliHLTTPCSpacePointData::GetSlice(clusterID);
      // UInt_t clusterPart  =AliHLTTPCSpacePointData::GetPatch(clusterID);
      // int row=AliHLTTPCGeometry::GetPadRow(pClusterData->fSpacePoints[i].fX);
      // if (row<AliHLTTPCGeometry::GetFirstRow(clusterPart) || row>AliHLTTPCGeometry::GetLastRow(clusterPart)) {
      // 	HLTError("row number %d calculated from x value %f is outside slice %d partition %d, expected row %d"
      // 		 , row, pClusterData->fSpacePoints[i].fX, clusterSlice, clusterPart, pClusterData->fSpacePoints[i].fPadRow);
      // }
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

bool AliHLTTPCSpacePointContainer::Check(AliHLTUInt32_t clusterID) const
{
  // check if the cluster is available
  return fClusters.find(clusterID)!=fClusters.end();
}

const vector<AliHLTUInt32_t>* AliHLTTPCSpacePointContainer::GetClusterIDs(AliHLTUInt32_t mask)
{
  // get array of cluster IDs filtered by mask
  if (fSelections.find(mask)!=fSelections.end()) {
    // return existing selection
    return fSelections.find(mask)->second;
  }
  // create new collection
  vector<AliHLTUInt32_t>* selected=new vector<AliHLTUInt32_t>;
  if (!selected) return NULL;
  UInt_t slice=AliHLTTPCSpacePointData::GetSlice(mask);
  UInt_t partition=AliHLTTPCSpacePointData::GetPatch(mask);
  HLTInfo("creating collection 0x%08x", mask);
  for (std::map<AliHLTUInt32_t, AliHLTTPCSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    UInt_t s=AliHLTTPCSpacePointData::GetSlice(cl->first);
    UInt_t p=AliHLTTPCSpacePointData::GetPatch(cl->first);
    if ((slice>=(unsigned)AliHLTTPCGeometry::GetNSlice() || s==slice) && 
	(partition>=(unsigned)AliHLTTPCGeometry::GetNumberOfPatches() || p==partition)) {
      selected->push_back(cl->first);
    }
  }
  HLTInfo("collection 0x%08x with %d spacepoints", mask, selected->size());
  fSelections[mask]=selected;
  return selected;
}

float AliHLTTPCSpacePointContainer::GetX(AliHLTUInt32_t clusterID) const
{
  // get X coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  // FIXME: understand deviation from the nominal x value
  // there is a small deviation in the x coordinate - padrow number correlation
  // in principle, the clusterfinder only uses the mapping to set the x parameter.
  // now extracting the x value from the padrow no.
  //return fClusters.find(clusterID)->second.Data()->fX;
  return AliHLTTPCGeometry::Row2X(fClusters.find(clusterID)->second.Data()->fPadRow);
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

float AliHLTTPCSpacePointContainer::GetMaxSignal(AliHLTUInt32_t clusterID) const
{
  // get charge
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->fQMax;
}

float AliHLTTPCSpacePointContainer::GetPhi(AliHLTUInt32_t clusterID) const
{
  // get charge

  // phi can be derived directly from the id, no need to search
  // for existing cluster
  int slice=AliHLTTPCSpacePointData::GetSlice(clusterID);
  return ( slice + 0.5 ) * TMath::Pi() / 9.0;
}

UChar_t AliHLTTPCSpacePointContainer::GetPadRow(AliHLTUInt32_t clusterID) const
{
  // get pad row
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0;
  return fClusters.find(clusterID)->second.Data()->fPadRow;
}

void AliHLTTPCSpacePointContainer::Clear(Option_t * option)
{
  // clear the object and reset pointer references
  fClusters.clear();

  for (std::map<AliHLTUInt32_t, vector<AliHLTUInt32_t>*>::iterator selection=fSelections.begin();
       selection!=fSelections.end(); selection++) {
    if (selection->second) delete selection->second;
  }
  fSelections.clear();

  AliHLTSpacePointContainer::Clear(option);
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

AliHLTSpacePointContainer* AliHLTTPCSpacePointContainer::SelectByMask(AliHLTUInt32_t mask, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a space point mask
  std::auto_ptr<AliHLTTPCSpacePointContainer> c(new AliHLTTPCSpacePointContainer);
  if (!c.get()) return NULL;

  UInt_t slice=AliHLTTPCSpacePointData::GetSlice(mask);
  UInt_t partition=AliHLTTPCSpacePointData::GetPatch(mask);
  for (std::map<AliHLTUInt32_t, AliHLTTPCSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    UInt_t s=AliHLTTPCSpacePointData::GetSlice(cl->first);
    UInt_t p=AliHLTTPCSpacePointData::GetPatch(cl->first);
    if ((slice>=(unsigned)AliHLTTPCGeometry::GetNSlice() || s==slice) && 
	(partition>=(unsigned)AliHLTTPCGeometry::GetNumberOfPatches() || p==partition)) {
      c->fClusters[cl->first]=cl->second;
    }
  }
  return c.release();
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

int AliHLTTPCSpacePointContainer::GetTrackID(AliHLTUInt32_t clusterID) const
{
  /// get track id for specified cluster
  map<AliHLTUInt32_t, AliHLTTPCSpacePointProperties>::const_iterator element=fClusters.find(clusterID);
  if (element==fClusters.end()) return -1;
  return element->second.GetTrackId();
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
      << " " << data->GetSigmaY2() << " "  << data->GetSigmaZ2()
      << " " << data->GetCharge() << " "  << data->GetQMax()
      << " " << fTrackId << " " << fMCId << " " << fUsed;
}

ostream& operator<<(ostream &out, const AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties& p)
{
  p.Print(out);
  return out;
}
