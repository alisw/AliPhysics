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

/// @file   AliHLTTPCRawSpacePointContainer.h
/// @author Matthias Richter
/// @date   2011-04-29
/// @brief  Helper class for handling of HLT TPC space point data blocks
///

#include "AliHLTTPCRawSpacePointContainer.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTComponent.h"
#include "AliHLTTemplates.h"
#include "TMath.h"
#include <memory>
#include <algorithm>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCRawSpacePointContainer)

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointContainer()
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

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointContainer(const AliHLTTPCRawSpacePointContainer& c)
  : AliHLTSpacePointContainer(c)
  , fClusters(c.fClusters.begin(), c.fClusters.end())
  , fSelections()
{
  /// copy constructor
}

AliHLTTPCRawSpacePointContainer& AliHLTTPCRawSpacePointContainer::operator=(const AliHLTTPCRawSpacePointContainer& c)
{
  /// assignment operator
  if (&c==this) return *this;
  AliHLTSpacePointContainer::operator=(c);
  fClusters=c.fClusters;

  return *this;
}

AliHLTTPCRawSpacePointContainer::~AliHLTTPCRawSpacePointContainer()
{
  // destructor
  Clear();
}

int AliHLTTPCRawSpacePointContainer::AddInputBlock(const AliHLTComponentBlockData* pDesc)
{
  // add input block to the collection
  if (!pDesc) return -EINVAL;
  int count=0;
  if (pDesc->fDataType!=AliHLTTPCDefinitions::fgkRawClustersDataType) {
    HLTWarning("ignoring data block of type %s", AliHLTComponent::DataType2Text(pDesc->fDataType).c_str());
    return 0;
  }
  if (!pDesc->fPtr || pDesc->fSize<sizeof(AliHLTTPCRawClusterData)) return -ENODATA;

  // consistency check of the input block
  const AliHLTTPCRawClusterData* pClusterData=reinterpret_cast<AliHLTTPCRawClusterData*>(pDesc->fPtr);
  if (pClusterData->fCount*sizeof(AliHLTTPCRawCluster)+sizeof(AliHLTTPCRawClusterData)!=pDesc->fSize) {
    HLTError("data block of type %s corrupted: number of entries %d is not consistent with block size %d",
	     AliHLTComponent::DataType2Text(pDesc->fDataType).c_str(), pClusterData->fCount, pDesc->fSize);
    return -EBADMSG;
  }

  AliHLTUInt8_t minslice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
  //AliHLTUInt8_t maxslice = AliHLTTPCDefinitions::GetMaxSliceNr( pDesc->fSpecification );
  AliHLTUInt8_t minpart  = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );
  //AliHLTUInt8_t maxpart  = AliHLTTPCDefinitions::GetMaxPatchNr( pDesc->fSpecification );
  //bool bIsSinglePartition=(pDesc->fSpecification==kAliHLTVoidDataSpec?false:minslice==maxslice && minpart==maxpart);

  for (UInt_t i=0; i<pClusterData->fCount; i++) {
    AliHLTUInt32_t clusterID=~(AliHLTUInt32_t)0;
    // cluster ID from slice, partition and index
    clusterID=AliHLTTPCSpacePointData::GetID(minslice, minpart, i);

    if (fClusters.find(clusterID)==fClusters.end()) {
      // new cluster
      fClusters[clusterID]=AliHLTTPCRawSpacePointProperties(&pClusterData->fClusters[i]);
      count++;
    } else {
      HLTError("cluster with ID 0x%08x already existing, skipping cluster %d of data block 0x%08x",
	       clusterID, i, pDesc->fSpecification);
    }
  }

  return count;
}

int AliHLTTPCRawSpacePointContainer::GetClusterIDs(vector<AliHLTUInt32_t>& tgt) const
{
  // get array of cluster IDs
  tgt.clear();
  transform(fClusters.begin(), fClusters.end(), back_inserter(tgt), HLT::AliGetKey());
  return tgt.size();
}

bool AliHLTTPCRawSpacePointContainer::Check(AliHLTUInt32_t clusterID) const
{
  // check if the cluster is available
  return fClusters.find(clusterID)!=fClusters.end();
}

const vector<AliHLTUInt32_t>* AliHLTTPCRawSpacePointContainer::GetClusterIDs(AliHLTUInt32_t mask)
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
  for (std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    UInt_t s=AliHLTTPCSpacePointData::GetSlice(cl->first);
    UInt_t p=AliHLTTPCSpacePointData::GetPatch(cl->first);
    if ((slice>=(unsigned)AliHLTTPCTransform::GetNSlice() || s==slice) && 
	(partition>=(unsigned)AliHLTTPCTransform::GetNumberOfPatches() || p==partition)) {
      selected->push_back(cl->first);
    }
  }
  HLTInfo("collection 0x%08x with %d spacepoints", mask, selected->size());
  fSelections[mask]=selected;
  return selected;
}

float AliHLTTPCRawSpacePointContainer::GetX(AliHLTUInt32_t clusterID) const
{
  // get X coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  // FIXME: understand deviation from the nominal x value
  // there is a small deviation in the x coordinate - padrow number correlation
  // in principle, the clusterfinder only uses the mapping to set the x parameter.
  // now extracting the x value from the padrow no.
  //return fClusters.find(clusterID)->second.Data()->fX;
  return AliHLTTPCTransform::Row2X(fClusters.find(clusterID)->second.Data()->fPadRow);
}

float AliHLTTPCRawSpacePointContainer::GetXWidth(AliHLTUInt32_t clusterID) const
{
  // get error for X coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return 0.0; // fixed in padrow number
}

float AliHLTTPCRawSpacePointContainer::GetY(AliHLTUInt32_t clusterID) const
{
  // get Y coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->GetPad();
}

float AliHLTTPCRawSpacePointContainer::GetYWidth(AliHLTUInt32_t clusterID) const
{
  // get error for Y coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->GetSigmaY2();
}

float AliHLTTPCRawSpacePointContainer::GetZ(AliHLTUInt32_t clusterID) const
{
  // get Z coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->GetTime();
}

float AliHLTTPCRawSpacePointContainer::GetZWidth(AliHLTUInt32_t clusterID) const
{
  // get error for Z coordinate
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->GetSigmaZ2();
}

float AliHLTTPCRawSpacePointContainer::GetCharge(AliHLTUInt32_t clusterID) const
{
  // get charge
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->GetCharge();
}

float AliHLTTPCRawSpacePointContainer::GetMaxSignal(AliHLTUInt32_t clusterID) const
{
  // get charge
  if (fClusters.find(clusterID)==fClusters.end() ||
      fClusters.find(clusterID)->second.Data()==NULL) return 0.0;
  return fClusters.find(clusterID)->second.Data()->GetQMax();
}

float AliHLTTPCRawSpacePointContainer::GetPhi(AliHLTUInt32_t clusterID) const
{
  // get charge

  // phi can be derived directly from the id, no need to search
  // for existing cluster
  int slice=AliHLTTPCSpacePointData::GetSlice(clusterID);
  return ( slice + 0.5 ) * TMath::Pi() / 9.0;
}

void AliHLTTPCRawSpacePointContainer::Clear(Option_t * option)
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

void AliHLTTPCRawSpacePointContainer::Print(ostream& out, Option_t */*option*/) const
{
  // print to stream
  out << "AliHLTTPCRawSpacePointContainer::Print" << endl;
  out << "n clusters: " << fClusters.size() << endl;
  for (std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    out << " " << cl->first << cl->second << endl;
  }
}

AliHLTSpacePointContainer* AliHLTTPCRawSpacePointContainer::SelectByMask(AliHLTUInt32_t mask, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a space point mask
  std::auto_ptr<AliHLTTPCRawSpacePointContainer> c(new AliHLTTPCRawSpacePointContainer);
  if (!c.get()) return NULL;

  UInt_t slice=AliHLTTPCSpacePointData::GetSlice(mask);
  UInt_t partition=AliHLTTPCSpacePointData::GetPatch(mask);
  for (std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    UInt_t s=AliHLTTPCSpacePointData::GetSlice(cl->first);
    UInt_t p=AliHLTTPCSpacePointData::GetPatch(cl->first);
    if ((slice>=(unsigned)AliHLTTPCTransform::GetNSlice() || s==slice) && 
	(partition>=(unsigned)AliHLTTPCTransform::GetNumberOfPatches() || p==partition)) {
      c->fClusters[cl->first]=cl->second;
    }
  }
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCRawSpacePointContainer::SelectByTrack(int trackId, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a specific track
  std::auto_ptr<AliHLTTPCRawSpacePointContainer> c(new AliHLTTPCRawSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties, int>(&AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::GetTrackId,trackId));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCRawSpacePointContainer::SelectByMC(int mcId, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a specific MC track
  std::auto_ptr<AliHLTTPCRawSpacePointContainer> c(new AliHLTTPCRawSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties, int>(&AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::GetMCId,mcId));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCRawSpacePointContainer::UsedClusters(bool /*bAlloc*/) const
{
  /// create a collection of all used clusters
  std::auto_ptr<AliHLTTPCRawSpacePointContainer> c(new AliHLTTPCRawSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties, bool>(&AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::IsUsed,true));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCRawSpacePointContainer::UnusedClusters(bool /*bAlloc*/) const
{
  /// create a collection of all unused clusters
  std::auto_ptr<AliHLTTPCRawSpacePointContainer> c(new AliHLTTPCRawSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties, bool>(&AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::IsUsed,false));
  return c.release();
}

int AliHLTTPCRawSpacePointContainer::MarkUsed(const AliHLTUInt32_t* clusterIDs, int arraySize)
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

int AliHLTTPCRawSpacePointContainer::SetTrackID(int trackID, const AliHLTUInt32_t* clusterIDs, int arraySize)
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

int AliHLTTPCRawSpacePointContainer::GetTrackID(AliHLTUInt32_t clusterID) const
{
  /// get track id for specified cluster
  map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator element=fClusters.find(clusterID);
  if (element==fClusters.end()) return -1;
  return element->second.GetTrackId();
}

int AliHLTTPCRawSpacePointContainer::SetMCID(int mcID, const AliHLTUInt32_t* clusterIDs, int arraySize)
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

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::AliHLTTPCRawSpacePointProperties()
  : fCluster(NULL)
  , fUsed(false)
  , fTrackId(-1)
  , fMCId(-1)
{
  // constructor
}

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::AliHLTTPCRawSpacePointProperties(const AliHLTTPCRawCluster* pCluster)
  : fCluster(pCluster)
  , fUsed(false)
  , fTrackId(-1)
  , fMCId(-1)
{
  // constructor
}

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::AliHLTTPCRawSpacePointProperties(const AliHLTTPCRawSpacePointProperties& src)
  : fCluster(src.fCluster)
  , fUsed(src.fUsed)
  , fTrackId(src.fTrackId)
  , fMCId(src.fMCId)
{
  // copy constructor
}

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties& AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::operator=(const AliHLTTPCRawSpacePointProperties& src)
{
  // assignment operator
  if (&src==this) return *this;
  fCluster=src.fCluster;
  fUsed=src.fUsed;
  fTrackId=src.fTrackId;
  fMCId=src.fMCId;

  return *this;
}

void AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::Print(ostream& out, Option_t */*option*/) const
{
  // print info
  if (!Data()) {
    out << "no data";
    return;
  }
  const AliHLTTPCRawCluster* data=Data();
  out << " " << data->GetPadRow() << " " << data->GetPad() << " " << data->GetTime()
      << " " << data->GetSigmaY2() << " "  << data->GetSigmaZ2()
      << " " << data->GetCharge() << " "  << data->GetQMax()
      << " " << fTrackId << " " << fMCId << " " << fUsed;
}

ostream& operator<<(ostream &out, const AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties& p)
{
  p.Print(out);
  return out;
}
