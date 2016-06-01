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

/// @file   AliHLTTPCRawSpacePointContainer.cxx
/// @author Matthias Richter, Sergey Gorbunov
/// @date   2011-08-08
/// @brief  Helper class for handling of HLT TPC raw cluster data blocks

#include "AliHLTTPCRawSpacePointContainer.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCClusterFlagsData.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTComponent.h"
#include "AliHLTTemplates.h"
#include "AliHLTDataDeflater.h"
#include "AliLog.h"
#include "TMath.h"
#include <memory>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCRawSpacePointContainer)

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointContainer(int mode, int createFlags)
  : AliHLTSpacePointContainer()
  , fClusters()
  , fSelections()
  , fBlocks()
  , fSingleBlock()
  , fMode(mode)
  , fCreateFlags(createFlags)
  , fWrittenClusterIds(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (fMode&kModeSingle) {
    fSingleBlock.SetDecoder(NULL);
    fSingleBlock.SetGrid(AllocateIndexGrid());
  }
}

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointContainer(const AliHLTTPCRawSpacePointContainer& c)
  : AliHLTSpacePointContainer(c)
  , fClusters(c.fClusters.begin(), c.fClusters.end())
  , fSelections()
  , fBlocks()
  , fSingleBlock()
  , fMode(c.fMode)
  , fCreateFlags(c.fCreateFlags)
  , fWrittenClusterIds(NULL)
{
  /// copy constructor
}

AliHLTTPCRawSpacePointContainer& AliHLTTPCRawSpacePointContainer::operator=(const AliHLTTPCRawSpacePointContainer& c)
{
  /// assignment operator
  if (&c==this) return *this;
  AliHLTSpacePointContainer::operator=(c);
  fClusters=c.fClusters;
  fMode=c.fMode;
  fCreateFlags = c.fCreateFlags;
  fWrittenClusterIds=NULL;

  return *this;
}

AliHLTTPCRawSpacePointContainer::~AliHLTTPCRawSpacePointContainer()
{
  // destructor
  Clear();
  if (fSingleBlock.GetGrid()) delete fSingleBlock.GetGrid();
  if (fWrittenClusterIds) delete fWrittenClusterIds;
}

int AliHLTTPCRawSpacePointContainer::AddInputBlock(const AliHLTComponentBlockData* pDesc)
{
  // add input block to the collection
  if (!pDesc) return -EINVAL;
  int iResult=0;
  int count=0;
  if (pDesc->fDataType!=AliHLTTPCDefinitions::fgkRawClustersDataType && pDesc->fDataType!=AliHLTTPCDefinitions::fgkRawClustersDataTypeNotCompressed) {
    HLTWarning("ignoring data block of type %s", AliHLTComponent::DataType2Text(pDesc->fDataType).c_str());
    return 0;
  }
  if (!pDesc->fPtr) return -ENODATA;
  if (pDesc->fSize<sizeof(AliHLTTPCRawClusterData)) return 0;

  AliHLTTPCRawClusterData* rawClusters = (AliHLTTPCRawClusterData*)(pDesc->fPtr);  

  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
  AliHLTUInt8_t part  = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );

  AliHLTUInt32_t decoderIndex=AliHLTTPCGeometry::CreateClusterID(slice, part, 0);

  AliHLTSpacePointPropertyGrid* pGrid=NULL;
  if (fMode&kModeSingle) {
    pGrid=fSingleBlock.GetGrid();
  } else {
    if (fBlocks.find(decoderIndex)!=fBlocks.end()) {
      HLTError("data block of slice %d partition %d already added, skipping data block", slice, part);
      return -EEXIST;
    }
  }

  if (fMode&kModeSingle && !pGrid) {
    pGrid=AllocateIndexGrid();
    if (!pGrid) {
      return -ENOMEM;
    }
  }

  if (fMode&kModeCreateMap) { // register immediately

    //UInt_t nofClusters=rawClusters->fCount;

    for( UInt_t icl=0; icl<rawClusters->fCount; icl++){      
      const AliHLTTPCRawCluster &cl = rawClusters->fClusters[icl];
      AliHLTUInt32_t clusterID=~(AliHLTUInt32_t)0;
      // cluster ID from slice, partition and index
      clusterID=AliHLTTPCGeometry::CreateClusterID(slice, part, icl);

      if (fClusters.find(clusterID)==fClusters.end()) {
	// new cluster
	fClusters[clusterID]=AliHLTTPCRawSpacePointProperties(&cl);
	count++;
      } else {
	HLTError("cluster with ID 0x%08x already existing, skipping cluster %d of data block 0x%08x",
		 clusterID, icl, pDesc->fSpecification);
      }
    }
  }

  if (pGrid && (iResult=PopulateAccessGrid(pGrid, rawClusters, slice, part))<0) {
    HLTError("failed to populate access grid for block %s 0x%09x: %d",
	     AliHLTComponent::DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification, iResult);
    return iResult;
  }
  
  if (fMode&kModeSingle) {
    fSingleBlock.SetDecoder(rawClusters);
    fSingleBlock.SetGrid(pGrid);
    fSingleBlock.SetId(decoderIndex);
  } else {
    fBlocks[decoderIndex]=AliHLTTPCRawSpacePointBlock(decoderIndex, rawClusters, pGrid);
  }
  return count;
}

AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* AliHLTTPCRawSpacePointContainer::AllocateIndexGrid()
{
  // allocate index grid, one single point to define the dimensions
  
  // max 33 padrows, step 1 padrow
  // max 140 pads, step 2x max delta pad
  // max 1024 time bins, step 2x max delta time
  return new AliHLTSpacePointPropertyGrid(33, 1.0,
					  140, 2*AliHLTTPCDefinitions::GetMaxClusterDeltaPad(),
					  1024, 2*AliHLTTPCDefinitions::GetMaxClusterDeltaTime()
					  );
}

int AliHLTTPCRawSpacePointContainer::PopulateAccessGrid(AliHLTSpacePointPropertyGrid* pGrid, AliHLTUInt32_t mask) const
{
  // populate an access grid
  if (!pGrid) return -EINVAL;

  pGrid->Clear();
  
  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(mask);
  AliHLTUInt8_t partition  = AliHLTTPCDefinitions::GetMinPatchNr(mask);
  AliHLTUInt32_t decoderIndex=AliHLTTPCGeometry::CreateClusterID(slice, partition, 0);
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointBlock>::const_iterator block=fBlocks.find(decoderIndex);
  if (block==fBlocks.end()) {
    HLTError("can not find data block of id 0x%08x", mask);
    return -ENOENT;
  }
  return PopulateAccessGrid(pGrid, block->second.GetDecoder(), slice, partition);
}

int AliHLTTPCRawSpacePointContainer::PopulateAccessGrid(AliHLTSpacePointPropertyGrid* pGrid, AliHLTTPCRawClusterData* pDecoder,
							 int slice, int partition) const
{
  // populate an access grid
  if (!pDecoder) return -EINVAL;
  int iResult=0;

  for( UInt_t icl=0; icl<pDecoder->fCount; icl++){
    const AliHLTTPCRawCluster &cl = pDecoder->fClusters[icl];
    iResult=pGrid->CountSpacePoint(cl.GetPadRow(), cl.GetPad(), cl.GetTime());
    if (iResult<0)
      HLTError("CountSpacePoint %d %f %f failed: %d", cl.GetPadRow(), cl.GetPad(), cl.GetTime(), iResult);
  }
  
  for( UInt_t icl=0; icl<pDecoder->fCount; icl++ ){
    const AliHLTTPCRawCluster &cl = pDecoder->fClusters[icl];
    AliHLTUInt32_t id=AliHLTTPCGeometry::CreateClusterID(slice, partition, icl);
    iResult=pGrid->AddSpacePoint(AliHLTSpacePointProperties(id), cl.GetPadRow(), cl.GetPad(), cl.GetTime());
    if (iResult<0)
      HLTError("AddSpacePoint 0x%08x %d %f %f failed: %d", id, cl.GetPadRow(), cl.GetPad(), cl.GetTime(), iResult);
  }

  return 0;
}

const AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* AliHLTTPCRawSpacePointContainer::GetSpacePointPropertyGrid(AliHLTUInt32_t mask) const
{
  // get the access grid for a data block
  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(mask);
  AliHLTUInt8_t part  = AliHLTTPCDefinitions::GetMinPatchNr(mask);
  AliHLTUInt32_t decoderIndex=AliHLTTPCGeometry::CreateClusterID(slice, part, 0);
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointBlock>::const_iterator block=fBlocks.find(decoderIndex);
  if (block==fBlocks.end()) {
    HLTError("can not find data block of id 0x%08x", mask);
    return NULL;
  }
  return block->second.GetGrid();
}

int AliHLTTPCRawSpacePointContainer::SetSpacePointPropertyGrid(AliHLTUInt32_t mask, AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pGrid)
{
  // set the access grid for a data block
  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(mask);
  AliHLTUInt8_t part  = AliHLTTPCDefinitions::GetMinPatchNr(mask);
  AliHLTUInt32_t decoderIndex=AliHLTTPCGeometry::CreateClusterID(slice, part, 0);
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointBlock>::iterator block=fBlocks.find(decoderIndex);
  if (block==fBlocks.end()) {
    HLTError("can not find data block of id 0x%08x", mask);
    return -ENOENT;
  }
  if (block->second.GetGrid()!=NULL && pGrid!=NULL && block->second.GetGrid()!=pGrid) {
    // there is trouble ahead because this will delete the index grid instance
    // but it might be an external pointer supposed to be deleted by someone else
    ALIHLTERRORGUARD(1, "overriding previous instance of index grid, potential memory leak or invalid deallocation ahead");
  }
  block->second.SetGrid(pGrid);
  return 0;
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
  UInt_t slice=AliHLTTPCGeometry::CluID2Slice(mask);
  UInt_t partition=AliHLTTPCGeometry::CluID2Partition(mask);
  //HLTInfo("creating collection 0x%08x", mask);

  // the first cluster with number 0 has equal ID to mask unless
  // the mask selects multiple partitions/slices
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(mask);
  bool bAll=false;
  if (slice>=(unsigned)AliHLTTPCGeometry::GetNSlice() ||
      partition>=(unsigned)AliHLTTPCGeometry::GetNumberOfPatches()) {
    cl=fClusters.begin();
    bAll=true;
  }
  for (; cl!=fClusters.end(); cl++) {
    UInt_t s=AliHLTTPCGeometry::CluID2Slice(cl->first);
    UInt_t p=AliHLTTPCGeometry::CluID2Partition(cl->first);
    if ((slice>=(unsigned)AliHLTTPCGeometry::GetNSlice() || s==slice) && 
	(partition>=(unsigned)AliHLTTPCGeometry::GetNumberOfPatches() || p==partition)) {
      selected->push_back(cl->first);
    } else if (!bAll) {
      // no need to continue, we are out of the range
      break;
    }
  }
  //HLTInfo("collection 0x%08x with %d spacepoints", mask, selected->size());
  fSelections[mask]=selected;
  return selected;
}

float AliHLTTPCRawSpacePointContainer::GetX(AliHLTUInt32_t clusterID) const
{
  // get X coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.GetCluster()==NULL) return 0.0;
  // FIXME: understand deviation from the nominal x value
  // there is a small deviation in the x coordinate - padrow number correlation
  // in principle, the clusterfinder only uses the mapping to set the x parameter.
  // now extracting the x value from the padrow no.
  //return cl->second.Decoder()->fX;
  return cl->second.GetCluster()->GetPadRow();
}

float AliHLTTPCRawSpacePointContainer::GetXWidth(AliHLTUInt32_t clusterID) const
{
  // get error for X coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.GetCluster()==NULL) return 0.0;
  return 0.0; // fixed in padrow number
}

float AliHLTTPCRawSpacePointContainer::GetY(AliHLTUInt32_t clusterID) const
{
  // get Y coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.GetCluster()==NULL) return 0.0;
  return cl->second.GetCluster()->GetPad();
}

float AliHLTTPCRawSpacePointContainer::GetYWidth(AliHLTUInt32_t clusterID) const
{
  // get error for Y coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.GetCluster()==NULL) return 0.0;  
  return cl->second.GetCluster()->GetSigmaPad2();
}

float AliHLTTPCRawSpacePointContainer::GetZ(AliHLTUInt32_t clusterID) const
{
  // get Z coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.GetCluster()==NULL) return 0.0;
  return cl->second.GetCluster()->GetTime();
}

float AliHLTTPCRawSpacePointContainer::GetZWidth(AliHLTUInt32_t clusterID) const
{
  // get error for Z coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.GetCluster()==NULL) return 0.0;
  return cl->second.GetCluster()->GetSigmaTime2();
}

float AliHLTTPCRawSpacePointContainer::GetCharge(AliHLTUInt32_t clusterID) const
{
  // get charge
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.GetCluster()==NULL) return 0.0;
  return cl->second.GetCluster()->GetCharge();
}

float AliHLTTPCRawSpacePointContainer::GetQMax(AliHLTUInt32_t clusterID) const
{
  // get charge
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.GetCluster()==NULL) return 0.0;
  return cl->second.GetCluster()->GetQMax();
}

float AliHLTTPCRawSpacePointContainer::GetPhi(AliHLTUInt32_t clusterID) const
{
  // get charge

  // phi can be derived directly from the id, no need to search
  // for existing cluster
  int slice=AliHLTTPCGeometry::CluID2Slice(clusterID);
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

  for (std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointBlock>::iterator block=fBlocks.begin();
       block!=fBlocks.end(); block++) {
    if (block->second.GetGrid()) delete block->second.GetGrid();
  }
  fBlocks.clear();

  fSingleBlock.SetDecoder(NULL);
  if (fSingleBlock.GetGrid()) fSingleBlock.GetGrid()->Clear();

  AliHLTSpacePointContainer::Clear(option);
}

void AliHLTTPCRawSpacePointContainer::Print(ostream& out, Option_t */*option*/) const
{
  // print to stream
  std::stringstream str;
  str << "AliHLTTPCRawSpacePointContainer::Print" << endl;
  str << "n clusters: " << fClusters.size() << endl;
  for (std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    str << " 0x" << hex << setw(8) << setfill('0') << cl->first << dec << cl->second << endl;
  }
  out << str.rdbuf();
}

AliHLTSpacePointContainer* AliHLTTPCRawSpacePointContainer::SelectByMask(AliHLTUInt32_t mask, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a space point mask
  std::auto_ptr<AliHLTTPCRawSpacePointContainer> c(new AliHLTTPCRawSpacePointContainer);
  if (!c.get()) return NULL;

  UInt_t slice=AliHLTTPCGeometry::CluID2Slice(mask);
  UInt_t partition=AliHLTTPCGeometry::CluID2Partition(mask);
  for (std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    UInt_t s=AliHLTTPCGeometry::CluID2Slice(cl->first);
    UInt_t p=AliHLTTPCGeometry::CluID2Partition(cl->first);
    if ((slice>=(unsigned)AliHLTTPCGeometry::GetNSlice() || s==slice) && 
	(partition>=(unsigned)AliHLTTPCGeometry::GetNumberOfPatches() || p==partition)) {
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

int AliHLTTPCRawSpacePointContainer::Write(AliHLTUInt8_t* outputPtr,
					    AliHLTUInt32_t size,
					    AliHLTComponentBlockDataList&
					    outputBlocks,
					    AliHLTDataDeflater* pDeflater,
					    const char* option) const
{
  /// write blocks to HLT component output
  AliHLTUInt32_t offset=0;
  if (outputBlocks.size()>0) {
    offset=outputBlocks.back().fOffset+outputBlocks.back().fSize;
  }
  return Write(outputPtr, size, offset, outputBlocks, pDeflater, option);
}

int AliHLTTPCRawSpacePointContainer::Write(AliHLTUInt8_t* outputPtr,
					    AliHLTUInt32_t size,
					    AliHLTUInt32_t offset,
					    AliHLTComponentBlockDataList&
					    outputBlocks,
					    AliHLTDataDeflater* pDeflater,
					    const char* option) const
{
  /// write blocks to HLT component output
  return WriteSorted(outputPtr, size, offset, outputBlocks, pDeflater, option);
}

int AliHLTTPCRawSpacePointContainer::WriteSorted(AliHLTUInt8_t* outputPtr,
						  AliHLTUInt32_t size,
						  AliHLTUInt32_t offset,
						  AliHLTComponentBlockDataList&
						  outputBlocks,
						  AliHLTDataDeflater* pDeflater,
						  const char* option) const
{
  /// write blocks to HLT component output
  int iResult=0;

  if (fMode&kModeSingle) {
    iResult=WriteSorted(outputPtr, size, offset, fSingleBlock.GetDecoder(), fSingleBlock.GetGrid(), fSingleBlock.GetId(), outputBlocks, pDeflater, option);
  } else {
    iResult=-ENOENT;
    for (std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointBlock>::const_iterator block=fBlocks.begin();
	 block!=fBlocks.end(); block++) {
      AliHLTTPCRawClusterData* pDecoder=block->second.GetDecoder();
      AliHLTSpacePointPropertyGrid* pGrid=block->second.GetGrid();
      AliHLTUInt32_t mask=block->first;
      // FIXME: have to propagate the parameter which block is currently to be written
      // for now the index grid is only set for that one
      if (!pGrid) continue;
      iResult=WriteSorted(outputPtr, size, offset, pDecoder, pGrid, mask, outputBlocks, pDeflater, option);
      break; // only one is supposed to be written
    }
    if (iResult==-ENOENT) {
      HLTError("could not find the index grid of the partition to be written");
    }
  }
  return iResult;
}

int AliHLTTPCRawSpacePointContainer::WriteSorted(AliHLTUInt8_t* outputPtr,
						  AliHLTUInt32_t size,
						  AliHLTUInt32_t offset,
						  AliHLTTPCRawClusterData* pDecoder,
						  AliHLTSpacePointPropertyGrid* pGrid,
						  AliHLTUInt32_t mask,
						  AliHLTComponentBlockDataList&
						  outputBlocks,
						  AliHLTDataDeflater* pDeflater,
						  const char* option) const
{
  /// write blocks to HLT component output

  if (!outputPtr || !pDecoder || !pGrid) return -EINVAL;
  if (pDecoder->fCount==0) return 0;
  if (option) {
    // this is only for sending mc labels in simulation and testing
    // no impact to real running
    if (!fWrittenClusterIds && strcmp(option, "write-cluster-ids")==0) {
      const_cast<AliHLTTPCRawSpacePointContainer*>(this)->fWrittenClusterIds=new vector<AliHLTUInt32_t>;
    }
  }
  int iResult=0;
  const AliHLTUInt32_t capacity=size;
  size=0;

  int slice=AliHLTTPCGeometry::CluID2Slice(mask);
  int part=AliHLTTPCGeometry::CluID2Partition(mask);

  // Note: the offset parameter is only for the block descriptors, output pointer and size
  // consider already the offset

  if( fWrittenClusterIds ) fWrittenClusterIds->clear();

  if( capacity < sizeof(AliHLTTPCRawClusterData) ) return -ENOSPC;
 
  AliHLTTPCRawClusterData* blockout=reinterpret_cast<AliHLTTPCRawClusterData*>(outputPtr);
  blockout->fVersion=0;
  blockout->fCount=0;
 
  size = sizeof(AliHLTTPCRawClusterData);

  if (pDeflater) {
    pDeflater->Clear();
    pDeflater->InitBitDataOutput(reinterpret_cast<AliHLTUInt8_t*>(blockout->fClusters), capacity-size);
    blockout->fVersion=pDeflater->GetDeflaterVersion();
    if (fMode&kModeDifferentialPadTime) blockout->fVersion+=2;
  }

  unsigned lastPadRow=0;
  AliHLTUInt64_t lastPad64=0;
  AliHLTUInt64_t lastTime64=0;
  pDeflater->StartEncoder();
  AliHLTSpacePointPropertyGrid::iterator clusterID=pGrid->begin();

  for (; clusterID!=pGrid->end(); clusterID++) {
    if (clusterID.Data().fTrackId>-1) {
      // this is an assigned cluster, skip
      // TODO: introduce selectors into AliHLTIndexGrid::begin to loop
      // consistently over entries, e.g. this check has to be done also
      // in the forwarding of MC labels in
      // AliHLTTPCDataCompressionComponent::ForwardMCLabels
      continue;
    }
    if ((unsigned)slice!=AliHLTTPCGeometry::CluID2Slice(clusterID.Data().fId) ||
	(unsigned)part!=AliHLTTPCGeometry::CluID2Partition(clusterID.Data().fId)) {
      HLTError("cluster id 0x%08x out of slice %d partition %d", clusterID.Data().fId, slice, part);
      iResult = -EBADMSG;
      break;
    }
    int index=AliHLTTPCGeometry::CluID2Index(clusterID.Data().fId);

    if( index<0 || index>= (int)pDecoder->fCount ){
      HLTError("cluster index %d id 0x%08x out of range %d", index, clusterID.Data().fId, (int)pDecoder->fCount);
      iResult = -EBADMSG;
      break;
      continue;
    }

    if( size + sizeof(AliHLTTPCRawCluster) > capacity ){
      iResult = -ENOSPC;
      break;
    }

    const AliHLTTPCRawCluster &input = pDecoder->fClusters[index];

    //const AliHLTTPCRawData::iterator& input=pDecoder->find(index);
    //if (!(input!=pDecoder->end())) continue;
    if (fWrittenClusterIds) fWrittenClusterIds->push_back(clusterID.Data().fId);

    int padrow=input.GetPadRow();
    float pad =input.GetPad();
    float time =input.GetTime();
    float sigmaY2=input.GetSigmaPad2();
    float sigmaZ2=input.GetSigmaTime2();

    if (padrow<0 || pad<0 || time<0 || sigmaY2<0 || sigmaZ2<0 ) {
      // something wrong here, padrow is stored in the cluster header
      // word which has bit pattern 0x3 in bits bit 30 and 31 which was
      // not recognized
      HLTError("wrong cluster data: padrow %d pad %f time %f sigmaY2 %f sigmaZ2 %f", padrow, pad, time, sigmaY2, sigmaZ2 );
      iResult = -EBADMSG;
      break;
    }

    if (!pDeflater) {
      AliHLTTPCRawCluster& c=blockout->fClusters[blockout->fCount];
      padrow+=AliHLTTPCGeometry::GetFirstRow(part);
      c.SetPadRow(padrow);
      c.SetCharge(input.GetCharge());
      c.SetPad(pad);  
      c.SetTime(time);
      c.SetSigmaPad2(sigmaY2);
      c.SetSigmaTime2(sigmaZ2);
      c.SetQMax(input.GetQMax());
      size += sizeof(AliHLTTPCRawCluster);
    } else {
      AliHLTUInt32_t oldSize = pDeflater->GetBitDataOutputSizeBytes(); 
      
      AliHLTUInt64_t padrow64=input.GetPadRow();
      if (padrow64==lastPadRow) {
	padrow64-=lastPadRow;
      } else if (padrow64>lastPadRow) {
	padrow64-=lastPadRow;
	lastPadRow+=padrow64;
      } else {
	HLTError("padrows are not ordered, skip data");
	iResult = -EBADMSG;
	break;
      }
      
      AliHLTUInt32_t padType=0;
      AliHLTUInt32_t signdPad=0;
      AliHLTUInt64_t pad64=0;
      if ((fMode&kModeDifferentialPadTime)!=0 && sigmaY2<.00001) {
	// single pad cluster
	// use twice the pad position to take account for the 0.5 offset
	// added in the AliHLTTPCRawData decoder in accordance with the
	// offline definition. Using the factor 2, this offset is not
	// cut off by rounding
	pad64=(AliHLTUInt64_t)round(2*pad);
	padType=1;
      } else {
	if (!isnan(pad)) pad64=(AliHLTUInt64_t)round(pad*AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kPad].fScale);
      }
      if (fMode&kModeDifferentialPadTime && padType==0) {
	AliHLTUInt64_t dpad64=0;
	if (pad64<lastPad64) {
	  dpad64=lastPad64-pad64;
	  signdPad=1;
	} else {
	  dpad64=pad64-lastPad64;
	  signdPad=0;
	}
	lastPad64=pad64;
	pad64=dpad64;
      }
      AliHLTUInt32_t signdTime=0;
      AliHLTUInt64_t time64=0;
      if (!isnan(time)) time64=(AliHLTUInt64_t)round(time*AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kTime].fScale);
      if (fMode&kModeDifferentialPadTime) {
	AliHLTUInt64_t dtime64=0;
	if (time64<lastTime64) {
	  dtime64=lastTime64-time64;
	  signdTime=1;
	} else {
	  dtime64=time64-lastTime64;
	  signdTime=0;
	}
	lastTime64=time64;
	time64=dtime64;
      }
      AliHLTUInt64_t sigmaY264=0;
      if (!isnan(sigmaY2)) sigmaY264=(AliHLTUInt64_t)round(sigmaY2*AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kSigmaY2].fScale);
      // we can safely use the upper limit as this is an unphysical cluster, no impact to physics
      if (sigmaY264 >= (unsigned)1<<AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kSigmaY2].fBitLength) {
        sigmaY264 = (1<<AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kSigmaY2].fBitLength)-1;
      }
      AliHLTUInt64_t sigmaZ264=0;
      if (!isnan(sigmaZ2)) sigmaZ264=(AliHLTUInt64_t)round(sigmaZ2*AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kSigmaZ2].fScale);
      // we can safely use the upper limit as this is an unphysical cluster, no impact to physics
      if (sigmaZ264 >= (unsigned)1<<AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kSigmaZ2].fBitLength) {
        sigmaZ264 = (1<<AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kSigmaZ2].fBitLength)-1;
      }
      pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kPadRow , padrow64);
      pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kPad    , pad64);
      if (fMode&kModeDifferentialPadTime) pDeflater->OutputBit(padType);
      if (fMode&kModeDifferentialPadTime && padType==0) pDeflater->OutputBit(signdPad);
      pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kTime   , time64);
      if (fMode&kModeDifferentialPadTime) pDeflater->OutputBit(signdTime);
      if (padType==0) pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kSigmaY2, sigmaY264);
      pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kSigmaZ2, sigmaZ264);
      pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kCharge , input.GetCharge());
      pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kQMax   , input.GetQMax());
      
      size += pDeflater->GetBitDataOutputSizeBytes() - oldSize;
    }
    blockout->fCount++;
  }
  pDeflater->StopEncoder();
  AliHLTComponent_BlockData bd;
  AliHLTComponent::FillBlockData(bd);
  bd.fOffset        = offset;
  if (!pDeflater) {
    bd.fSize        = sizeof(AliHLTTPCRawClusterData)+blockout->fCount*sizeof(AliHLTTPCRawCluster);
    bd.fDataType    = AliHLTTPCDefinitions::fgkRawClustersDataTypeNotCompressed;
  } else {
    pDeflater->Pad8Bits();
    bd.fSize        = sizeof(AliHLTTPCRawClusterData)+pDeflater->GetBitDataOutputSizeBytes();
    pDeflater->CloseBitDataOutput();
    bd.fDataType    = AliHLTTPCDefinitions::RemainingClustersCompressedDataType();
  }
  bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(slice, slice, part, part);

  if (iResult>=0 && bd.fSize+size<=capacity) {
    outputBlocks.push_back(bd);
    size = bd.fSize;
  } else {
    iResult=-ENOSPC;
  }
  
  if (iResult >= 0 && fCreateFlags)
  {
    if (size + sizeof(AliHLTTPCClusterFlagsData) > capacity)
    {
      iResult = -ENOSPC;
    }
    else
    {
      AliHLTComponent::FillBlockData(bd);
      bd.fOffset = size + offset;
      bd.fDataType = AliHLTTPCDefinitions::ClustersFlagsDataType();
      bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(slice, slice, part, part);
 
      AliHLTTPCClusterFlagsData* clusterFlagsData = (AliHLTTPCClusterFlagsData*) (outputPtr + size);
      clusterFlagsData->fVersion = 1;
      clusterFlagsData->fNumberOfFlags = 3;

      unsigned int tmpVal = 0;
      unsigned int tmppos = 0;
      unsigned int nEntries = 0;
      unsigned int nClusters = 0;
      for (AliHLTSpacePointPropertyGrid::iterator clusterID=pGrid->begin(); clusterID!=pGrid->end(); clusterID++)
      {
        if (clusterID.Data().fTrackId>-1) continue;  //Cannot handle clusters associated to tracks yet

        if (size + sizeof(AliHLTTPCClusterFlagsData) + (nEntries + 2) * sizeof(tmpVal) > capacity)
        {
          iResult = -ENOSPC;
          break;
        }

        int index=AliHLTTPCGeometry::CluID2Index(clusterID.Data().fId);
        const AliHLTTPCRawCluster &input = pDecoder->fClusters[index];

        unsigned int tmpFlags = input.fFlags & ((1 << clusterFlagsData->fNumberOfFlags) - 1);
        tmpVal |= tmpFlags << tmppos;
        tmppos += clusterFlagsData->fNumberOfFlags;
        nClusters++;
        if (tmppos >= sizeof(tmpVal) * 8)
        {
          tmppos -= sizeof(tmpVal) * 8;
          ((unsigned int*) clusterFlagsData->fData)[nEntries++] = tmpVal;
          tmpVal = 0;
          if (tmppos) tmpVal |= tmpFlags >> (clusterFlagsData->fNumberOfFlags - tmppos);
        }
      }
      if (iResult >= 0)
      {
        if (tmppos) ((unsigned int*) clusterFlagsData->fData)[nEntries++] = tmpVal;
        clusterFlagsData->fNumberOfClusters = nClusters;
        bd.fSize = sizeof(AliHLTTPCClusterFlagsData) + nEntries * sizeof(tmpVal);
        outputBlocks.push_back(bd);
        size += bd.fSize;
      }
    }
  }

  if ( (iResult>=0) && fWrittenClusterIds && (fWrittenClusterIds->size()>0) ) {
    AliHLTComponent::FillBlockData(bd);
    bd.fOffset        = size+offset;
    bd.fSize        = fWrittenClusterIds->size()*sizeof(vector<AliHLTUInt32_t>::value_type);
    if (bd.fSize+size<=capacity) {
      memcpy(outputPtr+size, &(*fWrittenClusterIds)[0], bd.fSize);
      bd.fDataType    = AliHLTTPCDefinitions::RemainingClusterIdsDataType();
      bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(slice, slice, part, part);
      outputBlocks.push_back(bd);    
      size += bd.fSize;
    } else {
      iResult=-ENOSPC;
    }
  }
  
  if (fWrittenClusterIds ) fWrittenClusterIds->clear();
 
  if (iResult<0) return iResult;
  return size;
}

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::AliHLTTPCRawSpacePointProperties()
  : fpCluster(NULL)
  , fUsed(false)
  , fTrackId(-1)
  , fMCId(-1)
{
  // constructor
}

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::AliHLTTPCRawSpacePointProperties(const AliHLTTPCRawCluster* pCluster)
  : fpCluster(pCluster)
  , fUsed(false)
  , fTrackId(-1)
  , fMCId(-1)
{
  // constructor
}

AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::AliHLTTPCRawSpacePointProperties(const AliHLTTPCRawSpacePointProperties& src)
  : fpCluster(src.fpCluster)
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
  fpCluster=src.fpCluster;
  fUsed=src.fUsed;
  fTrackId=src.fTrackId;
  fMCId=src.fMCId;

  return *this;
}

void AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties::Print(ostream& out, Option_t */*option*/) const
{
  // print info
  if (!fpCluster) {
    out << "no data";
    return;
  }
  std::stringstream str;

  str.setf(ios::fixed,ios::floatfield);
  str << " " << setfill(' ') << setw(3) << fpCluster->GetPadRow() 
      << " " << setw(8) << setprecision(3) << fpCluster->GetPad()
      << " " << setw(8) << setprecision(3) << fpCluster->GetTime()
      << " " << setw(8) << setprecision(1) << fpCluster->GetSigmaPad2() 
      << " " << setw(9) << setprecision(1) << fpCluster->GetSigmaTime2()
      << " " << setw(5) << fpCluster->GetCharge() 
      << " " << setw(5) << fpCluster->GetQMax()
      << " " << fTrackId << " " << fMCId << " " << fUsed;
  out << str.rdbuf();
}

ostream& operator<<(ostream &out, const AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties& p)
{
  p.Print(out);
  return out;
}
