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

/// @file   AliHLTTPCHWCFSpacePointContainer.cxx
/// @author Matthias Richter
/// @date   2011-08-08
/// @brief  Helper class for handling of HLT TPC cluster data blocks from the
///         HW ClusterFinder
/// @note   Class is a duplicate of AliHLTTPCHWCFSpacePointContainer and should
///         be merged with it in a generic way

#include "AliHLTTPCHWCFSpacePointContainer.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTComponent.h"
#include "AliHLTTemplates.h"
#include "AliHLTDataDeflater.h"
#include "AliRawDataHeader.h"
#include "AliLog.h"
#include "TMath.h"
#include <memory>
#include <algorithm>
#include <cmath>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCHWCFSpacePointContainer)

AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointContainer()
  : AliHLTSpacePointContainer()
  , fClusters()
  , fSelections()
  , fDecoders()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointContainer(const AliHLTTPCHWCFSpacePointContainer& c)
  : AliHLTSpacePointContainer(c)
  , fClusters(c.fClusters.begin(), c.fClusters.end())
  , fSelections()
  , fDecoders()
{
  /// copy constructor
}

AliHLTTPCHWCFSpacePointContainer& AliHLTTPCHWCFSpacePointContainer::operator=(const AliHLTTPCHWCFSpacePointContainer& c)
{
  /// assignment operator
  if (&c==this) return *this;
  AliHLTSpacePointContainer::operator=(c);
  fClusters=c.fClusters;

  return *this;
}

AliHLTTPCHWCFSpacePointContainer::~AliHLTTPCHWCFSpacePointContainer()
{
  // destructor
  Clear();
}

int AliHLTTPCHWCFSpacePointContainer::AddInputBlock(const AliHLTComponentBlockData* pDesc)
{
  // add input block to the collection
  if (!pDesc) return -EINVAL;
  int count=0;
  if (pDesc->fDataType!=AliHLTTPCDefinitions::fgkHWClustersDataType) {
    HLTWarning("ignoring data block of type %s", AliHLTComponent::DataType2Text(pDesc->fDataType).c_str());
    return 0;
  }
  if (!pDesc->fPtr) return -ENODATA;
  if (pDesc->fSize<=sizeof(AliRawDataHeader)) return 0;

  AliHLTUInt32_t *buffer=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);  
  // skip the first 8 32-bit CDH words
  buffer += 8;
  UInt_t bufferSize32 = ((Int_t)pDesc->fSize - sizeof(AliRawDataHeader) )/sizeof(AliHLTUInt32_t);

  std::auto_ptr<AliHLTTPCHWCFData> pDecoder(new AliHLTTPCHWCFData);
  if (!pDecoder.get()) return -ENOMEM;

  if (pDecoder->Init(reinterpret_cast<AliHLTUInt8_t*>(buffer), bufferSize32*sizeof(AliHLTUInt32_t))<0 ||
      (pDecoder->CheckVersion()<0 && (int)bufferSize32>pDecoder->GetRCUTrailerSize())) {
    HLTError("data block of type %s corrupted: can not decode format",
	     AliHLTComponent::DataType2Text(pDesc->fDataType).c_str());
    return -EBADMSG;
  }

  UInt_t nofClusters=pDecoder->GetNumberOfClusters();

  AliHLTUInt8_t minslice = AliHLTTPCDefinitions::GetMinSliceNr( pDesc->fSpecification );
  AliHLTUInt8_t minpart  = AliHLTTPCDefinitions::GetMinPatchNr( pDesc->fSpecification );

  for (UInt_t i=0; i<nofClusters; i++) {
    AliHLTUInt32_t clusterID=~(AliHLTUInt32_t)0;
    // cluster ID from slice, partition and index
    clusterID=AliHLTTPCSpacePointData::GetID(minslice, minpart, i);

    if (fClusters.find(clusterID)==fClusters.end()) {
      // new cluster
      fClusters[clusterID]=AliHLTTPCHWCFSpacePointProperties(pDecoder.get(), i);
      count++;
    } else {
      HLTError("cluster with ID 0x%08x already existing, skipping cluster %d of data block 0x%08x",
	       clusterID, i, pDesc->fSpecification);
    }
  }

  fDecoders.push_back(pDecoder.release());

  return count;
}

int AliHLTTPCHWCFSpacePointContainer::GetClusterIDs(vector<AliHLTUInt32_t>& tgt) const
{
  // get array of cluster IDs
  tgt.clear();
  transform(fClusters.begin(), fClusters.end(), back_inserter(tgt), HLT::AliGetKey());
  return tgt.size();
}

bool AliHLTTPCHWCFSpacePointContainer::Check(AliHLTUInt32_t clusterID) const
{
  // check if the cluster is available
  return fClusters.find(clusterID)!=fClusters.end();
}

const vector<AliHLTUInt32_t>* AliHLTTPCHWCFSpacePointContainer::GetClusterIDs(AliHLTUInt32_t mask)
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
  //HLTInfo("creating collection 0x%08x", mask);

  // the first cluster with number 0 has equal ID to mask unless
  // the mask selects multiple partitions/slices
  std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(mask);
  bool bAll=false;
  if (slice>=(unsigned)AliHLTTPCTransform::GetNSlice() ||
      partition>=(unsigned)AliHLTTPCTransform::GetNumberOfPatches()) {
    cl=fClusters.begin();
    bAll=true;
  }
  for (; cl!=fClusters.end(); cl++) {
    UInt_t s=AliHLTTPCSpacePointData::GetSlice(cl->first);
    UInt_t p=AliHLTTPCSpacePointData::GetPatch(cl->first);
    if ((slice>=(unsigned)AliHLTTPCTransform::GetNSlice() || s==slice) && 
	(partition>=(unsigned)AliHLTTPCTransform::GetNumberOfPatches() || p==partition)) {
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

float AliHLTTPCHWCFSpacePointContainer::GetX(AliHLTUInt32_t clusterID) const
{
  // get X coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.Decoder()==NULL) return 0.0;
  // FIXME: understand deviation from the nominal x value
  // there is a small deviation in the x coordinate - padrow number correlation
  // in principle, the clusterfinder only uses the mapping to set the x parameter.
  // now extracting the x value from the padrow no.
  //return cl->second.Decoder()->fX;
  int index=AliHLTTPCSpacePointData::GetNumber(cl->first);
  return AliHLTTPCTransform::Row2X(cl->second.Decoder()->GetPadRow(index));
}

float AliHLTTPCHWCFSpacePointContainer::GetXWidth(AliHLTUInt32_t clusterID) const
{
  // get error for X coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.Decoder()==NULL) return 0.0;
  return 0.0; // fixed in padrow number
}

float AliHLTTPCHWCFSpacePointContainer::GetY(AliHLTUInt32_t clusterID) const
{
  // get Y coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.Decoder()==NULL) return 0.0;
  int index=AliHLTTPCSpacePointData::GetNumber(cl->first);
  return cl->second.Decoder()->GetPad(index);
}

float AliHLTTPCHWCFSpacePointContainer::GetYWidth(AliHLTUInt32_t clusterID) const
{
  // get error for Y coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.Decoder()==NULL) return 0.0;
  int index=AliHLTTPCSpacePointData::GetNumber(cl->first);
  return cl->second.Decoder()->GetSigmaY2(index);
}

float AliHLTTPCHWCFSpacePointContainer::GetZ(AliHLTUInt32_t clusterID) const
{
  // get Z coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.Decoder()==NULL) return 0.0;
  int index=AliHLTTPCSpacePointData::GetNumber(cl->first);
  return cl->second.Decoder()->GetTime(index);
}

float AliHLTTPCHWCFSpacePointContainer::GetZWidth(AliHLTUInt32_t clusterID) const
{
  // get error for Z coordinate
  std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.Decoder()==NULL) return 0.0;
  int index=AliHLTTPCSpacePointData::GetNumber(cl->first);
  return cl->second.Decoder()->GetSigmaZ2(index);
}

float AliHLTTPCHWCFSpacePointContainer::GetCharge(AliHLTUInt32_t clusterID) const
{
  // get charge
  std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(clusterID);
  if (cl==fClusters.end() ||
      cl->second.Decoder()==NULL) return 0.0;
  int index=AliHLTTPCSpacePointData::GetNumber(cl->first);
  return cl->second.Decoder()->GetCharge(index);
}

float AliHLTTPCHWCFSpacePointContainer::GetPhi(AliHLTUInt32_t clusterID) const
{
  // get charge

  // phi can be derived directly from the id, no need to search
  // for existing cluster
  int slice=AliHLTTPCSpacePointData::GetSlice(clusterID);
  return ( slice + 0.5 ) * TMath::Pi() / 9.0;
}

void AliHLTTPCHWCFSpacePointContainer::Clear(Option_t * option)
{
  // clear the object and reset pointer references
  fClusters.clear();

  for (std::map<AliHLTUInt32_t, vector<AliHLTUInt32_t>*>::iterator selection=fSelections.begin();
       selection!=fSelections.end(); selection++) {
    if (selection->second) delete selection->second;
  }
  fSelections.clear();

  for (std::vector<AliHLTTPCHWCFData*>::iterator decoder=fDecoders.begin();
       decoder!=fDecoders.end(); decoder++) {
    if (*decoder) delete *decoder;
  }
  fDecoders.clear();

  AliHLTSpacePointContainer::Clear(option);
}

void AliHLTTPCHWCFSpacePointContainer::Print(ostream& out, Option_t */*option*/) const
{
  // print to stream
  out << "AliHLTTPCHWCFSpacePointContainer::Print" << endl;
  out << "n clusters: " << fClusters.size() << endl;
  for (std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.begin();
       cl!=fClusters.end(); cl++) {
    out << " " << cl->first << cl->second << endl;
  }
}

AliHLTSpacePointContainer* AliHLTTPCHWCFSpacePointContainer::SelectByMask(AliHLTUInt32_t mask, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a space point mask
  std::auto_ptr<AliHLTTPCHWCFSpacePointContainer> c(new AliHLTTPCHWCFSpacePointContainer);
  if (!c.get()) return NULL;

  UInt_t slice=AliHLTTPCSpacePointData::GetSlice(mask);
  UInt_t partition=AliHLTTPCSpacePointData::GetPatch(mask);
  for (std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.begin();
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

AliHLTSpacePointContainer* AliHLTTPCHWCFSpacePointContainer::SelectByTrack(int trackId, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a specific track
  std::auto_ptr<AliHLTTPCHWCFSpacePointContainer> c(new AliHLTTPCHWCFSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties, int>(&AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::GetTrackId,trackId));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCHWCFSpacePointContainer::SelectByMC(int mcId, bool /*bAlloc*/) const
{
  /// create a collection of clusters for a specific MC track
  std::auto_ptr<AliHLTTPCHWCFSpacePointContainer> c(new AliHLTTPCHWCFSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties, int>(&AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::GetMCId,mcId));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCHWCFSpacePointContainer::UsedClusters(bool /*bAlloc*/) const
{
  /// create a collection of all used clusters
  std::auto_ptr<AliHLTTPCHWCFSpacePointContainer> c(new AliHLTTPCHWCFSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties, bool>(&AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::IsUsed,true));
  return c.release();
}

AliHLTSpacePointContainer* AliHLTTPCHWCFSpacePointContainer::UnusedClusters(bool /*bAlloc*/) const
{
  /// create a collection of all unused clusters
  std::auto_ptr<AliHLTTPCHWCFSpacePointContainer> c(new AliHLTTPCHWCFSpacePointContainer);
  if (!c.get()) return NULL;

  HLT::copy_map_if(fClusters.begin(), fClusters.end(), c->fClusters, HLT::AliHLTUnaryPredicate<AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties, bool>(&AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::IsUsed,false));
  return c.release();
}

int AliHLTTPCHWCFSpacePointContainer::MarkUsed(const AliHLTUInt32_t* clusterIDs, int arraySize)
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

int AliHLTTPCHWCFSpacePointContainer::SetTrackID(int trackID, const AliHLTUInt32_t* clusterIDs, int arraySize)
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

int AliHLTTPCHWCFSpacePointContainer::GetTrackID(AliHLTUInt32_t clusterID) const
{
  /// get track id for specified cluster
  map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator element=fClusters.find(clusterID);
  if (element==fClusters.end()) return -1;
  return element->second.GetTrackId();
}

int AliHLTTPCHWCFSpacePointContainer::SetMCID(int mcID, const AliHLTUInt32_t* clusterIDs, int arraySize)
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

int AliHLTTPCHWCFSpacePointContainer::Write(AliHLTUInt8_t* outputPtr,
					    AliHLTUInt32_t size,
					    AliHLTComponentBlockDataList&
					    outputBlocks,
					    AliHLTDataDeflater* pDeflater,
					    const char* /*option*/) const
{
  /// write blocks to HLT component output
  if (!outputPtr) return -EINVAL;
  int iResult=0;
  AliHLTUInt32_t capacity=size;
  size=0;

  for (int slice=0; slice<AliHLTTPCTransform::GetNSlice() && iResult>=0; slice++) {
    for (int part=0; part<AliHLTTPCTransform::GetNPatches() && iResult>=0; part++) {
      AliHLTUInt32_t mask=AliHLTTPCSpacePointData::GetID(slice,part,0);
      // FIXME: make GetClusterIDs a const function and handle the cast there
      const vector<AliHLTUInt32_t>* collection=const_cast<AliHLTTPCHWCFSpacePointContainer*>(this)->GetClusterIDs(mask);
      if (!collection) continue;
      if (size+sizeof(AliHLTTPCRawClusterData)+collection->size()*sizeof(AliHLTTPCRawCluster)>capacity) {
	ALIHLTERRORGUARD(1,"too little space to write cluster output block");
	iResult=-ENOSPC;
	break;
      }
      AliHLTTPCRawClusterData* blockout=reinterpret_cast<AliHLTTPCRawClusterData*>(outputPtr+size);
      blockout->fVersion=0;
      blockout->fCount=0;

      if (pDeflater) {
	pDeflater->Clear();
	pDeflater->InitBitDataOutput(reinterpret_cast<AliHLTUInt8_t*>(blockout->fClusters), capacity-size-sizeof(AliHLTTPCRawClusterData));
	blockout->fVersion=pDeflater->GetDeflaterVersion();
      }

      unsigned lastPadRow=0;
      vector<AliHLTUInt32_t>::const_iterator clusterID=collection->begin();
      if (clusterID!=collection->end()) {
	std::map<AliHLTUInt32_t, AliHLTTPCHWCFSpacePointProperties>::const_iterator cl=fClusters.find(*clusterID);
	for (; clusterID!=collection->end(); clusterID++, (cl!=fClusters.end())?cl++:cl) {
	  if (cl!=fClusters.end() && cl->first!=*clusterID) cl=fClusters.find(*clusterID);
	  if (cl==fClusters.end() || cl->second.Decoder()==NULL) continue;
	  int index=AliHLTTPCSpacePointData::GetNumber(cl->first);
	  int padrow=cl->second.Decoder()->GetPadRow(index);
	  if (padrow<0) {
	    // something wrong here, padrow is stored in the cluster header
	    // word which has bit pattern 0x3 in bits bit 30 and 31 which was
	    // not recognized
	    ALIHLTERRORGUARD(1, "can not read cluster header word");
	    break;
	  }

	  // FIXME: the HW ClusterFinder returns only the sum
	  // sum(q_i*pad_i*pad_i)/sum(q_i)
	  // where the mean needs to be subtracted, not yet in the decoder
	  // but should be implemented there
	  float pad =cl->second.Decoder()->GetPad(index);
	  float time =cl->second.Decoder()->GetTime(index);
	  float sigmaY2=cl->second.Decoder()->GetSigmaY2(index);
	  float sigmaZ2=cl->second.Decoder()->GetSigmaZ2(index);
	  sigmaY2-=pad*pad;
	  sigmaZ2-=time*time;

	  if (!pDeflater) {
	    AliHLTTPCRawCluster& c=blockout->fClusters[blockout->fCount];
	    padrow+=AliHLTTPCTransform::GetFirstRow(part);
	    c.SetPadRow(padrow);
	    c.SetCharge(cl->second.Decoder()->GetCharge(index));
	    c.SetPad(pad);  
	    c.SetTime(time);
	    c.SetSigmaY2(sigmaY2);
	    c.SetSigmaZ2(sigmaZ2);
	    c.SetQMax(cl->second.Decoder()->GetQMax(index));
	  } else {
	    AliHLTUInt64_t padrow64=cl->second.Decoder()->GetPadRow(index);
	    // enable if padrows are ordered
	    if (padrow64>=lastPadRow) {
	    //   padrow64-=lastPadRow;
	    //   lastPadRow+=padrow64;
	    } else {
	    //   AliFatal("padrows not ordered");
	    }

	    AliHLTUInt64_t pad64
	      =(AliHLTUInt64_t)round(pad*AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kPad].fScale);
	    AliHLTUInt64_t time64
	      =(AliHLTUInt64_t)round(time*AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kTime].fScale);
	    AliHLTUInt64_t sigmaY264
	      =(AliHLTUInt64_t)round(sigmaY2*AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kSigmaY2].fScale);
	    AliHLTUInt64_t sigmaZ264
	      =(AliHLTUInt64_t)round(sigmaZ2*AliHLTTPCDefinitions::fgkClusterParameterDefinitions[AliHLTTPCDefinitions::kSigmaZ2].fScale);
	    pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kPadRow , padrow64);
	    pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kPad    , pad64);  
	    pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kTime   , time64);
	    pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kSigmaY2, sigmaY264);
	    pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kSigmaZ2, sigmaZ264);
	    pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kCharge , cl->second.Decoder()->GetCharge(index));
	    pDeflater->OutputParameterBits(AliHLTTPCDefinitions::kQMax   , cl->second.Decoder()->GetQMax(index));
	  }
	  blockout->fCount++;
	}
      }
      AliHLTComponent_BlockData bd;
      AliHLTComponent::FillBlockData(bd);
      bd.fOffset        = size;
      if (!pDeflater) {
	bd.fSize        = sizeof(AliHLTTPCRawClusterData)+blockout->fCount*sizeof(AliHLTTPCRawCluster);
      } else {
	pDeflater->Pad8Bits();
	bd.fSize        = sizeof(AliHLTTPCRawClusterData)+pDeflater->GetBitDataOutputSizeBytes();
	pDeflater->CloseBitDataOutput();
      }
      bd.fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(slice, slice, part, part);
      bd.fDataType      = AliHLTTPCDefinitions::fgkRawClustersDataType;
      outputBlocks.push_back(bd);
      
      size += bd.fSize;
    }
  }

  if (iResult<0) return iResult;
  return size;
}

AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::AliHLTTPCHWCFSpacePointProperties()
  : fDecoder(NULL)
  , fIndex(-1)
  , fUsed(false)
  , fTrackId(-1)
  , fMCId(-1)
{
  // constructor
}

AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::AliHLTTPCHWCFSpacePointProperties(const AliHLTTPCHWCFData* pDecoder, int index)
  : fDecoder(pDecoder)
  , fIndex(index)
  , fUsed(false)
  , fTrackId(-1)
  , fMCId(-1)
{
  // constructor
}

AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::AliHLTTPCHWCFSpacePointProperties(const AliHLTTPCHWCFSpacePointProperties& src)
  : fDecoder(src.fDecoder)
  , fIndex(src.fIndex)
  , fUsed(src.fUsed)
  , fTrackId(src.fTrackId)
  , fMCId(src.fMCId)
{
  // copy constructor
}

AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties& AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::operator=(const AliHLTTPCHWCFSpacePointProperties& src)
{
  // assignment operator
  if (&src==this) return *this;
  fDecoder=src.fDecoder;
  fIndex=src.fIndex;
  fUsed=src.fUsed;
  fTrackId=src.fTrackId;
  fMCId=src.fMCId;

  return *this;
}

void AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties::Print(ostream& out, Option_t */*option*/) const
{
  // print info
  if (!Decoder()) {
    out << "no data";
    return;
  }
  const AliHLTTPCHWCFData* decoder=Decoder();
  out << " " << decoder->GetPadRow(fIndex) << " " << decoder->GetPad(fIndex) << " " << decoder->GetTime(fIndex)
      << " " << decoder->GetSigmaY2(fIndex) << " "  << decoder->GetSigmaZ2(fIndex)
      << " " << decoder->GetCharge(fIndex) << " "  << decoder->GetQMax(fIndex)
      << " " << fTrackId << " " << fMCId << " " << fUsed;
}

ostream& operator<<(ostream &out, const AliHLTTPCHWCFSpacePointContainer::AliHLTTPCHWCFSpacePointProperties& p)
{
  p.Print(out);
  return out;
}
