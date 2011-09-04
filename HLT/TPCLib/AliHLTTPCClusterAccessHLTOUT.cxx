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

/// @file   AliHLTTPCClusterAccessHLTOUT.h
/// @author Matthias Richter
/// @date   2011-06-06
/// @brief  Interface to HLT TPC clusters
///

#include "AliHLTTPCClusterAccessHLTOUT.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTOUT.h"
#include "AliHLTComponent.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataInflater.h"
#include "AliHLTTPCDefinitions.h"
#include "AliLog.h"
#include "AliHLTSystem.h"
#include "AliHLTPluginBase.h"
#include "AliTPCclusterMI.h"
#include "TClonesArray.h"
#include <cstdlib>
#include <string>
#include <memory>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterAccessHLTOUT)

AliHLTTPCClusterAccessHLTOUT::AliHLTTPCClusterAccessHLTOUT()
  : TObject()
  , fVerbosity(0)
  , fClusters(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCClusterAccessHLTOUT::~AliHLTTPCClusterAccessHLTOUT()
{
  // destructor
  if (fClusters) {
    fClusters->Clear();
    delete fClusters;
    fClusters=NULL;
  }
}

void AliHLTTPCClusterAccessHLTOUT::Execute(const char *method,  const char *params, Int_t *error)
{
  /// inherited from TObject: abstract command interface
  if (strcmp(method, "read")==0) {
    int iResult=ProcessClusters(params);
    if (error) *error=iResult;
    return;
  }
  if (strcmp(method, "verbosity")==0) {
    int iResult=0;
    if (params) {
      char* dummy;
      int value=strtol(params, &dummy, 0);
      if (dummy==NULL) {
	fVerbosity=value;
      } else {
	AliError("invalid argument for command 'verbosity', expecting string with number");
      }
    } else {
      iResult=-EINVAL;
    }
    if (error) *error=iResult;
    return;
  }
}

TObject* AliHLTTPCClusterAccessHLTOUT::FindObject(const char *name) const
{
  /// inherited from TObject: return the cluster array if name id "clusterarray"
  if (strcmp(name, "clusterarray")==0) return fClusters;
  return TObject::FindObject(name);
}

void AliHLTTPCClusterAccessHLTOUT::Clear(Option_t * /*option*/)
{
  /// inherited from TObject: cleanup
  if (fClusters) fClusters->Clear();
}

void AliHLTTPCClusterAccessHLTOUT::Print(Option_t */*option*/) const
{
  /// inherited from TObject
  if (!fClusters) return;
  for (int i=0; i<fClusters->GetEntriesFast(); i++) {
    if (!fClusters->At(i)) continue;
    AliTPCclusterMI* pCluster=dynamic_cast<AliTPCclusterMI*>(fClusters->At(i));
    if (!pCluster) break;
    cout << "AliTPCclusterMI:"
	 << "  row="    << pCluster->GetRow() 
	 << "  pad="    << pCluster->GetPad()
	 << "  time="   << pCluster->GetTimeBin()
	 << "  charge=" << pCluster->GetQ()
	 << "  maxq="   << pCluster->GetMax()
	 << endl;
  }
}

int AliHLTTPCClusterAccessHLTOUT::ProcessClusters(const char* params)
{
  /// process the cluster data from HLTOUT and fill array
  /// the cluster data can be in many different formats, e.g.
  /// raw or compressed
  int iResult=0;
  TString strparams(params);
  int minSlice=0, maxSlice=35, minPart=0, maxPart=5;
  std::auto_ptr<TObjArray> tokens(strparams.Tokenize(" "));
  if (!tokens.get()) return -ENOMEM;
  for (int i=0; i< tokens->GetEntriesFast(); i++) {
    if (!tokens->At(i)) continue;
    TString argument=tokens->At(i)->GetName();
    if (argument.BeginsWith("sector=")) {
      argument.ReplaceAll("sector=", "");
      int sector=argument.Atoi();
      // the offline code enumerates first the 36 inner (partitions 0+1) and then 36 outer
      // sectors (partitions 2-5)
      if (fVerbosity>0) AliInfo(Form("processing HLT clusters for sector %d", sector));
      if (sector<36) { // inner sectors
	minSlice=maxSlice=sector;
	minPart=0; maxPart=1;
      } else { // outer sectors
	minSlice=maxSlice=sector-36;
	minPart=2; maxPart=5;
      }
    }
  }

  if (!fClusters) {
    fClusters=new TClonesArray("AliTPCclusterMI");
  }
  if (!fClusters) return -ENOMEM;

  AliHLTSystem* pSystem=AliHLTPluginBase::GetInstance();
  if (!pSystem) {
    AliError("can not access HLT system");
    return -ENODEV;
  }
  AliHLTOUT* pHLTOUT=pSystem->RequestHLTOUT();
  if (!pHLTOUT) {
    AliError("can not access HLTOUT");
    return -ENODEV;
  }

  for (int slice=minSlice; slice<=maxSlice; slice++) {
    for (int part=minPart; part<=maxPart; part++) {
      if (fVerbosity>0) AliInfo(Form("processing HLT clusters for slice %d partitions %d", slice, part));
      AliHLTUInt32_t spec=slice<<24 | slice<<16 | part<<8 | part;
      AliHLTTPCClusterMCDataList tpcClusterLabels;
      bool bHaveLabels=false;
      if (pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo, spec)>=0) {
	iResult=ReadAliHLTTPCClusterMCData(pHLTOUT, tpcClusterLabels);
	bHaveLabels=true;
      }

      if (pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::RemainingClustersCompressedDataType(), spec)>=0) {
	iResult=ReadRemainingClustersCompressed(pHLTOUT, fClusters, bHaveLabels?&tpcClusterLabels:NULL);
      } else if (pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::fgkRawClustersDataType, spec)>=0) {
	iResult=ReadAliHLTTPCRawClusterData(pHLTOUT, fClusters, bHaveLabels?&tpcClusterLabels:NULL);
      } else if (pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::fgkClustersDataType, spec)>=0) {
	ALIHLTERRORGUARD(1, "HLTOUT data contains tarnsformed TPC clusters instead of raw TPC clusters, can not create clusters for reconstruction");
      }
    }
  }

  pSystem->ReleaseHLTOUT(pHLTOUT);
  return iResult;
}

int AliHLTTPCClusterAccessHLTOUT::ReadAliHLTTPCClusterMCData(AliHLTOUT* pHLTOUT, AliHLTTPCClusterMCDataList &tpcClusterLabels) const
{
  // read cluster data from AliHLTTPCClusterData
  int iResult=0;
  if (!pHLTOUT) return -EINVAL;
  do {
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    if ((iResult=pHLTOUT->GetDataBuffer(pBuffer, size))<0) {
      continue;
    }
    if (pBuffer==NULL || size<4) {
      AliError("invalid cluster mc data block");
      continue;
    }
    const AliHLTTPCClusterMCData* clusterMCData = reinterpret_cast<const AliHLTTPCClusterMCData*>(pBuffer);
    Int_t nLabels = (Int_t) clusterMCData->fCount;
    if (nLabels*sizeof(AliHLTTPCClusterMCLabel) + sizeof(AliHLTTPCClusterMCData) != size) {
      AliError("inconsistent cluster mc data block size, skipping block");
      continue;
    }
    // id of the cluster is 
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t specification=kAliHLTVoidDataSpec;
    if (pHLTOUT->GetDataBlockDescription(dt, specification)<0) {
      AliError("failed to retrieve data block description, skipping mc cluster data block ...");
      continue;
    }
    AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(specification);
    AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(specification);
    if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr(specification) ||
	partition!=AliHLTTPCDefinitions::GetMaxPatchNr(specification)) {
      AliError(Form("can not read cluster mc data block with data of multiple partitions, skipping block %s %08x",
		    AliHLTComponent::DataType2Text(dt).c_str(), specification));
      continue;
    }
    const AliHLTTPCClusterMCLabel *labels = clusterMCData->fLabels;
    for (int i=0; i<nLabels; i++) {
      AliHLTUInt32_t id=AliHLTTPCSpacePointData::GetID(slice, partition, i);
      if (tpcClusterLabels.find(id)==tpcClusterLabels.end()) {
	// new cluster
	tpcClusterLabels[id]=labels[i];
      } else {
	AliError(Form("cluster with ID 0x%08x already existing, skipping cluster %d of data block 0x%08x",
		      id, i, specification));
      }
    }
  } while (pHLTOUT->SelectNextDataBlock()>=0);
  return iResult;
}

int AliHLTTPCClusterAccessHLTOUT::ReadAliHLTTPCClusterData(AliHLTOUT* pHLTOUT, TClonesArray* pClusters, const AliHLTTPCClusterMCDataList *tpcClusterLabels) const
{
  // read cluster data from AliHLTTPCClusterData
  int iResult=0;
  if (!pHLTOUT || !pClusters) return -EINVAL;
  do {
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    if ((iResult=pHLTOUT->GetDataBuffer(pBuffer, size))<0) {
      continue;
    }
    if (pBuffer==NULL || size<4) {
      AliError("invalid cluster data block");
      continue;
    }
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t specification=kAliHLTVoidDataSpec;
    if (pHLTOUT->GetDataBlockDescription(dt, specification)<0) {
      AliError("failed to retrieve data block description, skipping mc cluster data block ...");
      continue;
    }
    const AliHLTTPCClusterData* clusterData = reinterpret_cast<const AliHLTTPCClusterData*>(pBuffer);
    Int_t nSpacepoints = (Int_t) clusterData->fSpacePointCnt;
    if (nSpacepoints*sizeof(AliHLTTPCSpacePointData) + sizeof(AliHLTTPCClusterData) != size) {
      AliError("inconsistent cluster data block size, skipping block");
      continue;
    }
    const AliHLTTPCSpacePointData *clusters = clusterData->fSpacePoints;
    int offset=pClusters->GetEntries();
    pClusters->ExpandCreate(offset+nSpacepoints);
    AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(specification);
    AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(specification);
    // FIXME: get first row number of outer sectors from a common definition instead using number
    unsigned rowOffset=partition<2?0:63;
    for (int i=0; i<nSpacepoints; i++) {
      if (!pClusters->At(offset+i)) continue;
      AliTPCclusterMI* pCluster=dynamic_cast<AliTPCclusterMI*>(pClusters->At(offset+i));
      if (!pCluster) {
	AliError("invalid object type, expecting AliTPCclusterMI");
	break; // this is a problem of all objects
      }
      if (clusters[i].fPadRow<rowOffset) {
	AliError(Form("invalid row number %d, expecting minimum row number %d for slice %d partition %d", clusters[i].fPadRow, rowOffset, slice, partition));
      } else {
      pCluster->SetRow(clusters[i].fPadRow-rowOffset);
      }
      pCluster->SetPad(clusters[i].fY);
      pCluster->SetTimeBin(clusters[i].fZ);
      pCluster->SetSigmaY2(clusters[i].fSigmaY2);
      pCluster->SetSigmaZ2(clusters[i].fSigmaZ2);
      pCluster->SetQ(clusters[i].fCharge);
      pCluster->SetMax(clusters[i].fQMax);
      if (tpcClusterLabels) {
	if (tpcClusterLabels->find(clusters[i].fID)!=tpcClusterLabels->end()) {
	  const AliHLTTPCClusterMCWeight* mcWeights=tpcClusterLabels->find(clusters[i].fID)->second.fClusterID;
	  for (int k=0; k<3; k++) {
	    // TODO: sort the labels according to the weight in order to assign the most likely mc label
	    // to the first component 
	    pCluster->SetLabel(mcWeights[k].fMCID, k);
	  }
	} else {
	  AliError(Form("can not find mc label of cluster with id %0x08x", clusters[i].fID));
	}
      }
    }
    if (fVerbosity>0) AliInfo(Form("converted %d cluster(s) from block %s 0x%08x", nSpacepoints, AliHLTComponent::DataType2Text(dt).c_str(), specification));
  } while (pHLTOUT->SelectNextDataBlock()>=0);
  return iResult;
}

int AliHLTTPCClusterAccessHLTOUT::ReadAliHLTTPCRawClusterData(AliHLTOUT* pHLTOUT, TClonesArray* pClusters, const AliHLTTPCClusterMCDataList *tpcClusterLabels)
{
  // read cluster data from AliHLTTPCClusterData

  // FIXME: this is in large parts like ReadAliHLTTPCClusterData,
  // make a common method
  int iResult=0;
  if (!pHLTOUT || !pClusters) return -EINVAL;
  do {
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    if ((iResult=pHLTOUT->GetDataBuffer(pBuffer, size))<0) {
      continue;
    }
    if (pBuffer==NULL || size<4) {
      AliError("invalid cluster data block");
      continue;
    }
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t specification=kAliHLTVoidDataSpec;
    if (pHLTOUT->GetDataBlockDescription(dt, specification)<0) {
      AliError("failed to retrieve data block description, skipping mc cluster data block ...");
      continue;
    }
    const AliHLTTPCRawClusterData* clusterData = reinterpret_cast<const AliHLTTPCRawClusterData*>(pBuffer);
    Int_t nCount = (Int_t) clusterData->fCount;
    if (clusterData->fVersion!=0) {
      // this is encoded data of different formats
      switch (clusterData->fVersion) {
      case 1: 
	iResult=ReadAliHLTTPCRawClusterDataDeflateSimple(reinterpret_cast<const AliHLTUInt8_t*>(clusterData->fClusters),
							 size-sizeof(AliHLTTPCRawClusterData), nCount, specification,
							 pClusters, tpcClusterLabels);
	break;
      default:
	iResult=-EPROTO;
      }
      return iResult;
    }

    if (nCount*sizeof(AliHLTTPCRawCluster) + sizeof(AliHLTTPCRawClusterData) != size) {
      AliError("inconsistent cluster data block size, skipping block");
      continue;
    }
    const AliHLTTPCRawCluster *clusters = clusterData->fClusters;
    int offset=pClusters->GetEntries();
    pClusters->ExpandCreate(offset+nCount);
    AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(specification);
    AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(specification);
    // FIXME: get first row number of outer sectors from a common definition instead using number
    int rowOffset=partition<2?0:63;
    for (int i=0; i<nCount; i++) {
      if (!pClusters->At(offset+i)) continue;
      AliTPCclusterMI* pCluster=dynamic_cast<AliTPCclusterMI*>(pClusters->At(offset+i));
      if (!pCluster) {
	AliError("invalid object type, expecting AliTPCclusterMI");
	break; // this is a problem of all objects
      }
      if (fVerbosity>1) AliInfo(Form("cluster padrow %d (slice %d partition %d)", clusters[i].GetPadRow(), slice, partition));
      if (clusters[i].GetPadRow()<rowOffset) {
	AliError(Form("invalid row number %d, expecting minimum row number %d for slice %d partition %d", clusters[i].GetPadRow(), rowOffset, slice, partition));
      } else {
      pCluster->SetRow(clusters[i].GetPadRow()-rowOffset);
      }
      pCluster->SetPad(clusters[i].GetPad());
      pCluster->SetTimeBin(clusters[i].GetTime());
      pCluster->SetSigmaY2(clusters[i].GetSigmaY2());
      pCluster->SetSigmaZ2(clusters[i].GetSigmaZ2());
      pCluster->SetQ(clusters[i].GetCharge());
      pCluster->SetMax(clusters[i].GetQMax());
      if (tpcClusterLabels) {
	UInt_t clusterID=AliHLTTPCSpacePointData::GetID(slice, partition, i);
	if (tpcClusterLabels->find(clusterID)!=tpcClusterLabels->end()) {
	  const AliHLTTPCClusterMCWeight* mcWeights=tpcClusterLabels->find(clusterID)->second.fClusterID;
	  for (int k=0; k<3; k++) {
	    // TODO: sort the labels according to the weight in order to assign the most likely mc label
	    // to the first component 
	    pCluster->SetLabel(mcWeights[k].fMCID, k);
	  }
	} else {
	  AliError(Form("can not find mc label of cluster with id %0x08x", clusterID));
	}
      }
    }
    if (fVerbosity>0) AliInfo(Form("converted %d cluster(s) from block %s 0x%08x", nCount, AliHLTComponent::DataType2Text(dt).c_str(), specification));
  } while (pHLTOUT->SelectNextDataBlock()>=0);
  return iResult;
}

int AliHLTTPCClusterAccessHLTOUT::ReadRemainingClustersCompressed(AliHLTOUT* pHLTOUT, TClonesArray* pClusters, const AliHLTTPCClusterMCDataList *tpcClusterLabels)
{
  // read cluster data from AliHLTTPCClusterData
  int iResult=0;
  if (!pHLTOUT || !pClusters) return -EINVAL;
  do {
    const AliHLTUInt8_t* pBuffer=NULL;
    AliHLTUInt32_t size=0;
    if ((iResult=pHLTOUT->GetDataBuffer(pBuffer, size))<0) {
      continue;
    }
    if (pBuffer==NULL || size<4) {
      AliError("invalid cluster data block");
      continue;
    }
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t specification=kAliHLTVoidDataSpec;
    if (pHLTOUT->GetDataBlockDescription(dt, specification)<0) {
      AliError("failed to retrieve data block description, skipping mc cluster data block ...");
      continue;
    }
    const AliHLTTPCRawClusterData* clusterData = reinterpret_cast<const AliHLTTPCRawClusterData*>(pBuffer);
    Int_t nCount = (Int_t) clusterData->fCount;

    // this is encoded data of different formats
    switch (clusterData->fVersion) {
    case 1: 
      iResult=ReadAliHLTTPCRawClusterDataDeflateSimple(reinterpret_cast<const AliHLTUInt8_t*>(clusterData->fClusters),
						       size-sizeof(AliHLTTPCRawClusterData), nCount, specification,
						       pClusters, tpcClusterLabels);
      break;
    default:
      AliError(Form("invalid cluster format version %d", clusterData->fVersion));
      iResult=-EPROTO;
    }

    if (fVerbosity>0) AliInfo(Form("converted %d cluster(s) from block %s 0x%08x", nCount, AliHLTComponent::DataType2Text(dt).c_str(), specification));
  } while (pHLTOUT->SelectNextDataBlock()>=0 && iResult>=0);

  return iResult;
}

int AliHLTTPCClusterAccessHLTOUT::ReadAliHLTTPCRawClusterDataDeflateSimple(const AliHLTUInt8_t* pData, int dataSize,
									   int nofClusters, AliHLTUInt32_t specification,
									   TClonesArray* pClusters,
									   const AliHLTTPCClusterMCDataList *tpcClusterLabels)
{
  // read cluster data from AliHLTTPCClusterData

  // FIXME: quick implementation to read the compressed cluster data from HLTOUT
  // the data definition below is the same as in AliHLTTPCDataCompressionComponent
  // but needs to be moved to a common class (AliHLTTPCDefinitions?)
  // Think about a decoder class supporting iterator objects for various types
  // of cluster data
  int iResult=0;
  if (!pData || !pClusters) return -EINVAL;
  AliHLTDataInflater inflater;
  if ((iResult=inflater.InitBitDataInput(pData, dataSize))<0) {
    return iResult;
  }

  int offset=pClusters->GetEntries();
  pClusters->ExpandCreate(offset+nofClusters);
  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(specification);
  AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(specification);
  // the compressed format stores the difference of the local row number in
  // the partition to the row of the last cluster
  // add the first row in the partition to get global row number
  // offline uses row number in physical sector, inner sector consists of
  // partitions 0 and 1, outer sector of partition 2-5
  int rowOffset=AliHLTTPCTransform::GetFirstRow(partition)-(partition<2?0:AliHLTTPCTransform::GetFirstRow(2));

  int parameterId=0;
  int outClusterCnt=0;
  AliHLTUInt8_t switchBit=0;
  AliHLTUInt64_t value=0;
  AliTPCclusterMI* pCluster=NULL;
  AliHLTUInt32_t lastPadRow=0;
  while (outClusterCnt<nofClusters && inflater.InputBit(switchBit)) {
    const AliHLTTPCDefinitions::AliClusterParameter& parameter
      =AliHLTTPCDefinitions::fgkClusterParameterDefinitions[parameterId];
    // in mode DeflaterSimple, the optional parameter of the cluster parameter definition
    // corresponds to the number bits of the reduced format
    if (!inflater.InputBits(value, switchBit?parameter.fBitLength:parameter.fOptional)) {
      break;
    }

    if (!pCluster) {
      if (!pClusters->At(offset+outClusterCnt)) {
	// here we should not get anymore because of the condition outClusterCnt<nofClusters
	return -ENOSPC;
      }
      pCluster=dynamic_cast<AliTPCclusterMI*>(pClusters->At(offset+outClusterCnt));
      if (!pCluster) {
	AliError("invalid object type, expecting AliTPCclusterMI");
	iResult=-EBADF; // this is a problem of all objects
	break;
      }
    }
    switch (parameterId) {
    case AliHLTTPCDefinitions::kPadRow:
      {pCluster->SetRow(value+lastPadRow+rowOffset); lastPadRow+=value;break;}
    case AliHLTTPCDefinitions::kPad:
      {float pad=value; pad/=parameter.fScale; pCluster->SetPad(pad); break;}
    case AliHLTTPCDefinitions::kTime:
      {float time=value; time/=parameter.fScale; pCluster->SetTimeBin(time); break;}
    case AliHLTTPCDefinitions::kSigmaY2:
      {float sigmaY2=value; sigmaY2/=parameter.fScale; pCluster->SetSigmaY2(sigmaY2); break;}
    case AliHLTTPCDefinitions::kSigmaZ2:
      {float sigmaZ2=value; sigmaZ2/=parameter.fScale; pCluster->SetSigmaZ2(sigmaZ2); break;}
    case AliHLTTPCDefinitions::kCharge:
      {pCluster->SetQ(value); break;}
    case AliHLTTPCDefinitions::kQMax:
      {pCluster->SetMax(value); break;}
    }
    if (parameterId>=AliHLTTPCDefinitions::kLast) {
      // switch to next cluster
      if (tpcClusterLabels) {
	UInt_t clusterID=AliHLTTPCSpacePointData::GetID(slice, partition, outClusterCnt);
	if (tpcClusterLabels->find(clusterID)!=tpcClusterLabels->end()) {
	  const AliHLTTPCClusterMCWeight* mcWeights=tpcClusterLabels->find(clusterID)->second.fClusterID;
	  for (int k=0; k<3; k++) {
	    // TODO: sort the labels according to the weight in order to assign the most likely mc label
	    // to the first component 
	    pCluster->SetLabel(mcWeights[k].fMCID, k);
	  }
	} else {
	  AliError(Form("can not find mc label of cluster with id 0x%08x", clusterID));
	}
      }
      outClusterCnt++;
      pCluster=NULL;
      parameterId=-1;
    }
    parameterId++;
  }
  inflater.Pad8Bits();
  if (inflater.InputBit(switchBit)) {
    AliWarning("format error of compressed clusters, there is more data than expected");
  }
  inflater.CloseBitDataInput();
  if (iResult>=0 && nofClusters!=outClusterCnt) {
    // is this a Fatal?
    AliError(Form("error reading compressed cluster format: expected %d, read only %d cluster(s)", nofClusters, outClusterCnt));
    return -EPROTO;
  }
  return iResult;
}
