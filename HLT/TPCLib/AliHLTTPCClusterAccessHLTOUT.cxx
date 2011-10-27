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
#include "AliHLTTPCDataCompressionDecoder.h"
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
#include <iostream>
#include <iomanip>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterAccessHLTOUT)

AliHLTTPCClusterAccessHLTOUT::AliHLTTPCClusterAccessHLTOUT()
  : TObject()
  , fVerbosity(0)
  , fClusters(NULL)
  , fCurrentSector(-1)
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
  if (strcmp(name, "clusterarray")==0) {
    if (fCurrentSector<0) return NULL;
    return fClusters->GetSectorArray(fCurrentSector);
  }
  return TObject::FindObject(name);
}

void AliHLTTPCClusterAccessHLTOUT::Clear(Option_t * option)
{
  /// inherited from TObject: cleanup
  if (strcmp(option, "event")==0) {
    if (fClusters) fClusters->Clear();
    fCurrentSector=-1;
  }
}

void AliHLTTPCClusterAccessHLTOUT::Print(Option_t *option) const
{
  /// inherited from TObject
  if (fClusters) fClusters->Print(option);
}

int AliHLTTPCClusterAccessHLTOUT::ProcessClusters(const char* params)
{
  /// process the cluster data from HLTOUT and fill array
  /// the cluster data can be in many different formats, e.g.
  /// raw or compressed
  int iResult=0;
  TString strparams(params);
  int sector=-1;
  std::auto_ptr<TObjArray> tokens(strparams.Tokenize(" "));
  if (!tokens.get()) return -ENOMEM;
  for (int i=0; i< tokens->GetEntriesFast(); i++) {
    if (!tokens->At(i)) continue;
    TString argument=tokens->At(i)->GetName();
    // the offline code enumerates first the 36 inner (partitions 0+1) and then 36 outer
    // sectors (partitions 2-5)
    if (argument.BeginsWith("sector=")) {
      argument.ReplaceAll("sector=", "");
      sector=argument.Atoi();
    }
  }
  if (sector<0) {
    AliError("invalid argument, please specify \"sector=sectorno\"");
    return -EINVAL;
  }
  if (sector>=76) {
    AliError(Form("invalid sector number %d", sector));
    return -EINVAL;
  }

  if (!fClusters) {
    fClusters=new AliTPCclusterMIContainer;
  }
  if (!fClusters) return -ENOMEM;

  if (fCurrentSector>=0) {
    // cluster container already filled
    fCurrentSector=sector;
    TObjArray* pArray=fClusters->GetSectorArray(fCurrentSector);
    if (!pArray) {
      AliError(Form("can not get cluster array for sector %d", sector));
      return -ENOBUFS;
    }
    if (fVerbosity>0) AliInfo(Form("converted %d cluster(s) for sector %d", pArray->GetEntriesFast() ,sector));
    return pArray->GetEntriesFast();
  }

  // fill the cluster container
  AliHLTSystem* pSystem=AliHLTPluginBase::GetInstance();
  if (!pSystem) {
    AliError("can not access HLT system");
    return -ENODEV;
  }
  AliHLTOUT* pHLTOUT=pSystem->RequestHLTOUT();
  if (!pHLTOUT) {
    AliError("can not access HLTOUT");
    return -EACCES;
  }

  bool bNextBlock=false;
  bool bHaveLabels=false;
  bool bHaveIds=false;
  // add cluster id and mc information data blocks
  for (bNextBlock=(pHLTOUT->SelectFirstDataBlock()>=0);
       bNextBlock; bNextBlock=(pHLTOUT->SelectNextDataBlock()>=0)) {
    AliHLTComponentBlockData desc;
    // FIXME: extend HLTOUT to get the full descriptor
    const AliHLTUInt8_t* buffer=NULL;
    if ((iResult=pHLTOUT->GetDataBuffer(buffer, desc.fSize))<0) {
      continue;
    }
    desc.fPtr=(void*)buffer;
    if (pHLTOUT->GetDataBlockDescription(desc.fDataType, desc.fSpecification)<0) {
      continue;
    }
    if (desc.fDataType==AliHLTTPCDefinitions::AliHLTDataTypeClusterMCInfo()) {
      // add mc information
      if ((iResult=fClusters->AddClusterMCData(&desc))<0) {
	return iResult;
      }
      bHaveLabels=true;
    }
    if (desc.fDataType==AliHLTTPCDefinitions::RemainingClusterIdsDataType() ||
	desc.fDataType==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
      // add cluster ids
      if ((iResult=fClusters->AddClusterIds(&desc))<0) {
	return iResult;
      }
      bHaveIds=true;
    }
  }

  // read data
  iResult=-ENODATA;
  AliHLTTPCDataCompressionDecoder decoder;
  decoder.SetVerbosity(fVerbosity);
  int nExtractedClusters=0;
  for (bNextBlock=(pHLTOUT->SelectFirstDataBlock()>=0);
       bNextBlock; bNextBlock=(pHLTOUT->SelectNextDataBlock()>=0)) {
    AliHLTComponentBlockData desc;
    // FIXME: extend HLTOUT to get the full descriptor with one call
    const AliHLTUInt8_t* buffer=NULL;
    if ((iResult=pHLTOUT->GetDataBuffer(buffer, desc.fSize))<0) {
      continue;
    }
    desc.fPtr=(void*)buffer;
    if (pHLTOUT->GetDataBlockDescription(desc.fDataType, desc.fSpecification)<0) {
      continue;
    }
    if (!TestBit(kSkipPartitionClusters) &&
	(desc.fDataType==AliHLTTPCDefinitions::RemainingClustersCompressedDataType() ||
	 desc.fDataType==AliHLTTPCDefinitions::RawClustersDataType())) {
      iResult=decoder.ReadClustersPartition(fClusters->BeginRemainingClusterBlock(0, desc.fSpecification),
					    reinterpret_cast<AliHLTUInt8_t*>(desc.fPtr),
					    desc.fSize,
					    desc.fSpecification);
      if (iResult>0) nExtractedClusters+=iResult;
      continue;
    } else if (!TestBit(kSkipTrackClusters) &&
	       desc.fDataType==AliHLTTPCDefinitions::ClusterTracksCompressedDataType()) {
      iResult=decoder.ReadTrackModelClustersCompressed(fClusters->BeginTrackModelClusterBlock(0),
							reinterpret_cast<AliHLTUInt8_t*>(desc.fPtr),
							desc.fSize,
							desc.fSpecification);
      continue;
    }
  }

  pSystem->ReleaseHLTOUT(pHLTOUT);

  if (iResult<0) return iResult;
  if (fVerbosity>0) {
    int nConvertedClusters=0;
    for (int s=0; s<72; s++) {
      TObjArray* pArray=fClusters->GetSectorArray(s);
      if (!pArray) continue;
      nConvertedClusters+=pArray->GetEntriesFast();
    }
    AliInfo(Form("extracted HLT clusters: %d, converted HLT clusters: %d", nExtractedClusters, nConvertedClusters));
  }

  fCurrentSector=sector;
  TObjArray* pArray=fClusters->GetSectorArray(fCurrentSector);
  if (!pArray) {
    AliError(Form("can not get cluster array for sector %d", sector));
    return -ENOBUFS;
  }
  if (fVerbosity>0) AliInfo(Form("converted %d cluster(s) for sector %d", pArray->GetEntriesFast() ,sector));
  return pArray->GetEntriesFast();
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

AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::AliTPCclusterMIContainer()
  : fClusterArrays()
  , fRemainingClusterIds()
  , fTrackModelClusterIds()
  , fCurrentClusterIds(NULL)
  , fClusterMCData()
  , fIterator()

{
  /// constructor
  for (int i=0; i<72; i++) {
    fClusterArrays.push_back(new TClonesArray("AliTPCclusterMI"));
  }
}

AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::~AliTPCclusterMIContainer()
{
  /// dectructor
  for (vector<TClonesArray*>::iterator i=fClusterArrays.begin(); i!=fClusterArrays.end(); i++) {
    if (*i) {
      (*i)->Clear();
      delete *i;
    }
  }
}

AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::iterator& AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::BeginRemainingClusterBlock(int /*count*/, AliHLTUInt32_t specification)
{
  /// iterator of remaining clusters block of specification
  AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(specification);
  AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(specification);
  unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
  if (index<fRemainingClusterIds.size())
    fCurrentClusterIds=&fRemainingClusterIds[index];
  else
    fCurrentClusterIds=NULL;
  fIterator=iterator(this);
  return fIterator;
}

AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::iterator& AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::BeginTrackModelClusterBlock(int /*count*/)
{
  /// iterator of track model clusters
  if (fTrackModelClusterIds.fIds && fTrackModelClusterIds.fSize>0)
    fCurrentClusterIds=&fTrackModelClusterIds;
  else
    fCurrentClusterIds=NULL;
  fIterator=iterator(this);
  return fIterator;
}

int AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::AddClusterMCData(const AliHLTComponentBlockData* pDesc)
{
  /// add cluster mc data block
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType==AliHLTTPCDefinitions::AliHLTDataTypeClusterMCInfo()) {
    AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(pDesc->fSpecification);
    AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(pDesc->fSpecification);
    unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
    if (fClusterMCData.size()<=index) {
      if ((int)fClusterMCData.size()<AliHLTTPCTransform::GetNSlice()*AliHLTTPCTransform::GetNumberOfPatches()) {
	fClusterMCData.resize(AliHLTTPCTransform::GetNSlice()*AliHLTTPCTransform::GetNumberOfPatches(), NULL);
      } else {
	fClusterMCData.resize(index+1, NULL);
      }
    }
    if (pDesc->fSize<sizeof(AliHLTTPCClusterMCData)) return -EINVAL;
    const AliHLTTPCClusterMCData* pData=reinterpret_cast<const AliHLTTPCClusterMCData*>(pDesc->fPtr);
    unsigned nLabels = pData->fCount;
    if (nLabels*sizeof(AliHLTTPCClusterMCLabel) + sizeof(AliHLTTPCClusterMCData) != pDesc->fSize) {
      return -EINVAL;
    }
    fClusterMCData[index]=pData;
    return 0;
  }
  return -ENODATA;
}

int AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::AddClusterIds(const AliHLTComponentBlockData* pDesc)
{
  /// add cluster id block for remaining or track model clusters
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
    fTrackModelClusterIds.fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fTrackModelClusterIds.fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  if (pDesc->fDataType==AliHLTTPCDefinitions::RemainingClusterIdsDataType()) {
    AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(pDesc->fSpecification);
    AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(pDesc->fSpecification);
    unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
    if (fRemainingClusterIds.size()<=index) {
      if ((int)fRemainingClusterIds.size()<AliHLTTPCTransform::GetNSlice()*AliHLTTPCTransform::GetNumberOfPatches()) {
	fRemainingClusterIds.resize(AliHLTTPCTransform::GetNSlice()*AliHLTTPCTransform::GetNumberOfPatches());
      } else {
	fRemainingClusterIds.resize(index+1);
      }
    }
    fRemainingClusterIds[index].fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fRemainingClusterIds[index].fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  return -ENODATA;
}

AliHLTUInt32_t AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::GetClusterId(int clusterNo) const
{
  /// get the cluster id from the current cluster id block (optional)
  if (!fCurrentClusterIds ||
      (int)fCurrentClusterIds->fSize<=clusterNo ||
      clusterNo<0)
    return kAliHLTVoidDataSpec;
  return fCurrentClusterIds->fIds[clusterNo];
}

AliTPCclusterMI* AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::NextCluster(int slice, int partition)
{
  /// load next cluster from array of the sepcific sector
  unsigned sector=partition<2?slice:slice+36;
  if (fClusterArrays.size()<=sector ||
      fClusterArrays[sector]==NULL) {
    AliErrorClass(Form("no cluster array available for sector %d", sector));
    return NULL;
  }
  TClonesArray& array=*(fClusterArrays[sector]);
  int count=array.GetEntriesFast();
  return new (array[count]) AliTPCclusterMI;
}

int AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::SetMC(AliTPCclusterMI* pCluster, AliHLTUInt32_t clusterId)
{
  /// set MC data for the cluster
  if (!pCluster) return -EINVAL;
  if (clusterId==kAliHLTVoidDataSpec) return 0;

  unsigned slice=AliHLTTPCSpacePointData::GetSlice(clusterId);
  unsigned partition=AliHLTTPCSpacePointData::GetPatch(clusterId);
  unsigned number=AliHLTTPCSpacePointData::GetNumber(clusterId);
  if ((int)slice>=AliHLTTPCTransform::GetNSlice() ||
      (int)partition>=AliHLTTPCTransform::GetNumberOfPatches()) return -EDOM;
  unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
  if (fClusterMCData.size()<=index ||
      fClusterMCData[index]==NULL ||
      fClusterMCData[index]->fCount<=number) return 0;
  const AliHLTTPCClusterMCWeight* mcWeights=fClusterMCData[index]->fLabels[number].fClusterID;
  for (int k=0; k<3; k++) {
    // TODO: sort the labels according to the weight in order to assign the most likely mc label
    // to the first component 
    pCluster->SetLabel(mcWeights[k].fMCID, k);
  }

  return 0;
}

void  AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::Clear(Option_t* /*option*/)
{
  /// internal cleanup
  {
    for (vector<TClonesArray*>::iterator i=fClusterArrays.begin(); i!=fClusterArrays.end(); i++)
      if (*i) (*i)->Clear();
  }
  {
    for (vector<AliClusterIdBlock>::iterator i=fRemainingClusterIds.begin(); i!=fRemainingClusterIds.end(); i++)
      {i->fIds=NULL; i->fSize=0;}
  }
  fTrackModelClusterIds.fIds=NULL; fTrackModelClusterIds.fSize=0;
  fCurrentClusterIds=NULL;
  {
    for (vector<const AliHLTTPCClusterMCData*>::iterator i=fClusterMCData.begin(); i!=fClusterMCData.end(); i++)
      *i=NULL;
  }
}

TObjArray* AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::GetSectorArray(unsigned sector) const
{
  /// get the cluster array for a sector
  if (fClusterArrays.size()<=sector) return NULL;
  return fClusterArrays[sector];
}

void AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::Print(Option_t *option) const
{
  /// inherited from TObject
  cout << "AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer" << endl;
  ios::fmtflags coutflags=cout.flags(); // backup cout status flags
  bool bAll=false;
  if ((bAll=(strcmp(option, "full")==0)) ||
      strcmp(option, "short")==0) {
    for (unsigned iArray=0; iArray<fClusterArrays.size(); iArray++) {
      if (fClusterArrays[iArray]) {
	TClonesArray* pArray=fClusterArrays[iArray];
	cout << "  sector " << setfill(' ') << setw(2) << iArray << ": " << pArray->GetEntriesFast() << endl;
	if (bAll) {
	  for (int iCluster=0; iCluster<pArray->GetEntriesFast(); iCluster++) {
	    if (!pArray->At(iCluster)) continue;
	    AliTPCclusterMI* pCluster=dynamic_cast<AliTPCclusterMI*>(pArray->At(iCluster));
	    if (!pCluster) break;
	    cout << "    AliTPCclusterMI:"
		 << "  row="    << pCluster->GetRow() 
		 << "  pad="    << pCluster->GetPad()
		 << "  time="   << pCluster->GetTimeBin()
		 << "  charge=" << pCluster->GetQ()
		 << "  maxq="   << pCluster->GetMax()
		 << endl;
	  }
	}
      }
    }
  }
  cout.flags(coutflags); // restore the original flags
}

AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::iterator& AliHLTTPCClusterAccessHLTOUT::AliTPCclusterMIContainer::iterator::Next(int slice, int partition)
{
  // switch to next cluster
  if (!fData) {
    fCluster=NULL;
    fClusterId=kAliHLTVoidDataSpec;
    return *this;
  }
  if (fClusterNo>=0 && !fCluster) {
    // end was reached before
    return *this;
  }
  fCluster=fData->NextCluster(slice, partition);
  fClusterId=fData->GetClusterId(++fClusterNo);
  if (fCluster && fClusterId!=kAliHLTVoidDataSpec) {
    fData->SetMC(fCluster, fClusterId);
  }
  // offline uses row number in physical sector, inner sector consists of
  // partitions 0 and 1, outer sector of partition 2-5
  fRowOffset=partition<2?0:AliHLTTPCTransform::GetFirstRow(2);
  return *this;
}
