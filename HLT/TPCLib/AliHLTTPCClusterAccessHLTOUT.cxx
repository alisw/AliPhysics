// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
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
#include "AliHLTTPCGeometry.h"
#include "AliHLTOUT.h"
#include "AliHLTComponent.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataInflater.h"
#include "AliLog.h"
#include "AliHLTSystem.h"
#include "AliHLTPluginBase.h"
#include "AliTPCclusterMI.h"
#include "AliTPCClustersRow.h"
#include "AliTPCParam.h"
#include "TClonesArray.h"
#include "TString.h"
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
  , fCurrentRow(-1)
  , fPropagateSplitClusterFlag(0)
  , fPropagateEdgeClusterFlag(0)
  , fMarkEdgeClusters(0)
  , fpDecoder(NULL)
  , fTPCParam(NULL)
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
  if (fpDecoder) {
    fpDecoder->Clear();
    delete fpDecoder;
    fpDecoder=NULL;
  }
  if (fTPCParam) {
    // FIXME: a copy of the TPCParam object is not possible because there is
    // no appropriate copy constructor or assignment operator, using as
    // external pointer
    //delete fTPCParam;
    fTPCParam=NULL;
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
  if (strcmp(method, "prepare_copy")==0) {
    int iResult=ScanParameters(params);
    if (error) *error=iResult;
    return;
  }
  if (strcmp(method, "get_edge_flags_set")==0) {
    *error = GetPropagateEdgeClusterFlag() || GetMarkEdgeClusterFlag();
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
    return fClusters->GetSectorArray(fCurrentSector, fPropagateSplitClusterFlag, fPropagateEdgeClusterFlag, fMarkEdgeClusters);
  }
  return TObject::FindObject(name);
}

void AliHLTTPCClusterAccessHLTOUT::Copy(TObject &object) const
{
  /// inherited from TObject: supports writing of data to AliTPCClustersRow
  AliTPCClustersRow* rowcl=dynamic_cast<AliTPCClustersRow*>(&object);
  if (rowcl) {
    fClusters->FillSectorArray(rowcl->GetArray(), fCurrentSector, fCurrentRow, fPropagateSplitClusterFlag, fPropagateEdgeClusterFlag, fMarkEdgeClusters);
    return;
  }
  return TObject::Copy(object);
}


void AliHLTTPCClusterAccessHLTOUT::Clear(Option_t * option)
{
  /// inherited from TObject: cleanup
  if (fClusters) fClusters->Clear(option);
  fCurrentSector=-1;
  fCurrentRow=-1;
}

void AliHLTTPCClusterAccessHLTOUT::Print(Option_t *option) const
{
  /// inherited from TObject
  if (fClusters) fClusters->Print(option);
}

int AliHLTTPCClusterAccessHLTOUT::ScanParameters(const char* params)
{
  TString strparams(params);
  int sector=-1;
  int row=-1;
  fCurrentSector=-1;
  fCurrentRow=-1;
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
    if (argument.BeginsWith("row=")) {
      argument.ReplaceAll("row=", "");
      row=argument.Atoi();
    }
  }
  if (sector<0) {
    AliError("invalid argument, please specify \"sector=sectorno\"");
    return -EINVAL;
  }
  if (sector>=72) {
    AliError(Form("invalid sector number %d", sector));
    return -EINVAL;
  }

  fCurrentSector=sector;
  fCurrentRow=row;

  return 0;
}

int AliHLTTPCClusterAccessHLTOUT::ProcessClusters(const char* params)
{
  /// process the cluster data from HLTOUT and fill array
  /// the cluster data can be in many different formats, e.g.
  /// raw or compressed
  int iResult=0;
  iResult = ScanParameters(params);
  if (iResult<0) return iResult;

  if (!fClusters) {
    fClusters=new AliRawClusterContainer;
  }
  if (!fClusters) return -ENOMEM;

  if (fClusters->HaveData()) {
    // cluster container already filled
//     TObjArray* pArray=fClusters->GetSectorArray(fCurrentSector, fPropagateSplitClusterFlag, fPropagateEdgeClusterFlag, fMarkEdgeClusters);
//     if (!pArray) {
//       AliError(Form("can not get cluster array for sector %d", sector));
//       return -ENOBUFS;
//     }
//     if (fVerbosity>0) AliInfo(Form("converted %d cluster(s) for sector %d", pArray->GetEntriesFast() ,sector));
    return 0; //pArray->GetEntriesFast();
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

  if (!fpDecoder) {
    fpDecoder=new AliHLTTPCDataCompressionDecoder;
    
    //Also load config objects together with fpDecoder
    TString cdbPath("HLT/ConfigTPC/TPCClusterAccessHLTOUT");
    TObject* pConf=AliHLTMisc::Instance().ExtractObject(AliHLTMisc::Instance().LoadOCDBEntry(cdbPath));
    if (pConf)
    {
      TObjString* pStr = dynamic_cast<TObjString*>(pConf);
      if (pStr)
      {
        if (pStr->GetString().Contains("-propagate-split-cluster-flag")) fPropagateSplitClusterFlag = 1;
        if (pStr->GetString().Contains("-propagate-edge-cluster-flag")) fPropagateEdgeClusterFlag = 1;
        if (pStr->GetString().Contains("-mark-edge-clusters")) fMarkEdgeClusters = 1;
      }
    }
  }

  if (!fpDecoder) {
    AliError("failed to create decoder instance");
    return -ENODEV;
  }

  AliHLTTPCDataCompressionDecoder& decoder=*fpDecoder;
  decoder.Clear();
  decoder.SetVerbosity(fVerbosity);

  bool bHavePartitionRawData=false;
  bool bHavePartitionCompressedData=false;

  bool bNextBlock=false;
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
    if (desc.fDataType==AliHLTTPCDefinitions::RawClustersDescriptorDataType()) {
      // header      
      if ((iResult=decoder.AddRawClustersDescriptor(&desc))<0) {
	return iResult;
      }
      bHavePartitionRawData = kTRUE;
    }
    //CompressionDescriptor should have priority over rawcluster descriptor in case both are present, because this describes the actual compressed data.
    if (desc.fDataType==AliHLTTPCDefinitions::DataCompressionDescriptorDataType()) {
      // header      
      if ((iResult=decoder.AddCompressionDescriptor(&desc))<0) {
	return iResult;
      }
      bHavePartitionCompressedData = kTRUE;
    }
    if (desc.fDataType==AliHLTTPCDefinitions::AliHLTDataTypeClusterMCInfo()) {
      // add mc information
      if ((iResult=decoder.AddClusterMCData(&desc))<0) {
	return iResult;
      }
    }
    if (desc.fDataType==AliHLTTPCDefinitions::RemainingClusterIdsDataType() ||
	desc.fDataType==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
      // add cluster ids
      if ((iResult=decoder.AddClusterIds(&desc))<0) {
	return iResult;
      }
    }
    if (desc.fDataType==AliHLTTPCDefinitions::ClustersFlagsDataType()) {
      // add cluster flag information
      if ((iResult=decoder.AddClusterFlags(&desc))<0) {
	return iResult;
      }
    }
  }

  vector<bool> bHavePartitionData(216, false);

  // read data
  iResult=-ENODATA;
  int nExtractedClusters=0;
  for (bNextBlock=(pHLTOUT->SelectFirstDataBlock()>=0);
       bNextBlock; bNextBlock=(pHLTOUT->SelectNextDataBlock()>=0)) {
    decoder.SetPadShift(0.0);
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
	(desc.fDataType==AliHLTTPCDefinitions::RawClustersDataType())) {
      // This is a special handling of data blocks produced with v5-01-Release
      // The pad shift by 0.5 was not included in the data but was applied in the
      // unpacking in this class. Changed in r51306, the next tag containing this
      // change in the online system is v5-01-Rev-07. There are only very few runs
      // of Sep 2011 with recorded clusters not containing the 0.5 shift
      // There was also a change in the data type of the compressed partition
      // cluster blocks which helps to identify the blocks which need the pad shift
      // here
      if (desc.fSize<sizeof(AliHLTTPCRawClusterData)) continue;
      const AliHLTTPCRawClusterData* clusterData = reinterpret_cast<const AliHLTTPCRawClusterData*>(buffer);
      if (!clusterData) continue;
      if (clusterData->fVersion==1) {
	// compressed clusters without the pad shift
	// no raw clusters (version==0) have ever been recorded
	decoder.SetPadShift(0.5);
      }
      AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(desc.fSpecification);
      AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(desc.fSpecification);
      if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr(desc.fSpecification) ||
	  partition!=AliHLTTPCDefinitions::GetMaxPatchNr(desc.fSpecification)) {
	AliFatal(Form("inconsistent cluster data: can not handle blocks containing multiple partitions, "
		      "block specification 0x%08x", desc.fSpecification));
      }
      iResult=decoder.ReadClustersPartition(fClusters->BeginRemainingClusterBlock(0, desc.fSpecification),
					    reinterpret_cast<AliHLTUInt8_t*>(desc.fPtr),
					    desc.fSize,
					    desc.fSpecification);
      if (iResult>=0) nExtractedClusters+=iResult;
      else {
	AliFatal(Form("processing of cluster block 0x%08x failed with error code %d", desc.fSpecification, iResult));
      }
      unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
      if (index>=bHavePartitionData.size()) bHavePartitionData.resize(index, false);
      if (bHavePartitionData[index]) {
	AliFatal(Form("inconsistent cluster data: multiple data blocks of identical specification indicate a failure "
		      "in the production of the data. Probably an HLT emulation chain is executed in the reconstruction "
		      "and produces data in addition to HLTOUT. Option 'ignore-hltout' is required in that case; "
		      "block specification 0x%08x", desc.fSpecification));
      }
      bHavePartitionData[index]=true;
      if (bHavePartitionCompressedData) {
	AliFatal(Form("inconsistent cluster data: both compressed and raw cluster blocks present in HLTOUT, indicates a failure "
		      "in the production of the data. Probably an HLT emulation chain is executed in the reconstruction "
		      "and produces data in addition to HLTOUT. Option 'ignore-hltout' is required in that case; "
		      "block specification 0x%08x", desc.fSpecification));
      }
      bHavePartitionRawData=true;
      continue;
    } else if (!TestBit(kSkipPartitionClusters) &&
	       (desc.fDataType==AliHLTTPCDefinitions::RemainingClustersCompressedDataType())) {
      AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(desc.fSpecification);
      AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(desc.fSpecification);
      if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr(desc.fSpecification) ||
	  partition!=AliHLTTPCDefinitions::GetMaxPatchNr(desc.fSpecification)) {
	AliFatal(Form("inconsistent cluster data: can not handle blocks containing multiple partitions, "
		      "block specification 0x%08x", desc.fSpecification));
      }
      iResult=decoder.ReadClustersPartition(fClusters->BeginRemainingClusterBlock(0, desc.fSpecification),
					    reinterpret_cast<AliHLTUInt8_t*>(desc.fPtr),
					    desc.fSize,
					    desc.fSpecification);
      if (iResult>0) nExtractedClusters+=iResult;
      unsigned index=slice*AliHLTTPCGeometry::GetNumberOfPatches()+partition;
      if (index>=bHavePartitionData.size()) bHavePartitionData.resize(index, false);
      if (bHavePartitionData[index]) {
	AliFatal(Form("inconsistent cluster data: multiple data blocks of identical specification indicate a failure "
		      "in the production of the data. Probably an HLT emulation chain is executed in the reconstruction "
		      "and produces data in addition to HLTOUT. Option 'ignore-hltout' is required in that case; "
		      "block specification 0x%08x", desc.fSpecification));
      }
      bHavePartitionData[index]=true;
      bHavePartitionData[index]=true;
      if (bHavePartitionRawData) {
	AliFatal(Form("inconsistent cluster data: both compressed and raw cluster blocks present in HLTOUT, indicates a failure "
		      "in the production of the data. Probably an HLT emulation chain is executed in the reconstruction "
		      "and produces data in addition to HLTOUT. Option 'ignore-hltout' is required in that case; "
		      "block specification 0x%08x", desc.fSpecification));
      }
      bHavePartitionCompressedData=true;
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
//   if (fVerbosity>0) {
//     int nConvertedClusters=0;
//     for (int s=0; s<72; s++) {
//       TObjArray* pArray=fClusters->GetSectorArray(s, fPropagateSplitClusterFlag, fPropagateEdgeClusterFlag, fMarkEdgeClusters);
//       if (!pArray) continue;
//       nConvertedClusters+=pArray->GetEntriesFast();
//     }
//     AliInfo(Form("extracted HLT clusters: %d, converted HLT clusters: %d", nExtractedClusters, nConvertedClusters));
//   }

  fClusters->MarkValid();
//   TObjArray* pArray=fClusters->GetSectorArray(fCurrentSector, fPropagateSplitClusterFlag, fPropagateEdgeClusterFlag, fMarkEdgeClusters);
//   if (!pArray) {
//     AliError(Form("can not get cluster array for sector %d", sector));
//     return -ENOBUFS;
//   }
//   if (fVerbosity>0) AliInfo(Form("converted %d cluster(s) for sector %d", pArray->GetEntriesFast() ,sector));
  return 0; //pArray->GetEntriesFast();
}

AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::AliRawClusterContainer()
  : fClusterMaps()
  , fSectorArray(new TClonesArray(AliTPCclusterMI::Class()))
  , fIterator()
  , fHaveData(false)
{
  /// constructor
  for (int i=0; i<72; i++) {
    fClusterMaps.push_back(new AliRawClusterEntryVector);
    fClusterMaps.back()->reserve(30000);
  }
}

AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::~AliRawClusterContainer()
{
  /// dectructor
  {
    for (vector<AliRawClusterEntryVector*>::iterator i=fClusterMaps.begin(); i!=fClusterMaps.end(); i++) {
      if (*i) {
	delete *i;
      }
    }
  }
  if (fSectorArray) {
    fSectorArray->Clear();
    delete fSectorArray;
    fSectorArray=NULL;
  }
}

AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::iterator& AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::BeginPartitionClusterBlock(int count, AliHLTUInt32_t specification)
{
  /// iterator of remaining clusters block of specification

  // reserve space in the array of all clusters
  // reserve space in the map of the partition
  unsigned index=AliHLTTPCDefinitions::GetMinSliceNr(specification);
  AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(specification);
  if (partition>=2) index+=36;
  if (index<fClusterMaps.size() &&
      fClusterMaps[index]!=NULL &&
      fClusterMaps[index]->size()+count>fClusterMaps[index]->capacity()) {
    fClusterMaps[index]->reserve(fClusterMaps[index]->size()+count);
  }

  fIterator=iterator(this);
  return fIterator;
}

AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::iterator& AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::BeginTrackModelClusterBlock(int /*count*/)
{
  /// iterator of track model clusters
  fIterator=iterator(this);
  return fIterator;
}

AliHLTTPCClusterAccessHLTOUT::AliRawClusterEntry* AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::NextCluster(int slice, int partition)
{
  /// load next cluster from array of the sepcific sector
  unsigned sector=partition<2?slice:slice+36;
  if (fClusterMaps.size()<=sector || 
      fClusterMaps[sector]==NULL) {
    AliErrorClass(Form("no cluster array available for sector %d", sector));
    return NULL;
  }
  AliRawClusterEntryVector& map=*(fClusterMaps[sector]);
  map.push_back(AliRawClusterEntry());
  return &map.back();
}

void  AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::Clear(Option_t* option)
{
  /// internal cleanup
  if (strcmp(option, "event")==0) {
    for (vector<AliRawClusterEntryVector*>::iterator i=fClusterMaps.begin(); i!=fClusterMaps.end(); i++)
      if (*i) (*i)->clear();
    if (fSectorArray) fSectorArray->Clear("C");
    fHaveData=false;
  }
  if (strcmp(option, "sector")==0) {
    if (fSectorArray) fSectorArray->Clear("C");
  }
}

TObjArray* AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::GetSectorArray(unsigned sector, int propagateSplitClusterFlag, int propagateEdgeClusterFlag, int markEdgeClusters) const
{
  /// get the cluster array for a sector
  if (fClusterMaps.size()<=sector) return NULL;
  if (fSectorArray &&
      FillSectorArray(fSectorArray, sector, -1, propagateSplitClusterFlag, propagateEdgeClusterFlag, markEdgeClusters)<0) {
    fSectorArray->Clear();
  }
  return fSectorArray;
}

int AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::FillSectorArray(TClonesArray* pSectorArray, unsigned sector, int row, int propagateSplitClusterFlag, int propagateEdgeClusterFlag, int markEdgeClusters) const
{
  /// fill the cluster array for a sector and specific row if specified
  if (!pSectorArray) return -EINVAL;
  if (fClusterMaps.size()<=sector) return -ERANGE;
  pSectorArray->Clear();

  AliRawClusterEntryVector& map=*fClusterMaps[sector];
  unsigned nFilled=0;
  for (unsigned i=0; i<map.size(); i++) {    
    if (row>=0 && map[i].fCluster.GetPadRow()!=row) continue;
    AliTPCclusterMI* pCluster=new ((*pSectorArray)[nFilled]) AliTPCclusterMI;
    if (!pCluster) break;
    
    pCluster->SetRow(map[i].fCluster.GetPadRow());
    pCluster->SetPad(map[i].fCluster.GetPad());
    pCluster->SetTimeBin(map[i].fCluster.GetTime());
    pCluster->SetSigmaY2(map[i].fCluster.GetSigmaPad2());
    pCluster->SetSigmaZ2(map[i].fCluster.GetSigmaTime2());
    pCluster->SetQ(map[i].fCluster.GetCharge());
    pCluster->SetMax(map[i].fCluster.GetQMax());
    if (propagateSplitClusterFlag)
    {
      pCluster->SetType((int) map[i].fCluster.GetFlagSplitPad() + ((int) map[i].fCluster.GetFlagSplitTime() << 1));
    }
    else
    {
      pCluster->SetType(0);
    }
    
    if (propagateEdgeClusterFlag && map[i].fCluster.GetFlagEdge())
    {
        pCluster->SetType(-(pCluster->GetType() + 3));
    }
    else if (markEdgeClusters)
    {
      int fullRow = pCluster->GetRow();
      if (sector > 36) fullRow += AliHLTTPCGeometry::GetFirstRow(2);
      int maxPad = AliHLTTPCGeometry::GetNPads(fullRow);
      if (pCluster->GetPad() < 1.2 + pCluster->GetSigmaY2() || pCluster->GetPad() > maxPad - 1.2 - pCluster->GetSigmaY2())
      {
        //printf("Cluster Sector %d Row %d Pad %f Time %f Sigma %f MaxPad %d\n", sector, (int) pCluster->GetRow(), pCluster->GetPad(), pCluster->GetTimeBin(), (float) pCluster->GetSigmaY2(), maxPad);
        pCluster->SetType(-(pCluster->GetType() + 3));
      }
    }

    for (int k=0; k<3; k++) {
      // TODO: sort the labels according to the weight in order to assign the most likely mc label
      // to the first component 
      pCluster->SetLabel(map[i].fMC.fClusterID[k].fMCID, k);    
    }
    nFilled++;
  }

  return 0;
}

void AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::Print(Option_t *option) const
{
  /// inherited from TObject
  cout << "AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer" << endl;
  ios::fmtflags coutflags=cout.flags(); // backup cout status flags
  bool bAll=false;
  if ((bAll=(strcmp(option, "full")==0)) ||
      strcmp(option, "short")==0) {
    for (unsigned iArray=0; iArray<fClusterMaps.size(); iArray++) {
      if (fClusterMaps[iArray]) {
	AliRawClusterEntryVector& map=*fClusterMaps[iArray];
	cout << "  sector " << setfill(' ') << setw(2) << iArray << ": " << map.size() << endl;
	if (bAll) {
	  for (unsigned iCluster=0; iCluster<map.size(); iCluster++) {
	    AliHLTTPCRawCluster &cluster = map[iCluster].fCluster;
	    cout << "    AliTPCclusterMI:"
		 << "  row="    << cluster.GetPadRow() 
		 << "  pad="    << cluster.GetPad()
		 << "  time="   << cluster.GetTime()
		 << "  charge=" << cluster.GetCharge()
		 << "  maxq="   << cluster.GetQMax()
		 << endl;
	  }
	}
      }
    }
  }
  cout.flags(coutflags); // restore the original flags
}

AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::iterator& AliHLTTPCClusterAccessHLTOUT::AliRawClusterContainer::iterator::Next(int slice, int partition)
{
  // switch to next cluster
  if (!fData) {
    fEntry=NULL;
    return *this;
  }
  if (fClusterNo>=0 && !fEntry) {
    // end was reached before
    return *this;
  }
  fEntry=fData->NextCluster(slice, partition);

  // offline uses row number in physical sector, inner sector consists of
  // partitions 0 and 1, outer sector of partition 2-5
  fRowOffset=partition<2?0:AliHLTTPCGeometry::GetFirstRow(2);
  return *this;
}
