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
  , fpDecoder(NULL)
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

  if (!fpDecoder) {
    fpDecoder=new AliHLTTPCDataCompressionDecoder;
  }

  if (!fpDecoder) {
    AliError("failed to create decoder instance");
    return -ENODEV;
  }

  AliHLTTPCDataCompressionDecoder& decoder=*fpDecoder;
  decoder.Clear();
  decoder.SetVerbosity(fVerbosity);
  decoder.EnableClusterMerger();

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
  }

  bool bHavePartitionRawData=false;
  bool bHavePartitionCompressedData=false;
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
      // There was also a chenge in the data type of the compressed partition
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
      unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
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
      unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
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
    return *this;
  }
  if (fClusterNo>=0 && !fCluster) {
    // end was reached before
    return *this;
  }
  fCluster=fData->NextCluster(slice, partition);
  // offline uses row number in physical sector, inner sector consists of
  // partitions 0 and 1, outer sector of partition 2-5
  fRowOffset=partition<2?0:AliHLTTPCTransform::GetFirstRow(2);
  return *this;
}
