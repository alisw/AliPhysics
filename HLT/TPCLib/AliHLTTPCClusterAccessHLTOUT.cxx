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
#include "AliHLTOUT.h"
#include "AliHLTComponent.h"
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
      if (pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo, spec)>=0) {
	iResult=ReadAliHLTTPCClusterMCData(pHLTOUT, tpcClusterLabels);
      }

      if (pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::fgkClustersDataType, spec)>=0) {
	iResult=ReadAliHLTTPCClusterData(pHLTOUT, fClusters);
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
    for (int i=0; i<nSpacepoints; i++) {
      if (!pClusters->At(offset+i)) continue;
      AliTPCclusterMI* pCluster=dynamic_cast<AliTPCclusterMI*>(pClusters->At(offset+i));
      if (!pCluster) {
	AliError("invalid object type, expecting AliTPCclusterMI");
	break; // this is a problem of all objects
      }
      pCluster->SetRow(clusters[i].fPadRow);
      pCluster->SetTimeBin(clusters[i].fZ);
      pCluster->SetPad(clusters[i].fY);
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
    if (fVerbosity>0) AliInfo(Form("converted %d cluster(s) from block 0x%08x", nSpacepoints, specification));
  } while (pHLTOUT->SelectNextDataBlock()>=0);
  return iResult;
}
