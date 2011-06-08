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
#include "AliLog.h"
#include "AliHLTSystem.h"
#include "AliHLTPluginBase.h"
#include "AliTPCclusterMI.h"
#include "TClonesArray.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterAccessHLTOUT)

AliHLTTPCClusterAccessHLTOUT::AliHLTTPCClusterAccessHLTOUT()
  : TObject()
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

void AliHLTTPCClusterAccessHLTOUT::Execute(const char *method,  const char */*params*/, Int_t *error)
{
  /// inherited from TObject: abstract command interface
  if (strcmp(method, "read")==0) {
    int iResult=ProcessClusters();
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

int AliHLTTPCClusterAccessHLTOUT::ProcessClusters()
{
  /// process the cluster data from HLTOUT and fill array
  /// the cluster data can be in many different formats, e.g.
  /// raw or compressed
  int iResult=0;
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

  if (pHLTOUT->SelectFirstDataBlock(AliHLTTPCDefinitions::fgkClustersDataType)>=0) {
    iResult=ReadAliHLTTPCClusterData(pHLTOUT, fClusters);
  }

  pSystem->ReleaseHLTOUT(pHLTOUT);
  return iResult;
}

int AliHLTTPCClusterAccessHLTOUT::ReadAliHLTTPCClusterData(AliHLTOUT* pHLTOUT, TClonesArray* pClusters) const
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
    const AliHLTTPCClusterData* clusterData = reinterpret_cast<const AliHLTTPCClusterData*>(pBuffer);
    Int_t nSpacepoints = (Int_t) clusterData->fSpacePointCnt;
    if (nSpacepoints*sizeof(AliHLTTPCSpacePointData) + sizeof(AliHLTTPCClusterData) != size) {
      AliError("inconsistent cluster data block size, skipping block");
      continue;
    }
    const AliHLTTPCSpacePointData *clusters = clusterData->fSpacePoints;
    int offset=fClusters->GetEntries();
    fClusters->ExpandCreate(offset+nSpacepoints);
    for (int i=0; i<nSpacepoints; i++) {
      if (!fClusters->At(offset+i)) continue;
      AliTPCclusterMI* pCluster=dynamic_cast<AliTPCclusterMI*>(fClusters->At(offset+i));
      if (!pCluster) {
	AliError("invalid object type, expecting AliTPCclusterMI");
	break; // this is a problem of all objects
      }
      pCluster->SetRow(clusters[i].fPadRow);
      pCluster->SetTimeBin(clusters[i].fZ);
      pCluster->SetPad(clusters[i].fX);
      pCluster->SetQ(clusters[i].fCharge);
      pCluster->SetMax(clusters[i].fQMax);
    }
    AliInfo(Form("converted %d cluster(s)", nSpacepoints));
  } while (pHLTOUT->SelectNextDataBlock()>=0);
  return iResult;
}
