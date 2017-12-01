// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kenneth Aamodt <Kenneth.Aamodt@cern.ch>               *
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

/** @file   AliHLTTPCClusterDumpComponent.cxx
    @author Kenneth Aamodt
    @date   
    @brief  Special file writer converting TPC clusters input to readable
            ASCII format. 
*/

#include <cassert>
#include "AliHLTTPCClusterDumpComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCGeometry.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterDumpComponent)

  AliHLTTPCClusterDumpComponent::AliHLTTPCClusterDumpComponent()
    :
    AliHLTFileWriter()
//fSlice(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCClusterDumpComponent::~AliHLTTPCClusterDumpComponent()
{
  // see header file for class documentation
}

const char* AliHLTTPCClusterDumpComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCClusterDump";
}

void AliHLTTPCClusterDumpComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
}

AliHLTComponent* AliHLTTPCClusterDumpComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCClusterDumpComponent;
}

int AliHLTTPCClusterDumpComponent::InitWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCClusterDumpComponent::ScanArgument(int /*argc*/, const char** /*argv*/)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  bool bMissingParam=0;
  int i=0;
  
  if (bMissingParam) iResult=-EPROTO;
  else if (iResult>=0) iResult=i;
  
  return iResult;
}

int AliHLTTPCClusterDumpComponent::CloseWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCClusterDumpComponent::DumpEvent( const AliHLTComponentEventData& evtData,
					      AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation

  HLTDebug("Entering DumpEvent");

  int iResult=0;
  int blockno=0;
  const AliHLTComponentBlockData* pDesc=NULL;

  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  Int_t spacePointCounter=0;
  
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    TString filename;
    iResult=BuildFileName(evtData.fEventID, 0, pDesc->fDataType, 0, filename);
    ios::openmode filemode=(ios::openmode)0;
    if (fCurrentFileName.CompareTo(filename)==0) {
      filemode=ios::app;
    } else {
      fCurrentFileName=filename;
    }
    
    if (iResult>=0) {
      ofstream dump(fCurrentFileName.Data(), filemode);
      if (dump.good()) {
	
	if(pDesc->fDataType!=AliHLTTPCDefinitions::fgkClustersDataType){continue;}
	
	//if (dump.good() || 1) {//the || 1 is there since dump.good() will return false( EOF )
	iResult=1;
	const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*) pDesc->fPtr;
	Int_t nSpacepoints = (Int_t) clusterData->fSpacePointCnt;
	AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*) &clusterData->fSpacePoints;
	
	for(int i=0;i<nSpacepoints;i++){
	  UInt_t idCluster = clusters[i].GetRawID();
	  Int_t slice = AliHLTTPCGeometry::CluID2Slice(idCluster);
	  Int_t patch = AliHLTTPCGeometry::CluID2Partition(idCluster);
	  
	  dump << "" << endl;
	  dump << "ClusterNumber: " << spacePointCounter << endl;
	  dump << "Slice:         " << slice << endl;
	  dump << "Partition:     " << patch << endl;
	  dump << "[X,Y,Z]:       [" << clusters[i].fX<<" , "<<clusters[i].fY<<" , "<<clusters[i].fZ <<"]"<< endl;
	  //Float_t xyz[3]={clusters[i].fX,clusters[i].fY,clusters[i].fZ};
	  //AliHLTTPCGeometry::LocHLT2Raw(xyz,(Int_t)(clusters[i].fID/10),(Int_t)(clusters[i].fID%10));
	  //dump << "[R,P,T]:       [" << xyz[0]<<" , "<<xyz[1]<<" , "<<xyz[2] <<"]"<< endl;
	  dump << "Total Charge:  " << clusters[i].fCharge         << endl;
	  dump << "Q Max:         " << clusters[i].fQMax           << endl;
	  spacePointCounter++;
	}
      }
      else {
	HLTError("can not open file %s for writing", fCurrentFileName.Data());
	iResult=-EBADF;
      }
      dump.close();
    }
  }
    
  return iResult;
}
