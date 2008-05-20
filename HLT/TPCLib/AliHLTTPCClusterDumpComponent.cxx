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

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterDumpComponent)

  AliHLTTPCClusterDumpComponent::AliHLTTPCClusterDumpComponent()
    :
    AliHLTFileWriter(),
    fDirectory("")
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

int AliHLTTPCClusterDumpComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation

  TString argument="";
  bool bMissingParam=0;
  int i=0;
  do {
    if (i>=argc || (argument=argv[i]).IsNull()) continue;
    
    // -directory
    if (argument.CompareTo("-directory")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fDirectory=argv[i];
      break;
    }
  }while(0);

  HLTWarning("AliHLTTPCClusterDumpComponent does not have any arguments at this time");
  return 0;
}

int AliHLTTPCClusterDumpComponent::CloseWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCClusterDumpComponent::DumpEvent( const AliHLTComponentEventData& evtData,
					      const AliHLTComponentBlockData* /*blocks*/, 
					      AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation

  HLTDebug("Entering DumpEvent");

  int iResult=0;
  int blockno=0;
  const AliHLTComponentBlockData* pDesc=NULL;

  Int_t spacePointCounter=0;

  //building the filename
  fCurrentFileName="";
  ios::openmode filemode=(ios::openmode)0;
  if (!fDirectory.IsNull()) {
    fCurrentFileName+=fDirectory;
  }
  fCurrentFileName+="TPCClusterDump_Event";
  fCurrentFileName+=Form("_%d", GetEventCount());
  ofstream dump(fCurrentFileName.Data(), filemode);

  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    HLTDebug("event %Lu block %d: %s 0x%08x size %d", evtData.fEventID, blockno, DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification, pDesc->fSize);

    if(pDesc->fDataType!=AliHLTTPCDefinitions::fgkClustersDataType){continue;}
 
     if (dump.good()) {
       iResult=1;
       const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*) pDesc->fPtr;
       Int_t nSpacepoints = (Int_t) clusterData->fSpacePointCnt;
       AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*) &clusterData->fSpacePoints;
       
       for(int i=0;i<nSpacepoints;i++){
	 dump << "" << endl;
	 dump << "ClusterNumber: " << spacePointCounter << endl;
	 dump << "Slice:         " << clusters[i].fUsed    << endl;//quick fix to get the partiion and slice numbers to the clusterdump
	 dump << "Partition:     " << clusters[i].fTrackN  << endl;//quick fix to get the partiion and slice numbers to the clusterdump
	 dump << "X:             " << clusters[i].fX       << endl;
	 dump << "Y:             " << clusters[i].fY       << endl;
	 dump << "Z:             " << clusters[i].fZ       << endl;
	 dump << "ID:            " << clusters[i].fID      << endl;
	 dump << "Pad row:       " << clusters[i].fPadRow << endl;
	 dump << "fSigmaY2:      " << clusters[i].fSigmaY2 << endl;
	 dump << "fSigmaZ2:      " << clusters[i].fSigmaZ2 << endl;
	 dump << "Charge:        " << clusters[i].fCharge << endl;
	 dump << "Q Max:         " << clusters[i].fMaxQ << endl;
	 spacePointCounter++;
       }
       
     }
     else {
       HLTError("can not open file %s for writing", fCurrentFileName.Data());
       iResult=-EBADF;
     }

    dump.close();
  }
  return iResult;
}
