// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk <ovrebekk@ift.uib.no>                  *
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

/** @file   AliHLTTPCClusterHistoComponent.cxx
    @author Gaute Ovrebekk
    @brief  Component for ploting charge in clusters
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCClusterHistoComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDefinitions.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"

//#include "AliHLTTPC.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCClusterHistoComponent)

AliHLTTPCClusterHistoComponent::AliHLTTPCClusterHistoComponent()
:
fTotalClusterChargeOROCAll(NULL),                     
  fTotalClusterChargeIROCAll(NULL),                   
  fTotalClusterChargeROCSelection(NULL),               
  fTotalClusterChargePartitionSelection(NULL),         
  fQMaxPartitionAll(NULL),                             
  fQMaxROCAll(NULL),                              
  fNumberOfClusters(NULL),                        
  fPlotChargeOROCAll(kTRUE),   
  fPlotChargeIROCAll(kTRUE),
  fPlotChargeROCSel(kFALSE), 
  fPlotChargePartSel(kFALSE),
  fPlotQmaxPartAll(kTRUE),  
  fPlotQmaxROCAll(kTRUE),   
  fPlotNClusters(kTRUE)    
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTTPCClusterHistoComponent::~AliHLTTPCClusterHistoComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCClusterHistoComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "TPCClusterHisto";
}

void AliHLTTPCClusterHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
}

AliHLTComponentDataType AliHLTTPCClusterHistoComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram;

}

void AliHLTTPCClusterHistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTTPCClusterHistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCClusterHistoComponent;
}

int AliHLTTPCClusterHistoComponent::DoInit( int argc, const char** argv )
{
  //  HLTFatal("Initializing with log fatal");
  //  cout<<"Initializing with cout"<<endl;
  
  fPlotChargeOROCAll=kTRUE;   
  fPlotChargeIROCAll=kTRUE;
  fPlotChargeROCSel=kFALSE; 
  fPlotChargePartSel=kFALSE;
  fPlotQmaxPartAll=kTRUE;  
  fPlotQmaxROCAll=kTRUE;   
  fPlotNClusters=kTRUE;
   
  if(fPlotChargeOROCAll){fTotalClusterChargeOROCAll = new TH1F("fTotalClusterChargeOROCAll","Total Charge of clusters in all OROC",4000,0,4000);}
  if(fPlotChargeIROCAll){fTotalClusterChargeIROCAll = new TH1F("fTotalClusterChargeIROCAll","Total Charge of clusters in all IROC",4000,0,4000);}
  if(fPlotChargeROCSel){fTotalClusterChargeROCSelection = new TH1F("fTotalClusterChargeROCSelection","Total Charge of clusters in selection ROC",4000,0,4000);}
  if(fPlotChargePartSel){fTotalClusterChargePartitionSelection = new TH1F("fTotalClusterChargePartitionSelection","Total Charge of clusters in sel Part",4000,0,4000);}
  if(fPlotQmaxPartAll){fQMaxPartitionAll = new TH1F("fQMaxPartitionAll","QMax for All Partitions",216,0,216);}
  if(fPlotQmaxROCAll){fQMaxROCAll = new TH1F("fQMaxROCAll","QMax for All Partitions",72,0,72);}
  if(fPlotNClusters){fNumberOfClusters = new TH1F("fNumberOfClusters","Total Number of Clusters",1,0,1);}

  int iResult=0;
  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }
  
  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  }  

  return iResult; 
}
  
int AliHLTTPCClusterHistoComponent::DoDeinit()
{
  // see header file for class documentation
  if(fTotalClusterChargeOROCAll!=NULL) delete fTotalClusterChargeOROCAll;
  if(fTotalClusterChargeIROCAll!=NULL) delete fTotalClusterChargeIROCAll;     
  if(fQMaxPartitionAll!=NULL) delete fQMaxPartitionAll;
  if(fQMaxROCAll!=NULL) delete fQMaxROCAll;
  return 0;
}

int AliHLTTPCClusterHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  
  int TotalSpacePoint = 0;
  
  const AliHLTComponentBlockData* iter = NULL;
  
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;
  
  fQMaxPartitionAll->Reset();
  fQMaxROCAll->Reset();

  for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock() ) {
    
    Int_t thisrow=-1,thissector=-1,row=-1;
    
    AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
    AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
    row = AliHLTTPCTransform::GetFirstRow(patch); 
    AliHLTTPCTransform::Slice2Sector(slice,row,thissector,thisrow);
    
    HLTDebug ( "Input Data - TPC cluster - Slice/Patch: %d/%d.", slice, patch );
    
    const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*) iter->fPtr;
    Int_t nSpacepoint = (Int_t) clusterData->fSpacePointCnt;
    TotalSpacePoint += nSpacepoint;
    //HLTInfo("KryptonHisto found %d Spacepoints in slice %d patch %d", nSpacepoint, slice, patch);
    AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*) &clusterData->fSpacePoints;
    
    UInt_t tmpQPart = 0;//,tmpQROC = -1;
    
    for(int i=0;i<nSpacepoint;i++){
      if(fPlotChargeOROCAll){
	if(thissector>=36){
	  fTotalClusterChargeOROCAll->Fill(clusters[i].fCharge);
	}
      }
      if(fPlotChargeIROCAll){
	if(thissector<=35){
	  fTotalClusterChargeIROCAll->Fill(clusters[i].fCharge);
	}
      }
      if(fPlotChargeROCSel){
	
      }
      if(fPlotChargePartSel){
	
      }
      if(fPlotQmaxPartAll){
	if(clusters[i].fQMax>tmpQPart){
	  fQMaxPartitionAll->SetBinContent(patch+6*slice,clusters[i].fQMax);
	  tmpQPart=clusters[i].fQMax;
	}
      }
      if(fPlotQmaxROCAll){
	if(clusters[i].fQMax>fQMaxROCAll->GetBinContent(thissector)){
	  fQMaxROCAll->SetBinContent(thissector,clusters[i].fQMax);
	  //	  tmpQROC=clusters[i].fQMax;
	}
      }
    }
    if(fPlotNClusters){
      fNumberOfClusters->Fill(nSpacepoint);
    }
  }
  
  //delete til dodeinit
  if(fPlotChargeOROCAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,2,5);
    PushBack( (TObject*) fTotalClusterChargeOROCAll,kAliHLTDataTypeHistogram,fSpecification);
  }
  if(fPlotChargeIROCAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,1);
    PushBack( (TObject*) fTotalClusterChargeIROCAll,kAliHLTDataTypeHistogram,fSpecification);
  }
  if(fPlotChargeROCSel){
    
    
  }
  if(fPlotChargePartSel){
    
    
  }
  if(fPlotQmaxPartAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fQMaxPartitionAll,kAliHLTDataTypeHistogram,fSpecification);
    //delete fQMaxPartitionAll;
  }
  if(fPlotQmaxROCAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fQMaxROCAll,kAliHLTDataTypeHistogram,fSpecification);
    //delete fQMaxROCAll;
  }
  if(fPlotNClusters){
    
    
  }
  
  HLTInfo("KryptonHisto found %d Total Spacepoints", TotalSpacePoint);
  
  return 0;
}

int AliHLTTPCClusterHistoComponent::Configure(const char* arguments)
{
  
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  
  TObjArray* pTokens=allArgs.Tokenize(" ");
  
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
      
      if (argument.CompareTo("-plot-all")==0) {
	HLTInfo("Ploting charge of all clusters");
	//fPlotAll = kTRUE;
	continue;
      }
      
      else if (argument.CompareTo("-plot-trackclusters")==0) {
	HLTInfo("Ploting charge of clusters used on a track");
	//fPlotAll = kFALSE;
	continue;
      }
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  
  //if hvis det eksisterer
  if(fPlotChargeOROCAll){fTotalClusterChargeOROCAll = new TH1F("fTotalClusterChargeOROCAll","Total Charge of clusters in all OROC",4000,0,4000);}
  if(fPlotChargeIROCAll){fTotalClusterChargeIROCAll = new TH1F("fTotalClusterChargeIROCAll","Total Charge of clusters in all IROC",4000,0,4000);}
  if(fPlotChargeROCSel){fTotalClusterChargeROCSelection = new TH1F("fTotalClusterChargeROCSelection","Total Charge of clusters in selection ROC",4000,0,4000);}
  if(fPlotChargePartSel){fTotalClusterChargePartitionSelection = new TH1F("fTotalClusterChargePartitionSelection","Total Charge of clusters in sel Part",4000,0,4000);}
  if(fPlotQmaxPartAll){fQMaxPartitionAll = new TH1F("fQMaxPartitionAll","QMax for All Partitions",216,0,216);}
  if(fPlotQmaxROCAll){fQMaxROCAll = new TH1F("fQMaxROCAll","QMax for All Partitions",72,0,72);}
  if(fPlotNClusters){fNumberOfClusters = new TH1F("fNumberOfClusters","Total Number of Clusters",100,0,100);}

  return iResult;
}

int AliHLTTPCClusterHistoComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path="HLT/ConfigTPC/KryptonHistoComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	iResult=Configure(pString->GetString().Data());
      } else {
	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }

  return iResult;
}
