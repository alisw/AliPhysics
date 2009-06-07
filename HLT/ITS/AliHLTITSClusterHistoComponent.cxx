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

/** @file   AliHLTITSClusterHistoComponent.cxx
    @author Gaute Ovrebekk
    @brief  Component for ploting charge in clusters
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTITSClusterHistoComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliITSRecPoint.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"

//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSClusterHistoComponent)

AliHLTITSClusterHistoComponent::AliHLTITSClusterHistoComponent()
:
fXY(NULL),                     
  fXYZ(NULL),                   
  fCharge(NULL),   
  fPlotCharge(kFALSE),   
  fPlotXY(kTRUE),
  fPlotXYZ(kFALSE) 
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTITSClusterHistoComponent::~AliHLTITSClusterHistoComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTITSClusterHistoComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "ITSClusterHisto";
}

void AliHLTITSClusterHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeTObjArray );
}

AliHLTComponentDataType AliHLTITSClusterHistoComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram;

}

void AliHLTITSClusterHistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 10;
}

AliHLTComponent* AliHLTITSClusterHistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTITSClusterHistoComponent;
}

int AliHLTITSClusterHistoComponent::DoInit( int argc, const char** argv )
{
  fPlotCharge=kFALSE;   
  fPlotXY=kTRUE;
  fPlotXYZ=kFALSE; 
     
  if(fPlotCharge){fCharge = new TH1F("fCharge","Total Charge of clusters",4000,0,4000);}
  if(fPlotXY){fXY = new TH2F("fXY","Global XY of ITS clusters",1600,-80,80,1600,-80,80);}
  if(fPlotXYZ){fXYZ = new TH3F("fXYZ","Global XYZ of ITS clusters",1600,-80,80,1600,-80,80,2000,-100,100);}
  
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
  
int AliHLTITSClusterHistoComponent::DoDeinit()
{
  // see header file for class documentation
  if(fCharge!=NULL) delete fCharge;
  if(fXY!=NULL) delete fXY;     
  if(fXYZ!=NULL) delete fXYZ;
  return 0;
}

int AliHLTITSClusterHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  
  int TotalSpacePoint = 0;
  
  const AliHLTComponentBlockData* iter = NULL;
  
  if(!IsDataEvent())
    return 0;
  
  for ( iter = GetFirstInputBlock(kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD); iter != NULL; iter = GetNextInputBlock() ) {
    
    const AliHLTITSClusterData* clusterData = (const AliHLTITSClusterData*) iter->fPtr;
    Int_t nSpacepoint = (Int_t) clusterData->fSpacePointCnt;
    TotalSpacePoint += nSpacepoint;
    AliHLTITSSpacePointData *clusters = (AliHLTITSSpacePointData*) clusterData->fSpacePoints;

    for(int i=0;i<nSpacepoint;i++){
      Int_t lab[4]={0,0,0,0};
      Float_t hit[6]={0,0,0,0,0,0};
      Int_t info[3]={0,0,0};
      
      lab[0]=clusters[i].fTracks[0];
      lab[1]=clusters[i].fTracks[1];
      lab[2]=clusters[i].fTracks[2];
      lab[3]=clusters[i].fIndex;
      hit[0]=clusters[i].fY;
      hit[1]=clusters[i].fZ;
      hit[2]=clusters[i].fSigmaY2;
      hit[3]=clusters[i].fSigmaZ2;
      hit[4]=clusters[i].fQ;
      hit[5]=clusters[i].fSigmaYZ;
      info[0]=clusters[i].fNy;
      info[1]=clusters[i].fNz;
      info[2]=clusters[i].fLayer;

      Float_t xyz[3];
      AliITSRecPoint recpoint(lab,hit,info);
      recpoint.GetGlobalXYZ(xyz);
      if(fPlotXY){
	fXY->Fill(xyz[0],xyz[1]);
      }
      if(fPlotXYZ){
	fXYZ->Fill(xyz[0],xyz[1],xyz[2]);
      }
      if(fPlotCharge){
	fCharge->Fill(recpoint.GetQ());
      }
    }
  }
  
  if(fPlotCharge){
    AliHLTUInt32_t fSpecification = 0x0;
    PushBack( (TObject*) fCharge,kAliHLTDataTypeHistogram,fSpecification);
  }
  if(fPlotXY){
    AliHLTUInt32_t fSpecification = 0x0;
    PushBack( (TObject*) fXY,kAliHLTDataTypeHistogram,fSpecification);
  }
  if(fPlotXYZ){
    AliHLTUInt32_t fSpecification = 0x0;
    PushBack( (TObject*) fXYZ,kAliHLTDataTypeHistogram,fSpecification);
  }
  
  HLTInfo("ITSClusterHisto found %d Total Spacepoints", TotalSpacePoint);
  
  return 0;
}

int AliHLTITSClusterHistoComponent::Configure(const char* arguments)
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
	HLTInfo("Ploting all historgams");
	fPlotXY = kTRUE;
	fPlotXYZ = kTRUE;
	fPlotCharge = kTRUE;
	continue;
      }
      
      else if (argument.CompareTo("-plot-xy")==0) {
	HLTInfo("Ploting Global XY");
	fPlotXY = kTRUE;
	continue;
      }

      else if (argument.CompareTo("-plot-xyz")==0) {
	HLTInfo("Ploting Global XYZ");
	//fPlotXYZ = kTRUE;
	continue;
      }
      else if (argument.CompareTo("-plot-charge")==0) {
	HLTInfo("Ploting charge of clusters");
	fPlotCharge = kTRUE;
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
  
  if(!fCharge && fPlotCharge){fCharge = new TH1F("fCharge","Total Charge of clusters",4000,0,4000);}
  if(!fXY && fPlotXY){fXY = new TH2F("fXY","Global XY of ITS clusters",1600,-80,80,1600,-80,80);}
  if(!fXYZ && fPlotXYZ){fXYZ = new TH3F("fXYZ","Global XYZ of ITS clusters",1600,-80,80,1600,-80,80,2000,-100,100);}
  
  return iResult;
}

int AliHLTITSClusterHistoComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  
  const char* path="HLT/ConfigITS/HistoComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path);
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
