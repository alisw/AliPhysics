// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Piergiorgio Cerello <cerello@to.infn.it>              *
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

/// @file   AliHLTITSQHistoComponent.cxx
/// @author Piergiorgio Cerello cerello@to.infn.it
/// @date   2009-07-03
/// @brief  Interface component to the ITS QA

#include "AliHLTITSQAComponent.h"
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

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSQAComponent)

AliHLTITSQAComponent::AliHLTITSQAComponent()
:
fAliITSQADataMakerRec(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTITSQAComponent::~AliHLTITSQAComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTITSQAComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "ITSClusterQA";
}

void AliHLTITSQAComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeTObjArray );
}

AliHLTComponentDataType AliHLTITSQAComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram;

}

void AliHLTITSQAComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 10;
}

AliHLTComponent* AliHLTITSQAComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTITSQAComponent;
}

int AliHLTITSQAComponent::DoInit( int argc, const char** argv )
{
// add AliITSQADataMakerRec constructor
  fAliITSQADataMakerRec = new AliITSQADataMakerRec();
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
  
int AliHLTITSQAComponent::DoDeinit()
{
  // see header file for class documentation
// add AliITSQADataMakerRec destruction
  if(fAliITSQADataMakerRec) delete fAliITSQADataMakerRec;
  return 0;
}

int AliHLTITSQAComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  
  int TotalSpacePoint = 0;
  
  const AliHLTComponentBlockData* iter = NULL;
  
  if(!IsDataEvent())
    return 0;
  
  static Int_t rp = 1; 
  // Check id histograms already created for this Event Specie
  if ( rp ) { fAliITSQADataMakerRec->InitRecPoints(); rp = 0; }
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

      AliITSRecPoint recpoint(lab,hit,info);
      fAliITSQADataMakerRec->FillRecPoint(recpoint);
    }
  }
	for(Int_t i=0; i<6; i++) {
		//if(fPlotCharge){
		AliHLTUInt32_t specification = 0x0;
		PushBack( (TObject*) fAliITSQADataMakerRec->GetITSGlobalHisto(i),kAliHLTDataTypeHistogram,specification);
		//}
	}
	HLTInfo("ITSClusterHisto found %d Total Spacepoints", TotalSpacePoint);
  
  return 0;
}

int AliHLTITSQAComponent::Configure(const char* arguments)
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
/*      
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
*/
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  
  return iResult;
}

int AliHLTITSQAComponent::Reconfigure(const char* cdbEntry, const char* chainId)
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
