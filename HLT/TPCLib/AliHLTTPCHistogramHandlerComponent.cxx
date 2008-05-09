// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
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

/** @file   AliHLTTPCHistogramHandlerComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  The Histogram Handler component
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCHistogramHandlerComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <sys/time.h>
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TMath.h"


AliHLTTPCHistogramHandlerComponent gAliHLTTPCHistogramHandlerComponent;

ClassImp(AliHLTTPCHistogramHandlerComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCHistogramHandlerComponent::AliHLTTPCHistogramHandlerComponent()
    :    
    fNoiseHistograms(0),
    fKryptonHistograms(0),
    fSpecificationTPCA(0),
    fSpecificationTPCC(0),
    fSlice(-99),
    fHistTH1Tmp(NULL),
    fHistTH2Tmp(NULL),
    fHistTPCSideA(NULL),  
    fHistTPCSideC(NULL)    
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCHistogramHandlerComponent::~AliHLTTPCHistogramHandlerComponent() { 
// see header file for class documentation

}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCHistogramHandlerComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCHistogramHandler";
}

void AliHLTTPCHistogramHandlerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
// see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeHistogram );
}

AliHLTComponentDataType AliHLTTPCHistogramHandlerComponent::GetOutputDataType() { 
// see header file for class documentation

  return kAliHLTDataTypeHistogram;
}

int AliHLTTPCHistogramHandlerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList) { 
// see header file for class documentation

  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeHistogram);
  return tgtList.size();
}

void AliHLTTPCHistogramHandlerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
// see header file for class documentation

  constBase=0;
  inputMultiplier=2.0;
}

AliHLTComponent* AliHLTTPCHistogramHandlerComponent::Spawn() { 
// see header file for class documentation

  return new AliHLTTPCHistogramHandlerComponent();
}
	
int AliHLTTPCHistogramHandlerComponent::DoInit( int argc, const char** argv ) { 
// see header file for class documentation

  Int_t i = 0;
  Char_t* cpErr;
  
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
  } else {
    iResult=Reconfigure(NULL, NULL);
  }

 
  while ( i < argc ) {      
    if (!strcmp( argv[i], "-sum-noise-histograms")) {
        fNoiseHistograms = strtoul( argv[i+1], &cpErr ,0);
            
    if ( *cpErr ) {
        HLTError("Cannot convert sum-noise-histograms specifier '%s'.", argv[i+1]);
        return EINVAL;
    }
      i+=2;
      continue;
    }
    
    if (!strcmp( argv[i], "-sum-krypton-histograms")) {
        fKryptonHistograms = strtoul( argv[i+1], &cpErr ,0);
            
    if ( *cpErr ) {
        HLTError("Cannot convert sum-krypton-histograms specifier '%s'.", argv[i+1]);
        return EINVAL;
    }
      i+=2;
      continue;
    }    
                  
    Logging(kHLTLogError, "HLT::TPCHistogramHandler::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  } // end while
  
  return 0;
} // end DoInit()

int AliHLTTPCHistogramHandlerComponent::DoDeinit() { 
// see header file for class documentation 
   
   return 0;
}

int AliHLTTPCHistogramHandlerComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData){
// see header file for class documentation

  HLTInfo("--- Entering DoEvent() in TPCHistogramHandler ---");
  
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;
  
  fHistTH2Tmp = new TH2F("fHistTH2Tmp","fHistTH2Tmp",250,-250,250,250,-250,250);    
  fHistTPCSideA = new TH2F("fHistTPCSideA","TPC side A (max signal)",250,-250,250,250,-250,250);
  fHistTPCSideA->SetXTitle("global X (cm)"); fHistTPCSideA->SetYTitle("global Y (cm)");
  fHistTPCSideC = new TH2F("fHistTPCSideC","TPC side C (max signal)",250,-250,250,250,-250,250);
  fHistTPCSideC->SetXTitle("global X (cm)"); fHistTPCSideC->SetYTitle("global Y (cm)");
  
  const TObject *iter = NULL;  
        
  for(iter = GetFirstInputObject(kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC); iter != NULL; iter = GetNextInputObject()){
   
     HLTInfo("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", 
              evtData.fEventID, evtData.fEventID,
              DataType2Text(GetDataType(iter)).c_str(), 
              DataType2Text(kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC).c_str());
   
     if (GetDataType(iter) == (kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC) && GetEventCount()<2){
         HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeHistogram)!", 
	 DataType2Text(kAliHLTDataTypeHistogram).c_str(),
         DataType2Text(kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC).c_str());
     }      
     
     if (GetDataType(iter) != (kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC)) continue;
         
     // Summing the output histograms of the TPCNoiseMapComponent (from partition to TPC sides)
     if(fNoiseHistograms){  
       
       fHistTH2Tmp = (TH2F*)iter;
       UInt_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(GetSpecification(iter)); 
       UInt_t maxSlice     = AliHLTTPCDefinitions::GetMaxSliceNr(GetSpecification(iter)); 
       UInt_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(GetSpecification(iter)); 
       UInt_t maxPartition = AliHLTTPCDefinitions::GetMaxPatchNr(GetSpecification(iter)); 
       
       if((minSlice!=maxSlice) || (minPartition!=maxPartition)){
           HLTWarning("TPCHistogramHandler::The Noise Map component is not running on partition level!");
       }

       if(minSlice<18) fHistTPCSideA->Add(fHistTPCSideA,fHistTH2Tmp,1,1);
       else	       fHistTPCSideC->Add(fHistTPCSideC,fHistTH2Tmp,1,1);
       // minSlice=maxSlice, when the Noise Map component runs on partition level (as it should)
       
       fSpecificationTPCA = AliHLTTPCDefinitions::EncodeDataSpecification( 17,  0, 5, 0 ); 
       fSpecificationTPCC = AliHLTTPCDefinitions::EncodeDataSpecification( 35, 18, 5, 0 ); 

     } // endif fNoiseHistograms==kTRUE   
     
     
//  Summing the output of TPCKryptonClusterFinderComponent
//      if(fKryptonHistograms){
//         
//      } //endif fKryptonHistograms==kTRUE	   	         
  } // end for loop over histogram blocks
  
  MakeHistosPublic();

  return 0;
} // end DoEvent()

void AliHLTTPCHistogramHandlerComponent::MakeHistosPublic() {
// see header file for class documentation
  
 
  if(fNoiseHistograms){ 
    PushBack( (TObject*) fHistTPCSideA, kAliHLTDataTypeHistogram, fSpecificationTPCA);
    PushBack( (TObject*) fHistTPCSideC, kAliHLTDataTypeHistogram, fSpecificationTPCC);
    delete fHistTH2Tmp;
    delete fHistTPCSideA;
    delete fHistTPCSideC;
  }
 
//  TObjArray histos;
   
//   if(fPlotSideA) histos.Add(fHistSideA);
//   if(fPlotSideC) histos.Add(fHistSideC);
//   if(fApplyNoiseMap) histos.Add(fHistCDBMap);
//   
//   TIter iterator(&histos);
//   while(TObject *pObj=iterator.Next()){
//   
//         PushBack(pObj, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC, fSpecification);
//   }
//   
//   
//   //PushBack( (TObject*) &histos, kAliHLTDataTypeHistogram, fSpecification);    
}

int AliHLTTPCHistogramHandlerComponent::Configure(const char* arguments) { 
// see header file for class documentation
  
  int iResult=0;
  if (!arguments) return iResult;
  HLTInfo("parsing configuration string \'%s\'", arguments);

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
     
      if (argument.CompareTo("-sum-noise-histograms")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-sum-noise-histograms\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } else if (argument.CompareTo("-sum-krypton-histograms")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-sum-krypton-histograms\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    } // end for
  
    delete pTokens;
  
  } // end if pTokens
  
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTPCHistogramHandlerComponent::Reconfigure(const char* cdbEntry, const char* chainId) { 
// see header file for class documentation

  int iResult=0;
  const char* path="HLT/ConfigTPC/TPCHistogramHandlerComponent";
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
      HLTError("cannot fetch object \"%s\" from CDB", path);
    }
  }
  return iResult;
}
