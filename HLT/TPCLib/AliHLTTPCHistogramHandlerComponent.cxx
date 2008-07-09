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
#include "AliHLTTPCTransform.h"

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

ClassImp(AliHLTTPCHistogramHandlerComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCHistogramHandlerComponent::AliHLTTPCHistogramHandlerComponent()
    :    
    fSpecification(0),
    fNoiseHistograms(0),
    fKryptonHistograms(0),
    fSlice(-99),
    
    fHistTH1Tmp(NULL),
    fTotalClusterChargeIROCAll(NULL),
    fTotalClusterChargeOROCAll(NULL),
    fQMaxPartitionAll(NULL),
    fPlotQmaxROCAll(NULL),
    fNumberOfClusters(NULL),
    
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
  for (int j=0; j<argc && iResult>=0; j++) {
    
    argument=argv[j];
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

int AliHLTTPCHistogramHandlerComponent::DoEvent(const AliHLTComponentEventData&/* evtData*/, AliHLTComponentTriggerData& /*trigData*/){
// see header file for class documentation

  HLTInfo("--- Entering DoEvent() in TPCHistogramHandler ---");
  
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;
 
  fTotalClusterChargeIROCAll = new TH1F("fTotalClusterChargeIROCAll","Total Charge of clusters in all IROC",4000,0,4000);  
  fTotalClusterChargeOROCAll = new TH1F("fTotalClusterChargeOROCAll","Total Charge of clusters in all OROC",4000,0,4000);
  fQMaxPartitionAll          = new TH1F("fQMaxPartitionAll",         "QMax for All Partitions",             216,0,216);
  fPlotQmaxROCAll            = new TH1F("fQMaxROCAll",               "QMax for All ROC",                    72,0,72);
  fNumberOfClusters          = new TH1F("fNumberOfClusters",         "Total Number of Clusters",            1,0,1);
    
  fHistTH2Tmp = new TH2F("fHistTH2Tmp","fHistTH2Tmp",250,-250,250,250,-250,250);    
  fHistTPCSideA = new TH2F("fHistTPCSideA","TPC side A (max signal)",250,-250,250,250,-250,250);
  fHistTPCSideA->SetXTitle("global X (cm)"); fHistTPCSideA->SetYTitle("global Y (cm)");
  fHistTPCSideC = new TH2F("fHistTPCSideC","TPC side C (max signal)",250,-250,250,250,-250,250);
  fHistTPCSideC->SetXTitle("global X (cm)"); fHistTPCSideC->SetYTitle("global Y (cm)");
  
  const TObject *iter = NULL;  
        
  for(iter = GetFirstInputObject(kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC); iter != NULL; iter = GetNextInputObject()){
   
//      HLTInfo("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", 
//               evtData.fEventID, evtData.fEventID,
//               DataType2Text(GetDataType(iter)).c_str(), 
//               DataType2Text(kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC).c_str());
   
//      if (GetDataType(iter) == (kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC) && GetEventCount()<2){
//          HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeHistogram)!", 
// 	 DataType2Text(kAliHLTDataTypeHistogram).c_str(),
//          DataType2Text(kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC).c_str());
//      }      
     
     if (GetDataType(iter) != (kAliHLTDataTypeHistogram | kAliHLTDataOriginTPC)) continue;
         
     // Summing the output histograms of the AliHLTTPCNoiseMapComponent (from partition to TPC sides)
     if(fNoiseHistograms){  
       
        fHistTH2Tmp = (TH2F*)iter;
        UInt_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(GetSpecification(iter)); 
        UInt_t maxSlice     = AliHLTTPCDefinitions::GetMaxSliceNr(GetSpecification(iter)); 
        UInt_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(GetSpecification(iter)); 
        UInt_t maxPartition = AliHLTTPCDefinitions::GetMaxPatchNr(GetSpecification(iter)); 
       
        if((minSlice!=maxSlice) || (minPartition!=maxPartition)){
            HLTWarning("TPCHistogramHandler::The Noise Map component is not running on partition level!");
        }

        // minSlice=maxSlice, when the Noise Map component runs on partition level (as it should)
        if(minSlice<18) fHistTPCSideA->Add(fHistTPCSideA,fHistTH2Tmp,1,1);
        else		fHistTPCSideC->Add(fHistTPCSideC,fHistTH2Tmp,1,1);
     } // endif fNoiseHistograms==kTRUE   
     
     
     // Summing the output of AliHLTTPCClusterHistoComponent
     if(fKryptonHistograms){
       Int_t thisrow=-1,thissector=-1,row=-1;
       
       AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(GetSpecification(iter));
       AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr(GetSpecification(iter));
       row = AliHLTTPCTransform::GetFirstRow(patch); 
       AliHLTTPCTransform::Slice2Sector(slice,row,thissector,thisrow);
       
       fHistTH1Tmp = (TH1F*)iter;   	
       //cout << fHistTH1Tmp->GetName() << "\t" << fHistTH1Tmp->GetEntries() << endl;
       
	TString name = fHistTH1Tmp->GetName();
		
	if(name=="fTotalClusterChargeIROCAll"){
	   fTotalClusterChargeIROCAll->Add(fTotalClusterChargeIROCAll,fHistTH1Tmp,1,1);
	} 
	else if(name=="fTotalClusterChargeOROCAll"){
	   fTotalClusterChargeOROCAll->Add(fTotalClusterChargeOROCAll,fHistTH1Tmp,1,1);
	} 
	else if(name=="fQMaxPartitionAll"){
	  for(Int_t t=0;t<216;t++){
	    if(fHistTH1Tmp->GetBinContent(t)>fQMaxPartitionAll->GetBinContent(t)){
	      fQMaxPartitionAll->SetBinContent(t,fHistTH1Tmp->GetBinContent(t));
	    }
	  } 
	} 
	else if(name=="fQMaxROCAll"){
	  for(Int_t t=0;t<72;t++){
	    if(fHistTH1Tmp->GetBinContent(t)>fPlotQmaxROCAll->GetBinContent(t)){
	      fPlotQmaxROCAll->SetBinContent(t,fHistTH1Tmp->GetBinContent(t));
	    }
	  }
	} 
	else if(name=="fNumberOfClusters"){ 
	  fNumberOfClusters->Add(fNumberOfClusters,fHistTH1Tmp,1,1);
	} 
	else{
	  HLTWarning("No histogram names match. %s",name.Data());
	  continue;
	}     
     } //endif fKryptonHistograms==kTRUE	   	         
  } // end for loop over histogram blocks
  
  MakeHistosPublic();
  return 0;
} // end DoEvent()

void AliHLTTPCHistogramHandlerComponent::MakeHistosPublic() {
// see header file for class documentation
 
  if(fNoiseHistograms){ 
    PushBack((TObject*)fHistTPCSideA,kAliHLTDataTypeHistogram,AliHLTTPCDefinitions::EncodeDataSpecification( 0,17,0,5));
    PushBack((TObject*)fHistTPCSideC,kAliHLTDataTypeHistogram,AliHLTTPCDefinitions::EncodeDataSpecification(18,35,0,5));

    if(fHistTH2Tmp)  { delete fHistTH2Tmp;   fHistTH2Tmp=NULL;   }
    if(fHistTPCSideA){ delete fHistTPCSideA; fHistTPCSideA=NULL; }
    if(fHistTPCSideC){ delete fHistTPCSideC; fHistTPCSideC=NULL; }
  }  
  
  if(fKryptonHistograms){
     PushBack((TObject*)fTotalClusterChargeIROCAll,kAliHLTDataTypeHistogram,AliHLTTPCDefinitions::EncodeDataSpecification(0,17,0,1));
     PushBack((TObject*)fTotalClusterChargeOROCAll,kAliHLTDataTypeHistogram,AliHLTTPCDefinitions::EncodeDataSpecification(0,17,2,5));
     PushBack((TObject*)fQMaxPartitionAll,	   kAliHLTDataTypeHistogram,AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5));
     PushBack((TObject*)fPlotQmaxROCAll,	   kAliHLTDataTypeHistogram,AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5));
     PushBack((TObject*)fNumberOfClusters,	   kAliHLTDataTypeHistogram,AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5));
          
     if(fTotalClusterChargeIROCAll){ delete fTotalClusterChargeIROCAll; fTotalClusterChargeIROCAll=NULL; }
     if(fTotalClusterChargeOROCAll){ delete fTotalClusterChargeOROCAll; fTotalClusterChargeOROCAll=NULL; }
     if(fQMaxPartitionAll)         { delete fQMaxPartitionAll;          fQMaxPartitionAll=NULL;          }
     if(fPlotQmaxROCAll)           { delete fPlotQmaxROCAll;            fPlotQmaxROCAll=NULL;            }
     if(fNumberOfClusters)         { delete fNumberOfClusters;          fNumberOfClusters=NULL;          }
     if(fHistTH1Tmp)               { delete fHistTH1Tmp;                fHistTH1Tmp=NULL;                }
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
