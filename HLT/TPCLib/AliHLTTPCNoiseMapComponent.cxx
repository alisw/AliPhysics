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

/** @file   AliHLTTPCNoiseMapComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  The TPC Noise Map component
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCNoiseMapComponent.h"
#include "AliHLTTPCDigitReaderDecoder.h"
#include "AliHLTTPCDigitReaderPacked.h"
#include "AliHLTTPCTransform.h"
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
#include "TH2.h"


AliHLTTPCNoiseMapComponent gAliHLTTPCNoiseMapComponent;

ClassImp(AliHLTTPCNoiseMapComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCNoiseMapComponent::AliHLTTPCNoiseMapComponent()
    :    
    fNoiseMap(0),
    fIsPacked(1),
    fIsUnpacked(0),
    fCurrentPartition(0),
    fCurrentRow(0)   
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCNoiseMapComponent::~AliHLTTPCNoiseMapComponent() { 
// see header file for class documentation

}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCNoiseMapComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCNoiseMap";
}

void AliHLTTPCNoiseMapComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
// see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCNoiseMapComponent::GetOutputDataType() { 
// see header file for class documentation

  return AliHLTTPCDefinitions::fgkNoiseHistoDataType;
}

int AliHLTTPCNoiseMapComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList) { 
// see header file for class documentation

  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkNoiseHistoDataType);
  return tgtList.size();
}

void AliHLTTPCNoiseMapComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
// see header file for class documentation

  constBase=0;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTTPCNoiseMapComponent::Spawn() { 
// see header file for class documentation

  return new AliHLTTPCNoiseMapComponent();
}
	
int AliHLTTPCNoiseMapComponent::DoInit( int argc, const char** argv ) { 
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
    
   
    // -- fill noise histograms
    if (!strcmp( argv[i], "-noisemap")) {
        fNoiseMap = strtoul( argv[i+1], &cpErr ,0);
      
    //if(fNoiseMap) {
      
       fHistSideC = new TH2F("fHistSideC","",500,-250,250,500,-250,250);    
       HLTInfo("---- HAVE CREATED HISTOGRAM(S) ----");
       //fSideA, fSideC, fSliceA[18], fSliceC[18]
    
    //}
    if ( *cpErr ) {
        HLTError("Cannot convert noisemap specifier '%s'.", argv[i+1]);
        return EINVAL;
    }
      i+=2;
      continue;
    }
        
    Logging(kHLTLogError, "HLT::TPCNoiseMap::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  }

  //HLTDebug("using AliHLTTPCDigitReaderDecoder");
  //pDigitReader = new AliHLTTPCDigitReaderDecoder(); // double-loop
  //pDigitReader = new AliHLTTPCDigitReaderPacked();
  
//   if(fIsPacked){
//     //cout<<"Digit reader decoder is chosen"<<endl;
//     fDigitReader = new AliHLTTPCDigitReaderDecoder();
//   }
//   else if(fIsUnpacked){
//     //cout<<"Digit reader unpacked is chosen"<<endl;
//     fDigitReader = new AliHLTTPCDigitReaderDecoderUnpacked();
//   }
//   else{
//     HLTFatal("Neither of the two options of digit readers is set no data will be read.");
//   }


  return 0;
} // end DoInit()

int AliHLTTPCNoiseMapComponent::DoDeinit() { 
// see header file for class documentation 
      
    return 0;
}

int AliHLTTPCNoiseMapComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData&, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, 
					      vector<AliHLTComponentBlockData>& outputBlocks ) { 
// see header file for class documentation
 if(fNoiseMap) { 
 
  HLTInfo("--- Entering DoEvent() in TPCNoiseMap ---");

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;
  
  Float_t xyz[3]; 
  Int_t thissector, thisrow;
  

  for ( ndx=0; ndx<evtData.fBlockCnt; ndx++ ) {
      iter = blocks+ndx;
      
      HLTInfo("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", evtData.fEventID, evtData.fEventID, 
	       DataType2Text( iter->fDataType).c_str(), DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());

      if (iter->fDataType == AliHLTTPCDefinitions::fgkDDLPackedRawDataType && GetEventCount()<2) {
	  HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeDDLRaw)!", DataType2Text(AliHLTTPCDefinitions::fgkDDLPackedRawDataType).c_str(),
	 	      DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
      } //endif      
      
      if (iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC)) continue;
  
     
      UInt_t slice     = AliHLTTPCDefinitions::GetMinSliceNr( *iter ); 
      UInt_t partition = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
          
      //if ( partition < fMinPartition ) fMinPartition = partition;
      //if ( partition > fMaxPartition ) fMaxPartition = partition; // add a warning
      
      //fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( slice, slice, partition, partition );
      fSpecification = iter->fSpecification;
      
      AliHLTTPCDigitReader *pDigitReader = new AliHLTTPCDigitReaderDecoder;
      pDigitReader->InitBlock(iter->fPtr,iter->fSize,partition,slice);
      if(!pDigitReader) break; //AliHLTComponent.cxx, altrochannelselector rcu folder
            
      rowOffset = 0;      
      
      switch(partition){
      case 0:
	rowOffset=AliHLTTPCTransform::GetFirstRow(0);
	break;
      case 1:
	rowOffset=AliHLTTPCTransform::GetFirstRow(1);
	break;
      case 2:
	rowOffset=AliHLTTPCTransform::GetFirstRow(2);
	break;
      case 3:
	rowOffset=AliHLTTPCTransform::GetFirstRow(3);
	break;
      case 4:
	rowOffset=AliHLTTPCTransform::GetFirstRow(4);
	break;
      case 5:
	rowOffset=AliHLTTPCTransform::GetFirstRow(5);
	break;
      }
          	
      
        //while( pDigitReader->Next() ){ 
        while( pDigitReader->NextChannel()) { 
	 
	 fCurrentRow  = pDigitReader->GetRow();	 
	 fCurrentRow += rowOffset;
	
	 if(pDigitReader->GetPad()==1000 || fCurrentRow==1000) continue;

	 AliHLTTPCTransform::Slice2Sector(slice,fCurrentRow,thissector,thisrow);
	 AliHLTTPCTransform::Raw2Global(xyz,thissector,thisrow,pDigitReader->GetPad(),0);
	 //AliHLTTPCTransform::Raw2Local(xyz,thissector,thisrow,pDigitReader->GetPad(),0);
	 	 
	 //use AliTPCPad to fill the data there and use the functions to ask for max charge etc.

	 Int_t  maxSignal = 0;
	 while( pDigitReader->NextBunch()) {
     
	     //in case we want to fill the histograms with the signal value, an additional loop is necessary      	                              
             for(Int_t i=0;i<pDigitReader->GetBunchSize();i++) {
	     	
	     	 const UInt_t *bunchData = pDigitReader->GetSignals();
	     	 if(bunchData[i]>maxSignal) maxSignal = bunchData[i]; 
             	 //cout<<"Time: "<<pDigitReader->GetTime()+i<<"    Signal: "<<bunchData[i]<<endl;
             }  					   
	 } // end of inner while loop

cout<<slice<<" "<<partition<<" "<<fCurrentRow<<" "<<pDigitReader->GetPad()<<" "<<xyz[0]<<" "<<xyz[1]<<" "<<maxSignal<<endl;
	  
	  fHistSideC->Fill(xyz[0],xyz[1],maxSignal);
      
        } // end of outer while pad loop
  } // end of data block loop 
 
  SaveAndResetHistograms();
 
 } // end if noisemap==kTRUE  
 return 0;
} // end DoEvent()

void AliHLTTPCNoiseMapComponent::SaveAndResetHistograms() {
// see header file for class documentation
  
  TFile *outputfile = new TFile("test.root","RECREATE");
  fHistSideC->Write();
  outputfile->Save();
  outputfile->Close();
   
  PushBack( (TObject*) fHistSideC, AliHLTTPCDefinitions::fgkNoiseHistoDataType, fSpecification);
  //fill it with the right specification for every histogram
  //make a TObjArray and add all histos
  //check which pointers are empty and publish only the ones that hold something
    
  delete fHistSideC;
}

int AliHLTTPCNoiseMapComponent::Configure(const char* arguments) { 
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
     
      if (argument.CompareTo("-noisemap")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-noisemap\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    } // endfor
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;

}

int AliHLTTPCNoiseMapComponent::Reconfigure(const char* cdbEntry, const char* chainId) { 
// see header file for class documentation

  int iResult=0;
  const char* path="HLT/ConfigTPC/TPCNoiseMapComponent";
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
