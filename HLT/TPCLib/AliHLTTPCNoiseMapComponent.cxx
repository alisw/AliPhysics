/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Kenneth Aamodt <Kenneth.Aamodt@student.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCNoiseMapComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  The TPC Noise Map component
*/

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

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



AliHLTTPCNoiseMapComponent gAliHLTTPCNoiseMapComponent;


ClassImp(AliHLTTPCNoiseMapComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCNoiseMapComponent::AliHLTTPCNoiseMapComponent()
    :
    fFirstTimeBin(0),
    fLastTimeBin(AliHLTTPCTransform::GetNTimeBins()),
    fNSigmaThreshold(0),
    fSignalThreshold(0),
    fMinimumNumberOfSignals(0),
    fOldRCUFormat(0),
    fSortPads(0),
    fDigitReader(NULL),
    fNoiseMap(0),
    fSpecification(0),
    fMaxPatch(0),
    fMinPatch(0), 
    fCurrentPatch(0),
    fCurrentRow(0),
    fIsPacked(1),
    fIsUnpacked(0)


{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCNoiseMapComponent::~AliHLTTPCNoiseMapComponent() { // see header file for class documentation

}


// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCNoiseMapComponent::GetComponentID() { // see header file for class documentation

  return "TPCNoiseMap";
}

void AliHLTTPCNoiseMapComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { // see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCNoiseMapComponent::GetOutputDataType() { // see header file for class documentation

  //return AliHLTTPCDefinitions::fgkUnpackedRawDataType;
  return AliHLTTPCDefinitions::fgkNoiseHistoDataType;
}

int AliHLTTPCNoiseMapComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList) { // see header file for class documentation

  tgtList.clear();
  //tgtList.push_back(AliHLTTPCDefinitions::fgkUnpackedRawDataType);
  tgtList.push_back(AliHLTTPCDefinitions::fgkNoiseHistoDataType);
  return tgtList.size();
}

void AliHLTTPCNoiseMapComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { // see header file for class documentation

  constBase=0;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTTPCNoiseMapComponent::Spawn() { // see header file for class documentation

  return new AliHLTTPCNoiseMapComponent();
}
	
int AliHLTTPCNoiseMapComponent::DoInit( int argc, const char** argv ) { // see header file for class documentation

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
      
    if(fNoiseMap) {
      
       hPatch = new TH2F("hPatch","",500,-250,250,500,-250,250);    
       HLTInfo("---- HAVE CREATED HISTOGRAM(S) ----");
    
    }
    if ( *cpErr ) {
        HLTError("Cannot convert noisemap specifier '%s'.", argv[i+1]);
        return EINVAL;
    }
      i+=2;
      continue;
    }
    
    
//     // -- zero suppression threshold
//     if ( !strcmp( argv[i], "adc-threshold" ) ) {
//       fSignalThreshold = strtoul( argv[i+1], &cpErr ,0);
//       if ( *cpErr ) {
// 	HLTError("Cannot convert adc-threshold specifier '%s'.", argv[i+1]);
// 	return EINVAL;
//       }
//       i+=2;
//       continue;
//     }

      
    Logging(kHLTLogError, "HLT::TPCNoiseMap::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  }

  HLTDebug("using AliHLTTPCDigitReaderDecoder");
  fDigitReader = new AliHLTTPCDigitReaderDecoder(); // double-loop
  //fDigitReader = new AliHLTTPCDigitReaderPacked();
  
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

int AliHLTTPCNoiseMapComponent::DoDeinit() { // see header file for class documentation 
      
    return 0;
}

int AliHLTTPCNoiseMapComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData&, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, 
					      vector<AliHLTComponentBlockData>& outputBlocks ) { // see header file for class documentation
 if(fNoiseMap) { 
 
  HLTInfo("--- Entering DoEvent() in TPCNoiseMap ---");

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;
  
  //cout << "Number of blocks: " << evtData.fBlockCnt << endl;

  Float_t xyz[3]; 
  Int_t thissector, thisrow;
  Int_t rowOffset = 0;
  

  for ( ndx=0; ndx<evtData.fBlockCnt; ndx++ ) {

      iter = blocks+ndx;
      
      HLTInfo("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", evtData.fEventID, evtData.fEventID, 
	       DataType2Text( iter->fDataType).c_str(), DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());

      if (iter->fDataType == AliHLTTPCDefinitions::fgkDDLPackedRawDataType && GetEventCount()<2) {
	  HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeDDLRaw)!", DataType2Text(AliHLTTPCDefinitions::fgkDDLPackedRawDataType).c_str(),
	 	      DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
      } //endif      
      
//       if ( iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC) &&  iter->fDataType != AliHLTTPCDefinitions::fgkDDLPackedRawDataType ) {
// 	   continue;
//       } // endif
  
     
      UInt_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      UInt_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
      

      //if ( patch >= 2 ) // Outer sector, patches 2, 3, 4, 5
      //rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );      
      //fCurrentRow += rowOffset;
      

      if ( patch < fMinPatch ) fMinPatch =  patch;
      if ( patch > fMaxPatch ) fMaxPatch =  patch;
      
      
      fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( slice, slice, fMaxPatch, fMinPatch );
      fDigitReader->InitBlock(iter->fPtr,iter->fSize,patch,slice);
           	
      if(!fSortPads) { //the code within this if() does NOT sort the pads, and will only work with the UNSORTED (new) clusterfinding
      
        //while( fDigitReader->Next() ){ 
        while( fDigitReader->NextChannel()) { 
          while( fDigitReader->NextBunch()) {
          
             
             //AliHLTTPCTransform::Slice2Sector(slice,fCurrentRow,thissector,thisrow);            
	     AliHLTTPCTransform::Slice2Sector(slice,fDigitReader->GetRow(),thissector,thisrow);      
	     AliHLTTPCTransform::Raw2Global(xyz,slice,thisrow,fDigitReader->GetPad(),fDigitReader->GetTime());
             
	     // this is where the messing up takes place, I am not sure about the arguments of the above 2 functions

             
	     //hPatch->Fill(xyz[0],xyz[1]); // as a consequence, this looks strange
	     hPatch->Fill(fDigitReader->GetPad(), fDigitReader->GetRow());
	     
	     
	     //in case we want to fill the histograms with the signal value, an additional loop is necessary
//              for(Int_t i=0;i<fDigitReader->GetBunchSize();i++){ 
//           	
//                AliHLTTPCTransform::Raw2Global(xyz,slice,thisrow,fDigitReader->GetPad(),fDigitReader->GetTime()+i);              
//                //  const UInt_t *bunchData = fDigitReader->GetSignals();
//                //  for(Int_t i=0;i<fDigitReader->GetBunchSize();i++) {
//                //      cout<<"Time: "<<fDigitReader->GetTime()+i<<"    Signal: "<<bunchData[i]<<endl;
//                //  } 	                             
// 	     } // end for
           } // end of inner while loop
        } // end of outer while loop
      } // endif(!fSortPads)                   
  } // end for loop 
 
  SaveAndResetHistograms();
 
 }    
  
  return 0;
} // end DoEvent()

void AliHLTTPCNoiseMapComponent::SaveAndResetHistograms() {
  
  
  TFile *outputfile = new TFile("test.root","RECREATE");
  hPatch->Write();
  
  PushBack( (TObject*) hPatch, AliHLTTPCDefinitions::fgkNoiseHistoDataType, fSpecification);
  
  outputfile->Save();
  outputfile->Close();
  
  delete hPatch;
}

int AliHLTTPCNoiseMapComponent::Configure(const char* arguments) { // see header file for class documentation
  
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

int AliHLTTPCNoiseMapComponent::Reconfigure(const char* cdbEntry, const char* chainId) { // see header file for class documentation

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
