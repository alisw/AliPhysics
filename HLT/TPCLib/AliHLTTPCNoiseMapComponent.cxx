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
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDefinitions.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliHLTTPCNoiseMap.h"

#include "AliTPCCalPad.h"
#include "AliTPCROC.h"
#include "AliTPCCalROC.h"

#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <sys/time.h>
#include "TH2.h"
#include "TH3.h"

ClassImp(AliHLTTPCNoiseMapComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCNoiseMapComponent::AliHLTTPCNoiseMapComponent()
    :    
    fSpecification(0),
    fReadNoiseMap(0),
    fResetHistograms(0),
    fInitHist(kTRUE),    
    fCurrentRow(-99),
    fHistSignal(NULL),
    fHistSideAMaxSignal(NULL),
    fHistSideATotSignal(NULL),
    fHistSideAPadRMS(NULL),
    fHistSideCMaxSignal(NULL),
    fHistSideCTotSignal(NULL),
    fHistSideCPadRMS(NULL),
    fHistCDBMap(NULL)
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

  return kAliHLTDataTypeHistogram;
}

int AliHLTTPCNoiseMapComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList) { 
// see header file for class documentation

  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeHistogram);
  return tgtList.size();
}

void AliHLTTPCNoiseMapComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
// see header file for class documentation

  constBase=800000;
  inputMultiplier=0.0;
}

AliHLTComponent* AliHLTTPCNoiseMapComponent::Spawn() { 
// see header file for class documentation

  return new AliHLTTPCNoiseMapComponent();
}
	
int AliHLTTPCNoiseMapComponent::DoInit( int argc, const char** argv ) { 
// see header file for class documentation
 
  //Int_t i = 0;
  //Char_t* cpErr;
  
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
    //    iResult=Reconfigure(NULL, NULL);
  }

 
//   while ( i < argc ) {      
//     if (!strcmp( argv[i], "-read-noisemap")) {
//         fApplyNoiseMap = strtoul( argv[i+1], &cpErr ,0);
//             
//     if ( *cpErr ) {
//         HLTError("Cannot convert apply-noisemap specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
//     
//     if (!strcmp( argv[i], "-reset-histograms")) {
//         fResetHistograms = strtoul( argv[i+1], &cpErr ,0);
//     
//     if ( *cpErr ) {
//         HLTError("Cannot convert reset-histograms specifier '%s'.", argv[i+1]);
//         return EINVAL;
//     }
//       i+=2;
//       continue;
//     }
//                    
//     Logging(kHLTLogError, "HLT::TPCNoiseMap::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
//     return EINVAL;
// 
//   } // end while
  
  if(fReadNoiseMap){
     AliHLTTPCNoiseMap *nm = AliHLTTPCNoiseMap::Instance();
     if(!nm) { 
         HLTWarning("AliHLTTPCNoiseMap instance not existent."); 
     }
     else {
         AliTPCCalPad *noisePad = nm->ReadNoiseMap(GetRunNo());
         if(noisePad) {
	    fHistCDBMap = noisePad->MakeHisto2D(1);
	    
	 }
     }
  }

//   if(fApplyNoiseMap){
//     //TFile *f = TFile::Open("/scratch/noiseComp/Run3398_4000_v0_s72.root");
//     TFile *f = TFile::Open("/home/kanaki/noiseComp/Run3398_4000_v0_s72.root");
//     AliCDBEntry *pEntry = (AliCDBEntry*)f->Get("AliCDBEntry"); 
//     noisePad = (AliTPCCalPad*)pEntry->GetObject();
//     //fHistCDBMap = noisePad->MakeHisto2D(1); //side C
//   }
   
//   HLTDebug("using AliHLTTPCDigitReaderDecoder");
//   pDigitReader = new AliHLTTPCDigitReaderDecoder(); // double-loop
//   pDigitReader = new AliHLTTPCDigitReaderPacked();
  
  return 0;

} // end DoInit()

int AliHLTTPCNoiseMapComponent::DoDeinit() { 
  // see header file for class documentation  

  if(fHistSideAMaxSignal) delete fHistSideAMaxSignal; fHistSideAMaxSignal = NULL;
  if(fHistSideATotSignal) delete fHistSideATotSignal; fHistSideATotSignal = NULL;
  if(fHistSideAPadRMS)    delete fHistSideAPadRMS;    fHistSideAPadRMS    = NULL;
  
  if(fHistSideCMaxSignal) delete fHistSideCMaxSignal; fHistSideCMaxSignal = NULL;
  if(fHistSideCTotSignal) delete fHistSideCTotSignal; fHistSideCTotSignal = NULL;
  if(fHistSideCPadRMS)	  delete fHistSideCPadRMS;    fHistSideCPadRMS    = NULL;
         
  return 0;
}

void AliHLTTPCNoiseMapComponent::InitializeHistograms(UInt_t minSlice, UInt_t maxSlice, UInt_t minPartition, UInt_t maxPartition){
  // see header file for class documentation
   
  Char_t name1[50], name2[50], name3[50];
 
  if(minSlice<18){
     sprintf(name1, "fHistSideAMaxSignal_Slice_%.2d%.2d_Partition_%.2d%.2d", minSlice, maxSlice, minPartition, maxPartition);
     sprintf(name2, "fHistSideATotSignal_Slice_%.2d%.2d_Partition_%.2d%.2d", minSlice, maxSlice, minPartition, maxPartition);
     sprintf(name3, "fHistSideAPadRMS_Slice_%.2d%.2d_Partition_%.2d%.2d",    minSlice, maxSlice, minPartition, maxPartition);
     fHistSideAMaxSignal = new TH2F(name1,name1,250,-250,250,250,-250,250);
     fHistSideATotSignal = new TH2F(name2,name2,250,-250,250,250,-250,250);
     fHistSideAPadRMS	 = new TH2F(name3,name3,250,-250,250,250,-250,250);

  } else {
     sprintf(name1, "fHistSideCMaxSignal_Slice_%.2d%.2d_Partition_%.2d%.2d", minSlice, maxSlice, minPartition, maxPartition);
     sprintf(name2, "fHistSideCTotSignal_Slice_%.2d%.2d_Partition_%.2d%.2d", minSlice, maxSlice, minPartition, maxPartition);
     sprintf(name3, "fHistSideCPadRMS_Slice_%.2d%.2d_Partition_%.2d%.2d",    minSlice, maxSlice, minPartition, maxPartition);  
     fHistSideCMaxSignal = new TH2F(name1,name1,250,-250,250,250,-250,250);
     fHistSideCTotSignal = new TH2F(name2,name2,250,-250,250,250,-250,250);
     fHistSideCPadRMS	 = new TH2F(name3,name3,250,-250,250,250,-250,250);
  }
    
  fInitHist=kFALSE;
  
}

int AliHLTTPCNoiseMapComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/){
  // see header file for class documentation
 
 
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;
  //HLTInfo("--- Entering DoEvent() in TPCNoiseMap ---");
   
  const AliHLTComponentBlockData *iter = NULL;

  Float_t xyz[3]; 
  Int_t thissector, thisrow;
    
  for(iter = GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC); iter != NULL; iter = GetNextInputBlock()){
      
//     HLTInfo("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", 
// 	    evtData.fEventID, evtData.fEventID,
// 	    DataType2Text(iter->fDataType).c_str(), 
// 	    DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
// 
//     if (iter->fDataType == AliHLTTPCDefinitions::fgkDDLPackedRawDataType && GetEventCount()<2){
//       HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeDDLRaw)!", 
// 		 DataType2Text(AliHLTTPCDefinitions::fgkDDLPackedRawDataType).c_str(),
// 		 DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
//     }      
     
    if (iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC)) continue;
      
    UInt_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter); 
    UInt_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);
    UInt_t maxSlice     = AliHLTTPCDefinitions::GetMaxSliceNr(*iter); 
    UInt_t maxPartition = AliHLTTPCDefinitions::GetMaxPatchNr(*iter);
    
    if(fInitHist==kTRUE) InitializeHistograms(minSlice, maxSlice, minPartition, maxPartition);    
    fSpecification = iter->fSpecification;
     
    AliHLTTPCDigitReader *pDigitReader = new AliHLTTPCDigitReaderDecoder;

    pDigitReader->InitBlock(iter->fPtr,iter->fSize,minPartition,minSlice);
    if(!pDigitReader) break;
                  
    //  while(pDigitReader->Next()){ 

    Float_t maxSignal     = 0.;
    Float_t totalSignal   = 0.;
    Float_t squaredSignal = 0.;
    Float_t rms = 0.; 
    
    while(pDigitReader->NextChannel()) { // pad loop 
      
      fCurrentRow  = pDigitReader->GetRow();  
      fCurrentRow += pDigitReader->GetRowOffset();
      

      AliHLTTPCTransform::Slice2Sector(minSlice,fCurrentRow,thissector,thisrow);
      AliHLTTPCTransform::Raw2Local(xyz,thissector,thisrow,pDigitReader->GetPad(),0);
            
      if(minSlice>17) xyz[1] = (-1.0)*xyz[1];
      else continue;
      
      AliHLTTPCTransform::Local2Global(xyz,minSlice);
      // temporarily the transformation Raw2Global will be broken down to 2 steps,
      // as there is a correction necessary at the y coordinate of the local xyz.
      
      //AliHLTTPCTransform::Raw2Global(xyz,thissector,thisrow,pDigitReader->GetPad(),0);
      // transformation from pad-row coordinates to global ones
      // time info is not taken into account
      
      //       AliTPCCalROC *calRoc = noisePad->GetCalROC(thissector);
      //       calRoc->GetValue(thisrow,pDigitReader->GetPad());
      
      while(pDigitReader->NextBunch()) {
    
	const UInt_t *bunchData = pDigitReader->GetSignals();
      	//Int_t time = pDigitReader->GetTime();
     
	for(Int_t i=0;i<pDigitReader->GetBunchSize();i++){
          
	  if((Float_t)(bunchData[i])>maxSignal){ maxSignal = (Float_t)(bunchData[i]); }
	  totalSignal   += (Float_t)bunchData[i];
	  squaredSignal += (Float_t)bunchData[i]*(Float_t)bunchData[i];
	  //fHistSignal->Fill(time+i, bunchData[i]);
	} // end for loop

	rms = TMath::Sqrt(squaredSignal/pDigitReader->GetBunchSize());
            
      } // end of inner while loop           
     
     
      if(minSlice<18){
	if(maxSignal>0)
	  fHistSideAMaxSignal->Fill(xyz[0],xyz[1],maxSignal);	  
	if(totalSignal>0)
	  fHistSideATotSignal->Fill(xyz[0],xyz[1],totalSignal);     
	if(rms>0)
	  fHistSideAPadRMS->Fill(xyz[0],xyz[1],rms);
      } else if(minSlice>17){
	if(maxSignal>0)
	  fHistSideCMaxSignal->Fill(xyz[0],xyz[1],maxSignal);
	if(totalSignal>0)
	  fHistSideCTotSignal->Fill(xyz[0],xyz[1],totalSignal);	
	if(rms>0)
	  fHistSideCPadRMS->Fill(xyz[0],xyz[1],rms);
      } else continue;
      
      maxSignal     = 0.;
      totalSignal   = 0.;
      squaredSignal = 0.;
      rms = 0.; 
    } // end of while loop over pads
    
    pDigitReader->Reset();
    delete pDigitReader;
  } // end of for loop over data blocks
 
  if(fResetHistograms) ResetHistograms(); fResetHistograms = kFALSE;
  MakeHistosPublic();
 
  return 0;
} // end DoEvent()


void AliHLTTPCNoiseMapComponent::MakeHistosPublic() {
// see header file for class documentation
 
//   TFile *outputfile = new TFile("test.root","RECREATE");
//   fHistSignal->Write();
//   outputfile->Save();
//   outputfile->Close();

  TObjArray histos;
  histos.Add(fHistSideAMaxSignal);
  histos.Add(fHistSideATotSignal);
  histos.Add(fHistSideAPadRMS);
  
  histos.Add(fHistSideCMaxSignal);
  histos.Add(fHistSideCTotSignal);
  histos.Add(fHistSideCPadRMS);
  
  histos.Add(fHistCDBMap);
  //histos.Add(fHistSignal);
  
  TIter iterator(&histos);
  while(TObject *pObj=iterator.Next()){ PushBack(pObj, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC, fSpecification); }
    
  //PushBack( (TObject*) &histos, kAliHLTDataTypeHistogram, fSpecification);    
}

void AliHLTTPCNoiseMapComponent::ResetHistograms(){
// see header file for class documentation

  if(fHistSideAMaxSignal) fHistSideAMaxSignal->Reset();
  if(fHistSideATotSignal) fHistSideATotSignal->Reset();
  if(fHistSideAPadRMS)    fHistSideAPadRMS->Reset();
 
  if(fHistSideCMaxSignal) fHistSideCMaxSignal->Reset();
  if(fHistSideCTotSignal) fHistSideCTotSignal->Reset();
  if(fHistSideCPadRMS)	  fHistSideCPadRMS->Reset();
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
     
      if (argument.CompareTo("-read-noisemap")==0) {
	//if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	fReadNoiseMap = kTRUE;
	HLTInfo("got \'-read-noisemap\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if(argument.CompareTo("-reset-histograms")==0){
	//if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	fResetHistograms = kTRUE;
	HLTInfo("got \'-reset-histograms\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());	
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

int AliHLTTPCNoiseMapComponent::Reconfigure(const char* cdbEntry, const char* chainId) { 
// see header file for class documentation
  
  int iResult = 0;
  const char* path = "HLT/ConfigTPC/TPCNoiseMapComponent";
  const char* defaultNotify = "";
  if(cdbEntry){
      path          = cdbEntry;
      defaultNotify = "(manual operator entry)";
  }
  
  if(path){          
     HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>" );
     
     AliCDBPath argumentPath(path);
     AliCDBStorage *stor = AliCDBManager::Instance()->GetDefaultStorage();
        
     if(stor){
        Int_t version    = stor->GetLatestVersion(path, GetRunNo());
        Int_t subVersion = stor->GetLatestSubVersion(path, GetRunNo(), version);
        AliCDBEntry *pEntry = stor->Get(argumentPath,GetRunNo(), version, subVersion);

        if(pEntry){
            TObjString* pString = dynamic_cast<TObjString*>(pEntry->GetObject());
            if(pString){
               HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
               iResult = Configure(pString->GetString().Data());
            } // if pString is valid
	    else {
               HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
            }
        } // if pEntry is valid
	else {
           HLTError("cannot fetch object \"%s\" from CDB", path);
        }
     } // end if stor is valid
  } // end if path is valid
  
  return iResult;


}
