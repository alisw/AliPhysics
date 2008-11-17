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
//#include "AliHLTTPCDigitReaderPacked.h"
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
    fPlotSideA(0),
    fPlotSideC(0),    
    fApplyNoiseMap(0),
    fResetHistograms(0),
    fIsPacked(0),
    fIsUnpacked(0),
    fCurrentSlice(-99),
    fCurrentPartition(-99),
    fCurrentRow(-99),
    fHistSignal(NULL),
    fHistMaxSignal(NULL),
    fHistTotSignal(NULL),
    fHistPadRMS(NULL),
    fHistCDBMap(NULL),    
    fHistSideA(NULL),  
    fHistSideC(NULL)
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
    if (!strcmp( argv[i], "-apply-noisemap")) {
        fApplyNoiseMap = strtoul( argv[i+1], &cpErr ,0);
            
    if ( *cpErr ) {
        HLTError("Cannot convert apply-noisemap specifier '%s'.", argv[i+1]);
        return EINVAL;
    }
      i+=2;
      continue;
    }
    
    if (!strcmp( argv[i], "-plot-side-a")) {
        fPlotSideA = strtoul( argv[i+1], &cpErr ,0);
            
    if ( *cpErr ) {
        HLTError("Cannot convert plot-side-a specifier '%s'.", argv[i+1]);
        return EINVAL;
    }
      i+=2;
      continue;
    }
    
    if (!strcmp( argv[i], "-plot-side-c")) {
        fPlotSideC = strtoul( argv[i+1], &cpErr ,0);
    
    if ( *cpErr ) {
        HLTError("Cannot convert plot-side-c specifier '%s'.", argv[i+1]);
        return EINVAL;
    }
      i+=2;
      continue;
    }

    if (!strcmp( argv[i], "-reset-histograms")) {
        fResetHistograms = strtoul( argv[i+1], &cpErr ,0);
    
    if ( *cpErr ) {
        HLTError("Cannot convert reset-histograms specifier '%s'.", argv[i+1]);
        return EINVAL;
    }
      i+=2;
      continue;
    }
                   
    Logging(kHLTLogError, "HLT::TPCNoiseMap::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  } // end while
  
  if(fApplyNoiseMap){
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
   
  if(fPlotSideA){
     fHistSideA = new TH2F("fHistSideA","TPC Side A",250,-250,250,250,-250,250);		
     fHistSideA->SetXTitle("global X (cm)"); fHistSideA->SetYTitle("global Y (cm)");
  }   
  
  if(fPlotSideC){    
     fHistSideC = new TH2F("fHistSideC","TPC Side C",250,-250,250,250,-250,250);
     fHistSideC->SetXTitle("global X (cm)"); fHistSideC->SetYTitle("global Y (cm)");
  }
 
  fHistMaxSignal = new TH2F("fHistMaxSignal","maximum signal",   250,-250,250,250,-250,250);
  fHistTotSignal = new TH2F("fHistTotSignal","total signal",     250,-250,250,250,-250,250);
  fHistPadRMS    = new TH2F("fHistPadRMS",   "RMS",              250,-250,250,250,-250,250);
  //fHistSignal    = new TH1F("fHistSignal", "signal distribution per pad",1024,0,1024);
 
//   HLTDebug("using AliHLTTPCDigitReaderDecoder");
//   pDigitReader = new AliHLTTPCDigitReaderDecoder(); // double-loop
//   pDigitReader = new AliHLTTPCDigitReaderPacked();
  
  return 0;

} // end DoInit()

int AliHLTTPCNoiseMapComponent::DoDeinit() { 
  // see header file for class documentation  

  if(fHistMaxSignal) delete fHistMaxSignal; fHistMaxSignal = NULL;
  if(fHistTotSignal) delete fHistTotSignal; fHistTotSignal = NULL;
  if(fHistPadRMS)    delete fHistPadRMS;    fHistPadRMS    = NULL;
  if(fHistSideA)     delete fHistSideA;     fHistSideA     = NULL;
  if(fHistSideC)     delete fHistSideC;     fHistSideC     = NULL;
       
  return 0;
}

int AliHLTTPCNoiseMapComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/){
  // see header file for class documentation
 
  //HLTInfo("--- Entering DoEvent() in TPCNoiseMap ---");
 
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;
   
  const AliHLTComponentBlockData *iter = NULL;

  Float_t xyz[3]; 
  Int_t thissector, thisrow;
    
  for(iter = GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC); iter != NULL; iter = GetNextInputBlock()){
      
    HLTInfo("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s", 
	    evtData.fEventID, evtData.fEventID,
	    DataType2Text(iter->fDataType).c_str(), 
	    DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());

    if (iter->fDataType == AliHLTTPCDefinitions::fgkDDLPackedRawDataType && GetEventCount()<2){
      HLTWarning("data type %s is depricated, use %s (kAliHLTDataTypeDDLRaw)!", 
		 DataType2Text(AliHLTTPCDefinitions::fgkDDLPackedRawDataType).c_str(),
		 DataType2Text(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC).c_str());
    }      
     
    if (iter->fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC)) continue;
      
    UInt_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter); 
    UInt_t partition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);
     
    fSpecification = iter->fSpecification;
     
    AliHLTTPCDigitReader *pDigitReader = new AliHLTTPCDigitReaderDecoder;

    pDigitReader->InitBlock(iter->fPtr,iter->fSize,partition,slice);
    if(!pDigitReader) break;
       
    //sprintf(name,"hMaxSignal_slice%d_partition%d", slice, partition);
    //fHistMaxSignal = new TH2F(name,name,250,-250,250,250,-250,250);
            
    //  while(pDigitReader->Next()){ 


    Float_t maxSignal     = 0.;
    Float_t totalSignal   = 0.;
    Float_t squaredSignal = 0.;
    Float_t rms = 0.; 
    
    while( pDigitReader->NextChannel()) { // pad loop 
      
      fCurrentRow  = pDigitReader->GetRow();  
      fCurrentRow += pDigitReader->GetRowOffset();

      AliHLTTPCTransform::Slice2Sector(slice,fCurrentRow,thissector,thisrow);
      AliHLTTPCTransform::Raw2Local(xyz,thissector,thisrow,pDigitReader->GetPad(),0);
      
      if(slice>17) xyz[1] = (-1.0)*xyz[1];
      else continue;
      
      AliHLTTPCTransform::Local2Global(xyz,slice);
      // temporarily the transformation Raw2Global will be broken down to 2 steps,
      // as there is a correction necessary at the y coordinate of the local xyz.
      
      //AliHLTTPCTransform::Raw2Global(xyz,thissector,thisrow,pDigitReader->GetPad(),0);
      // transformation from pad-row coordinates to global ones
      // time info is not taken into account
      
      //       AliTPCCalROC *calRoc = noisePad->GetCalROC(thissector);
      //       calRoc->GetValue(thisrow,pDigitReader->GetPad());
      
      
      while( pDigitReader->NextBunch()) {
    
	const UInt_t *bunchData = pDigitReader->GetSignals();
      
	//fHistSignal = new TH1F("fHistSignal", "signal distribution per pad",1024,0,1024);
      
	//fHistSignal->Reset();
	//Int_t time = pDigitReader->GetTime();
     
	for(Int_t i=0;i<pDigitReader->GetBunchSize();i++){
          
	  if((Float_t)(bunchData[i])>maxSignal){ maxSignal = (Float_t)(bunchData[i]); }
	  totalSignal += (Float_t)bunchData[i];
	  squaredSignal += (Float_t)bunchData[i]*(Float_t)bunchData[i];
	  //fHistSignal->Fill(time+i, bunchData[i]);
	} // end for loop over bunches
	rms = TMath::Sqrt(squaredSignal/pDigitReader->GetBunchSize());
            
      } // end of inner while loop
           
    } // end of while loop over pads
     
    fHistMaxSignal->Fill(xyz[0],xyz[1],maxSignal);
    fHistTotSignal->Fill(xyz[0],xyz[1],totalSignal);
     
    fHistPadRMS->Fill(xyz[0],xyz[1],rms);
      
    //fHistPadRMS->Fill(xyz[0],xyz[1],fHistSignal->GetRMS());
    //delete fHistSignal; fHistSignal = NULL;
      
    if(fPlotSideA || fPlotSideC){
      if(slice<18) fHistSideA->Fill(xyz[0],xyz[1],maxSignal);
      else	      fHistSideC->Fill(xyz[0],xyz[1],maxSignal);			     
    } // end if plotting sides    
      
    
    maxSignal     = 0.;
    totalSignal   = 0.;
    squaredSignal = 0.;
    rms = 0.; 

    pDigitReader->Reset();
    delete pDigitReader;
  } // end of for loop over data blocks
 
  if(fResetHistograms) ResetHistograms();
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
  histos.Add(fHistMaxSignal);
  histos.Add(fHistTotSignal);
  histos.Add(fHistPadRMS);
  histos.Add(fHistCDBMap);
  //histos.Add(fHistSignal);
  if(fPlotSideA) histos.Add(fHistSideA);
  if(fPlotSideC) histos.Add(fHistSideC);
  
  TIter iterator(&histos);
  while(TObject *pObj=iterator.Next()){ PushBack(pObj, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC, fSpecification); }
    
  //PushBack( (TObject*) &histos, kAliHLTDataTypeHistogram, fSpecification);    
 
//   if(fHistMaxSignal) delete fHistMaxSignal; fHistMaxSignal = NULL;
//   if(fHistTotSignal) delete fHistTotSignal; fHistTotSignal = NULL;
//   if(fHistPadRMS)    delete fHistPadRMS;    fHistPadRMS    = NULL;
//   if(fHistSideA)     delete fHistSideA;     fHistSideA     = NULL;
//   if(fHistSideC)     delete fHistSideC;     fHistSideC     = NULL;
  
}

void AliHLTTPCNoiseMapComponent::ResetHistograms(){
// see header file for class documentation

  //if(fHistPartition) fHistPartition->Reset();  
  if(fHistMaxSignal) fHistMaxSignal->Reset();
  if(fHistTotSignal) fHistTotSignal->Reset();
  if(fHistPadRMS)    fHistPadRMS->Reset();

  if(fHistSideA) fHistSideA->Reset();
  if(fHistSideC) fHistSideC->Reset();
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
     
      if (argument.CompareTo("-apply-noisemap")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-apply-noisemap\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if (argument.CompareTo("-plot-side-c")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-plot-side-c\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      } 
      else if (argument.CompareTo("-plot-side-a")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("got \'-plot-side-a\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	
      }       
      else if(argument.CompareTo("-reset-histograms")==0){
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
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
  int iResult=0;
  const char* path="HLT/ConfigTPC/TPCNoiseMapComponent";
  const char* defaultNotify="";
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
     } // if stor is valid
  } // if path is valid
  
  return iResult;


}
