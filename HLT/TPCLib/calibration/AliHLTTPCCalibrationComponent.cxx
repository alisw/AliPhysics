// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
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

/** @file   AliHLTTPCCalibrationComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCCalibrationComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCAnalysisTaskcalib.h"
#include "AliHLTReadoutList.h"

#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
//#include "AliESDInputHandler.h"

#include "AliTPCcalibTime.h"
#include "AliTPCcalibTimeGain.h"
#include "AliTPCseed.h"

#include "TString.h"
#include "TObjArray.h"
#include "TTimeStamp.h"

#include <cstdlib>
#include <cerrno>

ClassImp(AliHLTTPCCalibrationComponent) // ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCCalibrationComponent::AliHLTTPCCalibrationComponent()
  :
  fCalibTask(NULL),
  fCalibTime(NULL),
  fCalibTimeGain(NULL),
  fESDEvent(NULL),
  fSeedArray(NULL),
  fMinPartition(5),
  fMaxPartition(0),
  fMinSlice(35),
  fMaxSlice(0),
  fSpecification(0) ,
  fEnableAnalysis(kTRUE)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCCalibrationComponent::~AliHLTTPCCalibrationComponent() {
// see header file for class documentation
}

const char* AliHLTTPCCalibrationComponent::GetComponentID() {
// see header file for class documentation

  return "TPCCalibration";
}

void AliHLTTPCCalibrationComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
// see header file for class documentation

  list.clear();     
  list.push_back( kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC );  // output of TPCCalibSeedMaker
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut ); // output of global esd converter 

}

AliHLTComponentDataType AliHLTTPCCalibrationComponent::GetOutputDataType() {
// see header file for class documentation

  return AliHLTTPCDefinitions::fgkCalibCEDataType;
}

void AliHLTTPCCalibrationComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
// see header file for class documentation

  constBase = 0;
  inputMultiplier = (2.0); // to be estimated
}

AliHLTComponent* AliHLTTPCCalibrationComponent::Spawn() {
// see header file for class documentation

  return new AliHLTTPCCalibrationComponent();
}  


Int_t AliHLTTPCCalibrationComponent::ScanArgument( Int_t argc, const char** argv ) {
// see header file for class documentation

  Int_t iResult = 0;
  TString argument = "";
  TString parameter = "";

  if(!argc) return -EINVAL;

  argument = argv[iResult];
  
  if(argument.IsNull()) return -EINVAL;

  if( argument.CompareTo("-enable-analysis") == 0 ){
    HLTInfo( "Analysis before shipping data to FXS enabled." );
    fEnableAnalysis = kTRUE;
  }
  else {
    iResult = -EINVAL;
  }
      
  return iResult;
}

Int_t AliHLTTPCCalibrationComponent::InitCalibration() {
// see header file for class documentation
  
  if(fCalibTask) return EINPROGRESS;
  fCalibTask = new AliHLTTPCAnalysisTaskcalib("TPC Calibration Task");
  
  if(fCalibTime) return EINPROGRESS;
  //fCalibTime = new AliTPCcalibTime();
  fCalibTime = new AliTPCcalibTime("calibTime","time dependent Vdrift calibration",-2, 2, 1);
  
  fCalibTime->SetDebugLevel(20);
  fCalibTime->SetStreamLevel(10);
  fCalibTime->SetTriggerMask(-1,-1,kFALSE); //accept everything 

  if(fCalibTimeGain) return EINPROGRESS;
  fCalibTimeGain = new AliTPCcalibTimeGain("calibTimeGain","time dependent gain calibration",-2, 2, 1);
  
  fCalibTask->AddJob(fCalibTime);
  fCalibTask->AddJob(fCalibTimeGain);
  fCalibTask->GetJobs();
 
  return 0;
}

Int_t AliHLTTPCCalibrationComponent::DeinitCalibration() {
// see header file for class documentation

  if(fCalibTask)     delete fCalibTask;     fCalibTask     = NULL;
  if(fCalibTime)     delete fCalibTime;     fCalibTime     = NULL;
  if(fCalibTimeGain) delete fCalibTimeGain; fCalibTimeGain = NULL;

  return 0;
}

Int_t AliHLTTPCCalibrationComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
// see header file for class documentation
  
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;

  TObject *iter = NULL;

  //--------------- output over TObjArray of AliTPCseed objects (output of TPCSeedMaker) -------------------//
  
  for(iter = (TObject*)GetFirstInputObject(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC); iter != NULL; iter = (TObject*)GetNextInputObject()){  
              
      if(GetDataType(iter) != (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTPC)) continue;
      fSeedArray = dynamic_cast<TObjArray*>(iter); 
  }
  
  //----------- loop over output of global esd converter ----------------//
  
  for(iter = (TObject*)GetFirstInputObject(kAliHLTDataTypeESDObject | kAliHLTDataOriginOut); iter != NULL; iter = (TObject*)GetNextInputObject()){   
      
      if(GetDataType(iter) != (kAliHLTDataTypeESDObject | kAliHLTDataOriginOut)) continue;
  
      AliHLTUInt8_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(GetSpecification(iter)); 
      AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(GetSpecification(iter));
      
      if( partition < fMinPartition ) fMinPartition = partition;
      if( partition > fMaxPartition ) fMaxPartition = partition;
      if( slice < fMinSlice ) fMinSlice = slice;
      if( slice > fMaxSlice ) fMaxSlice = slice;
      
      fESDEvent = dynamic_cast<AliESDEvent*>(iter);
      fESDEvent->CreateStdContent();
     
      HLTDebug("# Seeds: %i\n", fSeedArray->GetEntriesFast());
      
      for(Int_t i=0; i<fSeedArray->GetEntriesFast(); i++){
          
	  AliTPCseed *seed = (AliTPCseed*)fSeedArray->UncheckedAt(i);
          if(!seed) continue;
	  AliESDtrack *esd = fESDEvent->GetTrack(i);	
	  AliTPCseed *seedCopy = new AliTPCseed(*seed, kTRUE); 
	  esd->AddCalibObject(seedCopy);
      }      
 }   
  
  fCalibTask->Process(fESDEvent);

  fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( fMinSlice, fMaxSlice, fMinPartition, fMaxPartition );
  PushBack( (TObject*) fCalibTask, AliHLTTPCDefinitions::fgkCalibCEDataType | kAliHLTDataOriginOut, fSpecification);
  
  return 0;
}

Int_t AliHLTTPCCalibrationComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
  // see header file for class documentation
    
  if(fEnableAnalysis) fCalibTask->Analyze();  
  static AliHLTReadoutList rdList(AliHLTReadoutList::kTPC);
  PushToFXS( (TObject*) fCalibTask, "TPC", "CALIB", &rdList ) ;
  
  return 0;
} 


