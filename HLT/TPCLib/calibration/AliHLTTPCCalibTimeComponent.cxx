// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTTPCCalibTimeComponent.cxx
    @author Kalliopi Kanaki
    @date   2009-07-08
    @brief  A calibration component for interfacing the offline calculation of TPC drift velocity correction 
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCCalibTimeComponent.h"

#include "AliHLTTPCDefinitions.h"
#include "AliTPCcalibTime.h"
#include "AliESDEvent.h"

#include "TObjArray.h"
#include "TString.h"

#include <cstdlib>
#include <cerrno>


ClassImp(AliHLTTPCCalibTimeComponent) // ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCCalibTimeComponent::AliHLTTPCCalibTimeComponent()
  :
  fCalibTime(NULL),
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

AliHLTTPCCalibTimeComponent::~AliHLTTPCCalibTimeComponent() {
// see header file for class documentation
}


const char* AliHLTTPCCalibTimeComponent::GetComponentID() {
// see header file for class documentation

  return "TPCCalibTime";
}

void AliHLTTPCCalibTimeComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
// see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeTObjArray ); // output of TPCCalibSeedMaker
  list.push_back( kAliHLTDataTypeESDObject ); // output of TPCEsdConverter
}

AliHLTComponentDataType AliHLTTPCCalibTimeComponent::GetOutputDataType() {
// see header file for class documentation

  return AliHLTTPCDefinitions::fgkCalibCEDataType;
}

void AliHLTTPCCalibTimeComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
// see header file for class documentation

  constBase = 0;
  inputMultiplier = (2.0); // to be estimated
}

AliHLTComponent* AliHLTTPCCalibTimeComponent::Spawn() {
// see header file for class documentation

  return new AliHLTTPCCalibTimeComponent();
}  


Int_t AliHLTTPCCalibTimeComponent::ScanArgument( Int_t argc, const char** argv ) {
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

Int_t AliHLTTPCCalibTimeComponent::InitCalibration() {
// see header file for class documentation
    
  if(fCalibTime) return EINPROGRESS;
  //fCalibTime = new AliTPCcalibTime();

  fCalibTime = new AliTPCcalibTime("cosmicTime","cosmicTime",0, 1213.9e+06, 1213.96e+06);
  //AliTPCcalibTime::AliTPCcalibTime(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeVdrift)
 
  return 0;
}

Int_t AliHLTTPCCalibTimeComponent::DeinitCalibration() {
// see header file for class documentation

  if(fCalibTime) delete fCalibTime;
  fCalibTime = NULL;

  return 0;
}

Int_t AliHLTTPCCalibTimeComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
// see header file for class documentation
  
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;

  TObject *iterESD, *iterSEED = NULL;

  //----------- loop over output of TPCEsdConverter ----------------//

  for(iterESD = (TObject*)GetFirstInputObject(kAliHLTDataTypeESDObject|kAliHLTDataOriginTPC); iterESD != NULL; iterESD = (TObject*)GetNextInputObject()){   
       
      if(GetDataType(iterSEED) != (kAliHLTDataTypeESDObject | kAliHLTDataOriginTPC)) continue;
  
      AliHLTUInt8_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(GetSpecification(iterESD)); 
      AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(GetSpecification(iterESD));
      
      if( partition < fMinPartition ) fMinPartition = partition;
      if( partition > fMaxPartition ) fMaxPartition = partition;
      if( slice < fMinSlice ) fMinSlice = slice;
      if( slice > fMaxSlice ) fMaxSlice = slice;
      
      fESDEvent = dynamic_cast<AliESDEvent*>(iterESD);
      fESDEvent->CreateStdContent();
  } 
 
 
  //--------------- output over TObjArray of AliTPCseed objects (output of TPCSeedMaker) -------------------//
  
  for(iterSEED = (TObject*)GetFirstInputObject(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC); iterSEED != NULL; iterSEED = (TObject*)GetNextInputObject()){  
              
      if(GetDataType(iterSEED) != (kAliHLTDataTypeTObjArray | kAliHLTDataOriginTPC)) continue;
       
      fSeedArray = dynamic_cast<TObjArray*>(iterSEED); 
  }
 
  fCalibTime->Process(fESDEvent);
 
 
 
  fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( fMinSlice, fMaxSlice, fMinPartition, fMaxPartition );
  PushBack( (TObject*)fCalibTime, AliHLTTPCDefinitions::fgkCalibCEDataType, fSpecification);
  
  return 0;
}

Int_t AliHLTTPCCalibTimeComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
// see header file for class documentation
    
  if(fEnableAnalysis) fCalibTime->Analyze();  
  PushToFXS( (TObject*)fCalibTime, "TPC", "Time" ) ;
  
  return 0;
} 


