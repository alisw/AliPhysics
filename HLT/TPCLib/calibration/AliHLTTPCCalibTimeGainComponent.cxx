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

/** @file   AliHLTTPCCalibTimeGainComponent.cxx
    @author Kalliopi Kanaki
    @date   2009-07-08
    @brief  A calibration component for the TPC gain variation vs. time. 
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCCalibTimeGainComponent.h"

#include "AliHLTTPCDefinitions.h"
#include "AliTPCcalibTimeGain.h"
#include "AliESDEvent.h"
#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjArray.h"

ClassImp(AliHLTTPCCalibTimeGainComponent) // ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCCalibTimeGainComponent::AliHLTTPCCalibTimeGainComponent()
  :
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

AliHLTTPCCalibTimeGainComponent::~AliHLTTPCCalibTimeGainComponent() {
// see header file for class documentation
}

const char* AliHLTTPCCalibTimeGainComponent::GetComponentID() {
// see header file for class documentation

  return "TPCCalibTimeGain";
}

void AliHLTTPCCalibTimeGainComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
// see header file for class documentation

  list.clear();   
  list.push_back( kAliHLTDataTypeTObjArray ); // output of TPCCalibSeedMaker
  list.push_back( kAliHLTDataTypeESDObject ); // output of TPCEsdConverter 
}

AliHLTComponentDataType AliHLTTPCCalibTimeGainComponent::GetOutputDataType() {
// see header file for class documentation

  return AliHLTTPCDefinitions::fgkCalibCEDataType;
}

void AliHLTTPCCalibTimeGainComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
// see header file for class documentation

  constBase = 0;
  inputMultiplier = (2.0); // to be estimated
}

AliHLTComponent* AliHLTTPCCalibTimeGainComponent::Spawn() {
// see header file for class documentation

  return new AliHLTTPCCalibTimeGainComponent();
}  


Int_t AliHLTTPCCalibTimeGainComponent::ScanArgument( Int_t argc, const char** argv ) {
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

Int_t AliHLTTPCCalibTimeGainComponent::InitCalibration() {
// see header file for class documentation
    
  if(fCalibTimeGain) return EINPROGRESS;
  fCalibTimeGain = new AliTPCcalibTimeGain();
  //AliTPCcalibTimeGain::AliTPCcalibTimeGain(const Text_t *name, const Text_t *title, UInt_t StartTime, UInt_t EndTime, Int_t deltaIntegrationTimeGain)
  return 0;
}

Int_t AliHLTTPCCalibTimeGainComponent::DeinitCalibration() {
// see header file for class documentation

  if(fCalibTimeGain) delete fCalibTimeGain;
  fCalibTimeGain = NULL;

  return 0;
}

Int_t AliHLTTPCCalibTimeGainComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
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
 
  fCalibTimeGain->Process(fESDEvent);

  fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( fMinSlice, fMaxSlice, fMinPartition, fMaxPartition );
  PushBack( (TObject*) fCalibTimeGain, AliHLTTPCDefinitions::fgkCalibCEDataType, fSpecification);
  
  return 0;
}

Int_t AliHLTTPCCalibTimeGainComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
  // see header file for class documentation
    
  if(fEnableAnalysis) fCalibTimeGain->Analyze();  
  PushToFXS( (TObject*) fCalibTimeGain, "TPC", "TimeGain" ) ;
  
  return 0;
} 


