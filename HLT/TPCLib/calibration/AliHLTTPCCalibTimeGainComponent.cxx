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
#include "AliHLTReadoutList.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"

#include "AliTPCcalibTimeGain.h"
#include "AliTPCseed.h"

#include "TObjArray.h"
#include "TString.h"

#include <cstdlib>
#include <cerrno>


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
  list.push_back( kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC ); // output of TPCCalibSeedMaker
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOut ); // output of global esd converter 
}

AliHLTComponentDataType AliHLTTPCCalibTimeGainComponent::GetOutputDataType() {
// see header file for class documentation

  return AliHLTTPCDefinitions::fgkCalibCEDataType|kAliHLTDataOriginOut;
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
  //fCalibTimeGain = new AliTPCcalibTimeGain();
  fCalibTimeGain = new AliTPCcalibTimeGain("calibTimeGain","time dependent gain calibration",-2, 2, 1);
  return 0;
}

Int_t AliHLTTPCCalibTimeGainComponent::DeinitCalibration() {
// see header file for class documentation

  if(fCalibTimeGain) delete fCalibTimeGain; fCalibTimeGain = NULL;

  return 0;
}

Int_t AliHLTTPCCalibTimeGainComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
// see header file for class documentation
  
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR )) return 0;

  TObject *iter = NULL;

  //--------------- loop over the TObjArray output of AliTPCseed objects (output of TPCSeedMaker) -------------------//
  
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
   
  fCalibTimeGain->Process(fESDEvent);

  fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( fMinSlice, fMaxSlice, fMinPartition, fMaxPartition );
  PushBack( (TObject*) fCalibTimeGain, AliHLTTPCDefinitions::fgkCalibCEDataType| kAliHLTDataOriginOut, fSpecification);
  
  return 0;
}

Int_t AliHLTTPCCalibTimeGainComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,  AliHLTComponentTriggerData& /*trigData*/ ){
  // see header file for class documentation
    
  if(fEnableAnalysis) fCalibTimeGain->Analyze();  
  static AliHLTReadoutList rdList(AliHLTReadoutList::kTPC);
  PushToFXS( (TObject*) fCalibTimeGain, "TPC", "TimeGain", &rdList ) ;
  
  return 0;
} 


