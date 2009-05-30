//-*- Mode: C++ -*-
// $Id: $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTJETConeJetComponent.cxx
    @author Jochen Thaeder <thaeder@kip.uni-heidelberg.de>
    @date   
    @brief   Component to run the ConeJet jetfinder
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <cstdlib>
#include <cerrno>
#include <sys/time.h>

#include "AliHLTJETConeJetComponent.h" 

#include "TString.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETConeJetComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTJETConeJetComponent::AliHLTJETConeJetComponent() 
  :
  //fJetFinder(NULL),
  //  fJetHeader(NULL),*/
  fJetReader(NULL),
  fJetReaderHeader(NULL),
  fJetTrackCuts(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETConeJetComponent::~AliHLTJETConeJetComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTJETConeJetComponent::GetComponentID() {
  // see header file for class documentation
  return "JETConeJetFinder";
}

// #################################################################################
void AliHLTJETConeJetComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeMCObject|kAliHLTDataOriginOffline );
  list.push_back( kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT );
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOffline );
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT );
}

// #################################################################################
AliHLTComponentDataType AliHLTJETConeJetComponent::GetOutputDataType() {
  // see header file for class documentation
  return (kAliHLTDataTypeESDObject| kAliHLTDataOriginHLT);
}

// #################################################################################
void AliHLTJETConeJetComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation

  constBase = 0;
  inputMultiplier = 0.3;
}

// #################################################################################
AliHLTComponent* AliHLTJETConeJetComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTJETConeJetComponent();
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETConeJetComponent::DoInit( Int_t /*argc*/, const Char_t** /*argv*/ ) {
  // see header file for class documentation

  if ( /*fJetFinder || fJetHeader ||*/ fJetReader || fJetReader || fJetTrackCuts)
    return -EINPROGRESS;


  // -- Jet Track Cuts
  // -------------------------------------------
  if ( ! (fJetTrackCuts = new AliHLTJETTrackCuts()) ) {
    HLTError("Error initializing Track Cuts");
    return -EINPROGRESS;
  }

  // fJetTrackCuts->Set ...
  
  // -- Jet Reader Header
  // -------------------------------------------
  if ( ! (fJetReaderHeader = new AliHLTJETReaderHeader()) ) {
    HLTError("Error initializing Jet Reader Header");
    return -EINPROGRESS;
  }

  fJetReaderHeader->SetAnalysisCuts( dynamic_cast<AliAnalysisCuts*>(fJetTrackCuts) );

  // -- Jet Reader
  // -------------------------------------------
  if ( ! (fJetReader = new AliHLTJETReader()) ) {
    HLTError("Error initializing Jet Reader");
    return -EINPROGRESS;
  }

  fJetReader->SetReaderHeader(fJetReaderHeader);

#if 0
  // -- Jet Header
  // -------------------------------------------
  if ( ! (fJetHeader = new AliFastJetHeader()) ) {
    HLTError("Error initializing Jet Header");
    return -EINPROGRESS;
  }

  fJetHeader->SetRparam(0.7); 

  // -- Jet Finder
  // -------------------------------------------
  if ( ! (fJetFinder = new AliFastJetFinder()) ) {
    HLTError("Error initializing Jet Finder");
    return -EINPROGRESS;
  }

  fJetFinder->SetJetHeader(fJetHeader);
  fJetFinder->SetJetReader(fJetReader);
  fJetFinder->SetOutputFile("jets.root");

  // -- Initialize Jet Finder
  // -------------------------------------------
  fJetFinder->Init();
#endif


  return 0;
}

// #################################################################################
Int_t AliHLTJETConeJetComponent::DoDeinit() {
  // see header file for class documentation

  /*
  if ( fJetFinder )
    delete fJetFinder;
  fJetFinder = NULL;

  if ( fJetHeader )
    delete fJetHeader;
  fJetHeader = NULL;
  */
  /*
  if ( fJetReader )
    delete fJetReader;
  fJetReader = NULL;

  if ( fJetReaderHeader )
    delete fJetReaderHeader;
  fJetReaderHeader = NULL;

  if ( fJetTrackCuts )
    delete fJetTrackCuts;
  fJetTrackCuts = NULL;



  */
  return 0;
}

// #################################################################################
Int_t AliHLTJETConeJetComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
					  AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation

  Int_t iResult = 0;

  const TObject* iter = NULL;

  // -- Start-Of-Run
  // -----------------
  if ( GetFirstInputObject(kAliHLTDataTypeSOR) && !iResult ) {
    HLTInfo("On-line SOR Event");
  }

  // -- ADD MC Object -- Off-line
  // ------------------------------
  for ( iter=GetFirstInputObject(kAliHLTDataTypeMCObject|kAliHLTDataOriginOffline); iter != NULL && !iResult; iter=GetNextInputObject() ) {
    HLTInfo("Off-line MC Event");
  }

  // -- ADD MC Object -- On-line
  // ------------------------------
  for ( iter=GetFirstInputObject(kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT); iter != NULL && !iResult; iter=GetNextInputObject() ) {
    HLTInfo("On-line MC Event");
    
    // -- Set input event
    fJetReader->SetInputEvent( NULL, NULL, const_cast<TObject*>(iter) );    

    // -- FillMomentumArrayMC
    if ( ! (fJetReader->FillMomentumArrayMC()) )
      iResult = -1;  

    // -- Process Event
    // fJetFinder->ProcessEvent();

  }

  // -- ADD ESD Object -- Off-line
  // -------------------------------
  for ( iter=GetFirstInputObject(kAliHLTDataTypeESDObject|kAliHLTDataOriginOffline); iter != NULL && !iResult; iter=GetNextInputObject() ) {
    HLTInfo("Off-line ESD Event");
  }

  // -- ADD ESD Object -- On-line
  // ------------------------------
  for ( iter=GetFirstInputObject(kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT); iter != NULL && !iResult; iter=GetNextInputObject() ) {
    HLTInfo("On-line ESD Event");
  }

  // -- End-Of-Run
  // ---------------
  if ( GetFirstInputObject(kAliHLTDataTypeEOR) && !iResult ) {
    HLTInfo("On-line EOR Event");
    
    // -- Finish Event ?
    // fJetFinder->FinishRun();
  }
      
  // -- PushBack
  // -------------
  
  return 0;
}
