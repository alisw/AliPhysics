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

/** @file   AliHLTJETFastJetComponent.cxx
    @author Jochen Thaeder <thaeder@kip.uni-heidelberg.de>
    @date   
    @brief   Component to run the FastJet jetfinder
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include <cstdlib>
#include <cerrno>
#include <sys/time.h>

#include "AliHLTJETFastJetComponent.h" 

//#include "AliJetESDReader.h"
//#include "AliJetESDReaderHeader.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"

#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliStack.h"

#include "TString.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETFastJetComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTJETFastJetComponent::AliHLTJETFastJetComponent() 
  :
  fJetFinder(NULL),
  fJetHeader(NULL),
  fJetReader(NULL),
  fJetReaderHeader(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// #################################################################################
AliHLTJETFastJetComponent::~AliHLTJETFastJetComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTJETFastJetComponent::GetComponentID() {
  // see header file for class documentation
  return "JETFastJetFinder";
}

// #################################################################################
void AliHLTJETFastJetComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeMCObject|kAliHLTDataOriginOffline );
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginOffline );
  list.push_back( kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT );
}

// #################################################################################
AliHLTComponentDataType AliHLTJETFastJetComponent::GetOutputDataType() {
  // see header file for class documentation
  return (kAliHLTDataTypeESDObject| kAliHLTDataOriginHLT);
}

// #################################################################################
void AliHLTJETFastJetComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation

  constBase = 0;
  inputMultiplier = 0.3;
}

// #################################################################################
AliHLTComponent* AliHLTJETFastJetComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTJETFastJetComponent();
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETFastJetComponent::DoInit( Int_t /*argc*/, const Char_t** /*argv*/ ) {
  // see header file for class documentation

  if ( fJetFinder || fJetReader || fJetHeader || fJetReader )
    return -EINPROGRESS;


  // -- Jet Reader Header
  // -------------------------------------------
  if ( ! (fJetReaderHeader = new AliJetKineReaderHeader()) ) {
    HLTError("Error initializing Jet Reader Header");
    return -EINPROGRESS;
  }
  
  fJetReaderHeader->SetComment("MC full Kinematics");
  fJetReaderHeader->SetFastSimTPC(kFALSE);
  fJetReaderHeader->SetFastSimEMCAL(kFALSE);
  fJetReaderHeader->SetPtCut(0.);

  // -- Jet Reader
  // -------------------------------------------
  if ( ! (fJetReader = new AliJetKineReader()) ) {
    HLTError("Error initializing Jet Reader");
    return -EINPROGRESS;
  }

  fJetReader->SetReaderHeader(fJetReaderHeader);

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


#if 0
  fJetReaderHeader = new AliJetESDReaderHeader();
  fJetReaderHeader->SetComment("Testing");
  fJetReaderHeader->SetFirstEvent(0);
  fJetReaderHeader->SetLastEvent(4);
  fJetReader = new AliJetESDReader();
  fJetReader->SetReaderHeader(fJetReaderHeader);
#endif
  
  return 0;
}

// #################################################################################
Int_t AliHLTJETFastJetComponent::DoDeinit() {
  // see header file for class documentation

  if ( fJetFinder )
    delete fJetFinder;
  fJetFinder = NULL;

  if ( fJetHeader )
    delete fJetHeader;
  fJetHeader = NULL;
 
  if ( fJetReader )
    delete fJetReader;
  fJetReader = NULL;
 
  if ( fJetReaderHeader )
    delete fJetReaderHeader;
  fJetReaderHeader = NULL;

  return 0;
}

// #################################################################################
Int_t AliHLTJETFastJetComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
					  AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation

  const TObject* iter = NULL;
  
  for ( iter=GetFirstInputObject(kAliHLTDataTypeMCObject|kAliHLTDataOriginOffline); iter != NULL; iter=GetNextInputObject() ) {

    // ADD MC Object
    
    // -- Set Input Event
    fJetFinder->GetReader()->SetInputEvent( NULL, NULL, const_cast<TObject*>(iter) );    

    AliMCEvent* foo = ( AliMCEvent* ) iter;
    cout << foo->GetNumberOfTracks() << " -- "
         << foo->Stack() << " -- "
         << foo->Header()->Stack() << " -- "
         << foo->Header() << endl;

    AliStack* stack =  foo->Stack();
    cout << "N tracks" << stack->GetNtrack() << endl;

    for (Int_t iterStack = 0; iterStack < stack->GetNtrack(); iterStack++) {
      cout << iterStack << " -- " << stack->Particle(iterStack) << endl;

    }
    // -- Process Event
    //    fJetFinder->ProcessEvent();
  }

  
  if ( GetFirstInputObject(kAliHLTDataTypeEOR) ) {
    //fJetFinder->FinishRun();
  }



  for ( iter=GetFirstInputObject(kAliHLTDataTypeESDObject|kAliHLTDataOriginOffline); iter != NULL; iter=GetNextInputObject() ) {
    // ADD ESD Object -- Offline

    printf ("  ---   ESD-Offline  ---  \n");
  }

  for ( iter=GetFirstInputObject(kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT); iter != NULL; iter=GetNextInputObject() ) {
    // ADD ESD Object -- HLT
    printf ("  ---   ESD-HLT  ---  \n");
  }

 
  // ** PushBack ** \\


  return 0;
}
