//-*- Mode: C++ -*-
// $Id: AliHLTJETAnalysisComponent.cxx $
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

/** @file   AliHLTJETAnalysisComponent.cxx
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

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "AliHLTJETAnalysisComponent.h" 

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTJETAnalysisComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTJETAnalysisComponent::AliHLTJETAnalysisComponent() :
  fAnalysisJets(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

// #################################################################################
AliHLTJETAnalysisComponent::~AliHLTJETAnalysisComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTJETAnalysisComponent::GetComponentID() {
  // see header file for class documentation
  return "JETAnalysis";
}

// #################################################################################
void AliHLTJETAnalysisComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.clear(); 
  list.push_back( kAliHLTDataTypeJet|kAliHLTDataOriginHLT );
  list.push_back( kAliHLTDataTypeJet|kAliHLTDataOriginOffline );
  list.push_back( kAliHLTDataTypeMCObject|kAliHLTDataOriginOffline );
  list.push_back( kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT );
}

// #################################################################################
AliHLTComponentDataType AliHLTJETAnalysisComponent::GetOutputDataType() {
  // see header file for class documentation
  return (kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT); // XXX To Be CHANGED
}

// #################################################################################
void AliHLTJETAnalysisComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation

  constBase = 20000000;
  inputMultiplier = 0.3;
}

// #################################################################################
AliHLTComponent* AliHLTJETAnalysisComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTJETAnalysisComponent();
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTJETAnalysisComponent::DoInit( Int_t argc, const Char_t** argv ) {
  // see header file for class documentation

  Int_t iResult = 0;

  if ( fAnalysisJets )
    return -EINPROGRESS;

  Int_t bMissingParam=0;
  TString argument="";

  // -- Loop over all arguments
  for ( Int_t iter = 0; iter<argc && iResult>=0; iter++) {
    argument=argv[iter];

    if (argument.IsNull()) 
      continue;

#if 0
    /*
    // -- nEvents
    if ( !argument.CompareTo("-nevents") ) {
    if ((bMissingParam=(++iter>=argc))) break;
    
    TString parameter(argv[iter]);
    parameter.Remove(TString::kLeading, ' ');
    
    if ( parameter.IsDigit() ) {
    fNEventsMax = parameter.Atoi() - 1;
    }
    else {
    HLTError("Wrong parameter %s for argument %s.", parameter.Data(), argument.Data());
    iResult=-EINVAL;
    }
    } 
    
    // -- Argument not known
    else {
    HLTError("Unknown argument %s.", argument.Data());
    iResult = -EINVAL;
    }
    */
#endif
  } // for ( Int iter = 0; iter<argc && iResult>=0; iter++) {
  
  // -- Check if parameter is missing
  if ( bMissingParam ) {
    HLTError("Missing parameter for argument %s.", argument.Data());
    iResult=-EINVAL;
  }

  if ( !iResult ) {
    fAnalysisJets = new AliHLTJETAnalysisJets();
    fAnalysisJets->Initialize();
  }
  
  return iResult;
}

// #################################################################################
Int_t AliHLTJETAnalysisComponent::DoDeinit() {
  // see header file for class documentation

  if ( fAnalysisJets )
    delete fAnalysisJets;
  fAnalysisJets = NULL;

  return 0;
}

// #################################################################################
Int_t AliHLTJETAnalysisComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
					  AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation

  Int_t iResult = 0;

  const TObject* iter = NULL;

  // ------------------------------------------------
  // -- DATA Event 
  // ------------------------------------------------
  if ( IsDataEvent() ) {

    fAnalysisJets->ResetEvent();
    
    // -- Loop over jet objects
    // ------------------------------
    for ( iter=GetFirstInputObject(kAliHLTDataTypeJet|kAliHLTDataOriginHLT); 
	  iter != NULL && !iResult; iter=GetNextInputObject() ) {
      
      fAnalysisJets->SetJets(reinterpret_cast<AliHLTJETJets*>(const_cast<TObject*>(iter)));
    }
    
    // -- ADD MC Object -- On-line
    // ------------------------------
    for ( iter=GetFirstInputObject(kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT); 
	  iter != NULL && !iResult; iter=GetNextInputObject() ) {
      
      fAnalysisJets->SetHLTMC(reinterpret_cast<AliHLTMCEvent*>(const_cast<TObject*>(iter)));
    }
        
    // -- Process event
    // ------------------
    iResult = fAnalysisJets->Analyze();
  }

  // ------------------------------------------------
  // -- DATA Event 
  // ------------------------------------------------
  else {
    if ( GetFirstInputBlock(kAliHLTDataTypeEOR) ) {
      PushBack(fAnalysisJets, kAliHLTDataTypeJet|kAliHLTDataOriginHLT, GetSpecification());
    }
  }
   
  return iResult;
}
