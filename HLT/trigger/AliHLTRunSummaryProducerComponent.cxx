//-*- Mode: C++ -*-
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

/** @file   AliHLTRunSummaryProducerComponent.cxx
    @author Jochen Thaeder
    @date   
    @brief  Produces a run summary as @see AliHLTRunSummary
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTRunSummaryProducerComponent.h"
#include "AliHLTEventSummary.h"
#include "AliHLTDataTypes.h"
#include "AliRawDataHeader.h"

#include <cerrno>

// ** This is a global object used for automatic component registration, do not use this
AliHLTRunSummaryProducerComponent gAliHLTRunSummaryProducerComponent;

ClassImp(AliHLTRunSummaryProducerComponent)
    
// ------------------------------------------------------------------------------------------
AliHLTRunSummaryProducerComponent::AliHLTRunSummaryProducerComponent() : 
  fRunSummary(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// ------------------------------------------------------------------------------------------
AliHLTRunSummaryProducerComponent::~AliHLTRunSummaryProducerComponent() {
  // see header file for class documentation
}

// ------------------------------------------------------------------------------------------
const char* AliHLTRunSummaryProducerComponent::GetComponentID() {
  // see header file for class documentation
  return "TriggerRunSummaryProducer"; 
}

// ------------------------------------------------------------------------------------------
void AliHLTRunSummaryProducerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list) {
  // see header file for class documentationAliHLTRunSummary

  list.clear(); 
  list.push_back( kAliHLTDataTypeEventSummary );
  list.push_back( kAliHLTDataTypeRunStatistics );
}

// ------------------------------------------------------------------------------------------
AliHLTComponentDataType AliHLTRunSummaryProducerComponent::GetOutputDataType() {
  // see header file for class documentation

  return kAliHLTDataTypeRunSummary;
}

// ------------------------------------------------------------------------------------------
void AliHLTRunSummaryProducerComponent::GetOutputDataSize( unsigned long& constBase, 
							   double& inputMultiplier ) {
  // see header file for class documentation

  constBase = sizeof( AliHLTRunSummary );
  inputMultiplier = 0.0;
}

// ------------------------------------------------------------------------------------------
AliHLTComponent* AliHLTRunSummaryProducerComponent::Spawn() {
  // Spawn function, return new instance of this class
  // see header file for class documentation

  return new AliHLTRunSummaryProducerComponent;
}


// ------------------------------------------------------------------------------------------
Int_t AliHLTRunSummaryProducerComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation

  Int_t iResult = 0;

  if ( fRunSummary )
    return EINPROGRESS;
  
  fRunSummary = new AliHLTRunSummary();
  
  return iResult;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTRunSummaryProducerComponent::DoDeinit() {
  // see header file for class documentation

  if ( fRunSummary )
    delete fRunSummary;
  fRunSummary = NULL;

  return 0; 
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTRunSummaryProducerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData ) {
  // see header file for class documentation

  // ** Process EventSummary Block
  AliHLTEventSummary* eventSummary = (AliHLTEventSummary*) GetFirstInputObject ( kAliHLTDataTypeEventSummary, "AliHLTEventSummary" );
  
  if ( eventSummary ) 
    ProcessEventSummary ( eventSummary );   

  // ** Increase Number of Events
  fRunSummary->AddNEvents();

  // ** Process CTP trigger information
  ProcessTriggerData( trigData );

  // ** Set RunType / Run Number
  fRunSummary->SetRunNumber ( GetRunNo() );
  fRunSummary->SetRunType   ( GetRunType() );

  // ** PushBack Run Summary
  PushBack ( (TObject*) fRunSummary, sizeof(AliHLTRunSummary), kAliHLTDataTypeRunSummary, (AliHLTUInt32_t) 0 );

  return 0;
}

// -- **********************************************************************************************
// -- *******************************   Processing Functions  **************************************
// -- **********************************************************************************************

// -- **********************************************************************************************
void AliHLTRunSummaryProducerComponent::ProcessEventSummary( AliHLTEventSummary* eventSummary ) {
  // see header file for class documentation
  
  // ** Check if event was accepted or rejected 
  if ( ! eventSummary->IsAccepted() )
    fRunSummary->AddNEventsRejected();
}

// -- **********************************************************************************************
void AliHLTRunSummaryProducerComponent::ProcessTriggerData( AliHLTComponentTriggerData& trigData ) {
  // see header file for class documentation
  
  const AliRawDataHeader* cdh = NULL;
  if (AliHLTComponent::ExtractTriggerData(trigData, NULL, NULL, &cdh, NULL, true) != 0) return;
  AliHLTUInt64_t triggerClasses = cdh->GetTriggerClasses();

  for ( Int_t ndx = 0; ndx < gkNCTPTriggerClasses; ndx ++ ) {
    
    if ( triggerClasses & 0x1 )
      fRunSummary->AddTriggerClass( ndx );
    
    triggerClasses = triggerClasses >> 1;
  }

}

