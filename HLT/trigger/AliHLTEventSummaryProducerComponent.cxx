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

/** @file   AliHLTEventSummaryProducerComponent.cxx
    @author Jochen Thaeder
    @date   
    @brief  Produces a event summary as @see AliHLTEventSummary
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTEventSummaryProducerComponent.h"
#include "AliHLTTPCEventStatistics.h"
#include "AliRawDataHeader.h"

#include <cerrno>

// ** This is a global object used for automatic component registration, do not use this
AliHLTEventSummaryProducerComponent gAliHLTEventSummaryProducerComponent;

ClassImp(AliHLTEventSummaryProducerComponent)
    
// ------------------------------------------------------------------------------------------
AliHLTEventSummaryProducerComponent::AliHLTEventSummaryProducerComponent() : 
  fEventSummary(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// ------------------------------------------------------------------------------------------
AliHLTEventSummaryProducerComponent::~AliHLTEventSummaryProducerComponent() {
  // see header file for class documentation
}

// ------------------------------------------------------------------------------------------
const char* AliHLTEventSummaryProducerComponent::GetComponentID() {
  // see header file for class documentation
  return "TriggerEventSummaryProducer"; 
}

// ------------------------------------------------------------------------------------------
void AliHLTEventSummaryProducerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list) {
  // see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeEventStatistics );
}

// ------------------------------------------------------------------------------------------
AliHLTComponentDataType AliHLTEventSummaryProducerComponent::GetOutputDataType() {
  // see header file for class documentation

  return kAliHLTDataTypeEventSummary;
}

// ------------------------------------------------------------------------------------------
void AliHLTEventSummaryProducerComponent::GetOutputDataSize( unsigned long& constBase, 
							     double& inputMultiplier ) {
  // see header file for class documentation

  constBase = sizeof( AliHLTEventSummary );
  inputMultiplier = 0.0;
}

// ------------------------------------------------------------------------------------------
AliHLTComponent* AliHLTEventSummaryProducerComponent::Spawn() {
  // Spawn function, return new instance of this class
  // see header file for class documentation

  return new AliHLTEventSummaryProducerComponent;
}


// ------------------------------------------------------------------------------------------
Int_t AliHLTEventSummaryProducerComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation

  Int_t iResult = 0;
 
  return iResult;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTEventSummaryProducerComponent::DoDeinit() {
  // see header file for class documentation

  if ( fEventSummary )
    delete fEventSummary;
  fEventSummary = NULL;

  return 0;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTEventSummaryProducerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData ) {
  // see header file for class documentation

  if ( fEventSummary )
    delete fEventSummary;

  fEventSummary = new AliHLTEventSummary;

  const TObject* iter = NULL;

  // ** Loop over all event statistics input objs and add them to the event summary
  for ( iter = GetFirstInputObject(kAliHLTDataTypeEventStatistics|kAliHLTDataOriginAny); iter != NULL; iter = GetNextInputObject() ) {

    fEventSummary->AddDetector( (TObject*) iter );

  } // for ( iter = GetFirstInputObject(kAliHLTDataTypeEventEventStatistics|kAliHLTDataOriginAny); iter != NULL; iter = GetNextInputObject() ) {
  
  // ** Process CTP trigger information
  ProcessTriggerData( trigData );

  // ** Set RunType / Run Number
  fEventSummary->SetRunNumber ( GetRunNo() );
  fEventSummary->SetRunType   ( GetRunType() );

  // ** PushBack Event Summary
  PushBack ( (TObject*) fEventSummary, sizeof(AliHLTEventSummary), kAliHLTDataTypeEventSummary, (AliHLTUInt32_t) 0 );

  return 0;
}

// -- **********************************************************************************************
// -- *******************************   Processing Functions  **************************************
// -- **********************************************************************************************

// -- **********************************************************************************************
void AliHLTEventSummaryProducerComponent::ProcessTriggerData( AliHLTComponentTriggerData& trigData ) {
  // see header file for class documentation
  
  const AliRawDataHeader* cdh = NULL;
  if (AliHLTComponent::ExtractTriggerData(trigData, NULL, NULL, &cdh, NULL, true) != 0) return;
  AliHLTUInt64_t triggerClass = cdh->GetTriggerClasses();
  fEventSummary->SetTriggerClass( triggerClass );
}
