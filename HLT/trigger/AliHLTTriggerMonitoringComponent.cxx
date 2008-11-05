//-*- Mode: C++ -*-
// $Id: AliHLTTriggerMonitoringComponent.cxx 24328 2008-03-06 13:26:00Z richterm $
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

/** @file   AliHLTTriggerMonitoringComponent.cxx
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

#include "AliHLTTriggerMonitoringComponent.h"
#include "AliHLTTPCEventStatistics.h"
#include "AliHLTEventStatistics.h"
#include <cerrno>

// ** This is a global object used for automatic component registration, do not use this
AliHLTTriggerMonitoringComponent gAliHLTTriggerMonitoringComponent;

ClassImp(AliHLTTriggerMonitoringComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
   
// ------------------------------------------------------------------------------------------
AliHLTTriggerMonitoringComponent::AliHLTTriggerMonitoringComponent() : 
  fTotalTrackCut(-1),
  fLongTrackCut(-1),
  fEventSummary(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// ------------------------------------------------------------------------------------------
AliHLTTriggerMonitoringComponent::~AliHLTTriggerMonitoringComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// ------------------------------------------------------------------------------------------
const char* AliHLTTriggerMonitoringComponent::GetComponentID() {
  // see header file for class documentation
  return "TriggerMonitoring"; 
}

// ------------------------------------------------------------------------------------------
void AliHLTTriggerMonitoringComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list) {
  // see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeEventSummary );
}

// ------------------------------------------------------------------------------------------
AliHLTComponentDataType AliHLTTriggerMonitoringComponent::GetOutputDataType() {
  // see header file for class documentation

  return kAliHLTDataTypeEventSummary;
}

// ------------------------------------------------------------------------------------------
void AliHLTTriggerMonitoringComponent::GetOutputDataSize( unsigned long& constBase, 
							  double& inputMultiplier ) {
  // see header file for class documentation

  constBase = sizeof( AliHLTEventSummary );
  inputMultiplier = 0.0;
}

// ------------------------------------------------------------------------------------------
AliHLTComponent* AliHLTTriggerMonitoringComponent::Spawn() {
  // Spawn function, return new instance of this class
  // see header file for class documentation

  return new AliHLTTriggerMonitoringComponent;
}


// ------------------------------------------------------------------------------------------
Int_t AliHLTTriggerMonitoringComponent::DoInit( int argc, const char** argv ) {
  // see header file for class documentation

  Int_t iResult = 0;
 
  // ** Argument handling

  TString argument = "";
  TString parameter = "";
  Int_t bMissingParam=0;
  
  for ( Int_t ii=0; ii<argc && iResult>=0; ii++ ) {
  
    argument = argv[ii];

    if ( argument.IsNull() ) continue;

    // -triggerTotalTracks
    if ( ! argument.CompareTo("-triggerTotalTracks") ) {

      if ( ( bMissingParam=( ++ii >= argc ) ) ) break;
      parameter = argv[ii];  
      parameter.Remove( TString::kLeading, ' ' );
      
      if  (parameter.IsDigit() ) {
	fTotalTrackCut = (Int_t) parameter.Atoi();
	HLTInfo( "Trigger on total number of tracks is activated, and threshold is set to  %d.", fTotalTrackCut );
      } 
      else {
	HLTError( "Cannot convert triggerTotalTracks specifier '%s'.", parameter.Data() );
	iResult = -EINVAL;
	break;
      }
      
    } // if ( ! argument.CompareTo("-triggerTotalTracks") ) {

    // -triggerLongTracks
    else if ( ! argument.CompareTo("-triggerLongTracks") ) {

      if ( ( bMissingParam=( ++ii >= argc ) ) ) break;
      parameter = argv[ii];  
      parameter.Remove( TString::kLeading, ' ' );
      
      if  (parameter.IsDigit() ) {
	fLongTrackCut = (Int_t) parameter.Atoi();
	HLTInfo( "Trigger on number of long tracks is activated, and threshold is set to  %d.", fTotalTrackCut );
      } 
      else {
	HLTError( "Cannot convert triggerLongTracks specifier '%s'.", parameter.Data() );
	iResult = -EINVAL;
	break;
      }
    } // if ( ! argument.CompareTo("-triggerLongTracks") ) {
    // - unknow parameter
    else {
      iResult = -EINVAL;
      HLTError("Unknown argument '%s'", argument.Data() );
    }

  } // for ( Int_t ii=0; ii<argc && iResult>=0; ii++ ) {

  if ( bMissingParam ) {
    HLTError( "Missing parameter for argument '%s'.", argument.Data() );
    iResult = -EPROTO;
  }

  return iResult;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTTriggerMonitoringComponent::DoDeinit() {
  // see header file for class documentation

  if ( fEventSummary )
    delete fEventSummary;
  fEventSummary = NULL;

  return 0;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTTriggerMonitoringComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation

  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ))
    return 0;
 
  const TObject* iter = NULL;
  iter = GetFirstInputObject(kAliHLTDataTypeEventSummary|kAliHLTDataOriginAny);

  fEventSummary = (AliHLTEventSummary*) iter;
    
  // Trigger on eventSummary
  Trigger();
  
  if (fEventSummary->IsAccepted()) {
    HLTWarning( "Triggered" );
  } else {
    HLTWarning( "Discarded" );
  }

  HLTWarning( "RunNumber %d ", GetRunNo() );
  
  // ** PushBack Event Summary
  PushBack ( (TObject*) fEventSummary, sizeof(AliHLTEventSummary), kAliHLTDataTypeEventSummary, (AliHLTUInt32_t) 0 );

  return 0;
}

// -- **********************************************************************************************
// -- *******************************   Processing Functions  **************************************
// -- **********************************************************************************************

// -- **********************************************************************************************
void AliHLTTriggerMonitoringComponent::Trigger() {
  // see header file for class documentation
  
  TObjArray* detArray = fEventSummary->GetDetectorArray();

  for ( Int_t iter = 0; iter <= detArray->GetLast(); iter++ ) {
      
    AliHLTEventStatistics* evStat = (AliHLTEventStatistics*) (*detArray)[iter];
    TString detName = evStat->GetDetectorName();
    
    if ( detName.CompareTo("TPC") )
      continue;
    
    HLTWarning("TPC statistics found...");
    
    AliHLTTPCEventStatistics* tpcStat = (AliHLTTPCEventStatistics*) (*detArray)[iter];
    
    if ( fLongTrackCut != -1 && tpcStat->GetNTracksAboveClusterThreshold() < fLongTrackCut ) {
      fEventSummary->RejectEvent();
      return;
    }
    
    if ( fTotalTrackCut != -1 && tpcStat->GetNTotalTracks() < fTotalTrackCut ) {
      fEventSummary->RejectEvent();
      return;
    }

  } // for ( Int_t iter = 0; iter <= detArray->GetLast(); iter++ ) {

}
