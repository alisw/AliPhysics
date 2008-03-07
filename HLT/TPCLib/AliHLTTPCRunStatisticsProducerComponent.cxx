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

/** @file   AliHLTTPCRunStatisticsProducerComponent.cxx
    @author Jochen Thaeder
    @date   
    @brief  Component for the @see AliHLTTPCEventStatisticsProducer class
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt   

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCRunStatisticsProducerComponent.h"
#include "AliHLTTPCEventStatistics.h"
#include "TMath.h"

#include <cerrno>

// ** This is a global object used for automatic component registration, do not use this
AliHLTTPCRunStatisticsProducerComponent gAliHLTTPCRunStatisticsProducerComponent;

ClassImp(AliHLTTPCRunStatisticsProducerComponent)
    
// ------------------------------------------------------------------------------------------
AliHLTTPCRunStatisticsProducerComponent::AliHLTTPCRunStatisticsProducerComponent() : 
  fRunStat(NULL),
  fIsHeader(kFALSE) {

  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// ------------------------------------------------------------------------------------------
AliHLTTPCRunStatisticsProducerComponent::~AliHLTTPCRunStatisticsProducerComponent() {
  // see header file for class documentation
}

// ------------------------------------------------------------------------------------------
const char* AliHLTTPCRunStatisticsProducerComponent::GetComponentID() {
  // see header file for class documentation
  return "TPCRunStatisticsProducer"; 
}

// ------------------------------------------------------------------------------------------
void AliHLTTPCRunStatisticsProducerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list) {
  // see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeEventStatistics );
}

// ------------------------------------------------------------------------------------------
AliHLTComponentDataType AliHLTTPCRunStatisticsProducerComponent::GetOutputDataType() {
  // see header file for class documentation

  return kAliHLTDataTypeRunStatistics|kAliHLTDataOriginTPC;
}

// ------------------------------------------------------------------------------------------
void AliHLTTPCRunStatisticsProducerComponent::GetOutputDataSize( unsigned long& constBase, 
								 double& inputMultiplier ) {
  // see header file for class documentation

  constBase = sizeof( AliHLTTPCRunStatistics );
  inputMultiplier = 0.0;
}

// ------------------------------------------------------------------------------------------
AliHLTComponent* AliHLTTPCRunStatisticsProducerComponent::Spawn() {
  // Spawn function, return new instance of this class
  // see header file for class documentation

  return new AliHLTTPCRunStatisticsProducerComponent;
}
 
// ------------------------------------------------------------------------------------------
Int_t AliHLTTPCRunStatisticsProducerComponent::DoInit( int /*argc*/, const char** /*argv*/ ) {
  // see header file for class documentation

  Int_t iResult = 0;
  
  if ( fRunStat )
    return EINPROGRESS;

  fRunStat = new AliHLTTPCRunStatistics();

  return iResult;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTTPCRunStatisticsProducerComponent::DoDeinit() {
  // see header file for class documentation

  if ( fRunStat )
    delete fRunStat;
  fRunStat = NULL;

  return 0;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTTPCRunStatisticsProducerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
							AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation

  // ** Process EventStatistics Block
  AliHLTTPCEventStatistics* evStat = (AliHLTTPCEventStatistics*) GetFirstInputObject ( kAliHLTDataTypeEventStatistics|kAliHLTDataOriginTPC, 
										       "AliHLTTPCEventStatistics" );
  if ( evStat )
    ProcessEventStatistics ( evStat );
  
  if ( ! fIsHeader ) {
    fRunStat->SetDetectorName("TPC");
    
    if ( evStat )
      fRunStat->SetClusterThreshold( evStat->GetClusterThreshold() );
  }

  // ** increase number of events
  fRunStat->AddNEvents();


  PushBack ( (TObject*) GetRunStatistics(), sizeof(AliHLTTPCRunStatistics), kAliHLTDataTypeRunStatistics|kAliHLTDataOriginTPC, (AliHLTUInt32_t) 0 );

  return 0;
}

// -- **********************************************************************************************
// -- *******************************   Processing Functions  **************************************
// -- **********************************************************************************************

// -- **********************************************************************************************
void AliHLTTPCRunStatisticsProducerComponent::ProcessEventStatistics( AliHLTTPCEventStatistics* evStat ) {
  // see header file for class documentation
  
  // ** add number of total tracks
  fRunStat->AddNTotalTracks( evStat->GetNTotalTracks() );

  // ** add number of long tracks  ( tracks above cluster threshold )
  fRunStat->AddNTracksAboveClusterThreshold( evStat->GetNTracksAboveClusterThreshold() );

  // ** add number of total clusters
  fRunStat->AddNTotalCluster( evStat->GetNTotalCluster() );

  // ** add number of cluster used in tracks
  fRunStat->AddNUsedCluster( evStat->GetNUsedCluster() );

  // ** add number of average clusters per track
  fRunStat->AddAvgClusterPerTrack( evStat->GetAvgClusterPerTrack() );

}



