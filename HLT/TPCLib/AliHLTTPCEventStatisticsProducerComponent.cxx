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

/** @file   AliHLTTPCEventStatisticsProducerComponent.cxx
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

#include "AliHLTTPCEventStatisticsProducerComponent.h"

#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCTrack.h"

#include "TMath.h"

#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCEventStatisticsProducerComponent)
    
// ------------------------------------------------------------------------------------------
AliHLTTPCEventStatisticsProducerComponent::AliHLTTPCEventStatisticsProducerComponent() : 
  fClusterThreshold(0),
  fEvStat(NULL),
  fNTracksPerSlice(),
  fTracks(NULL),
  fGlobalTracks(kFALSE),
  fNSlice(0) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// ------------------------------------------------------------------------------------------
AliHLTTPCEventStatisticsProducerComponent::~AliHLTTPCEventStatisticsProducerComponent() {
  // see header file for class documentation
}

// ------------------------------------------------------------------------------------------
const char* AliHLTTPCEventStatisticsProducerComponent::GetComponentID() {
  // see header file for class documentation
  return "TPCEventStatisticsProducer"; 
}

// ------------------------------------------------------------------------------------------
void AliHLTTPCEventStatisticsProducerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list) {
  // see header file for class documentation

  list.clear(); 
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
  list.push_back( AliHLTTPCDefinitions::fgkTracksDataType );
  list.push_back( AliHLTTPCDefinitions::fgkTrackSegmentsDataType );
}

// ------------------------------------------------------------------------------------------
AliHLTComponentDataType AliHLTTPCEventStatisticsProducerComponent::GetOutputDataType() {
  // see header file for class documentation

  return kAliHLTDataTypeEventStatistics|kAliHLTDataOriginTPC;
}

// ------------------------------------------------------------------------------------------
void AliHLTTPCEventStatisticsProducerComponent::GetOutputDataSize( unsigned long& constBase, 
								   double& inputMultiplier ) {
  // see header file for class documentation

  constBase = sizeof( AliHLTTPCEventStatistics );
  inputMultiplier = 0.0;
}

// ------------------------------------------------------------------------------------------
// Spawn function, return new instance of this class
AliHLTComponent* AliHLTTPCEventStatisticsProducerComponent::Spawn() {
  // see header file for class documentation

  return new AliHLTTPCEventStatisticsProducerComponent;
}


// ------------------------------------------------------------------------------------------
Int_t AliHLTTPCEventStatisticsProducerComponent::DoInit( int argc, const char** argv ) {
  // see header file for class documentation

  Int_t iResult = 0;
  
  // ** Defaults

  Int_t threshold = 0;
  
  // ** Argument handling

  TString argument = "";
  TString parameter = "";
  Int_t bMissingParam=0;
  
  for ( Int_t ii=0; ii<argc && iResult>=0; ii++ ) {
  
    argument = argv[ii];

    if ( argument.IsNull() ) continue;

    // -detector
    if ( ! argument.CompareTo("-clusterThreshold") ) {

      if ( ( bMissingParam=( ++ii >= argc ) ) ) break;

      parameter = argv[ii];  
      parameter.Remove( TString::kLeading, ' ' );
      
      if  (parameter.IsDigit() ) {
	threshold = (Int_t) parameter.Atoi();
	HLTInfo( "Threshold of clusters for long tracks is set to %d.", threshold );
      } 
      else {
	HLTError( "Cannot convert clusterThreshold  specifier '%s'.", parameter.Data() );
	iResult = -EINVAL;
	break;
      }
      
    } // if ( ! argument.CompareTo("-clusterThreshold") ) {
    
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

  // set members
  
  fClusterThreshold =  threshold; 

  return iResult;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTTPCEventStatisticsProducerComponent::DoDeinit() {
  // see header file for class documentation

  if ( fEvStat )
    delete fEvStat;
  fEvStat = NULL;

  if ( fTracks ) 
    delete fTracks; 
  fTracks = NULL;

  return 0;
}

// ------------------------------------------------------------------------------------------
Int_t AliHLTTPCEventStatisticsProducerComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation

  const AliHLTComponentBlockData* iter = NULL;

  // ** No readout list for SOR and EOR event
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;
  
  //
  // ** Setup for new Event
  //
  InitializeEvent();
  
  //
  // ** Read in cluster from ClusterFinder
  //
  
  // ** Loop over all input blocks and specify which data format should be read
  for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock() ) {
    
    AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
    AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );

    HLTDebug ( "Input Data - TPC cluster - Slice/Patch: %d/%d.", slice, patch );

    AddClusters( iter->fPtr, (Int_t) slice, (Int_t) patch );
  } // for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock() ) {
  

  //
  // ** Read in track segments from Tracker
  //

  iter = NULL;

  // ** Loop over all input blocks and specify which data format should be read
  for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTrackSegmentsDataType); iter != NULL; iter = GetNextInputBlock() ) {
    
    AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
    HLTDebug ( "Input Data - TPC track segments - Slice: %d.", slice );

    AddTracks( iter->fPtr, (Int_t) slice );
  } //   for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTrackSegmentsDataType); iter != NULL; iter = GetNextInputBlock() ) {
  
  //
  // ** Read in tracks from GlobalMerger
  //

  iter = NULL;

  // ** Loop over all input blocks and specify which data format should be read
  for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTracksDataType); iter != NULL; iter = GetNextInputBlock() ) {

    HLTDebug ( "Input Data - TPC track segments." );

    AddTracks( iter->fPtr );
  } //   for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTracksDataType); iter != NULL; iter = GetNextInputBlock() ) {
  

  TObject* iterObject = NULL;

  // ** GetESDs
  iterObject = (TObject* )GetFirstInputObject(kAliHLTDataTypeESDTree|kAliHLTDataOriginTPC);

  if ( iterObject ) {

    HLTDebug ( "Input Data - TPC ESD." );

    AddESD( (TTree*) iterObject );
  } 
  
  //
  // ** Produce event statistics
  //
  
  ProcessEvent();
  
  PushBack ( (TObject*) GetEventStatistics(), sizeof(AliHLTTPCEventStatistics), kAliHLTDataTypeEventStatistics|kAliHLTDataOriginTPC, (AliHLTUInt32_t) 0 );

  return 0;
}

// -- **********************************************************************************************
// -- *******************************   Processing Functions  **************************************
// -- **********************************************************************************************

// -- **********************************************************************************************
void AliHLTTPCEventStatisticsProducerComponent::InitializeEvent() {
  // see header file for class documentation
  
  if ( fEvStat )
    delete fEvStat;

  fEvStat = new AliHLTTPCEventStatistics();
  fEvStat->SetClusterThreshold( fClusterThreshold );

  // ** Initialize arrays to 0 
  memset( fNTracksPerSlice, -1, 36*sizeof(Int_t) ); 

  // ** New AliHLTTPCTrackArray
  if ( fTracks ) 
    delete fTracks; 
  fTracks = new AliHLTTPCTrackArray;
  
  fGlobalTracks = kFALSE;
  fNSlice = 0;
}

// -- **********************************************************************************************
void AliHLTTPCEventStatisticsProducerComponent::AddClusters( void* ptr, Int_t slice, Int_t patch ) {
  // see header file for class documentation

  const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*) ptr;

  // these 2 variables have been introduced to avoid warning: unused variable in production compile
  Int_t sliceAntiWarning = slice;
  Int_t patchAntiWarning = patch;

  Int_t nSpacepoint = (Int_t) clusterData->fSpacePointCnt;
  HLTDebug( "%d Clusters found for slice %u - patch %u\n", nSpacepoint, sliceAntiWarning, patchAntiWarning );

  sliceAntiWarning = 0;
  patchAntiWarning = 0;

  // ** Add to event statistics
  fEvStat->AddNTotalCluster( nSpacepoint );

}

// -- **********************************************************************************************
void AliHLTTPCEventStatisticsProducerComponent::AddTracks( void* ptr, Int_t slice ) {
  // see header file for class documentation

  const AliHLTTPCTrackletData* trackData = (const AliHLTTPCTrackletData*) ptr;

  Int_t nTracks = (Int_t) trackData->fTrackletCnt;

  AliHLTTPCTrackSegmentData* segmentData = (AliHLTTPCTrackSegmentData*) &( trackData->fTracklets[0] );

  // ** Tracks which are not from GlobalMerger have to be transformed to global Coordinates ( bTransform=1, else 0 )
  Int_t bTransform;

  
  // ** Global tracks
  if ( slice == -1 ) {
    HLTDebug( "%d tracks found for event.", nTracks );
  
    bTransform = 0;
    fGlobalTracks = kTRUE;
  } 
  else { 
    HLTDebug( "%d tracks found for slice %u.", nTracks, slice );

    // ** Fill number if tracks per slice
    fNTracksPerSlice[slice] = (Int_t) trackData->fTrackletCnt;
    fNSlice++;

    bTransform = 1;
  }
  
  fTracks->FillTracksChecked( segmentData, nTracks, 0, slice, bTransform );

}

// -- **********************************************************************************************
void AliHLTTPCEventStatisticsProducerComponent::AddESD( TTree* esdTree ) {
  // see header file for class documentation

  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);
  esdTree->GetEntry(0);

  HLTWarning("Number of tracks found : %d", esd->GetNumberOfTracks() );
  
  if ( esd )
    delete esd;

}

// -- **********************************************************************************************
void AliHLTTPCEventStatisticsProducerComponent::ProcessEvent() { 
  // see header file for class documentation

  // ** Total number of tracks -- Add to event statistics
  fEvStat->SetNTotalTracks( fTracks->GetNTracks() ) ;

  // ** Loop over all tracks
  for( Int_t trackNdx=0; trackNdx < fTracks->GetNTracks(); trackNdx++ ) {

    AliHLTTPCTrack *track = fTracks->GetCheckedTrack( trackNdx ); 
    if( !track ) continue;

    // ** Used cluster -- Add to event statistics
    fEvStat->AddNUsedCluster( track->GetNHits() );

    // ** If is long track -- Add to event statistics
    if ( track->GetNHits() >= fClusterThreshold )
      fEvStat->AddNTracksAboveClusterThreshold();
  }

  // ** Average of ( number of clusters per track ), floored -- Add to event statistics
  //    N used cluster / N tracks
  if (  fEvStat->GetNTotalTracks() > 0 ) {
    fEvStat->SetAvgClusterPerTrack( (Int_t) TMath::FloorNint( (Double_t) fEvStat->GetNUsedCluster() / (Double_t) fEvStat->GetNTotalTracks()  ) );
  }
  fEvStat->SetAvgClusterPerTrack( 0 );

  if ( fGlobalTracks ) FillTracksPerSlice();

  // ** Average  of ( tracks per slice )
  if ( fNSlice > 0 )
    fEvStat->SetNAvgTracksPerSector( (Int_t) TMath::FloorNint( (Double_t) fEvStat->GetNTotalTracks() / (Double_t) fNSlice ) );
  else
    fEvStat->SetNAvgTracksPerSector( 0 );

  // ** Max and Min of Tracks per slice 
  for ( Int_t ii=0; ii < 36; ii++ ) {

    if ( fNTracksPerSlice[ii] == -1 ) continue;

    if ( fNTracksPerSlice[ii] > fEvStat->GetNMaxTracksPerSector() ) 
      fEvStat->SetNMaxTracksPerSector( fNTracksPerSlice[ii] );

    if ( fNTracksPerSlice[ii] < fEvStat->GetNMinTracksPerSector() ) 
      fEvStat->SetNMinTracksPerSector( fNTracksPerSlice[ii] );
  }

  // ** Info prints

  HLTInfo( "Total N Tracks : %d", fEvStat->GetNTotalTracks() );
  HLTInfo( "Total N Tracks with more than %d cluster : %d",  fEvStat->GetClusterThreshold(), fEvStat->GetNTracksAboveClusterThreshold() );
  HLTInfo( "Min N Tracks per slice : %d", fEvStat->GetNMinTracksPerSector() );
  HLTInfo( "Max N Tracks per slice : %d", fEvStat->GetNMaxTracksPerSector() );
  HLTInfo( "Avg N Tracks per slice : %d", fEvStat->GetNAvgTracksPerSector() );

  HLTInfo( "Total N Cluster : %d", fEvStat->GetNTotalCluster() );
  HLTInfo( "Used N Cluster : %d", fEvStat->GetNUsedCluster() );
  HLTInfo( "Average (Cluster per track) : %d", fEvStat->GetAvgClusterPerTrack() );
}

// -- **********************************************************************************************
void AliHLTTPCEventStatisticsProducerComponent::FillTracksPerSlice() { 
  // see header file for class documentation

  // fill fNTracksPerSlice[36];   
  // fill fNSlice
}
