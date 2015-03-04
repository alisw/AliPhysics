// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/// @file   AliHLTTPCCalibProcessor.cxx
/// @author Chiara Zampolli
/// @date   2013-11-26
/// @brief  A test HLT processing component
/// @ingroup alihlt_tutorial

#include "AliHLTTPCCalibProcessor.h"
#include "AliHLTDataTypes.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTReadoutList.h"
#include "AliLog.h"
#include "TMap.h"
#include "TObjString.h"
#include "TH2F.h"
#include "TH1F.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCalibProcessor)

AliHLTTPCCalibProcessor::AliHLTTPCCalibProcessor():
fGlobalDistr(kFALSE),
  fMinTracks(0),
  fMinEvents(0),
  fEvents(0),
  fMeanPt(0),
  fTotTracks(0),
  fHistoMeanPt(0x0),
  fHistoMult(0x0),
  fHistoPt(0x0),
  fHistoClusters(0x0),
  fHistoClustersPerTrack(0x0)

{
  // constructor
  // note: all internal allocation only happens in DoInit

  Printf("----------> TPCCalibProcessor: constructor");
}

AliHLTTPCCalibProcessor::~AliHLTTPCCalibProcessor()
{
  // destructor
}

const char* AliHLTTPCCalibProcessor::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "TPCCalibProcessor";
}

void AliHLTTPCCalibProcessor::GetInputDataTypes( AliHLTComponentDataTypeList& list)
{
  /// inherited from AliHLTComponent: list of data types in the vector reference
  list.clear();
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC );
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
  
}

AliHLTComponentDataType AliHLTTPCCalibProcessor::GetOutputDataType()
{
  /// inherited from AliHLTComponent: output data type of the component.
  return kAliHLTDataTypeHistogram;
}

void AliHLTTPCCalibProcessor::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  /// inherited from AliHLTComponent: output data size estimator
  /// estimated size = inputMultiplier*input + constBase
  constBase=100000;
  inputMultiplier=1.;
}

AliHLTComponent* AliHLTTPCCalibProcessor::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCCalibProcessor;
}

int AliHLTTPCCalibProcessor::ScanArgument(int argc, const char** argv)
{
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
	
  // Three arguments provided: 
  // 1. minimum number of tracks to accept an event
  // 2. minimum number of accepted events to start building average pt histo
  // 3. building distribution of number of tracks/event

  int i=0;
  TString argument=argv[i];

  if (argument.IsNull()) return 0;

  // minTracks
  if (argument.CompareTo("-minTracks")==0) {
    if (++i>=argc) return -EINVAL;
    argument=argv[i];
    fMinTracks = argument.Atoi();
    HLTInfo("got \'-minTracks\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }

  // minEvents
  if (argument.CompareTo("-minEvents")==0) {
    if (++i>=argc) return -EINVAL;
    argument=argv[i];
    fMinEvents = argument.Atoi();
    HLTInfo("got \'-minEvents\' argument: %s", argv[i]);
    return 2; // keyword + 1 argument
  }

  // -argument2: one keyword
  if (argument.CompareTo("-globalDistr")==0) {
    fGlobalDistr = kTRUE;
    HLTInfo("got \'-globalDistr\' argument");
    return 1; // only keyword
  }

  return 0;
}

int AliHLTTPCCalibProcessor::InitCalibration()
{
  // component initialization
  Printf("----------> TPCCalibProcessor: InitCalibration");

  int iResult=0;

  // init stage 1: default values for all data members
  // TODO: insert member initializtion here
  fGlobalDistr = kFALSE;
  fMinTracks = 0;
  fMinEvents = 0;
  fEvents = 0;
  fMeanPt = 0.;
  fTotTracks = 0;

  // implement the component initialization
  fHistoMeanPt = new TH2F("fHistoMeanPt", "fHistoMeanPt; #Sigma#_{mult}; <p#_T>", 1000, -0.5, 499.5, 200., 0., 20.);
  fHistoMult = new TH1F("fHistoMult", "fHistoMult; mult", 100, -0.5, 99.5);
  fHistoPt = new TH1F("fHistoPt", "fHistoPt; p#_{t}", 200, 0, 20.);
  fHistoClusters = new TH1F("fHistoClusters", "fHistoClusters; num of clusters per event", 100000, -0.5, 99999.5);
  fHistoClustersPerTrack = new TH1F("fHistoClustersPerTrack", "fHistoClustersPerTrack; num of clusters per track", 160, -0.5, 159.5);

  return 0;

}

int AliHLTTPCCalibProcessor::DeinitCalibration()
{
  // component cleanup, delete all instances of helper classes here

  if (fHistoMult) {
    delete fHistoMult; 
    fHistoMult = 0x0;
  }
  if (fHistoPt) {
    delete fHistoPt; 
    fHistoPt = 0x0;
  }
  if (fHistoClusters) {
    delete fHistoClusters; 
    fHistoClusters = 0x0;
  }
  if (fHistoClustersPerTrack) {
    delete fHistoClustersPerTrack; 
    fHistoClustersPerTrack = 0x0;
  }
  if (fHistoMeanPt) {
    delete fHistoMeanPt; 
    fHistoMeanPt = 0x0;
  }

  return 0;
}

int AliHLTTPCCalibProcessor::ProcessCalibration(
				     const AliHLTComponentEventData& evtData,   
				     const AliHLTComponent_BlockData* blocks,
				     AliHLTComponentTriggerData& /*trigData*/,
				     AliHLTUInt8_t* outputPtr,
				     AliHLTUInt32_t& size,
				     vector<AliHLTComponentBlockData>& outputBlocks)
{
  Printf("Chiara Log ==================== event %3d ================================", GetEventCount());

  // event processing function
  if (!IsDataEvent()) {
    printf("This is not data\n");
    return 0; // skip if it is not a normal data event
  }
  if ( evtData.fBlockCnt <= 0 ) {
    HLTWarning( "no blocks in event" );
    return 0;
  }

  int iResult=0;
  int nBlocks = (int)evtData.fBlockCnt;
  int nInputClusters = 0;
  int nInputTracks = 0;

  // first read all the clusters
	
  for (int ndx=0; ndx<nBlocks && iResult>=0; ndx++) {
    printf("ndx = %d\n", ndx);
    const AliHLTComponentBlockData* iter = blocks+ndx;
    if ( iter->fDataType != AliHLTTPCDefinitions::fgkClustersDataType ) continue;
		
    // checking the header of the data
    Int_t slice=AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
    Int_t patch=AliHLTTPCDefinitions::GetMinPatchNr(iter->fSpecification);
    Int_t slicepatch=slice*6+patch;
    //if( slicepatch >= fkNPatches ){
    //	HLTWarning("Wrong header of TPC cluster data, slice %d, patch %d",
    //			   slice, patch );
    //	continue;
    //}
    AliHLTTPCClusterData* inPtrSP = ( AliHLTTPCClusterData* )( iter->fPtr );
    nInputClusters += inPtrSP->fSpacePointCnt;
  }

  Printf("Chiara Log: nInputCluster = %d", nInputClusters);
  fHistoClusters->Fill(nInputClusters);

  // loop over the input tracks: calculate dEdx and write output

  unsigned int outSize = 0;
  int tempindex = 0;

  //for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrack);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {
    
    tempindex++;
    Printf("Chiara Log: I am reading the %d-th block", tempindex);

    AliHLTTracksData* dataPtr = ( AliHLTTracksData* ) pBlock->fPtr;  // array of tracks, with info on the number of tracks
    int nTracks = dataPtr->fCount;
    AliHLTExternalTrackParam* currTrack = dataPtr->fTracklets;  // array of tracks
    nInputTracks += nTracks;
		
    fHistoMult->Fill(nTracks);
    if (nTracks >= fMinTracks){
      fEvents++;
    }

    // looping over the tracks

    for( int itr=0; 
	 itr<nTracks && ( (AliHLTUInt8_t *)currTrack < ((AliHLTUInt8_t *) pBlock->fPtr)+pBlock->fSize); 
	 itr++ ){

      // get the clusters of the current track

      int nclusters = currTrack->fNPoints;
      Printf("Chiara Log: nclusters of current track = %d", nclusters);
      fHistoClustersPerTrack->Fill(nclusters);

      // getting all the clusters one by one, just as an exercise

      for( UInt_t ic=0; ic<currTrack->fNPoints; ic++){	    
	UInt_t id = currTrack->fPointIDs[ic];
	int iSlice = AliHLTTPCSpacePointData::GetSlice(id);
	int iPatch = AliHLTTPCSpacePointData::GetPatch(id);
	int iCluster = AliHLTTPCSpacePointData::GetNumber(id);
	if( iSlice<0 || iSlice>36 || iPatch<0 || iPatch>5 ){
	  HLTError("Corrupted TPC cluster Id: slice %d, patch %d, cluster %d",
		   iSlice, iPatch,iCluster );
	  continue;
	}
      }
			
      // Getting 1/pt

      Float_t pt = 1./currTrack->fq1Pt;
      Printf("Chiara Log: track %d has pt = %f", itr, pt);
      fHistoPt->Fill(pt);

      if (fEvents > fMinEvents){						
	fMeanPt += pt;
      }

      unsigned int step = sizeof( AliHLTExternalTrackParam ) + currTrack->fNPoints * sizeof( unsigned int );
      currTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currTrack) + step );  
    } // end loop on tracks

    // now filling the other histogram

    fHistoMeanPt->Fill(fTotTracks, fMeanPt/fTotTracks);
		
  }

  // publish the two histograms as output
  PushBack(fHistoMeanPt, kAliHLTDataTypeHistogram|kAliHLTDataOriginSample);
  PushBack(fHistoPt, kAliHLTDataTypeHistogram|kAliHLTDataOriginSample);
  PushBack(fHistoClusters, kAliHLTDataTypeHistogram|kAliHLTDataOriginSample);
  PushBack(fHistoClustersPerTrack, kAliHLTDataTypeHistogram|kAliHLTDataOriginSample);
  if (fGlobalDistr) PushBack(fHistoMult, kAliHLTDataTypeHistogram|kAliHLTDataOriginSample);

  return 0;
}


Int_t AliHLTTPCCalibProcessor::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, 
						      AliHLTComponentTriggerData& /*trigData*/ ) {
  // see header file for class documentation
    
  // ** PushBack data to FXS ...
  static AliHLTReadoutList rdList(AliHLTReadoutList::kTPC);
  PushToFXS( (TObject*) fHistoMeanPt, "TPC", "MeanPt", &rdList ) ;
  
  return 0;
} 

