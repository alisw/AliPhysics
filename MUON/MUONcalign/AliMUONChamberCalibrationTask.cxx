/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

//-----------------------------------------------------------------------------
/// \file   AliMUONChamberCalibrationTask.cxx
/// \brief  Implementation of the AliMUONChamberCalibrationTask 
/// \author Andry Rakotozafindrabe CEA/IRFU/SPhN
//-----------------------------------------------------------------------------

#include <Riostream.h>

#include <TBranch.h>
#include <TChain.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

#include "AliMUONChamberCalibrationTask.h"

// STEER includes
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliRecoParam.h"
#include "AliTracker.h"

// ANALYSIS includes
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"

// MUON includes
#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONClusterInfo.h"
#include "AliMUONESDInterface.h"
#include "AliMUONPadInfo.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVCluster.h"
//#include "AliMUONVClusterStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"

/// \cond CLASSIMP
ClassImp( AliMUONChamberCalibrationTask )
/// \endcond CLASSIMP

//______________________________________________________________
AliMUONChamberCalibrationTask::AliMUONChamberCalibrationTask():
  AliAnalysisTaskSE( "AliMUONChamberCalibrationTask" ),
  fOCDBPath( "local://$ALICE_ROOT/OCDB" ),
  fCalibChoice(kNOGAIN),
  fClusterInfoTree(0x0),
  fMuonRecoParam(0x0),
  fClusterInfo(0x0),
  fCalibData(0x0),
  fESDInterface(0x0),
  fDigitStore(0x0),
  fESDInputHandler(0x0),
  fESDInputEvent(0x0)
{
  //
  /// Default constructor
  //

}

//______________________________________________________________
AliMUONChamberCalibrationTask::AliMUONChamberCalibrationTask( const char* name,
							      char* ocdbpath,
							      const Int_t my_calib_option ):
  AliAnalysisTaskSE( name ),
  fOCDBPath( "local://$ALICE_ROOT/OCDB" ),
  fCalibChoice(kNOGAIN),
  fClusterInfoTree(0x0),
  fMuonRecoParam(0x0),
  fClusterInfo(0x0),
  fCalibData(0x0),
  fESDInterface(0x0),
  fDigitStore(0x0),
  fESDInputHandler(0x0),
  fESDInputEvent(0x0)
{
  //
  /// constructor
  //

  fOCDBPath = ocdbpath;
  if ( (my_calib_option >= ((Int_t)kNOGAIN)) && (my_calib_option <= ((Int_t)kINJECTIONGAIN)) ) 
    fCalibChoice = (Calibration_t)my_calib_option;
  else {
    AliWarning( Form("Wrong value of the calibration option %d not within [%d, %d] !!! Will use NOGAIN", 
		     my_calib_option, (Int_t)kNOGAIN, (Int_t)kINJECTIONGAIN ) );
    fCalibChoice = kNOGAIN;
  }
}

//______________________________________________________________
AliMUONChamberCalibrationTask::~AliMUONChamberCalibrationTask()
{
  //
  /// destructor
  //

  delete fMuonRecoParam;
  delete fClusterInfo;
  delete fESDInterface;

}
//______________________________________________________________
void AliMUONChamberCalibrationTask::CreateOutputObjects()
{
  //
  /// Creates the output TTree
  //

  AliDebug( 1, "" ); 

  TFile* clusterInfoFile = OpenFile( 0, "RECREATE" );
  if( clusterInfoFile ) clusterInfoFile->SetCompressionLevel(1);
  else AliError( "no output file created !!!" );

  if ( !fClusterInfoTree ) fClusterInfoTree = new TTree( "clusterInfoTree", "clusterInfoTree" ); 
  fClusterInfoTree->Branch( "clusterInfo" , &fClusterInfo, 32000, 99);
}

//______________________________________________________________
void AliMUONChamberCalibrationTask::LocalInit()
{
  //
  /// Initialization
  /// Initialize the cluster info and the ESD interface
  /// Set the magnetic field, the mapping and the reconstruction parameters
  //

  AliDebug( 1, "" );

  // initialize the cluster info and the ESD interface

  fClusterInfo = new AliMUONClusterInfo();
  fESDInterface = new AliMUONESDInterface();

  gRandom->SetSeed(0);
  
  // set mag field

  if ( !TGeoGlobalMagField::Instance()->GetField() ) {
    AliInfo( "Loading field map..." );
    AliGRPManager *grpMan = new AliGRPManager();
    grpMan->ReadGRPEntry();
    grpMan->SetMagField();
    delete grpMan;
  }
  
  // Load mapping

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage( fOCDBPath );
  man->SetSpecificStorage( "MUON/Calib/MappingData", fOCDBPath );
  man->SetSpecificStorage( "MUON/Calib/MappingRunData", fOCDBPath ); // for the manu serial numbers
  man->SetRun(0);
  man->Print();
  if ( ! AliMpCDB::LoadDDLStore() ) {
    AliFatal( "Could not access mapping from OCDB !" );
    exit(-1); 
  }

  // Set the reconstruction parameters for track refitting 
  // (needed when applying any of the with-gain options)

  fMuonRecoParam = AliMUONRecoParam::GetCosmicParam();

  TString caliboption1 = "NOGAIN";
  TString caliboption2 = "GAINCONSTANTCAPA";
  TString caliboption3 = "GAIN";
  TString caliboption4 = "INJECTIONGAIN";

  TString caliboption = caliboption1;
  if ( fCalibChoice == kGAINCONSTANTCAPA ) caliboption = caliboption2;
  if ( fCalibChoice == kGAIN ) caliboption = caliboption3;  
  if ( fCalibChoice == kINJECTIONGAIN ) caliboption = caliboption4;
  fMuonRecoParam->SetCalibrationMode(caliboption.Data());

  for (Int_t iCh=0; iCh<10; iCh++) {
    fMuonRecoParam->SetDefaultNonBendingReso( iCh, 0.152 ); // input ESD was aligned (default cosmic settings)
    fMuonRecoParam->SetDefaultBendingReso( iCh, 0.027 );
  }
  fMuonRecoParam->SetMaxNonBendingDistanceToTrack(5.); // was at 1. in default cosmic settings
  fMuonRecoParam->SetMaxBendingDistanceToTrack(5.);

  fMuonRecoParam->RequestStation(1, kTRUE); // only St 4 and 5 enabled in default cosmic settings
  fMuonRecoParam->ImproveTracks(kTRUE, 7.); // was 6. in default cosmic settings

  AliInfo( "reconstruction parameters initialized as follows :" );
  fMuonRecoParam->Print("FULL");

  AliMUONESDInterface::ResetTracker(fMuonRecoParam);
}

//______________________________________________________________
void AliMUONChamberCalibrationTask::ConnectInputData( Option_t* /*option*/ )
{
  //
  /// Connect to ESD here
  /// (called once)
  //

  AliDebug( 1, "" );

  fESDInputHandler = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  TTree* tree = NULL;

  if ( fESDInputHandler ) {

    // The properly initialized ESD input handler reads ESD tree 
    // and connect it to ESD event, so we only need to retrieve the later
    fESDInputEvent = (AliESDEvent*)fESDInputHandler->GetEvent();
    if ( !fESDInputEvent ) {

      AliFatal( "Could not get input ESD event !!! ");

    } 
  } else {

      AliError( "Could not get input ESD handler !!!" );
      // If no input event handler we need to get the tree once
      // from input slot 0 for the chain
      tree = dynamic_cast<TTree*> (GetInputData(0));
      if ( tree ) tree->GetReadEntry();
      else AliFatal( "Could not read tree from input slot 0 !!!" );
    }
}

//______________________________________________________________
void AliMUONChamberCalibrationTask::Exec( Option_t* /*option*/ )
{
  //
  /// Process the current event
  /// (called for each event)
  //

  static Bool_t first = kTRUE;

  if ( first ) { 
    AliDebug( 1, "" );
    first = kFALSE;
  }

  if ( !fESDInputEvent ) {
    AliError( "Input ESD event not available !!! " );
    return;
  }

  // load the current ESD event
  fESDInterface->LoadEvent( *fESDInputEvent );

  // get digit store
  fDigitStore = fESDInterface->GetDigits();

  // prepare access to calibration data 
  if ( !fCalibData ) fCalibData = new AliMUONCalibrationData( fESDInputEvent->GetESDRun()->GetRunNumber() );

  // --------------------------------------------------------------------
  // fill cluster info from clusters attached to each track of this event
  // --------------------------------------------------------------------

  Int_t nTracks = (Int_t)fESDInputEvent->GetNumberOfMuonTracks(); 
  if ( nTracks < 1 ) return; 

  TIter nextTrack( fESDInterface->CreateTrackIterator() );
  AliMUONTrack* track;

  while ( (track = static_cast<AliMUONTrack*>(nextTrack())) ) { // loop over tracks

    UInt_t muonClusterMap = BuildClusterMap( *track );
    
    AliMUONTrackParam* trackParam = 
      static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->First());

    while ( trackParam ) { // loop over clusters

      fClusterInfo->Clear("C");
      
      // fill cluster info
      AliMUONVCluster* cluster = trackParam->GetClusterPtr();
      fClusterInfo->SetRunId( fESDInputEvent->GetRunNumber() );
      fClusterInfo->SetEventId( fESDInputEvent->GetEventNumberInFile() );
      fClusterInfo->SetZ( cluster->GetZ() );
      fClusterInfo->SetClusterId( cluster->GetUniqueID() );
      fClusterInfo->SetClusterXY( cluster->GetX(), cluster->GetY() );
      fClusterInfo->SetClusterXYErr( cluster->GetErrX(), cluster->GetErrY() );
      fClusterInfo->SetClusterChi2( cluster->GetChi2() );
      fClusterInfo->SetClusterCharge( cluster->GetCharge() );
      
      // fill track info
      fClusterInfo->SetTrackId( track->GetUniqueID() );
      fClusterInfo->SetTrackXY( trackParam->GetNonBendingCoor(), trackParam->GetBendingCoor() );
      fClusterInfo->SetTrackThetaXY( TMath::ATan( trackParam->GetNonBendingSlope() ), 
				     TMath::ATan( trackParam->GetBendingSlope() ) );
      fClusterInfo->SetTrackP( trackParam->P() );
      const TMatrixD paramCov = trackParam->GetCovariances();
      fClusterInfo->SetTrackXYErr( TMath::Sqrt( paramCov(0,0) ), 
				   TMath::Sqrt( paramCov(2,2) ) );
      fClusterInfo->SetTrackChi2( track->GetNormalizedChi2() );
      fClusterInfo->SetTrackCharge( (Short_t)trackParam->GetCharge() );
      fClusterInfo->SetTrackNHits( track->GetNClusters() );
      fClusterInfo->SetTrackChamberHitMap( muonClusterMap );
      
      // fill pad info if available	  
      for ( Int_t i=0; i<cluster->GetNDigits(); i++ ) {

	AliMUONVDigit* digit = fDigitStore->FindObject( cluster->GetDigitId(i) );
	if ( !digit ) continue;
	
	// pad location
	const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->
	  GetMpSegmentation( digit->DetElemId(), AliMp::GetCathodType( digit->Cathode() ) );
	AliMpPad pad = seg->PadByIndices( digit->PadX(), digit->PadY() );
	
	// calibration parameters
	AliMUONVCalibParam* ped = fCalibData->Pedestals( digit->DetElemId(), digit->ManuId() );
	AliMUONVCalibParam* gain = fCalibData->Gains( digit->DetElemId(), digit->ManuId() );
	Int_t manuChannel = digit->ManuChannel();
	Int_t planeType = 0;
	if ( digit->ManuId() & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) ) {
	  planeType = 1;
	}
	
	// fill pad info
	AliMUONPadInfo padInfo;
	padInfo.SetPadId( digit->GetUniqueID() );
	padInfo.SetPadPlaneType( planeType );
	padInfo.SetPadXY( pad.GetPositionX(), pad.GetPositionY() );
	padInfo.SetPadDimXY( pad.GetDimensionX(), pad.GetDimensionY() );
	padInfo.SetPadCharge( (Double_t)digit->Charge() );
	padInfo.SetPadADC( digit->ADC() );
	padInfo.SetSaturated( digit->IsSaturated() );
	padInfo.SetCalibrated( digit->IsCalibrated() );
	padInfo.SetPedestal( ped->ValueAsFloatFast(manuChannel,0), // mean
			     ped->ValueAsFloatFast(manuChannel,1) ); // sigma
	padInfo.SetGain( gain->ValueAsFloatFast(manuChannel,0), // a0
			 gain->ValueAsFloatFast(manuChannel,1), // a1
			 (Int_t)gain->ValueAsFloatFast(manuChannel,2), // threshold
			 (Int_t)gain->ValueAsFloatFast(manuChannel,3) ); // fit quality
	
	fClusterInfo->AddPad( padInfo );
      }
            
      // fill cluster info tree
      fClusterInfoTree->Fill();
      
      trackParam = static_cast<AliMUONTrackParam*>(track->GetTrackParamAtCluster()->After(trackParam));
    }
    
  }
    // Added protection in case the derived task is not an AOD producer.
    AliAnalysisDataSlot *out0 = GetOutputSlot(0);
    if (out0 && out0->IsConnected()) PostData( 0, fClusterInfoTree );    

}


//______________________________________________________________
UInt_t AliMUONChamberCalibrationTask::BuildClusterMap( AliMUONTrack &track )
{
  //
  /// Build the map of clusters in tracking chambers
  //

  UInt_t muonClusterMap = 0;
  
  AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First());
  while ( trackParam ) {
    
    muonClusterMap |= BIT(trackParam->GetClusterPtr()->GetChamberId());
    
    trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
  }
  
  return muonClusterMap;
}

//______________________________________________________________
void AliMUONChamberCalibrationTask::Terminate( Option_t* /*option*/ )
{
  //
  /// Called once per task on the client machine at the end of the analysis.
  //

  AliDebug( 1, "" );
}
