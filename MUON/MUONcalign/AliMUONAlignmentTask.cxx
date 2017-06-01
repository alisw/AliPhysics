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
/// \class AliMUONAlignmentTask
/// AliAnalysisTask to align the MUON spectrometer.
/// The Task reads as input ESDs and feeds the MUONTracks to AliMUONAlignment.
/// The alignment itself is performed by AliMillePede2.
/// A OCDB entry is written with the alignment parameters.
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN
/// \author Hugo Pereira Da Costa, CEA/Saclay - Irfu/SPhN
//-----------------------------------------------------------------------------

// this class header must always come first
#include "AliMUONAlignmentTask.h"

// local header must come before system headers
#include "AliAnalysisManager.h"
#include "AliAlignObjMatrix.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliCDBEntry.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliGeomManager.h"

#include "AliMagF.h"
#include "AliMillePedeRecord.h"
#include "AliMpCDB.h"
#include "AliMUONAlignment.h"
#include "AliMUONCDB.h"
#include "AliMUONESDInterface.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONConstants.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVCluster.h"

// system headers
#include <cassert>
#include <fstream>

// root headers
#include <TString.h>
#include <TError.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <Riostream.h>

///\cond CLASSIMP
ClassImp(AliMUONAlignmentTask)
ClassImp(AliMUONAlignmentTrackParam)
///\endcond

//________________________________________________________________________
Int_t AliMUONAlignmentTrackParam::GetNChambers( void ) const
{
  Int_t out(0);
  for( Int_t iCh = 0; iCh < 10; ++iCh )
  { if( fHitPattern&(1<<iCh) ) out++; }

  return out;

}

//______________________________________
Int_t AliMUONAlignmentTrackParam::GetNStations( void ) const
{

  Int_t out(0);
  for( Int_t iSt = 0; iSt < 5; ++iSt )
  { if( fHitPattern&( (1<<(iSt*2)) | (1<<(iSt*2+1) ) ) ) out++; }

  return out;

}


//________________________________________________________________________
// Double_t Square( Double_t x ) { return x*x; }
// defined in AliMUONAlignent0
Double_t Square( Double_t x );

//________________________________________________________________________
AliMUONAlignmentTask::AliMUONAlignmentTask( const char *name ):
  AliAnalysisTaskSE(name),
  fReadRecords( kFALSE ),
  fWriteRecords( kFALSE ),
  fDoAlignment( kTRUE ),
  fMergeAlignmentCDBs( kTRUE ),
  fForceBField( kFALSE ),
  fBFieldOn( kFALSE ),
  fRefitStraightTracks( kFALSE ),
  fDoEvaluation( kFALSE ),
  fMinPt( 0 ),
  fMinPyPz( 0 ),
  fMinStations( 0 ),
  fAlign(0x0),
  fDefaultStorage(),
  fOldAlignStorage(),
  fNewAlignStorage( "local://ReAlignOCDB" ),
  fOldGeoTransformer(0x0),
  fNewGeoTransformer(0x0),
  fLoadOCDBOnce(kFALSE),
  fOCDBLoaded(kFALSE),
  fEvent(0),
  fTrackTot(0),
  fTrackOk(0),
  fRecordsTot(0),
  fRunNumberMin( 0 ),
  fRunNumberMax( AliCDBRunRange::Infinity() ),
  fTrackParams(0x0),
  fRecords(0x0),
  fRecordCount(0)
{
  /// Default Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

  // initialize parameters ...
  for(Int_t k=0;k<AliMUONAlignment::fNGlobal;k++)
  {
    fParameters[k]=0.;
    fErrors[k]=0.;
    fPulls[k]=0.;
  }

  // create alignment object
  fAlign = new AliMUONAlignment();

  // create old geometry transformer
  fOldGeoTransformer = new AliMUONGeometryTransformer();

}

//________________________________________________________________________
AliMUONAlignmentTask::~AliMUONAlignmentTask()
{
  /// destructor
  delete fAlign;
  delete fOldGeoTransformer;
  delete fNewGeoTransformer;
}

//________________________________________________________________________
void AliMUONAlignmentTask::LocalInit()
{
  /**
  Local initialization, called once per task on the client machine
  where the analysis train is assembled
  **/

  /* must run alignment when reading records */
  if( fReadRecords && !fDoAlignment )
  {

    std::cout << "I-AliMUONAlignmentTask::LocalInit: " << "Automatically setting fDoAlignment to kTRUE because fReadRecords is kTRUE" << std::endl;
    SetDoAlignment( kTRUE );

  }

  // print configuration
  if( fRunNumberMin > 0 || fRunNumberMax > 0 )
  { std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "run range: %i - %i", fRunNumberMin, fRunNumberMax ) << std::endl; }

  std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fReadRecords: %s", (fReadRecords ? "kTRUE":"kFALSE" ) ) << std::endl;
  std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fWriteRecords: %s", (fWriteRecords ? "kTRUE":"kFALSE" ) ) << std::endl;
  std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fDoAlignment: %s", (fDoAlignment ? "kTRUE":"kFALSE" ) ) << std::endl;

  if( fDoAlignment )
  {
    // merge alignment DB flag is irrelevant if not actually performing the alignemnt
    std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fMergeAlignmentCDBs: %s", (fMergeAlignmentCDBs ? "kTRUE":"kFALSE" ) ) << std::endl;
  }

  // storage elements
  if( !fDefaultStorage.IsNull() ) std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fDefaultStorage: %s", fDefaultStorage.Data() ) << std::endl;
  if( !fOldAlignStorage.IsNull() ) std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fOldAlignStorage: %s", fOldAlignStorage.Data() ) << std::endl;

  if( fDoAlignment )
  {
    // new alignment storage is irrelevant if not actually performing the alignemnt
    if( !fNewAlignStorage.IsNull() ) std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fNewAlignStorage: %s", fNewAlignStorage.Data() ) << std::endl;
    else AliFatal( "Invalid new alignment storage path" );
  }

  if( !fReadRecords )
  {
    // following flags are only relevant if not reading records
    if( fForceBField ) std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fBFieldOn: %s", (fBFieldOn ? "kTRUE":"kFALSE" ) ) << std::endl;
    else std::cout << "I-AliMUONAlignmentTask::LocalInit: " << "fBFieldOn: from GRP" << std::endl;

    // refit straight tracks
    std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fRefitStraightTracks: %s", (fRefitStraightTracks ? "kTRUE":"kFALSE" ) ) << std::endl;
    std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fDoEvaluation: %s", (fDoEvaluation ? "kTRUE":"kFALSE" ) ) << std::endl;

    // min px/pz cut
    if( fMinPt > 0 ) std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fMinPt: %f", fMinPt ) << std::endl;

    // min px/pz cut
    if( fMinPyPz > 0 ) std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fMinPyPz: %f", fMinPyPz ) << std::endl;

    // min number of stations
    if( fMinStations > 0 ) std::cout << "I-AliMUONAlignmentTask::LocalInit: " << Form( "fMinStations: %i", fMinStations ) << std::endl;

  }

  // consistency checks between flags
  /* need at least one of the flags set to true */
  if( !( fReadRecords || fWriteRecords || fDoAlignment ) )
  { AliFatal( "Need at least one of the three flags fReadRecords, fWriteRecords and fDoAlignment set to kTRUE" ); }

  /* cannot read and write records at the same time */
  if( fReadRecords && fWriteRecords )
  { AliFatal( "Cannot set both fReadRecords and fWriteRecords to kTRUE at the same time" ); }

  /*
  fix detectors
  warning, counting chambers from 1.
  this must be done before calling the Init method
  */
  if( fDoAlignment )
  {

//     // Note: grouping must be done before fixing chambers
//     AliInfo( "Grouping all detectors by half chambers" );
//     for( Int_t iCh = 5; iCh <= 10; ++iCh )
//     {
//       fAlign->GroupHalfChamber( iCh, 0 );
//       fAlign->GroupHalfChamber( iCh, 1 );
//     }

    // fix chambers
    const Int_t chambers[] = { 1, 10, 0 };
    for( Int_t i = 0; chambers[i] > 0; ++i )
    {
      AliInfo( Form( "Fixing chamber %i", chambers[i] ) );
      fAlign->FixChamber( chambers[i] );
    }

  } else {

    AliInfo( "Not fixing detector elements, since alignment is not required" );

  }

  // propagate magnetic field flag
  fAlign->SetBFieldOn( fBFieldOn );
  fAlign->SetRefitStraightTracks( fRefitStraightTracks );
  fAlign->SetDoEvaluation( fDoEvaluation );

  // initialize
  fAlign->Init();

  // Set expected resolution (see AliMUONAlignment)
  fAlign->SetSigmaXY(0.15,0.10);

  // initialize global parameters to provided values
  fAlign->InitGlobalParameters(fParameters);

}

//________________________________________________________________________
void AliMUONAlignmentTask::UserCreateOutputObjects()
{

  // connect AOD output
  if( fWriteRecords )
  {
    // get AOD output handler and add Branch
    AliAODHandler* handler = dynamic_cast<AliAODHandler*>( AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler() );
    if( handler )
    {

      // track params
      fTrackParams = new TClonesArray( "AliMUONAlignmentTrackParam", 0 );
      fTrackParams->SetName( "trackParams" );
      handler->AddBranch("TClonesArray", &fTrackParams);

      // records
      fRecords = new TClonesArray( "AliMillePedeRecord", 0 );
      fRecords->SetName( "records" );
      handler->AddBranch("TClonesArray", &fRecords);

      fRecordCount = 0;

    } else AliFatal( "Error: invalid output event handler" );

  }

}

//________________________________________________________________________
void AliMUONAlignmentTask::UserExec(Option_t *)
{

  // do nothing if run number not in range
  if( fRunNumberMin > 0 && fCurrentRunNumber < fRunNumberMin ) return;
  if( fRunNumberMax > 0 && fCurrentRunNumber > fRunNumberMax ) return;

  // increase event count
  ++fEvent;

  // clear array
  if( fWriteRecords )
  {
    if( fTrackParams ) fTrackParams->Clear();
    if( fRecords ) fRecords->Clear();
    fRecordCount = 0;
  }

  if( fReadRecords )
  {

    AliAODEvent* lAOD( dynamic_cast<AliAODEvent*>(InputEvent() ) );

    // check AOD
    if( !lAOD )
    {
      AliInfo("Error: AOD event not available");
      return;
    }

    // read track parameters
    TClonesArray* trackParams = static_cast<TClonesArray*>( lAOD->FindListObject( "trackParams" ) );
    if( !trackParams )
    {
      AliInfo( "Unable to read object name \"trackParams\"" );
      return;
    }

    // read records
    TClonesArray* records = static_cast<TClonesArray*>( lAOD->FindListObject( "records" ) );
    if( !records )
    {
      AliInfo( "Unable to read object name \"records\"" );
      return;
    }

    // get number of records and print
    const Int_t lRecordsCount( records->GetEntriesFast() );
    fTrackTot += lRecordsCount;

    // loop over records
    for( Int_t index = 0; index < lRecordsCount; ++index )
    {

      // read track param
      AliMUONAlignmentTrackParam* trackParam = 0;
      if( !( trackParam = dynamic_cast<AliMUONAlignmentTrackParam*>( trackParams->UncheckedAt(index) ) ) )
      {
        AliInfo( Form( "Invalid trackParam at %i", index ) );
        continue;
      }

      // apply minPyPz cut
      const Bool_t pypzOk( fMinPyPz <= 0 || TMath::Sqrt( Square(trackParam->fPy) + Square(trackParam->fPz) ) >= fMinPyPz );
      if( !pypzOk ) continue;

      // apply minPt cut
      const Bool_t ptOk( fMinPt <= 0 || TMath::Sqrt( Square(trackParam->fPx) + Square(trackParam->fPy) ) >= fMinPt );
      if( !ptOk ) continue;

      // check number of stations
      if( fMinStations > 0 && trackParam->GetNStations() < fMinStations )
      {
        AliInfo( "Skipping track because no enough stations hit" );
        continue;
      }

      // read record
      AliMillePedeRecord* record = 0;
      if( !( record = dynamic_cast<AliMillePedeRecord*>( records->UncheckedAt(index) ) ) )
      {
        AliInfo( Form( "Invalid record at %i", index ) );
        continue;
      }

      // process track
      fAlign->ProcessTrack( record );
      ++fTrackOk;

      if(!(fTrackTot%100))
      { AliInfo( Form( "Processed %i Tracks and %i were fitted.", fTrackTot, fTrackOk ) ); }

    }

  } else {

    /// Main loop, called for each event
    AliESDEvent* lESD( dynamic_cast<AliESDEvent*>(InputEvent()) );

    // check ESD
    if( !lESD )
    {
      AliInfo("Error: ESD event not available");
      return;
    }

    Int_t nTracks = Int_t(lESD->GetNumberOfMuonTracks());
    for( Int_t iTrack = 0; iTrack < nTracks; iTrack++ )
    {

      AliESDMuonTrack* esdTrack = lESD->GetMuonTrack(iTrack);
      if (!esdTrack->ContainTrackerData()) continue;
      if (!esdTrack->ContainTriggerData()) continue;

      // apply minPyPz and pt cut
      const Bool_t pypzOk( fMinPyPz <= 0 || TMath::Sqrt( Square(esdTrack->Py()) + Square(esdTrack->Pz()) ) >= fMinPyPz );
      const Bool_t ptOk( fMinPt <= 0 || esdTrack->Pt() >= fMinPt );
      const Bool_t invBendMomentumOk( TMath::Abs( esdTrack->GetInverseBendingMomentum() ) <= 1.04 );
      if( pypzOk && ptOk && invBendMomentumOk )
      {

        AliMUONTrack track;
        AliMUONESDInterface::ESDToMUON(*esdTrack, track);

        // process track and retrieve corresponding records, for storage
        const AliMillePedeRecord* lRecords( fAlign->ProcessTrack( &track, fDoAlignment ) );
        ++fTrackOk;

        // store in array
        if( fWriteRecords )
        {

          if( fTrackParams )
          {

            // create local track parameters
            AliMUONAlignmentTrackParam trackParams( esdTrack->Px(), esdTrack->Py(), esdTrack->Pz() );

            // fill hit pattern
            AliMUONTrack track;
            AliMUONESDInterface::ESDToMUON(*esdTrack, track);
            TObjArray* trackParamAtCluster( track.GetTrackParamAtCluster() );
            for( Int_t iCluster = 0; iCluster < trackParamAtCluster->GetEntries(); ++iCluster )
            {

              AliMUONTrackParam* trackParam( (AliMUONTrackParam*) trackParamAtCluster->At(iCluster) );
              if (!trackParam) continue;

              // cluster
              AliMUONVCluster* cluster( trackParam->GetClusterPtr() );
              if( !cluster ) continue;

              trackParams.fHitPattern |= (1<<cluster->GetChamberId());

            }

            // append to TClonesArray
            new((*fTrackParams)[fRecordCount]) AliMUONAlignmentTrackParam( trackParams );

          }
          if( fRecords ) new((*fRecords)[fRecordCount]) AliMillePedeRecord( *lRecords );
          ++fRecordCount;
        }

      }

      ++fTrackTot;

      if(!(fTrackTot%100))
      { AliInfo( Form( "Processed %i Tracks and %i were fitted.", fTrackTot, fTrackOk ) ); }

      // Post final data. Write histo list to a file with option "RECREATE"
      // PostData(1,fList);

    }

    // save AOD
    fRecordsTot += fRecords->GetEntriesFast();
    if( fWriteRecords && fRecordCount > 0 )
    {
      AliAODHandler* handler = dynamic_cast<AliAODHandler*>( AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler() );
      if( handler )
      {
        AliAODEvent* aod = handler->GetAOD();
        AliAODHeader* header = dynamic_cast<AliAODHeader*>(aod->GetHeader());
        header->SetRunNumber(lESD->GetRunNumber());
        AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE);
        AliInfo( Form( "Fitted Tracks: %i, records: %i", fTrackOk, fRecordsTot ) );

      } else AliInfo( "Error: invalid output event handler" );

    }

  }

}

//________________________________________________________________________
void AliMUONAlignmentTask::FinishTaskOutput()
{

  /// Called once per task on the client machine at the end of the analysis.
  AliInfo( Form( "Processed %i tracks.", fTrackTot ) );
  AliInfo( Form( "Accepted %i tracks.", fTrackOk ) );

  // terminate output
  fAlign->Terminate();

  // stop here if no alignment is to be performed
  if( !fDoAlignment ) return;

  AliLog::SetGlobalLogLevel(AliLog::kInfo);

  // Perform global fit
  fAlign->GlobalFit(fParameters,fErrors,fPulls);

  // store misalignments from OCDB into old transformers
  if( fMergeAlignmentCDBs )
  { SaveMisAlignmentData( fOldGeoTransformer ); }

  // Re Align
  fNewGeoTransformer = fAlign->ReAlign( fOldGeoTransformer, fParameters, true );

  // Generate realigned data in local cdb
  const TClonesArray* array = fNewGeoTransformer->GetMisAlignmentData();

  // 100 mum residual resolution for chamber misalignments?
  fAlign->SetAlignmentResolution( array, -1, 0.01, 0.01, 0.004, 0.003 );

  // CDB manager
  AliLog::SetGlobalDebugLevel(2);

  // recover default storage full name (raw:// cannot be used to set specific storage)
  AliCDBManager* cdbManager = AliCDBManager::Instance();

  // unload old alignment path
  if( cdbManager->GetEntryCache()->Contains("MUON/Align/Data") )
  { cdbManager->UnloadFromCache("MUON/Align/Data"); }

  // load new alignment path
  if( !fNewAlignStorage.IsNull() ) cdbManager->SetSpecificStorage("MUON/Align/Data",fNewAlignStorage.Data());
  else {

    TString defaultStorage(cdbManager->GetDefaultStorage()->GetType());
    if (defaultStorage == "alien") defaultStorage += Form("://folder=%s", cdbManager->GetDefaultStorage()->GetBaseFolder().Data());
    else defaultStorage += Form("://%s", cdbManager->GetDefaultStorage()->GetBaseFolder().Data());
    cdbManager->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());

  }

  // create new DB entry
  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON alignment objects with residual misalignment");
  AliCDBId id("MUON/Align/Data", fRunNumberMin, fRunNumberMax );
  cdbManager->Put(const_cast<TClonesArray*>(array), id, cdbData);

}

//________________________________________________________________________
void AliMUONAlignmentTask::NotifyRun()
{

  /// run number (re)initialization

  // propagate run number to fAlign
  if( fAlign ) fAlign->SetRunNumber( fCurrentRunNumber );

  // update ocdb
  if (fOCDBLoaded && fLoadOCDBOnce)
  {
    AliError(Form("OCDB already loaded %d %d",fOCDBLoaded,fLoadOCDBOnce));
    return;
  }

  AliCDBManager* cdbManager = AliCDBManager::Instance();

  // Initialize default storage
  if( !( cdbManager->IsDefaultStorageSet() || fDefaultStorage.IsNull() ) )
  {

    AliInfo( Form( "Initializing default storage: %s", fDefaultStorage.Data() ) );
    cdbManager->SetDefaultStorage(fDefaultStorage.Data());

  } else if( !fDefaultStorage.IsNull() ) {

    AliInfo( "Default storage already set. Ignoring fDefaultStorage" );

  } else {

    AliInfo( "Default storage already set" );

  }

  // Initialize run number
  if( cdbManager->GetRun() < 0 )
  {
    AliInfo( Form( "Setting run number: %i", fCurrentRunNumber ) );
    cdbManager->SetRun(fCurrentRunNumber);
  }

  // following initialization is not needed when reading records
  if( !fReadRecords )
  {

    // load magnetic field if needed
    if( !TGeoGlobalMagField::Instance()->IsLocked() )
    {

      AliInfo( "Loading magnetic field" );
      if( !AliMUONCDB::LoadField() )
      {
        AliError( "Failed to load magnetic field" );
        return;
      }

    } else { AliInfo( "Magnetic field already locked" ); }

    // checking magnitic field
    if( !fForceBField )
    {
      AliInfo( "Reading magnetic field setup from GRP" );

      // decide bFieldOn value base on dipole current
      // propagete to Alignment class
      // and printout
      AliMagF* magF = dynamic_cast<AliMagF*>( TGeoGlobalMagField::Instance()->GetField() );
      fBFieldOn = TMath::Abs( magF->GetFactorDip() ) > 1e-5;
      fAlign->SetBFieldOn( fBFieldOn );
      AliInfo( Form( "Dipole magnetic field factor: %.2f", magF->GetFactorDip() ) );
      AliInfo( Form( "fBFieldOn = %s", (fBFieldOn ? "kTRUE":"kFALSE") ) );
    }

    AliInfo( "Loading muon mapping" );
    if( !AliMUONCDB::LoadMapping(kFALSE) )
    {
      AliError( "Failed to load muon mapping" );
      return;
    }

    AliInfo( "Assigning field to Track extrapolator" );
    AliMUONTrackExtrap::SetField();

  }

  // load geometry if needed
  if( !AliGeomManager::GetGeometry() )
  {

    // reset existing geometry/alignment if any
    if( cdbManager->GetEntryCache()->Contains("GRP/Geometry/Data") )
    {
      AliInfo( "Unloading GRP/Geometry/Data" );
      cdbManager->UnloadFromCache("GRP/Geometry/Data");
    }

    if( cdbManager->GetEntryCache()->Contains("MUON/Align/Data") )
    {
      AliInfo( Form( "Unloading MUON/Align/Data from %s", cdbManager->GetSpecificStorage( "MUON/Align/Data" )->GetBaseFolder().Data() ) );
      cdbManager->UnloadFromCache("MUON/Align/Data");
    }

    // get original geometry transformer
    AliInfo( "Loading geometry" );
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry())
    {
      AliError("Failed to load geometry");
      return;
    }

    if (!fOldAlignStorage.IsNull())
    {

      AliInfo( Form( "Initializing MUON/Align/Data using: %s", fOldAlignStorage.Data() ) );
      cdbManager->SetSpecificStorage("MUON/Align/Data",fOldAlignStorage.Data());

    } else {

      // recover default storage full name (raw:// cannot be used to set specific storage)
      TString defaultStorage(cdbManager->GetDefaultStorage()->GetType());
      if (defaultStorage == "alien") defaultStorage += Form("://folder=%s", cdbManager->GetDefaultStorage()->GetBaseFolder().Data());
      else defaultStorage += Form("://%s", cdbManager->GetDefaultStorage()->GetBaseFolder().Data());

      AliInfo( Form( "Re-initializing MUON/Align/Data using: %s", defaultStorage.Data() ) );
      cdbManager->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());

    }

    AliInfo("Loading muon Alignment objects");
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");

  } else { AliInfo( "Geometry already initialized - fOldAlignStorage ignored" ); }

  // update geometry transformer and pass to Alignment object
  AliInfo("Loading muon geometry in transformer");
  fOldGeoTransformer->LoadGeometryData();
  fAlign->SetGeometryTransformer(fOldGeoTransformer);

  fOCDBLoaded = kTRUE;

}

//_____________________________________________________________________________________
void AliMUONAlignmentTask::SaveMisAlignmentData( AliMUONGeometryTransformer* transformer ) const
{

  // clear transformer
  transformer->ClearMisAlignmentData();

  // load MUON/Align/Data from OCDB
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  AliCDBEntry* cdbEntry = cdbManager->Get( "MUON/Align/Data" );
  if (!cdbEntry)
  {

    AliError( "unable to load entry for path MUON/Align/Data" );
    return;

  }

  // get TClonesArray and check
  TClonesArray* misArray = (TClonesArray*)cdbEntry->GetObject();
  if( !misArray )
  {

    AliError( "unable to load old misalignment array" );
    return;

  }

  // loop over stored entries
  for (Int_t index=0; index<misArray->GetEntriesFast(); ++index )
  {

    // load matrix and check
    AliAlignObjMatrix* matrix = (AliAlignObjMatrix*) misArray->At( index );
    if( !matrix )
    {
      AliError( Form( "unable to load matrix for index %i", index ) );
      continue;
    }

    // get volume ID
    AliGeomManager::ELayerID layerId;
    Int_t moduleId;
    matrix->GetVolUID( layerId, moduleId);

    AliInfo( Form( "layer: %i, module: %i name: %s", layerId, moduleId, matrix->GetSymName() ) );

    // make sure ELayerID is correct
    assert( layerId == AliGeomManager::kMUON );

    // printout
    // AliInfo( Form( "Found matrix for %s %i", matrix->GetSymName(), moduleId ) );

    // get matrix
    TGeoHMatrix misMatrix;
    matrix->GetMatrix(misMatrix);

    // add to geometry transformer
    // need the detector element
    // "detElement->GetId()"
    // see fOldGeoTransformer->GetMisAlignmentData( ... )

    if( TString( matrix->GetSymName() ).Contains( "DE" ) ) transformer->AddMisAlignDetElement( moduleId,  misMatrix);
    else transformer->AddMisAlignModule( moduleId,  misMatrix);
  }

  return;
}

