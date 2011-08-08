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
/// The alignment itself is performed by AliMillepede.
/// A OCDB entry is written with the alignment parameters.
/// A list of graph are written to monitor the alignment parameters.
///
/// \author Javier Castillo, CEA/Saclay - Irfu/SPhN
//-----------------------------------------------------------------------------

// this class header must always come first
#include "AliMUONAlignmentTask.h"

// local header must come before system headers
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliMagF.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliGeomManager.h"

#include "AliMpCDB.h"
#include "AliMUONCDB.h"
#include "AliMUONAlignment.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONESDInterface.h"

// system headers
#include <fstream>
#include <TString.h>
#include <TError.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <Riostream.h>

///\cond CLASSIMP
ClassImp(AliMUONAlignmentTask)
///\endcond

//________________________________________________________________________
AliMUONAlignmentTask::AliMUONAlignmentTask(const char *name, const char *newalignocdb, const char *oldalignocdb, const char *defaultocdb, const char *geofilename):
  AliAnalysisTaskSE(name),
  fReadRecords( kFALSE ),
  fWriteRecords( kFALSE ),
  fDoAlignment( kTRUE ),
  fAlign(0x0),
  fGeoFilename(geofilename),
  fDefaultStorage(defaultocdb),
  fOldAlignStorage(oldalignocdb),
  fNewAlignStorage(newalignocdb),
  fOldGeoTransformer(NULL),
  fNewGeoTransformer(NULL),
  fLoadOCDBOnce(kFALSE),
  fOCDBLoaded(kFALSE),
  fTrackTot(0),
  fTrackOk(0),
  fLastRunNumber(-1),
  fMSDEx(0x0),
  fMSDEy(0x0),
  fMSDEz(0x0),
  fMSDEp(0x0),
  fList(0x0),
  fRecords(0x0),
  fRecordCount(0)
{
  /// Default Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

  // Output slot #0 writes NTuple/histos into a TList
  DefineOutput(1, TList::Class());

  // initialize parameters ...
  for(Int_t k=0;k<4*156;k++)
  {
    fParameters[k]=0.;
    fErrors[k]=0.;
    fPulls[k]=0.;
  }

  fOldGeoTransformer = new AliMUONGeometryTransformer();

}

//________________________________________________________________________
AliMUONAlignmentTask::AliMUONAlignmentTask(const AliMUONAlignmentTask& other):
  AliAnalysisTaskSE(other),
  fReadRecords( other.fReadRecords ),
  fWriteRecords( other.fWriteRecords ),
  fDoAlignment( other.fDoAlignment ),
  fAlign( other.fAlign ),
  fGeoFilename( other.fGeoFilename ),
  fDefaultStorage( other.fDefaultStorage ),
  fOldAlignStorage( other.fOldAlignStorage ),
  fNewAlignStorage( other.fNewAlignStorage ),
  fOldGeoTransformer( other.fOldGeoTransformer ),
  fNewGeoTransformer( other.fNewGeoTransformer ),
  fLoadOCDBOnce( other.fLoadOCDBOnce ),
  fOCDBLoaded( other.fOCDBLoaded ),
  fTrackTot( other.fTrackTot ),
  fTrackOk( other.fTrackOk ),
  fLastRunNumber( other.fLastRunNumber ),
  fMSDEx( other.fMSDEx ),
  fMSDEy( other.fMSDEy ),
  fMSDEz( other.fMSDEz ),
  fMSDEp( other.fMSDEp ),
  fList( other.fList ),
  fRecords( other.fRecords ),
  fRecordCount( other.fRecordCount )
{

  // initialize parameters
  for(Int_t k=0;k<4*156;k++)
  {

    fParameters[k]=other.fParameters[k];
    fErrors[k]=other.fErrors[k];
    fPulls[k]=other.fPulls[k];

  }

}

//________________________________________________________________________
AliMUONAlignmentTask& AliMUONAlignmentTask::operator=(const AliMUONAlignmentTask& other)
{
  /// Assignment
  AliAnalysisTaskSE::operator=(other);

  fReadRecords = other.fReadRecords;
  fWriteRecords = other.fWriteRecords;
  fDoAlignment = other.fDoAlignment;

  // this breaks in destructor
  fAlign = other.fAlign;

  fGeoFilename = other.fGeoFilename;
  fDefaultStorage = other.fDefaultStorage;
  fOldAlignStorage = other.fOldAlignStorage;
  fNewAlignStorage = other.fNewAlignStorage;

  // this breaks in destructor
  fOldGeoTransformer = other.fOldGeoTransformer;
  fNewGeoTransformer = other.fNewGeoTransformer;

  fLoadOCDBOnce = other.fLoadOCDBOnce;
  fOCDBLoaded = other.fOCDBLoaded;
  fTrackTot = other.fTrackTot;
  fTrackOk = other.fTrackOk;
  fLastRunNumber = other.fLastRunNumber;
  fMSDEx = other.fMSDEx;
  fMSDEy = other.fMSDEy;
  fMSDEz = other.fMSDEz;
  fMSDEp = other.fMSDEp;
  fList = other.fList;
  fRecords = other.fRecords;
  fRecordCount = other.fRecordCount;

  // initialize parameters ...
  for( Int_t k=0; k<4*156; ++k)
  {

    fParameters[k]=other.fParameters[k];
    fErrors[k]=other.fErrors[k];
    fPulls[k]=other.fPulls[k];

  }

  return *this;
}

//________________________________________________________________________
AliMUONAlignmentTask::~AliMUONAlignmentTask()
{
  /*
  it is safe to delete NULL pointers, so
  there is no need to test their validity here.
  However, it crashes here, because of incorrect assignment operator
  and copy constructor, resulting in double deletion of objects.
  Would require deep copy instead.
  */
  delete fAlign;
  delete fOldGeoTransformer;
  delete fNewGeoTransformer;
}

//________________________________________________________________________
void AliMUONAlignmentTask::LocalInit()
{
  /// Local initialization, called once per task on the client machine
  /// where the analysis train is assembled

  // print configuration
  AliInfo( Form( "fReadRecords: %s", (fReadRecords ? "kTRUE":"kFALSE" ) ) );
  AliInfo( Form( "fWriteRecords: %s", (fWriteRecords ? "kTRUE":"kFALSE" ) ) );
  AliInfo( Form( "fDoAlignment: %s", (fDoAlignment ? "kTRUE":"kFALSE" ) ) );

  // consistency checks between flags
  /* need at least one of the flags set to true */
  if( !( fReadRecords || fWriteRecords || fDoAlignment ) )
  { AliFatal( "Need at least one of the three flags fReadRecords, fWriteRecords and fDoAlignment set to kTRUE" ); }

  /* cannot read and write records at the same time */
  if( fReadRecords && fWriteRecords )
  { AliFatal( "Cannot set both fReadRecords and fWriteRecords to kTRUE at the same time" ); }

  /* must run alignment when reading records */
  if( fReadRecords && !fDoAlignment )
  {

    AliInfo( "Automatically setting fDoAlignment to kTRUE because fReadRecords is kTRUE" );
    SetDoAlignment( kTRUE );

  }

  // Set initial values here, good guess may help convergence
  // St 1
  //  Int_t iPar = 0;
  //  fParameters[iPar++] =  0.010300 ;  fParameters[iPar++] =  0.010600 ;  fParameters[iPar++] =  0.000396 ;

  fAlign = new AliMUONAlignment();
  fAlign->InitGlobalParameters(fParameters);

//  AliCDBManager::Instance()->Print();
//
//  fAlign->SetGeometryTransformer(fOldGeoTransformer);

  // Do alignment with magnetic field off
  fAlign->SetBFieldOn(kFALSE);

  // Set tracking station to use
  //  Bool_t bStOnOff[5] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  Bool_t bChOnOff[10] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};

  // Set degrees of freedom to align (see AliMUONAlignment)
  fAlign->AllowVariations(bChOnOff);

  // Set expected resolution (see AliMUONAlignment)
  fAlign->SetSigmaXY(0.15,0.10);

  // Fix parameters or add constraints here
  //   for (Int_t iSt=0; iSt<5; iSt++)
  //     if (!bStOnOff[iSt]) fAlign->FixStation(iSt+1);
  for (Int_t iCh=0; iCh<10; iCh++)
    if (!bChOnOff[iCh]) fAlign->FixChamber(iCh+1);

  // Left and right sides of the detector are independent, one can choose to align
  // only one side
  Bool_t bSpecLROnOff[2] = {kTRUE,kTRUE};
  fAlign->FixHalfSpectrometer(bChOnOff,bSpecLROnOff);

  fAlign->SetChOnOff(bChOnOff);
  fAlign->SetSpecLROnOff(bChOnOff);

    // Here we can fix some detection elements
  //  fAlign->FixDetElem(908);
	//  fAlign->FixDetElem(1012);
	fAlign->FixDetElem(608);

  // Set predifined global constrains: X, Y, P, XvsZ, YvsZ, PvsZ, XvsY, YvsY, PvsY
//   Bool_t bVarXYT[9] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
//   Bool_t bDetTLBR[4] = {kFALSE,kTRUE,kFALSE,kTRUE};
  //  fAlign->AddConstraints(bChOnOff,bVarXYT,bDetTLBR,bSpecLROnOff);

}

//________________________________________________________________________
void AliMUONAlignmentTask::UserCreateOutputObjects()
{

  if( fDoAlignment )
  {

    // Creating graphs
    fMSDEx = new TGraphErrors(156);
    fMSDEy = new TGraphErrors(156);
    fMSDEz = new TGraphErrors(156);
    fMSDEp = new TGraphErrors(156);

    // Add Ntuples to the list
    fList = new TList();
    fList->Add(fMSDEx);
    fList->Add(fMSDEy);
    fList->Add(fMSDEz);
    fList->Add(fMSDEp);

    fList->SetOwner(kTRUE);
    PostData(1, fList);
  }

  // connect AOD output
  if( fWriteRecords )
  {
    AliAODHandler* handler = dynamic_cast<AliAODHandler*>( AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler() );
    if( handler )
    {

      // get AOD output handler and add Branch
      fRecords = new TClonesArray( "AliMUONAlignmentTrackRecord", 0 );
      fRecords->SetName( "records" );
      handler->AddBranch("TClonesArray", &fRecords);

      fRecordCount = 0;

    } else AliInfo( "Error: invalid output event handler" );

  }

}

//________________________________________________________________________
void AliMUONAlignmentTask::UserExec(Option_t *)
{
  // print configuration
//  AliInfo( Form( "fReadRecords: %s", (fReadRecords ? "kTRUE":"kFALSE" ) ) );
//  AliInfo( Form( "fWriteRecords: %s", (fWriteRecords ? "kTRUE":"kFALSE" ) ) );
//  AliInfo( Form( "fDoAlignment: %s", (fDoAlignment ? "kTRUE":"kFALSE" ) ) );
  // clear array
  if( fWriteRecords && fRecords )
  {
    fRecords->Clear();
    fRecordCount = 0;
  }

  // local track parameters
  Double_t trackParams[8] = {0.,0.,0.,0.,0.,0.,0.,0.};

  if( fReadRecords ) {

    AliAODEvent* lAOD( dynamic_cast<AliAODEvent*>(InputEvent() ) );

    // check AOD
    if( !lAOD )
    {
      AliInfo("Error: AOD event not available");
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

      if( AliMUONAlignmentTrackRecord* record = dynamic_cast<AliMUONAlignmentTrackRecord*>( records->UncheckedAt(index) ) )
      {

        fAlign->ProcessTrack( record, fDoAlignment );
        if( fDoAlignment )
        { fAlign->LocalFit( fTrackOk++, trackParams, 0 ); }

      } else AliInfo( Form( "Invalid record at %i", index ) );

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

// HUGO: Comment out check on run number, to be able to run on MC
//     if (!lESD->GetRunNumber())
//     {
//       AliInfo( Form( "Current Run Number: %i", lESD->GetRunNumber() ) );
//       return;
//     }

    //  if (lESD->GetRunNumber()!=fLastRunNumber){
    //    fLastRunNumber = lESD->GetRunNumber();
    //    Prepare(fGeoFilename.Data(),fDefaultOCDB.Data(),fMisAlignOCDB.Data());
    //  }

    Int_t nTracks = Int_t(lESD->GetNumberOfMuonTracks());
//    cout << " there are " << nTracks << " tracks" << endl;
    for( Int_t iTrack = 0; iTrack < nTracks; iTrack++ )
    {

      AliESDMuonTrack* esdTrack = lESD->GetMuonTrack(iTrack);
      if (!esdTrack->ClustersStored()) continue;
      if (!esdTrack->ContainTriggerData()) continue;

      Double_t invBenMom = esdTrack->GetInverseBendingMomentum();
      //     fInvBenMom->Fill(invBenMom);
      //     fBenMom->Fill(1./invBenMom);
      if (TMath::Abs(invBenMom)<=1.04)
      {

        AliMUONTrack track;
        AliMUONESDInterface::ESDToMUON(*esdTrack, track);

        // process track and retrieve corresponding records, for storage
        const AliMUONAlignmentTrackRecord* lRecords( fAlign->ProcessTrack( &track, fDoAlignment ) );

        // do the fit, if required
        if( fDoAlignment ) fAlign->LocalFit(fTrackOk++,trackParams,0);
        else fTrackOk++;

        // store in array
        if( fWriteRecords && fRecords )
        {					
          new((*fRecords)[fRecordCount]) AliMUONAlignmentTrackRecord( *lRecords );
          ++fRecordCount;
        }

      }			
			
      ++fTrackTot;

      if(!(fTrackTot%1000))
      { AliInfo( Form( "Processed %i Tracks and %i were fitted.", fTrackTot, fTrackOk ) ); }

//      cout << "Processed " << fTrackTot << " Tracks." << endl;
      // Post final data. Write histo list to a file with option "RECREATE"
      PostData(1,fList);

    }

    // save AOD
    if( fWriteRecords && fRecordCount > 0 ) { 
			AliAODHandler* handler = dynamic_cast<AliAODHandler*>( AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler() );
//			printf("handler: %p\n",handler);
			AliAODEvent* aod = handler->GetAOD();
//			printf("aod: %p\n",aod);
			AliAODHeader* header = aod->GetHeader();
//			printf("header: %p\n",header);
			header->SetRunNumber(lESD->GetRunNumber());
//			printf("RunNumber: %d\n",lESD->GetRunNumber());
			AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()->SetFillAOD(kTRUE); }

  }

}

//________________________________________________________________________
void AliMUONAlignmentTask::Terminate(const Option_t*)
{ return; }

//________________________________________________________________________
void AliMUONAlignmentTask::FinishTaskOutput()
{

  /// Called once per task on the client machine at the end of the analysis.
  AliInfo( Form( "Processed %i tracks.", fTrackTot ) );
  AliInfo( Form( "Accepted %i tracks.", fTrackOk ) );

  // stop here if no alignment is to be performed
  if( !fDoAlignment ) return;

  AliLog::SetGlobalLogLevel(AliLog::kInfo);

  // Perform global fit
  fAlign->GlobalFit(fParameters,fErrors,fPulls);

//   // Update pointers reading them from the output slot
//   fList = (TList*)GetOutputData(0);
//   fMSDEx = (TGraphErrors*)fList->At(0);
//   fMSDEy = (TGraphErrors*)fList->At(1);
//   fMSDEz = (TGraphErrors*)fList->At(2);
//   fMSDEp = (TGraphErrors*)fList->At(3);

  // Store results
  Double_t DEid[156] = {0};
  Double_t MSDEx[156] = {0};
  Double_t MSDEy[156] = {0};
  Double_t MSDEz[156] = {0};
  Double_t MSDEp[156] = {0};
  Double_t DEidErr[156] = {0};
  Double_t MSDExErr[156] = {0};
  Double_t MSDEyErr[156] = {0};
  Double_t MSDEzErr[156] = {0};
  Double_t MSDEpErr[156] = {0};
  Int_t lNDetElem = 4*2+4*2+18*2+26*2+26*2;
  Int_t lNDetElemCh[10] = {4,4,4,4,18,18,26,26,26,26};
  // Int_t lSNDetElemCh[10] = {4,8,12,16,34,52,78,104,130,156};
  Int_t idOffset = 0; // 400
  Int_t lSDetElemCh = 0;
  for(Int_t iDE=0; iDE<lNDetElem; iDE++)
  {

    DEidErr[iDE] = 0.;
    DEid[iDE] = idOffset+100;
    DEid[iDE] += iDE;
    lSDetElemCh = 0;
    for(Int_t iCh=0; iCh<9; iCh++)
    {
      lSDetElemCh += lNDetElemCh[iCh];

      if (iDE>=lSDetElemCh)
      {
        DEid[iDE] += 100;
        DEid[iDE] -= lNDetElemCh[iCh];
      }

    }

    MSDEx[iDE]=fParameters[4*iDE+0];
    MSDEy[iDE]=fParameters[4*iDE+1];
    MSDEz[iDE]=fParameters[4*iDE+3];
    MSDEp[iDE]=fParameters[4*iDE+2];
    MSDExErr[iDE]=(Double_t)fAlign->GetParError(4*iDE+0);
    MSDEyErr[iDE]=(Double_t)fAlign->GetParError(4*iDE+1);
    MSDEzErr[iDE]=(Double_t)fAlign->GetParError(4*iDE+3);
    MSDEpErr[iDE]=(Double_t)fAlign->GetParError(4*iDE+2);
    fMSDEx->SetPoint(iDE,DEid[iDE],fParameters[4*iDE+0]);
    fMSDEx->SetPointError(iDE,DEidErr[iDE],(Double_t)fAlign->GetParError(4*iDE+0));
    fMSDEy->SetPoint(iDE,DEid[iDE],fParameters[4*iDE+1]);
    fMSDEy->SetPointError(iDE,DEidErr[iDE],(Double_t)fAlign->GetParError(4*iDE+1));
    fMSDEz->SetPoint(iDE,DEid[iDE],fParameters[4*iDE+3]);
    fMSDEz->SetPointError(iDE,DEidErr[iDE],(Double_t)fAlign->GetParError(4*iDE+3));
    fMSDEp->SetPoint(iDE,DEid[iDE],fParameters[4*iDE+2]);
    fMSDEp->SetPointError(iDE,DEidErr[iDE],(Double_t)fAlign->GetParError(4*iDE+2));

  }

  // Post final data. Write histo list to a file with option "RECREATE"
  PostData(1,fList);

  // HUGO: stop here to test reproducibility
	//  return;

  // Re Align
  fNewGeoTransformer = fAlign->ReAlign(fOldGeoTransformer,fParameters,true);
  //  newTransform->WriteTransformations("transform2ReAlign.dat");

  // Generate realigned data in local cdb
  const TClonesArray* array = fNewGeoTransformer->GetMisAlignmentData();

  // 100 mum residual resolution for chamber misalignments?
  fAlign->SetAlignmentResolution(array,-1,0.01,0.01,0.004,0.003);

  // CDB manager
  AliLog::SetGlobalDebugLevel(2);
  AliCDBManager* cdbm = AliCDBManager::Instance();
  // cdbManager->SetDefaultStorage(fDefaultOCDB.Data());

	// recover default storage full name (raw:// cannot be used to set specific storage)
	TString defaultStorage(cdbm->GetDefaultStorage()->GetType());
	if (defaultStorage == "alien") defaultStorage += Form("://folder=%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
	else defaultStorage += Form("://%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());

	if (fOldAlignStorage != "none") cdbm->UnloadFromCache("MUON/Align/Data");
	if (!fNewAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fNewAlignStorage.Data());
	//	else cdbm->SetSpecificStorage("MUON/Align/Data",fDefaultStorage.Data());
	else cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());


  AliCDBMetaData* cdbData = new AliCDBMetaData();
  cdbData->SetResponsible("Dimuon Offline project");
  cdbData->SetComment("MUON alignment objects with residual misalignment");
  AliCDBId id("MUON/Align/Data", 0, AliCDBRunRange::Infinity());
  cdbm->Put(const_cast<TClonesArray*>(array), id, cdbData);

	gSystem->Exec("cp -a ReAlignOCDB/MUON/Align/Data/Run0_999999999_v0_s0.root Run0_999999999_v0_s0.root");
	gSystem->Exec("ls -l");

}

//________________________________________________________________________
void AliMUONAlignmentTask::NotifyRun()
{
  
  if (fOCDBLoaded && fLoadOCDBOnce) { AliError(Form("OCDB already loaded %d %d",fOCDBLoaded,fLoadOCDBOnce));
    return;
  }

  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage(fDefaultStorage.Data());
//	cdbm->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  // HUGO: add to comment out the specific settings below in order to be able to run.
  //cdbm->SetSpecificStorage("GRP/GRP/Data",fDefaultStorage.Data());
	//cdbm->SetSpecificStorage("GRP/Geometry/data",fDefaultStorage.Data());
	//cdbm->SetSpecificStorage("MUON/Calib/MappingData",fDefaultStorage.Data());
	cdbm->SetRun(fCurrentRunNumber);

  if (!AliMUONCDB::LoadField()) { AliError("Problem field"); return;}

  // set the magnetic field for track extrapolations
  AliMUONTrackExtrap::SetField();

  if (!AliMUONCDB::LoadMapping(kTRUE)) { AliError("Problem mapping"); return;}

	// recover default storage full name (raw:// cannot be used to set specific storage)
	TString defaultStorage(cdbm->GetDefaultStorage()->GetType());
	if (defaultStorage == "alien") defaultStorage += Form("://folder=%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
	else defaultStorage += Form("://%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());

	// reset existing geometry/alignment if any
	if (cdbm->GetEntryCache()->Contains("GRP/Geometry/Data")) cdbm->UnloadFromCache("GRP/Geometry/Data");
	if (cdbm->GetEntryCache()->Contains("MUON/Align/Data")) cdbm->UnloadFromCache("MUON/Align/Data");
	if (AliGeomManager::GetGeometry()) AliGeomManager::GetGeometry()->UnlockGeometry();

	// get original geometry transformer
	AliGeomManager::LoadGeometry();
	if (!AliGeomManager::GetGeometry()) { AliError("Problem geometry"); return;}
	if (fOldAlignStorage != "none")
  {

		if (!fOldAlignStorage.IsNull()) cdbm->SetSpecificStorage("MUON/Align/Data",fOldAlignStorage.Data());
		else cdbm->SetSpecificStorage("MUON/Align/Data",defaultStorage.Data());
		//		else cdbm->SetSpecificStorage("MUON/Align/Data",fDefaultStorage.Data());
		
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");

  }

  // fOldGeoTransformer = new AliMUONGeometryTransformer();
	fOldGeoTransformer->LoadGeometryData();
	fAlign->SetGeometryTransformer(fOldGeoTransformer);

	fOCDBLoaded = kTRUE;

}
