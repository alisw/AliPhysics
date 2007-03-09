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
/* $Id$ */

//-----------------------------
// Class AliMUONReconstructor  
//-----------------------------
// Class for the 
// MUON track reconstruction

#include "AliMUONReconstructor.h"

#include "AliMUONConstants.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONClusterFinderAZ.h"
#include "AliMUONClusterReconstructor.h"
#include "AliMUONData.h"
#include "AliMUONDigitCalibrator.h"
#include "AliMUONEventRecoCombi.h" 
#include "AliMUONDigitMaker.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTracker.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrackReconstructor.h"
#include "AliMUONTrackReconstructorK.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONSegFactory.h"
#include "AliMUONSegmentation.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMUONClusterFinderCOG.h"
#include "AliMUONClusterFinderSimpleFit.h"
#include "AliMUONClusterFinderMLEM.h"

#include "AliESD.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"

#include "TTask.h"
#include "TStopwatch.h"

/// \cond CLASSIMP
ClassImp(AliMUONReconstructor)
/// \endcond

//_____________________________________________________________________________
AliMUONReconstructor::AliMUONReconstructor()
  : AliReconstructor(), 
    fDigitMaker(new AliMUONDigitMaker()), 
    fCalibrationData(0x0),
    fCrateManager(new AliMUONTriggerCrateStore()),
    fTriggerCircuit(new TClonesArray("AliMUONTriggerCircuit", 234)),
    fTransformer(new AliMUONGeometryTransformer(kTRUE)),
    fSegmentation(0x0),
    fMUONData(new AliMUONData(0x0,"MUON","MUON"))
{
/// Default constructor

    AliDebug(1,"");
    // Crate manager
    fCrateManager->ReadFromFile();

    // set to digit maker
    fDigitMaker->SetCrateManager(fCrateManager);

    // transformater
    fTransformer->ReadGeometryData("volpath.dat", "geometry.root");
    
    // create segmentation and pass it to EventRecoCombi
    AliMUONSegFactory factory(fTransformer);
    fSegmentation = factory.CreateSegmentation();
    AliMUONEventRecoCombi::Instance(fSegmentation); 

    // trigger circuit
    for (Int_t i = 0; i < AliMUONConstants::NTriggerCircuit(); i++)  {
      AliMUONTriggerCircuit* c = new AliMUONTriggerCircuit();
      c->SetTransformer(fTransformer);
      c->Init(i,*fCrateManager);
      TClonesArray& circuit = *fTriggerCircuit;
      new(circuit[circuit.GetEntriesFast()])AliMUONTriggerCircuit(*c);
      delete c;
    }

  
}

//_____________________________________________________________________________
AliMUONReconstructor::~AliMUONReconstructor()
{
/// Destructor

  AliDebug(1,"");
  delete fCalibrationData;
  delete fDigitMaker;
  delete fCrateManager;
  delete fTriggerCircuit;
  delete fTransformer;
  delete fSegmentation;
  delete fMUONData;
}

//_____________________________________________________________________________
TTask* 
AliMUONReconstructor::GetCalibrationTask() const
{
/// Create the calibration task(s). 
  
  //const AliRun* run = fRunLoader->GetAliRun();
  //Int_t runNumber = run->GetRunNumber();     
  Int_t runNumber = AliCDBManager::Instance()->GetRun();
  AliInfo("Calibration will occur.");
 
  fCalibrationData = new AliMUONCalibrationData(runNumber);
  if ( !fCalibrationData->IsValid() )
    {
      AliError("Could not retrieve calibrations !");
      delete fCalibrationData;
      fCalibrationData = 0x0;
      return 0x0;
    }    
  TTask* calibration = new TTask("MUONCalibrator","MUON Digit calibrator");
  
  TString opt(GetOption());
  opt.ToUpper();
  Bool_t statusMap(kTRUE);
  
  if ( strstr(opt,"NOSTATUSMAP") )
  {
    AliWarning("Disconnecting status map : SHOULD BE USED FOR DEBUG ONLY. NOT FOR PRODUCTION !!!");
    statusMap = kFALSE; 
  }
  calibration->Add(new AliMUONDigitCalibrator(fMUONData,fCalibrationData,statusMap));
  return calibration;
}

//_____________________________________________________________________________
AliMUONClusterReconstructor*
AliMUONReconstructor::CreateClusterReconstructor() const
{
/// Create cluster reconstructor

  AliMUONVClusterFinder* clusterFinder(0x0);
  
  TString opt(GetOption());
  opt.ToUpper();
  
  if ( strstr(opt,"PRECLUSTER") )
  {
    clusterFinder = new AliMUONPreClusterFinder;
  }  
  else if ( strstr(opt,"COG") )
  {
    clusterFinder = new AliMUONClusterFinderCOG;
  }  
  else if ( strstr(opt,"SIMPLEFIT") )
  {
    clusterFinder = new AliMUONClusterFinderSimpleFit;
  }
  else if ( strstr(opt,"MLEM:DRAW") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kTRUE);
  }
  else if ( strstr(opt,"MLEM") )
  {
    clusterFinder = new AliMUONClusterFinderMLEM(kFALSE);
  } 
  
  if ( clusterFinder) 
  {
    AliInfo(Form("Will use %s for clusterizing",clusterFinder->ClassName()));
  }
  
  AliMUONClusterReconstructor* clusterReco = 
    new AliMUONClusterReconstructor(fMUONData,clusterFinder,fTransformer);
  return clusterReco;
}

//_____________________________________________________________________________
void AliMUONReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
/// Reconstruct
/// \todo add more

  AliLoader* loader = runLoader->GetLoader("MUONLoader");
  Int_t nEvents     = runLoader->GetNumberOfEvents();
  Int_t evtNumber   = runLoader->GetEventNumber();

  fMUONData->SetLoader(loader);

// passing loader as argument.
  AliMUONVTrackReconstructor* recoEvent;
  if (strstr(GetOption(),"Original")) recoEvent = new AliMUONTrackReconstructor(fMUONData);
  else if (strstr(GetOption(),"Combi")) recoEvent = new AliMUONTrackReconstructorK(fMUONData,"Combi");
  else recoEvent = new AliMUONTrackReconstructorK(fMUONData,"Kalman");
  
  recoEvent->SetTriggerCircuit(fTriggerCircuit);

  AliMUONClusterReconstructor* recoCluster = CreateClusterReconstructor();
  
  AliMUONClusterFinderVS *recModel = recoCluster->GetRecoModel();

  if (!strstr(GetOption(),"VS")) {
    recModel = (AliMUONClusterFinderVS*) new AliMUONClusterFinderAZ();
    recoCluster->SetRecoModel(recModel);
  }
  recModel->SetGhostChi2Cut(10);
  recModel->SetEventNumber(evtNumber);

  loader->LoadDigits("READ");
  loader->LoadRecPoints("RECREATE");
  loader->LoadTracks("RECREATE");
  
  TTask* calibration = GetCalibrationTask();
  
  Int_t chBeg = (strstr(GetOption(),"Combi") ? 6 : 0);

  //   Loop over events              
  for(Int_t ievent = 0; ievent < nEvents; ievent++) {

    AliDebug(1,Form("Event %d",ievent));
    
    runLoader->GetEvent(ievent);

    //----------------------- digit2cluster & Trigger2Trigger -------------------
    if (!loader->TreeR()) loader->MakeRecPointsContainer();
     
    // tracking branch
    if (!strstr(GetOption(),"Combi")) {
      fMUONData->MakeBranch("RC");
      fMUONData->SetTreeAddress("D,RC");
    } else {
      fMUONData->SetTreeAddress("D");
      fMUONData->SetTreeAddress("RCC");
    }
    // Important for avoiding a memory leak when reading digits ( to be investigated more in detail)
    // In any case the reading of GLT is needed for the Trigger2Tigger method below
    fMUONData->SetTreeAddress("GLT");

    fMUONData->GetDigits();
    
    if ( calibration ) 
    {
      calibration->ExecuteTask();
    }
    
    recoCluster->Digits2Clusters(chBeg); 
    
    if (strstr(GetOption(),"Combi")) {
      // Combined cluster / track finder
      AliMUONEventRecoCombi::Instance()->FillEvent(fMUONData, (AliMUONClusterFinderAZ*)recModel);
      ((AliMUONClusterFinderAZ*) recModel)->SetReco(2); 
    }
    else fMUONData->Fill("RC"); 

    // trigger branch
    fMUONData->MakeBranch("TC");
    fMUONData->SetTreeAddress("TC");
    recoCluster->Trigger2Trigger(); 
    fMUONData->Fill("TC");

    //AZ loader->WriteRecPoints("OVERWRITE");

    //---------------------------- Track & TriggerTrack ---------------------
    if (!loader->TreeT()) loader->MakeTracksContainer();

    // trigger branch
    fMUONData->MakeBranch("RL"); //trigger track
    fMUONData->SetTreeAddress("RL");
    recoEvent->EventReconstructTrigger();
    fMUONData->Fill("RL");

    // tracking branch
    fMUONData->MakeBranch("RT"); //track
    fMUONData->SetTreeAddress("RT");
    recoEvent->EventReconstruct();
    fMUONData->Fill("RT");

    loader->WriteTracks("OVERWRITE"); 
  
    if (strstr(GetOption(),"Combi")) { 
      // Combined cluster / track
      ((AliMUONClusterFinderAZ*) recModel)->SetReco(1);
      fMUONData->MakeBranch("RC");
      fMUONData->SetTreeAddress("RC");
      AliMUONEventRecoCombi::Instance()->FillRecP(fMUONData, (AliMUONTrackReconstructorK*)recoEvent); 
      fMUONData->Fill("RC"); 
    }
    loader->WriteRecPoints("OVERWRITE"); 

    //--------------------------- Resetting branches -----------------------
    fMUONData->ResetDigits();
    fMUONData->ResetRawClusters();
    fMUONData->ResetTrigger();
    fMUONData->ResetRecTracks();  
    fMUONData->ResetRecTriggerTracks();

  }
  loader->UnloadDigits();
  loader->UnloadRecPoints();
  loader->UnloadTracks();

  delete recoCluster;
  delete recoEvent;
  delete calibration;
}

//_____________________________________________________________________________
void AliMUONReconstructor::Reconstruct(AliRunLoader* runLoader, 
                                       AliRawReader* rawReader) const
{
/// Recontruct
/// \todo add more

  //  AliLoader
  AliLoader* loader = runLoader->GetLoader("MUONLoader");
  Int_t evtNumber   = runLoader->GetEventNumber();
 
  fMUONData->SetLoader(loader);

  // passing loader as argument.
  fDigitMaker->SetMUONData(fMUONData);

  // disable trigger rawdata reading
  if (strstr(GetOption(),"TriggerDisable"))
      fDigitMaker->DisableTrigger();

  AliMUONClusterReconstructor* recoCluster = CreateClusterReconstructor();

  AliMUONClusterFinderVS *recModel = recoCluster->GetRecoModel();

  if (!strstr(GetOption(),"VS")) 
  {
    recModel = (AliMUONClusterFinderVS*) new AliMUONClusterFinderAZ();
    recoCluster->SetRecoModel(recModel);
  }
  recModel->SetGhostChi2Cut(10);
  recModel->SetEventNumber(evtNumber);
  
  TTask* calibration = GetCalibrationTask();
  
  loader->LoadRecPoints("RECREATE");
 
  // Digits are not stored on disk and created on flight from rawdata.
  // In order to write digits on disk the following line should be uncommented
  // loader->LoadDigits("RECREATE"); 

  //   Loop over events  
  Int_t iEvent = 0;
           
  TStopwatch totalTimer;
  TStopwatch rawTimer;
  TStopwatch calibTimer;
  TStopwatch clusterTimer;
  
  rawTimer.Start(kTRUE); rawTimer.Stop();
  calibTimer.Start(kTRUE); calibTimer.Stop();
  clusterTimer.Start(kTRUE); clusterTimer.Stop();
  
  totalTimer.Start(kTRUE);
  
  while (rawReader->NextEvent()) 
  {
    AliDebug(1,Form("Event %d",iEvent));
    
    runLoader->GetEvent(iEvent++);
    
    //----------------------- raw2digits & raw2trigger-------------------
    //  if (!loader->TreeD()) 
    //  {
    //    AliDebug(1,Form("Making Digit Container for event %d",iEvent));
    //    loader->MakeDigitsContainer();
    //  }
    //  Digits are not stored on disk and created on flight from rawdata.
    //  In order to write digits on disk the following lines should be uncommented
    //  fMUONData->MakeBranch("D,GLT");
    //  fMUONData->SetTreeAddress("D,GLT");
    fMUONData->SetDataContainer("D, GLT");
    rawTimer.Start(kFALSE);
    fDigitMaker->Raw2Digits(rawReader);
    rawTimer.Stop();
    
    if ( calibration )
    {
      calibTimer.Start(kFALSE);
      calibration->ExecuteTask();
      calibTimer.Stop();
    }
    // Digits are not stored on disk and created on flight from rawdata.
    // In order to write digits on disk the following lines should be uncommented
    // fMUONData->Fill("D,GLT");
    // loader->WriteDigits("OVERWRITE");
    //----------------------- digit2cluster & Trigger2Trigger -------------------
    clusterTimer.Start(kFALSE);

    if (!loader->TreeR()) loader->MakeRecPointsContainer();
    
    // tracking branch
    fMUONData->MakeBranch("RC");
    fMUONData->SetTreeAddress("RC");
    recoCluster->Digits2Clusters(); 
    fMUONData->Fill("RC"); 

    // trigger branch
    fMUONData->MakeBranch("TC");
    fMUONData->SetTreeAddress("TC");
    fMUONData->Fill("TC");
    
    loader->WriteRecPoints("OVERWRITE");

    clusterTimer.Stop();
    
    
    //--------------------------- Resetting branches -----------------------
    fMUONData->ResetDigits();
    fMUONData->ResetRawClusters();
    fMUONData->ResetTrigger();
  }
  
  totalTimer.Stop();
  
  loader->UnloadRecPoints();
  loader->UnloadDigits();
  
  delete recoCluster;
  
  AliInfo(Form("Execution time for converting RAW data to digits in MUON : R:%.2fs C:%.2fs",
               rawTimer.RealTime(),rawTimer.CpuTime()));
  AliInfo(Form("Execution time for calibrating MUON : R:%.2fs C:%.2fs",
               calibTimer.RealTime(),calibTimer.CpuTime()));
  AliInfo(Form("Execution time for clusterizing MUON : R:%.2fs C:%.2fs",
               clusterTimer.RealTime(),clusterTimer.CpuTime()));
  AliInfo(Form("Total Execution time for Reconstruct(from raw) MUON : R:%.2fs C:%.2fs",
               totalTimer.RealTime(),totalTimer.CpuTime()));
}

//_____________________________________________________________________________
void AliMUONReconstructor::FillESD(AliRunLoader* runLoader, AliESD* esd) const
{
/// Fill ESD
/// \todo add more

  TClonesArray* recTracksArray = 0;
  TClonesArray* recTrigTracksArray = 0;
  
  AliLoader* loader = runLoader->GetLoader("MUONLoader");
  loader->LoadTracks("READ");
  fMUONData->SetLoader(loader);

   // declaration  
  Int_t iEvent;// nPart;
  Int_t nTrackHits;// nPrimary;
  Double_t fitFmin;

  Double_t bendingSlope, nonBendingSlope, inverseBendingMomentum;
  Double_t xRec, yRec, zRec, chi2MatchTrigger;
  Bool_t matchTrigger;

  // setting pointer for tracks, triggertracks & trackparam at vertex
  AliMUONTrack* recTrack = 0;
  AliMUONTrackParam trackParam;
  AliMUONTriggerTrack* recTriggerTrack = 0;

  iEvent = runLoader->GetEventNumber(); 
  runLoader->GetEvent(iEvent);

  // Get vertex 
  Double_t vertex[3] = {0};
  const AliESDVertex *esdVert = esd->GetVertex(); 
  if (esdVert->GetNContributors()) {
    esdVert->GetXYZ(vertex);
    AliDebug(1, "find vertex\n");
  }
  // setting ESD MUON class
  AliESDMuonTrack* theESDTrack = new  AliESDMuonTrack() ;

  //-------------------- trigger tracks-------------
  Long_t trigPat = 0;
  fMUONData->SetTreeAddress("RL");
  fMUONData->GetRecTriggerTracks();
  recTrigTracksArray = fMUONData->RecTriggerTracks();

  // ready global trigger pattern from first track
  if (recTrigTracksArray) 
    recTriggerTrack = (AliMUONTriggerTrack*) recTrigTracksArray->First();
  if (recTriggerTrack) trigPat = recTriggerTrack->GetGTPattern();

  //printf(">>> Event %d Number of Recconstructed tracks %d \n",iEvent, nrectracks);
 
  // -------------------- tracks-------------
  fMUONData->SetTreeAddress("RT");
  fMUONData->GetRecTracks();
  recTracksArray = fMUONData->RecTracks();
        
  Int_t nRecTracks = 0;
  if (recTracksArray)
    nRecTracks = (Int_t) recTracksArray->GetEntriesFast(); //
  
  // loop over tracks
  for (Int_t iRecTracks = 0; iRecTracks <  nRecTracks;  iRecTracks++) {

    // reading info from tracks
    recTrack = (AliMUONTrack*) recTracksArray->At(iRecTracks);

    trackParam = *((AliMUONTrackParam*) (recTrack->GetTrackParamAtHit())->First());
   
    // extrapolate to the vertex if required
    //   if the vertex is not available, extrapolate to (0,0,0)
    if (!strstr(GetOption(),"NoExtrapToVtx")) {
      if (esdVert->GetNContributors())
        AliMUONTrackExtrap::ExtrapToVertex(&trackParam, vertex[0],vertex[1],vertex[2]);
      else
        AliMUONTrackExtrap::ExtrapToVertex(&trackParam, 0.,0.,0.);
    }
    
    bendingSlope            = trackParam.GetBendingSlope();
    nonBendingSlope         = trackParam.GetNonBendingSlope();
    inverseBendingMomentum  = trackParam.GetInverseBendingMomentum();
    xRec  = trackParam.GetNonBendingCoor();
    yRec  = trackParam.GetBendingCoor();
    zRec  = trackParam.GetZ();

    nTrackHits       = recTrack->GetNTrackHits();
    fitFmin          = recTrack->GetFitFMin();
    matchTrigger     = recTrack->GetMatchTrigger();
    chi2MatchTrigger = recTrack->GetChi2MatchTrigger();

    // setting data member of ESD MUON
    theESDTrack->SetInverseBendingMomentum(inverseBendingMomentum);
    theESDTrack->SetThetaX(TMath::ATan(nonBendingSlope));
    theESDTrack->SetThetaY(TMath::ATan(bendingSlope));
    theESDTrack->SetZ(zRec);
    theESDTrack->SetBendingCoor(yRec); // calculate vertex at ESD or Tracking level ?
    theESDTrack->SetNonBendingCoor(xRec);
    theESDTrack->SetChi2(fitFmin);
    theESDTrack->SetNHit(nTrackHits);
    theESDTrack->SetMatchTrigger(matchTrigger);
    theESDTrack->SetChi2MatchTrigger(chi2MatchTrigger);

    // storing ESD MUON Track into ESD Event 
    if (nRecTracks != 0)  
      esd->AddMuonTrack(theESDTrack);
  } // end loop tracks

  // reset muondata
  fMUONData->ResetRecTracks();
  fMUONData->ResetRecTriggerTracks();

  //} // end loop on event  
  loader->UnloadTracks(); 

  delete theESDTrack;
}

//_____________________________________________________________________________
void AliMUONReconstructor::FillESD(AliRunLoader* runLoader, AliRawReader* /*rawReader*/, AliESD* esd) const
{
/// Fill ESD
/// \todo add more

  // don't need rawReader ???
  FillESD(runLoader, esd);
}

//_____________________________________________________________________________
AliTracker* AliMUONReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
  /// create tracker for MUON
  /// go into the tracking framework finally (Ch.F)
 
  AliLoader* loader = runLoader->GetLoader("MUONLoader");

  fMUONData->SetLoader(loader);

  AliMUONTracker* tracker = new AliMUONTracker();
  tracker->SetMUONData(fMUONData);
  tracker->SetTriggerCircuit(fTriggerCircuit);
  tracker->SetOption(GetOption());

  return tracker;

}
