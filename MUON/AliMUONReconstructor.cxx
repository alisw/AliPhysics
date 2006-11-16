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

#include "AliESD.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliMUON.h"
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
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrackReconstructor.h"
#include "AliMUONTrackReconstructorK.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerCrateStore.h"

#include "AliMpSegmentation.h"

#include "AliMUONPreClusterFinder.h"
#include "AliMUONClusterFinderCOG.h"
#include "AliMUONClusterFinderSimpleFit.h"
#include "AliMUONClusterFinderMLEM.h"
  
#include "AliRawReader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "TTask.h"
#include "TStopwatch.h"

/// \cond CLASSIMP
ClassImp(AliMUONReconstructor)
/// \endcond

//_____________________________________________________________________________
AliMUONReconstructor::AliMUONReconstructor()
  : AliReconstructor(), 
    fRunLoader(0x0),
    fDigitMaker(new AliMUONDigitMaker()), 
    fCalibrationData(0x0),
    fCrateManager(new AliMUONTriggerCrateStore()),
    fTriggerCircuit(new TClonesArray("AliMUONTriggerCircuit", 234)),
    fTransformer(new AliMUONGeometryTransformer(kTRUE))

{
/// Default constructor

    AliDebug(1,"");
    // Crate manager
    fCrateManager->ReadFromFile();

    // set to digit maker
    fDigitMaker->SetCrateManager(fCrateManager);

    // transformater
    fTransformer->ReadGeometryData("volpath.dat", "geometry.root");

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
}

//_____________________________________________________________________________
TTask* 
AliMUONReconstructor::GetCalibrationTask(AliMUONData* data) const
{
/// Create the calibration task(s). 
  
  const AliRun* run = fRunLoader->GetAliRun();

  AliInfo("Calibration will occur.");
  Int_t runNumber = run->GetRunNumber();     
  fCalibrationData = new AliMUONCalibrationData(runNumber);
  if ( !fCalibrationData->IsValid() )
    {
      AliError("Could not retrieve calibrations !");
      delete fCalibrationData;
      fCalibrationData = 0x0;
      return 0x0;
    }    
  TTask* calibration = new TTask("MUONCalibrator","MUON Digit calibrator");
  calibration->Add(new AliMUONDigitCalibrator(data,fCalibrationData));
  //FIXME: calibration->Add(something about dead channels should go here).
  return calibration;

}

//_____________________________________________________________________________
void
AliMUONReconstructor::Init(AliRunLoader* runLoader)
{
/// Initialize

  fRunLoader = runLoader;
}

//_____________________________________________________________________________
AliMUONClusterReconstructor*
AliMUONReconstructor::CreateClusterReconstructor(AliMUONData* data) const
{
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
    new AliMUONClusterReconstructor(data,clusterFinder,fTransformer);
  return clusterReco;
}

//_____________________________________________________________________________
void AliMUONReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
/// Reconstruct
/// \todo add more

  AliLoader* loader = runLoader->GetLoader("MUONLoader");
  Int_t nEvents = runLoader->GetNumberOfEvents();

  AliMUONData* data = new AliMUONData(loader,"MUON","MUON");

// passing loader as argument.
  AliMUONVTrackReconstructor* recoEvent;
  if (strstr(GetOption(),"Original")) recoEvent = new AliMUONTrackReconstructor(data);
  else if (strstr(GetOption(),"Combi")) recoEvent = new AliMUONTrackReconstructorK(data,"Combi");
  else recoEvent = new AliMUONTrackReconstructorK(data,"Kalman");
  
  recoEvent->SetTriggerCircuit(fTriggerCircuit);

  AliMUONClusterReconstructor* recoCluster = CreateClusterReconstructor(data);
  
  AliMUONClusterFinderVS *recModel = recoCluster->GetRecoModel();

  if (!strstr(GetOption(),"VS")) {
    recModel = (AliMUONClusterFinderVS*) new AliMUONClusterFinderAZ();
    recoCluster->SetRecoModel(recModel);
  }
  recModel->SetGhostChi2Cut(10);

  loader->LoadDigits("READ");
  loader->LoadRecPoints("RECREATE");
  loader->LoadTracks("RECREATE");
  
  TTask* calibration = GetCalibrationTask(data);
  
  Int_t chBeg = (strstr(GetOption(),"Combi") ? 6 : 0);

  //   Loop over events              
  for(Int_t ievent = 0; ievent < nEvents; ievent++) {

    AliDebug(1,Form("Event %d",ievent));
    
    runLoader->GetEvent(ievent);

    //----------------------- digit2cluster & Trigger2Trigger -------------------
    if (!loader->TreeR()) loader->MakeRecPointsContainer();
     
    // tracking branch
    if (!strstr(GetOption(),"Combi")) {
      data->MakeBranch("RC");
      data->SetTreeAddress("D,RC");
    } else {
      data->SetTreeAddress("D");
      data->SetTreeAddress("RCC");
    }
    // Important for avoiding a memory leak when reading digits ( to be investigated more in detail)
    // In any case the reading of GLT is needed for the Trigger2Tigger method below
    data->SetTreeAddress("GLT");

    data->GetDigits();
    
    if ( calibration ) 
    {
      calibration->ExecuteTask();
    }
    
    recoCluster->Digits2Clusters(chBeg); 
    
    if (strstr(GetOption(),"Combi")) {
      // Combined cluster / track finder
      AliMUONEventRecoCombi::Instance()->FillEvent(data, (AliMUONClusterFinderAZ*)recModel);
      ((AliMUONClusterFinderAZ*) recModel)->SetReco(2); 
    }
    else data->Fill("RC"); 

    // trigger branch
    data->MakeBranch("TC");
    data->SetTreeAddress("TC");
    recoCluster->Trigger2Trigger(); 
    data->Fill("TC");

    //AZ loader->WriteRecPoints("OVERWRITE");

    //---------------------------- Track & TriggerTrack ---------------------
    if (!loader->TreeT()) loader->MakeTracksContainer();

    // trigger branch
    data->MakeBranch("RL"); //trigger track
    data->SetTreeAddress("RL");
    recoEvent->EventReconstructTrigger();
    data->Fill("RL");

    // tracking branch
    data->MakeBranch("RT"); //track
    data->SetTreeAddress("RT");
    recoEvent->EventReconstruct();
    data->Fill("RT");

    loader->WriteTracks("OVERWRITE"); 
  
    if (strstr(GetOption(),"Combi")) { 
      // Combined cluster / track
      ((AliMUONClusterFinderAZ*) recModel)->SetReco(1);
      data->MakeBranch("RC");
      data->SetTreeAddress("RC");
      AliMUONEventRecoCombi::Instance()->FillRecP(data, (AliMUONTrackReconstructorK*)recoEvent); 
      data->Fill("RC"); 
    }
    loader->WriteRecPoints("OVERWRITE"); 

    //--------------------------- Resetting branches -----------------------
    data->ResetDigits();
    data->ResetRawClusters();
    data->ResetTrigger();
    data->ResetRecTracks();  
    data->ResetRecTriggerTracks();

  }
  loader->UnloadDigits();
  loader->UnloadRecPoints();
  loader->UnloadTracks();

  delete recoCluster;
  delete recoEvent;
  delete data;
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
  AliMUONData data(loader,"MUON","MUON");

  // passing loader as argument.
  fDigitMaker->SetMUONData(&data);

  AliMUONClusterReconstructor* recoCluster = CreateClusterReconstructor(&data);

  AliMUONVTrackReconstructor *recoEvent;
  if (strstr(GetOption(),"Original")) recoEvent = new AliMUONTrackReconstructor(&data);
  else if (strstr(GetOption(),"Combi")) recoEvent = new AliMUONTrackReconstructorK(&data,"Combi");
  else recoEvent = new AliMUONTrackReconstructorK(&data,"Kalman");

  recoEvent->SetTriggerCircuit(fTriggerCircuit);

  AliMUONClusterFinderVS *recModel = recoCluster->GetRecoModel();

  if (!strstr(GetOption(),"VS")) 
  {
    recModel = (AliMUONClusterFinderVS*) new AliMUONClusterFinderAZ();
    recoCluster->SetRecoModel(recModel);
  }
  recModel->SetGhostChi2Cut(10);

  TTask* calibration = GetCalibrationTask(&data);
  
  loader->LoadRecPoints("RECREATE");
  loader->LoadTracks("RECREATE");
  loader->LoadDigits("READ");
  
  //   Loop over events  
  Int_t iEvent = 0;
           
  TStopwatch totalTimer;
  TStopwatch rawTimer;
  TStopwatch calibTimer;
  TStopwatch clusterTimer;
  TStopwatch trackingTimer;
  
  rawTimer.Start(kTRUE); rawTimer.Stop();
  calibTimer.Start(kTRUE); calibTimer.Stop();
  clusterTimer.Start(kTRUE); clusterTimer.Stop();
  trackingTimer.Start(kTRUE); trackingTimer.Stop();
  
  totalTimer.Start(kTRUE);
  
  while (rawReader->NextEvent()) 
  {
    AliDebug(1,Form("Event %d",iEvent));
    
    runLoader->GetEvent(iEvent++);

    //----------------------- raw2digits & raw2trigger-------------------
    if (!loader->TreeD()) 
    {
      AliDebug(1,Form("Making Digit Container for event %d",iEvent));
      loader->MakeDigitsContainer();
    }
    
    data.SetTreeAddress("D,GLT");
    rawTimer.Start(kFALSE);
    fDigitMaker->Raw2Digits(rawReader);
    rawTimer.Stop();
    
    if ( calibration )
    {
      calibTimer.Start(kFALSE);
      calibration->ExecuteTask();
      calibTimer.Stop();
    }
  
    //----------------------- digit2cluster & Trigger2Trigger -------------------
    clusterTimer.Start(kFALSE);

    if (!loader->TreeR()) loader->MakeRecPointsContainer();
     
    // tracking branch
    data.MakeBranch("RC");
    data.SetTreeAddress("RC");
    recoCluster->Digits2Clusters(); 
    data.Fill("RC"); 

    // trigger branch
    data.MakeBranch("TC");
    data.SetTreeAddress("TC");
    recoCluster->Trigger2Trigger();
    data.Fill("TC");
    
    loader->WriteRecPoints("OVERWRITE");

    clusterTimer.Stop();

    //---------------------------- Track & TriggerTrack ---------------------
    trackingTimer.Start(kFALSE);
    if (!loader->TreeT()) loader->MakeTracksContainer();

    // trigger branch
    data.MakeBranch("RL"); //trigger track
    data.SetTreeAddress("RL");
    recoEvent->EventReconstructTrigger();
    data.Fill("RL");

    // tracking branch
    data.MakeBranch("RT"); //track
    data.SetTreeAddress("RT");
    recoEvent->EventReconstruct();
    data.Fill("RT");

    loader->WriteTracks("OVERWRITE");  
    trackingTimer.Stop();
    
    //--------------------------- Resetting branches -----------------------
    data.ResetDigits();
    data.ResetRawClusters();
    data.ResetTrigger();
    data.ResetRecTracks();
    data.ResetRecTriggerTracks();
  
  }
  
  totalTimer.Stop();
  
  loader->UnloadRecPoints();
  loader->UnloadTracks();
  loader->UnloadDigits();
  
  delete recoEvent;

  delete recoCluster;
  
  AliInfo(Form("Execution time for converting RAW data to digits in MUON : R:%.2fs C:%.2fs",
               rawTimer.RealTime(),rawTimer.CpuTime()));
  AliInfo(Form("Execution time for calibrating MUON : R:%.2fs C:%.2fs",
               calibTimer.RealTime(),calibTimer.CpuTime()));
  AliInfo(Form("Execution time for clusterizing MUON : R:%.2fs C:%.2fs",
               clusterTimer.RealTime(),clusterTimer.CpuTime()));
  AliInfo(Form("Execution time for tracking MUON : R:%.2fs C:%.2fs",
               trackingTimer.RealTime(),trackingTimer.CpuTime()));
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
  AliMUONData* muonData = new AliMUONData(loader,"MUON","MUON");

   // declaration  
  Int_t iEvent;// nPart;
  Int_t nTrackHits;// nPrimary;
  Double_t fitFmin;

  Double_t bendingSlope, nonBendingSlope, inverseBendingMomentum;
  Double_t xRec, yRec, zRec, chi2MatchTrigger;
  Bool_t matchTrigger;

  // setting pointer for tracks, triggertracks & trackparam at vertex
  AliMUONTrack* recTrack = 0;
  AliMUONTrackParam* trackParam = 0;
  AliMUONTriggerTrack* recTriggerTrack = 0;

  iEvent = runLoader->GetEventNumber(); 
  runLoader->GetEvent(iEvent);

  // Get vertex 
  Double_t vertex[3] = {0};
  const AliESDVertex *esdVert = esd->GetVertex(); 
  if (esdVert->GetNContributors()) {
    esdVert->GetXYZ(vertex);
    printf("find vertex\n");
  }
  // setting ESD MUON class
  AliESDMuonTrack* theESDTrack = new  AliESDMuonTrack() ;

  //-------------------- trigger tracks-------------
  Long_t trigPat = 0;
  muonData->SetTreeAddress("RL");
  muonData->GetRecTriggerTracks();
  recTrigTracksArray = muonData->RecTriggerTracks();

  // ready global trigger pattern from first track
  if (recTrigTracksArray) 
    recTriggerTrack = (AliMUONTriggerTrack*) recTrigTracksArray->First();
  if (recTriggerTrack) trigPat = recTriggerTrack->GetGTPattern();

  //printf(">>> Event %d Number of Recconstructed tracks %d \n",iEvent, nrectracks);
 
  // -------------------- tracks-------------
  muonData->SetTreeAddress("RT");
  muonData->GetRecTracks();
  recTracksArray = muonData->RecTracks();
        
  Int_t nRecTracks = 0;
  if (recTracksArray)
    nRecTracks = (Int_t) recTracksArray->GetEntriesFast(); //
  
  // loop over tracks
  for (Int_t iRecTracks = 0; iRecTracks <  nRecTracks;  iRecTracks++) {

    // reading info from tracks
    recTrack = (AliMUONTrack*) recTracksArray->At(iRecTracks);

    trackParam = (AliMUONTrackParam*) (recTrack->GetTrackParamAtHit())->First();
   
    if (esdVert->GetNContributors())
      AliMUONTrackExtrap::ExtrapToVertex(trackParam, vertex[0],vertex[1],vertex[2]);

    bendingSlope            = trackParam->GetBendingSlope();
    nonBendingSlope         = trackParam->GetNonBendingSlope();
    inverseBendingMomentum  = trackParam->GetInverseBendingMomentum();
    xRec  = trackParam->GetNonBendingCoor();
    yRec  = trackParam->GetBendingCoor();
    zRec  = trackParam->GetZ();

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
  muonData->ResetRecTracks();
  muonData->ResetRecTriggerTracks();

  //} // end loop on event  
  loader->UnloadTracks(); 

  delete theESDTrack;
  delete muonData;
}//_____________________________________________________________________________
void AliMUONReconstructor::FillESD(AliRunLoader* runLoader, AliRawReader* /*rawReader*/, AliESD* esd) const
{
/// Fill ESD
/// \todo add more

  // don't need rawReader ???
  FillESD(runLoader, esd);
}
