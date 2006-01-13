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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for MUON reconstruction                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <TParticle.h>
#include <TArrayF.h>

#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliESD.h"
#include "AliMUONReconstructor.h"
 
#include "AliMUONData.h"
#include "AliMUONTrackReconstructor.h"
#include "AliMUONClusterReconstructor.h"
#include "AliMUONClusterFinderVS.h"
#include "AliMUONClusterFinderAZ.h"
#include "AliMUONEventRecoCombi.h" 
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTriggerTrack.h"
#include "AliESDMuonTrack.h"
#include "AliMUONRawData.h"

#include "AliRawReader.h"


ClassImp(AliMUONReconstructor)
//_____________________________________________________________________________
AliMUONReconstructor::AliMUONReconstructor()
  : AliReconstructor()
{
}
//_____________________________________________________________________________
AliMUONReconstructor::~AliMUONReconstructor()
{
}
//_____________________________________________________________________________
void AliMUONReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
//  AliLoader
  AliLoader* loader = runLoader->GetLoader("MUONLoader");
  Int_t nEvents = runLoader->GetNumberOfEvents();

  AliMUONData* data = new AliMUONData(loader,"MUON","MUON");

// passing loader as argument.
  AliMUONTrackReconstructor* recoEvent = new AliMUONTrackReconstructor(loader, data);

  if (strstr(GetOption(),"Original")) 
    recoEvent->SetTrackMethod(1); // Original tracking
  else if (strstr(GetOption(),"Combi")) 
    recoEvent->SetTrackMethod(3); // Combined cluster / track
  else
    recoEvent->SetTrackMethod(2); // Kalman

  AliMUONClusterReconstructor* recoCluster = new AliMUONClusterReconstructor(loader, data);
  
  AliMUONClusterFinderVS *recModel = recoCluster->GetRecoModel();

  if (!strstr(GetOption(),"VS")) {
    recModel = (AliMUONClusterFinderVS*) new AliMUONClusterFinderAZ();
    recoCluster->SetRecoModel(recModel);
  }
  recModel->SetGhostChi2Cut(10);

  loader->LoadDigits("READ");
  loader->LoadRecPoints("RECREATE");
  loader->LoadTracks("RECREATE");
  
  Int_t chBeg = recoEvent->GetTrackMethod() == 3 ? 6 : 0; 
  //   Loop over events              
  for(Int_t ievent = 0; ievent < nEvents; ievent++) {
    printf("Event %d\n",ievent);
    runLoader->GetEvent(ievent);

    //----------------------- digit2cluster & Trigger2Trigger -------------------
    if (!loader->TreeR()) loader->MakeRecPointsContainer();
     
    // tracking branch
    if (recoEvent->GetTrackMethod() != 3) {
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
    recoCluster->Digits2Clusters(chBeg); 
    if (recoEvent->GetTrackMethod() == 3) {
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
  
    if (recoEvent->GetTrackMethod() == 3) { 
      // Combined cluster / track
      ((AliMUONClusterFinderAZ*) recModel)->SetReco(1);
      data->MakeBranch("RC");
      data->SetTreeAddress("RC");
      AliMUONEventRecoCombi::Instance()->FillRecP(data, recoEvent); 
      data->Fill("RC"); 
    }
    loader->WriteRecPoints("OVERWRITE"); 

    //--------------------------- Resetting branches -----------------------
    data->ResetDigits();
    data->ResetRawClusters();
    data->ResetTrigger();

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
}

//_____________________________________________________________________________
void AliMUONReconstructor::Reconstruct(AliRunLoader* runLoader, AliRawReader* rawReader) const
{
//  AliLoader
  AliLoader* loader = runLoader->GetLoader("MUONLoader");
  AliMUONData* data = new AliMUONData(loader,"MUON","MUON");

// passing loader as argument.
  AliMUONTrackReconstructor* recoEvent = new AliMUONTrackReconstructor(loader, data);

  AliMUONRawData* rawData = new AliMUONRawData(loader, data);

  AliMUONClusterReconstructor* recoCluster = new AliMUONClusterReconstructor(loader, data);
  AliMUONClusterFinderVS *recModel = recoCluster->GetRecoModel();
  recModel->SetGhostChi2Cut(10);

  loader->LoadRecPoints("RECREATE");
  loader->LoadTracks("RECREATE");
  loader->LoadDigits("RECREATE");


  //   Loop over events  
  Int_t iEvent = 0;
            
  while (rawReader->NextEvent()) {
    printf("Event %d\n",iEvent);
    runLoader->GetEvent(iEvent++);

    //----------------------- raw2digits & raw2trigger-------------------
    if (!loader->TreeD()) loader->MakeDigitsContainer();

    // tracking branch
    data->MakeBranch("D");
    data->SetTreeAddress("D");
    rawData->ReadTrackerDDL(rawReader);
    data->Fill("D"); 

    // trigger branch
    data->MakeBranch("GLT");
    data->SetTreeAddress("GLT");
    rawData->ReadTriggerDDL(rawReader);
    data->Fill("GLT"); 

    loader->WriteDigits("OVERWRITE");

    //----------------------- digit2cluster & Trigger2Trigger -------------------
    if (!loader->TreeR()) loader->MakeRecPointsContainer();
     
    // tracking branch
    data->MakeBranch("RC");
    data->SetTreeAddress("RC");
    recoCluster->Digits2Clusters(); 
    data->Fill("RC"); 

    // trigger branch
    data->MakeBranch("TC");
    data->SetTreeAddress("TC");
    recoCluster->Trigger2Trigger(); 
    data->Fill("TC");

    loader->WriteRecPoints("OVERWRITE");

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
  
    //--------------------------- Resetting branches -----------------------
    data->ResetDigits();
    data->ResetRawClusters();
    data->ResetTrigger();

    data->ResetRawClusters();
    data->ResetTrigger();
    data->ResetRecTracks();
    data->ResetRecTriggerTracks();
  
  }
  loader->UnloadRecPoints();
  loader->UnloadTracks();
  loader->UnloadDigits();

  delete recoCluster;
  delete recoEvent;
  delete data;
}

//_____________________________________________________________________________
void AliMUONReconstructor::FillESD(AliRunLoader* runLoader, AliESD* esd) const
{
  TClonesArray* recTracksArray = 0;
  TClonesArray* recTrigTracksArray = 0;
  
  AliLoader* loader = runLoader->GetLoader("MUONLoader");
  loader->LoadTracks("READ");
  AliMUONData* muonData = new AliMUONData(loader,"MUON","MUON");

   // declaration  
  Int_t iEvent;// nPart;
  Int_t nTrackHits;// nPrimary;
  Double_t fitFmin;
  TArrayF vertex(3);

  Double_t bendingSlope, nonBendingSlope, inverseBendingMomentum;
  Double_t xRec, yRec, zRec, chi2MatchTrigger;
  Bool_t matchTrigger;

  // setting pointer for tracks, triggertracks & trackparam at vertex
  AliMUONTrack* recTrack = 0;
  AliMUONTrackParam* trackParam = 0;
  AliMUONTriggerTrack* recTriggerTrack = 0;
//   TParticle* particle = new TParticle();
//   AliGenEventHeader* header = 0;
  iEvent = runLoader->GetEventNumber(); 
  runLoader->GetEvent(iEvent);

  // vertex calculation (maybe it exists already somewhere else)
  vertex[0] = vertex[1] = vertex[2] = 0.;
 //  nPrimary = 0;
//   if ( (header = runLoader->GetHeader()->GenEventHeader()) ) {
//     header->PrimaryVertex(vertex);
//   } else {
//     runLoader->LoadKinematics("READ");
//     runLoader->TreeK()->GetBranch("Particles")->SetAddress(&particle);
//     nPart = (Int_t)runLoader->TreeK()->GetEntries();
//     for(Int_t iPart = 0; iPart < nPart; iPart++) {
//       runLoader->TreeK()->GetEvent(iPart);
//       if (particle->GetFirstMother() == -1) {
// 	vertex[0] += particle->Vx();
// 	vertex[1] += particle->Vy();
// 	vertex[2] += particle->Vz();
// 	nPrimary++;
//       }
//       if (nPrimary) {
// 	vertex[0] /= (double)nPrimary;
// 	vertex[1] /= (double)nPrimary;
// 	vertex[2] /= (double)nPrimary;
//       }
//     }
//   }
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
    trackParam->ExtrapToVertex(vertex[0],vertex[1],vertex[2]);

    bendingSlope            = trackParam->GetBendingSlope();
    nonBendingSlope         = trackParam->GetNonBendingSlope();
    inverseBendingMomentum = trackParam->GetInverseBendingMomentum();
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

  // add global trigger pattern
  if (nRecTracks != 0)  
    esd->SetTrigger(trigPat);

  // reset muondata
  muonData->ResetRecTracks();
  muonData->ResetRecTriggerTracks();

  //} // end loop on event  
  loader->UnloadTracks(); 
 //  if (!header)
//     runLoader->UnloadKinematics();
  delete theESDTrack;
  delete muonData;
  // delete particle;
}//_____________________________________________________________________________
void AliMUONReconstructor::FillESD(AliRunLoader* runLoader, AliRawReader* /*rawReader*/, AliESD* esd) const
{
  // don't need rawReader ???
  FillESD(runLoader, esd);
}
