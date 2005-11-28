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

// used local container for each method
// passing fLoader as argument, could be avoided ???
  AliMUONTrackReconstructor* recoEvent = new AliMUONTrackReconstructor(loader);
  AliMUONData* dataEvent = recoEvent->GetMUONData();
  if (strstr(GetOption(),"Kalman")) recoEvent->SetTrackMethod(2); // Kalman
  else if (strstr(GetOption(),"Combi")) recoEvent->SetTrackMethod(3); // Combined cluster / track
  else recoEvent->SetTrackMethod(1); // original

  AliMUONClusterReconstructor* recoCluster = new AliMUONClusterReconstructor(loader);
  AliMUONData* dataCluster = recoCluster->GetMUONData();
  AliMUONClusterFinderVS *recModel = recoCluster->GetRecoModel();
  if (strstr(GetOption(),"AZ") || strstr(GetOption(),"Combi")) {
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
      dataCluster->MakeBranch("RC");
      dataCluster->SetTreeAddress("D,RC");
    } else {
      dataCluster->SetTreeAddress("D");
      dataCluster->SetTreeAddress("RCC");
    }
    // Important for avoiding a memory leak when reading digits ( to be investigated more in detail)
    // In any case the reading of GLT is needed for the Trigger2Tigger method below
    dataCluster->SetTreeAddress("GLT");

    recoCluster->Digits2Clusters(chBeg); 
    if (recoEvent->GetTrackMethod() == 3) {
      // Combined cluster / track finder
      AliMUONEventRecoCombi::Instance()->FillEvent(dataCluster, dataEvent, (AliMUONClusterFinderAZ*)recModel);
      ((AliMUONClusterFinderAZ*) recModel)->SetReco(2); 
    }
    else dataCluster->Fill("RC"); 

    // trigger branch
    dataCluster->MakeBranch("TC");
    dataCluster->SetTreeAddress("TC");
    recoCluster->Trigger2Trigger(); 
    dataCluster->Fill("TC");

    //AZ loader->WriteRecPoints("OVERWRITE");

    //---------------------------- Track & TriggerTrack ---------------------
    if (!loader->TreeT()) loader->MakeTracksContainer();

    // trigger branch
    dataEvent->MakeBranch("RL"); //trigger track
    dataEvent->SetTreeAddress("RL");
    recoEvent->EventReconstructTrigger();
    dataEvent->Fill("RL");

    // tracking branch
    dataEvent->MakeBranch("RT"); //track
    dataEvent->SetTreeAddress("RT");
    recoEvent->EventReconstruct();
    dataEvent->Fill("RT");

    loader->WriteTracks("OVERWRITE"); 
  
    if (recoEvent->GetTrackMethod() == 3) { 
      // Combined cluster / track
      ((AliMUONClusterFinderAZ*) recModel)->SetReco(1);
      dataCluster->MakeBranch("RC");
      dataCluster->SetTreeAddress("RC");
      AliMUONEventRecoCombi::Instance()->FillRecP(dataCluster, recoEvent); 
      dataCluster->Fill("RC"); 
    }
    loader->WriteRecPoints("OVERWRITE"); 

    //--------------------------- Resetting branches -----------------------
    dataCluster->ResetDigits();
    dataCluster->ResetRawClusters();
    dataCluster->ResetTrigger();

    dataEvent->ResetRawClusters();
    dataEvent->ResetTrigger();
    dataEvent->ResetRecTracks();  
    dataEvent->ResetRecTriggerTracks();

  }
  loader->UnloadDigits();
  loader->UnloadRecPoints();
  loader->UnloadTracks();

  delete recoCluster;
  delete recoEvent;
}

//_____________________________________________________________________________
void AliMUONReconstructor::Reconstruct(AliRunLoader* runLoader, AliRawReader* rawReader) const
{
//  AliLoader
  AliLoader* loader = runLoader->GetLoader("MUONLoader");

// used local container for each method
// passing fLoader as argument, could be avoided ???
  AliMUONTrackReconstructor* recoEvent = new AliMUONTrackReconstructor(loader);
  AliMUONData* dataEvent = recoEvent->GetMUONData();

  AliMUONRawData* rawData = new AliMUONRawData(loader);
  AliMUONData* dataCluster = rawData->GetMUONData();

  AliMUONClusterReconstructor* recoCluster = new AliMUONClusterReconstructor(loader, dataCluster);
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
    dataCluster->MakeBranch("D");
    dataCluster->SetTreeAddress("D");
    rawData->ReadTrackerDDL(rawReader);
    dataCluster->Fill("D"); 

    // trigger branch
    dataCluster->MakeBranch("GLT");
    dataCluster->SetTreeAddress("GLT");
    rawData->ReadTriggerDDL(rawReader);
    dataCluster->Fill("GLT"); 

    loader->WriteDigits("OVERWRITE");

    //----------------------- digit2cluster & Trigger2Trigger -------------------
    if (!loader->TreeR()) loader->MakeRecPointsContainer();
     
    // tracking branch
    dataCluster->MakeBranch("RC");
    dataCluster->SetTreeAddress("RC");
    recoCluster->Digits2Clusters(); 
    dataCluster->Fill("RC"); 

    // trigger branch
    dataCluster->MakeBranch("TC");
    dataCluster->SetTreeAddress("TC");
    recoCluster->Trigger2Trigger(); 
    dataCluster->Fill("TC");

    loader->WriteRecPoints("OVERWRITE");

    //---------------------------- Track & TriggerTrack ---------------------
    if (!loader->TreeT()) loader->MakeTracksContainer();

    // trigger branch
    dataEvent->MakeBranch("RL"); //trigger track
    dataEvent->SetTreeAddress("RL");
    recoEvent->EventReconstructTrigger();
    dataEvent->Fill("RL");

    // tracking branch
    dataEvent->MakeBranch("RT"); //track
    dataEvent->SetTreeAddress("RT");
    recoEvent->EventReconstruct();
    dataEvent->Fill("RT");

    loader->WriteTracks("OVERWRITE");  
  
    //--------------------------- Resetting branches -----------------------
    dataCluster->ResetDigits();
    dataCluster->ResetRawClusters();
    dataCluster->ResetTrigger();

    dataEvent->ResetRawClusters();
    dataEvent->ResetTrigger();
    dataEvent->ResetRecTracks();
    dataEvent->ResetRecTriggerTracks();
  
  }
  loader->UnloadRecPoints();
  loader->UnloadTracks();
  loader->UnloadDigits();

  delete recoCluster;
  delete recoEvent;
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
