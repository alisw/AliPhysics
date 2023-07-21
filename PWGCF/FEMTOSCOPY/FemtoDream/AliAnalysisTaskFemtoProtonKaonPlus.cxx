/*
 * AliAnalysisTaskFemtoProtonKaonPlus.cxx
 *
 *  Created on: 2 May 2023
 *  Author: Italiano Luca
 */

#include "AliAnalysisTaskFemtoProtonKaonPlus.h"
#include "AliAnalysisTaskNanoPt.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliNanoAODTrack.h"
#include "AliAODTrack.h"
#include <sstream>

ClassImp(AliAnalysisTaskFemtoProtonKaonPlus)
AliAnalysisTaskFemtoProtonKaonPlus::AliAnalysisTaskFemtoProtonKaonPlus()
  : AliAnalysisTaskSE(),
    fisLightWeight(false),
    fTrackBufferSize(),
    fIsMC(false),
    fDoPairCleaning(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsKaon(nullptr),
    fTrackCutsAntiKaon(nullptr),
    fTrackCutsProton(nullptr),
    fTrackCutsAntiProton(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fKaonList(nullptr),
    fKaonMCList(nullptr),
    fAntiKaonList(nullptr),
    fAntiKaonMCList(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr){
}

AliAnalysisTaskFemtoProtonKaonPlus::AliAnalysisTaskFemtoProtonKaonPlus(
  const char *name, bool isMC)
  : AliAnalysisTaskSE(name),
    fisLightWeight(false),
    fTrackBufferSize(2000),
    fIsMC(isMC),
    fDoPairCleaning(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsKaon(nullptr),
    fTrackCutsAntiKaon(nullptr),
    fTrackCutsProton(nullptr),
    fTrackCutsAntiProton(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fKaonList(nullptr),
    fKaonMCList(nullptr),
    fAntiKaonList(nullptr),
    fAntiKaonMCList(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr){
    
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Kaon Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiKaon Cuts
  DefineOutput(6, TList::Class());  //Output for the Results
  DefineOutput(7, TList::Class());  //Output for the Results QA 
  if (fIsMC) {
    DefineOutput(8, TList::Class());  //Output for the Proton MC
    DefineOutput(9, TList::Class());  //Output for the AntiProton MC
    DefineOutput(10, TList::Class());  //Output for the Kaon MC
    DefineOutput(11, TList::Class());  //Output for the AntiKaon MC
  }
}

AliAnalysisTaskFemtoProtonKaonPlus::~AliAnalysisTaskFemtoProtonKaonPlus() {
  delete fEvent;
  delete fTrack;
  delete fTrackCutsKaon; 
  delete fTrackCutsAntiKaon;
  delete fTrackCutsProton;
  delete fTrackCutsAntiProton;
  delete fPairCleaner;
  delete fPartColl;
}

//==================================================================================================================================================

void AliAnalysisTaskFemtoProtonKaonPlus::UserCreateOutputObjects() {

  fGTI = new AliAODTrack*[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }

  if (!fTrackCutsProton) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsProton->Init();
    fProtonList = fTrackCutsProton->GetQAHists();
    if (fIsMC) {
    fProtonMCList = fTrackCutsProton->GetMCQAHists();
    }
  }

  if (!fTrackCutsAntiProton) {
    AliError("No AntiProton cuts \n");
  } else {
    fTrackCutsAntiProton->Init();
    fAntiProtonList = fTrackCutsAntiProton->GetQAHists();
    if (fIsMC) {
      fAntiProtonMCList = fTrackCutsAntiProton->GetMCQAHists();
    }
  }

  if (!fTrackCutsKaon) {
    AliError("No Kaon cuts \n");
  } else {
    fTrackCutsKaon->Init();
    fKaonList = fTrackCutsKaon->GetQAHists();
    if (fIsMC) {
      fKaonMCList = fTrackCutsKaon->GetMCQAHists();
    }
  }

  if (!fTrackCutsAntiKaon) {
    AliError("No AntiKaon cuts \n");
  } else {
    fTrackCutsAntiKaon->Init();
    fAntiKaonList = fTrackCutsAntiKaon->GetQAHists();
    if (fIsMC) {
      fAntiKaonMCList = fTrackCutsAntiKaon->GetMCQAHists();
    }
  }
  //////////////////////////////////////////////////////////////////////////
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
        fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 0,
        fConfig->GetMinimalBookingME()); 
  }

  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());
  //fEvent->SetCalcSpherocity(fEventCuts->GetDoSpherocityCuts()); 
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");

  if (fConfig->GetUseEventMixing()) {
    fResults = fPartColl->GetHistList();
    if (!fConfig->GetMinimalBookingME()) {
      fResultsQA->Add(fPartColl->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  } else {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }

//Deleted section on ancestor studies here

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fKaonList);
  PostData(5, fAntiKaonList);
  PostData(6, fResults);
  PostData(7, fResultsQA); 

  if (fTrackCutsProton->GetIsMonteCarlo()) {
    PostData(8, fProtonMCList);
  }
  if (fTrackCutsAntiProton->GetIsMonteCarlo()) {
    PostData(9, fAntiProtonMCList);
  }

  if (fTrackCutsKaon->GetIsMonteCarlo()) {
    PostData(10, fKaonMCList);
  }
  if (fTrackCutsAntiKaon->GetIsMonteCarlo()) {
    PostData(11, fAntiKaonMCList);
  }

} //void AliAnalysisTaskFemtoProtonKaonPlus::UserCreateOutputObjects()

//==================================================================================================================================================

void AliAnalysisTaskFemtoProtonKaonPlus::UserExec(Option_t*) {
  AliAODEvent *Event = static_cast<AliAODEvent*>(InputEvent());

  if (!Event) {
    AliWarning("No Input Event");
    return;
  } 

  fEvent->SetEvent(Event);
  if (!fEventCuts->isSelected(fEvent)) {
       return;
  }

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard NanoAOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }

  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  static std::vector<AliFemtoDreamBasePart> SelectedKaons; 
  SelectedKaons.clear();
  static std::vector<AliFemtoDreamBasePart> SelectedAntiKaons; 
  SelectedAntiKaons.clear();
  static std::vector<AliFemtoDreamBasePart> SelectedProtons;
  SelectedProtons.clear();
  static std::vector<AliFemtoDreamBasePart> SelectedAntiProtons;
  SelectedAntiProtons.clear();

  //Now we loop over all the tracks in the reconstructed event.
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
    if (!track) {
      continue;
    }

    fTrack->SetTrack(track);

    /*bool isProton = fTrackCutsProton->isSelected(fTrack);
    bool isAntiProton = fTrackCutsAntiProton->isSelected(fTrack);
    bool isKaon = fTrackCutsKaon->isSelected(fTrack);
    bool isAntiKaon = fTrackCutsAntiKaon->isSelected(fTrack);*/

    if (fTrackCutsProton->isSelected(fTrack) && fTrackCutsKaon->isSelected(fTrack)){
      continue;
    }
    if (fTrackCutsAntiProton->isSelected(fTrack) && fTrackCutsAntiKaon->isSelected(fTrack)){
      continue;
    }

    if (fTrackCutsProton->isSelected(fTrack)) {
      SelectedProtons.push_back(*fTrack);
    }
    if (fTrackCutsAntiProton->isSelected(fTrack)) {
      SelectedAntiProtons.push_back(*fTrack);
    }
    if (fTrackCutsKaon->isSelected(fTrack)){ 
      SelectedKaons.push_back(*fTrack);
    }
    if (fTrackCutsAntiKaon->isSelected(fTrack)){
      SelectedAntiKaons.push_back(*fTrack);
    }
  }

  //loop once over the MC stack to calculate Efficiency/Purity
  if (fIsMC) {
  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliMCEvent* fMC = eventHandler->MCEvent();

  for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
      AliAODMCParticle *mcPart = (AliAODMCParticle*) fMC->GetTrack(iPart);
      if (mcPart->IsPhysicalPrimary()) {
        if (mcPart->GetPdgCode() == fTrackCutsProton->GetPDGCode()) {
          fTrackCutsProton->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fTrackCutsAntiProton->GetPDGCode()) {
          fTrackCutsAntiProton->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fTrackCutsKaon->GetPDGCode()) {
          fTrackCutsKaon->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fTrackCutsAntiKaon->GetPDGCode()) {
          fTrackCutsAntiKaon->FillGenerated(mcPart->Pt());
        }
      }
    }
  }

  fPairCleaner->ResetArray();

  if(fDoPairCleaning){
    fPairCleaner->CleanTrackAndDecay(&SelectedProtons, &SelectedKaons, 0); 
    fPairCleaner->CleanTrackAndDecay(&SelectedAntiProtons, &SelectedAntiKaons, 1); 
  }

  fPairCleaner->StoreParticle(SelectedProtons);
  fPairCleaner->StoreParticle(SelectedAntiProtons);
  fPairCleaner->StoreParticle(SelectedKaons);
  fPairCleaner->StoreParticle(SelectedAntiKaons);

  //Official FemtoDream Two-Body Calculations
  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fKaonList);
  PostData(5, fAntiKaonList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
//-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fTrackCutsProton->GetIsMonteCarlo()) {
    PostData(8, fProtonMCList);
  }
  if (fTrackCutsAntiProton->GetIsMonteCarlo()) {
    PostData(9, fAntiProtonMCList);
  }
  if (fTrackCutsKaon->GetIsMonteCarlo()) {
    PostData(10, fKaonMCList);
  }
  if (fTrackCutsAntiKaon->GetIsMonteCarlo()) {
    PostData(11, fAntiKaonMCList);
  }
}//void AliAnalysisTaskFemtoProtonKaonPlus::UserExec(Option_t*)

//==================================================================================================================================================

void AliAnalysisTaskFemtoProtonKaonPlus::ResetGlobalTrackReference() {

  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//==================================================================================================================================================

void AliAnalysisTaskFemtoProtonKaonPlus::StoreGlobalTrackReference(AliAODTrack *track){
  // see AliFemtoDreamAnalysis for details
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  if (fGTI[trackID]) {
    if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if ((fGTI[trackID])->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}
