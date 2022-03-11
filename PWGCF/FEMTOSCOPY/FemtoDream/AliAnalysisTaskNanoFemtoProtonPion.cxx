/*
 * AliAnalysisTaskNanoFemtoProtonPion.cxx
 *
 *  Created on: 11 Mar 2022
 *  Author: Lesch Marcel
 */

#include "AliAnalysisTaskNanoFemtoProtonPion.h"
#include "AliAnalysisTaskNanoPt.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisTaskNanoFemtoProtonPion)
AliAnalysisTaskNanoFemtoProtonPion::AliAnalysisTaskNanoFemtoProtonPion()
  : AliAnalysisTaskSE(),
    fisLightWeight(false),
    fTrackBufferSize(),
    fIsMC(false),
    fDoPairCleaning(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsPion(nullptr),
    fTrackCutsAntiPion(nullptr),
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
    fPionList(nullptr),
    fPionMCList(nullptr),
    fAntiPionList(nullptr),
    fAntiPionMCList(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr){
}

AliAnalysisTaskNanoFemtoProtonPion::AliAnalysisTaskNanoFemtoProtonPion(
  const char *name, bool isMC)
  : AliAnalysisTaskSE(name),
    fisLightWeight(false),
    fTrackBufferSize(2000),
    fIsMC(isMC),
    fDoPairCleaning(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsPion(nullptr),
    fTrackCutsAntiPion(nullptr),
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
    fPionList(nullptr),
    fPionMCList(nullptr),
    fAntiPionList(nullptr),
    fAntiPionMCList(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr){
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Pion Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiPion Cuts
  DefineOutput(6, TList::Class());  //Output for the Results
  DefineOutput(7, TList::Class());  //Output for the Results QA
  if (fIsMC) {
    DefineOutput(8, TList::Class());  //Output for the Proton MC
    DefineOutput(9, TList::Class());  //Output for the AntiProton MC
    DefineOutput(10, TList::Class());  //Output for the Pion MC
    DefineOutput(11, TList::Class());  //Output for the AntiPion MC
  }
}

AliAnalysisTaskNanoFemtoProtonPion::~AliAnalysisTaskNanoFemtoProtonPion() {
  delete fEvent;
  delete fTrack;
  delete fTrackCutsPion; 
  delete fTrackCutsAntiPion;
  delete fTrackCutsProton;
  delete fTrackCutsAntiProton;
  delete fPairCleaner;
  delete fPartColl;
}

void AliAnalysisTaskNanoFemtoProtonPion::UserCreateOutputObjects() {

  fGTI = new AliVTrack*[fTrackBufferSize];

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

  if (!fTrackCutsPion) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsPion->Init();
    fPionList = fTrackCutsPion->GetQAHists();
    if (fIsMC) {
      fPionMCList = fTrackCutsPion->GetMCQAHists();
    }
  }

  if (!fTrackCutsAntiPion) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsAntiPion->Init();
    fAntiPionList = fTrackCutsAntiPion->GetQAHists();
    if (fIsMC) {
      fAntiPionMCList = fTrackCutsAntiPion->GetMCQAHists();
    }
  }
//--------------------------------------------------------------------------------------------------------------------
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

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPionList);
  PostData(5, fAntiPionList);
  PostData(6, fResults);
  PostData(7, fResultsQA);

  if (fTrackCutsProton->GetIsMonteCarlo()) {
    PostData(8, fProtonMCList);
  }
  if (fTrackCutsAntiProton->GetIsMonteCarlo()) {
    PostData(9, fAntiProtonMCList);
  }

  if (fTrackCutsPion->GetIsMonteCarlo()) {
    PostData(10, fPionMCList);
  }
  if (fTrackCutsAntiPion->GetIsMonteCarlo()) {
    PostData(11, fAntiPionMCList);
  }
}

void AliAnalysisTaskNanoFemtoProtonPion::UserExec(Option_t*) {
  AliVEvent *Event= fInputEvent;

  if (!Event) {
    AliWarning("No Input Event");
  } else {
    fEvent->SetEvent(Event);
    if (fEventCuts->isSelected(fEvent)) {
      ResetGlobalTrackReference();
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliVTrack *track = static_cast<AliVTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard NanoAOD");
          return;
        }
        StoreGlobalTrackReference(track);
      }

      fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
      static std::vector<AliFemtoDreamBasePart> SelectedPions; 
      SelectedPions.clear();
      static std::vector<AliFemtoDreamBasePart> SelectedAntiPions; 
      SelectedAntiPions.clear();
      static std::vector<AliFemtoDreamBasePart> SelectedProtons;
      SelectedProtons.clear();
      static std::vector<AliFemtoDreamBasePart> SelectedAntiProtons;
      SelectedAntiProtons.clear();

      //Now we loop over all the tracks in the reconstructed event.
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliVTrack *track = static_cast<AliVTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }

        fTrack->SetTrack(track, fInputEvent);
        if (fTrackCutsProton->isSelected(fTrack)) {
          SelectedProtons.push_back(*fTrack);
        }
        if (fTrackCutsAntiProton->isSelected(fTrack)) {
          SelectedAntiProtons.push_back(*fTrack);
        }
        if (fTrackCutsPion->isSelected(fTrack)){ 
            SelectedPions.push_back(*fTrack);
        }
        if (fTrackCutsAntiPion->isSelected(fTrack)){
            SelectedAntiPions.push_back(*fTrack);
        }
      }
      //loop once over the MC stack to calculate Efficiency/Purity
      if (fIsMC) {
        AliAODInputHandler *eventHandler =
          dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()
                                            ->GetInputEventHandler());
        AliMCEvent* fMC = eventHandler->MCEvent();

        for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
          AliAODMCParticle *mcPart = (AliAODMCParticle*) fMC->GetTrack(iPart);
          if (mcPart->IsPhysicalPrimary()) {
            if (mcPart->GetPdgCode() == fTrackCutsProton->GetPDGCode()) {
              fTrackCutsProton->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fTrackCutsAntiProton->GetPDGCode()) {
              fTrackCutsAntiProton->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fTrackCutsPion->GetPDGCode()) {
              fTrackCutsPion->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fTrackCutsAntiPion->GetPDGCode()) {
              fTrackCutsAntiPion->FillGenerated(mcPart->Pt());
            }
          }
        }
      }

      if(fDoPairCleaning){
        fPairCleaner->CleanTrackAndDecay(&SelectedProtons, &SelectedPions, 0); 
        fPairCleaner->CleanTrackAndDecay(&SelectedAntiProtons, &SelectedAntiPions, 0); 
      }

      fPairCleaner->ResetArray();
      fPairCleaner->StoreParticle(SelectedProtons);
      fPairCleaner->StoreParticle(SelectedAntiProtons);
      fPairCleaner->StoreParticle(SelectedPions);
      fPairCleaner->StoreParticle(SelectedAntiPions);
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                          fEvent->GetZVertex(), fEvent->GetRefMult08(),
                          fEvent->GetV0MCentrality());
    }
  }

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPionList);
  PostData(5, fAntiPionList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
//-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fTrackCutsProton->GetIsMonteCarlo()) {
    PostData(8, fProtonMCList);
  }
  if (fTrackCutsAntiProton->GetIsMonteCarlo()) {
    PostData(9, fAntiProtonMCList);
  }
  if (fTrackCutsPion->GetIsMonteCarlo()) {
    PostData(10, fPionMCList);
  }
  if (fTrackCutsAntiPion->GetIsMonteCarlo()) {
    PostData(11, fAntiPionMCList);
  }
}

void AliAnalysisTaskNanoFemtoProtonPion::ResetGlobalTrackReference() {

  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskNanoFemtoProtonPion::StoreGlobalTrackReference(AliVTrack *track) {
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack*>(track);
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
    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}

