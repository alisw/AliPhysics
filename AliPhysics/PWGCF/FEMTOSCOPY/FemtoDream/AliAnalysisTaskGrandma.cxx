/*
 * AliAnalysisTaskGrandma.cxx
 *
 *  Created on: Sep 11, 2018
 *      Author: hohlweger
 */
#include <vector>
#include "AliAnalysisTaskGrandma.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "TClonesArray.h"
#include "AliLog.h"
ClassImp(AliAnalysisTaskGrandma)

AliAnalysisTaskGrandma::AliAnalysisTaskGrandma()
    : AliAnalysisTaskSE(),
      fTrackBufferSize(2000),
      fIsMC(false),
      fQA(nullptr),
      fMinBookingME(false),
      fMinBookingSample(false),
      fMVPileUp(false),
      fEvtCutQA(false),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fEvtHistList(nullptr),
      fFemtoTrack(nullptr),
      fTrackCuts(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCuts(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr),
      fFemtov0(nullptr),
      fv0Cuts(nullptr),
      fv0CutHistList(nullptr),
      fv0CutHistMCList(nullptr),
      fAntiv0Cuts(nullptr),
      fAntiv0CutHistList(nullptr),
      fAntiv0CutHistMCList(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fConfig(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr),
      fGTI(nullptr) {
}

AliAnalysisTaskGrandma::AliAnalysisTaskGrandma(const char* name, bool isMC)
    : AliAnalysisTaskSE(name),
      fTrackBufferSize(2000),
      fIsMC(isMC),
      fQA(nullptr),
      fMinBookingME(false),
      fMinBookingSample(false),
      fMVPileUp(false),
      fEvtCutQA(false),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fEvtHistList(nullptr),
      fFemtoTrack(nullptr),
      fTrackCuts(nullptr),
      fTrackCutHistList(nullptr),
      fTrackCutHistMCList(nullptr),
      fAntiTrackCuts(nullptr),
      fAntiTrackCutHistList(nullptr),
      fAntiTrackCutHistMCList(nullptr),
      fFemtov0(nullptr),
      fv0Cuts(nullptr),
      fv0CutHistList(nullptr),
      fv0CutHistMCList(nullptr),
      fAntiv0Cuts(nullptr),
      fAntiv0CutHistList(nullptr),
      fAntiv0CutHistMCList(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fConfig(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Class and Pair Cleaner
  DefineOutput(2, TList::Class());  //Output for the Event Cuts
  DefineOutput(3, TList::Class());  //Output for the Track Cuts
  DefineOutput(4, TList::Class());  //Output for the Anti Track Cuts
  DefineOutput(5, TList::Class());  //Output for the V0 Cuts
  DefineOutput(6, TList::Class());  //Output for the Anti V0 Cuts
  DefineOutput(7, TList::Class());  //Output for the Results
  DefineOutput(8, TList::Class());  //Output for the Results QA
  if (fIsMC){
    DefineOutput(9, TList::Class());  //Output for the MC Track Cuts
    DefineOutput(10, TList::Class());  //Output for the MC V0 Cuts
    DefineOutput(11, TList::Class());  //Output for the MC Anti Track Cuts
    DefineOutput(12, TList::Class());  //Output for the MC Anti V0 Cuts
  }
}

AliAnalysisTaskGrandma::~AliAnalysisTaskGrandma() {
}

void AliAnalysisTaskGrandma::UserCreateOutputObjects() {
  //first flag: Pile up rejection based on multiple vertices
  // (not used anymore since we event selection is done by AliEventCuts)
  //second flag: for the QA output of the AliEventCuts
  // might want to turn this off for systematics

  fEvent = new AliFemtoDreamEvent(true, fEvtCutQA, GetCollisionCandidates());
  fEvent->SetCalcSpherocity(true);

  fFemtoTrack = new AliFemtoDreamTrack();
  fFemtoTrack->SetUseMCInfo(fIsMC);

  fFemtov0 = new AliFemtoDreamv0();
  fFemtov0->SetPDGCode(fv0Cuts->GetPDGv0());
  fFemtov0->SetUseMCInfo(fIsMC);
  fFemtov0->SetPDGDaughterPos(fv0Cuts->GetPDGPosDaug());  //order +sign doesnt play a role
  fFemtov0->GetPosDaughter()->SetUseMCInfo(fIsMC);
  fFemtov0->SetPDGDaughterNeg(fv0Cuts->GetPDGNegDaug());  //only used for MC Matching
  fFemtov0->GetNegDaughter()->SetUseMCInfo(fIsMC);

  //the pair cleaner you have to setup yourself depending on the pairs you want to investigate
  fPairCleaner = new AliFemtoDreamPairCleaner(2, 3, false);  //false - full booking, true - minimal booking

  fQA = new TList();
  fQA->SetOwner();
  fQA->SetName("QA");
  fQA->Add(fEvent->GetEvtCutList());
  fQA->Add(fPairCleaner->GetHistList());

  // if (!fEvtCuts->GetMinimalBooking()) {
  //   if (fAnalysis->GetEventCutHists()) {
  //     fEvtHistList = fAnalysis->GetEventCutHists();
  //   }
  // } else {
  //   fEvtHistList = new TList();
  //   fEvtHistList->SetName("EventCuts");
  //   fEvtHistList->SetOwner();
  // }

  if (fEvtCuts) {
    fEvtCuts->InitQA();
    if (fEvtCuts->GetHistList()) {
      fEvtHistList = fEvtCuts->GetHistList();
    }
  } else {
    AliWarning("Event cuts are missing! \n");
  }

  fGTI = new AliAODTrack*[fTrackBufferSize];
  if (fTrackCuts) {
    fTrackCuts->Init();
    if (fTrackCuts->GetQAHists()) {
      fTrackCutHistList = fTrackCuts->GetQAHists();
    }
    if (fTrackCuts->GetMCQAHists()) {
      fTrackCutHistMCList = fTrackCuts->GetMCQAHists();
    }
  } else {
    AliWarning("Track cuts are missing! \n");
  }

  if (fAntiTrackCuts) {
    fAntiTrackCuts->Init();
    if (fAntiTrackCuts->GetQAHists()) {
      fAntiTrackCutHistList = fAntiTrackCuts->GetQAHists();
    }
    if (fAntiTrackCuts->GetMCQAHists()) {
      fAntiTrackCutHistMCList = fAntiTrackCuts->GetMCQAHists();
    }
  } else {
    AliWarning("Anti Track cuts are missing! \n");
  }

if (fv0Cuts) {
    fv0Cuts->Init();
    if (fv0Cuts->GetQAHists()) {
      fv0CutHistList = fv0Cuts->GetQAHists();
    }
    if (fv0Cuts->GetMCQAHists()) {
      fv0CutHistMCList = fv0Cuts->GetMCQAHists();
    }
  } else {
    AliWarning("V0 cuts are missing! \n");
  }
  if (fAntiv0Cuts) {
    fAntiv0Cuts->Init();
    if (fAntiv0Cuts->GetQAHists()) {
      fAntiv0CutHistList = fAntiv0Cuts->GetQAHists();
    }
    if (fAntiv0Cuts->GetMCQAHists()) {
      fAntiv0CutHistMCList = fAntiv0Cuts->GetMCQAHists();
    }
  } else {
    AliWarning("Anti V0 cuts are missing! \n");
  }

  if (fConfig->GetUseEventMixing()) {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    fResultList = fPartColl->GetHistList();
    fResultQAList = fPartColl->GetQAList();
  }
  PostData(1, fQA);
  PostData(2, fEvtHistList);
  PostData(3, fTrackCutHistList);
  PostData(4, fAntiTrackCutHistList);
  PostData(5, fv0CutHistList);
  PostData(6, fAntiv0CutHistList);
  PostData(7, fResultList);
  PostData(8, fResultQAList);
  if (fIsMC){
    PostData(9, fTrackCutHistMCList);
    PostData(10, fv0CutHistMCList);
    PostData(11, fAntiTrackCutHistMCList);
    PostData(12, fAntiv0CutHistMCList);
  }
}
void AliAnalysisTaskGrandma::UserExec(Option_t *) {
  AliAODEvent *Event = static_cast<AliAODEvent*>(fInputEvent);
  if (!Event) {
    AliWarning("No Input Event");
  } else {

    fEvent->SetEvent(Event);
    if (fEvtCuts->isSelected(fEvent)) {
      ResetGlobalTrackReference();
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
        StoreGlobalTrackReference(track);
      }


      std::vector<AliFemtoDreamBasePart> Particles;
      std::vector<AliFemtoDreamBasePart> AntiParticles;
      fFemtoTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
      for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
        AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
        fFemtoTrack->SetTrack(track, fEvent->GetMultiplicity());
        if (fTrackCuts->isSelected(fFemtoTrack)) {
          Particles.push_back(*fFemtoTrack);
        }
        if (fAntiTrackCuts->isSelected(fFemtoTrack)) {
          AntiParticles.push_back(*fFemtoTrack);
        }
      }

     std::vector<AliFemtoDreamBasePart> Decays;
     std::vector<AliFemtoDreamBasePart> AntiDecays;
      //  Look for the lambda, store it in an event
      //  Get a V0 from the event:
      TClonesArray *v01 = static_cast<TClonesArray*>(Event->GetV0s());
      //number of V0s:

      fFemtov0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
      int entriesV0 = v01->GetEntriesFast();
      for (int iv0 = 0; iv0 < entriesV0; iv0++) {
        AliAODv0 *v0 = Event->GetV0(iv0);
        fFemtov0->Setv0(Event, v0, fEvent->GetMultiplicity());
        if (fv0Cuts->isSelected(fFemtov0)) {
          Decays.push_back(*fFemtov0);
        }
        if (fAntiv0Cuts->isSelected(fFemtov0)) {
          AntiDecays.push_back(*fFemtov0);
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
          if (TMath::Abs(mcPart->Eta()) < 0.8 && mcPart->IsPhysicalPrimary()) {
            if (mcPart->GetPdgCode() == fTrackCuts->GetPDGCode()) {
              fTrackCuts->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fAntiTrackCuts->GetPDGCode()) {
              fAntiTrackCuts->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fv0Cuts->GetPDGv0()) {
              fv0Cuts->FillGenerated(mcPart->Pt());
            } else if (mcPart->GetPdgCode() == fAntiv0Cuts->GetPDGv0()) {
              fAntiv0Cuts->FillGenerated(mcPart->Pt());
            }
          }
        }
      }


      fPairCleaner->ResetArray();
      fPairCleaner->CleanTrackAndDecay(&Particles, &AntiDecays, 0);
      //  fPairCleaner->CleanTrackAndDecay(&Particles, &XiDecays, 2);
      fPairCleaner->CleanTrackAndDecay(&AntiParticles, &Decays, 1);

      //  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiXiDecays, 3);
      //
        fPairCleaner->CleanDecay(&Decays, 0);
        fPairCleaner->CleanDecay(&AntiDecays, 1);
        fPairCleaner->CleanDecayAndDecay(&Decays, &AntiDecays, 2);

      //  fPairCleaner->CleanDecay(&XiDecays, 2);
      //  fPairCleaner->CleanDecay(&AntiXiDecays, 3);

      fPairCleaner->StoreParticle(Particles);
      fPairCleaner->StoreParticle(AntiParticles);
      fPairCleaner->StoreParticle(Decays);
      fPairCleaner->StoreParticle(AntiDecays);
      //  fPairCleaner->StoreParticle(XiDecays);
      //  fPairCleaner->StoreParticle(AntiXiDecays);
      if (fConfig->GetUseEventMixing()) {
        fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                            fEvent->GetZVertex(), fEvent->GetMultiplicity(),
                            fEvent->GetV0MCentrality());
      }
    }
  }
  PostData(1, fQA);
  PostData(2, fEvtHistList);
  PostData(3, fTrackCutHistList);
  PostData(4, fAntiTrackCutHistList);
  PostData(5, fv0CutHistList);
  PostData(6, fAntiv0CutHistList);
  PostData(7, fResultList);
  PostData(8, fResultQAList);
  if (fIsMC){
    PostData(9, fTrackCutHistMCList);
    PostData(10, fv0CutHistMCList);
    PostData(11, fAntiTrackCutHistMCList);
    PostData(12, fAntiv0CutHistMCList);
  }
  return;
}
void AliAnalysisTaskGrandma::StoreGlobalTrackReference(AliAODTrack *track) {
  //This method was inherited form H. Beck analysis

  //bhohlweg@cern.ch: We ask for the Unique Track ID that points back to the
  //ESD. Seems like global tracks have a positive ID, Tracks with Filterbit
  //128 only have negative ID, this is used to match the Tracks later to their
  //global counterparts

  // Stores the pointer to the global track

  // This was AOD073
  // // Don't use the filter bits 2 (ITS standalone) and 128 TPC only
  // // Remove this return statement and you'll see
  // // any TPC signal
  // if(track->TestFilterBit(128) || track->TestFilterBit(2))
  //   return;
  // This is AOD086
  // Another set of tracks was introduced: Global constrained.
  // We only want filter bit 1 <-- NO! we also want no
  // filter bit at all, which are the v0 tracks
  //  if(!track->TestFilterBit(1))
  //    return;

  // There are also tracks without any filter bit, i.e. filter map 0,
  // at the beginning of the event: they have ~id 1 to 5, 1 to 12
  // This are tracks that didn't survive the primary track filter but
  // got written cause they are V0 daughters

  // Check whether the track has some info
  // I don't know: there are tracks with filter bit 0
  // and no TPC signal. ITS standalone V0 daughters?
  // if(!track->GetTPCsignal()){
  //   printf("Warning: track has no TPC signal, "
  //     //    "not adding it's info! "
  //     "ID: %d FilterMap: %d\n"
  //     ,track->GetID(),track->GetFilterMap());
  //   //    return;
  // }

  // Check that the id is positive
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }

  // Check id is not too big for buffer
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  // Warn if we overwrite a track
  if (fGTI[trackID]) {
    // Seems like there are FilterMap 0 tracks
    // that have zero TPCNcls, don't store these!
    if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    // Imagine the other way around, the zero map zero clusters track
    // is stored and the good one wants to be added. We ommit the warning
    // and just overwrite the 'bad' track
    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }  // Two tracks same id

  // // There are tracks with filter bit 0,
  // // do they have TPCNcls stored?
  // if(!track->GetFilterMap()){
  //   printf("Filter map is zero, TPCNcls: %u\n"
  //     ,track->GetTPCNcls());
  // }

  // Assign the pointer
  (fGTI[trackID]) = track;
}

void AliAnalysisTaskGrandma::ResetGlobalTrackReference() {
  //This method was inherited form H. Beck analysis

  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}
