/*
 * AliAnalysisTaskGrandma.cxx
 *
 *  Created on: Sep 11, 2018
 *      Author: hohlweger
 */

#include "AliAnalysisTaskGrandma.h"
#include "AliLog.h"
ClassImp(AliAnalysisTaskGrandma)

AliAnalysisTaskGrandma::AliAnalysisTaskGrandma()
    : AliAnalysisTaskSE(),
      fIsMC(false),
      fQA(nullptr),
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
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fConfig(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskGrandma::AliAnalysisTaskGrandma(const char* name, bool isMC)
    : AliAnalysisTaskSE(name),
      fIsMC(isMC),
      fQA(nullptr),
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
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fConfig(nullptr),
      fResultList(nullptr),
      fResultQAList(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Class and Pair Cleaner
  DefineOutput(2, TList::Class());  //Output for the Event Cuts
  DefineOutput(3, TList::Class());  //Output for the Track Cuts
  DefineOutput(4, TList::Class());  //Output for the Anti Track Cuts
  DefineOutput(5, TList::Class());  //Output for the Results
  DefineOutput(6, TList::Class());  //Output for the Results QA
  if (fIsMC)
    DefineOutput(7, TList::Class());  //Output for the MC Track Cuts
  if (fIsMC)
    DefineOutput(8, TList::Class());  //Output for the MC Anti Track Cuts
}

AliAnalysisTaskGrandma::~AliAnalysisTaskGrandma() {

}

void AliAnalysisTaskGrandma::UserCreateOutputObjects() {
  //first flag: Pile up rejection based on multiple vertices
  // (not used anymore since we event selection is done by AliEventCuts)
  //second flag: for the QA output of the AliEventCuts
  // might want to turn this off for systematics
  fEvent = new AliFemtoDreamEvent(true, true, GetCollisionCandidates());
  fFemtoTrack = new AliFemtoDreamTrack();
  fFemtoTrack->SetUseMCInfo(fIsMC);
  //the pair cleaner you have to setup yourself depending on the pairs you want to investigate
  fPairCleaner = new AliFemtoDreamPairCleaner(0, 0, false);  //false - full booking, true - minimal booking

  fQA = new TList();
  fQA->SetOwner();
  fQA->SetName("QA");
  fQA->Add(fEvent->GetEvtCutList());
  fQA->Add(fPairCleaner->GetHistList());

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
  PostData(5, fResultList);
  PostData(6, fResultQAList);
  if (fIsMC)
    PostData(7, fTrackCutHistMCList);
  if (fIsMC)
    PostData(8, fAntiTrackCutHistMCList);
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
      fPairCleaner->ResetArray();
      //  fPairCleaner->CleanTrackAndDecay(&Particles, &Decays, 0);
      //  fPairCleaner->CleanTrackAndDecay(&Particles, &XiDecays, 2);
      //  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiDecays, 1);
      //  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiXiDecays, 3);
      //
      //  fPairCleaner->CleanDecay(&Decays, 0);
      //  fPairCleaner->CleanDecay(&AntiDecays, 1);
      //  fPairCleaner->CleanDecay(&XiDecays, 2);
      //  fPairCleaner->CleanDecay(&AntiXiDecays, 3);

      fPairCleaner->StoreParticle(Particles);
      fPairCleaner->StoreParticle(AntiParticles);
//        fPairCleaner->StoreParticle(Decays);
//        fPairCleaner->StoreParticle(AntiDecays);
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
  PostData(5, fResultList);
  PostData(6, fResultQAList);
  if (fIsMC)
    PostData(7, fTrackCutHistMCList);
  if (fIsMC)
    PostData(8, fAntiTrackCutHistMCList);
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
