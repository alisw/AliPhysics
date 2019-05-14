/*
 * AliAnalysisTaskNanoXioton.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */
#include "AliAnalysisTaskNanoXioton.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisTaskNanoXioton)
AliAnalysisTaskNanoXioton::AliAnalysisTaskNanoXioton()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskNanoXioton::AliAnalysisTaskNanoXioton(const char* name)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Xi Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiXi Cuts
}

AliAnalysisTaskNanoXioton::~AliAnalysisTaskNanoXioton() {
  if (fEvent) {
    delete fEvent;
  }
  if (fEventCuts) {
    delete fEventCuts;
  }
  if (fTrack) {
    delete fTrack;
  }
  if (fProton) {
    delete fProton;
  }
  if (fAntiProton) {
    delete fAntiProton;
  }
  if (fCascade) {
    delete fCascade;
  }
  if (fXi) {
    delete fXi;
  }
  if (fAntiXi) {
    delete fAntiXi;
  }
}

void AliAnalysisTaskNanoXioton::UserCreateOutputObjects() {
  fGTI = new AliVTrack *[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }
  if (!fProton) {
    AliError("No Proton cuts \n");
  } else {
    fProton->Init();
  }
  if (!fAntiProton) {
    AliError("No AntiProton cuts \n");
  } else {
    fAntiProton->Init();
  }
  if (!fXi) {
    AliError("No Xi cuts \n");
  } else {
    fXi->Init();
  }
  if (!fAntiXi) {
    AliError("No AntiXi cuts \n");
  } else {
    fAntiXi->Init();
  }
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), false);
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(false);

  fCascade = new AliFemtoDreamCascade();
  fCascade->SetUseMCInfo(false);
  //PDG Codes should be set assuming Xi- to also work for Xi+
  fCascade->SetPDGCode(3312);
  fCascade->SetPDGDaugPos(2212);
  fCascade->GetPosDaug()->SetUseMCInfo(false);
  fCascade->SetPDGDaugNeg(211);
  fCascade->GetNegDaug()->SetUseMCInfo(false);
  fCascade->SetPDGDaugBach(211);
  fCascade->GetBach()->SetUseMCInfo(false);
  fCascade->Setv0PDGCode(3122);

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }
  if (!fProton->GetMinimalBooking()) {
    fProtonList = fProton->GetQAHists();
  } else {
    fProtonList = new TList();
    fProtonList->SetName("TrackCuts");
    fProtonList->SetOwner();
  }
  if (!fAntiProton->GetMinimalBooking()) {
    fAntiProtonList = fAntiProton->GetQAHists();
  } else {
    fAntiProtonList = new TList();
    fAntiProtonList->SetName("AntiTrackCuts");
    fAntiProtonList->SetOwner();
  }
  if (!fXi->GetMinimalBooking()) {
    fXiList = fXi->GetQAHists();
  } else {
    fXiList = new TList();
    fXiList->SetName("XiCuts");
    fXiList->SetOwner();
  }
  if (!fAntiXi->GetMinimalBooking()) {
    fAntiXiList = fAntiXi->GetQAHists();
  } else {
    fAntiXiList = new TList();
    fAntiXiList->SetName("AntiXiCuts");
    fAntiXiList->SetOwner();
  }
  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
}

void AliAnalysisTaskNanoXioton::UserExec(Option_t *option) {
  AliVEvent *fInputEvent = InputEvent();
  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }

  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  std::vector<AliFemtoDreamBasePart> Protons;
  std::vector<AliFemtoDreamBasePart> AntiProtons;
  const int multiplicity = fEvent->GetMultiplicity();
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    fTrack->SetTrack(track, fInputEvent, multiplicity);
    if (fProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
    }
    if (fAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
    }
  }
  std::vector<AliFemtoDreamBasePart> Xis;
  std::vector<AliFemtoDreamBasePart> AntiXis;
  if (fInputEvent->GetNumberOfCascades() != 0) {
    std::cout << " +=================================================== \n";
    std::cout << " +=================================================== \n";
    std::cout << "fInput->GetCascades(): " << fInputEvent->GetNumberOfCascades() << std::endl;
    std::cout << " +=================================================== \n";
    std::cout << " +=================================================== \n";
  }
  AliAODEvent* aodEvt = dynamic_cast<AliAODEvent*>(fInputEvent);
  for (int iCasc = 0;
      iCasc
          < static_cast<TClonesArray *>(aodEvt->GetCascades())->GetEntriesFast();
      ++iCasc) {
    AliAODcascade* casc = aodEvt->GetCascade(iCasc);
    fCascade->SetCascade(fInputEvent, casc);
    if (fXi->isSelected(fCascade)) {
      Xis.push_back(*fCascade);
    }
    if (fAntiXi->isSelected(fCascade)) {
      AntiXis.push_back(*fCascade);
    }
  }

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoXioton::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoXioton::StoreGlobalTrackReference(AliVTrack *track) {
  // see AliFemtoDreamAnalysis for details
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
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
    if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}
