/*
 * AliAnalysisTaskNanoXioton.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */
#include "AliAnalysisTaskAODXioton.h"

ClassImp(AliAnalysisTaskAODXioton)

AliAnalysisTaskAODXioton::AliAnalysisTaskAODXioton()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fProtonMCList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fAntiProtonMCList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fXiMCList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fAntiXiMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskAODXioton::AliAnalysisTaskAODXioton(const char* name, bool isMC)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fProtonMCList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fAntiProtonMCList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fXiMCList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fAntiXiMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Xi Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiXi Cuts
  DefineOutput(6, TList::Class());  //Output for the Results
  DefineOutput(7, TList::Class());  //Output for the Results QA
  if (isMC) {
    DefineOutput(8, TList::Class());  //Output for the Track MC
    DefineOutput(9, TList::Class());  //Output for the Anti Track MC
    DefineOutput(10, TList::Class());  //Output for the V0 MC
    DefineOutput(11, TList::Class());  //Output for the Anti V0 MC
  }
}

AliAnalysisTaskAODXioton::~AliAnalysisTaskAODXioton() {
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
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
}

void AliAnalysisTaskAODXioton::UserCreateOutputObjects() {
  fGTI = new AliAODTrack *[fTrackBufferSize];

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
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 2,
                                                fConfig->GetMinimalBookingME());
  }
  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                  GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(
      fProton->GetIsMonteCarlo() || fAntiProton->GetIsMonteCarlo());

  fCascade = new AliFemtoDreamCascade();
  fCascade->SetUseMCInfo(fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  //PDG Codes should be set assuming Xi- to also work for Xi+
  fCascade->SetPDGCode(3312);
  fCascade->SetPDGDaugPos(2212);
  fCascade->GetPosDaug()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGDaugNeg(211);
  fCascade->GetNegDaug()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->SetPDGDaugBach(211);
  fCascade->GetBach()->SetUseMCInfo(
      fXi->GetIsMonteCarlo() || fAntiXi->GetIsMonteCarlo());
  fCascade->Setv0PDGCode(3122);

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  fProtonList = fProton->GetQAHists();
  fAntiProtonList = fAntiProton->GetQAHists();
  fXiList = fXi->GetQAHists();
  fAntiXiList = fAntiXi->GetQAHists();

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
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
  PostData(6, fResults);
  PostData(7, fResultsQA);

  if (fProton->GetIsMonteCarlo()) {
    if (!fProton->GetMinimalBooking()) {
      fProtonMCList = fProton->GetMCQAHists();
    } else {
      fProtonMCList = new TList();
      fProtonMCList->SetName("MCTrkCuts");
      fProtonMCList->SetOwner();
    }
    PostData(8, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    if (!fAntiProton->GetMinimalBooking()) {
      fAntiProtonMCList = fAntiProton->GetMCQAHists();
    } else {
      fAntiProtonMCList = new TList();
      fAntiProtonMCList->SetName("MCAntiTrkCuts");
      fAntiProtonMCList->SetOwner();
    }
    PostData(9, fAntiProtonMCList);
  }

  if (fXi->GetIsMonteCarlo()) {
    if (!fXi->GetMinimalBooking()) {
      fXiMCList = fXi->GetMCQAHists();
    } else {
      fXiMCList = new TList();
      fXiMCList->SetName("MCXiCuts");
      fXiMCList->SetOwner();
    }
    PostData(10, fXiMCList);
  }
  if (fAntiXi->GetIsMonteCarlo()) {
    if (!fAntiXi->GetMinimalBooking()) {
      fAntiXiMCList = fAntiXi->GetMCQAHists();
    } else {
      fAntiXiMCList = new TList();
      fAntiXiMCList->SetName("MCAntiv0Cuts");
      fAntiXiMCList->SetOwner();
    }
    PostData(11, fAntiXiMCList);
  }
}

void AliAnalysisTaskAODXioton::UserExec(Option_t *option) {
  AliAODEvent *evt = static_cast<AliAODEvent*>(InputEvent());
  if (!evt) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(evt);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }

  // PROTON SELECTION
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(evt->GetTrack(iTrack));
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
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(evt->GetTrack(iTrack));
    fTrack->SetTrack(track, multiplicity);
    if (fProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
    }
    if (fAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
    }
  }

  std::vector<AliFemtoDreamBasePart> Xis;
  std::vector<AliFemtoDreamBasePart> AntiXis;
  for (int iCasc = 0;
      iCasc < static_cast<TClonesArray *>(evt->GetCascades())->GetEntriesFast();
      ++iCasc) {
    AliAODcascade* casc = evt->GetCascade(iCasc);
    fCascade->SetCascade(evt, casc);
    if (fXi->isSelected(fCascade)) {
      Xis.push_back(*fCascade);
    }
    if (fAntiXi->isSelected(fCascade)) {
      AntiXis.push_back(*fCascade);
    }
  }
  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Protons, &Xis, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiXis, 1);

  fPairCleaner->CleanDecay(&Xis, 0);
  fPairCleaner->CleanDecay(&AntiXis, 1);

  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Xis);
  fPairCleaner->StoreParticle(AntiXis);

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  if (fProton->GetIsMonteCarlo()) {
    PostData(8, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    PostData(9, fAntiProtonMCList);
  }
  if (fXi->GetIsMonteCarlo()) {
    PostData(10, fXiMCList);
  }
  if (fAntiXi->GetIsMonteCarlo()) {
    PostData(11, fAntiXiMCList);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskAODXioton::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskAODXioton::StoreGlobalTrackReference(AliAODTrack *track) {
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
