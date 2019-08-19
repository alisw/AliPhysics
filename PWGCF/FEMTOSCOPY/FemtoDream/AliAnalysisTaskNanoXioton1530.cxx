/*
 * AliAnalysisTaskNanoXioton1530.cxx
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */
#include "AliAnalysisTaskNanoXioton1530.h"
#include "AliNanoAODTrack.h"
ClassImp(AliAnalysisTaskNanoXioton1530)
AliAnalysisTaskNanoXioton1530::AliAnalysisTaskNanoXioton1530()
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
      fPion(nullptr),
      fPionList(nullptr),
      fAntiPion(nullptr),
      fAntiPionList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fXi1530(nullptr),
      fXi1530Cuts(nullptr),
      fXi1530List(nullptr),
      fAntiXi1530Cuts(nullptr),
      fAntiXi1530List(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}



AliAnalysisTaskNanoXioton1530::AliAnalysisTaskNanoXioton1530(const char* name)
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
      fPion(nullptr),
      fPionList(nullptr),
      fAntiPion(nullptr),
      fAntiPionList(nullptr),
      fCascade(nullptr),
      fXi(nullptr),
      fXiList(nullptr),
      fAntiXi(nullptr),
      fAntiXiList(nullptr),
      fXi1530(nullptr),
      fXi1530Cuts(nullptr),
      fXi1530List(nullptr),
      fAntiXi1530Cuts(nullptr),
      fAntiXi1530List(nullptr),
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
  DefineOutput(4, TList::Class());  //Output for the Pion Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiPion Cuts
  DefineOutput(6, TList::Class());  //Output for the Xi Cuts
  DefineOutput(7, TList::Class());  //Output for the AntiXi Cuts
  DefineOutput(8, TList::Class());  //Output for the Xi1530 Cuts
  DefineOutput(9, TList::Class());  //Output for the AntiXi1530 Cuts
  DefineOutput(10, TList::Class()); //Output for the Results QA
  DefineOutput(11, TList::Class()); //Output for the Results QA
}

AliAnalysisTaskNanoXioton1530::~AliAnalysisTaskNanoXioton1530() {
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
  if (fPion) {
    delete fPion;
  }
  if (fAntiPion) {
    delete fAntiPion;
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
  if (fXi1530) {
    delete fXi1530;
  }
  if (fXi1530Cuts) {
    delete fXi1530Cuts;
  }
  if (fAntiXi1530Cuts) {
    delete fAntiXi1530Cuts;
  }
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
}

void AliAnalysisTaskNanoXioton1530::UserCreateOutputObjects() {
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
  if (!fPion) {
    AliError("No Pion cuts \n");
  } else {
    fPion->Init();
  }
  if (!fAntiPion) {
    AliError("No AntiPion cuts \n");
  } else {
    fAntiPion->Init();
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
  if (!fXi1530Cuts) {
    AliError("No Xi1530 cuts \n");
  } else {
    fXi1530Cuts->Init();
    fXi1530Cuts->SetName("Xi1530Cuts");
  }
  if (!fAntiXi1530Cuts) {
    AliError("No AntiXi1530 cuts \n");
  } else {
    fAntiXi1530Cuts->Init();
    fAntiXi1530Cuts->SetName("AntiXi1530Cuts");
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
                                  GetCollisionCandidates(), false);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(false);

  fXi1530 = new AliFemtoDreamv0();
  fXi1530->SetPDGCode(fXi1530Cuts->GetPDGv0());
  fXi1530->SetUseMCInfo(false);
  fXi1530->SetPDGDaughterPos(fXi1530Cuts->GetPDGPosDaug());  // order +sign doesnt play a role
  fXi1530->GetPosDaughter()->SetUseMCInfo(false);
  fXi1530->SetPDGDaughterNeg(fXi1530Cuts->GetPDGNegDaug());  // only used for MC Matching
  fXi1530->GetNegDaughter()->SetUseMCInfo(false);

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
  if (!fPion->GetMinimalBooking()) {
    fPionList = fPion->GetQAHists();
  } else {
    fPionList = new TList();
    fPionList->SetName("PionCuts");
    fPionList->SetOwner();
  }
  if (!fAntiPion->GetMinimalBooking()) {
    fAntiPionList = fAntiPion->GetQAHists();
  } else {
    fAntiPionList = new TList();
    fAntiPionList->SetName("AntiPionCuts");
    fAntiPionList->SetOwner();
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
  if (!fXi1530Cuts->GetMinimalBooking()) {
    fXi1530List = fXi1530Cuts->GetQAHists();
  } else {
    fXi1530List = new TList();
    fXi1530List->SetName("Xi1530Cuts");
    fXi1530List->SetOwner();
  }
  if (!fAntiXi1530Cuts->GetMinimalBooking()) {
    fAntiXi1530List = fAntiXi1530Cuts->GetQAHists();
  } else {
    fAntiXi1530List = new TList();
    fAntiXi1530List->SetName("AntiXi1530Cuts");
    fAntiXi1530List->SetOwner();
  }
  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");
  if (!fConfig->GetMinimalBookingME()) {
    fResults = fPartColl->GetHistList();
    fResultsQA->Add(fPartColl->GetQAList());
    fResultsQA->Add(fPairCleaner->GetHistList());
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
  PostData(6, fXiList);
  PostData(7, fAntiXiList);
  PostData(8, fXi1530List);
  PostData(9, fAntiXi1530List);
  PostData(10, fResults);
  PostData(11, fResultsQA);
}

void AliAnalysisTaskNanoXioton1530::UserExec(Option_t *option) {
//  AliVEvent *fInputEvent = InputEvent();
  if (!fInputEvent) {
    AliError("No input event");
    return;
  }
  fEvent->SetEvent(fInputEvent);
  if (!fEventCuts->isSelected(fEvent)) {
    return;
  }
  static const float piMass = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  static const float XiMass = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
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
  std::vector<AliFemtoDreamBasePart> PionPlus;
  std::vector<AliFemtoDreamBasePart> PionMinus;
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
    fTrack->SetInvMass(piMass);
    if (fPion->isSelected(fTrack)) {
      PionPlus.push_back(*fTrack);
    }
    if (fAntiPion->isSelected(fTrack)) {
      PionMinus.push_back(*fTrack);
    }
  }
  std::vector<AliFemtoDreamBasePart> Xis;
  std::vector<AliFemtoDreamBasePart> AntiXis;
  AliAODEvent* aodEvt = dynamic_cast<AliAODEvent*>(fInputEvent);
  for (int iCasc = 0;
      iCasc
          < static_cast<TClonesArray *>(aodEvt->GetCascades())->GetEntriesFast();
      ++iCasc) {
    AliAODcascade* casc = aodEvt->GetCascade(iCasc);
    fCascade->SetCascade(fInputEvent, casc);
    fCascade->SetInvMass(XiMass);
    if (fXi->isSelected(fCascade)) {
      Xis.push_back(*fCascade);
    }
    if (fAntiXi->isSelected(fCascade)) {
      AntiXis.push_back(*fCascade);
    }
  }
  std::vector<AliFemtoDreamBasePart> Xi1530s;
  std::vector<AliFemtoDreamBasePart> AntiXi1530s;
  fXi1530->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (const auto &pi : PionPlus) {
    for (const auto &xi : Xis) {
      fXi1530->Setv0(pi, xi, fInputEvent, false, true, fAntiXi1530Cuts->GetCutDaughters());
      if (fXi1530Cuts->isSelected(fXi1530)) {
        fXi1530->SetCPA(gRandom->Uniform());  //cpacode needed for CleanDecay v0;
        Xi1530s.push_back(*fXi1530);
      }
    }
  }
  for (const auto &pi : PionMinus) {
    for (const auto &xi : AntiXis) {
      fXi1530->Setv0(xi, pi, fInputEvent, true, false, fAntiXi1530Cuts->GetCutDaughters());
      if (fAntiXi1530Cuts->isSelected(fXi1530)) {
        fXi1530->SetCPA(gRandom->Uniform());  //cpacode needed for CleanDecay v0;
        AntiXi1530s.push_back(*fXi1530);
      }
    }
  }
  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Protons, &Xi1530s, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiXi1530s, 1);
  fPairCleaner->CleanDecay(&Xi1530s, 0);
  fPairCleaner->CleanDecay(&AntiXi1530s, 1);

  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Xi1530s);
  fPairCleaner->StoreParticle(AntiXi1530s);

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPionList);
  PostData(5, fAntiPionList);
  PostData(6, fXiList);
  PostData(7, fAntiXiList);
  PostData(8, fXi1530List);
  PostData(9, fAntiXi1530List);
  PostData(10, fResults);
  PostData(11, fResultsQA);
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoXioton1530::ResetGlobalTrackReference() {
  // see AliFemtoDreamAnalysis for details
  for (int i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskNanoXioton1530::StoreGlobalTrackReference(
    AliVTrack *track) {
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
