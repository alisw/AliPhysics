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
      fProtonProtonDump(nullptr),
      fAntiProtonAntiProtonDump(nullptr),
      fDumpster(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskNanoXioton::AliAnalysisTaskNanoXioton(const char* name,
                                                     bool isMC)
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
      fProtonProtonDump(nullptr),
      fAntiProtonAntiProtonDump(nullptr),
      fDumpster(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Xi Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiXi Cuts
  DefineOutput(6, TList::Class());  //Output for the Results
  DefineOutput(7, TList::Class());  //Output for the Results QA
  DefineOutput(8, TList::Class());  //Output for the Dumpster
  if (isMC) {
    DefineOutput(9, TList::Class());  //Output for the Track MC
    DefineOutput(10, TList::Class());  //Output for the Anti Track MC
    DefineOutput(11, TList::Class());  //Output for the V0 MC
    DefineOutput(12, TList::Class());  //Output for the Anti V0 MC
  }
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
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }

  if (fProtonProtonDump) {
    delete fProtonProtonDump;
  }
  if (fAntiProtonAntiProtonDump) {
    delete fAntiProtonAntiProtonDump;
  }
  if (fDumpster) {
    delete fDumpster;
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

  fDumpster = new TList();
  fDumpster->SetName("Dumpster");
  fDumpster->SetOwner(kTRUE);

  if (fConfig->GetUseEventMixing()) {
    fResults = fPartColl->GetHistList();
    if (!fConfig->GetMinimalBookingME()) {
      fResultsQA->Add(fPartColl->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());

      fProtonProtonDump = new AliFemtoDreamDump("pp");
      fProtonProtonDump->SetkstarThreshold(0.1);
      fDumpster->Add(fProtonProtonDump->GetOutput());

      fAntiProtonAntiProtonDump = new AliFemtoDreamDump("apap");
      fAntiProtonAntiProtonDump->SetkstarThreshold(0.1);
      fDumpster->Add(fAntiProtonAntiProtonDump->GetOutput());
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
  PostData(8, fDumpster);

  if (fProton->GetIsMonteCarlo()) {
    if (!fProton->GetMinimalBooking()) {
      fProtonMCList = fProton->GetMCQAHists();
    } else {
      fProtonMCList = new TList();
      fProtonMCList->SetName("MCTrkCuts");
      fProtonMCList->SetOwner();
    }
    PostData(9, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    if (!fAntiProton->GetMinimalBooking()) {
      fAntiProtonMCList = fAntiProton->GetMCQAHists();
    } else {
      fAntiProtonMCList = new TList();
      fAntiProtonMCList->SetName("MCAntiTrkCuts");
      fAntiProtonMCList->SetOwner();
    }
    PostData(10, fAntiProtonMCList);
  }

  if (fXi->GetIsMonteCarlo()) {
    if (!fXi->GetMinimalBooking()) {
      fXiMCList = fXi->GetMCQAHists();
    } else {
      fXiMCList = new TList();
      fXiMCList->SetName("MCXiCuts");
      fXiMCList->SetOwner();
    }
    PostData(11, fXiMCList);
  }
  if (fAntiXi->GetIsMonteCarlo()) {
    if (!fAntiXi->GetMinimalBooking()) {
      fAntiXiMCList = fAntiXi->GetMCQAHists();
    } else {
      fAntiXiMCList = new TList();
      fAntiXiMCList->SetName("MCAntiv0Cuts");
      fAntiXiMCList->SetOwner();
    }
    PostData(12, fAntiXiMCList);
  }
}

void AliAnalysisTaskNanoXioton::UserExec(Option_t *option) {
//  AliVEvent *fInputEvent = InputEvent();
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
  fPairCleaner->ResetArray();
  fPairCleaner->CleanTrackAndDecay(&Protons, &Xis, 0);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiXis, 1);

  fPairCleaner->CleanDecay(&Xis, 0);
  fPairCleaner->CleanDecay(&AntiXis, 1);

  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Xis);
  fPairCleaner->StoreParticle(AntiXis);

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());

  if(fProtonProtonDump) {
    fProtonProtonDump->SetEvent(Protons, fEvent, 2212);
  }
  if (fAntiProtonAntiProtonDump) {
    fAntiProtonAntiProtonDump->SetEvent(AntiProtons, fEvent, 2212);
  }

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fXiList);
  PostData(5, fAntiXiList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  PostData(8, fDumpster);
  if (fProton->GetIsMonteCarlo()) {
    PostData(9, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    PostData(10, fAntiProtonMCList);
  }
  if (fXi->GetIsMonteCarlo()) {
    PostData(11, fXiMCList);
  }
  if (fAntiXi->GetIsMonteCarlo()) {
    PostData(12, fAntiXiMCList);
  }
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
