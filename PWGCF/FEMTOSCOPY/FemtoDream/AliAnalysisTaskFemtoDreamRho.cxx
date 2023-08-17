/*
 * AliAnalysisTaskFemtoDreamRho.cxx
 *
 *  Created on: 19 Jul 2023
 *      Author: M. Korwieser
 */

#include "AliAnalysisTaskFemtoDreamRho.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"

ClassImp(AliAnalysisTaskFemtoDreamRho)
    AliAnalysisTaskFemtoDreamRho::AliAnalysisTaskFemtoDreamRho()
    : AliAnalysisTaskSE(),
      fTrigger(AliVEvent::kINT7),
      fIsMC(false),
      fDoCleaning(false),
      fOutput(),
      fEvent(),
      fTrack(),
      fRhoParticle(),
      fEventCuts(),
      fPosPionCuts(),
      fNegPionCuts(),
      fRhoCuts(),
      fPosProtonCuts(),
      fNegProtonCuts(),
      fConfig(),
      fPairCleaner(),
      fPartColl(),
      fGTI(),
      fTrackBufferSize() {}

AliAnalysisTaskFemtoDreamRho::AliAnalysisTaskFemtoDreamRho(const char *name,
                                                           bool isMC, bool doCleaning)
    : AliAnalysisTaskSE(name),
      fTrigger(AliVEvent::kINT7),
      fIsMC(isMC),
      fDoCleaning(doCleaning),
      fOutput(),
      fEvent(),
      fTrack(),
      fRhoParticle(),
      fEventCuts(),
      fPosPionCuts(),
      fNegPionCuts(),
      fRhoCuts(),
      fPosProtonCuts(),
      fNegProtonCuts(),
      fConfig(),
      fPairCleaner(),
      fPartColl(),
      fGTI(),
      fTrackBufferSize(2000) {
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskFemtoDreamRho::~AliAnalysisTaskFemtoDreamRho() {}

void AliAnalysisTaskFemtoDreamRho::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetName("Output");
  fOutput->SetOwner();

  fEvent = new AliFemtoDreamEvent(false, true, fTrigger);
  fOutput->Add(fEvent->GetEvtCutList());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fRhoParticle = new AliFemtoDreamv0();
  fRhoParticle->SetPDGCode(fRhoCuts->GetPDGv0());
  fRhoParticle->SetUseMCInfo(fIsMC);
  fRhoParticle->SetPDGDaughterPos(
      fRhoCuts->GetPDGPosDaug());  // order +sign doesnt play a role
  fRhoParticle->GetPosDaughter()->SetUseMCInfo(fIsMC);
  fRhoParticle->SetPDGDaughterNeg(
      fRhoCuts->GetPDGNegDaug());  // only used for MC Matching
  fRhoParticle->GetNegDaughter()->SetUseMCInfo(fIsMC);

  fGTI = new AliAODTrack *[fTrackBufferSize];

  if (!fEventCuts) {
    AliFatal("Event Cuts not set!");
  }
  fEventCuts->InitQA();
  fOutput->Add(fEventCuts->GetHistList());

  if (!fPosPionCuts) {
    AliFatal("Track Cuts for positive pion not set!");
  }
  fPosPionCuts->Init();
  fPosPionCuts->SetName("PosPion");
  fOutput->Add(fPosPionCuts->GetQAHists());

  if (fPosPionCuts->GetIsMonteCarlo()) {
    fPosPionCuts->SetMCName("MCPosPion");
    fOutput->Add(fPosPionCuts->GetMCQAHists());
  }

  if (!fNegPionCuts) {
    AliFatal("Track Cuts for negative pion not set!");
  }
  fNegPionCuts->Init();
  fNegPionCuts->SetName("NegPion");
  fOutput->Add(fNegPionCuts->GetQAHists());
  if (fNegPionCuts->GetIsMonteCarlo()) {
    fNegPionCuts->SetMCName("MCNegPion");
    fOutput->Add(fNegPionCuts->GetMCQAHists());
  }

  if (!fRhoCuts) {
    AliFatal("Cuts for the Rho not set!");
  }
  fRhoCuts->Init();
  fRhoCuts->SetName("Rho");
  fOutput->Add(fRhoCuts->GetQAHists());
  if (fRhoCuts->GetIsMonteCarlo()) {
    fRhoCuts->SetMCName("MCRho");
    fOutput->Add(fRhoCuts->GetMCQAHists());
  }

  if (!fPosProtonCuts) {
    AliFatal("Track Cuts for Proton not set!");
  }
  fPosProtonCuts->Init();
  fPosProtonCuts->SetName("Proton");
  fOutput->Add(fPosProtonCuts->GetQAHists());
  if (fPosProtonCuts->GetIsMonteCarlo()) {
    fPosProtonCuts->SetMCName("MCProton");
    fOutput->Add(fPosProtonCuts->GetMCQAHists());
  }

  if (!fNegProtonCuts) {
    AliFatal("Track Cuts for AntiProton not set!");
  }
  fNegProtonCuts->Init();
  fNegProtonCuts->SetName("AntiProton");
  fOutput->Add(fNegProtonCuts->GetQAHists());
  if (fNegProtonCuts->GetIsMonteCarlo()) {
    fNegProtonCuts->SetMCName("MCAntiProton");
    fOutput->Add(fNegProtonCuts->GetMCQAHists());
  }

  fPairCleaner =
      new AliFemtoDreamPairCleaner(3, 0, fConfig->GetMinimalBookingME());
  fOutput->Add(fPairCleaner->GetHistList());

  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());

  PostData(1, fOutput);
}

void AliAnalysisTaskFemtoDreamRho::UserExec(Option_t *) {
  AliAODEvent *Event = static_cast<AliAODEvent *>(fInputEvent);
  if (!Event) {
    AliWarning("No Input Event");
  }

  fEvent->SetEvent(Event);
  if (!fEventCuts->isSelected(fEvent)) return;

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
    if (!track) continue;
    StoreGlobalTrackReference(track);
  }
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  // First we want to combine all charged pions with each other in the SE in order to find Rhos
  static std::vector<AliFemtoDreamBasePart> Particles; //pi+ candidates
  Particles.clear();
  static std::vector<AliFemtoDreamBasePart> AntiParticles; //pi- candidates
  AntiParticles.clear();
  static std::vector<AliFemtoDreamBasePart> V0Particles; //Rho candidates
  V0Particles.clear();
  static std::vector<AliFemtoDreamBasePart> Protons; //proton candidates
  Protons.clear();
  static std::vector<AliFemtoDreamBasePart> AntiProtons; //anti-proton candidates
  AntiProtons.clear();

  static float massChargedPion =
      TDatabasePDG::Instance()->GetParticle(fPosPionCuts->GetPDGCode())->Mass(); // as usual to minimize uncert.
  fRhoParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  // Loop to identify all charged pions & protons
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
    if (!track) continue;
    fTrack->SetTrack(track);
    if (fPosPionCuts->isSelected(fTrack)) {
      fTrack->SetInvMass(massChargedPion);//Since we combine these later we set the inv. mass to the PDG value
      Particles.push_back(*fTrack);
    }
    if (fNegPionCuts->isSelected(fTrack)) {
      fTrack->SetInvMass(massChargedPion);//Since we combine these later we set the inv. mass to the PDG value
      AntiParticles.push_back(*fTrack);
    }
    if (fPosProtonCuts->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
    }
    if (fNegProtonCuts->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
    }
  }

  // Construct the v0 for the Rho decay
  for (const auto &posPion : Particles) { // Build charged pion pairs!
    for (const auto &negPion : AntiParticles) {
      fRhoParticle->Setv0(posPion, negPion, Event);
      if (fRhoCuts->isSelected(fRhoParticle)) { // Check for proper Rho candidates
        V0Particles.push_back(*fRhoParticle); 
      }
    }
  }
 
  // For now no default pair cleaning but may change to CleanTrackandTrack of same charge pion and proton  
  if (fDoCleaning) {
    fPairCleaner->CleanTrackAndDecay(&Particles, &Protons, 0);
    fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiProtons, 1);
    fPairCleaner->CleanDecay(&V0Particles, 0);
  }

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  fPairCleaner->StoreParticle(V0Particles);
  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);

  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

  PostData(1, fOutput);
}

void AliAnalysisTaskFemtoDreamRho::ResetGlobalTrackReference() {
  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskFemtoDreamRho::StoreGlobalTrackReference(
    AliAODTrack *track) {
  // for documentation see AliFemtoDreamAnalysis

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
    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}