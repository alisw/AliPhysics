#include "AliAnalysisTaskFemtoDreamPhi.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"

ClassImp(AliAnalysisTaskFemtoDreamPhi)
    AliAnalysisTaskFemtoDreamPhi::AliAnalysisTaskFemtoDreamPhi()
    : AliAnalysisTaskSE(),
      fIsMC(false),
      fOutput(),
      fEvent(),
      fTrack(),
      fPhiParticle(),
      fEventCuts(),
      fPosKaonCuts(),
      fNegKaonCuts(),
      fPhiCuts(),
      fConfig(),
      fPairCleaner(),
      fPartColl(),
      fGTI(),
      fTrackBufferSize() {}

AliAnalysisTaskFemtoDreamPhi::AliAnalysisTaskFemtoDreamPhi(const char *name,
                                                           bool isMC)
    : AliAnalysisTaskSE(name),
      fIsMC(isMC),
      fOutput(),
      fEvent(),
      fTrack(),
      fPhiParticle(),
      fEventCuts(),
      fPosKaonCuts(),
      fNegKaonCuts(),
      fPhiCuts(),
      fConfig(),
      fPairCleaner(),
      fPartColl(),
      fGTI(),
      fTrackBufferSize(2000) {
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskFemtoDreamPhi::~AliAnalysisTaskFemtoDreamPhi() {}

void AliAnalysisTaskFemtoDreamPhi::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetName("Output");
  fOutput->SetOwner();

  fEvent = new AliFemtoDreamEvent(false, true, AliVEvent::kINT7);
  fOutput->Add(fEvent->GetEvtCutList());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fPhiParticle = new AliFemtoDreamv0();
  fPhiParticle->SetPDGCode(fPhiCuts->GetPDGv0());
  fPhiParticle->SetUseMCInfo(fIsMC);
  fPhiParticle->SetPDGDaughterPos(
      fPhiCuts->GetPDGPosDaug());  // order +sign doesnt play a role
  fPhiParticle->GetPosDaughter()->SetUseMCInfo(fIsMC);
  fPhiParticle->SetPDGDaughterNeg(
      fPhiCuts->GetPDGNegDaug());  // only used for MC Matching
  fPhiParticle->GetNegDaughter()->SetUseMCInfo(fIsMC);

  fGTI = new AliAODTrack *[fTrackBufferSize];

  if (!fEventCuts) {
    AliFatal("Event Cuts not set!");
  }
  fEventCuts->InitQA();
  fOutput->Add(fEventCuts->GetHistList());

  if (!fPosKaonCuts) {
    AliFatal("Track Cuts for Particle One not set!");
  }
  fPosKaonCuts->Init();
  fPosKaonCuts->SetName("Particle1");
  fOutput->Add(fPosKaonCuts->GetQAHists());

  if (fPosKaonCuts->GetIsMonteCarlo()) {
    fPosKaonCuts->SetMCName("MCParticle1");
    fOutput->Add(fPosKaonCuts->GetMCQAHists());
  }

  if (!fNegKaonCuts) {
    AliFatal("Track Cuts for Particle One not set!");
  }
  fNegKaonCuts->Init();
  fNegKaonCuts->SetName("Particle2");
  fOutput->Add(fNegKaonCuts->GetQAHists());
  if (fNegKaonCuts->GetIsMonteCarlo()) {
    fNegKaonCuts->SetMCName("MCParticle2");
    fOutput->Add(fNegKaonCuts->GetMCQAHists());
  }

  if (!fPhiCuts) {
    AliFatal("Cuts for the phi not set!");
  }
  fPhiCuts->Init();
  fPhiCuts->SetName("Phi");
  fOutput->Add(fPhiCuts->GetQAHists());
  if (fPhiCuts->GetIsMonteCarlo()) {
    fPhiCuts->SetMCName("MCPhi");
    fOutput->Add(fPhiCuts->GetMCQAHists());
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

void AliAnalysisTaskFemtoDreamPhi::UserExec(Option_t *) {
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

  static std::vector<AliFemtoDreamBasePart> Particles;
  Particles.clear();
  static std::vector<AliFemtoDreamBasePart> AntiParticles;
  AntiParticles.clear();
  static std::vector<AliFemtoDreamBasePart> V0Particles;
  V0Particles.clear();

  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));
    if (!track) continue;
    fTrack->SetTrack(track);
    if (fPosKaonCuts->isSelected(fTrack)) {
      Particles.push_back(*fTrack);
    }
    if (fNegKaonCuts->isSelected(fTrack)) {
      AntiParticles.push_back(*fTrack);
    }
  }

  static float massKaon =
      TDatabasePDG::Instance()->GetParticle(fPosKaonCuts->GetPDGCode())->Mass();
  fPhiParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (const auto &posK : Particles) {
    for (const auto &negK : AntiParticles) {
      fPhiParticle->Setv0(posK, massKaon, negK, massKaon);
      if (fPhiCuts->isSelected(fPhiParticle)) {
        V0Particles.push_back(*fPhiParticle);
      }
    }
  }

  fPairCleaner->CleanTrackAndDecay(&Particles, &AntiParticles, 0);
  fPairCleaner->CleanTrackAndDecay(&Particles, &V0Particles, 1);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &V0Particles, 2);
  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  fPairCleaner->StoreParticle(V0Particles);
  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

  PostData(1, fOutput);
}

void AliAnalysisTaskFemtoDreamPhi::ResetGlobalTrackReference() {
  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskFemtoDreamPhi::StoreGlobalTrackReference(
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
