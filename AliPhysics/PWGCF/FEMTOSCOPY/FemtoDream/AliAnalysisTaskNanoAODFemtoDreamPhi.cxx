#include "AliAnalysisTaskNanoAODFemtoDreamPhi.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisTaskNanoAODFemtoDreamPhi)
    AliAnalysisTaskNanoAODFemtoDreamPhi::AliAnalysisTaskNanoAODFemtoDreamPhi()
    : AliAnalysisTaskSE(),
      fIsMC(false),
      fTrigger(AliVEvent::kINT7),
      fOutput(),
      fInputEvent(nullptr),
      fEvent(),
      fTrack(),
      fPhiParticle(),
      fEventCuts(),
      fPosKaonCuts(),
      fNegKaonCuts(),
      fPhiCuts(),
      fTrackCutsPartProton(),
      fTrackCutsPartAntiProton(),
      fConfig(),
      fPairCleaner(),
      fPartColl(),
      fGTI(),
      fTrackBufferSize() {}

AliAnalysisTaskNanoAODFemtoDreamPhi::AliAnalysisTaskNanoAODFemtoDreamPhi(const char *name,
                                                           bool isMC)
    : AliAnalysisTaskSE(name),
      fIsMC(isMC),
      fTrigger(AliVEvent::kINT7),
      fOutput(),
      fInputEvent(nullptr),
      fEvent(),
      fTrack(),
      fPhiParticle(),
      fEventCuts(),
      fPosKaonCuts(),
      fNegKaonCuts(),
      fPhiCuts(),
      fTrackCutsPartProton(),
      fTrackCutsPartAntiProton(),
      fConfig(),
      fPairCleaner(),
      fPartColl(),
      fGTI(),
      fTrackBufferSize(2000) {
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskNanoAODFemtoDreamPhi::~AliAnalysisTaskNanoAODFemtoDreamPhi() {}

void AliAnalysisTaskNanoAODFemtoDreamPhi::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetName("Output");
  fOutput->SetOwner();

  fEvent = new AliFemtoDreamEvent(false, true, fTrigger);
  fEvent->SetCalcSpherocity(true);
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

  fGTI = new AliVTrack *[fTrackBufferSize];


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

  if (!fTrackCutsPartProton) {
    AliFatal("Track Cuts for Proton not set!");
  }
  fTrackCutsPartProton->Init();
  fTrackCutsPartProton->SetName("Proton");
  fOutput->Add(fTrackCutsPartProton->GetQAHists());
  if (fTrackCutsPartProton->GetIsMonteCarlo()) {
    fTrackCutsPartProton->SetMCName("MCProton");
    fOutput->Add(fTrackCutsPartProton->GetMCQAHists());
  }

  if (!fTrackCutsPartAntiProton) {
    AliFatal("Track Cuts for AntiProton not set!");
  }
  fTrackCutsPartAntiProton->Init();
  fTrackCutsPartAntiProton->SetName("AntiProton");
  fOutput->Add(fTrackCutsPartAntiProton->GetQAHists());
  if (fTrackCutsPartAntiProton->GetIsMonteCarlo()) {
    fTrackCutsPartAntiProton->SetMCName("MCAntiProton");
    fOutput->Add(fTrackCutsPartAntiProton->GetMCQAHists());
  }

  fPairCleaner =
      new AliFemtoDreamPairCleaner(3, 1, fConfig->GetMinimalBookingME());
  fOutput->Add(fPairCleaner->GetHistList());

  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());

  PostData(1, fOutput);
}

void AliAnalysisTaskNanoAODFemtoDreamPhi::UserExec(Option_t *) {
    AliVEvent *fInputEvent = InputEvent();

    // PREAMBLE - CHECK EVERYTHING IS THERE
    if (!fInputEvent) {
      AliError("No Input event");
      return;
    }

    // EVENT SELECTION
    fEvent->SetEvent(fInputEvent);
    if (!fEventCuts->isSelected(fEvent)) return;

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }
  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

  static std::vector<AliFemtoDreamBasePart> Particles;
  Particles.clear();
  static std::vector<AliFemtoDreamBasePart> AntiParticles;
  AntiParticles.clear();
  static std::vector<AliFemtoDreamBasePart> V0Particles;
  V0Particles.clear();
  static std::vector<AliFemtoDreamBasePart> Protons;
  Protons.clear();
  static std::vector<AliFemtoDreamBasePart> AntiProtons;
  AntiProtons.clear();

  static float massKaon =
      TDatabasePDG::Instance()->GetParticle(fPosKaonCuts->GetPDGCode())->Mass();

  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));

    if (!track) continue;
    fTrack->SetTrack(track, fInputEvent);
    //find mothers of MC Kaons (if phi->stop loop, else loop until arriving to g,q,p) and set MotherPDG
    if (fIsMC) {
        TClonesArray *mcarray = dynamic_cast<TClonesArray*> (fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
        if (!mcarray) {
          AliError("SPTrack: MC Array not found");
        }
        if (fTrack->GetID() >= 0) {
          AliAODMCParticle * mcPart = (AliAODMCParticle*) mcarray->At(fTrack->GetID());
          if (!(mcPart)) {
            break;
          }
          int motherID = mcPart->GetMother();
          int lastMother = motherID;
          AliAODMCParticle *mcMother = nullptr;
          while (motherID != -1) {
            lastMother = motherID;
            mcMother = (AliAODMCParticle *) mcarray->At(motherID);
            motherID = mcMother->GetMother();
            if (mcMother->GetPdgCode()==333) {
              break;
            }
          }
          if (lastMother != -1) {
            mcMother = (AliAODMCParticle *) mcarray->At(lastMother);
          }
          if (mcMother) {
            fTrack->SetMotherPDG(mcMother->GetPdgCode());
//          std::cout<<"Track Mother: "<<fTrack->GetMotherPDG()<<endl;
            fTrack->SetMotherID(lastMother);
         }
        }
        else {
          break;  //if we don't have MC Information, don't use that track
        }
    }
    fTrack->SetInvMass(massKaon);
    if (fPosKaonCuts->isSelected(fTrack)) {
      Particles.push_back(*fTrack);
    }
    if (fNegKaonCuts->isSelected(fTrack)) {
      AntiParticles.push_back(*fTrack);
    }
    if (fTrackCutsPartProton->isSelected(fTrack)) {
      Protons.push_back(*fTrack);
    }
    if (fTrackCutsPartAntiProton->isSelected(fTrack)) {
      AntiProtons.push_back(*fTrack);
    }
  }

  fPhiParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (const auto &posK : Particles) {
    for (const auto &negK : AntiParticles) {
      fPhiParticle->Setv0(posK, negK, fInputEvent, false, false, true);
      fPhiParticle->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);
//      std::cout<<"ID mother kp: "<<posK.GetMotherID()<<endl;
//      std::cout<<"ID mother km: "<<negK.GetMotherID()<<endl;
//      std::cout<<"PDG kp: "<<posK.GetMCPDGCode()<<endl;
//      std::cout<<"PDG km: "<<negK.GetMCPDGCode()<<endl;
//      std::cout<<"PDG phi: "<<fPhiParticle->GetMCPDGCode()<<endl;
      if (fPhiCuts->isSelected(fPhiParticle)) {
        fPhiParticle->SetCPA(gRandom->Uniform()); //cpacode needed for CleanDecay v0;
        V0Particles.push_back(*fPhiParticle);
//      std::cout<<"PDG phi cut: "<<fPhiParticle->GetMCPDGCode()<<endl;
      }
    }
  }

  fPairCleaner->CleanTrackAndDecay(&Protons, &AntiProtons, 0);
  fPairCleaner->CleanTrackAndDecay(&Protons, &V0Particles, 1);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &V0Particles, 2);
  fPairCleaner->CleanDecay(&V0Particles, 0);
  fPairCleaner->ResetArray();

  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(V0Particles);
  fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                      fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

  PostData(1, fOutput);
}

void AliAnalysisTaskNanoAODFemtoDreamPhi::ResetGlobalTrackReference() {
  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

void AliAnalysisTaskNanoAODFemtoDreamPhi::StoreGlobalTrackReference(
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
        if (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID])->GetFilterMap() ||
            fGTI[trackID]->GetTPCNcls()) {
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
