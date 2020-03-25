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
      fUseOMixing(false),
      fTrigger(AliVEvent::kINT7),
      fOutput(),
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

AliAnalysisTaskFemtoDreamPhi::AliAnalysisTaskFemtoDreamPhi(const char *name,
                                                           bool isMC)
    : AliAnalysisTaskSE(name),
      fIsMC(isMC),
      fUseOMixing(false),
      fTrigger(AliVEvent::kINT7),
      fOutput(),
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

AliAnalysisTaskFemtoDreamPhi::~AliAnalysisTaskFemtoDreamPhi() {}

void AliAnalysisTaskFemtoDreamPhi::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetName("Output");
  fOutput->SetOwner();

  fEvent = new AliFemtoDreamEvent(false, true, fTrigger);
  // fEvent->SetCalcSpherocity(true);

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
  static std::vector<AliFemtoDreamBasePart> Protons;
  Protons.clear();
  static std::vector<AliFemtoDreamBasePart> AntiProtons;
  AntiProtons.clear();
  static std::vector<AliFemtoDreamBasePart> PhiTRUE;
  PhiTRUE.clear();
  static std::vector<AliFemtoDreamBasePart> ProtonTRUE;
  ProtonTRUE.clear();
  static std::vector<AliFemtoDreamBasePart> AProtonTRUE;
  AProtonTRUE.clear();
  static std::vector<AliFemtoDreamBasePart> PhiALL;
  PhiALL.clear();

  static float massKaon =
      TDatabasePDG::Instance()->GetParticle(fPosKaonCuts->GetPDGCode())->Mass();

  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(Event->GetTrack(iTrack));

    if (!track) continue;
    fTrack->SetTrack(track);
    // find mothers of MC Kaons (if phi->stop loop, else loop until arriving to
    // g,q,p) and set MotherPDG
    if (fIsMC) {
      TClonesArray *mcarray = dynamic_cast<TClonesArray *>(
          Event->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcarray) {
        AliError("SPTrack: MC Array not found");
      }
      if (fTrack->GetID() >= 0) {
        AliAODMCParticle *mcPart =
            (AliAODMCParticle *)mcarray->At(fTrack->GetID());
        if (!(mcPart)) {
          break;
        }
        int motherID = mcPart->GetMother();
        int lastMother = motherID;
        AliAODMCParticle *mcMother = nullptr;
        while (motherID != -1) {
          lastMother = motherID;
          mcMother = (AliAODMCParticle *)mcarray->At(motherID);
          motherID = mcMother->GetMother();
          if (mcMother->GetPdgCode() == 333) {
            break;
          }
        }
        if (lastMother != -1) {
          mcMother = (AliAODMCParticle *)mcarray->At(lastMother);
        }
        if (mcMother) {
          fTrack->SetMotherPDG(mcMother->GetPdgCode());
          //          std::cout<<"Track Mother: "<<fTrack->GetMotherPDG()<<endl;
          fTrack->SetMotherID(lastMother);
        }
      } else {
        break;  // if we don't have MC Information, don't use that track
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


  if (fIsMC) {
    TClonesArray *fArrayMCAOD = dynamic_cast<TClonesArray *>(
        Event->FindListObject(AliAODMCParticle::StdBranchName()));
    int noPart = fArrayMCAOD->GetEntriesFast();
    int mcpdg;
    AliFemtoDreamBasePart part;
    AliFemtoDreamBasePart part2;
    for (int iPart = 1; iPart < noPart; iPart++) {
      AliAODMCParticle *mcPart = (AliAODMCParticle *)fArrayMCAOD->At(iPart);
      if (!(mcPart)) {
        std::cout << "NO MC particle" << std::endl;
        continue;
      }
      if (mcPart->GetLabel() < 0) {
        continue;
      }
      mcpdg = mcPart->GetPdgCode();
      if (mcpdg == 333) {
        int firstdaughter = mcPart->GetDaughterFirst();
        if (firstdaughter <= noPart) {
          AliAODMCParticle *mcDaughter =
              (AliAODMCParticle *)fArrayMCAOD->At(firstdaughter);

          if (mcDaughter) {
            int dpdg = mcDaughter->GetPdgCode();
            double dpt = mcDaughter->Pt();
            double deta = mcDaughter->Eta();
            if (std::abs(dpdg) == 321) {
              if ((dpt < 999 && dpt > 0.15) && (deta > -0.8 && deta < 0.8)) {
                part.SetMCParticleRePart(mcPart);
                PhiTRUE.push_back(part);
              }
            }
          }
        }
      }

      if (mcpdg == 333) {
        part.SetMCParticleRePart(mcPart);
        PhiALL.push_back(part);
      }


      if (mcpdg == 2212) {
        double pt = mcPart->Pt();
        double eta = mcPart->Eta();

        if ((pt < 4.05 && pt > 0.5) && (eta > -0.8 && eta < 0.8)) {
          part.SetMCParticleRePart(mcPart);
          ProtonTRUE.push_back(part);
        }
      }

      if (mcpdg == -2212) {
        double pt = mcPart->Pt();
        double eta = mcPart->Eta();
        if ((pt < 4.05 && pt > 0.5) && (eta > -0.8 && eta < 0.8)) {
          part.SetMCParticleRePart(mcPart);
          AProtonTRUE.push_back(part);
        }
      }

    }
  }

  fPhiParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (const auto &posK : Particles) {
    for (const auto &negK : AntiParticles) {
      fPhiParticle->Setv0(posK, negK, Event, false, false, true);
      fPhiParticle->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);
      //      std::cout<<"ID mother kp: "<<posK.GetMotherID()<<endl;
      //      std::cout<<"ID mother km: "<<negK.GetMotherID()<<endl;
      //      std::cout<<"PDG kp: "<<posK.GetMCPDGCode()<<endl;
      //      std::cout<<"PDG km: "<<negK.GetMCPDGCode()<<endl;
      //      std::cout<<"PDG phi: "<<fPhiParticle->GetMCPDGCode()<<endl;
      if (fPhiCuts->isSelected(fPhiParticle)) {
        fPhiParticle->SetCPA(
            gRandom->Uniform());  // cpacode needed for CleanDecay v0;
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
  fPairCleaner->StoreParticle(ProtonTRUE);
  fPairCleaner->StoreParticle(AProtonTRUE);
  fPairCleaner->StoreParticle(PhiTRUE);
  fPairCleaner->StoreParticle(PhiALL);



  if (fConfig->GetUseEventMixing()) {
    if (fUseOMixing) {
      std::vector<std::vector<AliFemtoDreamBasePart>> &Particles =
          fPairCleaner->GetCleanParticles();
      int size = Particles.size();
      if (size == 7) {
        if ((Particles.at(2)).size() > 0) {
          if (((Particles.at(0)).size() > 0) ||
              ((Particles.at(1)).size() > 0)) {
              fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                                  fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

          }
        }

        if (fIsMC) {
          if ((Particles.at(5)).size() > 0 ||
              ((Particles.at(6)).size() > 0)) {
            if (((Particles.at(3)).size() > 0) ||
                ((Particles.at(4)).size() > 0)) {
                fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                                    fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

            }
          }
        }
      }
    } else {
        fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                            fEvent->GetRefMult08(), fEvent->GetV0MCentrality());

    }
  }



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
