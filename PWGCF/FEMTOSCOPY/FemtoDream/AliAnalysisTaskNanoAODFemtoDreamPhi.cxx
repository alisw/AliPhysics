#include "AliAnalysisTaskNanoAODFemtoDreamPhi.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliNanoAODTrack.h"
#include "AliVEvent.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskNanoAODFemtoDreamPhi)
    AliAnalysisTaskNanoAODFemtoDreamPhi::AliAnalysisTaskNanoAODFemtoDreamPhi()
    : AliAnalysisTaskSE(),
      fIsMC(false),
      fIsMCTruth(false),
      fIsmixTRUTHREAL(false),
      fIsmixTRUTHFAKE(false),
      fIsmixREC(false),
      fUseOMixing(false),
      fInvMassCutSBdown(0.0),
      fInvMassCutSBup(0.0),
      fprotonpT(0.5),
      fprotoneta(0.8),
      fkaonpT(0.15),
      fkaoneta(0.8),
      fTrigger(AliVEvent::kINT7),
      fOutput(nullptr),
      fInputEvent(nullptr),
      fEvent(nullptr),
      fTrack(nullptr),
      fPhiParticle(nullptr),
      fEventCuts(nullptr),
      fPosKaonCuts(nullptr),
      fNegKaonCuts(nullptr),
      fPhiCuts(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fSample(nullptr),
      fGTI(nullptr),
      fTrackBufferSize() {}

AliAnalysisTaskNanoAODFemtoDreamPhi::AliAnalysisTaskNanoAODFemtoDreamPhi(
    const char *name, bool isMC)
    : AliAnalysisTaskSE(name),
      fIsMC(isMC),
      fIsMCTruth(false),
      fIsmixTRUTHREAL(false),
      fIsmixTRUTHFAKE(false),
      fIsmixREC(false),
      fUseOMixing(false),
      fInvMassCutSBdown(0.0),
      fInvMassCutSBup(0.0),
      fprotonpT(0.5),
      fprotoneta(0.8),
      fkaonpT(0.15),
      fkaoneta(0.8),
      fTrigger(AliVEvent::kINT7),
      fOutput(nullptr),
      fInputEvent(nullptr),
      fEvent(nullptr),
      fTrack(nullptr),
      fPhiParticle(nullptr),
      fEventCuts(nullptr),
      fPosKaonCuts(nullptr),
      fNegKaonCuts(nullptr),
      fPhiCuts(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fSample(nullptr),
      fGTI(nullptr),
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
      new AliFemtoDreamPairCleaner(13, 1, fConfig->GetMinimalBookingME());

  fOutput->Add(fPairCleaner->GetHistList());

  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  if (fConfig->GetUsePhiSpinning()) {
    fSample = new AliFemtoDreamControlSample(fConfig);
    fOutput->Add(fSample->GetHistList());
    fOutput->Add(fSample->GetQAList());
  }

  if (fConfig->GetUseEventMixing()) {
    fOutput->Add(fPartColl->GetHistList());
    fOutput->Add(fPartColl->GetQAList());
  }

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
  static std::vector<AliFemtoDreamBasePart> ProtonTRUE;
  ProtonTRUE.clear();
  static std::vector<AliFemtoDreamBasePart> AProtonTRUE;
  AProtonTRUE.clear();
  static std::vector<AliFemtoDreamBasePart> ParticlesTRUE;
  ParticlesTRUE.clear();
  static std::vector<AliFemtoDreamBasePart> AntiParticlesTRUE;
  AntiParticlesTRUE.clear();
  static std::vector<AliFemtoDreamBasePart> PhiTRUEinvmass;
  PhiTRUEinvmass.clear();
  static std::vector<AliFemtoDreamBasePart> PhiFAKEinvmass;
  PhiFAKEinvmass.clear();
  static std::vector<AliFemtoDreamBasePart> Kaonsinvmass;
  Kaonsinvmass.clear();
  static std::vector<AliFemtoDreamBasePart> AntiKaonsinvmass;
  AntiKaonsinvmass.clear();

  static float massKaon =
      TDatabasePDG::Instance()->GetParticle(fPosKaonCuts->GetPDGCode())->Mass();

  for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack *>(fInputEvent->GetTrack(iTrack));

    if (!track) continue;
    fTrack->SetTrack(track, fInputEvent);
    // find mothers of MC Kaons (if phi->stop loop, else loop until arriving to
    // g,q,p) and set MotherPDG
    if (fIsMC) {
      TClonesArray *mcarray = dynamic_cast<TClonesArray *>(
          fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
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
        int a = 0;
        while (motherID != -1) {
          lastMother = motherID;
          mcMother = (AliAODMCParticle *)mcarray->At(motherID);
          motherID = mcMother->GetMother();
          a++;
          if (mcMother->GetPdgCode() == 333) {
            //              cout<<a<< " motherid: "<<lastMother<<" Pdg:
            //              "<<mcMother->GetPdgCode()<<endl;
            break;
          }
        }
        if ((lastMother != -1)) {
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

  if (fIsMC && fIsMCTruth) {
    TClonesArray *fArrayMCAOD = dynamic_cast<TClonesArray *>(
        fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    int noPart = fArrayMCAOD->GetEntriesFast();
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
      int mcpdg = mcPart->GetPdgCode();

      if (mcpdg == 2212) {
        double pt = mcPart->Pt();
        double eta = mcPart->Eta();
        if ((pt < 4.05 && pt > fprotonpT) && (eta > -fprotoneta && eta < fprotoneta)) {
          part.SetMCParticleRePart(mcPart);
          ProtonTRUE.push_back(part);
        }
      }

      if (mcpdg == -2212) {
        double pt = mcPart->Pt();
        double eta = mcPart->Eta();
        if ((pt < 4.05 && pt > fprotonpT) && (eta > -fprotoneta && eta < fprotoneta)) {
          part.SetMCParticleRePart(mcPart);
          AProtonTRUE.push_back(part);
        }
      }

      if (mcpdg == 321) {
        double pt = mcPart->Pt();
        double eta = mcPart->Eta();
        if ((pt < 999 && pt > fkaonpT) && (eta > -fkaoneta && eta < fkaoneta)) {
          part.SetMCParticleRePart(mcPart);
          part.SetPDGCode(mcPart->GetPdgCode());

          int motherID = mcPart->GetMother();

          if (motherID <= noPart) {
            AliAODMCParticle *mcMother = nullptr;

            mcMother = (AliAODMCParticle *)fArrayMCAOD->At(motherID);

            if (mcMother) {
              if (mcMother->GetLabel() < 0) {
                continue;
              }
              part.SetMotherID(motherID);
              part.SetMotherPDG(mcMother->GetPdgCode());

            } else {
              part.SetMotherID(0);
              part.SetMotherPDG(0);
            }

            ParticlesTRUE.push_back(part);
          }
        }
      }

      if (mcpdg == -321) {
        double pt = mcPart->Pt();
        double eta = mcPart->Eta();
        if ((pt < 999 && pt > fkaonpT) && (eta > -fkaoneta && eta < fkaoneta)) {
          part.SetMCParticleRePart(mcPart);
          part.SetPDGCode(mcPart->GetPdgCode());
          int motherID = mcPart->GetMother();

          if (motherID <= noPart) {
            AliAODMCParticle *mcMother = nullptr;

            mcMother = (AliAODMCParticle *)fArrayMCAOD->At(motherID);

            if (mcMother) {
              if (mcMother->GetLabel() < 0) {
                continue;
              }
              part.SetMotherID(motherID);
              part.SetMotherPDG(mcMother->GetPdgCode());

            } else {
              part.SetMotherID(0);
              part.SetMotherPDG(0);
            }

            AntiParticlesTRUE.push_back(part);
          }
        }
      }
    }

    for (const auto &posDaughter : ParticlesTRUE) {
      int mcpdgm1 = posDaughter.GetMotherPDG();
      int motherIDKp = posDaughter.GetMotherID();

      for (const auto &negDaughter : AntiParticlesTRUE) {
        int mcpdgm2 = negDaughter.GetMotherPDG();
        int motherIDKm = negDaughter.GetMotherID();

        float posP[3], negP[3];
        posDaughter.GetMomentum().GetXYZ(posP);
        negDaughter.GetMomentum().GetXYZ(negP);
        TLorentzVector trackPos, trackNeg;
        float Kaonmass = TDatabasePDG::Instance()->GetParticle(321)->Mass();
        trackPos.SetXYZM(posP[0], posP[1], posP[2], Kaonmass);
        trackNeg.SetXYZM(negP[0], negP[1], negP[2], Kaonmass);
        TLorentzVector trackSum = trackPos + trackNeg;
        if ((trackSum.M() < fInvMassCutSBup) &&
            (trackSum.M() > fInvMassCutSBdown)) {
          // cout << trackSum.M() << endl;
          if ((motherIDKp == motherIDKm)) {
            if ((mcpdgm1 == 333) && (mcpdgm2 == 333)) {
                if(fIsmixTRUTHREAL){
              AliAODMCParticle *mcMother = nullptr;
              mcMother = (AliAODMCParticle *)fArrayMCAOD->At(motherIDKp);
              part.SetMCParticleRePart(mcMother);
              PhiTRUEinvmass.push_back(part);
              //            cout << "invmassphireal" << endl;
                }
            }
          } else {
              if(fIsmixTRUTHFAKE){
            part.SetInvMass(trackSum.M());
            part.SetMomentum(0, trackSum.Px(), trackSum.Py(), trackSum.Pz());
            part.SetEta(trackSum.Eta());
            part.SetPhi(trackSum.Phi());
            part.SetTheta(trackSum.Theta());
            part.SetPt(posDaughter.GetPt() + posDaughter.GetPt());

            PhiFAKEinvmass.push_back(part);
            //            cout << "invmassphifake" << endl;
              }
          }
        }
      }
    }
  }

  fPhiParticle->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  for (const auto &posK : Particles) {
    for (const auto &negK : AntiParticles) {
      fPhiParticle->Setv0(posK, negK, fInputEvent, false, false, true);
      fPhiParticle->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);

      if (fPhiCuts->isSelected(fPhiParticle)) {
        fPhiParticle->SetCPA(
            gRandom->Uniform());  // cpacode needed for CleanDecay v0;
        V0Particles.push_back(*fPhiParticle);
        Kaonsinvmass.push_back(posK);
        AntiKaonsinvmass.push_back(negK);

      }
    }
  }

  fPairCleaner->CleanTrackAndDecay(&Protons, &AntiProtons, 0);
  fPairCleaner->CleanTrackAndDecay(&Protons, &V0Particles, 1);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &V0Particles, 2);

  fPairCleaner->CleanTrackAndDecay(&Particles, &AntiParticles, 3);
  fPairCleaner->CleanTrackAndDecay(&Protons, &Particles, 4);
  fPairCleaner->CleanTrackAndDecay(&Protons, &AntiParticles, 5);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiParticles, 6);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &Particles, 7);

  fPairCleaner->CleanTrackAndDecay(&Kaonsinvmass, &AntiKaonsinvmass, 8);
  fPairCleaner->CleanTrackAndDecay(&Protons, &Kaonsinvmass, 9);
  fPairCleaner->CleanTrackAndDecay(&Protons, &AntiKaonsinvmass, 10);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &AntiKaonsinvmass, 11);
  fPairCleaner->CleanTrackAndDecay(&AntiProtons, &Kaonsinvmass, 12);

  fPairCleaner->CleanDecay(&V0Particles, 0);

  fPairCleaner->ResetArray();

  fPairCleaner->StoreParticle(Protons);         // 0
  fPairCleaner->StoreParticle(AntiProtons);     // 1
  fPairCleaner->StoreParticle(V0Particles);     // 2
  fPairCleaner->StoreParticle(ProtonTRUE);      // 3
  fPairCleaner->StoreParticle(AProtonTRUE);     // 4
  fPairCleaner->StoreParticle(PhiTRUEinvmass);  // 5
  fPairCleaner->StoreParticle(PhiFAKEinvmass);  // 6
  fPairCleaner->StoreParticle(Particles);  // 7
  fPairCleaner->StoreParticle(AntiParticles);  // 8
  fPairCleaner->StoreParticle(Kaonsinvmass);  // 9
  fPairCleaner->StoreParticle(AntiKaonsinvmass);  // 10

  if (fPairCleaner->GetCounter() > 0) {
    if (fConfig->GetUseEventMixing()) {
      if (fUseOMixing) {
          std::vector<std::vector<AliFemtoDreamBasePart>> &Particlevec =
              fPairCleaner->GetCleanParticles();
          int size = Particlevec.size();
          if (size == 7) {
          if (fIsMC) {
            if (fIsMCTruth) {
              if (fIsmixTRUTHREAL) {
                if (((Particlevec.at(5)).size() > 0 &&
                     (((Particlevec.at(3)).size() > 0) ||
                      ((Particlevec.at(4)).size() > 0)))) {
                    (Particlevec.at(6)).clear();
                    (Particlevec.at(0)).clear();
                    (Particlevec.at(1)).clear();
                    (Particlevec.at(2)).clear();
                  fPartColl->SetEvent(
                      Particlevec, fEvent->GetZVertex(),
                      fEvent->GetRefMult08(), fEvent->GetV0MCentrality());
                }
              }
              if (fIsmixTRUTHFAKE) {
                if (((Particlevec.at(6)).size() > 0 &&
                     (((Particlevec.at(3)).size() > 0) ||
                      ((Particlevec.at(4)).size() > 0)))) {
                    (Particlevec.at(5)).clear();
                    (Particlevec.at(0)).clear();
                    (Particlevec.at(1)).clear();
                    (Particlevec.at(2)).clear();
                  fPartColl->SetEvent(
                      Particlevec, fEvent->GetZVertex(),
                      fEvent->GetRefMult08(), fEvent->GetV0MCentrality());
                }
              }
            }
            if (fIsmixREC) {
              if (((Particlevec.at(2)).size() > 0 &&
                   (((Particlevec.at(0)).size() > 0) ||
                    ((Particlevec.at(1)).size() > 0)))) {
                  (Particlevec.at(3)).clear();
                  (Particlevec.at(4)).clear();
                  (Particlevec.at(5)).clear();
                  (Particlevec.at(6)).clear();
                fPartColl->SetEvent(
                    Particlevec, fEvent->GetZVertex(),
                    fEvent->GetRefMult08(), fEvent->GetV0MCentrality());
              }
            }

        } else {
          if (((Particlevec.at(2)).size() > 0 &&
               (((Particlevec.at(0)).size() > 0) ||
                ((Particlevec.at(1)).size() > 0)))) {
              (Particlevec.at(3)).clear();
              (Particlevec.at(4)).clear();
              (Particlevec.at(5)).clear();
              (Particlevec.at(6)).clear();
            fPartColl->SetEvent(Particlevec,
                                fEvent->GetZVertex(), fEvent->GetRefMult08(),
                                fEvent->GetV0MCentrality());
          }
        }
          }

      } else {
        fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                            fEvent->GetZVertex(), fEvent->GetRefMult08(),
                            fEvent->GetV0MCentrality());
      }
    }

    if (fConfig->GetUsePhiSpinning()) {
      fSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
    }
  }

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
