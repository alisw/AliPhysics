#include "AliSigma0PhotonMotherCuts.h"
#include <iostream>
#include "AliMultSelection.h"

ClassImp(AliSigma0PhotonMotherCuts)

    //____________________________________________________________________________________________________
    AliSigma0PhotonMotherCuts::AliSigma0PhotonMotherCuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fIsMC(false),
      fIsLightweight(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fSigma(),
      fSidebandUp(),
      fSidebandDown(),
      fLambdaMixed(),
      fPhotonMixed(),
      fLambdaMixedBinned(),
      fPhotonMixedBinned(),
      fTreeVariables(),
      fLambdaCuts(nullptr),
      fPhotonCuts(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fMixingDepth(10),
      fPDG(0),
      fPDGDaughter1(0),
      fPDGDaughter2(0),
      fMassSigma(0),
      fSigmaMassCut(0.01),
      fSidebandCutUp(0.05),
      fSidebandCutDown(0.01),
      fPhotonPtMin(0),
      fPhotonPtMax(999.f),
      fRapidityMax(0.5),
      fArmenterosCut(false),
      fArmenterosQtLow(0.f),
      fArmenterosQtUp(0.f),
      fArmenterosAlphaLow(0.f),
      fArmenterosAlphaUp(0.f),
      fMCHighMultThreshold(5.f),
      fHistCutBooking(nullptr),
      fHistNSigma(nullptr),
      fHistMassCutPt(nullptr),
      fHistInvMass(nullptr),
      fHistInvMassBeforeArmenteros(nullptr),
      fHistInvMassRecPhoton(nullptr),
      fHistInvMassRecLambda(nullptr),
      fHistInvMassRec(nullptr),
      fHistInvMassPt(nullptr),
      fHistInvMassEta(nullptr),
      fHistEtaPhi(nullptr),
      fHistPtY(),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistMixedPtY(),
      fHistMixedInvMassPt(nullptr),
      fHistMixedInvMassBinnedPt(nullptr),
      fHistLambdaPtPhi(nullptr),
      fHistLambdaPtEta(nullptr),
      fHistLambdaMassPt(nullptr),
      fHistPhotonPtPhi(nullptr),
      fHistPhotonPtEta(nullptr),
      fHistPhotonMassPt(nullptr),
      fHistMCTruthPtY(nullptr),
      fHistMCTruthPtEta(nullptr),
      fHistMCTruthDaughterPtY(nullptr),
      fHistMCTruthDaughterPtEta(nullptr),
      fHistMCTruthDaughterPtYAccept(nullptr),
      fHistMCTruthDaughterPtEtaAccept(nullptr),
      fHistMCTruthPtYHighMult(nullptr),
      fHistMCTruthPtEtaHighMult(nullptr),
      fHistMCTruthDaughterPtYHighMult(nullptr),
      fHistMCTruthDaughterPtEtaHighMult(nullptr),
      fHistMCTruthDaughterPtYAcceptHighMult(nullptr),
      fHistMCTruthDaughterPtEtaAcceptHighMult(nullptr),
      fHistMCV0Pt(nullptr),
      fHistMCV0Mass(nullptr),
      fHistMCV0Mother(nullptr),
      fHistMCV0Check(nullptr),
      fHistMCV0MotherCheck(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0PhotonMotherCuts::AliSigma0PhotonMotherCuts(
    const AliSigma0PhotonMotherCuts &ref)
    : TObject(ref),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fIsMC(false),
      fIsLightweight(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fSigma(),
      fSidebandUp(),
      fSidebandDown(),
      fLambdaMixed(),
      fPhotonMixed(),
      fTreeVariables(),
      fLambdaCuts(nullptr),
      fPhotonCuts(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fMixingDepth(10),
      fPDG(0),
      fPDGDaughter1(0),
      fPDGDaughter2(0),
      fMassSigma(0),
      fSigmaMassCut(0.01),
      fSidebandCutUp(0.05),
      fSidebandCutDown(0.01),
      fPhotonPtMin(0),
      fPhotonPtMax(999.f),
      fRapidityMax(0.5),
      fArmenterosCut(false),
      fArmenterosQtLow(0.f),
      fArmenterosQtUp(0.f),
      fArmenterosAlphaLow(0.f),
      fArmenterosAlphaUp(0.f),
      fMCHighMultThreshold(5.f),
      fHistCutBooking(nullptr),
      fHistNSigma(nullptr),
      fHistMassCutPt(nullptr),
      fHistInvMass(nullptr),
      fHistInvMassBeforeArmenteros(nullptr),
      fHistInvMassRecPhoton(nullptr),
      fHistInvMassRecLambda(nullptr),
      fHistInvMassRec(nullptr),
      fHistInvMassPt(nullptr),
      fHistInvMassEta(nullptr),
      fHistEtaPhi(nullptr),
      fHistPtY(),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistMixedPtY(),
      fHistMixedInvMassPt(nullptr),
      fHistMixedInvMassBinnedPt(nullptr),
      fHistLambdaPtPhi(nullptr),
      fHistLambdaPtEta(nullptr),
      fHistLambdaMassPt(nullptr),
      fHistPhotonPtPhi(nullptr),
      fHistPhotonPtEta(nullptr),
      fHistPhotonMassPt(nullptr),
      fHistMCTruthPtY(nullptr),
      fHistMCTruthPtEta(nullptr),
      fHistMCTruthDaughterPtY(nullptr),
      fHistMCTruthDaughterPtEta(nullptr),
      fHistMCTruthDaughterPtYAccept(nullptr),
      fHistMCTruthDaughterPtEtaAccept(nullptr),
      fHistMCTruthPtYHighMult(nullptr),
      fHistMCTruthPtEtaHighMult(nullptr),
      fHistMCTruthDaughterPtYHighMult(nullptr),
      fHistMCTruthDaughterPtEtaHighMult(nullptr),
      fHistMCTruthDaughterPtYAcceptHighMult(nullptr),
      fHistMCTruthDaughterPtEtaAcceptHighMult(nullptr),
      fHistMCV0Pt(nullptr),
      fHistMCV0Mass(nullptr),
      fHistMCV0Mother(nullptr),
      fHistMCV0Check(nullptr),
      fHistMCV0MotherCheck(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0PhotonMotherCuts &AliSigma0PhotonMotherCuts::operator=(
    const AliSigma0PhotonMotherCuts &ref) {
  // Assignment operator
  if (this == &ref) return *this;

  return (*this);
}

//____________________________________________________________________________________________________
AliSigma0PhotonMotherCuts *AliSigma0PhotonMotherCuts::DefaultCuts() {
  AliSigma0PhotonMotherCuts *photonMotherCuts = new AliSigma0PhotonMotherCuts();
  photonMotherCuts->SetPhotonMaxPt(2);
  return photonMotherCuts;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SelectPhotonMother(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    const std::vector<AliSigma0ParticleV0> &photonCandidates,
    const std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  fInputEvent = inputEvent;
  fMCEvent = mcEvent;

  if (fIsMC) {
    ProcessMC();
    fV0Reader =
        (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
            fV0ReaderName.Data());
  }
  if (!fIsLightweight) SingleV0QA(photonCandidates, lambdaCandidates);

  // Particle pairing
  SigmaToLambdaGamma(photonCandidates, lambdaCandidates);
  SigmaToLambdaGammaMixedEvent(photonCandidates, lambdaCandidates);
  SigmaToLambdaGammaMixedEventBinned(photonCandidates, lambdaCandidates);
  FillEventBuffer(photonCandidates, lambdaCandidates);
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SingleV0QA(
    const std::vector<AliSigma0ParticleV0> &photonCandidates,
    const std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  // Photon QA - keep track of kinematics
  for (const auto &photon : photonCandidates) {
    fHistPhotonPtPhi->Fill(photon.GetPt(), photon.GetPhi());
    fHistPhotonPtEta->Fill(photon.GetPt(), photon.GetEta());
    fHistPhotonMassPt->Fill(photon.GetPt(), photon.GetMass());
  }

  for (const auto &lambda : lambdaCandidates) {
    // Lambda QA - keep track of kinematics
    fHistLambdaPtPhi->Fill(lambda.GetPt(), lambda.GetPhi());
    fHistLambdaPtEta->Fill(lambda.GetPt(), lambda.GetEta());
    fHistLambdaMassPt->Fill(lambda.GetPt(), lambda.GetMass());
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGamma(
    const std::vector<AliSigma0ParticleV0> &photonCandidates,
    const std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  fSigma.clear();
  fSidebandUp.clear();
  fSidebandDown.clear();

  // Mulitplicity estimator: V0M
  Float_t lPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection =
      (AliMultSelection *)fInputEvent->FindListObject("MultSelection");
  if (MultSelection) {
    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }
  int nSigma = 0;

  // SAME EVENT
  for (const auto &photon : photonCandidates) {
    if (photon.GetPt() > fPhotonPtMax || photon.GetPt() < fPhotonPtMin)
      continue;

    for (const auto &lambda : lambdaCandidates) {
      AliSigma0ParticlePhotonMother sigma(lambda, photon, fInputEvent);
      const float invMass = sigma.GetMass();
      const float armAlpha = sigma.GetArmenterosAlpha();
      const float armQt = sigma.GetArmenterosQt();
      const float pT = sigma.GetPt();
      if (!fIsLightweight) {
        fHistArmenterosBefore->Fill(armAlpha, armQt);
        fHistInvMassBeforeArmenteros->Fill(invMass);
      }
      // Armenteros cut
      if (fArmenterosCut) {
        if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
        if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
          continue;
      }
      const float rap = sigma.GetRapidity();

      // Now write out the stuff to the Femto containers
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        fSigma.push_back(sigma);
        if (!fIsLightweight) {
          fHistMassCutPt->Fill(pT);
          fHistEtaPhi->Fill(sigma.GetEta(), sigma.GetPhi());
        }
        ++nSigma;
      }
      if (invMass < fMassSigma + fSidebandCutUp &&
          invMass > fMassSigma + fSidebandCutDown) {
        fSidebandUp.push_back(sigma);
      }
      if (invMass > fMassSigma - fSidebandCutUp &&
          invMass < fMassSigma - fSidebandCutDown) {
        fSidebandDown.push_back(sigma);
      }

      fHistInvMass->Fill(invMass);

      if (std::abs(rap) > fRapidityMax) continue;
      const int rapBin = GetRapidityBin(rap);
      if (!fIsLightweight) {
        fHistArmenterosAfter->Fill(armAlpha, armQt);
        fHistInvMassRecPhoton->Fill(pT, sigma.GetRecMassPhoton());
        fHistInvMassRecLambda->Fill(pT, sigma.GetRecMassLambda());
        fHistInvMassRec->Fill(pT, sigma.GetRecMass());
        if (rapBin > -1) fHistPtY[rapBin]->Fill(pT, invMass);
        fHistInvMassEta->Fill(sigma.GetEta(), invMass);
      }
      fHistInvMassPt->Fill(pT, invMass);

      if (fIsMC) {
        int pdgLambdaMother = 0;
        int pdgPhotonMother = 0;
        const int label =
            sigma.MatchToMC(fMCEvent, fPDG, {{fPDGDaughter1, fPDGDaughter2}},
                            pdgLambdaMother, pdgPhotonMother);
        if (label > 0) {
          fHistMCV0Pt->Fill(sigma.GetPt());
          fHistMCV0Mass->Fill(invMass);
        }

        // let's where the other particle comes from if one of them stems from
        // a Sigma0
        if (std::abs(pdgLambdaMother) == 3212 &&
            std::abs(pdgLambdaMother) != 3212) {
          fHistMCV0Mother->Fill(invMass, std::abs(pdgPhotonMother));
        }
        if (std::abs(pdgLambdaMother) == 3212 &&
            std::abs(pdgLambdaMother) != 3212) {
          fHistMCV0Mother->Fill(invMass, std::abs(pdgLambdaMother));
        }
        fHistMCV0MotherCheck->Fill(std::abs(pdgLambdaMother),
                                   std::abs(pdgPhotonMother));

        const int labV0 = photon.GetMCLabelV0();
        const int labPhoton = lambda.GetMCLabelV0();
        if (labV0 < 0 || labPhoton < 0) continue;

        AliMCParticle *partV0 =
            static_cast<AliMCParticle *>(fMCEvent->GetTrack(labV0));
        AliMCParticle *partPhoton =
            static_cast<AliMCParticle *>(fMCEvent->GetTrack(labPhoton));
        if (!partV0 || !partPhoton) continue;

        fHistMCV0Check->Fill(std::abs(partV0->PdgCode()),
                             std::abs(partPhoton->PdgCode()));
      }
    }
  }
  if (!fIsLightweight) fHistNSigma->Fill(nSigma);
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGammaMixedEvent(
    const std::vector<AliSigma0ParticleV0> &photonCandidates,
    const std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fLambdaMixed) {
    for (const auto &Lambda : LambdaContainer) {
      if (!Lambda.GetIsUse()) continue;
      for (auto Photon : photonCandidates) {
        if (Photon.GetPt() > fPhotonPtMax || Photon.GetPt() < fPhotonPtMin)
          continue;
        AliSigma0ParticlePhotonMother sigma(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma.GetArmenterosAlpha();
        const float armQt = sigma.GetArmenterosQt();
        const float pT = sigma.GetPt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
            continue;
        }
        const float invMass = sigma.GetMass();
        const float rap = sigma.GetRapidity();
        if (std::abs(rap) > fRapidityMax) continue;
        fHistMixedInvMassPt->Fill(pT, invMass);
        if (!fIsLightweight) {
          const int rapBin = GetRapidityBin(rap);
          if (rapBin > -1) fHistMixedPtY[rapBin]->Fill(pT, invMass);
        }
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto &PhotonContainer : fPhotonMixed) {
    for (const auto &Photon : PhotonContainer) {
      if (Photon.GetPt() > fPhotonPtMax || Photon.GetPt() < fPhotonPtMin)
        continue;
      for (const auto &Lambda : lambdaCandidates) {
        if (!Lambda.GetIsUse()) continue;
        AliSigma0ParticlePhotonMother sigma(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma.GetArmenterosAlpha();
        const float armQt = sigma.GetArmenterosQt();
        const float pT = sigma.GetPt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
            continue;
        }
        const float invMass = sigma.GetMass();
        const float rap = sigma.GetRapidity();
        if (std::abs(rap) > fRapidityMax) continue;
        fHistMixedInvMassPt->Fill(pT, invMass);
        if (!fIsLightweight) {
          const int rapBin = GetRapidityBin(rap);
          if (rapBin > -1) fHistMixedPtY[rapBin]->Fill(pT, invMass);
        }
      }
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGammaMixedEventBinned(
    const std::vector<AliSigma0ParticleV0> &photonCandidates,
    const std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  // Mulitplicity estimator: V0M
  Float_t lPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection =
      (AliMultSelection *)fInputEvent->FindListObject("MultSelection");
  if (MultSelection) {
    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }

  const float zVertex = fInputEvent->GetPrimaryVertex()->GetZ();

  const int zVertexBin = GetZvertexBin(zVertex);
  const int multBin = GetMultiplicityBin(lPercentile);

  if (zVertexBin < 0 || multBin < 0) return;

  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fLambdaMixedBinned[zVertexBin][multBin]) {
    for (const auto &Lambda : LambdaContainer) {
      if (!Lambda.GetIsUse()) continue;
      for (auto Photon : photonCandidates) {
        if (Photon.GetPt() > fPhotonPtMax || Photon.GetPt() < fPhotonPtMin)
          continue;
        AliSigma0ParticlePhotonMother sigma(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma.GetArmenterosAlpha();
        const float armQt = sigma.GetArmenterosQt();
        const float pT = sigma.GetPt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
            continue;
        }
        const float invMass = sigma.GetMass();
        const float rap = sigma.GetRapidity();
        if (std::abs(rap) > fRapidityMax) continue;
        fHistMixedInvMassBinnedPt->Fill(pT, invMass);
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto &PhotonContainer : fPhotonMixedBinned[zVertexBin][multBin]) {
    for (const auto &Photon : PhotonContainer) {
      if (Photon.GetPt() > fPhotonPtMax || Photon.GetPt() < fPhotonPtMin)
        continue;
      for (const auto &Lambda : lambdaCandidates) {
        if (!Lambda.GetIsUse()) continue;
        AliSigma0ParticlePhotonMother sigma(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma.GetArmenterosAlpha();
        const float armQt = sigma.GetArmenterosQt();
        const float pT = sigma.GetPt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
            continue;
        }
        const float invMass = sigma.GetMass();
        const float rap = sigma.GetRapidity();
        if (std::abs(rap) > fRapidityMax) continue;
        fHistMixedInvMassBinnedPt->Fill(pT, invMass);
      }
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::FillEventBuffer(
    const std::vector<AliSigma0ParticleV0> &photonCandidates,
    const std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  // ++++++++++++++
  // Photon
  if (static_cast<int>(photonCandidates.size()) > 0) {
    if (static_cast<int>(fPhotonMixed.size()) < fMixingDepth) {
      fPhotonMixed.push_back(photonCandidates);
    } else {
      fPhotonMixed.pop_front();
      fPhotonMixed.push_back(photonCandidates);
    }
  }

  // ++++++++++++++
  // Lambda
  if (static_cast<int>(lambdaCandidates.size()) > 0) {
    if (static_cast<int>(fLambdaMixed.size()) < fMixingDepth) {
      fLambdaMixed.push_back(lambdaCandidates);
    } else {
      fLambdaMixed.pop_front();
      fLambdaMixed.push_back(lambdaCandidates);
    }
  }

  // +++++++++++++++
  // Binned mixed event in multiplicity and z-vertex position

  // Mulitplicity estimator: V0M
  Float_t lPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection =
      (AliMultSelection *)fInputEvent->FindListObject("MultSelection");
  if (MultSelection) {
    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }

  const float zVertex = fInputEvent->GetPrimaryVertex()->GetZ();

  const int zVertexBin = GetZvertexBin(zVertex);
  const int multBin = GetMultiplicityBin(lPercentile);

  if (zVertexBin < 0 || multBin < 0) return;

  // Photon
  if (static_cast<int>(photonCandidates.size()) > 0) {
    if (static_cast<int>(fPhotonMixedBinned[zVertexBin][multBin].size()) <
        fMixingDepth) {
      fPhotonMixedBinned[zVertexBin][multBin].push_back(photonCandidates);
    } else {
      fPhotonMixedBinned[zVertexBin][multBin].pop_front();
      fPhotonMixedBinned[zVertexBin][multBin].push_back(photonCandidates);
    }
  }

  // ++++++++++++++
  // Lambda
  if (static_cast<int>(lambdaCandidates.size()) > 0) {
    if (static_cast<int>(fLambdaMixedBinned[zVertexBin][multBin].size()) <
        fMixingDepth) {
      fLambdaMixedBinned[zVertexBin][multBin].push_back(lambdaCandidates);
    } else {
      fLambdaMixedBinned[zVertexBin][multBin].pop_front();
      fLambdaMixedBinned[zVertexBin][multBin].push_back(lambdaCandidates);
    }
  }
}

//____________________________________________________________________________________________________
float AliSigma0PhotonMotherCuts::ComputeRapidity(float pt, float pz,
                                                 float m) const {
  // calculates rapidity keeping the sign in case E == pz

  float energy = std::sqrt(pt * pt + pz * pz + m * m);
  if (energy != std::fabs(pz))
    return 0.5 * std::log((energy + pz) / (energy - pz));
  return (pz >= 0) ? 1.e30 : -1.e30;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::ProcessMC() const {
  // Mulitplicity estimator: V0M
  Float_t lPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection =
      (AliMultSelection *)fInputEvent->FindListObject("MultSelection");
  if (MultSelection) {
    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }

  // Loop over the MC tracks
  for (int iPart = 1; iPart < (fMCEvent->GetNumberOfTracks()); iPart++) {
    AliMCParticle *mcParticle =
        static_cast<AliMCParticle *>(fMCEvent->GetTrack(iPart));
    if (!mcParticle) continue;
    //    if (!mcParticle->IsPhysicalPrimary()) continue;
    if (mcParticle->GetNDaughters() != 2) continue;
    if (mcParticle->PdgCode() != fPDG) continue;
    fHistMCTruthPtY->Fill(mcParticle->Y(), mcParticle->Pt());
    fHistMCTruthPtEta->Fill(mcParticle->Eta(), mcParticle->Pt());
    if (lPercentile < fMCHighMultThreshold) {
      fHistMCTruthPtYHighMult->Fill(mcParticle->Y(), mcParticle->Pt());
      fHistMCTruthPtEtaHighMult->Fill(mcParticle->Eta(), mcParticle->Pt());
    }

    if (!CheckDaughters(mcParticle)) continue;
    fHistMCTruthDaughterPtY->Fill(mcParticle->Y(), mcParticle->Pt());
    fHistMCTruthDaughterPtEta->Fill(mcParticle->Eta(), mcParticle->Pt());
    if (lPercentile < fMCHighMultThreshold) {
      fHistMCTruthDaughterPtYHighMult->Fill(mcParticle->Y(), mcParticle->Pt());
      fHistMCTruthDaughterPtEtaHighMult->Fill(mcParticle->Eta(),
                                              mcParticle->Pt());
    }

    if (!CheckDaughtersInAcceptance(mcParticle)) continue;
    fHistMCTruthDaughterPtYAccept->Fill(mcParticle->Y(), mcParticle->Pt());
    fHistMCTruthDaughterPtEtaAccept->Fill(mcParticle->Eta(), mcParticle->Pt());
    if (lPercentile < fMCHighMultThreshold) {
      fHistMCTruthDaughterPtYAcceptHighMult->Fill(mcParticle->Y(),
                                                  mcParticle->Pt());
      fHistMCTruthDaughterPtEtaAcceptHighMult->Fill(mcParticle->Eta(),
                                                    mcParticle->Pt());
    }
  }
}

//____________________________________________________________________________________________________
bool AliSigma0PhotonMotherCuts::CheckDaughters(
    const AliMCParticle *particle) const {
  AliMCParticle *posDaughter = nullptr;
  AliMCParticle *negDaughter = nullptr;

  if (particle->GetNDaughters() != 2) return false;

  for (int daughterIndex = particle->GetFirstDaughter();
       daughterIndex <= particle->GetLastDaughter(); ++daughterIndex) {
    if (daughterIndex < 0) continue;
    AliMCParticle *tmpDaughter =
        static_cast<AliMCParticle *>(fMCEvent->GetTrack(daughterIndex));
    if (!tmpDaughter) continue;
    const int pdgCode = tmpDaughter->PdgCode();
    if (pdgCode == fPDGDaughter1)
      posDaughter = tmpDaughter;
    else if (pdgCode == fPDGDaughter2)
      negDaughter = tmpDaughter;
  }

  if (!posDaughter || !negDaughter) return false;

  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0PhotonMotherCuts::CheckDaughtersInAcceptance(
    const AliMCParticle *particle) const {
  AliMCParticle *lambdaDaughter = nullptr;
  AliMCParticle *photonDaughter = nullptr;

  if (std::abs(particle->Y()) > fRapidityMax) return false;

  if (particle->GetNDaughters() != 2) return false;

  for (int daughterIndex = particle->GetFirstDaughter();
       daughterIndex <= particle->GetLastDaughter(); ++daughterIndex) {
    if (daughterIndex < 0) continue;
    AliMCParticle *tmpDaughter =
        static_cast<AliMCParticle *>(fMCEvent->GetTrack(daughterIndex));
    if (!tmpDaughter) continue;
    const int pdgCode = tmpDaughter->PdgCode();
    if (pdgCode == fPDGDaughter1) {
      lambdaDaughter = tmpDaughter;
    } else if (pdgCode == fPDGDaughter2) {
      photonDaughter = tmpDaughter;
    }
  }

  if (!lambdaDaughter || !photonDaughter) {
    return false;
  }

  if (photonDaughter->Pt() > fPhotonPtMax ||
      photonDaughter->Pt() < fPhotonPtMin)
    return false;

  // Check Daughters in acceptance
  if (fLambdaCuts) {
    if (!fLambdaCuts->CheckDaughtersInAcceptance(lambdaDaughter)) return false;
  }
  if (fPhotonCuts) {
    if (!fPhotonCuts->CheckDaughtersInAcceptance(photonDaughter)) return false;
  } else if (fV0Reader) {
    AliConversionPhotonCuts *photonCuts = fV0Reader->GetConversionCuts();
    if (photonCuts) {
      if (!photonCuts->PhotonIsSelectedMC(photonDaughter->Particle(), fMCEvent))
        return false;
    }
  }

  return true;
}

//____________________________________________________________________________________________________
int AliSigma0PhotonMotherCuts::GetRapidityBin(float rapidity) const {
  if (-10 < rapidity && rapidity <= -1)
    return 0;
  else if (-1 < rapidity && rapidity <= -0.5)
    return 1;
  else if (-0.5 < rapidity && rapidity <= -0.4)
    return 2;
  else if (-0.4 < rapidity && rapidity <= -0.3)
    return 3;
  else if (-0.3 < rapidity && rapidity <= -0.2)
    return 4;
  else if (-0.2 < rapidity && rapidity <= -0.1)
    return 5;
  else if (-0.1 < rapidity && rapidity <= 0.f)
    return 6;
  else if (0.f < rapidity && rapidity <= 0.1)
    return 7;
  else if (0.1 < rapidity && rapidity <= 0.2)
    return 8;
  else if (0.2 < rapidity && rapidity <= 0.3)
    return 9;
  else if (0.3 < rapidity && rapidity <= 0.4)
    return 10;
  else if (0.4 < rapidity && rapidity <= 0.5)
    return 11;
  else if (0.5 < rapidity && rapidity <= 1.f)
    return 12;
  else if (1.0 < rapidity && rapidity < 10.f)
    return 13;
  else
    return -1;
}

//____________________________________________________________________________________________________
int AliSigma0PhotonMotherCuts::GetMultiplicityBin(float percentile) const {
  if (50 < percentile && percentile <= 100)
    return 0;
  else if (30 < percentile && percentile <= 50)
    return 1;
  else if (10 < percentile && percentile <= 30)
    return 2;
  else if (5 < percentile && percentile <= 10)
    return 3;
  else if (1 < percentile && percentile <= 5)
    return 4;
  else if (percentile <= 1)
    return 5;
  else
    return -1;
}

//____________________________________________________________________________________________________
int AliSigma0PhotonMotherCuts::GetZvertexBin(float zVertex) const {
  if (-10.f < zVertex && zVertex <= -8.f)
    return 0;
  else if (-8.f < zVertex && zVertex <= -6.f)
    return 1;
  else if (-6.f < zVertex && zVertex <= -4.f)
    return 2;
  else if (-4.f < zVertex && zVertex <= -2.f)
    return 3;
  else if (-2.f < zVertex && zVertex <= 0)
    return 4;
  else if (0.f < zVertex && zVertex <= 2.f)
    return 5;
  else if (2.f < zVertex && zVertex <= 4.f)
    return 6;
  else if (4.f < zVertex && zVertex <= 6.f)
    return 7;
  else if (6.f < zVertex && zVertex <= 8.f)
    return 8;
  else if (8.f < zVertex && zVertex <= 10.f)
    return 9;
  else
    return -1;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::InitCutHistograms(TString appendix) {
  fMassSigma = fDataBasePDG.GetParticle(fPDG)->Mass();

  std::cout << "============================\n"
            << " PHOTON MOTHER CUT CONFIGURATION \n"
            << " Sigma0 mass      " << fMassSigma << "\n"
            << " Sigma0 selection " << fSigmaMassCut << "\n"
            << " Sigma0 sb up     " << fSidebandCutUp << "\n"
            << " Sigma0 sb down   " << fSidebandCutDown << "\n"
            << " Photon pT min    " << fPhotonPtMin << "\n"
            << " Photon pT max    " << fPhotonPtMax << "\n"
            << "============================\n";

  std::vector<float> rapBins = {{-10., -1.f, -0.5, -0.4, -0.3, -0.2, -0.1, 0.f,
                                 0.1, 0.2, 0.3, 0.4, 0.5, 1.f, 10.f}};

  TH1::AddDirectory(kFALSE);

  TString name;
  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    fHistograms->SetName(appendix);
  }

  fHistCutBooking = new TProfile("fHistCutBooking", ";;Cut value", 12, 0, 12);
  fHistCutBooking->GetXaxis()->SetBinLabel(1, "#Sigma^{0} selection");
  fHistCutBooking->GetXaxis()->SetBinLabel(2, "#Sigma^{0} sb down");
  fHistCutBooking->GetXaxis()->SetBinLabel(3, "#Sigma^{0} sb up");
  fHistCutBooking->GetXaxis()->SetBinLabel(4, "#gamma #it{p}_{T} min");
  fHistCutBooking->GetXaxis()->SetBinLabel(5, "#gamma #it{p}_{T} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(6, "#Sigma^{0} mixing depth");
  fHistCutBooking->GetXaxis()->SetBinLabel(7, "Armenteros q_{T} low");
  fHistCutBooking->GetXaxis()->SetBinLabel(8, "Armenteros q_{T} up");
  fHistCutBooking->GetXaxis()->SetBinLabel(9, "Armenteros #alpha low");
  fHistCutBooking->GetXaxis()->SetBinLabel(10, "Armenteros #alpha up");
  fHistCutBooking->GetXaxis()->SetBinLabel(11, "Rapidity y max");
  fHistCutBooking->GetXaxis()->SetBinLabel(12, "MC Mult for efficiency");
  fHistograms->Add(fHistCutBooking);

  fHistCutBooking->Fill(0.f, fSigmaMassCut);
  fHistCutBooking->Fill(1.f, fSidebandCutDown);
  fHistCutBooking->Fill(2.f, fSidebandCutUp);
  fHistCutBooking->Fill(3.f, fPhotonPtMin);
  fHistCutBooking->Fill(4.f, fPhotonPtMax);
  fHistCutBooking->Fill(5.f, fMixingDepth);
  fHistCutBooking->Fill(6.f, fArmenterosQtLow);
  fHistCutBooking->Fill(7.f, fArmenterosQtUp);
  fHistCutBooking->Fill(8.f, fArmenterosAlphaLow);
  fHistCutBooking->Fill(9.f, fArmenterosAlphaUp);
  fHistCutBooking->Fill(10.f, fRapidityMax);
  fHistCutBooking->Fill(11.f, fMCHighMultThreshold);

  fHistInvMassPt = new TH2F("fHistInvMassPt",
                            "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                            "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                            100, 0, 10, 500, 1., 1.5);
  fHistograms->Add(fHistInvMassPt);
  fHistMixedInvMassPt = new TH2F("fHistMixedInvMassPt",
                                 "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                                 "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                                 100, 0, 10, 500, 1., 1.5);
  fHistograms->Add(fHistMixedInvMassPt);
  fHistMixedInvMassBinnedPt =
      new TH2F("fHistMixedInvMassBinnedPt",
               "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
               "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
               100, 0, 10, 500, 1., 1.5);
  fHistograms->Add(fHistMixedInvMassBinnedPt);

  fHistInvMass =
      new TH1F("fHistInvMass", "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries",
               150, 1.15, 1.3);
  fHistograms->Add(fHistInvMass);

  if (!fIsLightweight) {
    fHistNSigma =
        new TH1F("fHistNSigma", ";# #Sigma candidates; Entries", 10, 0, 10);
    fHistMassCutPt = new TH1F(
        "fHistMassCutPt", "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries",
        100, 0, 10);

    fHistInvMassBeforeArmenteros = new TH1F(
        "fHistInvMassBeforeArmenteros",
        "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 150, 1.15, 1.3);

    fHistInvMassRecPhoton = new TH2F("fHistInvMassRecPhoton",
                                     "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                                     "M_{#Lambda#gamma_{rec}} (GeV/#it{c}^{2})",
                                     50, 0, 10, 150, 1.15, 1.3);
    fHistInvMassRecLambda = new TH2F("fHistInvMassRecLambda",
                                     "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                                     "M_{#Lambda_{rec}#gamma} (GeV/#it{c}^{2})",
                                     50, 0, 10, 150, 1.15, 1.3);
    fHistInvMassRec = new TH2F("fHistInvMassRec",
                               "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                               "M_{#Lambda_{rec}#gamma_{rec}} (GeV/#it{c}^{2})",
                               50, 0, 10, 150, 1.15, 1.3);
    fHistInvMassEta = new TH2F("fHistInvMassEta",
                               "; #eta; M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                               250, -1, 1, 150, 1.15, 1.3);
    fHistArmenterosBefore =
        new TH2F("fHistArmenterosBefore", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 200, -1, 1, 100, 0, 0.5);
    fHistArmenterosAfter =
        new TH2F("fHistArmenterosAfter", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 200, -1, 1, 100, 0, 0.5);
    fHistEtaPhi = new TH2F("fHistEtaPhi", "; #eta; #phi", 100, -1, 1, 100,
                           -TMath::Pi(), TMath::Pi());

    fHistograms->Add(fHistNSigma);
    fHistograms->Add(fHistMassCutPt);
    fHistograms->Add(fHistInvMassBeforeArmenteros);
    fHistograms->Add(fHistInvMassRecPhoton);
    fHistograms->Add(fHistInvMassRecLambda);
    fHistograms->Add(fHistInvMassRec);
    fHistograms->Add(fHistInvMassEta);
    fHistograms->Add(fHistEtaPhi);
    fHistograms->Add(fHistArmenterosBefore);
    fHistograms->Add(fHistArmenterosAfter);

    for (int i = 0; i < static_cast<int>(rapBins.size() - 1); ++i) {
      fHistPtY[i] =
          new TH2F(Form("fHistPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
                   Form("%.2f < y < %.2f ; #it{p}_{T} (GeV/#it{c}); "
                        "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                        rapBins[i], rapBins[i + 1]),
                   100, 0, 10, 150, 1.15, 1.3);
      fHistMixedPtY[i] =
          new TH2F(Form("fHistMixedPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
                   Form("%.2f < y < %.2f ; #it{p}_{T} (GeV/#it{c}); "
                        "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                        rapBins[i], rapBins[i + 1]),
                   100, 0, 10, 150, 1.15, 1.3);
      fHistograms->Add(fHistPtY[i]);
      fHistograms->Add(fHistMixedPtY[i]);
    }

    fHistLambdaPtPhi =
        new TH2F("fHistLambdaPtPhi", "; #it{p}_{T} (GeV/#it{c}); #phi (rad)",
                 100, 0, 10, 100, 0, 2.f * TMath::Pi());
    fHistLambdaPtEta =
        new TH2F("fHistLambdaPtEta", "; #it{p}_{T} (GeV/#it{c}); #eta", 100, 0,
                 10, 100, -1, 1);
    fHistLambdaMassPt =
        new TH2F("fHistLambdaMassPt", "; #it{p}_{T} (GeV/#it{c}); M_{p#pi}",
                 100, 0, 10, 250, 1., 1.5);
    fHistPhotonPtPhi =
        new TH2F("fHistPhotonPtPhi", "; #it{p}_{T} (GeV/#it{c}); #phi (rad)",
                 100, 0, 10, 100, 0, 2.f * TMath::Pi());
    fHistPhotonPtEta =
        new TH2F("fHistPhotonPtEta", "; #it{p}_{T} (GeV/#it{c}); #eta", 100, 0,
                 10, 100, -1, 1);
    fHistPhotonMassPt =
        new TH2F("fHistPhotonMassPt", "; #it{p}_{T} (GeV/#it{c}); M_{p#pi}",
                 100, 0, 10, 250, 0., 0.1);

    fHistograms->Add(fHistLambdaPtPhi);
    fHistograms->Add(fHistLambdaPtEta);
    fHistograms->Add(fHistLambdaMassPt);
    fHistograms->Add(fHistPhotonPtPhi);
    fHistograms->Add(fHistPhotonPtEta);
    fHistograms->Add(fHistPhotonMassPt);
  }

  if (fIsMC) {
    if (fHistogramsMC != nullptr) {
      delete fHistogramsMC;
      fHistogramsMC = nullptr;
    }
    if (fHistogramsMC == nullptr) {
      fHistogramsMC = new TList();
      fHistogramsMC->SetOwner(kTRUE);
      fHistogramsMC->SetName("MC");
    }

    fHistMCTruthPtY =
        new TH2F("fHistMCTruthPtY", "; y; #it{p}_{T} (GeV/#it{c})", 100, -10,
                 10, 100, 0, 10);
    fHistMCTruthPtEta =
        new TH2F("fHistMCTruthPtEta", "; #eta; #it{p}_{T} (GeV/#it{c})", 100,
                 -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtY =
        new TH2F("fHistMCTruthDaughterPtY", "; y; #it{p}_{T} (GeV/#it{c})", 100,
                 -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtEta =
        new TH2F("fHistMCTruthDaughterPtEta", "; #eta; #it{p}_{T} (GeV/#it{c})",
                 100, -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtYAccept =
        new TH2F("fHistMCTruthDaughterPtYAccept",
                 "; y; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtEtaAccept =
        new TH2F("fHistMCTruthDaughterPtEtaAccept",
                 "; #eta; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);
    fHistMCTruthPtYHighMult =
        new TH2F("fHistMCTruthPtYHighMult", "; y; #it{p}_{T} (GeV/#it{c})", 100,
                 -10, 10, 100, 0, 10);
    fHistMCTruthPtEtaHighMult =
        new TH2F("fHistMCTruthPtEtaHighMult", "; #eta; #it{p}_{T} (GeV/#it{c})",
                 100, -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtYHighMult =
        new TH2F("fHistMCTruthDaughterPtYHighMult",
                 "; y; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtEtaHighMult =
        new TH2F("fHistMCTruthDaughterPtEtaHighMult",
                 "; #eta; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtYAcceptHighMult =
        new TH2F("fHistMCTruthDaughterPtYAcceptHighMult",
                 "; y; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtEtaAcceptHighMult =
        new TH2F("fHistMCTruthDaughterPtEtaAcceptHighMult",
                 "; #eta; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);

    fHistMCV0Pt = new TH1F("fHistMCV0Pt",
                           "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries",
                           100, 0, 10);
    fHistMCV0Mass =
        new TH1F("fHistMCV0Mass",
                 "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 500, 1., 1.5);
    fHistMCV0Mother =
        new TH2F("fHistMCV0Mother",
                 "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); PDG code mother", 250,
                 1., 2., 4000, 0, 4000);
    fHistMCV0Check =
        new TH2F("fHistMCV0Check",
                 "; PDG code #Lambda candidate; PDG code #gamma candidate",
                 4000, 0, 4000, 4000, 0, 4000);
    fHistMCV0MotherCheck =
        new TH2F("fHistMCV0MotherCheck",
                 "; PDG code #Lambda mother; PDG code #gamma mother", 4000, 0,
                 4000, 4000, 0, 4000);

    fHistogramsMC->Add(fHistMCTruthPtY);
    fHistogramsMC->Add(fHistMCTruthPtEta);
    fHistogramsMC->Add(fHistMCTruthDaughterPtY);
    fHistogramsMC->Add(fHistMCTruthDaughterPtEta);
    fHistogramsMC->Add(fHistMCTruthDaughterPtYAccept);
    fHistogramsMC->Add(fHistMCTruthDaughterPtEtaAccept);
    fHistogramsMC->Add(fHistMCTruthPtYHighMult);
    fHistogramsMC->Add(fHistMCTruthPtEtaHighMult);
    fHistogramsMC->Add(fHistMCTruthDaughterPtYHighMult);
    fHistogramsMC->Add(fHistMCTruthDaughterPtEtaHighMult);
    fHistogramsMC->Add(fHistMCTruthDaughterPtYAcceptHighMult);
    fHistogramsMC->Add(fHistMCTruthDaughterPtEtaAcceptHighMult);
    fHistogramsMC->Add(fHistMCV0Pt);
    fHistogramsMC->Add(fHistMCV0Mass);
    fHistogramsMC->Add(fHistMCV0Mother);
    fHistogramsMC->Add(fHistMCV0Check);
    fHistogramsMC->Add(fHistMCV0MotherCheck);

    fHistograms->Add(fHistogramsMC);
  }
}
