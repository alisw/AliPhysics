#include "AliSigma0PhotonMotherCuts.h"
#include <iostream>
#include "AliMultSelection.h"
#include "TMath.h"

ClassImp(AliSigma0PhotonMotherCuts)

    //____________________________________________________________________________________________________
    AliSigma0PhotonMotherCuts::AliSigma0PhotonMotherCuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fIsMC(false),
      fDoCleanUp(true),
      fIsLightweight(false),
      fIsSpectrumAnalysis(true),
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
      fLambdaCuts(nullptr),
      fPhotonCuts(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fMultMode(AliVEvent::kINT7),
      fMixingDepth(10),
      fPDG(0),
      fPDGDaughter1(0),
      fPDGDaughter2(0),
      fMassWindowPt(false),
      fMassWindowP0(1.19257f),
      fMassWindowP1(3.66143e-03),
      fMassWindowP2(-0.772265f),
      fMassSigma(0),
      fSigmaMassCut(0.003),
      fSidebandCutUp(0.05),
      fSidebandCutDown(0.01),
      fPhotonPtMin(0),
      fPhotonPtMax(999.f),
      fPtMin(0.f),
      fRapidityMax(0.5),
      fArmenterosCut(false),
      fArmenterosQtLow(0.f),
      fArmenterosQtUp(0.f),
      fArmenterosAlphaLow(0.f),
      fArmenterosAlphaUp(0.f),
      fMCHighMultThreshold(5.f),
      fHistCutBooking(nullptr),
      fHistNSigma(nullptr),
      fHistNPhotonBefore(nullptr),
      fHistNPhotonAfter(nullptr),
      fHistNLambdaBefore(nullptr),
      fHistNLambdaAfter(nullptr),
      fHistNPhotonLabel(nullptr),
      fHistNLambdaLabel(nullptr),
      fHistNLambdaGammaLabel(nullptr),
      fHistMassCutPt(nullptr),
      fHistInvMass(nullptr),
      fHistInvMassK0Gamma(nullptr),
      fHistInvMassSelected(nullptr),
      fHistInvMassRecPhoton(nullptr),
      fHistInvMassRecLambda(nullptr),
      fHistInvMassRec(nullptr),
      fHistInvMassPt(nullptr),
      fHistInvMassPtRaw(nullptr),
      fHistEtaPhi(nullptr),
      fHistPtRapidity(nullptr),
      fHistPtMult(),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistMixedInvMassPt(nullptr),
      fHistMixedInvMassBinnedMultPt(),
      fHistDeltaEtaDeltaPhiGammaNegBefore(nullptr),
      fHistDeltaEtaDeltaPhiGammaPosBefore(nullptr),
      fHistDeltaEtaDeltaPhiLambdaNegBefore(nullptr),
      fHistDeltaEtaDeltaPhiLambdaPosBefore(nullptr),
      fHistDeltaEtaDeltaPhiLambdaGammaNegBefore(nullptr),
      fHistDeltaEtaDeltaPhiLambdaGammaPosBefore(nullptr),
      fHistDeltaEtaDeltaPhiGammaNegAfter(nullptr),
      fHistDeltaEtaDeltaPhiGammaPosAfter(nullptr),
      fHistDeltaEtaDeltaPhiLambdaNegAfter(nullptr),
      fHistDeltaEtaDeltaPhiLambdaPosAfter(nullptr),
      fHistDeltaEtaDeltaPhiLambdaGammaNegAfter(nullptr),
      fHistDeltaEtaDeltaPhiLambdaGammaPosAfter(nullptr),
      fHistLambdaPtPhi(nullptr),
      fHistLambdaPtEta(nullptr),
      fHistLambdaMassPt(nullptr),
      fHistPhotonPtPhi(nullptr),
      fHistPhotonPtEta(nullptr),
      fHistPhotonMassPt(nullptr),
      fHistSigmaLambdaPtCorr(nullptr),
      fHistSigmaPhotonPtCorr(nullptr),
      fHistSigmaLambdaPCorr(nullptr),
      fHistSigmaPhotonPCorr(nullptr),
      fHistMCTruthPt(nullptr),
      fHistMCTruthPtMult(),
      fHistMCTruthPtY(nullptr),
      fHistMCTruthDaughterPtY(nullptr),
      fHistMCTruthDaughterPtYAccept(nullptr),
      fHistMCTruthPtYHighMult(nullptr),
      fHistMCTruthDaughterPtYHighMult(nullptr),
      fHistMCTruthDaughterPtYAcceptHighMult(nullptr),
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
      fDoCleanUp(true),
      fIsLightweight(false),
      fIsSpectrumAnalysis(true),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fSigma(),
      fSidebandUp(),
      fSidebandDown(),
      fLambdaMixed(),
      fPhotonMixed(),
      fLambdaCuts(nullptr),
      fPhotonCuts(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fMultMode(AliVEvent::kINT7),
      fMixingDepth(10),
      fPDG(0),
      fPDGDaughter1(0),
      fPDGDaughter2(0),
      fMassWindowPt(false),
      fMassWindowP0(1.19257f),
      fMassWindowP1(3.66143e-03),
      fMassWindowP2(-0.772265f),
      fMassSigma(0),
      fSigmaMassCut(0.003),
      fSidebandCutUp(0.05),
      fSidebandCutDown(0.01),
      fPhotonPtMin(0),
      fPhotonPtMax(999.f),
      fPtMin(0.f),
      fRapidityMax(0.5),
      fArmenterosCut(false),
      fArmenterosQtLow(0.f),
      fArmenterosQtUp(0.f),
      fArmenterosAlphaLow(0.f),
      fArmenterosAlphaUp(0.f),
      fMCHighMultThreshold(5.f),
      fHistCutBooking(nullptr),
      fHistNSigma(nullptr),
      fHistNPhotonBefore(nullptr),
      fHistNPhotonAfter(nullptr),
      fHistNLambdaBefore(nullptr),
      fHistNLambdaAfter(nullptr),
      fHistNPhotonLabel(nullptr),
      fHistNLambdaLabel(nullptr),
      fHistNLambdaGammaLabel(nullptr),
      fHistMassCutPt(nullptr),
      fHistInvMass(nullptr),
      fHistInvMassK0Gamma(nullptr),
      fHistInvMassSelected(nullptr),
      fHistInvMassRecPhoton(nullptr),
      fHistInvMassRecLambda(nullptr),
      fHistInvMassRec(nullptr),
      fHistInvMassPt(nullptr),
      fHistInvMassPtRaw(nullptr),
      fHistEtaPhi(nullptr),
      fHistPtRapidity(nullptr),
      fHistPtMult(),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistMixedInvMassPt(nullptr),
      fHistMixedInvMassBinnedMultPt(),
      fHistDeltaEtaDeltaPhiGammaNegBefore(nullptr),
      fHistDeltaEtaDeltaPhiGammaPosBefore(nullptr),
      fHistDeltaEtaDeltaPhiLambdaNegBefore(nullptr),
      fHistDeltaEtaDeltaPhiLambdaPosBefore(nullptr),
      fHistDeltaEtaDeltaPhiLambdaGammaNegBefore(nullptr),
      fHistDeltaEtaDeltaPhiLambdaGammaPosBefore(nullptr),
      fHistDeltaEtaDeltaPhiGammaNegAfter(nullptr),
      fHistDeltaEtaDeltaPhiGammaPosAfter(nullptr),
      fHistDeltaEtaDeltaPhiLambdaNegAfter(nullptr),
      fHistDeltaEtaDeltaPhiLambdaPosAfter(nullptr),
      fHistDeltaEtaDeltaPhiLambdaGammaNegAfter(nullptr),
      fHistDeltaEtaDeltaPhiLambdaGammaPosAfter(nullptr),
      fHistLambdaPtPhi(nullptr),
      fHistLambdaPtEta(nullptr),
      fHistLambdaMassPt(nullptr),
      fHistPhotonPtPhi(nullptr),
      fHistPhotonPtEta(nullptr),
      fHistPhotonMassPt(nullptr),
      fHistSigmaLambdaPtCorr(nullptr),
      fHistSigmaPhotonPtCorr(nullptr),
      fHistSigmaLambdaPCorr(nullptr),
      fHistSigmaPhotonPCorr(nullptr),
      fHistMCTruthPt(nullptr),
      fHistMCTruthPtMult(),
      fHistMCTruthPtY(nullptr),
      fHistMCTruthDaughterPtY(nullptr),
      fHistMCTruthDaughterPtYAccept(nullptr),
      fHistMCTruthPtYHighMult(nullptr),
      fHistMCTruthDaughterPtYHighMult(nullptr),
      fHistMCTruthDaughterPtYAcceptHighMult(nullptr),
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
  photonMotherCuts->SetSigmaMassPt(true);
  return photonMotherCuts;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SelectPhotonMother(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    std::vector<AliSigma0ParticleV0> &photonCandidates,
    std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  fInputEvent = inputEvent;
  fMCEvent = mcEvent;

  fSigma.clear();
  fSidebandUp.clear();
  fSidebandDown.clear();

  if (fIsMC) {
    ProcessMC();
    fV0Reader =
        (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
            fV0ReaderName.Data());
  }
  
  if (fDoCleanUp) {
  CleanUpClones(photonCandidates, lambdaCandidates);
  }

  if (!fIsLightweight) SingleV0QA(photonCandidates, lambdaCandidates);

  // Particle pairing
  SigmaToLambdaGamma(photonCandidates, lambdaCandidates);
  if (fIsSpectrumAnalysis) {
    SigmaToLambdaGammaMixedEvent(photonCandidates, lambdaCandidates);
    SigmaToLambdaGammaMixedEventBinned(photonCandidates, lambdaCandidates);
    FillEventBuffer(photonCandidates, lambdaCandidates);
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::CleanUpClones(
    std::vector<AliSigma0ParticleV0> &photonCandidates,
    std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  // Keep track of the fantastic situation before applying the cuts

  if (!fIsLightweight) {
    fHistNPhotonBefore->Fill(photonCandidates.size());
    fHistNLambdaBefore->Fill(lambdaCandidates.size());
  }

  // Do the checks for the photons
  int nPhotonKilledLabel = 0;
  for (auto photon1 = photonCandidates.begin();
       photon1 < photonCandidates.end(); ++photon1) {
    if (!photon1->GetIsUse()) continue;
    const auto posDaughter1 = photon1->GetPosDaughter();
    const auto negDaughter1 = photon1->GetNegDaughter();

    for (auto photon2 = photon1 + 1; photon2 < photonCandidates.end();
         ++photon2) {
      if (!photon1->GetIsUse() || !photon2->GetIsUse()) continue;
      const auto posDaughter2 = photon2->GetPosDaughter();
      const auto negDaughter2 = photon2->GetNegDaughter();

      if (!fIsLightweight) {
        float deltaPhistarPos =
            posDaughter1.GetAveragePhiStar() - posDaughter2.GetAveragePhiStar();
        float deltaPhistarNeg =
            negDaughter1.GetAveragePhiStar() - negDaughter2.GetAveragePhiStar();
        float deltaEtaPos = posDaughter1.GetEta() - posDaughter2.GetEta();
        float deltaEtaNeg = negDaughter1.GetEta() - negDaughter2.GetEta();

        fHistDeltaEtaDeltaPhiGammaNegBefore->Fill(deltaEtaNeg, deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiGammaPosBefore->Fill(deltaEtaPos, deltaPhistarPos);
      }

      bool hasSameLabels =
          (photon1->GetTrackLabelPos() == photon2->GetTrackLabelPos() ||
           photon1->GetTrackLabelNeg() == photon2->GetTrackLabelNeg() ||
           photon1->GetTrackLabelNeg() == photon2->GetTrackLabelPos() ||
           photon1->GetTrackLabelPos() == photon2->GetTrackLabelNeg());

      // do the check for both daughters
      if (hasSameLabels) {
        const float cpa1 = photon1->GetCosineAlpha();
        const float cpa2 = photon2->GetCosineAlpha();
        ++nPhotonKilledLabel;
        if (cpa1 > cpa2) {
          photon2->SetUse(false);
        } else {
          photon1->SetUse(false);
        }
      }
      if (!photon1->GetIsUse()) break;
    }
  }

  // Do the checks for the Lambdas
  int nLambdaKilledLabel = 0;
  for (auto lambda1 = lambdaCandidates.begin();
       lambda1 < lambdaCandidates.end(); ++lambda1) {
    if (!lambda1->GetIsUse()) continue;

    const auto posDaughter1 = lambda1->GetPosDaughter();
    const auto negDaughter1 = lambda1->GetNegDaughter();

    for (auto lambda2 = lambda1 + 1; lambda2 < lambdaCandidates.end();
         ++lambda2) {
      if (!lambda1->GetIsUse() || !lambda2->GetIsUse()) continue;

      const auto posDaughter2 = lambda2->GetPosDaughter();
      const auto negDaughter2 = lambda2->GetNegDaughter();

      if (!fIsLightweight) {
        float deltaPhistarPos =
            posDaughter1.GetAveragePhiStar() - posDaughter2.GetAveragePhiStar();
        float deltaPhistarNeg =
            negDaughter1.GetAveragePhiStar() - negDaughter2.GetAveragePhiStar();
        float deltaEtaPos = posDaughter1.GetEta() - posDaughter2.GetEta();
        float deltaEtaNeg = negDaughter1.GetEta() - negDaughter2.GetEta();
        fHistDeltaEtaDeltaPhiLambdaNegBefore->Fill(deltaEtaNeg,
                                                   deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiLambdaPosBefore->Fill(deltaEtaPos,
                                                   deltaPhistarPos);
      }

      bool hasSameLabels =
          (lambda1->GetTrackLabelPos() == lambda2->GetTrackLabelPos() ||
           lambda1->GetTrackLabelNeg() == lambda2->GetTrackLabelNeg() ||
           lambda1->GetTrackLabelNeg() == lambda2->GetTrackLabelPos() ||
           lambda1->GetTrackLabelPos() == lambda2->GetTrackLabelNeg());

      // do the check for both daughters
      if (hasSameLabels) {
        ++nLambdaKilledLabel;
        const float cpa1 = lambda1->GetCosineAlpha();
        const float cpa2 = lambda2->GetCosineAlpha();
        if (cpa1 > cpa2) {
          lambda2->SetUse(false);
        } else {
          lambda1->SetUse(false);
        }
      }
      if (!lambda1->GetIsUse()) break;
    }
  }

  // Now do the exercise for photon-Lambda combination
  int nPhotonLambdaKilledLabel = 0;
  for (auto photon = photonCandidates.begin(); photon < photonCandidates.end();
       ++photon) {
    if (!photon->GetIsUse()) continue;
    const auto posDaughterPhoton = photon->GetPosDaughter();
    const auto negDaughterPhoton = photon->GetNegDaughter();

    for (auto lambda = lambdaCandidates.begin();
         lambda < lambdaCandidates.end(); ++lambda) {
      if (!lambda->GetIsUse() || !photon->GetIsUse()) continue;
      const auto posDaughterLambda = lambda->GetPosDaughter();
      const auto negDaughterLambda = lambda->GetNegDaughter();

      if (!fIsLightweight) {
        float deltaPhistarPos = posDaughterPhoton.GetAveragePhiStar() -
                                posDaughterLambda.GetAveragePhiStar();
        float deltaPhistarNeg = negDaughterPhoton.GetAveragePhiStar() -
                                negDaughterLambda.GetAveragePhiStar();
        float deltaEtaPos =
            posDaughterPhoton.GetEta() - posDaughterLambda.GetEta();
        float deltaEtaNeg =
            negDaughterPhoton.GetEta() - negDaughterLambda.GetEta();

        fHistDeltaEtaDeltaPhiLambdaGammaNegBefore->Fill(deltaEtaNeg,
                                                        deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiLambdaGammaPosBefore->Fill(deltaEtaPos,
                                                        deltaPhistarPos);
      }

      bool hasSameLabels =
          (photon->GetTrackLabelPos() == lambda->GetTrackLabelPos() ||
           photon->GetTrackLabelNeg() == lambda->GetTrackLabelNeg() ||
           photon->GetTrackLabelNeg() == lambda->GetTrackLabelPos() ||
           photon->GetTrackLabelPos() == lambda->GetTrackLabelNeg());

      // do the check for both daughters
      if (hasSameLabels) {
        const float cpaPhoton = photon->GetCosineAlpha();
        const float cpaLambda = lambda->GetCosineAlpha();
        ++nPhotonLambdaKilledLabel;
        if (cpaPhoton > cpaLambda) {
          lambda->SetUse(false);
        } else {
          photon->SetUse(false);
        }
      }
      if (!photon->GetIsUse()) break;
    }
  }

  // We're done here, now comes the QA

  // Now take a look what changed for the Photons
  int nPhotonAfter = 0;
  for (auto photon1 = photonCandidates.begin();
       photon1 < photonCandidates.end(); ++photon1) {
    if (!photon1->GetIsUse()) continue;
    const auto posDaughter1 = photon1->GetPosDaughter();
    const auto negDaughter1 = photon1->GetNegDaughter();
    ++nPhotonAfter;

    for (auto photon2 = photon1 + 1; photon2 < photonCandidates.end();
         ++photon2) {
      if (!photon2->GetIsUse()) continue;
      const auto posDaughter2 = photon2->GetPosDaughter();
      const auto negDaughter2 = photon2->GetNegDaughter();

      if (!fIsLightweight) {
        float deltaPhistarPos =
            posDaughter1.GetAveragePhiStar() - posDaughter2.GetAveragePhiStar();
        float deltaPhistarNeg =
            negDaughter1.GetAveragePhiStar() - negDaughter2.GetAveragePhiStar();
        float deltaEtaPos = posDaughter1.GetEta() - posDaughter2.GetEta();
        float deltaEtaNeg = negDaughter1.GetEta() - negDaughter2.GetEta();
        fHistDeltaEtaDeltaPhiGammaNegAfter->Fill(deltaEtaNeg, deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiGammaPosAfter->Fill(deltaEtaPos, deltaPhistarPos);
      }
    }
  }

  // Now take a look what changed for the Lambdas
  int nLambdaAfter = 0;
  for (auto photon = photonCandidates.begin(); photon < photonCandidates.end();
       ++photon) {
    if (!photon->GetIsUse()) continue;
    const auto posDaughter1 = photon->GetPosDaughter();
    const auto negDaughter1 = photon->GetNegDaughter();

    for (auto lambda = lambdaCandidates.begin();
         lambda < lambdaCandidates.end(); ++lambda) {
      if (!lambda->GetIsUse()) continue;

      const auto posDaughter2 = lambda->GetPosDaughter();
      const auto negDaughter2 = lambda->GetNegDaughter();

      if (!fIsLightweight) {
        float deltaPhistarPos =
            posDaughter1.GetAveragePhiStar() - posDaughter2.GetAveragePhiStar();
        float deltaPhistarNeg =
            negDaughter1.GetAveragePhiStar() - negDaughter2.GetAveragePhiStar();
        float deltaEtaPos = posDaughter1.GetEta() - posDaughter2.GetEta();
        float deltaEtaNeg = negDaughter1.GetEta() - negDaughter2.GetEta();
        fHistDeltaEtaDeltaPhiLambdaGammaNegAfter->Fill(deltaEtaNeg,
                                                       deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiLambdaGammaPosAfter->Fill(deltaEtaPos,
                                                       deltaPhistarPos);
      }
    }
  }

  for (auto lambda1 = lambdaCandidates.begin();
       lambda1 < lambdaCandidates.end(); ++lambda1) {
    if (!lambda1->GetIsUse()) continue;
    ++nLambdaAfter;

    const auto posDaughter1 = lambda1->GetPosDaughter();
    const auto negDaughter1 = lambda1->GetNegDaughter();

    for (auto lambda2 = lambda1 + 1; lambda2 < lambdaCandidates.end();
         ++lambda2) {
      if (!lambda2->GetIsUse()) continue;

      const auto posDaughter2 = lambda2->GetPosDaughter();
      const auto negDaughter2 = lambda2->GetNegDaughter();

      if (!fIsLightweight) {
        float deltaPhistarPos =
            posDaughter1.GetAveragePhiStar() - posDaughter2.GetAveragePhiStar();
        float deltaPhistarNeg =
            negDaughter1.GetAveragePhiStar() - negDaughter2.GetAveragePhiStar();
        float deltaEtaPos = posDaughter1.GetEta() - posDaughter2.GetEta();
        float deltaEtaNeg = negDaughter1.GetEta() - negDaughter2.GetEta();
        fHistDeltaEtaDeltaPhiLambdaNegAfter->Fill(deltaEtaNeg, deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiLambdaPosAfter->Fill(deltaEtaPos, deltaPhistarPos);
      }
    }
  }

  if (!fIsLightweight) {
    fHistNPhotonAfter->Fill(nPhotonAfter);
    fHistNLambdaAfter->Fill(nLambdaAfter);
    fHistNPhotonLabel->Fill(nPhotonKilledLabel);
    fHistNLambdaLabel->Fill(nLambdaKilledLabel);
    fHistNLambdaGammaLabel->Fill(nPhotonLambdaKilledLabel);
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SingleV0QA(
    const std::vector<AliSigma0ParticleV0> &photonCandidates,
    const std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  // Photon QA - keep track of kinematics
  for (const auto &photon : photonCandidates) {
    if (!photon.GetIsUse()) continue;
    fHistPhotonPtPhi->Fill(photon.GetPt(), photon.GetPhi());
    fHistPhotonPtEta->Fill(photon.GetPt(), photon.GetEta());
    fHistPhotonMassPt->Fill(photon.GetPt(), photon.GetMass());
  }

  // Lambda QA - keep track of kinematics
  for (const auto &lambda : lambdaCandidates) {
    if (!lambda.GetIsUse()) continue;
    fHistLambdaPtPhi->Fill(lambda.GetPt(), lambda.GetPhi());
    fHistLambdaPtEta->Fill(lambda.GetPt(), lambda.GetEta());
    fHistLambdaMassPt->Fill(lambda.GetPt(), lambda.GetMass());
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGamma(
    const std::vector<AliSigma0ParticleV0> &photonCandidates,
    const std::vector<AliSigma0ParticleV0> &lambdaCandidates) {
  static float massK0 = TDatabasePDG::Instance()->GetParticle(311)->Mass();

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
    if (!photon.GetIsUse()) continue;

    for (const auto &lambda : lambdaCandidates) {
      if (!lambda.GetIsUse()) continue;
      AliSigma0ParticlePhotonMother sigma(lambda, photon, fInputEvent);
      const float invMass = sigma.GetMass();
      const float armAlpha = sigma.GetArmenterosAlpha();
      const float armQt = sigma.GetArmenterosQt();
      const float pT = sigma.GetPt();

      if (pT < fPtMin) continue;

      if (!fIsLightweight) {
        fHistArmenterosBefore->Fill(armAlpha, armQt);
      }
      // Armenteros cut
      if (fArmenterosCut) {
        if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
        if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
          continue;
      }
      const float rap = sigma.GetRapidity();
      const int multBin = GetMultiplicityBin(lPercentile, fMultMode);

      if (!fIsLightweight) {
        fHistArmenterosAfter->Fill(armAlpha, armQt);
        fHistInvMassRecPhoton->Fill(pT, sigma.GetRecMassPhoton());
        fHistInvMassRecLambda->Fill(pT, sigma.GetRecMassLambda());
        fHistInvMassRec->Fill(pT, sigma.GetRecMass());

        TLorentzVector track1, track2;
        track1.SetXYZM(lambda.GetPx(), lambda.GetPy(), lambda.GetPz(), massK0);
        track2.SetXYZM(photon.GetPx(), photon.GetPy(), photon.GetPz(), 0.f);
        TLorentzVector trackSum = track1 + track2;
        fHistInvMassK0Gamma->Fill(trackSum.M());
      }

      int label = -10;
      int pdgLambdaMother = 0;
      int pdgPhotonMother = 0;
      if (fIsMC) {
        label =
            sigma.MatchToMC(fMCEvent, fPDG, {{fPDGDaughter1, fPDGDaughter2}},
                            pdgLambdaMother, pdgPhotonMother);
      }

      // Now write out the stuff to the Femto containers
      if (invMass < GetMassSigmaPt(pT) + fSigmaMassCut &&
          invMass > GetMassSigmaPt(pT) - fSigmaMassCut) {
        fSigma.push_back(sigma);
        if (!fIsLightweight) {
          fHistInvMassSelected->Fill(pT, invMass);
          fHistMassCutPt->Fill(pT);
          fHistEtaPhi->Fill(sigma.GetEta(), sigma.GetPhi());
          fHistPtRapidity->Fill(pT, rap);
          fHistSigmaLambdaPtCorr->Fill(pT, lambda.GetPt());
          fHistSigmaPhotonPtCorr->Fill(pT, photon.GetPt());
          fHistSigmaLambdaPCorr->Fill(sigma.GetP(), lambda.GetP());
          fHistSigmaPhotonPCorr->Fill(sigma.GetP(), photon.GetP());
        }
        ++nSigma;
      }
      if (invMass < GetMassSigmaPt(pT) + fSidebandCutUp &&
          invMass > GetMassSigmaPt(pT) + fSidebandCutDown) {
        fSidebandUp.push_back(sigma);
        if (!fIsLightweight) {
          fHistInvMassSelected->Fill(pT, invMass);
        }
      }
      if (invMass > GetMassSigmaPt(pT) - fSidebandCutUp &&
          invMass < GetMassSigmaPt(pT) - fSidebandCutDown) {
        fSidebandDown.push_back(sigma);
        if (!fIsLightweight) {
          fHistInvMassSelected->Fill(pT, invMass);
        }
      }

      fHistInvMass->Fill(invMass);
      fHistInvMassPtRaw->Fill(pT, invMass);

      if (TMath::Abs(rap) < fRapidityMax) {
        fHistInvMassPt->Fill(pT, invMass);

        if (multBin >= 0 && fIsSpectrumAnalysis) {
          fHistPtMult[multBin]->Fill(pT, invMass);
        }
      }

      if (fIsMC) {
        if (label > 0) {
          fHistMCV0Pt->Fill(sigma.GetPt());
          fHistMCV0Mass->Fill(invMass);
        }
        if (!fIsLightweight) {
          // let's where the other particle comes from if one of them stems
          // from
          // a Sigma0
          if (TMath::Abs(pdgLambdaMother) == 3212 &&
              TMath::Abs(pdgLambdaMother) != 3212) {
            fHistMCV0Mother->Fill(invMass, TMath::Abs(pdgPhotonMother));
          }
          if (TMath::Abs(pdgLambdaMother) == 3212 &&
              TMath::Abs(pdgLambdaMother) != 3212) {
            fHistMCV0Mother->Fill(invMass, TMath::Abs(pdgLambdaMother));
          }
          fHistMCV0MotherCheck->Fill(TMath::Abs(pdgLambdaMother),
                                     TMath::Abs(pdgPhotonMother));

          const int labV0 = photon.GetMCLabelV0();
          const int labPhoton = lambda.GetMCLabelV0();
          if (labV0 < 0 || labPhoton < 0) continue;

          AliMCParticle *partV0 =
              static_cast<AliMCParticle *>(fMCEvent->GetTrack(labV0));
          AliMCParticle *partPhoton =
              static_cast<AliMCParticle *>(fMCEvent->GetTrack(labPhoton));
          if (!partV0 || !partPhoton) continue;

          fHistMCV0Check->Fill(TMath::Abs(partV0->PdgCode()),
                               TMath::Abs(partPhoton->PdgCode()));
        }
      }
    }
  }
  fHistNSigma->Fill(nSigma);
}

//____________________________________________________________________________________________________
float AliSigma0PhotonMotherCuts::GetMassSigmaPt(float pt) const {
  float massSigma = -1.f;
  if (fMassWindowPt) {
    massSigma = fMassWindowP0 + fMassWindowP1 * std::exp(fMassWindowP2 * pt);
  } else {
    massSigma = fMassSigma;
  }
  return massSigma;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGammaMixedEvent(
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
        const int multBin = GetMultiplicityBin(lPercentile, fMultMode);
        if (TMath::Abs(rap) > fRapidityMax || multBin < 0) continue;
        if (fIsSpectrumAnalysis) fHistMixedInvMassPt->Fill(pT, invMass);
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
        const int multBin = GetMultiplicityBin(lPercentile, fMultMode);
        if (TMath::Abs(rap) > fRapidityMax || multBin < 0) continue;
        if (fIsSpectrumAnalysis) fHistMixedInvMassPt->Fill(pT, invMass);
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

  const int multBin = GetMultiplicityBin(lPercentile, fMultMode);
  if (multBin < 0) return;

  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fLambdaMixedBinned[multBin]) {
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
        if (TMath::Abs(rap) > fRapidityMax) continue;
        if (fIsSpectrumAnalysis)
          fHistMixedInvMassBinnedMultPt[multBin]->Fill(pT, invMass);
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto &PhotonContainer : fPhotonMixedBinned[multBin]) {
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
        if (TMath::Abs(rap) > fRapidityMax) continue;
        if (fIsSpectrumAnalysis)
          fHistMixedInvMassBinnedMultPt[multBin]->Fill(pT, invMass);
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
  const int multBin = GetMultiplicityBin(lPercentile, fMultMode);
  if (multBin < 0) return;

  // Photon
  if (static_cast<int>(photonCandidates.size()) > 0) {
    if (static_cast<int>(fPhotonMixedBinned[multBin].size()) < fMixingDepth) {
      fPhotonMixedBinned[multBin].push_back(photonCandidates);
    } else {
      fPhotonMixedBinned[multBin].pop_front();
      fPhotonMixedBinned[multBin].push_back(photonCandidates);
    }
  }

  // ++++++++++++++
  // Lambda
  if (static_cast<int>(lambdaCandidates.size()) > 0) {
    if (static_cast<int>(fLambdaMixedBinned[multBin].size()) < fMixingDepth) {
      fLambdaMixedBinned[multBin].push_back(lambdaCandidates);
    } else {
      fLambdaMixedBinned[multBin].pop_front();
      fLambdaMixedBinned[multBin].push_back(lambdaCandidates);
    }
  }
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
  const int multBin = GetMultiplicityBin(lPercentile, fMultMode);

  // Loop over the MC tracks
  for (int iPart = 1; iPart < (fMCEvent->GetNumberOfTracks()); iPart++) {
    AliMCParticle *mcParticle =
        static_cast<AliMCParticle *>(fMCEvent->GetTrack(iPart));
    if (!mcParticle) continue;
    //    if (!mcParticle->IsPhysicalPrimary()) continue;
    if (mcParticle->GetNDaughters() != 2) continue;
    if (mcParticle->PdgCode() != fPDG) continue;

    if (TMath::Abs(mcParticle->Y()) <= fRapidityMax) {
      fHistMCTruthPt->Fill(mcParticle->Pt());
    }
    if (fIsSpectrumAnalysis) {
      if (multBin >= 0 && TMath::Abs(mcParticle->Y()) <= fRapidityMax) {
        fHistMCTruthPtMult[multBin]->Fill(mcParticle->Pt());
      }

      fHistMCTruthPtY->Fill(mcParticle->Y(), mcParticle->Pt());
      if (lPercentile < fMCHighMultThreshold) {
        fHistMCTruthPtYHighMult->Fill(mcParticle->Y(), mcParticle->Pt());
      }

      if (!CheckDaughters(mcParticle)) continue;
      fHistMCTruthDaughterPtY->Fill(mcParticle->Y(), mcParticle->Pt());
      if (lPercentile < fMCHighMultThreshold) {
        fHistMCTruthDaughterPtYHighMult->Fill(mcParticle->Y(),
                                              mcParticle->Pt());
      }

      if (!CheckDaughtersInAcceptance(mcParticle)) continue;
      fHistMCTruthDaughterPtYAccept->Fill(mcParticle->Y(), mcParticle->Pt());
      if (lPercentile < fMCHighMultThreshold) {
        fHistMCTruthDaughterPtYAcceptHighMult->Fill(mcParticle->Y(),
                                                    mcParticle->Pt());
      }
    }
  }
}

//____________________________________________________________________________________________________
bool AliSigma0PhotonMotherCuts::CheckDaughters(
    const AliMCParticle *particle) const {
  AliMCParticle *posDaughter = nullptr;
  AliMCParticle *negDaughter = nullptr;

  if (particle->GetNDaughters() != 2) return false;

  for (int daughterIndex = particle->GetDaughterFirst();
       daughterIndex <= particle->GetDaughterLast(); ++daughterIndex) {
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

  if (TMath::Abs(particle->Y()) > fRapidityMax) return false;

  if (particle->GetNDaughters() != 2) return false;

  for (int daughterIndex = particle->GetDaughterFirst();
       daughterIndex <= particle->GetDaughterLast(); ++daughterIndex) {
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
int AliSigma0PhotonMotherCuts::GetMultiplicityBin(float percentile,
                                                  UInt_t multMode) {
  if (multMode == AliVEvent::kINT7) {
    if (0 < percentile && percentile <= 5)
      return 0;
    else if (5. < percentile && percentile <= 15.)
      return 1;
    else if (15. < percentile && percentile <= 30.)
      return 2;
    else if (30. < percentile && percentile <= 50.)
      return 3;
    else if (50. < percentile && percentile <= 100.)
      return 4;
    else
      return -1;
  } else if (multMode == AliVEvent::kHighMultV0) {
    if (0 < percentile && percentile <= 0.01)
      return 0;
    else if (0.01 < percentile && percentile <= 0.05)
      return 1;
    else if (0.05 < percentile && percentile <= 0.1)
      return 2;
    else if (0.1 < percentile && percentile <= 1.)
      return 3;
    else if (1. < percentile && percentile <= 100.)
      return 4;
    else
      return -1;
  } else {
    std::cerr << "Multiplicity mode not defined \n";
    return -1;
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::InitCutHistograms(TString appendix) {
  fMassSigma = fDataBasePDG.GetParticle(fPDG)->Mass();

  std::cout << "============================\n"
            << " PHOTON MOTHER CUT CONFIGURATION \n"
            << " Sigma0 mass       " << fMassSigma << "\n"
            << " Sigma0 selection  " << fSigmaMassCut << "\n"
            << " Sigma0 sb up      " << fSidebandCutUp << "\n"
            << " Sigma0 sb down    " << fSidebandCutDown << "\n"
            << " Photon pT min     " << fPhotonPtMin << "\n"
            << " Photon pT max     " << fPhotonPtMax << "\n"
            << " Armenteros Qt low " << fArmenterosQtLow << "\n"
            << " Armenteros Qt up  " << fArmenterosQtUp << "\n"
            << " Armenteros a low  " << fArmenterosAlphaLow << "\n"
            << " Armenteros a up   " << fArmenterosAlphaUp << "\n"
            << " Rapidity max      " << fRapidityMax << "\n"
            << "============================\n";

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

  fHistCutBooking = new TProfile("fHistCutBooking", ";;Cut value", 16, 0, 16);
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
  fHistCutBooking->GetXaxis()->SetBinLabel(13,
                                           "Use M_{#Sigma^{0}}(#it{p}_{T})");
  fHistCutBooking->GetXaxis()->SetBinLabel(14,
                                           "p_{0} M_{#Sigma^{0}}(#it{p}_{T})");
  fHistCutBooking->GetXaxis()->SetBinLabel(15,
                                           "p_{1} M_{#Sigma^{0}}(#it{p}_{T})");
  fHistCutBooking->GetXaxis()->SetBinLabel(16,
                                           "p_{2} M_{#Sigma^{0}}(#it{p}_{T})");
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
  fHistCutBooking->Fill(12.f, static_cast<double>(fMassWindowPt));
  fHistCutBooking->Fill(13.f, fMassWindowP0);
  fHistCutBooking->Fill(14.f, fMassWindowP1);
  fHistCutBooking->Fill(15.f, fMassWindowP2);

  fHistInvMass =
      new TH1F("fHistInvMass", "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries",
               150, 1.15, 1.3);
  fHistograms->Add(fHistInvMass);

  fHistInvMassPt = new TH2F("fHistInvMassPt",
                            "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                            "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                            100, 0, 10, 300, 1., 1.3);
  fHistograms->Add(fHistInvMassPt);

  fHistInvMassPtRaw =
      new TH2F("fHistInvMassPtRaw",
               "before y cut; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
               "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
               100, 0, 10, 300, 1., 1.3);
  fHistograms->Add(fHistInvMassPtRaw);

  fHistNSigma =
      new TH1F("fHistNSigma", ";# #Sigma candidates; Entries", 10, 0, 10);
  fHistograms->Add(fHistNSigma);

  std::vector<float> multBinsLow, multBinsUp;
  if (fMultMode == AliVEvent::kINT7) {
    multBinsLow = {{0, 5., 15., 30., 50.}};
    multBinsUp = {{5., 15., 30., 50., 100.}};
  } else if (fMultMode == AliVEvent::kHighMultV0) {
    multBinsLow = {{0, 0.01, 0.05, 0.1, 1}};
    multBinsUp = {{0.01, 0.05, 0.1, 1, 100.}};
  }

  if (fIsSpectrumAnalysis) {
    fHistMixedInvMassPt = new TH2F("fHistMixedInvMassPt",
                                   "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                                   "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                                   100, 0, 10, 300, 1., 1.3);
    fHistograms->Add(fHistMixedInvMassPt);

    for (int i = 0; i < static_cast<int>(multBinsUp.size()); ++i) {
      fHistPtMult[i] = new TH2F(
          TString::Format("fHistPtMult_%i", i),
          TString::Format("V0M: %.2f - %.2f %%; #it{p}_{T} (GeV/#it{c}); "
                          "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                          multBinsLow[i], multBinsUp[i]),
          100, 0, 10, 300, 1.15, 1.3);
      fHistograms->Add(fHistPtMult[i]);

      fHistMixedInvMassBinnedMultPt[i] = new TH2F(
          TString::Format("fHistMixedInvMassBinnedMultPt_%i", i),
          TString::Format("V0M: %.2f - %.2f %%; #it{p}_{T} (GeV/#it{c}); "
                          "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                          multBinsLow[i], multBinsUp[i]),
          100, 0, 10, 300, 1.15, 1.3);
      fHistograms->Add(fHistMixedInvMassBinnedMultPt[i]);
    }
  }

  if (!fIsLightweight) {
    fHistNPhotonBefore =
        new TH1F("fHistNPhotonBefore",
                 ";# #gamma candidates before clean-up; Entries", 15, 0, 15);
    fHistNPhotonAfter = new TH1F(
        "fHistNPhotonAfter",
        ";# #gamma candidates after clean-up candidates; Entries", 15, 0, 15);
    fHistNLambdaBefore =
        new TH1F("fHistNLambdaBefore",
                 ";# #Lambda candidates before clean-up; Entries", 15, 0, 15);
    fHistNLambdaAfter =
        new TH1F("fHistNLambdaAfter",
                 ";# #Lambda candidates after clean-up; Entries", 15, 0, 15);
    fHistMassCutPt = new TH1F(
        "fHistMassCutPt", "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries",
        100, 0, 10);

    fHistInvMassK0Gamma = new TH1F("fHistInvMassK0Gamma",
                                   "; M_{K^{0}#gamma} (GeV/#it{c}); Entries",
                                   1000, 0.5, 3);

    fHistInvMassSelected = new TH2F("fHistInvMassSelected",
                                    "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                                    "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                                    50, 0, 10, 150, 1.15, 1.3);
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
    fHistArmenterosBefore =
        new TH2F("fHistArmenterosBefore", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 250, -1, 1, 250, 0, 0.5);
    fHistArmenterosAfter =
        new TH2F("fHistArmenterosAfter", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 250, -1, 1, 250, 0, 0.5);
    fHistEtaPhi = new TH2F("fHistEtaPhi", "; #eta; #phi", 100, -1, 1, 100,
                           -TMath::Pi(), TMath::Pi());
    fHistPtRapidity =
        new TH2F("fHistPtRapidity", "; #it{p}_{T} (GeV/#it{c}); y", 100, 0, 10,
                 100, -2, 2);

    fHistograms->Add(fHistInvMassK0Gamma);
    fHistograms->Add(fHistNPhotonBefore);
    fHistograms->Add(fHistNPhotonAfter);
    fHistograms->Add(fHistNLambdaBefore);
    fHistograms->Add(fHistNLambdaAfter);
    fHistograms->Add(fHistMassCutPt);
    fHistograms->Add(fHistInvMassSelected);
    fHistograms->Add(fHistInvMassRecPhoton);
    fHistograms->Add(fHistInvMassRecLambda);
    fHistograms->Add(fHistInvMassRec);
    fHistograms->Add(fHistEtaPhi);
    fHistograms->Add(fHistPtRapidity);
    fHistograms->Add(fHistArmenterosBefore);
    fHistograms->Add(fHistArmenterosAfter);

    fHistDeltaEtaDeltaPhiGammaNegBefore =
        new TH2F("fHistDeltaEtaDeltaPhiGammaNegBefore",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiGammaPosBefore =
        new TH2F("fHistDeltaEtaDeltaPhiGammaPosBefore",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiLambdaNegBefore =
        new TH2F("fHistDeltaEtaDeltaPhiLambdaNegBefore",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiLambdaPosBefore =
        new TH2F("fHistDeltaEtaDeltaPhiLambdaPosBefore",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiLambdaGammaNegBefore =
        new TH2F("fHistDeltaEtaDeltaPhiLambdaGammaNegBefore",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiLambdaGammaPosBefore =
        new TH2F("fHistDeltaEtaDeltaPhiLambdaGammaPosBefore",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiGammaNegAfter =
        new TH2F("fHistDeltaEtaDeltaPhiGammaNegAfter",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiGammaPosAfter =
        new TH2F("fHistDeltaEtaDeltaPhiGammaPosAfter",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiLambdaNegAfter =
        new TH2F("fHistDeltaEtaDeltaPhiLambdaNegAfter",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiLambdaPosAfter =
        new TH2F("fHistDeltaEtaDeltaPhiLambdaPosAfter",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiLambdaGammaNegAfter =
        new TH2F("fHistDeltaEtaDeltaPhiLambdaGammaNegAfter",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);
    fHistDeltaEtaDeltaPhiLambdaGammaPosAfter =
        new TH2F("fHistDeltaEtaDeltaPhiLambdaGammaPosAfter",
                 "; #Delta #eta; #Delta #phi*", 201, -0.1, 0.1, 201, -0.1, 0.1);

    fHistograms->Add(fHistDeltaEtaDeltaPhiGammaNegBefore);
    fHistograms->Add(fHistDeltaEtaDeltaPhiGammaPosBefore);
    fHistograms->Add(fHistDeltaEtaDeltaPhiLambdaNegBefore);
    fHistograms->Add(fHistDeltaEtaDeltaPhiLambdaPosBefore);
    fHistograms->Add(fHistDeltaEtaDeltaPhiLambdaGammaNegBefore);
    fHistograms->Add(fHistDeltaEtaDeltaPhiLambdaGammaPosBefore);
    fHistograms->Add(fHistDeltaEtaDeltaPhiGammaNegAfter);
    fHistograms->Add(fHistDeltaEtaDeltaPhiGammaPosAfter);
    fHistograms->Add(fHistDeltaEtaDeltaPhiLambdaNegAfter);
    fHistograms->Add(fHistDeltaEtaDeltaPhiLambdaPosAfter);
    fHistograms->Add(fHistDeltaEtaDeltaPhiLambdaGammaNegAfter);
    fHistograms->Add(fHistDeltaEtaDeltaPhiLambdaGammaPosAfter);

    fHistNPhotonLabel = new TH1F(
        "fHistNPhotonLabel",
        ";# #gamma candidates eliminated due label duplicates; Entries", 15, 0,
        15);
    fHistNLambdaLabel = new TH1F(
        "fHistNLambdaLabel",
        ";# #Lambda candidates eliminated due label duplicates; Entries", 15, 0,
        15);
    fHistNLambdaGammaLabel = new TH1F(
        "fHistNLambdaGammaLabel",
        ";#Candidates eliminated due #Lambda-#gamma label duplicates; Entries",
        15, 0, 15);
    fHistograms->Add(fHistNPhotonLabel);
    fHistograms->Add(fHistNLambdaLabel);
    fHistograms->Add(fHistNLambdaGammaLabel);

    fHistLambdaPtPhi =
        new TH2F("fHistLambdaPtPhi", "; #it{p}_{T} (GeV/#it{c}); #phi (rad)",
                 100, 0, 10, 100, 0, 2.f * TMath::Pi());
    fHistLambdaPtEta =
        new TH2F("fHistLambdaPtEta", "; #it{p}_{T} (GeV/#it{c}); #eta", 100, 0,
                 10, 100, -1, 1);
    fHistLambdaMassPt =
        new TH2F("fHistLambdaMassPt", "; #it{p}_{T} (GeV/#it{c}); M_{p#pi}",
                 100, 0, 10, 100, 1., 1.3);
    fHistPhotonPtPhi =
        new TH2F("fHistPhotonPtPhi", "; #it{p}_{T} (GeV/#it{c}); #phi (rad)",
                 100, 0, 10, 100, 0, 2.f * TMath::Pi());
    fHistPhotonPtEta =
        new TH2F("fHistPhotonPtEta", "; #it{p}_{T} (GeV/#it{c}); #eta", 100, 0,
                 10, 100, -1, 1);
    fHistPhotonMassPt =
        new TH2F("fHistPhotonMassPt", "; #it{p}_{T} (GeV/#it{c}); M_{p#pi}",
                 100, 0, 10, 100, 0., 0.1);

    fHistograms->Add(fHistLambdaPtPhi);
    fHistograms->Add(fHistLambdaPtEta);
    fHistograms->Add(fHistLambdaMassPt);
    fHistograms->Add(fHistPhotonPtPhi);
    fHistograms->Add(fHistPhotonPtEta);
    fHistograms->Add(fHistPhotonMassPt);

    fHistSigmaLambdaPtCorr = new TH2F("fHistSigmaLambdaPtCorr",
                                      "; #Sigma^{0} #it{p}_{T} (GeV/#it{c}; "
                                      "#Lambda #it{p}_{T} (GeV/#it{c}",
                                      100, 0, 10, 100, 0, 10);
    fHistSigmaPhotonPtCorr = new TH2F("fHistSigmaPhotonPtCorr",
                                      "; #Sigma^{0} #it{p}_{T} (GeV/#it{c}; "
                                      "#gamma #it{p}_{T} (GeV/#it{c}",
                                      100, 0, 10, 100, 0, 10);
    fHistSigmaLambdaPCorr = new TH2F("fHistSigmaLambdaPCorr",
                                     "; #Sigma^{0} #it{p} (GeV/#it{c}; "
                                     "#Lambda #it{p} (GeV/#it{c}",
                                     100, 0, 10, 100, 0, 10);
    fHistSigmaPhotonPCorr = new TH2F("fHistSigmaPhotonPCorr",
                                     "; #Sigma^{0} #it{p} (GeV/#it{c}; "
                                     "#gamma #it{p} (GeV/#it{c}",
                                     100, 0, 10, 100, 0, 10);

    fHistograms->Add(fHistSigmaLambdaPtCorr);
    fHistograms->Add(fHistSigmaPhotonPtCorr);
    fHistograms->Add(fHistSigmaLambdaPCorr);
    fHistograms->Add(fHistSigmaPhotonPCorr);
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

    fHistMCTruthPt = new TH1F("fHistMCTruthPt",
                              "; #it{p}_{T} (GeV/#it{c}); Entries", 100, 0, 10);
    fHistogramsMC->Add(fHistMCTruthPt);

    if (fIsSpectrumAnalysis) {
      for (int i = 0; i < static_cast<int>(multBinsUp.size()); ++i) {
        fHistMCTruthPtMult[i] = new TH1F(
            TString::Format("fHistMCTruthPtMult%i", i),
            TString::Format(
                "V0M: %.2f - %.2f %%; #it{p}_{T} (GeV/#it{c}); Entries",
                multBinsLow[i], multBinsUp[i]),
            100, 0, 10);
        fHistogramsMC->Add(fHistMCTruthPtMult[i]);
      }

      fHistMCTruthPtY =
          new TH2F("fHistMCTruthPtY", "; y; #it{p}_{T} (GeV/#it{c})", 100, -5,
                   5, 100, 0, 10);
      fHistMCTruthDaughterPtY =
          new TH2F("fHistMCTruthDaughterPtY", "; y; #it{p}_{T} (GeV/#it{c})",
                   200, -10, 10, 100, 0, 10);
      fHistMCTruthDaughterPtYAccept =
          new TH2F("fHistMCTruthDaughterPtYAccept",
                   "; y; #it{p}_{T} (GeV/#it{c})", 200, -10, 10, 100, 0, 10);
      fHistMCTruthPtYHighMult =
          new TH2F("fHistMCTruthPtYHighMult", "; y; #it{p}_{T} (GeV/#it{c})",
                   200, -10, 10, 100, 0, 10);
      fHistMCTruthDaughterPtYHighMult =
          new TH2F("fHistMCTruthDaughterPtYHighMult",
                   "; y; #it{p}_{T} (GeV/#it{c})", 200, -10, 10, 100, 0, 10);
      fHistMCTruthDaughterPtYAcceptHighMult =
          new TH2F("fHistMCTruthDaughterPtYAcceptHighMult",
                   "; y; #it{p}_{T} (GeV/#it{c})", 200, -10, 10, 100, 0, 10);

      fHistogramsMC->Add(fHistMCTruthPtY);
      fHistogramsMC->Add(fHistMCTruthDaughterPtY);
      fHistogramsMC->Add(fHistMCTruthDaughterPtYAccept);
      fHistogramsMC->Add(fHistMCTruthPtYHighMult);
      fHistogramsMC->Add(fHistMCTruthDaughterPtYHighMult);
      fHistogramsMC->Add(fHistMCTruthDaughterPtYAcceptHighMult);
    }

    fHistMCV0Pt = new TH1F("fHistMCV0Pt",
                           "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries",
                           100, 0, 10);
    fHistMCV0Mass =
        new TH1F("fHistMCV0Mass",
                 "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 500, 1., 1.5);
    fHistogramsMC->Add(fHistMCV0Pt);
    fHistogramsMC->Add(fHistMCV0Mass);

    if (!fIsLightweight) {
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

      fHistogramsMC->Add(fHistMCV0Mother);
      fHistogramsMC->Add(fHistMCV0Check);
      fHistogramsMC->Add(fHistMCV0MotherCheck);
    }

    fHistograms->Add(fHistogramsMC);
  }
}
