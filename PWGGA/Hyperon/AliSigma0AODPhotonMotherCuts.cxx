#include "AliSigma0AODPhotonMotherCuts.h"
#include <iostream>
#include "TMath.h"
#include "AliFemtoDreamv0.h"

ClassImp(AliSigma0AODPhotonMotherCuts)

    //____________________________________________________________________________________________________
    AliSigma0AODPhotonMotherCuts::AliSigma0AODPhotonMotherCuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fIsMC(false),
      fDoCleanUp(true),
      fIsLightweight(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fSigma(),
      fLambda(),
      fPhoton(),
      fSidebandUp(),
      fSidebandDown(),
      fLambdaCandidates(),
      fPhotonCandidates(),
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
      fDoDeltaPhiEtaCut(false),
      fDeltaPhiEtaMax(-1.f),
      fHistCutBooking(nullptr),
      fHistNSigma(nullptr),
      fHistNCandidates(nullptr),
      fHistNPhotonBefore(nullptr),
      fHistNPhotonAfter(nullptr),
      fHistNLambdaBefore(nullptr),
      fHistNLambdaAfter(nullptr),
      fHistNPhotonLabel(nullptr),
      fHistNLambdaLabel(nullptr),
      fHistNLambdaGammaLabel(nullptr),
      fHistNPhotonSplit(nullptr),
      fHistNLambdaSplit(nullptr),
      fHistNLambdaGammaSplit(nullptr),
      fHistMassCutPt(nullptr),
      fHistInvMass(nullptr),
      fHistInvMassOtherChilren(),
      fHistInvMassSelected(nullptr),
      fHistInvMassRecPhoton(nullptr),
      fHistInvMassRecLambda(nullptr),
      fHistInvMassRec(nullptr),
      fHistInvMassPtRaw(nullptr),
      fHistEtaPhi(nullptr),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
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
      fHistMCV0Pt(nullptr),
      fHistMCV0Mass(nullptr),
      fHistMCV0Mother(nullptr),
      fHistMCV0Check(nullptr),
      fHistMCV0MotherCheck(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0AODPhotonMotherCuts::AliSigma0AODPhotonMotherCuts(
    const AliSigma0AODPhotonMotherCuts &ref)
    : TObject(ref),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fIsMC(false),
      fDoCleanUp(true),
      fIsLightweight(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fSigma(),
      fLambda(),
      fPhoton(),
      fSidebandUp(),
      fSidebandDown(),
      fLambdaCandidates(),
      fPhotonCandidates(),
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
      fDoDeltaPhiEtaCut(false),
      fDeltaPhiEtaMax(-1.f),
      fHistCutBooking(nullptr),
      fHistNSigma(nullptr),
      fHistNCandidates(nullptr),
      fHistNPhotonBefore(nullptr),
      fHistNPhotonAfter(nullptr),
      fHistNLambdaBefore(nullptr),
      fHistNLambdaAfter(nullptr),
      fHistNPhotonLabel(nullptr),
      fHistNLambdaLabel(nullptr),
      fHistNLambdaGammaLabel(nullptr),
      fHistNPhotonSplit(nullptr),
      fHistNLambdaSplit(nullptr),
      fHistNLambdaGammaSplit(nullptr),
      fHistMassCutPt(nullptr),
      fHistInvMass(nullptr),
      fHistInvMassOtherChilren(),
      fHistInvMassSelected(nullptr),
      fHistInvMassRecPhoton(nullptr),
      fHistInvMassRecLambda(nullptr),
      fHistInvMassRec(nullptr),
      fHistInvMassPtRaw(nullptr),
      fHistEtaPhi(nullptr),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
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
      fHistMCV0Pt(nullptr),
      fHistMCV0Mass(nullptr),
      fHistMCV0Mother(nullptr),
      fHistMCV0Check(nullptr),
      fHistMCV0MotherCheck(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0AODPhotonMotherCuts &AliSigma0AODPhotonMotherCuts::operator=(
    const AliSigma0AODPhotonMotherCuts &ref) {
  // Assignment operator
  if (this == &ref) return *this;

  return (*this);
}

//____________________________________________________________________________________________________
AliSigma0AODPhotonMotherCuts *AliSigma0AODPhotonMotherCuts::DefaultCuts() {
  AliSigma0AODPhotonMotherCuts *photonMotherCuts =
      new AliSigma0AODPhotonMotherCuts();
  photonMotherCuts->SetSigmaMassPt(true);
  photonMotherCuts->SetMinPt(1.0);
  return photonMotherCuts;
}

//____________________________________________________________________________________________________
void AliSigma0AODPhotonMotherCuts::SelectPhotonMother(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    const std::vector<AliFemtoDreamBasePart> &photonCandidates,
    const std::vector<AliFemtoDreamBasePart> &lambdaCandidates) {
  fInputEvent = inputEvent;
  fMCEvent = mcEvent;

  fSigma.clear();
  fSidebandUp.clear();
  fSidebandDown.clear();
  fLambda.clear();
  fPhoton.clear();
  fLambdaCandidates.clear();
  fPhotonCandidates.clear();

  if(photonCandidates.size() == 0 || lambdaCandidates.size() == 0) {
    return;
  }

  fLambdaCandidates.assign( lambdaCandidates.begin(), lambdaCandidates.end() );
  fPhotonCandidates.assign( photonCandidates.begin(), photonCandidates.end() );

  if (fDoCleanUp) {
    CleanUpClones();
  }

  if (!fIsLightweight) SingleV0QA();

  // Particle pairing
  SigmaToLambdaGamma();
}

//____________________________________________________________________________________________________
void AliSigma0AODPhotonMotherCuts::CleanUpClones() {
  // Keep track of the fantastic situation before applying the cuts

  if (!fIsLightweight) {
    fHistNPhotonBefore->Fill(fPhotonCandidates.size());
    fHistNLambdaBefore->Fill(fLambdaCandidates.size());
  }

  // Do the checks for the photons
  int nPhotonKilledLabel = 0;
  int nPhotonKilledSplit = 0;
  for (auto photon1 = fPhotonCandidates.begin();
       photon1 < fPhotonCandidates.end(); ++photon1) {
    if (!photon1->UseParticle()) continue;
    std::vector<int> photon1TrackLabels = photon1->GetIDTracks();
    std::vector<float> photon1Eta = photon1->GetEta();
    float photon1PosPhiStar = photon1->GetAveragePhiAtRadius(0);
    float photon1NegPhiStar = photon1->GetAveragePhiAtRadius(1);

    for (auto photon2 = photon1 + 1; photon2 < fPhotonCandidates.end();
         ++photon2) {
      if (!photon1->UseParticle() || !photon2->UseParticle()) continue;
      std::vector<int> photon2TrackLabels = photon2->GetIDTracks();
      std::vector<float> photon2Eta = photon2->GetEta();
      float photon2PosPhiStar = photon2->GetAveragePhiAtRadius(0);
      float photon2NegPhiStar = photon2->GetAveragePhiAtRadius(1);

      float deltaPhistarPos = photon1PosPhiStar - photon2PosPhiStar;
      float deltaPhistarNeg = photon1NegPhiStar - photon2NegPhiStar;
      float deltaEtaPos = photon1Eta[1] - photon2Eta[1];
      float deltaEtaNeg = photon1Eta[2] - photon2Eta[2];

      if (!fIsLightweight) {
        fHistDeltaEtaDeltaPhiGammaNegBefore->Fill(deltaEtaNeg, deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiGammaPosBefore->Fill(deltaEtaPos, deltaPhistarPos);
      }

      bool hasSameLabels = (photon1TrackLabels[0] == photon2TrackLabels[0] ||
                            photon1TrackLabels[1] == photon2TrackLabels[0] ||
                            photon1TrackLabels[1] == photon2TrackLabels[1] ||
                            photon1TrackLabels[0] == photon2TrackLabels[1]);
      bool hasTrackSplitting =
          (fDoDeltaPhiEtaCut
              && (((deltaPhistarPos * deltaPhistarPos
                  + deltaEtaPos * deltaEtaPos) < fDeltaPhiEtaMax)
                  || ((deltaPhistarNeg * deltaPhistarNeg
                      + deltaEtaNeg * deltaEtaNeg) < fDeltaPhiEtaMax)));

      // do the check for both daughters
      if (hasSameLabels) ++nPhotonKilledLabel;
      if (hasTrackSplitting) ++nPhotonKilledSplit;
      if (hasSameLabels || hasTrackSplitting) {
        const float cpa1 = photon1->GetCPA();
        const float cpa2 = photon2->GetCPA();
        if (cpa1 > cpa2) {
          photon2->SetUse(false);
        } else {
          photon1->SetUse(false);
        }
      }
      if (!photon1->UseParticle()) break;
    }
  }

  // Do the checks for the Lambdas
  int nLambdaKilledLabel = 0;
  int nLambdaKilledSplit = 0;
  for (auto lambda1 = fLambdaCandidates.begin();
       lambda1 < fLambdaCandidates.end(); ++lambda1) {
    if (!lambda1->UseParticle()) continue;
    std::vector<int> lambda1TrackLabels = lambda1->GetIDTracks();
    std::vector<float> lambda1Eta = lambda1->GetEta();
    float lambda1PosPhiStar = lambda1->GetAveragePhiAtRadius(0);
    float lambda1NegPhiStar = lambda1->GetAveragePhiAtRadius(1);

    for (auto lambda2 = lambda1 + 1; lambda2 < fLambdaCandidates.end();
         ++lambda2) {
      if (!lambda1->UseParticle() || !lambda2->UseParticle()) continue;
      std::vector<int> lambda2TrackLabels = lambda2->GetIDTracks();
      std::vector<float> lambda2Eta = lambda2->GetEta();
      float lambda2PosPhiStar = lambda2->GetAveragePhiAtRadius(0);
      float lambda2NegPhiStar = lambda2->GetAveragePhiAtRadius(1);

      float deltaPhistarPos = lambda1PosPhiStar - lambda2PosPhiStar;
      float deltaPhistarNeg = lambda1NegPhiStar - lambda2NegPhiStar;
      float deltaEtaPos = lambda1Eta[1] - lambda2Eta[1];
      float deltaEtaNeg = lambda1Eta[2] - lambda2Eta[2];

      if (!fIsLightweight) {
        fHistDeltaEtaDeltaPhiLambdaNegBefore->Fill(deltaEtaNeg,
                                                   deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiLambdaPosBefore->Fill(deltaEtaPos,
                                                   deltaPhistarPos);
      }

      bool hasSameLabels = (lambda1TrackLabels[0] == lambda2TrackLabels[0] ||
                            lambda1TrackLabels[1] == lambda2TrackLabels[0] ||
                            lambda1TrackLabels[1] == lambda2TrackLabels[1] ||
                            lambda1TrackLabels[0] == lambda2TrackLabels[1]);
      bool hasTrackSplitting =
          (fDoDeltaPhiEtaCut
              && (((deltaPhistarPos * deltaPhistarPos
                  + deltaEtaPos * deltaEtaPos) < fDeltaPhiEtaMax)
                  || ((deltaPhistarNeg * deltaPhistarNeg
                      + deltaEtaNeg * deltaEtaNeg) < fDeltaPhiEtaMax)));

      // do the check for both daughters
      if (hasSameLabels) ++nLambdaKilledLabel;
      if (hasTrackSplitting) ++nLambdaKilledSplit;
      if (hasSameLabels || hasTrackSplitting) {
        const float cpa1 = lambda1->GetCPA();
        const float cpa2 = lambda2->GetCPA();
        if (cpa1 > cpa2) {
          lambda2->SetUse(false);
        } else {
          lambda1->SetUse(false);
        }
      }
      if (!lambda1->UseParticle()) break;
    }
  }

  // Now do the exercise for photon-Lambda combination
  int nPhotonLambdaKilledLabel = 0;
  int nPhotonLambdaKilledSplit = 0;
  for (auto photon = fPhotonCandidates.begin(); photon < fPhotonCandidates.end();
       ++photon) {
    if (!photon->UseParticle()) continue;
    std::vector<int> photonTrackLabels = photon->GetIDTracks();
    std::vector<float> photonEta = photon->GetEta();
    float photonPosPhiStar = photon->GetAveragePhiAtRadius(0);
    float photonNegPhiStar = photon->GetAveragePhiAtRadius(1);

    for (auto lambda = fLambdaCandidates.begin();
         lambda < fLambdaCandidates.end(); ++lambda) {
      if (!lambda->UseParticle() || !photon->UseParticle()) continue;
      std::vector<int> lambdaTrackLabels = lambda->GetIDTracks();
      std::vector<float> lambdaEta = lambda->GetEta();
      float lambdaPosPhiStar = lambda->GetAveragePhiAtRadius(0);
      float lambdaNegPhiStar = lambda->GetAveragePhiAtRadius(1);

      float deltaPhistarPos = photonPosPhiStar - lambdaPosPhiStar;
      float deltaPhistarNeg = photonNegPhiStar - lambdaNegPhiStar;
      float deltaEtaPos = photonEta[1] - lambdaEta[1];
      float deltaEtaNeg = photonEta[2] - lambdaEta[2];

      if (!fIsLightweight) {
        fHistDeltaEtaDeltaPhiLambdaGammaNegBefore->Fill(deltaEtaNeg,
                                                        deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiLambdaGammaPosBefore->Fill(deltaEtaPos,
                                                        deltaPhistarPos);
      }

      bool hasSameLabels = (photonTrackLabels[0] == lambdaTrackLabels[0] ||
                            photonTrackLabels[1] == lambdaTrackLabels[0] ||
                            photonTrackLabels[1] == lambdaTrackLabels[1] ||
                            photonTrackLabels[0] == lambdaTrackLabels[1]);
      bool hasTrackSplitting =
          (fDoDeltaPhiEtaCut
              && (((deltaPhistarPos * deltaPhistarPos
                  + deltaEtaPos * deltaEtaPos) < fDeltaPhiEtaMax)
                  || ((deltaPhistarNeg * deltaPhistarNeg
                      + deltaEtaNeg * deltaEtaNeg) < fDeltaPhiEtaMax)));

      // do the check for both daughters
      if (hasSameLabels) ++nPhotonLambdaKilledLabel;
      if (hasTrackSplitting) ++nPhotonLambdaKilledSplit;
      if (hasSameLabels || hasTrackSplitting) {
        const float cpaPhoton = photon->GetCPA();
        const float cpaLambda = lambda->GetCPA();
        if (cpaPhoton > cpaLambda) {
          lambda->SetUse(false);
        } else {
          photon->SetUse(false);
        }
      }
      if (!photon->UseParticle()) break;
    }
  }

  // We're done here, now comes the QA

  // Now take a look what changed for the Photons
  int nPhotonAfter = 0;
  for (auto photon1 = fPhotonCandidates.begin();
       photon1 < fPhotonCandidates.end(); ++photon1) {
    if (!photon1->UseParticle()) continue;
    std::vector<int> photon1TrackLabels = photon1->GetIDTracks();
    std::vector<float> photon1Eta = photon1->GetEta();
    float photon1PosPhiStar = photon1->GetAveragePhiAtRadius(0);
    float photon1NegPhiStar = photon1->GetAveragePhiAtRadius(1);
    ++nPhotonAfter;

    for (auto photon2 = photon1 + 1; photon2 < fPhotonCandidates.end();
         ++photon2) {
      if (!photon2->UseParticle()) continue;
      std::vector<int> photon2TrackLabels = photon2->GetIDTracks();
      std::vector<float> photon2Eta = photon2->GetEta();
      float photon2PosPhiStar = photon2->GetAveragePhiAtRadius(0);
      float photon2NegPhiStar = photon2->GetAveragePhiAtRadius(1);

      if (!fIsLightweight) {
        float deltaPhistarPos = photon1PosPhiStar - photon2PosPhiStar;
        float deltaPhistarNeg = photon1NegPhiStar - photon2NegPhiStar;
        float deltaEtaPos = photon1Eta[1] - photon2Eta[1];
        float deltaEtaNeg = photon1Eta[2] - photon2Eta[2];
        fHistDeltaEtaDeltaPhiGammaNegAfter->Fill(deltaEtaNeg, deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiGammaPosAfter->Fill(deltaEtaPos, deltaPhistarPos);
      }
    }
  }

  // Now take a look what changed for the Lambdas
  int nLambdaAfter = 0;
  for (auto photon = fPhotonCandidates.begin(); photon < fPhotonCandidates.end();
       ++photon) {
    if (!photon->UseParticle()) continue;
    std::vector<int> photonTrackLabels = photon->GetIDTracks();
    std::vector<float> photonEta = photon->GetEta();
    float photonPosPhiStar = photon->GetAveragePhiAtRadius(0);
    float photonNegPhiStar = photon->GetAveragePhiAtRadius(1);

    for (auto lambda = fLambdaCandidates.begin();
         lambda < fLambdaCandidates.end(); ++lambda) {
      if (!lambda->UseParticle()) continue;
      std::vector<int> lambdaTrackLabels = lambda->GetIDTracks();
      std::vector<float> lambdaEta = lambda->GetEta();
      float lambdaPosPhiStar = lambda->GetAveragePhiAtRadius(0);
      float lambdaNegPhiStar = lambda->GetAveragePhiAtRadius(1);

      if (!fIsLightweight) {
        float deltaPhistarPos = photonPosPhiStar - lambdaPosPhiStar;
        float deltaPhistarNeg = photonNegPhiStar - lambdaNegPhiStar;
        float deltaEtaPos = photonEta[1] - lambdaEta[1];
        float deltaEtaNeg = photonEta[2] - lambdaEta[2];
        fHistDeltaEtaDeltaPhiLambdaGammaNegAfter->Fill(deltaEtaNeg,
                                                       deltaPhistarNeg);
        fHistDeltaEtaDeltaPhiLambdaGammaPosAfter->Fill(deltaEtaPos,
                                                       deltaPhistarPos);
      }
    }
  }

  for (auto lambda1 = fLambdaCandidates.begin();
       lambda1 < fLambdaCandidates.end(); ++lambda1) {
    if (!lambda1->UseParticle()) continue;
    ++nLambdaAfter;
    std::vector<int> lambda1TrackLabels = lambda1->GetIDTracks();
    std::vector<float> lambda1Eta = lambda1->GetEta();
    float lambda1PosPhiStar = lambda1->GetAveragePhiAtRadius(0);
    float lambda1NegPhiStar = lambda1->GetAveragePhiAtRadius(1);

    for (auto lambda2 = lambda1 + 1; lambda2 < fLambdaCandidates.end();
         ++lambda2) {
      if (!lambda2->UseParticle()) continue;
      std::vector<int> lambda2TrackLabels = lambda2->GetIDTracks();
      std::vector<float> lambda2Eta = lambda2->GetEta();
      float lambda2PosPhiStar = lambda2->GetAveragePhiAtRadius(0);
      float lambda2NegPhiStar = lambda2->GetAveragePhiAtRadius(1);

      if (!fIsLightweight) {
        float deltaPhistarPos = lambda1PosPhiStar - lambda2PosPhiStar;
        float deltaPhistarNeg = lambda1NegPhiStar - lambda2NegPhiStar;
        float deltaEtaPos = lambda1Eta[1] - lambda2Eta[1];
        float deltaEtaNeg = lambda1Eta[2] - lambda2Eta[2];
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
    fHistNPhotonSplit->Fill(nPhotonKilledSplit);
    fHistNLambdaSplit->Fill(nLambdaKilledSplit);
    fHistNLambdaGammaSplit->Fill(nPhotonLambdaKilledSplit);
  }
}

//____________________________________________________________________________________________________
void AliSigma0AODPhotonMotherCuts::SingleV0QA() {
  // Photon QA - keep track of kinematics
  for (const auto &photon : fPhotonCandidates) {
    if (!photon.UseParticle()) continue;
    fHistPhotonPtPhi->Fill(photon.GetPt(), photon.GetPhi()[0]);
    fHistPhotonPtEta->Fill(photon.GetPt(), photon.GetEta()[0]);
    fHistPhotonMassPt->Fill(photon.GetPt(), photon.GetInvMass());
  }

  // Lambda QA - keep track of kinematics
  for (const auto &lambda : fLambdaCandidates) {
    if (!lambda.UseParticle()) continue;
    fHistLambdaPtPhi->Fill(lambda.GetPt(), lambda.GetPhi()[0]);
    fHistLambdaPtEta->Fill(lambda.GetPt(), lambda.GetEta()[0]);
    fHistLambdaMassPt->Fill(lambda.GetPt(), lambda.GetInvMass());
  }
}

//____________________________________________________________________________________________________
void AliSigma0AODPhotonMotherCuts::SigmaToLambdaGamma() {
  const float massProton   = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const float massPion     = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const float massElectron = TDatabasePDG::Instance()->GetParticle(11)->Mass();

  const float massK0 = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  const float massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const float massAntiLambda = TDatabasePDG::Instance()->GetParticle(-3122)->Mass();

  std::vector<float> massHypD1 = {{ massProton, massPion, massPion, massElectron}};
  std::vector<float> massHypD2 = {{ massPion, massProton, massPion, massElectron}};
  std::vector<float> massHypTrack = {{massLambda, massAntiLambda, massK0, 0.f}};

  auto *dummy = static_cast<AliAODEvent*>(fInputEvent);

  int nSigma = 0;
  const float lambdaMass = fDataBasePDG.GetParticle(fPDGDaughter1)->Mass();

  for (const auto &lambda : fLambdaCandidates) {
    fHistNCandidates->Fill(0);
    if (!lambda.UseParticle()) continue;
    fHistNCandidates->Fill(1);
  }

  for (const auto &photon : fPhotonCandidates) {
    fHistNCandidates->Fill(2);
    if (!photon.UseParticle()) continue;
    fHistNCandidates->Fill(3);
  }

  // SAME EVENT
  AliFemtoDreamv0 sigma(7);
  for (const auto &photon : fPhotonCandidates) {
    if (photon.GetPt() > fPhotonPtMax || photon.GetPt() < fPhotonPtMin)
      continue;
    if (!photon.UseParticle()) continue;

    for (const auto &lambda : fLambdaCandidates) {
      if (!lambda.UseParticle()) continue;
      sigma.Setv0(lambda, photon, dummy, true, true, false);
      sigma.Setv0Mass(sigma.GetInvMass() - lambda.GetInvMass() + lambdaMass -
                      photon.GetInvMass());
      const float invMass = sigma.GetInvMass();
      const float armAlpha = GetArmenterosAlpha(photon, lambda, sigma);
      const float armQt = GetArmenterosQt(photon, lambda, sigma);
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
      if (!fIsLightweight) {
        fHistArmenterosAfter->Fill(armAlpha, armQt);
      }

      if (fIsMC) {
        // Get the mother of the particles
        AliAODMCParticle *partV0 =
            static_cast<AliAODMCParticle *>(fMCEvent->GetTrack(
                lambda.GetID()));
        // we need the first mother of the Lambda, not after the full cascade to match this to a Sigma0
        AliAODMCParticle *partMotherV0 = (partV0) ?
            static_cast<AliAODMCParticle *>(fMCEvent->GetTrack(
                partV0->GetMother())) : nullptr;
        AliAODMCParticle *partMotherPhoton =
            static_cast<AliAODMCParticle *>(fMCEvent->GetTrack(
                photon.GetMotherID()));

        fHistMCV0Check->Fill(TMath::Abs(lambda.GetMCPDGCode()),
                             TMath::Abs(photon.GetMCPDGCode()));

        if (partMotherV0 && partMotherPhoton) {
          fHistMCV0MotherCheck->Fill(
              TMath::Abs(partMotherV0->GetPdgCode()),
              TMath::Abs(partMotherPhoton->GetPdgCode()));

          if (partMotherV0 == partMotherPhoton && TMath::Abs(partMotherV0->GetPdgCode() == fPDG)) {
            sigma.SetMCParticle(partMotherV0, fMCEvent);
            sigma.SetMCPDGCode(fPDG);
            fHistMCV0Mother->Fill(invMass, TMath::Abs(partMotherV0->GetPdgCode()));

            fHistMCV0Pt->Fill(pT);
            fHistMCV0Mass->Fill(invMass);
          }
        }
      }

      // Now write out the stuff to the Femto containers
      if (invMass < GetMassSigmaPt(pT) + fSigmaMassCut &&
          invMass > GetMassSigmaPt(pT) - fSigmaMassCut) {
        fHistNCandidates->Fill(4);

        fSigma.push_back(sigma);
        fLambda.push_back(lambda);
        fPhoton.push_back(photon);
        if (!fIsLightweight) {
          fHistInvMassSelected->Fill(pT, invMass);
          fHistMassCutPt->Fill(pT);
          fHistEtaPhi->Fill(sigma.GetEta()[0], sigma.GetPhi()[0]);
          fHistSigmaLambdaPtCorr->Fill(pT, lambda.GetPt());
          fHistSigmaPhotonPtCorr->Fill(pT, photon.GetPt());
          fHistSigmaLambdaPCorr->Fill(sigma.GetP(), lambda.GetP());
          fHistSigmaPhotonPCorr->Fill(sigma.GetP(), photon.GetP());

          TLorentzVector track1D1, track1D2, track2D1, track2D2, track1, track2,
              trackSum;
          int count = 0;
          for (size_t i = 0; i < massHypTrack.size(); ++i) {
            track1D1.SetXYZM(lambda.GetMomentum(1).Px(),
                             lambda.GetMomentum(1).Py(),
                             lambda.GetMomentum(1).Pz(), massHypD1[i]);
            track1D2.SetXYZM(lambda.GetMomentum(2).Px(),
                             lambda.GetMomentum(2).Py(),
                             lambda.GetMomentum(2).Pz(), massHypD2[i]);
            track1 = track1D1 + track1D2;

            for (size_t j = 0; j < massHypTrack.size(); ++j) {
              track2D1.SetXYZM(photon.GetMomentum(1).Px(),
                               photon.GetMomentum(1).Py(),
                               photon.GetMomentum(1).Pz(), massHypD1[j]);
              track2D2.SetXYZM(photon.GetMomentum(2).Px(),
                               photon.GetMomentum(2).Py(),
                               photon.GetMomentum(2).Pz(), massHypD2[j]);
              track2 = track2D1 + track2D2;
              TLorentzVector trackSum = track1 + track2;
              fHistInvMassOtherChilren[count++]->Fill(trackSum.M());
            }
          }
        }
        ++nSigma;
      }
      if (invMass < GetMassSigmaPt(pT) + fSidebandCutUp &&
          invMass > GetMassSigmaPt(pT) + fSidebandCutDown) {
        fHistNCandidates->Fill(5);
        fSidebandUp.push_back(sigma);
        if (!fIsLightweight) {
          fHistInvMassSelected->Fill(pT, invMass);
        }
      }
      if (invMass > GetMassSigmaPt(pT) - fSidebandCutUp &&
          invMass < GetMassSigmaPt(pT) - fSidebandCutDown) {
        fHistNCandidates->Fill(6);
        fSidebandDown.push_back(sigma);
        if (!fIsLightweight) {
          fHistInvMassSelected->Fill(pT, invMass);
        }
      }

      fHistInvMass->Fill(invMass);
      fHistInvMassPtRaw->Fill(pT, invMass);
    }
  }
  fHistNSigma->Fill(nSigma);
}

//____________________________________________________________________________________________________
float AliSigma0AODPhotonMotherCuts::GetMassSigmaPt(float pt) const {
  float massSigma = -1.f;
  if (fMassWindowPt) {
    massSigma = fMassWindowP0 + fMassWindowP1 * std::exp(fMassWindowP2 * pt);
  } else {
    massSigma = fMassSigma;
  }
  return massSigma;
}

//____________________________________________________________________________________________________
float AliSigma0AODPhotonMotherCuts::GetArmenterosAlpha(
    const AliFemtoDreamBasePart &gamma, const AliFemtoDreamBasePart &lambda,
    const AliFemtoDreamBasePart &sigma) const {
  TVector3 gammaP = gamma.GetMomentum();
  TVector3 lambdaP = lambda.GetMomentum();
  TVector3 sigmaP = sigma.GetMomentum();
  return 1. - 2. / (1. + gammaP.Dot(sigmaP) / lambdaP.Dot(sigmaP));
}

//____________________________________________________________________________________________________
float AliSigma0AODPhotonMotherCuts::GetArmenterosQt(
    const AliFemtoDreamBasePart &gamma, const AliFemtoDreamBasePart &lambda,
    const AliFemtoDreamBasePart &sigma) const {
  TVector3 lambdaP = lambda.GetMomentum();
  TVector3 sigmaP = sigma.GetMomentum();
  return lambdaP.Perp(sigmaP);
}

//____________________________________________________________________________________________________
void AliSigma0AODPhotonMotherCuts::InitCutHistograms(TString appendix) {
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

  fHistCutBooking = new TProfile("fHistCutBooking", ";;Cut value", 14, 0, 14);
  fHistCutBooking->GetXaxis()->SetBinLabel(1, "#Sigma^{0} selection");
  fHistCutBooking->GetXaxis()->SetBinLabel(2, "#Sigma^{0} sb down");
  fHistCutBooking->GetXaxis()->SetBinLabel(3, "#Sigma^{0} sb up");
  fHistCutBooking->GetXaxis()->SetBinLabel(4, "#gamma #it{p}_{T} min");
  fHistCutBooking->GetXaxis()->SetBinLabel(5, "#gamma #it{p}_{T} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(6, "Armenteros q_{T} low");
  fHistCutBooking->GetXaxis()->SetBinLabel(7, "Armenteros q_{T} up");
  fHistCutBooking->GetXaxis()->SetBinLabel(8, "Armenteros #alpha low");
  fHistCutBooking->GetXaxis()->SetBinLabel(9, "Armenteros #alpha up");
  fHistCutBooking->GetXaxis()->SetBinLabel(10, "Rapidity y max");
  fHistCutBooking->GetXaxis()->SetBinLabel(11,
                                           "Use M_{#Sigma^{0}}(#it{p}_{T})");
  fHistCutBooking->GetXaxis()->SetBinLabel(12,
                                           "p_{0} M_{#Sigma^{0}}(#it{p}_{T})");
  fHistCutBooking->GetXaxis()->SetBinLabel(13,
                                           "p_{1} M_{#Sigma^{0}}(#it{p}_{T})");
  fHistCutBooking->GetXaxis()->SetBinLabel(14,
                                           "p_{2} M_{#Sigma^{0}}(#it{p}_{T})");
  fHistograms->Add(fHistCutBooking);

  fHistCutBooking->Fill(0.f, fSigmaMassCut);
  fHistCutBooking->Fill(1.f, fSidebandCutDown);
  fHistCutBooking->Fill(2.f, fSidebandCutUp);
  fHistCutBooking->Fill(3.f, fPhotonPtMin);
  fHistCutBooking->Fill(4.f, fPhotonPtMax);
  fHistCutBooking->Fill(5.f, fArmenterosQtLow);
  fHistCutBooking->Fill(6.f, fArmenterosQtUp);
  fHistCutBooking->Fill(7.f, fArmenterosAlphaLow);
  fHistCutBooking->Fill(8.f, fArmenterosAlphaUp);
  fHistCutBooking->Fill(9.f, fRapidityMax);
  fHistCutBooking->Fill(10.f, static_cast<double>(fMassWindowPt));
  fHistCutBooking->Fill(11.f, fMassWindowP0);
  fHistCutBooking->Fill(12.f, fMassWindowP1);
  fHistCutBooking->Fill(13.f, fMassWindowP2);

  fHistInvMass =
      new TH1F("fHistInvMass", "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries",
               600, 1.15, 1.75);
  fHistograms->Add(fHistInvMass);

  fHistInvMassPtRaw =
      new TH2F("fHistInvMassPtRaw",
               "before y cut; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
               "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
               100, 0, 10, 300, 1., 1.3);
  fHistograms->Add(fHistInvMassPtRaw);

  fHistNSigma =
      new TH1F("fHistNSigma", ";# #Sigma candidates; Entries", 10, 0, 10);
  fHistograms->Add(fHistNSigma);

  fHistNCandidates = new TH1F("fHistNCandidates", ";;Counts", 10, 0, 10);
  fHistNCandidates->GetXaxis()->SetBinLabel(1, "#Lambda in");
  fHistNCandidates->GetXaxis()->SetBinLabel(2, "#Lambda ok");
  fHistNCandidates->GetXaxis()->SetBinLabel(3, "#gamma in");
  fHistNCandidates->GetXaxis()->SetBinLabel(4, "#gamma ok");
  fHistNCandidates->GetXaxis()->SetBinLabel(5, "#Sigma^{0}");
  fHistNCandidates->GetXaxis()->SetBinLabel(6, "SB up");
  fHistNCandidates->GetXaxis()->SetBinLabel(7, "SB low");
  fHistograms->Add(fHistNCandidates);

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

    std::vector<TString> partNames = {{ "#Lambda", "#bar{#Lambda}", "K^{0}_{S}", "#gamma"}};
    int count = 0;
    for (size_t i = 0; i < partNames.size(); ++i) {
      for (size_t j = 0; j < partNames.size(); ++j) {
        fHistInvMassOtherChilren[count] = new TH1F(
            Form("fHistInvMassOtherChilren_%i", count),
            Form("; M_{%s%s} (GeV/#it{c}); Entries", partNames[i].Data(),
                 partNames[j].Data()),
            5000, 0., 10);
        fHistograms->Add(fHistInvMassOtherChilren[count]);
        count++;
      }
    }

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

    fHistNPhotonSplit= new TH1F(
        "fHistNPhotonSplit",
        ";# #gamma candidates eliminated due to #Delta#eta#Delta#phi*; Entries", 15, 0,
        15);
    fHistNLambdaSplit = new TH1F(
        "fHistNLambdaSplit",
        ";# #Lambda candidates eliminated due to #Delta#eta#Delta#phi*; Entries", 15, 0,
        15);
    fHistNLambdaGammaSplit = new TH1F(
        "fHistNLambdaGammaSplit",
        ";#Candidates eliminated due to #Lambda-#gamma #Delta#eta#Delta#phi*; Entries",
        15, 0, 15);
    fHistograms->Add(fHistNPhotonSplit);
    fHistograms->Add(fHistNLambdaSplit);
    fHistograms->Add(fHistNLambdaGammaSplit);

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
