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
      fMassSigma(0),
      fSigmaMassCut(0.005),
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
      fHistInvMassRecPhoton(nullptr),
      fHistInvMassRecLambda(nullptr),
      fHistInvMassRec(nullptr),
      fHistInvMassPt(nullptr),
      fHistEtaPhi(nullptr),
      fHistPtRapidity(nullptr),
      fHistPtMult(),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistMixedInvMassPt(nullptr),
      fHistMixedInvMassBinnedMultPt(),
	  fHistDiffPGamma(),
	  fHistDiffPLambda(),
	  fHistDiffPGammaPos(),
	  fHistDiffPLambdaPos(),
	  fHistDiffPGammaNeg(),
	  fHistDiffPLambdaNeg(),
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
      fHistMCTrueSigmaLambdaPtCorr(nullptr),
      fHistMCTrueSigmaPhotonPtCorr(nullptr),
      fHistMCTrueSigmaLambdaPCorr(nullptr),
      fHistMCTrueSigmaPhotonPCorr(nullptr),
      fHistMCBkgSigmaLambdaPtCorr(nullptr),
      fHistMCBkgSigmaPhotonPtCorr(nullptr),
      fHistMCBkgSigmaLambdaPCorr(nullptr),
      fHistMCBkgSigmaPhotonPCorr(nullptr),
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
      fMassSigma(0),
      fSigmaMassCut(0.005),
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
      fHistInvMassRecPhoton(nullptr),
      fHistInvMassRecLambda(nullptr),
      fHistInvMassRec(nullptr),
      fHistInvMassPt(nullptr),
      fHistEtaPhi(nullptr),
      fHistPtRapidity(nullptr),
      fHistPtMult(),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistMixedInvMassPt(nullptr),
      fHistMixedInvMassBinnedMultPt(),
	  fHistDiffPGamma(),
	  fHistDiffPLambda(),
	  fHistDiffPGammaPos(),
	  fHistDiffPLambdaPos(),
	  fHistDiffPGammaNeg(),
	  fHistDiffPLambdaNeg(),
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
      fHistMCTrueSigmaLambdaPtCorr(nullptr),
      fHistMCTrueSigmaPhotonPtCorr(nullptr),
      fHistMCTrueSigmaLambdaPCorr(nullptr),
      fHistMCTrueSigmaPhotonPCorr(nullptr),
      fHistMCBkgSigmaLambdaPtCorr(nullptr),
      fHistMCBkgSigmaPhotonPtCorr(nullptr),
      fHistMCBkgSigmaLambdaPCorr(nullptr),
      fHistMCBkgSigmaPhotonPCorr(nullptr),
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

  for (auto photon1 = photonCandidates.begin(); photon1 < photonCandidates.end(); ++photon1) {
    for (auto photon2 = photon1 + 1; photon2 < photonCandidates.end(); ++photon2) {
      fHistDiffPGamma[0]->Fill(photon1->GetPx() - photon2->GetPx());
      fHistDiffPGamma[1]->Fill(photon1->GetPy() - photon2->GetPy());
      fHistDiffPGamma[2]->Fill(photon1->GetPz() - photon2->GetPz());
      fHistDiffPGamma[3]->Fill(photon1->GetP() - photon2->GetP());

      fHistDiffPGammaPos[0]->Fill(photon1->GetPosDaughter().GetPx() - photon2->GetPosDaughter().GetPx());
      fHistDiffPGammaPos[1]->Fill(photon1->GetPosDaughter().GetPy() - photon2->GetPosDaughter().GetPy());
      fHistDiffPGammaPos[2]->Fill(photon1->GetPosDaughter().GetPz() - photon2->GetPosDaughter().GetPz());
      fHistDiffPGammaPos[3]->Fill(photon1->GetPosDaughter().GetP() - photon2->GetPosDaughter().GetP());

      fHistDiffPGammaNeg[0]->Fill(photon1->GetNegDaughter().GetPx() - photon2->GetNegDaughter().GetPx());
      fHistDiffPGammaNeg[1]->Fill(photon1->GetNegDaughter().GetPy() - photon2->GetNegDaughter().GetPy());
      fHistDiffPGammaNeg[2]->Fill(photon1->GetNegDaughter().GetPz() - photon2->GetNegDaughter().GetPz());
      fHistDiffPGammaNeg[3]->Fill(photon1->GetNegDaughter().GetP() - photon2->GetNegDaughter().GetP());
    }
  }

  for (const auto &lambda : lambdaCandidates) {
    // Lambda QA - keep track of kinematics
    fHistLambdaPtPhi->Fill(lambda.GetPt(), lambda.GetPhi());
    fHistLambdaPtEta->Fill(lambda.GetPt(), lambda.GetEta());
    fHistLambdaMassPt->Fill(lambda.GetPt(), lambda.GetMass());
  }

  for (auto lambda1 = lambdaCandidates.begin(); lambda1 < lambdaCandidates.end(); ++lambda1) {
    for (auto lambda2 = lambda1 + 1; lambda2 < lambdaCandidates.end(); ++lambda2) {
      fHistDiffPLambda[0]->Fill(lambda1->GetPx() - lambda2->GetPx());
      fHistDiffPLambda[1]->Fill(lambda1->GetPy() - lambda2->GetPy());
      fHistDiffPLambda[2]->Fill(lambda1->GetPz() - lambda2->GetPz());
      fHistDiffPLambda[3]->Fill(lambda1->GetP() - lambda2->GetP());

      fHistDiffPLambdaPos[0]->Fill(lambda1->GetPosDaughter().GetPx() - lambda2->GetPosDaughter().GetPx());
      fHistDiffPLambdaPos[1]->Fill(lambda1->GetPosDaughter().GetPy() - lambda2->GetPosDaughter().GetPy());
      fHistDiffPLambdaPos[2]->Fill(lambda1->GetPosDaughter().GetPz() - lambda2->GetPosDaughter().GetPz());
      fHistDiffPLambdaPos[3]->Fill(lambda1->GetPosDaughter().GetP() - lambda2->GetPosDaughter().GetP());

      fHistDiffPLambdaNeg[0]->Fill(lambda1->GetNegDaughter().GetPx() - lambda2->GetNegDaughter().GetPx());
      fHistDiffPLambdaNeg[1]->Fill(lambda1->GetNegDaughter().GetPy() - lambda2->GetNegDaughter().GetPy());
      fHistDiffPLambdaNeg[2]->Fill(lambda1->GetNegDaughter().GetPz() - lambda2->GetNegDaughter().GetPz());
      fHistDiffPLambdaNeg[3]->Fill(lambda1->GetNegDaughter().GetP() - lambda2->GetNegDaughter().GetP());
    }
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
      }
      // Armenteros cut
      if (fArmenterosCut) {
        if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
        if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
          continue;
      }
      const float rap = sigma.GetRapidity();
      const int multBin = GetMultiplicityBin(lPercentile, fMultMode);

      int label = -10;
      int pdgLambdaMother = 0;
      int pdgPhotonMother = 0;
      if (fIsMC) {
        label =
            sigma.MatchToMC(fMCEvent, fPDG, {{fPDGDaughter1, fPDGDaughter2}},
                            pdgLambdaMother, pdgPhotonMother);
      }

      // Now write out the stuff to the Femto containers
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        fSigma.push_back(sigma);
        if (!fIsLightweight) {
          fHistMassCutPt->Fill(pT);
          fHistEtaPhi->Fill(sigma.GetEta(), sigma.GetPhi());
          fHistPtRapidity->Fill(pT, rap);
          fHistSigmaLambdaPtCorr->Fill(pT, lambda.GetPt());
          fHistSigmaPhotonPtCorr->Fill(pT, photon.GetPt());
          fHistSigmaLambdaPCorr->Fill(sigma.GetP(), lambda.GetP());
          fHistSigmaPhotonPCorr->Fill(sigma.GetP(), photon.GetP());

          if (fIsMC && label > 0) {
            fHistMCTrueSigmaLambdaPtCorr->Fill(pT, lambda.GetPt());
            fHistMCTrueSigmaPhotonPtCorr->Fill(pT, photon.GetPt());
            fHistMCTrueSigmaLambdaPCorr->Fill(sigma.GetP(), lambda.GetP());
            fHistMCTrueSigmaPhotonPCorr->Fill(sigma.GetP(), photon.GetP());
          } else if (fIsMC && label <= 0) {
            fHistMCBkgSigmaLambdaPtCorr->Fill(pT, lambda.GetPt());
            fHistMCBkgSigmaPhotonPtCorr->Fill(pT, photon.GetPt());
            fHistMCBkgSigmaLambdaPCorr->Fill(sigma.GetP(), lambda.GetP());
            fHistMCBkgSigmaPhotonPCorr->Fill(sigma.GetP(), photon.GetP());
          }
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

      if (TMath::Abs(rap) > fRapidityMax) continue;
      fHistInvMassPt->Fill(pT, invMass);

      if (!fIsLightweight) {
        fHistArmenterosAfter->Fill(armAlpha, armQt);
        fHistInvMassRecPhoton->Fill(pT, sigma.GetRecMassPhoton());
        fHistInvMassRecLambda->Fill(pT, sigma.GetRecMassLambda());
        fHistInvMassRec->Fill(pT, sigma.GetRecMass());
      }
      if (multBin >= 0 && fIsSpectrumAnalysis)
        fHistPtMult[multBin]->Fill(pT, invMass);

      if (fIsMC) {
        if (label > 0) {
          fHistMCV0Pt->Fill(sigma.GetPt());
          fHistMCV0Mass->Fill(invMass);
        }
        if (!fIsLightweight) {
          // let's where the other particle comes from if one of them stems from
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
  if (!fIsLightweight) fHistNSigma->Fill(nSigma);
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
        fHistMixedInvMassPt->Fill(pT, invMass);
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
        fHistMixedInvMassPt->Fill(pT, invMass);
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
      fHistMCTruthDaughterPtYHighMult->Fill(mcParticle->Y(), mcParticle->Pt());
    }

    if (!CheckDaughtersInAcceptance(mcParticle)) continue;
    fHistMCTruthDaughterPtYAccept->Fill(mcParticle->Y(), mcParticle->Pt());
    if (lPercentile < fMCHighMultThreshold) {
      fHistMCTruthDaughterPtYAcceptHighMult->Fill(mcParticle->Y(),
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

  if (TMath::Abs(particle->Y()) > fRapidityMax) return false;

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
            << " Sigma0 mass      " << fMassSigma << "\n"
            << " Sigma0 selection " << fSigmaMassCut << "\n"
            << " Sigma0 sb up     " << fSidebandCutUp << "\n"
            << " Sigma0 sb down   " << fSidebandCutDown << "\n"
            << " Photon pT min    " << fPhotonPtMin << "\n"
            << " Photon pT max    " << fPhotonPtMax << "\n"
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

  fHistInvMass =
      new TH1F("fHistInvMass", "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries",
               150, 1.15, 1.3);
  fHistograms->Add(fHistInvMass);

  fHistInvMassPt = new TH2F("fHistInvMassPt",
                            "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                            "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                            100, 0, 10, 300, 1., 1.3);
  fHistograms->Add(fHistInvMassPt);
  fHistMixedInvMassPt = new TH2F("fHistMixedInvMassPt",
                                 "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                                 "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                                 100, 0, 10, 300, 1., 1.3);
  fHistograms->Add(fHistMixedInvMassPt);

  std::vector<float> multBinsLow, multBinsUp;
  if (fMultMode == AliVEvent::kINT7) {
    multBinsLow = {{0, 5., 15., 30., 50.}};
    multBinsUp = {{5., 15., 30., 50., 100.}};
  } else if (fMultMode == AliVEvent::kHighMultV0) {
    multBinsLow = {{0, 0.01, 0.05, 0.1, 1}};
    multBinsUp = {{0.01, 0.05, 0.1, 1, 100.}};
  }

  if (fIsSpectrumAnalysis) {
    for (int i = 0; i < static_cast<int>(multBinsUp.size()); ++i) {
      fHistPtMult[i] =
          new TH2F(Form("fHistPtMult_%i", i),
                   Form("V0M: %.2f - %.2f %%; #it{p}_{T} (GeV/#it{c}); "
                        "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                        multBinsLow[i], multBinsUp[i]),
                   100, 0, 10, 300, 1.15, 1.3);
      fHistograms->Add(fHistPtMult[i]);

      fHistMixedInvMassBinnedMultPt[i] =
          new TH2F(Form("fHistMixedInvMassBinnedMultPt_%i", i),
                   Form("V0M: %.2f - %.2f %%; #it{p}_{T} (GeV/#it{c}); "
                        "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                        multBinsLow[i], multBinsUp[i]),
                   100, 0, 10, 300, 1.15, 1.3);
      fHistograms->Add(fHistMixedInvMassBinnedMultPt[i]);
    }
  }

  if (!fIsLightweight) {
    fHistNSigma =
        new TH1F("fHistNSigma", ";# #Sigma candidates; Entries", 10, 0, 10);
    fHistMassCutPt = new TH1F(
        "fHistMassCutPt", "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries",
        100, 0, 10);

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
                 100, -1, 1, 100, 0, 0.5);
    fHistArmenterosAfter =
        new TH2F("fHistArmenterosAfter", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 100, -1, 1, 100, 0, 0.5);
    fHistEtaPhi = new TH2F("fHistEtaPhi", "; #eta; #phi", 100, -1, 1, 100,
                           -TMath::Pi(), TMath::Pi());
    fHistPtRapidity =
        new TH2F("fHistPtRapidity", "; #it{p}_{T} (GeV/#it{c}); y", 100, 0, 10,
                 50, -5, 5);

    fHistograms->Add(fHistNSigma);
    fHistograms->Add(fHistMassCutPt);
    fHistograms->Add(fHistInvMassRecPhoton);
    fHistograms->Add(fHistInvMassRecLambda);
    fHistograms->Add(fHistInvMassRec);
    fHistograms->Add(fHistEtaPhi);
    fHistograms->Add(fHistPtRapidity);
    fHistograms->Add(fHistArmenterosBefore);
    fHistograms->Add(fHistArmenterosAfter);

    TString coordinate[4] = {"x", "y", "z", "tot"};
    for (int i = 0; i < 4; ++i) {
      fHistDiffPGamma[i] = new TH1F(
              Form("fHistDiffPGamma_%s", coordinate[i].Data()),
              Form("; diff p_{%s,#gamma} (GeV/c); Entries", coordinate[i].Data()),
              1001, -1, 1);
      fHistDiffPLambda[i] = new TH1F(
              Form("fHistDiffPLambda_%s", coordinate[i].Data()),
              Form("; diff p_{%s,#Lambda} (GeV/c); Entries", coordinate[i].Data()),
              1001, -1, 1);
      fHistograms->Add(fHistDiffPGamma[i]);
      fHistograms->Add(fHistDiffPLambda[i]);

      fHistDiffPGammaPos[i] = new TH1F(
              Form("fHistDiffPGammaPos_%s", coordinate[i].Data()),
              Form("; diff p_{%s,#gamma pos} (GeV/c); Entries", coordinate[i].Data()),
              1001, -1, 1);
      fHistDiffPLambdaPos[i] = new TH1F(
              Form("fHistDiffPLambdaPos_%s", coordinate[i].Data()),
              Form("; diff p_{%s,#Lambda pos} (GeV/c); Entries", coordinate[i].Data()),
              1001, -1, 1);
      fHistograms->Add(fHistDiffPGammaPos[i]);
      fHistograms->Add(fHistDiffPLambdaPos[i]);

      fHistDiffPGammaNeg[i] = new TH1F(
              Form("fHistDiffPGammaNeg_%s", coordinate[i].Data()),
              Form("; diff p_{%s,#gamma neg} (GeV/c); Entries", coordinate[i].Data()),
              1001, -1, 1);
      fHistDiffPLambdaNeg[i] = new TH1F(
              Form("fHistDiffPLambdaNeg_%s", coordinate[i].Data()),
              Form("; diff p_{%s,#Lambda neg} (GeV/c); Entries", coordinate[i].Data()),
              1001, -1, 1);
      fHistograms->Add(fHistDiffPGammaNeg[i]);
      fHistograms->Add(fHistDiffPLambdaNeg[i]);
    }

    fHistLambdaPtPhi =
        new TH2F("fHistLambdaPtPhi", "; #it{p}_{T} (GeV/#it{c}); #phi (rad)",
                 50, 0, 10, 50, 0, 2.f * TMath::Pi());
    fHistLambdaPtEta =
        new TH2F("fHistLambdaPtEta", "; #it{p}_{T} (GeV/#it{c}); #eta", 50, 0,
                 10, 50, -1, 1);
    fHistLambdaMassPt =
        new TH2F("fHistLambdaMassPt", "; #it{p}_{T} (GeV/#it{c}); M_{p#pi}", 50,
                 0, 10, 50, 1., 1.3);
    fHistPhotonPtPhi =
        new TH2F("fHistPhotonPtPhi", "; #it{p}_{T} (GeV/#it{c}); #phi (rad)",
                 50, 0, 10, 50, 0, 2.f * TMath::Pi());
    fHistPhotonPtEta =
        new TH2F("fHistPhotonPtEta", "; #it{p}_{T} (GeV/#it{c}); #eta", 50, 0,
                 10, 50, -1, 1);
    fHistPhotonMassPt =
        new TH2F("fHistPhotonMassPt", "; #it{p}_{T} (GeV/#it{c}); M_{p#pi}", 50,
                 0, 10, 50, 0., 0.1);

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

    for (int i = 0; i < static_cast<int>(multBinsUp.size()); ++i) {
      fHistMCTruthPtMult[i] =
          new TH1F(Form("fHistMCTruthPtMult%i", i),
                   Form("V0M: %.2f - %.2f %%; #it{p}_{T} (GeV/#it{c}); Entries",
                        multBinsLow[i], multBinsUp[i]),
                   100, 0, 10);
      fHistogramsMC->Add(fHistMCTruthPtMult[i]);
    }

    fHistMCTruthPtY =
        new TH2F("fHistMCTruthPtY", "; y; #it{p}_{T} (GeV/#it{c})", 100, -5, 5,
                 100, 0, 10);
    fHistMCTruthDaughterPtY =
        new TH2F("fHistMCTruthDaughterPtY", "; y; #it{p}_{T} (GeV/#it{c})", 200,
                 -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtYAccept =
        new TH2F("fHistMCTruthDaughterPtYAccept",
                 "; y; #it{p}_{T} (GeV/#it{c})", 200, -10, 10, 100, 0, 10);
    fHistMCTruthPtYHighMult =
        new TH2F("fHistMCTruthPtYHighMult", "; y; #it{p}_{T} (GeV/#it{c})", 200,
                 -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtYHighMult =
        new TH2F("fHistMCTruthDaughterPtYHighMult",
                 "; y; #it{p}_{T} (GeV/#it{c})", 200, -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtYAcceptHighMult =
        new TH2F("fHistMCTruthDaughterPtYAcceptHighMult",
                 "; y; #it{p}_{T} (GeV/#it{c})", 200, -10, 10, 100, 0, 10);

    fHistMCV0Pt = new TH1F("fHistMCV0Pt",
                           "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries",
                           100, 0, 10);
    fHistMCV0Mass =
        new TH1F("fHistMCV0Mass",
                 "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 500, 1., 1.5);
    fHistogramsMC->Add(fHistMCV0Pt);
    fHistogramsMC->Add(fHistMCV0Mass);

    if (!fIsLightweight) {
      fHistMCTrueSigmaLambdaPtCorr =
          new TH2F("fHistMCTrueSigmaLambdaPtCorr",
                   "; #Sigma^{0} #it{p}_{T} (GeV/#it{c}; "
                   "#Lambda #it{p}_{T} (GeV/#it{c}",
                   250, 0, 10, 250, 0, 10);
      fHistMCTrueSigmaPhotonPtCorr =
          new TH2F("fHistMCTrueSigmaPhotonPtCorr",
                   "; #Sigma^{0} #it{p}_{T} (GeV/#it{c}; "
                   "#gamma #it{p}_{T} (GeV/#it{c}",
                   250, 0, 10, 250, 0, 10);
      fHistMCTrueSigmaLambdaPCorr = new TH2F("fHistMCTrueSigmaLambdaPCorr",
                                             "; #Sigma^{0} #it{p} (GeV/#it{c}; "
                                             "#Lambda #it{p} (GeV/#it{c}",
                                             250, 0, 10, 250, 0, 10);
      fHistMCTrueSigmaPhotonPCorr = new TH2F("fHistMCTrueSigmaPhotonPCorr",
                                             "; #Sigma^{0} #it{p} (GeV/#it{c}; "
                                             "#gamma #it{p} (GeV/#it{c}",
                                             250, 0, 10, 250, 0, 10);
      fHistMCBkgSigmaLambdaPtCorr =
          new TH2F("fHistMCBkgSigmaLambdaPtCorr",
                   "; #Sigma^{0} #it{p}_{T} (GeV/#it{c}; "
                   "#Lambda #it{p}_{T} (GeV/#it{c}",
                   250, 0, 10, 250, 0, 10);
      fHistMCBkgSigmaPhotonPtCorr =
          new TH2F("fHistMCBkgSigmaPhotonPtCorr",
                   "; #Sigma^{0} #it{p}_{T} (GeV/#it{c}; "
                   "#gamma #it{p}_{T} (GeV/#it{c}",
                   250, 0, 10, 250, 0, 10);
      fHistMCBkgSigmaLambdaPCorr = new TH2F("fHistMCBkgSigmaLambdaPCorr",
                                            "; #Sigma^{0} #it{p} (GeV/#it{c}; "
                                            "#Lambda #it{p} (GeV/#it{c}",
                                            250, 0, 10, 250, 0, 10);
      fHistMCBkgSigmaPhotonPCorr = new TH2F("fHistMCBkgSigmaPhotonPCorr",
                                            "; #Sigma^{0} #it{p} (GeV/#it{c}; "
                                            "#gamma #it{p} (GeV/#it{c}",
                                            250, 0, 10, 250, 0, 10);

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
      fHistogramsMC->Add(fHistMCTruthDaughterPtY);
      fHistogramsMC->Add(fHistMCTruthDaughterPtYAccept);
      fHistogramsMC->Add(fHistMCTruthPtYHighMult);
      fHistogramsMC->Add(fHistMCTruthDaughterPtYHighMult);
      fHistogramsMC->Add(fHistMCTruthDaughterPtYAcceptHighMult);
      fHistogramsMC->Add(fHistMCTrueSigmaLambdaPtCorr);
      fHistogramsMC->Add(fHistMCTrueSigmaPhotonPtCorr);
      fHistogramsMC->Add(fHistMCTrueSigmaLambdaPCorr);
      fHistogramsMC->Add(fHistMCTrueSigmaPhotonPCorr);
      fHistogramsMC->Add(fHistMCBkgSigmaLambdaPtCorr);
      fHistogramsMC->Add(fHistMCBkgSigmaPhotonPtCorr);
      fHistogramsMC->Add(fHistMCBkgSigmaLambdaPCorr);
      fHistogramsMC->Add(fHistMCBkgSigmaPhotonPCorr);
      fHistogramsMC->Add(fHistMCV0Mother);
      fHistogramsMC->Add(fHistMCV0Check);
      fHistogramsMC->Add(fHistMCV0MotherCheck);
    }

    fHistograms->Add(fHistogramsMC);
  }
}
