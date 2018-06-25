#include "AliSigma0PhotonMotherCuts.h"
#include <iostream>

ClassImp(AliSigma0PhotonMotherCuts)

    //____________________________________________________________________________________________________
    AliSigma0PhotonMotherCuts::AliSigma0PhotonMotherCuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsSigma(nullptr),
      fHistogramsSigmaMC(nullptr),
      fHistogramsAntiSigma(nullptr),
      fHistogramsAntiSigmaMC(nullptr),
      fIsMC(false),
      fIsLightweight(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fLambdaMixed(),
      fAntiLambdaMixed(),
      fPhotonMixed(),
      fMixingDepth(10),
      fMassSigma(0),
      fSigmaMassCut(0),
      fSigmaSidebandLow(0),
      fSigmaSidebandUp(0),
      fArmenterosCut(true),
      fArmenterosQtLow(0.f),
      fArmenterosQtUp(0.f),
      fArmenterosAlphaLow(0.f),
      fArmenterosAlphaUp(0.f),
      fHistCuts(nullptr),
      fHistNSigma(nullptr),
      fHistSigmaPt(nullptr),
      fHistSigmaMassCutPt(nullptr),
      fHistSigmaInvMass(nullptr),
      fHistSigmaInvMassRec(nullptr),
      fHistSigmaInvMassPt(nullptr),
      fHistSigmaInvMassEta(nullptr),
      fHistSigmaEtaPhi(nullptr),
      fHistSigmaRapidity(nullptr),
      fHistSigmaPtY(),
      fHistSigmaArmenterosBefore(nullptr),
      fHistSigmaArmenterosAfter(nullptr),
      fHistSigmaMixedPt(nullptr),
      fHistSigmaMixedInvMass(nullptr),
      fHistSigmaMixedPtY(),
      fHistSigmaMixedInvMassPt(nullptr),
      fHistSigmaMixedInvMassEta(nullptr),
      fHistNAntiSigma(nullptr),
      fHistAntiSigmaPt(nullptr),
      fHistAntiSigmaMassCutPt(nullptr),
      fHistAntiSigmaInvMass(nullptr),
      fHistAntiSigmaInvMassRec(nullptr),
      fHistAntiSigmaInvMassPt(nullptr),
      fHistAntiSigmaInvMassEta(nullptr),
      fHistAntiSigmaEtaPhi(nullptr),
      fHistAntiSigmaRapidity(nullptr),
      fHistAntiSigmaPtY(),
      fHistAntiSigmaArmenterosBefore(nullptr),
      fHistAntiSigmaArmenterosAfter(nullptr),
      fHistAntiSigmaMixedPt(nullptr),
      fHistAntiSigmaMixedInvMass(nullptr),
      fHistAntiSigmaMixedPtY(),
      fHistAntiSigmaMixedInvMassPt(nullptr),
      fHistAntiSigmaMixedInvMassEta(nullptr),
      fHistMCSigmaMassCutPt(nullptr),
      fHistMCTruthSigma0PhotonConvPt(nullptr),
      fHistMCTruthSigma0PhotonConvP(nullptr),
      fHistMCTruthSigma0PhotonConvInvMass(nullptr),
      fHistMCTruthSigma0PhotonConvInvMassPt(nullptr),
      fHistMCTruthSigma0PhotonConvPtEta(nullptr),
      fHistMCTruthSigma0PhotonConvR(nullptr),
      fHistMCTruthSigma0PhotonConvConvPointX(nullptr),
      fHistMCTruthSigma0PhotonConvConvPointY(nullptr),
      fHistMCTruthSigma0PhotonConvConvPointZ(nullptr),
      fHistMCTruthSigma0PhotonConvEleP(nullptr),
      fHistMCTruthSigma0PhotonConvElePt(nullptr),
      fHistMCTruthSigma0PhotonConvPtY(nullptr),
      fHistMCTruthSigma0Pt(nullptr),
      fHistMCTruthSigma0PtY(nullptr),
      fHistMCTruthSigma0PtEta(nullptr),
      fHistMCTruthSigma0PhotonPt(nullptr),
      fHistMCTruthSigma0PhotonPtY(nullptr),
      fHistMCTruthSigma0PhotonPtEta(nullptr),
      fHistMCAntiSigmaMassCutPt(nullptr),
      fHistMCTruthAntiSigma0PhotonConvPt(nullptr),
      fHistMCTruthAntiSigma0PhotonConvP(nullptr),
      fHistMCTruthAntiSigma0PhotonConvInvMass(nullptr),
      fHistMCTruthAntiSigma0PhotonConvInvMassPt(nullptr),
      fHistMCTruthAntiSigma0PhotonConvPtEta(nullptr),
      fHistMCTruthAntiSigma0PhotonConvR(nullptr),
      fHistMCTruthAntiSigma0PhotonConvConvPointX(nullptr),
      fHistMCTruthAntiSigma0PhotonConvConvPointY(nullptr),
      fHistMCTruthAntiSigma0PhotonConvConvPointZ(nullptr),
      fHistMCTruthAntiSigma0PhotonConvEleP(nullptr),
      fHistMCTruthAntiSigma0PhotonConvElePt(nullptr),
      fHistMCTruthAntiSigma0PhotonConvPtY(nullptr),
      fHistMCTruthAntiSigma0Pt(nullptr),
      fHistMCTruthAntiSigma0PtY(nullptr),
      fHistMCTruthAntiSigma0PtEta(nullptr),
      fHistMCTruthAntiSigma0PhotonPt(nullptr),
      fHistMCTruthAntiSigma0PhotonPtY(nullptr),
      fHistMCTruthAntiSigma0PhotonPtEta(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0PhotonMotherCuts::AliSigma0PhotonMotherCuts(
    const AliSigma0PhotonMotherCuts &ref)
    : TObject(ref),
      fHistograms(nullptr),
      fHistogramsSigma(nullptr),
      fHistogramsSigmaMC(nullptr),
      fHistogramsAntiSigma(nullptr),
      fHistogramsAntiSigmaMC(nullptr),
      fIsMC(false),
      fIsLightweight(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fLambdaMixed(),
      fAntiLambdaMixed(),
      fPhotonMixed(),
      fMixingDepth(10),
      fMassSigma(0),
      fSigmaMassCut(0),
      fSigmaSidebandLow(0),
      fSigmaSidebandUp(0),
      fArmenterosCut(true),
      fArmenterosQtLow(0.f),
      fArmenterosQtUp(0.f),
      fArmenterosAlphaLow(0.f),
      fArmenterosAlphaUp(0.f),
      fHistCuts(nullptr),
      fHistNSigma(nullptr),
      fHistSigmaPt(nullptr),
      fHistSigmaMassCutPt(nullptr),
      fHistSigmaInvMass(nullptr),
      fHistSigmaInvMassRec(nullptr),
      fHistSigmaInvMassPt(nullptr),
      fHistSigmaInvMassEta(nullptr),
      fHistSigmaEtaPhi(nullptr),
      fHistSigmaRapidity(nullptr),
      fHistSigmaPtY(),
      fHistSigmaArmenterosBefore(nullptr),
      fHistSigmaArmenterosAfter(nullptr),
      fHistSigmaMixedPt(nullptr),
      fHistSigmaMixedInvMass(nullptr),
      fHistSigmaMixedPtY(),
      fHistSigmaMixedInvMassPt(nullptr),
      fHistSigmaMixedInvMassEta(nullptr),
      fHistNAntiSigma(nullptr),
      fHistAntiSigmaPt(nullptr),
      fHistAntiSigmaMassCutPt(nullptr),
      fHistAntiSigmaInvMass(nullptr),
      fHistAntiSigmaInvMassRec(nullptr),
      fHistAntiSigmaInvMassPt(nullptr),
      fHistAntiSigmaInvMassEta(nullptr),
      fHistAntiSigmaEtaPhi(nullptr),
      fHistAntiSigmaRapidity(nullptr),
      fHistAntiSigmaPtY(),
      fHistAntiSigmaArmenterosBefore(nullptr),
      fHistAntiSigmaArmenterosAfter(nullptr),
      fHistAntiSigmaMixedPt(nullptr),
      fHistAntiSigmaMixedInvMass(nullptr),
      fHistAntiSigmaMixedPtY(),
      fHistAntiSigmaMixedInvMassPt(nullptr),
      fHistAntiSigmaMixedInvMassEta(nullptr),
      fHistMCSigmaMassCutPt(nullptr),
      fHistMCTruthSigma0PhotonConvPt(nullptr),
      fHistMCTruthSigma0PhotonConvP(nullptr),
      fHistMCTruthSigma0PhotonConvInvMass(nullptr),
      fHistMCTruthSigma0PhotonConvInvMassPt(nullptr),
      fHistMCTruthSigma0PhotonConvPtEta(nullptr),
      fHistMCTruthSigma0PhotonConvR(nullptr),
      fHistMCTruthSigma0PhotonConvConvPointX(nullptr),
      fHistMCTruthSigma0PhotonConvConvPointY(nullptr),
      fHistMCTruthSigma0PhotonConvConvPointZ(nullptr),
      fHistMCTruthSigma0PhotonConvEleP(nullptr),
      fHistMCTruthSigma0PhotonConvElePt(nullptr),
      fHistMCTruthSigma0PhotonConvPtY(nullptr),
      fHistMCTruthSigma0Pt(nullptr),
      fHistMCTruthSigma0PtY(nullptr),
      fHistMCTruthSigma0PtEta(nullptr),
      fHistMCTruthSigma0PhotonPt(nullptr),
      fHistMCTruthSigma0PhotonPtY(nullptr),
      fHistMCTruthSigma0PhotonPtEta(nullptr),
      fHistMCAntiSigmaMassCutPt(nullptr),
      fHistMCTruthAntiSigma0PhotonConvPt(nullptr),
      fHistMCTruthAntiSigma0PhotonConvP(nullptr),
      fHistMCTruthAntiSigma0PhotonConvInvMass(nullptr),
      fHistMCTruthAntiSigma0PhotonConvInvMassPt(nullptr),
      fHistMCTruthAntiSigma0PhotonConvPtEta(nullptr),
      fHistMCTruthAntiSigma0PhotonConvR(nullptr),
      fHistMCTruthAntiSigma0PhotonConvConvPointX(nullptr),
      fHistMCTruthAntiSigma0PhotonConvConvPointY(nullptr),
      fHistMCTruthAntiSigma0PhotonConvConvPointZ(nullptr),
      fHistMCTruthAntiSigma0PhotonConvEleP(nullptr),
      fHistMCTruthAntiSigma0PhotonConvElePt(nullptr),
      fHistMCTruthAntiSigma0PhotonConvPtY(nullptr),
      fHistMCTruthAntiSigma0Pt(nullptr),
      fHistMCTruthAntiSigma0PtY(nullptr),
      fHistMCTruthAntiSigma0PtEta(nullptr),
      fHistMCTruthAntiSigma0PhotonPt(nullptr),
      fHistMCTruthAntiSigma0PhotonPtY(nullptr),
      fHistMCTruthAntiSigma0PhotonPtEta(nullptr) {}

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
  photonMotherCuts->SetArmenterosCut(0., 0.12, -1., -0.6);
  return photonMotherCuts;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SelectPhotonMother(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    std::vector<AliAODConversionPhoton> &photonCandidates,
    std::vector<AliSigma0ParticleV0> &lambdaCandidates,
    std::vector<AliSigma0ParticleV0> &antiLambdaCandidates) {
  fMCEvent = mcEvent;

  fInputEvent = inputEvent;

  // Particle pairing
  SigmaToLambdaGamma(photonCandidates, lambdaCandidates, antiLambdaCandidates);
  SigmaToLambdaGammaMixedEvent(photonCandidates, lambdaCandidates,
                               antiLambdaCandidates);
  FillEventBuffer(photonCandidates, lambdaCandidates, antiLambdaCandidates);
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGamma(
    std::vector<AliAODConversionPhoton> &photonCandidates,
    std::vector<AliSigma0ParticleV0> &lambdaCandidates,
    std::vector<AliSigma0ParticleV0> &antiLambdaCandidates) {
  // SAME EVENT
  for (const auto &photon : photonCandidates) {
    int nAntiSigma = 0;
    for (const auto &antiLambda : antiLambdaCandidates) {
      /// Candidates with lambdas with shared daughter tracks are not used
      if (!antiLambda.GetIsUse()) continue;
      const AliSigma0ParticlePhotonMother antiSigma(antiLambda, photon,
                                                    fInputEvent);
      const float armAlpha = antiSigma.GetArmenterosAlpha();
      const float armQt = antiSigma.GetArmenterosQt();
      if (!fIsLightweight)
        fHistAntiSigmaArmenterosBefore->Fill(armAlpha, armQt);
      const float invMass = antiSigma.GetMass();
      // Armenteros cut
      if (fArmenterosCut) {
        if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
        if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
          continue;
      }
      const float rap = antiSigma.GetRapidity();
      const int rapBin = GetRapidityBin(rap);
      if (!fIsLightweight) {
        fHistAntiSigmaArmenterosAfter->Fill(armAlpha, armQt);
        fHistAntiSigmaPt->Fill(antiSigma.GetPt());
        fHistAntiSigmaInvMass->Fill(invMass);
        fHistAntiSigmaInvMassRec->Fill(antiSigma.GetRecMass());
        fHistAntiSigmaRapidity->Fill(rap);
        if (rapBin > -1)
          fHistAntiSigmaPtY[rapBin]->Fill(antiSigma.GetPt(), invMass);
        fHistAntiSigmaInvMassEta->Fill(antiSigma.GetEta(), invMass);
      }
      fHistAntiSigmaInvMassPt->Fill(antiSigma.GetPt(), invMass);
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        if (!fIsLightweight) {
          fHistAntiSigmaMassCutPt->Fill(antiSigma.GetPt());
          fHistAntiSigmaEtaPhi->Fill(antiSigma.GetEta(), antiSigma.GetPhi());
        }
        ++nAntiSigma;
        if (fIsMC) {
          if (antiSigma.IsTrueSigma(fMCEvent)) {
            fHistMCAntiSigmaMassCutPt->Fill(antiSigma.GetPt());
          }
        }
      }
    }
    if (!fIsLightweight) fHistNAntiSigma->Fill(nAntiSigma);

    int nSigma = 0;
    // Sigmas
    for (const auto &lambda : lambdaCandidates) {
      /// Candidates with lambdas with shared daughter tracks are not used
      if (!lambda.GetIsUse()) continue;
      const AliSigma0ParticlePhotonMother sigma(lambda, photon, fInputEvent);
      const float armAlpha = sigma.GetArmenterosAlpha();
      const float armQt = sigma.GetArmenterosQt();
      if (!fIsLightweight) fHistSigmaArmenterosBefore->Fill(armAlpha, armQt);
      const float invMass = sigma.GetMass();
      // Armenteros cut
      if (fArmenterosCut) {
        if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
        if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
          continue;
      }
      const float rap = sigma.GetRapidity();
      const int rapBin = GetRapidityBin(rap);
      if (!fIsLightweight) {
        fHistSigmaArmenterosAfter->Fill(armAlpha, armQt);
        fHistSigmaPt->Fill(sigma.GetPt());
        fHistSigmaInvMass->Fill(invMass);
        fHistSigmaInvMassRec->Fill(sigma.GetRecMass());
        fHistSigmaRapidity->Fill(rap);
        if (rapBin > -1) fHistSigmaPtY[rapBin]->Fill(sigma.GetPt(), invMass);
        fHistSigmaInvMassEta->Fill(sigma.GetEta(), invMass);
      }
      fHistSigmaInvMassPt->Fill(sigma.GetPt(), invMass);
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        if (!fIsLightweight) {
          fHistSigmaMassCutPt->Fill(sigma.GetPt());
          fHistSigmaEtaPhi->Fill(sigma.GetEta(), sigma.GetPhi());
        }
        ++nSigma;
        if (fIsMC) {
          if (sigma.IsTrueSigma(fMCEvent)) {
            fHistMCSigmaMassCutPt->Fill(sigma.GetPt());
          }
        }
      }
    }
    if (!fIsLightweight) fHistNSigma->Fill(nSigma);
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGammaMixedEvent(
    std::vector<AliAODConversionPhoton> &photonCandidates,
    std::vector<AliSigma0ParticleV0> &lambdaCandidates,
    std::vector<AliSigma0ParticleV0> &antiLambdaCandidates) {
  // Anti - Sigma0
  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fAntiLambdaMixed) {
    for (const auto &Lambda : LambdaContainer) {
      if (!Lambda.GetIsUse()) continue;
      for (auto Photon : photonCandidates) {
        const AliSigma0ParticlePhotonMother sigma(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma.GetArmenterosAlpha();
        const float armQt = sigma.GetArmenterosQt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
            continue;
        }
        const float invMass = sigma.GetMass();
        fHistAntiSigmaMixedInvMassPt->Fill(sigma.GetPt(), invMass);
        if (!fIsLightweight) {
          fHistAntiSigmaMixedPt->Fill(sigma.GetPt());
          fHistAntiSigmaMixedInvMass->Fill(invMass);
          const float rap = sigma.GetRapidity();
          const int rapBin = GetRapidityBin(rap);
          if (rapBin > -1)
            fHistAntiSigmaMixedPtY[rapBin]->Fill(sigma.GetPt(), invMass);
          fHistAntiSigmaMixedInvMassEta->Fill(sigma.GetEta(), invMass);
        }
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto &PhotonContainer : fPhotonMixed) {
    for (const auto &Photon : PhotonContainer) {
      for (const auto &Lambda : antiLambdaCandidates) {
        if (!Lambda.GetIsUse()) continue;
        const AliSigma0ParticlePhotonMother sigma(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma.GetArmenterosAlpha();
        const float armQt = sigma.GetArmenterosQt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
            continue;
        }
        const float invMass = sigma.GetMass();
        fHistAntiSigmaMixedInvMassPt->Fill(sigma.GetPt(), invMass);
        if (!fIsLightweight) {
          fHistAntiSigmaMixedPt->Fill(sigma.GetPt());
          fHistAntiSigmaMixedInvMass->Fill(invMass);
          const float rap = sigma.GetRapidity();
          const int rapBin = GetRapidityBin(rap);
          if (rapBin > -1)
            fHistAntiSigmaMixedPtY[rapBin]->Fill(sigma.GetPt(), invMass);
          fHistAntiSigmaMixedInvMassEta->Fill(sigma.GetEta(), invMass);
        }
      }
    }
  }

  // Sigma0
  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fLambdaMixed) {
    for (const auto &Lambda : LambdaContainer) {
      if (!Lambda.GetIsUse()) continue;
      for (auto Photon : photonCandidates) {
        const AliSigma0ParticlePhotonMother sigma(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma.GetArmenterosAlpha();
        const float armQt = sigma.GetArmenterosQt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
            continue;
        }
        const float invMass = sigma.GetMass();
        fHistSigmaMixedInvMassPt->Fill(sigma.GetPt(), invMass);
        if (!fIsLightweight) {
          fHistSigmaMixedPt->Fill(sigma.GetPt());
          fHistSigmaMixedInvMass->Fill(invMass);
          const float rap = sigma.GetRapidity();
          const int rapBin = GetRapidityBin(rap);
          if (rapBin > -1)
            fHistSigmaMixedPtY[rapBin]->Fill(sigma.GetPt(), invMass);
          fHistSigmaMixedInvMassEta->Fill(sigma.GetEta(), invMass);
        }
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto &PhotonContainer : fPhotonMixed) {
    for (const auto &Photon : PhotonContainer) {
      for (const auto &Lambda : lambdaCandidates) {
        if (!Lambda.GetIsUse()) continue;
        const AliSigma0ParticlePhotonMother sigma(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma.GetArmenterosAlpha();
        const float armQt = sigma.GetArmenterosQt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) continue;
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow)
            continue;
        }
        const float invMass = sigma.GetMass();
        fHistSigmaMixedInvMassPt->Fill(sigma.GetPt(), invMass);
        if (!fIsLightweight) {
          fHistSigmaMixedPt->Fill(sigma.GetPt());
          fHistSigmaMixedInvMass->Fill(invMass);
          const float rap = sigma.GetRapidity();
          const int rapBin = GetRapidityBin(rap);
          if (rapBin > -1)
            fHistSigmaMixedPtY[rapBin]->Fill(sigma.GetPt(), invMass);
          fHistSigmaMixedInvMassEta->Fill(sigma.GetEta(), invMass);
        }
      }
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::FillEventBuffer(
    std::vector<AliAODConversionPhoton> &photonCandidates,
    std::vector<AliSigma0ParticleV0> &lambdaCandidates,
    std::vector<AliSigma0ParticleV0> &antiLambdaCandidates) {
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

  // ++++++++++++++
  // Anti - Lambda
  if (static_cast<int>(antiLambdaCandidates.size()) > 0) {
    if (static_cast<int>(fAntiLambdaMixed.size()) < fMixingDepth) {
      fAntiLambdaMixed.push_back(antiLambdaCandidates);
    } else {
      fAntiLambdaMixed.pop_front();
      fAntiLambdaMixed.push_back(antiLambdaCandidates);
    }
  }
}

//____________________________________________________________________________________________________
int AliSigma0PhotonMotherCuts::GetRapidityBin(float rapidity) const {
  if (rapidity <= -1.0)
    return 0;
  else if (-1.0 < rapidity && rapidity <= -0.9)
    return 1;
  else if (-0.9 < rapidity && rapidity <= -0.8)
    return 2;
  else if (-0.8 < rapidity && rapidity <= -0.7)
    return 3;
  else if (-0.7 < rapidity && rapidity <= -0.6)
    return 4;
  else if (-0.6 < rapidity && rapidity <= -0.5)
    return 5;
  else if (-0.5 < rapidity && rapidity <= -0.4)
    return 6;
  else if (-0.4 < rapidity && rapidity <= -0.3)
    return 7;
  else if (-0.3 < rapidity && rapidity <= -0.2)
    return 8;
  else if (-0.2 < rapidity && rapidity <= -0.1)
    return 9;
  else if (-0.1 < rapidity && rapidity <= 0.f)
    return 10;
  else if (0.f < rapidity && rapidity <= 0.1)
    return 11;
  else if (0.1 < rapidity && rapidity <= 0.2)
    return 12;
  else if (0.2 < rapidity && rapidity <= 0.3)
    return 13;
  else if (0.3 < rapidity && rapidity <= 0.4)
    return 14;
  else if (0.4 < rapidity && rapidity <= 0.5)
    return 15;
  else if (0.5 < rapidity && rapidity <= 0.6)
    return 16;
  else if (0.6 < rapidity && rapidity <= 0.7)
    return 17;
  else if (0.7 < rapidity && rapidity <= 0.8)
    return 18;
  else if (0.8 < rapidity && rapidity <= 0.9)
    return 19;
  else if (0.9 < rapidity && rapidity <= 1.f)
    return 20;
  else if (1.0 < rapidity)
    return 21;
  else
    return -1;
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::InitCutHistograms(TString appendix) {
  std::cout << "============================\n"
            << " PHOTON MOTHER CUT CONFIGURATION \n"
            << " Sigma0 mass     " << fMassSigma << "\n"
            << " Sigma0 select   " << fSigmaMassCut << "\n"
            << " Sigma0 sideb    " << fSigmaSidebandLow << " "
            << fSigmaSidebandUp << "\n"
            << " Mixing depth    " << fMixingDepth << "\n"
            << "============================\n";

  std::vector<float> rapBins = {
      {-10.f, -1.f, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.f,
       0.1,   0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.f,  10.f}};

  TH1::AddDirectory(kFALSE);

  TString name;
  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    name = "PhotonMotherCut" + appendix;
    fHistograms->SetName(name);
  }

  if (fHistogramsSigma != nullptr) {
    delete fHistogramsSigma;
    fHistogramsSigma = nullptr;
  }
  if (fHistogramsSigma == nullptr) {
    fHistogramsSigma = new TList();
    fHistogramsSigma->SetOwner(kTRUE);
    name = "SigmaCut_QA" + appendix;
    fHistogramsSigma->SetName(name);
  }

  fHistCuts = new TProfile("fHistCuts", ";;Cut value", 10, 0, 10);
  fHistCuts->GetXaxis()->SetBinLabel(1, "#Sigma^{0} selection");
  fHistCuts->GetXaxis()->SetBinLabel(2, "#Sigma^{0} sideband low");
  fHistCuts->GetXaxis()->SetBinLabel(3, "#Sigma^{0} sideband up");
  fHistCuts->GetXaxis()->SetBinLabel(4, "#Sigma^{0} mixing depth");
  fHistCuts->GetXaxis()->SetBinLabel(5, "Armenteros q_{T} low");
  fHistCuts->GetXaxis()->SetBinLabel(6, "Armenteros q_{T} up");
  fHistCuts->GetXaxis()->SetBinLabel(7, "Armenteros #alpha low");
  fHistCuts->GetXaxis()->SetBinLabel(8, "Armenteros #alpha up");
  fHistograms->Add(fHistCuts);

  fHistCuts->Fill(0.f, fSigmaMassCut);
  fHistCuts->Fill(1.f, fSigmaSidebandLow);
  fHistCuts->Fill(2.f, fSigmaSidebandUp);
  fHistCuts->Fill(3.f, fMixingDepth);
  fHistCuts->Fill(4.f, fArmenterosQtLow);
  fHistCuts->Fill(5.f, fArmenterosQtUp);
  fHistCuts->Fill(6.f, fArmenterosAlphaLow);
  fHistCuts->Fill(7.f, fArmenterosAlphaUp);

  fHistSigmaInvMassPt = new TH2F("fHistSigmaInvMassPt",
                                 "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                                 "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                                 1000, 0, 20, 2000, 1., 2.);
  fHistogramsSigma->Add(fHistSigmaInvMassPt);
  fHistSigmaMixedInvMassPt =
      new TH2F("fHistSigmaMixedInvMassPt",
               "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
               "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassPt);

  if (!fIsLightweight) {
    fHistNSigma =
        new TH1F("fHistNSigma", ";# #Sigma candidates; Entries", 20, 0, 20);
    fHistSigmaPt = new TH1F("fHistSigmaPt",
                            "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries",
                            500, 0, 20);
    fHistSigmaMassCutPt = new TH1F(
        "fHistSigmaMassCutPt",
        "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries", 500, 0, 20);
    fHistSigmaInvMass =
        new TH1F("fHistSigmaInvMass",
                 "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
    fHistSigmaInvMassRec =
        new TH1F("fHistSigmaInvMassRec",
                 "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
    fHistSigmaInvMassEta = new TH2F(
        "fHistSigmaInvMassEta", "; #eta; M_{#Lambda#gamma} (GeV/#it{c}^{2})",
        1000, 0, 1, 2000, 1., 2.);
    fHistSigmaArmenterosBefore =
        new TH2F("fHistSigmaArmenterosBefore",
                 " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
    fHistSigmaArmenterosAfter =
        new TH2F("fHistSigmaArmenterosAfter",
                 " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
    fHistSigmaEtaPhi = new TH2F("fHistSigmaEtaPhi", "; #eta; #phi", 200, -1, 1,
                                200, 0, 2 * TMath::Pi());
    fHistSigmaRapidity =
        new TH1F("fHistSigmaRapidity", "; y; Entries", 200, -1, 1);
    fHistSigmaMixedPt = new TH1F(
        "fHistSigmaMixedPt", "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries",
        500, 0, 20);
    fHistSigmaMixedInvMass =
        new TH1F("fHistSigmaMixedInvMass",
                 "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
    fHistSigmaMixedInvMassEta = new TH2F(
        "fHistSigmaMixedInvMassEta",
        "; #eta; M_{#Lambda#gamma} (GeV/#it{c}^{2})", 1000, 0, 1, 2000, 1., 2.);

    fHistogramsSigma->Add(fHistNSigma);
    fHistogramsSigma->Add(fHistSigmaPt);
    fHistogramsSigma->Add(fHistSigmaMassCutPt);
    fHistogramsSigma->Add(fHistSigmaInvMass);
    fHistogramsSigma->Add(fHistSigmaInvMassRec);
    fHistogramsSigma->Add(fHistSigmaInvMassEta);
    fHistogramsSigma->Add(fHistSigmaEtaPhi);
    fHistogramsSigma->Add(fHistSigmaRapidity);
    fHistogramsSigma->Add(fHistSigmaArmenterosBefore);
    fHistogramsSigma->Add(fHistSigmaArmenterosAfter);
    fHistogramsSigma->Add(fHistSigmaMixedPt);
    fHistogramsSigma->Add(fHistSigmaMixedInvMass);
    fHistogramsSigma->Add(fHistSigmaMixedInvMassEta);

    for (int i = 0; i < static_cast<int>(rapBins.size() - 1); ++i) {
      fHistSigmaPtY[i] =
          new TH2F(Form("fHistSigmaPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
                   Form("%.2f < y < %.2f ; #it{p}_{T} [GeV/#it{c}]; "
                        "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                        rapBins[i], rapBins[i + 1]),
                   500, 0, 10, 1000, 1., 1.5);
      fHistSigmaMixedPtY[i] = new TH2F(
          Form("fHistSigmaMixedPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
          Form("%.2f < y < %.2f ; #it{p}_{T} [GeV/#it{c}]; "
               "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
               rapBins[i], rapBins[i + 1]),
          500, 0, 10, 1000, 1., 1.5);
      fHistogramsSigma->Add(fHistSigmaPtY[i]);
      fHistogramsSigma->Add(fHistSigmaMixedPtY[i]);
    }
  }
  fHistograms->Add(fHistogramsSigma);

  if (fHistogramsAntiSigma != nullptr) {
    delete fHistogramsAntiSigma;
    fHistogramsAntiSigma = nullptr;
  }
  if (fHistogramsAntiSigma == nullptr) {
    fHistogramsAntiSigma = new TList();
    fHistogramsAntiSigma->SetOwner(kTRUE);
    name = "AntiSigmaCut_QA" + appendix;
    fHistogramsAntiSigma->SetName(name);
  }

  fHistAntiSigmaInvMassPt =
      new TH2F("fHistAntiSigmaInvMassPt",
               "; #it{p}_{T}  #bar{#Lambda}#gamma (GeV/#it{c}); "
               "M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);

  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassPt);
  fHistAntiSigmaMixedInvMassPt =
      new TH2F("fHistAntiSigmaMixedInvMassPt",
               "; #it{p}_{T} #bar{#Lambda}#gamma (GeV/#it{c}); "
               "M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassPt);

  if (!fIsLightweight) {
    fHistNAntiSigma =
        new TH1F("fHistNAntiSigma", ";# #Sigma candidates; Entries", 20, 0, 20);
    fHistAntiSigmaPt = new TH1F(
        "fHistAntiSigmaPt",
        "; #it{p}_{T} #bar{#Lambda}#gamma (GeV/#it{c}); Entries", 1000, 0, 20);
    fHistAntiSigmaMassCutPt = new TH1F(
        "fHistAntiSigmaMassCutPt",
        "; #it{p}_{T} #bar{#Lambda}#gamma (GeV/#it{c}); Entries", 1000, 0, 20);
    fHistAntiSigmaInvMass = new TH1F(
        "fHistAntiSigmaInvMass",
        "; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
    fHistAntiSigmaInvMassRec = new TH1F(
        "fHistAntiSigmaInvMassRec",
        "; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
    fHistAntiSigmaInvMassEta =
        new TH2F("fHistAntiSigmaInvMassEta",
                 "; #eta; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})", 1000, 0, 1,
                 2000, 1., 2.);
    fHistAntiSigmaArmenterosBefore =
        new TH2F("fHistAntiSigmaArmenterosBefore",
                 " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
    fHistAntiSigmaArmenterosAfter =
        new TH2F("fHistAntiSigmaArmenterosAfter",
                 " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
    fHistAntiSigmaEtaPhi = new TH2F("fHistAntiSigmaEtaPhi", "; #eta; #phi", 200,
                                    -1, 1, 200, 0, 2 * TMath::Pi());
    fHistAntiSigmaRapidity =
        new TH1F("fHistAntiSigmaRapidity", "; y; Entries", 200, -1, 1);
    fHistAntiSigmaMixedPt = new TH1F(
        "fHistAntiSigmaMixedPt",
        "; #it{p}_{T} #bar{#Lambda}#gamma (GeV/#it{c}); Entries", 500, 0, 20);
    fHistAntiSigmaMixedInvMass = new TH1F(
        "fHistAntiSigmaMixedInvMass",
        "; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
    fHistAntiSigmaMixedInvMassEta =
        new TH2F("fHistAntiSigmaMixedInvMassEta",
                 "; #eta; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})", 1000, 0, 1,
                 2000, 1., 2.);

    fHistogramsAntiSigma->Add(fHistNAntiSigma);
    fHistogramsAntiSigma->Add(fHistAntiSigmaPt);
    fHistogramsAntiSigma->Add(fHistAntiSigmaMassCutPt);
    fHistogramsAntiSigma->Add(fHistAntiSigmaInvMass);
    fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassRec);
    fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassEta);
    fHistogramsAntiSigma->Add(fHistAntiSigmaRapidity);
    fHistogramsAntiSigma->Add(fHistAntiSigmaEtaPhi);
    fHistogramsAntiSigma->Add(fHistAntiSigmaArmenterosBefore);
    fHistogramsAntiSigma->Add(fHistAntiSigmaArmenterosAfter);
    fHistogramsAntiSigma->Add(fHistAntiSigmaMixedPt);
    fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMass);
    fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassEta);

    for (int i = 0; i < static_cast<int>(rapBins.size() - 1); ++i) {
      fHistAntiSigmaPtY[i] = new TH2F(
          Form("fHistAntiSigmaPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
          Form("%.2f < y < %.2f ; #it{p}_{T} [GeV/#it{c}]; "
               "M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})",
               rapBins[i], rapBins[i + 1]),
          500, 0, 10, 1000, 1., 1.5);
      fHistAntiSigmaMixedPtY[i] = new TH2F(
          Form("fHistAntiSigmaMixedPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
          Form("%.2f < y < %.2f ; #it{p}_{T} [GeV/#it{c}]; "
               "M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})",
               rapBins[i], rapBins[i + 1]),
          500, 0, 10, 1000, 1., 1.5);
      fHistogramsAntiSigma->Add(fHistAntiSigmaPtY[i]);
      fHistogramsAntiSigma->Add(fHistAntiSigmaMixedPtY[i]);
    }
  }

  fHistograms->Add(fHistogramsAntiSigma);

  if (fIsMC) {
    if (fHistogramsSigmaMC != nullptr) {
      delete fHistogramsSigmaMC;
      fHistogramsSigmaMC = nullptr;
    }
    if (fHistogramsSigmaMC == nullptr) {
      fHistogramsSigmaMC = new TList();
      fHistogramsSigmaMC->SetOwner(kTRUE);
      name = "SigmaCut_MC" + appendix;
      fHistogramsSigmaMC->SetName(name);
    }

    fHistMCTruthSigma0PhotonConvPt =
        new TH1F("fHistMCTruthSigma0PhotonConvPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0PhotonConvP =
        new TH1F("fHistMCTruthSigma0PhotonConvP",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0PhotonConvInvMass =
        new TH1F("fHistMCTruthSigma0PhotonConvInvMass",
                 "; Invariant mass (GeV/#it{c}^{2}); Entries", 500, 0, 1);
    fHistMCTruthSigma0PhotonConvInvMassPt =
        new TH2F("fHistMCTruthSigma0PhotonConvInvMassPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Invariant mass (GeV/#it{c}^{2})",
                 500, 0, 10, 500, 0, 1);
    fHistMCTruthSigma0PhotonConvPtEta =
        new TH2F("fHistMCTruthSigma0PhotonConvInvMassEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistMCTruthSigma0PhotonConvR =
        new TH2F("fHistMCTruthSigma0PhotonConvR",
                 "; #it{p}_{T} [GeV/#it{c}]; Conversion radius", 500, 0, 10,
                 500, 0, 200);
    fHistMCTruthSigma0PhotonConvConvPointX =
        new TH1F("fHistMCTruthSigma0PhotonConvConvPointX",
                 "; Conv. point x [cm]; Entries", 500, -200, 200);
    fHistMCTruthSigma0PhotonConvConvPointY =
        new TH1F("fHistMCTruthSigma0PhotonConvConvPointY",
                 "; Conv. point y [cm]; Entries", 500, -200, 200);
    fHistMCTruthSigma0PhotonConvConvPointZ =
        new TH1F("fHistMCTruthSigma0PhotonConvConvPointZ",
                 "; Conv. point z [cm]; Entries", 500, -200, 200);
    fHistMCTruthSigma0PhotonConvEleP =
        new TH1F("fHistMCTruthSigma0PhotonConvEleP",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0PhotonConvElePt =
        new TH1F("fHistMCTruthSigma0PhotonConvElePt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0PhotonConvPtY =
        new TH2F("fHistMCTruthSigma0PhotonConvPtY",
                 "; y; #it{p}_{T} [GeV/#it{c}]", 1000, -10, 10, 500, 0, 10);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvPt);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvP);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvInvMass);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvInvMassPt);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvPtEta);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvR);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvConvPointX);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvConvPointY);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvConvPointZ);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvEleP);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvElePt);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonConvPtY);

    fHistMCTruthSigma0Pt =
        new TH1F("fHistMCTruthSigma0Pt", "; #it{p}_{T} [GeV/#it{c}]; Entries",
                 500, 0, 10);
    fHistMCTruthSigma0PtY =
        new TH2F("fHistMCTruthSigma0PtY", "; y; #it{p}_{T} [GeV/#it{c}]", 1000,
                 -10, 10, 500, 0, 10);
    fHistMCTruthSigma0PtEta =
        new TH2F("fHistMCTruthSigma0PtEta", "; #eta; #it{p}_{T} [GeV/#it{c}]",
                 500, -10, 10, 500, 0, 10);
    fHistMCTruthSigma0PhotonPt =
        new TH1F("fHistMCTruthSigma0PhotonPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0PhotonPtY =
        new TH2F("fHistMCTruthSigma0PhotonPtY", "; y; #it{p}_{T} [GeV/#it{c}]",
                 1000, -10, 10, 500, 0, 10);
    fHistMCTruthSigma0PhotonPtEta =
        new TH2F("fHistMCTruthSigma0PhotonPtEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0Pt);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PtY);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PtEta);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonPt);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonPtY);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0PhotonPtEta);

    fHistMCSigmaMassCutPt = new TH1F(
        "fHistMCSigmaMassCutPt",
        "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries", 500, 0, 10);
    fHistogramsSigmaMC->Add(fHistMCSigmaMassCutPt);
    fHistogramsSigma->Add(fHistogramsSigmaMC);

    if (fHistogramsAntiSigmaMC != nullptr) {
      delete fHistogramsAntiSigmaMC;
      fHistogramsAntiSigmaMC = nullptr;
    }
    if (fHistogramsAntiSigmaMC == nullptr) {
      fHistogramsAntiSigmaMC = new TList();
      fHistogramsAntiSigmaMC->SetOwner(kTRUE);
      name = "AntiSigmaCut_MC" + appendix;
      fHistogramsAntiSigmaMC->SetName(name);
    }
    fHistMCTruthAntiSigma0PhotonConvPt =
        new TH1F("fHistMCTruthAntiSigma0PhotonConvPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0PhotonConvP =
        new TH1F("fHistMCTruthAntiSigma0PhotonConvP",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0PhotonConvInvMass =
        new TH1F("fHistMCTruthAntiSigma0PhotonConvInvMass",
                 "; Invariant mass (GeV/#it{c}^{2}); Entries", 500, 0, 1);
    fHistMCTruthAntiSigma0PhotonConvInvMassPt =
        new TH2F("fHistMCTruthAntiSigma0PhotonConvInvMassPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Invariant mass (GeV/#it{c}^{2})",
                 500, 0, 10, 500, 0, 1);
    fHistMCTruthAntiSigma0PhotonConvPtEta =
        new TH2F("fHistMCTruthAntiSigma0PhotonConvPtEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistMCTruthAntiSigma0PhotonConvR =
        new TH2F("fHistMCTruthAntiSigma0PhotonConvR",
                 "; #it{p}_{T} [GeV/#it{c}]; Conversion radius", 500, 0, 10,
                 500, 0, 200);
    fHistMCTruthAntiSigma0PhotonConvConvPointX =
        new TH1F("fHistMCTruthAntiSigma0PhotonConvConvPointX",
                 "; Conv. point x [cm]; Entries", 500, -200, 200);
    fHistMCTruthAntiSigma0PhotonConvConvPointY =
        new TH1F("fHistMCTruthAntiSigma0PhotonConvConvPointY",
                 "; Conv. point y [cm]; Entries", 500, -200, 200);
    fHistMCTruthAntiSigma0PhotonConvConvPointZ =
        new TH1F("fHistMCTruthAntiSigma0PhotonConvConvPointZ",
                 "; Conv. point z [cm]; Entries", 500, -200, 200);
    fHistMCTruthAntiSigma0PhotonConvEleP =
        new TH1F("fHistMCTruthAntiSigma0PhotonConvEleP",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0PhotonConvElePt =
        new TH1F("fHistMCTruthAntiSigma0PhotonConvElePt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0PhotonConvPtY =
        new TH2F("fHistMCTruthAntiSigma0PhotonConvPtY",
                 "; y; #it{p}_{T} [GeV/#it{c}]", 1000, -10, 10, 500, 0, 10);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvPt);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvP);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvInvMass);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvInvMassPt);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvPtEta);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvR);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvConvPointX);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvConvPointY);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvConvPointZ);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvEleP);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvElePt);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonConvPtY);

    fHistMCTruthAntiSigma0Pt =
        new TH1F("fHistMCTruthAntiSigma0Pt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0PtY =
        new TH2F("fHistMCTruthAntiSigma0PtY", "; y; #it{p}_{T} [GeV/#it{c}]",
                 1000, -10, 10, 500, 0, 10);
    fHistMCTruthAntiSigma0PtEta =
        new TH2F("fHistMCTruthAntiSigma0PtEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistMCTruthAntiSigma0PhotonPt =
        new TH1F("fHistMCTruthAntiSigma0PhotonPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0PhotonPtY =
        new TH2F("fHistMCTruthAntiSigma0PhotonPtY",
                 "; y; #it{p}_{T} [GeV/#it{c}]", 1000, -10, 10, 500, 0, 10);
    fHistMCTruthAntiSigma0PhotonPtEta =
        new TH2F("fHistMCTruthAntiSigma0PhotonPtEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0Pt);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PtY);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PtEta);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonPt);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonPtY);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0PhotonPtEta);

    fHistMCAntiSigmaMassCutPt = new TH1F(
        "fHistMCAntiSigmaMassCutPt",
        "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries", 500, 0, 10);
    fHistogramsAntiSigmaMC->Add(fHistMCAntiSigmaMassCutPt);
    fHistogramsAntiSigma->Add(fHistogramsAntiSigmaMC);
  }
}
