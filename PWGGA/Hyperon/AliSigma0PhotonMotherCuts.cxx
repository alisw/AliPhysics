#include "AliSigma0PhotonMotherCuts.h"

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliSigma0PhotonCuts.h"

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
      fHistSigmaInvMassFinal(nullptr),
      fHistSigmaInvMassRec(nullptr),
      fHistSigmaInvMassPt(nullptr),
      fHistSigmaInvMassEta(nullptr),
      fHistSigmaEtaPhi(nullptr),
      fHistSigmaRapidity(nullptr),
      fHistSigmaPtY(),
      fHistSigmaNShared(nullptr),
      fHistSigmaMassDiff(nullptr),
      fHistSigmaSharedInvMass(nullptr),
      fHistSigmaArmenterosBefore(nullptr),
      fHistSigmaArmenterosAfter(nullptr),
      fHistSigmaMixedPt(nullptr),
      fHistSigmaMixedInvMass(nullptr),
      fHistSigmaMixedPtY(),
      fHistSigmaMixedInvMassPt(nullptr),
      fHistSigmaMixedInvMassEta(nullptr),
      fHistSigmaPtTwoGamma(nullptr),
      fHistSigmaInvMassTwoGamma(nullptr),
      fHistSigmaInvMassPtTwoGamma(nullptr),
      fHistSigmaInvMassEtaTwoGamma(nullptr),
      fHistSigmaArmenterosBeforeTwoGamma(nullptr),
      fHistSigmaArmenterosAfterTwoGamma(nullptr),
      fHistSigmaMixedPtTwoGamma(nullptr),
      fHistSigmaMixedInvMassTwoGamma(nullptr),
      fHistSigmaMixedInvMassPtTwoGamma(nullptr),
      fHistSigmaMixedInvMassEtaTwoGamma(nullptr),
      fHistSigmaPtDiElectron(nullptr),
      fHistSigmaInvMassDiElectron(nullptr),
      fHistSigmaInvMassPtDiElectron(nullptr),
      fHistSigmaInvMassEtaDiElectron(nullptr),
      fHistSigmaArmenterosBeforeDiElectron(nullptr),
      fHistSigmaArmenterosAfterDiElectron(nullptr),
      fHistSigmaMixedPtDiElectron(nullptr),
      fHistSigmaMixedInvMassDiElectron(nullptr),
      fHistSigmaMixedInvMassPtDiElectron(nullptr),
      fHistSigmaMixedInvMassEtaDiElectron(nullptr),
      fHistNAntiSigma(nullptr),
      fHistAntiSigmaPt(nullptr),
      fHistAntiSigmaMassCutPt(nullptr),
      fHistAntiSigmaInvMass(nullptr),
      fHistAntiSigmaInvMassFinal(nullptr),
      fHistAntiSigmaInvMassRec(nullptr),
      fHistAntiSigmaInvMassPt(nullptr),
      fHistAntiSigmaInvMassEta(nullptr),
      fHistAntiSigmaEtaPhi(nullptr),
      fHistAntiSigmaRapidity(nullptr),
      fHistAntiSigmaPtY(),
      fHistAntiSigmaNShared(nullptr),
      fHistAntiSigmaMassDiff(nullptr),
      fHistAntiSigmaSharedInvMass(nullptr),
      fHistAntiSigmaArmenterosBefore(nullptr),
      fHistAntiSigmaArmenterosAfter(nullptr),
      fHistAntiSigmaMixedPt(nullptr),
      fHistAntiSigmaMixedInvMass(nullptr),
      fHistAntiSigmaMixedPtY(),
      fHistAntiSigmaMixedInvMassPt(nullptr),
      fHistAntiSigmaMixedInvMassEta(nullptr),
      fHistAntiSigmaPtTwoGamma(nullptr),
      fHistAntiSigmaInvMassTwoGamma(nullptr),
      fHistAntiSigmaInvMassPtTwoGamma(nullptr),
      fHistAntiSigmaInvMassEtaTwoGamma(nullptr),
      fHistAntiSigmaArmenterosBeforeTwoGamma(nullptr),
      fHistAntiSigmaArmenterosAfterTwoGamma(nullptr),
      fHistAntiSigmaMixedPtTwoGamma(nullptr),
      fHistAntiSigmaMixedInvMassTwoGamma(nullptr),
      fHistAntiSigmaMixedInvMassPtTwoGamma(nullptr),
      fHistAntiSigmaMixedInvMassEtaTwoGamma(nullptr),
      fHistAntiSigmaPtDiElectron(nullptr),
      fHistAntiSigmaInvMassDiElectron(nullptr),
      fHistAntiSigmaInvMassPtDiElectron(nullptr),
      fHistAntiSigmaInvMassEtaDiElectron(nullptr),
      fHistAntiSigmaArmenterosBeforeDiElectron(nullptr),
      fHistAntiSigmaArmenterosAfterDiElectron(nullptr),
      fHistAntiSigmaMixedPtDiElectron(nullptr),
      fHistAntiSigmaMixedInvMassDiElectron(nullptr),
      fHistAntiSigmaMixedInvMassPtDiElectron(nullptr),
      fHistAntiSigmaMixedInvMassEtaDiElectron(nullptr),
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
      fHistMCTruthSigma0TwoElectronPtPos(nullptr),
      fHistMCTruthSigma0TwoElectronPPos(nullptr),
      fHistMCTruthSigma0TwoElectronPtEtaPos(nullptr),
      fHistMCTruthSigma0TwoElectronPtNeg(nullptr),
      fHistMCTruthSigma0TwoElectronPNeg(nullptr),
      fHistMCTruthSigma0TwoElectronPtEtaNeg(nullptr),
      fHistMCTruthSigma0TwoGammaPt(nullptr),
      fHistMCTruthSigma0TwoGammaP(nullptr),
      fHistMCTruthSigma0TwoGammaPtEta(nullptr),
      fHistMCTruthSigma0TwoGammaPtEle(nullptr),
      fHistMCTruthSigma0TwoGammaPEle(nullptr),
      fHistMCTruthSigma0TwoGammaPtEtaEle(nullptr),
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
      fHistMCTruthAntiSigma0PhotonPtEta(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPtPos(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPPos(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPtEtaPos(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPtNeg(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPNeg(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPtEtaNeg(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPt(nullptr),
      fHistMCTruthAntiSigma0TwoGammaP(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPtEta(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPtEle(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPEle(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPtEtaEle(nullptr) {}

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
      fHistSigmaInvMassFinal(nullptr),
      fHistSigmaInvMassRec(nullptr),
      fHistSigmaInvMassPt(nullptr),
      fHistSigmaInvMassEta(nullptr),
      fHistSigmaEtaPhi(nullptr),
      fHistSigmaRapidity(nullptr),
      fHistSigmaPtY(),
      fHistSigmaNShared(nullptr),
      fHistSigmaMassDiff(nullptr),
      fHistSigmaSharedInvMass(nullptr),
      fHistSigmaArmenterosBefore(nullptr),
      fHistSigmaArmenterosAfter(nullptr),
      fHistSigmaMixedPt(nullptr),
      fHistSigmaMixedInvMass(nullptr),
      fHistSigmaMixedPtY(),
      fHistSigmaMixedInvMassPt(nullptr),
      fHistSigmaMixedInvMassEta(nullptr),
      fHistSigmaPtTwoGamma(nullptr),
      fHistSigmaInvMassTwoGamma(nullptr),
      fHistSigmaInvMassPtTwoGamma(nullptr),
      fHistSigmaInvMassEtaTwoGamma(nullptr),
      fHistSigmaArmenterosBeforeTwoGamma(nullptr),
      fHistSigmaArmenterosAfterTwoGamma(nullptr),
      fHistSigmaMixedPtTwoGamma(nullptr),
      fHistSigmaMixedInvMassTwoGamma(nullptr),
      fHistSigmaMixedInvMassPtTwoGamma(nullptr),
      fHistSigmaMixedInvMassEtaTwoGamma(nullptr),
      fHistSigmaPtDiElectron(nullptr),
      fHistSigmaInvMassDiElectron(nullptr),
      fHistSigmaInvMassPtDiElectron(nullptr),
      fHistSigmaInvMassEtaDiElectron(nullptr),
      fHistSigmaArmenterosBeforeDiElectron(nullptr),
      fHistSigmaArmenterosAfterDiElectron(nullptr),
      fHistSigmaMixedPtDiElectron(nullptr),
      fHistSigmaMixedInvMassDiElectron(nullptr),
      fHistSigmaMixedInvMassPtDiElectron(nullptr),
      fHistSigmaMixedInvMassEtaDiElectron(nullptr),
      fHistNAntiSigma(nullptr),
      fHistAntiSigmaPt(nullptr),
      fHistAntiSigmaMassCutPt(nullptr),
      fHistAntiSigmaInvMass(nullptr),
      fHistAntiSigmaInvMassFinal(nullptr),
      fHistAntiSigmaInvMassRec(nullptr),
      fHistAntiSigmaInvMassPt(nullptr),
      fHistAntiSigmaInvMassEta(nullptr),
      fHistAntiSigmaEtaPhi(nullptr),
      fHistAntiSigmaRapidity(nullptr),
      fHistAntiSigmaPtY(),
      fHistAntiSigmaNShared(nullptr),
      fHistAntiSigmaMassDiff(nullptr),
      fHistAntiSigmaSharedInvMass(nullptr),
      fHistAntiSigmaArmenterosBefore(nullptr),
      fHistAntiSigmaArmenterosAfter(nullptr),
      fHistAntiSigmaMixedPt(nullptr),
      fHistAntiSigmaMixedInvMass(nullptr),
      fHistAntiSigmaMixedPtY(),
      fHistAntiSigmaMixedInvMassPt(nullptr),
      fHistAntiSigmaMixedInvMassEta(nullptr),
      fHistAntiSigmaPtTwoGamma(nullptr),
      fHistAntiSigmaInvMassTwoGamma(nullptr),
      fHistAntiSigmaInvMassPtTwoGamma(nullptr),
      fHistAntiSigmaInvMassEtaTwoGamma(nullptr),
      fHistAntiSigmaArmenterosBeforeTwoGamma(nullptr),
      fHistAntiSigmaArmenterosAfterTwoGamma(nullptr),
      fHistAntiSigmaMixedPtTwoGamma(nullptr),
      fHistAntiSigmaMixedInvMassTwoGamma(nullptr),
      fHistAntiSigmaMixedInvMassPtTwoGamma(nullptr),
      fHistAntiSigmaMixedInvMassEtaTwoGamma(nullptr),
      fHistAntiSigmaPtDiElectron(nullptr),
      fHistAntiSigmaInvMassDiElectron(nullptr),
      fHistAntiSigmaInvMassPtDiElectron(nullptr),
      fHistAntiSigmaInvMassEtaDiElectron(nullptr),
      fHistAntiSigmaArmenterosBeforeDiElectron(nullptr),
      fHistAntiSigmaArmenterosAfterDiElectron(nullptr),
      fHistAntiSigmaMixedPtDiElectron(nullptr),
      fHistAntiSigmaMixedInvMassDiElectron(nullptr),
      fHistAntiSigmaMixedInvMassPtDiElectron(nullptr),
      fHistAntiSigmaMixedInvMassEtaDiElectron(nullptr),
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
      fHistMCTruthSigma0TwoElectronPtPos(nullptr),
      fHistMCTruthSigma0TwoElectronPPos(nullptr),
      fHistMCTruthSigma0TwoElectronPtEtaPos(nullptr),
      fHistMCTruthSigma0TwoElectronPtNeg(nullptr),
      fHistMCTruthSigma0TwoElectronPNeg(nullptr),
      fHistMCTruthSigma0TwoElectronPtEtaNeg(nullptr),
      fHistMCTruthSigma0TwoGammaPt(nullptr),
      fHistMCTruthSigma0TwoGammaP(nullptr),
      fHistMCTruthSigma0TwoGammaPtEta(nullptr),
      fHistMCTruthSigma0TwoGammaPtEle(nullptr),
      fHistMCTruthSigma0TwoGammaPEle(nullptr),
      fHistMCTruthSigma0TwoGammaPtEtaEle(nullptr),
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
      fHistMCTruthAntiSigma0PhotonPtEta(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPtPos(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPPos(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPtEtaPos(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPtNeg(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPNeg(nullptr),
      fHistMCTruthAntiSigma0TwoElectronPtEtaNeg(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPt(nullptr),
      fHistMCTruthAntiSigma0TwoGammaP(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPtEta(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPtEle(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPEle(nullptr),
      fHistMCTruthAntiSigma0TwoGammaPtEtaEle(nullptr) {}

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
    std::vector<AliSigma0ParticleV0> &antiLambdaCandidates,
    std::vector<AliSigma0ParticlePhotonMother> &sigmaCandidates,
    std::vector<AliSigma0ParticlePhotonMother> &antiSigmaCandidates) {
  sigmaCandidates.clear();
  antiSigmaCandidates.clear();

  fMCEvent = mcEvent;

  if (fIsMC) ProcessMC();

  fInputEvent = inputEvent;

  // Particle pairing
  SigmaToLambdaGamma(photonCandidates, lambdaCandidates, antiLambdaCandidates,
                     sigmaCandidates, antiSigmaCandidates);
  SigmaToLambdaGammaGamma(photonCandidates, lambdaCandidates,
                          antiLambdaCandidates);

  FillEventBuffer(photonCandidates, lambdaCandidates,
                  antiLambdaCandidates);

  // Check for pile-up
  for (auto start = sigmaCandidates.begin(); start != sigmaCandidates.end();
       ++start) {
    AliSigma0ParticlePhotonMother *iSigmaCandidate =
        static_cast<AliSigma0ParticlePhotonMother *>(&(*(start)));
    if (!iSigmaCandidate->GetIsUse()) continue;
    const float invMassI = iSigmaCandidate->GetMass();
    if (invMassI > fMassSigma + fSigmaMassCut ||
        invMassI < fMassSigma - fSigmaMassCut)
      continue;
    for (auto match = start + 1; match != sigmaCandidates.end(); ++match) {
      if (start == match) continue;
      auto *jSigmaCandidate =
          static_cast<AliSigma0ParticlePhotonMother *>(&(*(match)));
      if (!iSigmaCandidate->GetIsUse() || !jSigmaCandidate->GetIsUse())
        continue;
      const float invMassJ = jSigmaCandidate->GetMass();
      if (invMassJ > fMassSigma + fSigmaMassCut ||
          invMassJ < fMassSigma - fSigmaMassCut)
        continue;
      fHistSigmaMassDiff->Fill(invMassI - invMassJ);
    }
  }

  for (auto start = antiSigmaCandidates.begin();
       start != antiSigmaCandidates.end(); ++start) {
    AliSigma0ParticlePhotonMother *iSigmaCandidate =
        static_cast<AliSigma0ParticlePhotonMother *>(&(*(start)));
    if (!iSigmaCandidate->GetIsUse()) continue;
    const float invMassI = iSigmaCandidate->GetMass();
    if (invMassI > fMassSigma + fSigmaMassCut ||
        invMassI < fMassSigma - fSigmaMassCut)
      continue;
    for (auto match = start + 1; match != antiSigmaCandidates.end(); ++match) {
      if (start == match) continue;
      auto *jSigmaCandidate =
          static_cast<AliSigma0ParticlePhotonMother *>(&(*(match)));
      if (!iSigmaCandidate->GetIsUse() || !jSigmaCandidate->GetIsUse())
        continue;
      const float invMassJ = jSigmaCandidate->GetMass();
      if (invMassJ > fMassSigma + fSigmaMassCut ||
          invMassJ < fMassSigma - fSigmaMassCut)
        continue;
      fHistAntiSigmaMassDiff->Fill(invMassI - invMassJ);
    }
  }

  // Track cleaning for Femto
  TrackCleaner(sigmaCandidates, antiSigmaCandidates, fMassSigma - fSigmaMassCut,
               fMassSigma + fSigmaMassCut);
  TrackCleaner(sigmaCandidates, antiSigmaCandidates,
               fMassSigma + fSigmaSidebandLow, fMassSigma + fSigmaSidebandUp);
  TrackCleaner(sigmaCandidates, antiSigmaCandidates,
               fMassSigma - fSigmaSidebandUp, fMassSigma - fSigmaSidebandLow);

  // Book-keeping after track cleaner
  for (const auto &sigma : sigmaCandidates) {
    if (!sigma.GetIsUse()) continue;
    fHistSigmaInvMassFinal->Fill(sigma.GetMass());
  }

  for (const auto &sigma : antiSigmaCandidates) {
    if (!sigma.GetIsUse()) continue;
    fHistAntiSigmaInvMassFinal->Fill(sigma.GetMass());
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGamma(
    std::vector<AliAODConversionPhoton> &photonCandidates,
    std::vector<AliSigma0ParticleV0> &lambdaCandidates,
    std::vector<AliSigma0ParticleV0> &antiLambdaCandidates,
    std::vector<AliSigma0ParticlePhotonMother> &sigmaCandidates,
    std::vector<AliSigma0ParticlePhotonMother> &antiSigmaCandidates) {
  // SAME EVENT
  for (const auto &photon : photonCandidates) {
    int nAntiSigma = 0;
    for (const auto &antiLambda : antiLambdaCandidates) {
      /// Candidates with lambdas with shared daughter tracks are not used
      if (!antiLambda.GetIsUse()) continue;
      auto *antiSigma =
          new AliSigma0ParticlePhotonMother(antiLambda, photon, fInputEvent);
      const float armAlpha = antiSigma->GetArmenterosAlpha();
      const float armQt = antiSigma->GetArmenterosQt();
      fHistAntiSigmaArmenterosBefore->Fill(armAlpha, armQt);
      const float invMass = antiSigma->GetMass();
      // Armenteros cut
      if (fArmenterosCut) {
        if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
          delete antiSigma;
          continue;
        }
        if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow) {
          delete antiSigma;
          continue;
        }
      }
      fHistAntiSigmaArmenterosAfter->Fill(armAlpha, armQt);
      antiSigmaCandidates.push_back(*antiSigma);
      fHistAntiSigmaPt->Fill(antiSigma->GetPt());
      fHistAntiSigmaInvMass->Fill(invMass);
      fHistAntiSigmaInvMassRec->Fill(antiSigma->GetRecMass());
      fHistAntiSigmaInvMassPt->Fill(antiSigma->GetPt(), invMass);
      const float rap =
          ComputeRapidity(antiSigma->GetPt(), antiSigma->GetPz(), invMass);
      fHistAntiSigmaRapidity->Fill(rap);
      const int rapBin = GetRapidityBin(rap);
      if (rapBin > -1)
        fHistAntiSigmaPtY[rapBin]->Fill(antiSigma->GetPt(), invMass);
      fHistAntiSigmaInvMassEta->Fill(antiSigma->GetEta(), invMass);
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        fHistAntiSigmaMassCutPt->Fill(antiSigma->GetPt());
        fHistAntiSigmaEtaPhi->Fill(antiSigma->GetEta(), antiSigma->GetPhi());
        ++nAntiSigma;
        if (fIsMC) {
          if (antiSigma->IsTrueSigma(fMCEvent)) {
            fHistMCAntiSigmaMassCutPt->Fill(antiSigma->GetPt());
          }
        }
      }
      delete antiSigma;
    }
    fHistNAntiSigma->Fill(nAntiSigma);

    int nSigma = 0;
    // Sigmas
    for (const auto &lambda : lambdaCandidates) {
      /// Candidates with lambdas with shared daughter tracks are not used
      if (!lambda.GetIsUse()) continue;
      auto *sigma =
          new AliSigma0ParticlePhotonMother(lambda, photon, fInputEvent);
      const float armAlpha = sigma->GetArmenterosAlpha();
      const float armQt = sigma->GetArmenterosQt();
      fHistSigmaArmenterosBefore->Fill(armAlpha, armQt);
      const float invMass = sigma->GetMass();
      // Armenteros cut
      if (fArmenterosCut) {
        if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
          delete sigma;
          continue;
        }
        if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow) {
          delete sigma;
          continue;
        }
      }
      fHistSigmaArmenterosAfter->Fill(armAlpha, armQt);
      sigmaCandidates.push_back(*sigma);
      fHistSigmaPt->Fill(sigma->GetPt());
      fHistSigmaInvMass->Fill(invMass);
      fHistSigmaInvMassRec->Fill(sigma->GetRecMass());
      fHistSigmaInvMassPt->Fill(sigma->GetPt(), invMass);
      const float rap =
          ComputeRapidity(sigma->GetPt(), sigma->GetPz(), invMass);
      fHistSigmaRapidity->Fill(rap);
      const int rapBin = GetRapidityBin(rap);
      if (rapBin > -1) fHistSigmaPtY[rapBin]->Fill(sigma->GetPt(), invMass);
      fHistSigmaInvMassEta->Fill(sigma->GetEta(), invMass);
      if (invMass < fMassSigma + fSigmaMassCut &&
          invMass > fMassSigma - fSigmaMassCut) {
        fHistSigmaMassCutPt->Fill(sigma->GetPt());
        fHistSigmaEtaPhi->Fill(sigma->GetEta(), sigma->GetPhi());
        ++nSigma;
        if (fIsMC) {
          if (sigma->IsTrueSigma(fMCEvent)) {
            fHistMCSigmaMassCutPt->Fill(sigma->GetPt());
          }
        }
      }
      delete sigma;
    }
    fHistNSigma->Fill(nSigma);
  }

  // MIXED EVENT

  // Anti - Sigma0
  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fAntiLambdaMixed) {
    for (const auto &Lambda : LambdaContainer) {
      if (!Lambda.GetIsUse()) continue;
      for (auto Photon : photonCandidates) {
        AliSigma0ParticlePhotonMother *sigma =
            new AliSigma0ParticlePhotonMother(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma->GetArmenterosAlpha();
        const float armQt = sigma->GetArmenterosQt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
            delete sigma;
            continue;
          }
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow) {
            delete sigma;
            continue;
          }
        }
        const float invMass = sigma->GetMass();
        fHistAntiSigmaMixedPt->Fill(sigma->GetPt());
        fHistAntiSigmaMixedInvMass->Fill(invMass);
        const float rap =
            ComputeRapidity(sigma->GetPt(), sigma->GetPz(), invMass);
        const int rapBin = GetRapidityBin(rap);
        if (rapBin > -1)
          fHistAntiSigmaMixedPtY[rapBin]->Fill(sigma->GetPt(), invMass);
        fHistAntiSigmaMixedInvMassPt->Fill(sigma->GetPt(), invMass);
        fHistAntiSigmaMixedInvMassEta->Fill(sigma->GetEta(), invMass);
        delete sigma;
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto &PhotonContainer : fPhotonMixed) {
    for (const auto &Photon : PhotonContainer) {
      for (const auto &Lambda : antiLambdaCandidates) {
        if (!Lambda.GetIsUse()) continue;
        AliSigma0ParticlePhotonMother *sigma =
            new AliSigma0ParticlePhotonMother(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma->GetArmenterosAlpha();
        const float armQt = sigma->GetArmenterosQt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
            delete sigma;
            continue;
          }
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow) {
            delete sigma;
            continue;
          }
        }
        const float invMass = sigma->GetMass();
        fHistAntiSigmaMixedPt->Fill(sigma->GetPt());
        fHistAntiSigmaMixedInvMass->Fill(invMass);
        const float rap =
            ComputeRapidity(sigma->GetPt(), sigma->GetPz(), invMass);
        const int rapBin = GetRapidityBin(rap);
        if (rapBin > -1)
          fHistAntiSigmaMixedPtY[rapBin]->Fill(sigma->GetPt(), invMass);
        fHistAntiSigmaMixedInvMassPt->Fill(sigma->GetPt(), invMass);
        fHistAntiSigmaMixedInvMassEta->Fill(sigma->GetEta(), invMass);
        delete sigma;
      }
    }
  }

  // Sigma0
  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fLambdaMixed) {
    for (const auto &Lambda : LambdaContainer) {
      if (!Lambda.GetIsUse()) continue;
      for (auto Photon : photonCandidates) {
        AliSigma0ParticlePhotonMother *sigma =
            new AliSigma0ParticlePhotonMother(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma->GetArmenterosAlpha();
        const float armQt = sigma->GetArmenterosQt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
            delete sigma;
            continue;
          }
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow) {
            delete sigma;
            continue;
          }
        }
        const float invMass = sigma->GetMass();
        fHistSigmaMixedPt->Fill(sigma->GetPt());
        fHistSigmaMixedInvMass->Fill(invMass);
        const float rap =
            ComputeRapidity(sigma->GetPt(), sigma->GetPz(), invMass);
        const int rapBin = GetRapidityBin(rap);
        if (rapBin > -1)
          fHistSigmaMixedPtY[rapBin]->Fill(sigma->GetPt(), invMass);
        fHistSigmaMixedInvMassPt->Fill(sigma->GetPt(), invMass);
        fHistSigmaMixedInvMassEta->Fill(sigma->GetEta(), invMass);
        delete sigma;
      }
    }
  }

  // lambdas from this event with mixed photons
  for (const auto &PhotonContainer : fPhotonMixed) {
    for (const auto &Photon : PhotonContainer) {
      for (const auto &Lambda : lambdaCandidates) {
        if (!Lambda.GetIsUse()) continue;
        AliSigma0ParticlePhotonMother *sigma =
            new AliSigma0ParticlePhotonMother(Lambda, Photon, fInputEvent);
        // Armenteros cut
        const float armAlpha = sigma->GetArmenterosAlpha();
        const float armQt = sigma->GetArmenterosQt();
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
            delete sigma;
            continue;
          }
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow) {
            delete sigma;
            continue;
          }
        }
        const float invMass = sigma->GetMass();
        fHistSigmaMixedPt->Fill(sigma->GetPt());
        fHistSigmaMixedInvMass->Fill(invMass);
        const float rap =
            ComputeRapidity(sigma->GetPt(), sigma->GetPz(), invMass);
        const int rapBin = GetRapidityBin(rap);
        if (rapBin > -1)
          fHistSigmaMixedPtY[rapBin]->Fill(sigma->GetPt(), invMass);
        fHistSigmaMixedInvMassPt->Fill(sigma->GetPt(), invMass);
        fHistSigmaMixedInvMassEta->Fill(sigma->GetEta(), invMass);
        delete sigma;
      }
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::SigmaToLambdaGammaGamma(
    std::vector<AliAODConversionPhoton> &photonCandidates,
    std::vector<AliSigma0ParticleV0> &lambdaCandidates,
    std::vector<AliSigma0ParticleV0> &antiLambdaCandidates) {
  // SAME EVENT
  for (auto photon1 = photonCandidates.begin();
       photon1 != photonCandidates.end(); ++photon1) {
    for (auto photon2 = photon1 + 1; photon2 != photonCandidates.end();
         ++photon2) {
      // Anti - Sigma0
      for (const auto &antiLambda : antiLambdaCandidates) {
        /// Candidates with lambdas with shared daughter tracks are not used
        if (!antiLambda.GetIsUse()) continue;
        auto *antiSigma = new AliSigma0ParticlePhotonMother(
            antiLambda, *photon1, *photon2, fInputEvent);
        const float armAlpha = antiSigma->GetArmenterosAlpha();
        const float armQt = antiSigma->GetArmenterosQt();
        fHistAntiSigmaArmenterosBeforeTwoGamma->Fill(armAlpha, armQt);
        // Armenteros cut
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
            delete antiSigma;
            continue;
          }
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow) {
            delete antiSigma;
            continue;
          }
        }
        fHistAntiSigmaArmenterosAfterTwoGamma->Fill(armAlpha, armQt);
        const float invMass = antiSigma->GetMass();
        fHistAntiSigmaPtTwoGamma->Fill(antiSigma->GetPt());
        fHistAntiSigmaInvMassTwoGamma->Fill(invMass);
        fHistAntiSigmaInvMassPtTwoGamma->Fill(antiSigma->GetPt(), invMass);
        fHistAntiSigmaInvMassEtaTwoGamma->Fill(antiSigma->GetEta(), invMass);
        delete antiSigma;
      }

      // Sigma0
      for (const auto &lambda : lambdaCandidates) {
        /// Candidates with lambdas with shared daughter tracks are not used
        if (!lambda.GetIsUse()) continue;
        auto *sigma = new AliSigma0ParticlePhotonMother(lambda, *photon1,
                                                        *photon2, fInputEvent);
        const float armAlpha = sigma->GetArmenterosAlpha();
        const float armQt = sigma->GetArmenterosQt();
        fHistSigmaArmenterosBeforeTwoGamma->Fill(armAlpha, armQt);
        // Armenteros cut
        if (fArmenterosCut) {
          if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
            delete sigma;
            continue;
          }
          if (armAlpha > fArmenterosAlphaUp || armAlpha < fArmenterosAlphaLow) {
            delete sigma;
            continue;
          }
        }
        fHistSigmaArmenterosAfterTwoGamma->Fill(armAlpha, armQt);
        float invMass = sigma->GetMass();
        fHistSigmaPtTwoGamma->Fill(sigma->GetPt());
        fHistSigmaInvMassTwoGamma->Fill(invMass);
        fHistSigmaInvMassPtTwoGamma->Fill(sigma->GetPt(), invMass);
        fHistSigmaInvMassEtaTwoGamma->Fill(sigma->GetEta(), invMass);
        delete sigma;
      }
    }
  }

  // MIXED EVENT

  // Anti - Sigma0
  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fAntiLambdaMixed) {
    for (const auto &Lambda : LambdaContainer) {
      if (!Lambda.GetIsUse()) continue;
      for (auto photon1 = photonCandidates.begin();
           photon1 != photonCandidates.end(); ++photon1) {
        for (auto photon2 = photon1 + 1; photon2 != photonCandidates.end();
             ++photon2) {
          AliSigma0ParticlePhotonMother *sigma =
              new AliSigma0ParticlePhotonMother(Lambda, *photon1, *photon2,
                                                fInputEvent);
          // Armenteros cut
          const float armAlpha = sigma->GetArmenterosAlpha();
          const float armQt = sigma->GetArmenterosQt();
          if (fArmenterosCut) {
            if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
              delete sigma;
              continue;
            }
            if (armAlpha > fArmenterosAlphaUp ||
                armAlpha < fArmenterosAlphaLow) {
              delete sigma;
              continue;
            }
          }
          const float invMass = sigma->GetMass();
          fHistAntiSigmaMixedPt->Fill(sigma->GetPt());
          fHistAntiSigmaMixedInvMass->Fill(invMass);
          fHistAntiSigmaMixedInvMassPt->Fill(sigma->GetPt(), invMass);
          fHistAntiSigmaMixedInvMassEta->Fill(sigma->GetEta(), invMass);
          delete sigma;
        }
      }
    }
  }

  // lambdas from this event with one mixed photon and one photon from this
  // event
  for (const auto &PhotonContainer : fPhotonMixed) {
    for (const auto &photon1 : PhotonContainer) {
      for (const auto &photon2 : photonCandidates) {
        for (const auto &Lambda : antiLambdaCandidates) {
          if (!Lambda.GetIsUse()) continue;
          AliSigma0ParticlePhotonMother *sigma =
              new AliSigma0ParticlePhotonMother(Lambda, photon1, photon2,
                                                fInputEvent);
          // Armenteros cut
          const float armAlpha = sigma->GetArmenterosAlpha();
          const float armQt = sigma->GetArmenterosQt();
          if (fArmenterosCut) {
            if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
              delete sigma;
              continue;
            }
            if (armAlpha > fArmenterosAlphaUp ||
                armAlpha < fArmenterosAlphaLow) {
              delete sigma;
              continue;
            }
          }
          const float invMass = sigma->GetMass();
          fHistAntiSigmaMixedPt->Fill(sigma->GetPt());
          fHistAntiSigmaMixedInvMass->Fill(invMass);
          fHistAntiSigmaMixedInvMassPt->Fill(sigma->GetPt(), invMass);
          fHistAntiSigmaMixedInvMassEta->Fill(sigma->GetEta(), invMass);
          delete sigma;
        }
      }
    }
  }

  // Sigma0
  // photons from this event with mixed lambdas
  for (const auto &LambdaContainer : fLambdaMixed) {
    for (const auto &Lambda : LambdaContainer) {
      if (!Lambda.GetIsUse()) continue;
      for (auto photon1 = photonCandidates.begin();
           photon1 != photonCandidates.end(); ++photon1) {
        for (auto photon2 = photon1 + 1; photon2 != photonCandidates.end();
             ++photon2) {
          AliSigma0ParticlePhotonMother *sigma =
              new AliSigma0ParticlePhotonMother(Lambda, *photon1, *photon2,
                                                fInputEvent);
          // Armenteros cut
          const float armAlpha = sigma->GetArmenterosAlpha();
          const float armQt = sigma->GetArmenterosQt();
          if (fArmenterosCut) {
            if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
              delete sigma;
              continue;
            }
            if (armAlpha > fArmenterosAlphaUp ||
                armAlpha < fArmenterosAlphaLow) {
              delete sigma;
              continue;
            }
          }
          const float invMass = sigma->GetMass();
          fHistSigmaMixedPtTwoGamma->Fill(sigma->GetPt());
          fHistSigmaMixedInvMassTwoGamma->Fill(invMass);
          fHistSigmaMixedInvMassPtTwoGamma->Fill(sigma->GetPt(), invMass);
          fHistSigmaMixedInvMassEtaTwoGamma->Fill(sigma->GetEta(), invMass);
          delete sigma;
        }
      }
    }
  }

  // lambdas from this event with one mixed photon and one photon from this
  // event
  for (const auto &PhotonContainer : fPhotonMixed) {
    for (const auto &photon1 : PhotonContainer) {
      for (const auto &photon2 : photonCandidates) {
        for (const auto &Lambda : lambdaCandidates) {
          if (!Lambda.GetIsUse()) continue;
          AliSigma0ParticlePhotonMother *sigma =
              new AliSigma0ParticlePhotonMother(Lambda, photon1, photon2,
                                                fInputEvent);
          // Armenteros cut
          const float armAlpha = sigma->GetArmenterosAlpha();
          const float armQt = sigma->GetArmenterosQt();
          if (fArmenterosCut) {
            if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) {
              delete sigma;
              continue;
            }
            if (armAlpha > fArmenterosAlphaUp ||
                armAlpha < fArmenterosAlphaLow) {
              delete sigma;
              continue;
            }
          }
          const float invMass = sigma->GetMass();
          fHistSigmaMixedPtTwoGamma->Fill(sigma->GetPt());
          fHistSigmaMixedInvMassTwoGamma->Fill(invMass);
          fHistSigmaMixedInvMassPtTwoGamma->Fill(sigma->GetPt(), invMass);
          fHistSigmaMixedInvMassEtaTwoGamma->Fill(sigma->GetEta(), invMass);
          delete sigma;
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
void AliSigma0PhotonMotherCuts::TrackCleaner(
    std::vector<AliSigma0ParticlePhotonMother> &sigmaCandidates,
    std::vector<AliSigma0ParticlePhotonMother> &antiSigmaCandidates,
    float massLow, float massUp) {
  /// Remove Sigma0 candidates with shared tracks
  /// Only tracks in the proper mass interval are considered
  /// Compare the track labels of the photon (+/+, -/-, +/-, -/+) and of the
  /// lambda (+/+, -/-, +/-, -/+) and mixed photon lambda
  /// o If the lambda shares tracks, the Sigma0 candidate with the Lambda with
  /// worse cosPA is removed
  /// o If the photon shares tracks, the Sigma0 candidate with the photon with
  /// worse cosPA is removed
  /// o If both daughters have shared tracks, the Sigma0 candidates with the
  /// overall worse cosPA is removed
  /// Track sharing of lambda decay protons with the singleTrack container is
  /// circumvented by removing those lambdas before

  // +++++++++++++++++++++++++++++++++++++++++++
  // clean SIGMA - SIGMA
  int nShared = 0;
  for (auto start = sigmaCandidates.begin(); start != sigmaCandidates.end();
       ++start) {
    AliSigma0ParticlePhotonMother *iSigmaCandidate =
        static_cast<AliSigma0ParticlePhotonMother *>(&(*(start)));
    if (!iSigmaCandidate->GetIsUse()) continue;
    if (iSigmaCandidate->GetMass() > massUp ||
        iSigmaCandidate->GetMass() < massLow)
      continue;
    auto iSigmaLambda = iSigmaCandidate->GetV0();
    auto iSigmaPhoton = iSigmaCandidate->GetPhoton();
    fHistSigmaSharedInvMass->Fill(0.f, iSigmaCandidate->GetMass());
    /// Candidates with shared tracks between daughters are immediately not used
    if (iSigmaLambda.GetTrackLabelNeg() == iSigmaPhoton.GetTrackLabelNeg() ||
        iSigmaLambda.GetTrackLabelNeg() == iSigmaPhoton.GetTrackLabelPos() ||
        iSigmaLambda.GetTrackLabelPos() == iSigmaPhoton.GetTrackLabelPos() ||
        iSigmaLambda.GetTrackLabelPos() == iSigmaPhoton.GetTrackLabelNeg()) {
      fHistSigmaSharedInvMass->Fill(1.f, iSigmaCandidate->GetMass());
      iSigmaCandidate->SetUse(false);
    }
    for (auto match = start + 1; match != sigmaCandidates.end(); ++match) {
      if (start == match) continue;
      auto *jSigmaCandidate =
          static_cast<AliSigma0ParticlePhotonMother *>(&(*(match)));
      if (!iSigmaCandidate->GetIsUse() || !jSigmaCandidate->GetIsUse())
        continue;
      if (jSigmaCandidate->GetMass() > massUp ||
          jSigmaCandidate->GetMass() < massLow)
        continue;
      auto jSigmaLambda = jSigmaCandidate->GetV0();
      auto jSigmaPhoton = jSigmaCandidate->GetPhoton();

      // Check for daughter track labels
      if (iSigmaLambda.GetTrackLabelNeg() != jSigmaLambda.GetTrackLabelNeg() &&
          iSigmaLambda.GetTrackLabelNeg() != jSigmaLambda.GetTrackLabelPos() &&
          iSigmaLambda.GetTrackLabelPos() != jSigmaLambda.GetTrackLabelNeg() &&
          iSigmaLambda.GetTrackLabelPos() != jSigmaLambda.GetTrackLabelPos() &&
          iSigmaLambda.GetTrackLabelNeg() != jSigmaPhoton.GetTrackLabelNeg() &&
          iSigmaLambda.GetTrackLabelNeg() != jSigmaPhoton.GetTrackLabelPos() &&
          iSigmaLambda.GetTrackLabelPos() != jSigmaPhoton.GetTrackLabelNeg() &&
          iSigmaLambda.GetTrackLabelPos() != jSigmaPhoton.GetTrackLabelPos() &&
          iSigmaPhoton.GetTrackLabelNeg() != jSigmaPhoton.GetTrackLabelNeg() &&
          iSigmaPhoton.GetTrackLabelNeg() != jSigmaPhoton.GetTrackLabelPos() &&
          iSigmaPhoton.GetTrackLabelPos() != jSigmaPhoton.GetTrackLabelNeg() &&
          iSigmaPhoton.GetTrackLabelPos() != jSigmaPhoton.GetTrackLabelPos() &&
          iSigmaPhoton.GetTrackLabelNeg() != jSigmaLambda.GetTrackLabelNeg() &&
          iSigmaPhoton.GetTrackLabelNeg() != jSigmaLambda.GetTrackLabelPos() &&
          iSigmaPhoton.GetTrackLabelPos() != jSigmaLambda.GetTrackLabelNeg() &&
          iSigmaPhoton.GetTrackLabelPos() != jSigmaLambda.GetTrackLabelPos())
        continue;
      else {
        ++nShared;
        double sharedCase = -1;

        float iLambdaCosPA = 0.f;
        float iPhotonCosPA = 0.f;
        float jLambdaCosPA = 0.f;
        float jPhotonCosPA = 0.f;

        // lambdas share daughter, but photon does not - take Sigma candidate
        // with lambda with better CPA
        if ((iSigmaLambda.GetTrackLabelNeg() ==
                 jSigmaLambda.GetTrackLabelNeg() ||
             iSigmaLambda.GetTrackLabelNeg() ==
                 jSigmaLambda.GetTrackLabelPos() ||
             iSigmaLambda.GetTrackLabelPos() ==
                 jSigmaLambda.GetTrackLabelNeg() ||
             iSigmaLambda.GetTrackLabelPos() ==
                 jSigmaLambda.GetTrackLabelPos()) &&
            (iSigmaPhoton.GetTrackLabelNeg() !=
                 jSigmaPhoton.GetTrackLabelNeg() &&
             iSigmaPhoton.GetTrackLabelNeg() !=
                 jSigmaPhoton.GetTrackLabelPos() &&
             iSigmaPhoton.GetTrackLabelPos() !=
                 jSigmaPhoton.GetTrackLabelNeg() &&
             iSigmaPhoton.GetTrackLabelPos() !=
                 jSigmaPhoton.GetTrackLabelPos())) {
          iLambdaCosPA = iSigmaLambda.GetCosineAlpha();
          jLambdaCosPA = jSigmaLambda.GetCosineAlpha();
          sharedCase = 2;
        }
        // photons share daughter, but lambda does not - take Sigma candidate
        // with photon with better CPA
        else if ((iSigmaPhoton.GetTrackLabelNeg() ==
                      jSigmaPhoton.GetTrackLabelNeg() ||
                  iSigmaPhoton.GetTrackLabelNeg() ==
                      jSigmaPhoton.GetTrackLabelPos() ||
                  iSigmaPhoton.GetTrackLabelPos() ==
                      jSigmaPhoton.GetTrackLabelNeg() ||
                  iSigmaPhoton.GetTrackLabelPos() ==
                      jSigmaPhoton.GetTrackLabelPos()) &&
                 (iSigmaLambda.GetTrackLabelNeg() !=
                      jSigmaLambda.GetTrackLabelNeg() &&
                  iSigmaLambda.GetTrackLabelNeg() !=
                      jSigmaLambda.GetTrackLabelPos() &&
                  iSigmaLambda.GetTrackLabelPos() !=
                      jSigmaLambda.GetTrackLabelNeg() &&
                  iSigmaLambda.GetTrackLabelPos() !=
                      jSigmaLambda.GetTrackLabelPos())) {
          iPhotonCosPA = iSigmaPhoton.GetCosineAlpha();
          jPhotonCosPA = jSigmaPhoton.GetCosineAlpha();
          sharedCase = 3;
        }
        // both share daughter - take the sqrt of the respective CPAs
        else if ((iSigmaPhoton.GetTrackLabelNeg() ==
                      jSigmaPhoton.GetTrackLabelNeg() ||
                  iSigmaPhoton.GetTrackLabelNeg() ==
                      jSigmaPhoton.GetTrackLabelPos() ||
                  iSigmaPhoton.GetTrackLabelPos() ==
                      jSigmaPhoton.GetTrackLabelNeg() ||
                  iSigmaPhoton.GetTrackLabelPos() ==
                      jSigmaPhoton.GetTrackLabelPos()) &&
                 (iSigmaLambda.GetTrackLabelNeg() ==
                      jSigmaLambda.GetTrackLabelNeg() ||
                  iSigmaLambda.GetTrackLabelNeg() ==
                      jSigmaLambda.GetTrackLabelPos() ||
                  iSigmaLambda.GetTrackLabelPos() ==
                      jSigmaLambda.GetTrackLabelNeg() ||
                  iSigmaLambda.GetTrackLabelPos() ==
                      jSigmaLambda.GetTrackLabelPos())) {
          iLambdaCosPA = iSigmaLambda.GetCosineAlpha();
          jLambdaCosPA = jSigmaLambda.GetCosineAlpha();
          iPhotonCosPA = iSigmaPhoton.GetCosineAlpha();
          jPhotonCosPA = jSigmaPhoton.GetCosineAlpha();
          sharedCase = 4;
        }
        // mixed track sharing between any lambda & photon
        else {
          iLambdaCosPA = iSigmaLambda.GetCosineAlpha();
          jLambdaCosPA = jSigmaLambda.GetCosineAlpha();
          iPhotonCosPA = iSigmaPhoton.GetCosineAlpha();
          jPhotonCosPA = jSigmaPhoton.GetCosineAlpha();
          sharedCase = 5;
        }

        const float iSigmaSq = std::sqrt(iLambdaCosPA * iLambdaCosPA +
                                         iPhotonCosPA * iPhotonCosPA);
        const float jSigmaSq = std::sqrt(jLambdaCosPA * jLambdaCosPA +
                                         jPhotonCosPA * jPhotonCosPA);
        if (iSigmaSq >= jSigmaSq) {
          // keep i, remove j
          jSigmaCandidate->SetUse(false);
          fHistSigmaSharedInvMass->Fill(sharedCase, jSigmaCandidate->GetMass());
        } else {
          // remove i, keep j
          iSigmaCandidate->SetUse(false);
          fHistSigmaSharedInvMass->Fill(sharedCase, iSigmaCandidate->GetMass());
          break;
        }
      }
    }
  }
  fHistSigmaNShared->Fill(nShared);

  // +++++++++++++++++++++++++++++++++++++++++++
  // clean ANTI-SIGMA - ANTI-SIGMA
  nShared = 0;
  for (auto start = antiSigmaCandidates.begin();
       start != antiSigmaCandidates.end(); ++start) {
    auto *iSigmaCandidate =
        static_cast<AliSigma0ParticlePhotonMother *>(&(*(start)));
    if (!iSigmaCandidate->GetIsUse()) continue;
    if (iSigmaCandidate->GetMass() > massUp ||
        iSigmaCandidate->GetMass() < massLow)
      continue;
    auto iSigmaLambda = iSigmaCandidate->GetV0();
    auto iSigmaPhoton = iSigmaCandidate->GetPhoton();
    fHistAntiSigmaSharedInvMass->Fill(0.f, iSigmaCandidate->GetMass());
    /// Candidates with shared tracks between daughters are immediately not used
    if (iSigmaLambda.GetTrackLabelNeg() == iSigmaPhoton.GetTrackLabelNeg() ||
        iSigmaLambda.GetTrackLabelNeg() == iSigmaPhoton.GetTrackLabelPos() ||
        iSigmaLambda.GetTrackLabelPos() == iSigmaPhoton.GetTrackLabelPos() ||
        iSigmaLambda.GetTrackLabelPos() == iSigmaPhoton.GetTrackLabelNeg()) {
      fHistAntiSigmaSharedInvMass->Fill(1.f, iSigmaCandidate->GetMass());
      iSigmaCandidate->SetUse(false);
    }
    for (auto match = start + 1; match != antiSigmaCandidates.end(); ++match) {
      if (start == match) continue;
      auto *jSigmaCandidate =
          static_cast<AliSigma0ParticlePhotonMother *>(&(*(match)));
      if (!iSigmaCandidate->GetIsUse() || !jSigmaCandidate->GetIsUse())
        continue;
      if (jSigmaCandidate->GetMass() > massUp ||
          jSigmaCandidate->GetMass() < massLow)
        continue;
      auto jSigmaLambda = jSigmaCandidate->GetV0();
      auto jSigmaPhoton = jSigmaCandidate->GetPhoton();

      // Check for daughter track labels
      if (iSigmaLambda.GetTrackLabelNeg() != jSigmaLambda.GetTrackLabelNeg() &&
          iSigmaLambda.GetTrackLabelNeg() != jSigmaLambda.GetTrackLabelPos() &&
          iSigmaLambda.GetTrackLabelPos() != jSigmaLambda.GetTrackLabelNeg() &&
          iSigmaLambda.GetTrackLabelPos() != jSigmaLambda.GetTrackLabelPos() &&
          iSigmaLambda.GetTrackLabelNeg() != jSigmaPhoton.GetTrackLabelNeg() &&
          iSigmaLambda.GetTrackLabelNeg() != jSigmaPhoton.GetTrackLabelPos() &&
          iSigmaLambda.GetTrackLabelPos() != jSigmaPhoton.GetTrackLabelNeg() &&
          iSigmaLambda.GetTrackLabelPos() != jSigmaPhoton.GetTrackLabelPos() &&
          iSigmaPhoton.GetTrackLabelNeg() != jSigmaPhoton.GetTrackLabelNeg() &&
          iSigmaPhoton.GetTrackLabelNeg() != jSigmaPhoton.GetTrackLabelPos() &&
          iSigmaPhoton.GetTrackLabelPos() != jSigmaPhoton.GetTrackLabelNeg() &&
          iSigmaPhoton.GetTrackLabelPos() != jSigmaPhoton.GetTrackLabelPos() &&
          iSigmaPhoton.GetTrackLabelNeg() != jSigmaLambda.GetTrackLabelNeg() &&
          iSigmaPhoton.GetTrackLabelNeg() != jSigmaLambda.GetTrackLabelPos() &&
          iSigmaPhoton.GetTrackLabelPos() != jSigmaLambda.GetTrackLabelNeg() &&
          iSigmaPhoton.GetTrackLabelPos() != jSigmaLambda.GetTrackLabelPos())
        continue;
      else {
        ++nShared;
        double sharedCase = -1;

        float iLambdaCosPA = 0.f;
        float iPhotonCosPA = 0.f;
        float jLambdaCosPA = 0.f;
        float jPhotonCosPA = 0.f;

        // lambdas share daughter, but photon does not - take Sigma candidate
        // with lambda with better CPA
        if ((iSigmaLambda.GetTrackLabelNeg() ==
                 jSigmaLambda.GetTrackLabelNeg() ||
             iSigmaLambda.GetTrackLabelNeg() ==
                 jSigmaLambda.GetTrackLabelPos() ||
             iSigmaLambda.GetTrackLabelPos() ==
                 jSigmaLambda.GetTrackLabelNeg() ||
             iSigmaLambda.GetTrackLabelPos() ==
                 jSigmaLambda.GetTrackLabelPos()) &&
            (iSigmaPhoton.GetTrackLabelNeg() !=
                 jSigmaPhoton.GetTrackLabelNeg() &&
             iSigmaPhoton.GetTrackLabelNeg() !=
                 jSigmaPhoton.GetTrackLabelPos() &&
             iSigmaPhoton.GetTrackLabelPos() !=
                 jSigmaPhoton.GetTrackLabelNeg() &&
             iSigmaPhoton.GetTrackLabelPos() !=
                 jSigmaPhoton.GetTrackLabelPos())) {
          iLambdaCosPA = iSigmaLambda.GetCosineAlpha();
          jLambdaCosPA = jSigmaLambda.GetCosineAlpha();
          sharedCase = 2;
        }
        // photons share daughter, but lambda does not - take Sigma candidate
        // with photon with better CPA
        else if ((iSigmaPhoton.GetTrackLabelNeg() ==
                      jSigmaPhoton.GetTrackLabelNeg() ||
                  iSigmaPhoton.GetTrackLabelNeg() ==
                      jSigmaPhoton.GetTrackLabelPos() ||
                  iSigmaPhoton.GetTrackLabelPos() ==
                      jSigmaPhoton.GetTrackLabelNeg() ||
                  iSigmaPhoton.GetTrackLabelPos() ==
                      jSigmaPhoton.GetTrackLabelPos()) &&
                 (iSigmaLambda.GetTrackLabelNeg() !=
                      jSigmaLambda.GetTrackLabelNeg() &&
                  iSigmaLambda.GetTrackLabelNeg() !=
                      jSigmaLambda.GetTrackLabelPos() &&
                  iSigmaLambda.GetTrackLabelPos() !=
                      jSigmaLambda.GetTrackLabelNeg() &&
                  iSigmaLambda.GetTrackLabelPos() !=
                      jSigmaLambda.GetTrackLabelPos())) {
          iPhotonCosPA = iSigmaPhoton.GetCosineAlpha();
          jPhotonCosPA = jSigmaPhoton.GetCosineAlpha();
          sharedCase = 3;
        }
        // both share daughter - take the sqrt of the respective CPAs
        else if ((iSigmaPhoton.GetTrackLabelNeg() ==
                      jSigmaPhoton.GetTrackLabelNeg() ||
                  iSigmaPhoton.GetTrackLabelNeg() ==
                      jSigmaPhoton.GetTrackLabelPos() ||
                  iSigmaPhoton.GetTrackLabelPos() ==
                      jSigmaPhoton.GetTrackLabelNeg() ||
                  iSigmaPhoton.GetTrackLabelPos() ==
                      jSigmaPhoton.GetTrackLabelPos()) &&
                 (iSigmaLambda.GetTrackLabelNeg() ==
                      jSigmaLambda.GetTrackLabelNeg() ||
                  iSigmaLambda.GetTrackLabelNeg() ==
                      jSigmaLambda.GetTrackLabelPos() ||
                  iSigmaLambda.GetTrackLabelPos() ==
                      jSigmaLambda.GetTrackLabelNeg() ||
                  iSigmaLambda.GetTrackLabelPos() ==
                      jSigmaLambda.GetTrackLabelPos())) {
          iLambdaCosPA = iSigmaLambda.GetCosineAlpha();
          jLambdaCosPA = jSigmaLambda.GetCosineAlpha();
          iPhotonCosPA = iSigmaPhoton.GetCosineAlpha();
          jPhotonCosPA = jSigmaPhoton.GetCosineAlpha();
          sharedCase = 4;
        }
        // mixed track sharing between any lambda & photon
        else {
          iLambdaCosPA = iSigmaLambda.GetCosineAlpha();
          jLambdaCosPA = jSigmaLambda.GetCosineAlpha();
          iPhotonCosPA = iSigmaPhoton.GetCosineAlpha();
          jPhotonCosPA = jSigmaPhoton.GetCosineAlpha();
          sharedCase = 5;
        }

        const float iSigmaSq = std::sqrt(iLambdaCosPA * iLambdaCosPA +
                                         iPhotonCosPA * iPhotonCosPA);
        const float jSigmaSq = std::sqrt(jLambdaCosPA * jLambdaCosPA +
                                         jPhotonCosPA * jPhotonCosPA);
        if (iSigmaSq >= jSigmaSq) {
          // keep i, remove j
          jSigmaCandidate->SetUse(false);
          fHistAntiSigmaSharedInvMass->Fill(sharedCase,
                                            jSigmaCandidate->GetMass());
        } else {
          // remove i, keep j
          iSigmaCandidate->SetUse(false);
          fHistAntiSigmaSharedInvMass->Fill(sharedCase,
                                            iSigmaCandidate->GetMass());
          break;
        }
      }
    }
  }
  fHistAntiSigmaNShared->Fill(nShared);
}

//____________________________________________________________________________________________________
void AliSigma0PhotonMotherCuts::ProcessMC() {
  // Loop over the MC tracks
  for (int iPart = 1; iPart < (fMCEvent->GetNumberOfTracks()); iPart++) {
    TParticle *mcParticle =
        static_cast<AliMCParticle *>(fMCEvent->GetTrack(iPart))->Particle();
    if (!mcParticle) continue;
    const int nDaughter = mcParticle->GetNDaughters();
    const int pdgCode = mcParticle->GetPdgCode();

    // Sigma0 decay
    if (std::abs(pdgCode) == 3212) {
      if (pdgCode > 0) {
        fHistMCTruthSigma0Pt->Fill(mcParticle->Pt());
        fHistMCTruthSigma0PtY->Fill(
            ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(),
                            mcParticle->GetCalcMass()),
            mcParticle->Pt());
        fHistMCTruthSigma0PtEta->Fill(mcParticle->Eta(), mcParticle->Pt());
      } else {
        fHistMCTruthAntiSigma0Pt->Fill(mcParticle->Pt());
        fHistMCTruthAntiSigma0PtY->Fill(
            ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(),
                            mcParticle->GetCalcMass()),
            mcParticle->Pt());
        fHistMCTruthAntiSigma0PtEta->Fill(mcParticle->Eta(), mcParticle->Pt());
      }

      // Sigma -> Lambda gamma
      if (nDaughter == 2) {
        TParticle *daughter1 =
            static_cast<AliMCParticle *>(
                fMCEvent->GetTrack(mcParticle->GetFirstDaughter()))
                ->Particle();
        TParticle *daughter2 =
            static_cast<AliMCParticle *>(
                fMCEvent->GetTrack(mcParticle->GetLastDaughter()))
                ->Particle();

        // decays to Lambda + gamma
        if( (daughter1->GetPdgCode() == 22 && std::abs(daughter2->GetPdgCode()) == 3122) ||
            (daughter2->GetPdgCode() == 22 && std::abs(daughter1->GetPdgCode()) == 3122)) {
          if (pdgCode > 0) {
            fHistMCTruthSigma0PhotonPt->Fill(mcParticle->Pt());
            fHistMCTruthSigma0PhotonPtY->Fill(
                ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(),
                                mcParticle->GetCalcMass()),
                mcParticle->Pt());
            fHistMCTruthSigma0PhotonPtEta->Fill(mcParticle->Eta(), mcParticle->Pt());
          } else {
            fHistMCTruthAntiSigma0PhotonPt->Fill(mcParticle->Pt());
            fHistMCTruthAntiSigma0PhotonPtY->Fill(
                ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(),
                                mcParticle->GetCalcMass()),
                mcParticle->Pt());
            fHistMCTruthAntiSigma0PhotonPtEta->Fill(mcParticle->Eta(), mcParticle->Pt());
          }
        }

        // check that the photon converts
        bool firstDaugherPhoton = IsConvertedPhoton(daughter1);
        bool secondDaugherPhoton = IsConvertedPhoton(daughter2);

        if (!firstDaugherPhoton && !secondDaugherPhoton) continue;

        TParticle *photon = daughter1;
        TParticle *lambda = daughter2;
        if (secondDaugherPhoton) {
          photon = daughter2;
          lambda = daughter1;
        }

        if (std::abs(lambda->GetPdgCode()) != 3122) continue;

        if (pdgCode > 0) {
          fHistMCTruthSigma0PhotonConvP->Fill(photon->P());
          fHistMCTruthSigma0PhotonConvPt->Fill(photon->Pt());
          fHistMCTruthSigma0PhotonConvInvMass->Fill(photon->GetCalcMass());
          fHistMCTruthSigma0PhotonConvInvMassPt->Fill(photon->Pt(),
                                                      photon->GetCalcMass());
          fHistMCTruthSigma0PhotonConvPtEta->Fill(photon->Eta(), photon->Pt());
          fHistMCTruthSigma0PhotonConvPtY->Fill(
              ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(),
                              mcParticle->GetCalcMass()),
              mcParticle->Pt());
        } else {
          fHistMCTruthAntiSigma0PhotonConvP->Fill(photon->P());
          fHistMCTruthAntiSigma0PhotonConvPt->Fill(photon->Pt());
          fHistMCTruthAntiSigma0PhotonConvInvMass->Fill(photon->GetCalcMass());
          fHistMCTruthAntiSigma0PhotonConvInvMassPt->Fill(
              photon->Pt(), photon->GetCalcMass());
          fHistMCTruthAntiSigma0PhotonConvPtEta->Fill(photon->Eta(),
                                                      photon->Pt());
          fHistMCTruthAntiSigma0PhotonConvPtY->Fill(
              ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(),
                              mcParticle->GetCalcMass()),
              mcParticle->Pt());
        }

        TParticle *ePos = nullptr;
        TParticle *eNeg = nullptr;
        AliMCParticle *partEle = nullptr;
        if (photon->GetNDaughters() >= 2) {
          for (int daughterIndex = photon->GetFirstDaughter();
               daughterIndex <= photon->GetLastDaughter(); ++daughterIndex) {
            if (daughterIndex < 0) continue;
            TParticle *tmpDaughter =
                static_cast<AliMCParticle *>(fMCEvent->GetTrack(daughterIndex))
                    ->Particle();
            if (tmpDaughter->GetUniqueID() == 5) {
              if (tmpDaughter->GetPdgCode() == 11) {
                eNeg = tmpDaughter;
                partEle = static_cast<AliMCParticle *>(
                    fMCEvent->GetTrack(daughterIndex));
              } else if (tmpDaughter->GetPdgCode() == -11)
                ePos = tmpDaughter;
            }
          }
        }

        if (ePos != nullptr) {
          if (pdgCode > 0) {
            fHistMCTruthSigma0PhotonConvEleP->Fill(ePos->P());
            fHistMCTruthSigma0PhotonConvElePt->Fill(ePos->Pt());
          } else {
            fHistMCTruthAntiSigma0PhotonConvEleP->Fill(ePos->P());
            fHistMCTruthAntiSigma0PhotonConvElePt->Fill(ePos->Pt());
          }
        }
        if (eNeg != nullptr) {
          if (pdgCode > 0) {
            fHistMCTruthSigma0PhotonConvEleP->Fill(ePos->P());
            fHistMCTruthSigma0PhotonConvElePt->Fill(ePos->Pt());
          } else {
            fHistMCTruthAntiSigma0PhotonConvEleP->Fill(ePos->P());
            fHistMCTruthAntiSigma0PhotonConvElePt->Fill(ePos->Pt());
          }
        }

        if (eNeg != nullptr && ePos != nullptr) {
          float conversionVertexX = partEle->Xv();
          float conversionVertexY = partEle->Yv();
          float conversionVertexZ = partEle->Zv();
          if (pdgCode > 0) {
            fHistMCTruthSigma0PhotonConvR->Fill(
                photon->Pt(), std::sqrt(conversionVertexX * conversionVertexX +
                                        conversionVertexY * conversionVertexY));
            fHistMCTruthSigma0PhotonConvConvPointX->Fill(conversionVertexX);
            fHistMCTruthSigma0PhotonConvConvPointY->Fill(conversionVertexY);
            fHistMCTruthSigma0PhotonConvConvPointZ->Fill(conversionVertexZ);
          } else {
            fHistMCTruthAntiSigma0PhotonConvR->Fill(
                photon->Pt(), std::sqrt(conversionVertexX * conversionVertexX +
                                        conversionVertexY * conversionVertexY));
            fHistMCTruthAntiSigma0PhotonConvConvPointX->Fill(conversionVertexX);
            fHistMCTruthAntiSigma0PhotonConvConvPointY->Fill(conversionVertexY);
            fHistMCTruthAntiSigma0PhotonConvConvPointZ->Fill(conversionVertexZ);
          }
        }
      }

      if (nDaughter == 3) {
        // Sigma -> Lambda e+ e-
        TParticle *ePos = nullptr;
        TParticle *eNeg = nullptr;
        TParticle *lambda = nullptr;
        for (int daughterIndex = mcParticle->GetFirstDaughter();
             daughterIndex <= mcParticle->GetLastDaughter(); ++daughterIndex) {
          if (daughterIndex < 0) continue;
          TParticle *tmpDaughter =
              static_cast<AliMCParticle *>(fMCEvent->GetTrack(daughterIndex))
                  ->Particle();
          std::cout << daughterIndex << " " << tmpDaughter->GetPdgCode()
                    << "\n";
          const int pdgCode = tmpDaughter->GetPdgCode();
          if (pdgCode == 11)
            eNeg = tmpDaughter;
          else if (pdgCode == -11)
            ePos = tmpDaughter;
          else if (std::abs(pdgCode) == 3122)
            lambda = tmpDaughter;
        }

        if (ePos != nullptr && eNeg != nullptr && lambda != nullptr) {
          if (pdgCode > 0) {
            fHistMCTruthSigma0TwoElectronPtPos->Fill(ePos->Pt());
            fHistMCTruthSigma0TwoElectronPPos->Fill(ePos->P());
            fHistMCTruthSigma0TwoElectronPtEtaPos->Fill(ePos->Eta(),
                                                        ePos->Pt());
            fHistMCTruthSigma0TwoElectronPtNeg->Fill(eNeg->Pt());
            fHistMCTruthSigma0TwoElectronPNeg->Fill(eNeg->P());
            fHistMCTruthSigma0TwoElectronPtEtaNeg->Fill(eNeg->Eta(),
                                                        eNeg->Pt());
          } else if (pdgCode < 0) {
            fHistMCTruthAntiSigma0TwoElectronPtPos->Fill(ePos->Pt());
            fHistMCTruthAntiSigma0TwoElectronPPos->Fill(ePos->P());
            fHistMCTruthAntiSigma0TwoElectronPtEtaPos->Fill(ePos->Eta(),
                                                            ePos->Pt());
            fHistMCTruthAntiSigma0TwoElectronPtNeg->Fill(eNeg->Pt());
            fHistMCTruthAntiSigma0TwoElectronPNeg->Fill(eNeg->P());
            fHistMCTruthAntiSigma0TwoElectronPtEtaNeg->Fill(eNeg->Eta(),
                                                            eNeg->Pt());
          }

          // Sigma -> Lambda Gamma Gamma
          TParticle *photon1 = nullptr;
          TParticle *photon2 = nullptr;
          TParticle *lambdaPhoton = nullptr;
          for (int daughterIndex = mcParticle->GetFirstDaughter();
               daughterIndex <= mcParticle->GetLastDaughter();
               ++daughterIndex) {
            if (daughterIndex < 0) continue;
            TParticle *tmpDaughter =
                static_cast<AliMCParticle *>(fMCEvent->GetTrack(daughterIndex))
                    ->Particle();
            if (IsConvertedPhoton(tmpDaughter)) {
              if (photon1 == nullptr)
                photon1 = tmpDaughter;
              else if (photon1 != nullptr && photon2 == nullptr)
                photon2 = tmpDaughter;
            } else if (std::abs(tmpDaughter->GetPdgCode()) == 3122)
              lambdaPhoton = tmpDaughter;
          }

          if (photon1 != nullptr && photon2 != nullptr &&
              lambdaPhoton != nullptr) {
            TParticle *ePosPhoton1 = nullptr;
            TParticle *eNegPhoton1 = nullptr;
            if (photon1->GetNDaughters() >= 2) {
              for (int daughterIndex = photon1->GetFirstDaughter();
                   daughterIndex <= photon1->GetLastDaughter();
                   ++daughterIndex) {
                if (daughterIndex < 0) continue;
                TParticle *tmpDaughter = static_cast<AliMCParticle *>(
                                             fMCEvent->GetTrack(daughterIndex))
                                             ->Particle();
                if (tmpDaughter->GetUniqueID() == 5) {
                  if (tmpDaughter->GetPdgCode() == 11)
                    eNegPhoton1 = tmpDaughter;
                  else if (tmpDaughter->GetPdgCode() == -11)
                    ePosPhoton1 = tmpDaughter;
                }
              }
            }
            TParticle *ePosPhoton2 = nullptr;
            TParticle *eNegPhoton2 = nullptr;
            if (photon2->GetNDaughters() >= 2) {
              for (int daughterIndex = photon2->GetFirstDaughter();
                   daughterIndex <= photon2->GetLastDaughter();
                   ++daughterIndex) {
                if (daughterIndex < 0) continue;
                TParticle *tmpDaughter = static_cast<AliMCParticle *>(
                                             fMCEvent->GetTrack(daughterIndex))
                                             ->Particle();
                if (tmpDaughter->GetUniqueID() == 5) {
                  if (tmpDaughter->GetPdgCode() == 11)
                    eNegPhoton2 = tmpDaughter;
                  else if (tmpDaughter->GetPdgCode() == -11)
                    ePosPhoton2 = tmpDaughter;
                }
              }
            }

            if (pdgCode > 0) {
              fHistMCTruthSigma0TwoGammaPt->Fill(photon1->Pt());
              fHistMCTruthSigma0TwoGammaP->Fill(photon1->P());
              fHistMCTruthSigma0TwoGammaPtEta->Fill(photon1->Eta(),
                                                    photon1->Pt());
              fHistMCTruthSigma0TwoGammaPt->Fill(photon2->Pt());
              fHistMCTruthSigma0TwoGammaP->Fill(photon2->P());
              fHistMCTruthSigma0TwoGammaPtEta->Fill(photon2->Eta(),
                                                    photon2->Pt());

              fHistMCTruthSigma0TwoGammaPtEle->Fill(ePosPhoton1->Pt());
              fHistMCTruthSigma0TwoGammaPEle->Fill(ePosPhoton1->P());
              fHistMCTruthSigma0TwoGammaPtEtaEle->Fill(ePosPhoton1->Eta(),
                                                       ePosPhoton1->Pt());
              fHistMCTruthSigma0TwoGammaPtEle->Fill(eNegPhoton1->Pt());
              fHistMCTruthSigma0TwoGammaPEle->Fill(eNegPhoton1->P());
              fHistMCTruthSigma0TwoGammaPtEtaEle->Fill(eNegPhoton1->Eta(),
                                                       eNegPhoton1->Pt());
              fHistMCTruthSigma0TwoGammaPtEle->Fill(ePosPhoton2->Pt());
              fHistMCTruthSigma0TwoGammaPEle->Fill(ePosPhoton2->P());
              fHistMCTruthSigma0TwoGammaPtEtaEle->Fill(ePosPhoton2->Pt(),
                                                       ePosPhoton2->Eta());
              fHistMCTruthSigma0TwoGammaPtEle->Fill(eNegPhoton2->Pt());
              fHistMCTruthSigma0TwoGammaPEle->Fill(eNegPhoton2->P());
              fHistMCTruthSigma0TwoGammaPtEtaEle->Fill(eNegPhoton2->Eta(),
                                                       eNegPhoton2->Pt());
            } else if (pdgCode < 0) {
              fHistMCTruthAntiSigma0TwoGammaPt->Fill(photon1->Pt());
              fHistMCTruthAntiSigma0TwoGammaP->Fill(photon1->P());
              fHistMCTruthAntiSigma0TwoGammaPtEta->Fill(photon1->Eta(),
                                                        photon1->Pt());
              fHistMCTruthAntiSigma0TwoGammaPt->Fill(photon2->Pt());
              fHistMCTruthAntiSigma0TwoGammaP->Fill(photon2->P());
              fHistMCTruthAntiSigma0TwoGammaPtEta->Fill(photon2->Eta(),
                                                        photon2->Pt());

              fHistMCTruthAntiSigma0TwoGammaPtEle->Fill(ePosPhoton1->Pt());
              fHistMCTruthAntiSigma0TwoGammaPEle->Fill(ePosPhoton1->P());
              fHistMCTruthAntiSigma0TwoGammaPtEtaEle->Fill(ePosPhoton1->Eta(),
                                                           ePosPhoton1->Pt());
              fHistMCTruthAntiSigma0TwoGammaPtEle->Fill(eNegPhoton1->Pt());
              fHistMCTruthAntiSigma0TwoGammaPEle->Fill(eNegPhoton1->P());
              fHistMCTruthAntiSigma0TwoGammaPtEtaEle->Fill(eNegPhoton1->Eta(),
                                                           eNegPhoton1->Pt());
              fHistMCTruthAntiSigma0TwoGammaPtEle->Fill(ePosPhoton2->Pt());
              fHistMCTruthAntiSigma0TwoGammaPEle->Fill(ePosPhoton2->P());
              fHistMCTruthAntiSigma0TwoGammaPtEtaEle->Fill(ePosPhoton2->Eta(),
                                                           ePosPhoton2->Pt());
              fHistMCTruthAntiSigma0TwoGammaPtEle->Fill(eNegPhoton2->Pt());
              fHistMCTruthAntiSigma0TwoGammaPEle->Fill(eNegPhoton2->P());
              fHistMCTruthAntiSigma0TwoGammaPtEtaEle->Fill(eNegPhoton2->Eta(),
                                                           eNegPhoton2->Pt());
            }
          }
        }
      }
    }
  }
}

//____________________________________________________________________________________________________
bool AliSigma0PhotonMotherCuts::IsConvertedPhoton(TParticle *particle) const {
  if (particle->GetPdgCode() != 22) return false;

  // check if particle doesn't have a photon as mother
  TParticle *mcMother =
      static_cast<AliMCParticle *>(fMCEvent->GetTrack(particle->GetMother(0)))
          ->Particle();
  if (mcMother->GetPdgCode() == 22) return false;

  // looking for conversion gammas (electron + positron from pairbuilding (= 5)
  // )
  TParticle *ePos = nullptr;
  TParticle *eNeg = nullptr;
  if (particle->GetNDaughters() >= 2) {
    for (int daughterIndex = particle->GetFirstDaughter();
         daughterIndex <= particle->GetLastDaughter(); ++daughterIndex) {
      if (daughterIndex < 0) continue;
      TParticle *tmpDaughter =
          static_cast<AliMCParticle *>(fMCEvent->GetTrack(daughterIndex))
              ->Particle();
      if (tmpDaughter->GetUniqueID() == 5) {
        if (tmpDaughter->GetPdgCode() == 11)
          eNeg = tmpDaughter;
        else if (tmpDaughter->GetPdgCode() == -11)
          ePos = tmpDaughter;
      }
    }
  }
  // Two daughters from pair production?
  if (ePos == nullptr || eNeg == nullptr) return false;
  return true;
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
void AliSigma0PhotonMotherCuts::InitCutHistograms() {
  std::cout << "============================\n"
            << " PHOTON MOTHER CUT CONFIGURATION \n"
            << " Sigma0 mass     " << fMassSigma << "\n"
            << " Sigma0 select   " << fSigmaMassCut << "\n"
            << " Sigma0 sideb    " << fSigmaSidebandLow << " "
            << fSigmaSidebandUp << "\n"
            << " Mixing depth    " << fMixingDepth << "\n"
            << "============================\n";

  TH1::AddDirectory(kFALSE);

  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    fHistograms->SetName("PhotonMotherCut_QA");
  }

  if (fHistogramsSigma != nullptr) {
    delete fHistogramsSigma;
    fHistogramsSigma = nullptr;
  }
  if (fHistogramsSigma == nullptr) {
    fHistogramsSigma = new TList();
    fHistogramsSigma->SetOwner(kTRUE);
    fHistogramsSigma->SetName("SigmaCut_QA");
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

  fHistNSigma =
      new TH1F("fHistNSigma", ";# #Sigma candidates; Entries", 20, 0, 20);
  fHistSigmaPt =
      new TH1F("fHistSigmaPt",
               "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries", 500, 0, 20);
  fHistSigmaMassCutPt =
      new TH1F("fHistSigmaMassCutPt",
               "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries", 500, 0, 20);
  fHistSigmaInvMass =
      new TH1F("fHistSigmaInvMass",
               "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistSigmaInvMassFinal = new TH1F(
      "fHistSigmaInvMassFinal",
      "after track cleaning; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 500,
      1., 2.);
  fHistSigmaInvMassRec =
      new TH1F("fHistSigmaInvMassRec",
               "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistSigmaInvMassPt = new TH2F("fHistSigmaInvMassPt",
                                 "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
                                 "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                                 1000, 0, 20, 2000, 1., 2.);
  fHistSigmaInvMassEta = new TH2F("fHistSigmaInvMassEta",
                                  "; #eta; M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                                  1000, 0, 1, 2000, 1., 2.);
  fHistSigmaNShared = new TH1F(
      "fHistSigmaNShared",
      "; # shared V0#gamma-V0#gamma pairs per event; Entries", 50, 0, 50);
  fHistSigmaMassDiff =
      new TH1F("fHistSigmaMassDiff", "; #DeltaM #Sigma^{0} candidates; Entries",
               1000, -0.1, 0.1);
  fHistSigmaSharedInvMass = new TH2F("fHistSigmaSharedInvMass",
                                     ";; M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                                     10, 0, 10, 1000, 1., 1.5);
  fHistSigmaSharedInvMass->GetXaxis()->SetBinLabel(1, "#Sigma^{0}");
  fHistSigmaSharedInvMass->GetXaxis()->SetBinLabel(2, "direct");
  fHistSigmaSharedInvMass->GetXaxis()->SetBinLabel(3, "#Lambda_{1}#Lambda_{2}");
  fHistSigmaSharedInvMass->GetXaxis()->SetBinLabel(4, "#gamma_{1}#gamma_{2}");
  fHistSigmaSharedInvMass->GetXaxis()->SetBinLabel(
      5, "#Lambda_{1}#Lambda_{2} && #gamma_{1}#gamma_{2}");
  fHistSigmaSharedInvMass->GetXaxis()->SetBinLabel(
      6, "#Lambda_{1}#gamma_{2} || #Lambda_{2}#gamma_{1}");
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
  fHistSigmaMixedPt =
      new TH1F("fHistSigmaMixedPt",
               "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); Entries", 500, 0, 20);
  fHistSigmaMixedInvMass =
      new TH1F("fHistSigmaMixedInvMass",
               "; M_{#Lambda#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistSigmaMixedInvMassPt =
      new TH2F("fHistSigmaMixedInvMassPt",
               "; #it{p}_{T} #Lambda#gamma (GeV/#it{c}); "
               "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistSigmaMixedInvMassEta = new TH2F(
      "fHistSigmaMixedInvMassEta", "; #eta; M_{#Lambda#gamma} (GeV/#it{c}^{2})",
      1000, 0, 1, 2000, 1., 2.);

  fHistSigmaPtTwoGamma = new TH1F(
      "fHistSigmaPtTwoGamma",
      "; #it{p}_{T} #Lambda#gamma#gamma (GeV/#it{c}); Entries", 500, 0, 20);
  fHistSigmaInvMassTwoGamma = new TH1F(
      "fHistSigmaInvMassTwoGamma",
      "; M_{#Lambda#gamma#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistSigmaInvMassPtTwoGamma =
      new TH2F("fHistSigmaInvMassPtTwoGamma",
               "; #it{p}_{T} #Lambda#gamma#gamma (GeV/#it{c}); "
               "M_{#Lambda#gamma#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistSigmaInvMassEtaTwoGamma =
      new TH2F("fHistSigmaInvMassEtaTwoGamma",
               "; #eta; M_{#Lambda#gamma#gamma} (GeV/#it{c}^{2})", 1000, 0, 1,
               2000, 1., 2.);
  fHistSigmaArmenterosBeforeTwoGamma =
      new TH2F("fHistSigmaArmenterosBeforeTwoGamma",
               " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
  fHistSigmaArmenterosAfterTwoGamma =
      new TH2F("fHistSigmaArmenterosAfterTwoGamma",
               " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
  fHistSigmaMixedPtTwoGamma = new TH1F(
      "fHistSigmaMixedPtTwoGamma",
      "; #it{p}_{T} #Lambda#gamma#gamma (GeV/#it{c}); Entries", 500, 0, 20);
  fHistSigmaMixedInvMassTwoGamma = new TH1F(
      "fHistSigmaMixedInvMassTwoGamma",
      "; M_{#Lambda#gamma#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistSigmaMixedInvMassPtTwoGamma =
      new TH2F("fHistSigmaMixedInvMassPtTwoGamma",
               "; #it{p}_{T} #Lambda#gamma#gamma (GeV/#it{c}); "
               "M_{#Lambda#gamma#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistSigmaMixedInvMassEtaTwoGamma =
      new TH2F("fHistSigmaMixedInvMassEtaTwoGamma",
               "; #eta; M_{#Lambda#gamma#gamma} (GeV/#it{c}^{2})", 1000, 0, 1,
               2000, 1., 2.);

  fHistSigmaPtDiElectron = new TH1F(
      "fHistSigmaPtDiElectron",
      "; #it{p}_{T} #Lambdae^{+}e^{-} (GeV/#it{c}); Entries", 500, 0, 20);
  fHistSigmaInvMassDiElectron = new TH1F(
      "fHistSigmaInvMassDiElectron",
      "; M_{#Lambdae^{+}e^{-}} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistSigmaInvMassPtDiElectron =
      new TH2F("fHistSigmaInvMassPtDiElectron",
               "; #it{p}_{T} #Lambdae^{+}e^{-} (GeV/#it{c}); "
               "M_{#Lambdae^{+}e^{-}} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistSigmaInvMassEtaDiElectron =
      new TH2F("fHistSigmaInvMassEtaDiElectron",
               "; #eta; M_{#Lambdae^{+}e^{-}} (GeV/#it{c}^{2})", 1000, 0, 1,
               2000, 1., 2.);
  fHistSigmaArmenterosBeforeDiElectron =
      new TH2F("fHistSigmaArmenterosBeforeDiElectron",
               " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
  fHistSigmaArmenterosAfterDiElectron =
      new TH2F("fHistSigmaArmenterosAfterDiElectron",
               " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
  fHistSigmaMixedPtDiElectron = new TH1F(
      "fHistSigmaMixedPtDiElectron",
      "; #it{p}_{T} #Lambdae^{+}e^{-} (GeV/#it{c}); Entries", 500, 0, 20);
  fHistSigmaMixedInvMassDiElectron = new TH1F(
      "fHistSigmaMixedInvMassDiElectron",
      "; M_{#Lambdae^{+}e^{-}} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistSigmaMixedInvMassPtDiElectron =
      new TH2F("fHistSigmaMixedInvMassPtDiElectron",
               "; #it{p}_{T} #Lambdae^{+}e^{-} (GeV/#it{c}); "
               "M_{#Lambdae^{+}e^{-}} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistSigmaMixedInvMassEtaDiElectron =
      new TH2F("fHistSigmaMixedInvMassEtaDiElectron",
               "; #eta; M_{#Lambdae^{+}e^{-}} (GeV/#it{c}^{2})", 1000, 0, 1,
               2000, 1., 2.);

  fHistogramsSigma->Add(fHistNSigma);
  fHistogramsSigma->Add(fHistSigmaPt);
  fHistogramsSigma->Add(fHistSigmaMassCutPt);
  fHistogramsSigma->Add(fHistSigmaInvMass);
  fHistogramsSigma->Add(fHistSigmaInvMassFinal);
  fHistogramsSigma->Add(fHistSigmaInvMassRec);
  fHistogramsSigma->Add(fHistSigmaInvMassPt);
  fHistogramsSigma->Add(fHistSigmaInvMassEta);
  fHistogramsSigma->Add(fHistSigmaEtaPhi);
  fHistogramsSigma->Add(fHistSigmaRapidity);
  fHistogramsSigma->Add(fHistSigmaNShared);
  fHistogramsSigma->Add(fHistSigmaMassDiff);
  fHistogramsSigma->Add(fHistSigmaSharedInvMass);
  fHistogramsSigma->Add(fHistSigmaArmenterosBefore);
  fHistogramsSigma->Add(fHistSigmaArmenterosAfter);
  fHistogramsSigma->Add(fHistSigmaMixedPt);
  fHistogramsSigma->Add(fHistSigmaMixedInvMass);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassPt);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassEta);
  fHistogramsSigma->Add(fHistSigmaPtTwoGamma);
  fHistogramsSigma->Add(fHistSigmaInvMassTwoGamma);
  fHistogramsSigma->Add(fHistSigmaInvMassPtTwoGamma);
  fHistogramsSigma->Add(fHistSigmaInvMassEtaTwoGamma);
  fHistogramsSigma->Add(fHistSigmaArmenterosBeforeTwoGamma);
  fHistogramsSigma->Add(fHistSigmaArmenterosAfterTwoGamma);
  fHistogramsSigma->Add(fHistSigmaMixedPtTwoGamma);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassTwoGamma);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassPtTwoGamma);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassEtaTwoGamma);
  fHistogramsSigma->Add(fHistSigmaPtDiElectron);
  fHistogramsSigma->Add(fHistSigmaInvMassDiElectron);
  fHistogramsSigma->Add(fHistSigmaInvMassPtDiElectron);
  fHistogramsSigma->Add(fHistSigmaInvMassEtaDiElectron);
  fHistogramsSigma->Add(fHistSigmaArmenterosBeforeDiElectron);
  fHistogramsSigma->Add(fHistSigmaArmenterosAfterDiElectron);
  fHistogramsSigma->Add(fHistSigmaMixedPtDiElectron);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassDiElectron);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassPtDiElectron);
  fHistogramsSigma->Add(fHistSigmaMixedInvMassEtaDiElectron);

  std::vector<float> rapBins = {
      {-10.f, -1.f, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.f,
       0.1,   0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.f,  10.f}};
  for (int i = 0; i < static_cast<int>(rapBins.size() - 1); ++i) {
    fHistSigmaPtY[i] =
        new TH2F(Form("fHistSigmaPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
                 Form("%.2f < y < %.2f ; #it{p}_{T} [GeV/#it{c}]; "
                      "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                      rapBins[i], rapBins[i + 1]),
                 500, 0, 10, 1000, 1., 1.5);
    fHistSigmaMixedPtY[i] =
        new TH2F(Form("fHistSigmaMixedPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
                 Form("%.2f < y < %.2f ; #it{p}_{T} [GeV/#it{c}]; "
                      "M_{#Lambda#gamma} (GeV/#it{c}^{2})",
                      rapBins[i], rapBins[i + 1]),
                 500, 0, 10, 1000, 1., 1.5);
    fHistogramsSigma->Add(fHistSigmaPtY[i]);
    fHistogramsSigma->Add(fHistSigmaMixedPtY[i]);
  }
  fHistograms->Add(fHistogramsSigma);

  if (fHistogramsAntiSigma != nullptr) {
    delete fHistogramsAntiSigma;
    fHistogramsAntiSigma = nullptr;
  }
  if (fHistogramsAntiSigma == nullptr) {
    fHistogramsAntiSigma = new TList();
    fHistogramsAntiSigma->SetOwner(kTRUE);
    fHistogramsAntiSigma->SetName("AntiSigmaCut_QA");
  }

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
  fHistAntiSigmaInvMassFinal = new TH1F(
      "fHistAntiSigmaInvMassFinal",
      "after track cleaning; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2}); Entries",
      2000, 1., 2.);
  fHistAntiSigmaInvMassRec = new TH1F(
      "fHistAntiSigmaInvMassRec",
      "; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistAntiSigmaInvMassPt =
      new TH2F("fHistAntiSigmaInvMassPt",
               "; #it{p}_{T}  #bar{#Lambda}#gamma (GeV/#it{c}); "
               "M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistAntiSigmaInvMassEta =
      new TH2F("fHistAntiSigmaInvMassEta",
               "; #eta; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})", 1000, 0, 1,
               2000, 1., 2.);
  fHistAntiSigmaNShared = new TH1F(
      "fHistAntiSigmaNShared",
      "; # shared #bar{V0}#gamma-#bar{V0}#gamma pairs per event; Entries", 50,
      0, 50);
  fHistAntiSigmaMassDiff = new TH1F(
      "fHistAntiSigmaMassDiff",
      "; #DeltaM #bar{#Sigma^{0}} candidates; Entries", 1000, -0.1, 0.1);
  fHistAntiSigmaSharedInvMass = new TH2F(
      "fHistAntiSigmaSharedInvMass",
      ";; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})", 10, 0, 10, 1000, 1., 1.5);
  fHistAntiSigmaSharedInvMass->GetXaxis()->SetBinLabel(1, "#Sigma^{0}");
  fHistAntiSigmaSharedInvMass->GetXaxis()->SetBinLabel(2, "direct");
  fHistAntiSigmaSharedInvMass->GetXaxis()->SetBinLabel(
      3, "#bar{#Lambda}_{1}#bar{#Lambda}_{2}");
  fHistAntiSigmaSharedInvMass->GetXaxis()->SetBinLabel(4,
                                                       "#gamma_{1}#gamma_{2}");
  fHistAntiSigmaSharedInvMass->GetXaxis()->SetBinLabel(
      5, "#bar{#Lambda}_{1}#bar{#Lambda}_{2} && #gamma_{1}#gamma_{2}");
  fHistAntiSigmaSharedInvMass->GetXaxis()->SetBinLabel(
      6, "#bar{#Lambda}_{1}#gamma_{2} || #bar{#Lambda}_{2}#gamma_{1}");
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
  fHistAntiSigmaMixedInvMassPt =
      new TH2F("fHistAntiSigmaMixedInvMassPt",
               "; #it{p}_{T} #bar{#Lambda}#gamma (GeV/#it{c}); "
               "M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistAntiSigmaMixedInvMassEta =
      new TH2F("fHistAntiSigmaMixedInvMassEta",
               "; #eta; M_{#bar{#Lambda}#gamma} (GeV/#it{c}^{2})", 1000, 0, 1,
               2000, 1., 2.);

  fHistAntiSigmaPtTwoGamma =
      new TH1F("fHistAntiSigmaPtTwoGamma",
               "; #it{p}_{T} #bar{#Lambda}#gamma#gamma (GeV/#it{c}); Entries",
               500, 0, 20);
  fHistAntiSigmaInvMassTwoGamma =
      new TH1F("fHistAntiSigmaInvMassTwoGamma",
               "; M_{#bar{#Lambda}#gamma#gamma} (GeV/#it{c}^{2}); Entries",
               2000, 1., 2.);
  fHistAntiSigmaInvMassPtTwoGamma =
      new TH2F("fHistAntiSigmaInvMassPtTwoGamma",
               "; #it{p}_{T} #bar{#Lambda}#gamma#gamma (GeV/#it{c}); "
               "M_{#bar{#Lambda}#gamma#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistAntiSigmaInvMassEtaTwoGamma =
      new TH2F("fHistAntiSigmaInvMassEtaTwoGamma",
               "; #eta; M_{#bar{#Lambda}#gamma#gamma} (GeV/#it{c}^{2})", 1000,
               0, 1, 2000, 1., 2.);
  fHistAntiSigmaArmenterosBeforeTwoGamma =
      new TH2F("fHistAntiSigmaArmenterosBeforeTwoGamma",
               " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
  fHistAntiSigmaArmenterosAfterTwoGamma =
      new TH2F("fHistAntiSigmaArmenterosAfterTwoGamma",
               " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
  fHistAntiSigmaMixedPtTwoGamma =
      new TH1F("fHistAntiSigmaMixedPtTwoGamma",
               "; #it{p}_{T} #bar{#Lambda}#gamma#gamma (GeV/#it{c}); Entries",
               500, 0, 20);
  fHistAntiSigmaMixedInvMassTwoGamma =
      new TH1F("fHistAntiSigmaMixedInvMassTwoGamma",
               "; M_{#bar{#Lambda}#gamma#gamma} (GeV/#it{c}^{2}); Entries",
               2000, 1., 2.);
  fHistAntiSigmaMixedInvMassPtTwoGamma =
      new TH2F("fHistAntiSigmaMixedInvMassPtTwoGamma",
               "; #it{p}_{T} #bar{#Lambda}#gamma#gamma (GeV/#it{c}); "
               "M_{#bar{#Lambda}#gamma#gamma} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistAntiSigmaMixedInvMassEtaTwoGamma =
      new TH2F("fHistAntiSigmaMixedInvMassEtaTwoGamma",
               "; #eta; M_{#bar{#Lambda}#gamma#gamma} (GeV/#it{c}^{2})", 1000,
               0, 1, 2000, 1., 2.);

  fHistAntiSigmaPtDiElectron = new TH1F(
      "fHistAntiSigmaPtDiElectron",
      "; #it{p}_{T} #bar{#Lambda}e^{+}e^{-} (GeV/#it{c}); Entries", 500, 0, 20);
  fHistAntiSigmaInvMassDiElectron = new TH1F(
      "fHistAntiSigmaInvMassDiElectron",
      "; M_{#bar{#Lambda}e^{+}e^{-}} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistAntiSigmaInvMassPtDiElectron =
      new TH2F("fHistAntiSigmaInvMassPtDiElectron",
               "; #it{p}_{T} #bar{#Lambda}e^{+}e^{-} (GeV/#it{c}); "
               "M_{#bar{#Lambda}e^{+}e^{-}} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistAntiSigmaInvMassEtaDiElectron =
      new TH2F("fHistAntiSigmaInvMassEtaDiElectron",
               "; #eta; M_{#bar{#Lambda}e^{+}e^{-}} (GeV/#it{c}^{2})", 1000, 0,
               1, 2000, 1., 2.);
  fHistAntiSigmaArmenterosBeforeDiElectron =
      new TH2F("fHistAntiSigmaArmenterosBeforeDiElectron",
               " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
  fHistAntiSigmaArmenterosAfterDiElectron =
      new TH2F("fHistAntiSigmaArmenterosAfterDiElectron",
               " ; #alpha; #it{q}_{T} [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
  fHistAntiSigmaMixedPtDiElectron = new TH1F(
      "fHistAntiSigmaMixedPtDiElectron",
      "; #it{p}_{T} #bar{#Lambda}e^{+}e^{-} (GeV/#it{c}); Entries", 500, 0, 20);
  fHistAntiSigmaMixedInvMassDiElectron = new TH1F(
      "fHistAntiSigmaMixedInvMassDiElectron",
      "; M_{#bar{#Lambda}e^{+}e^{-}} (GeV/#it{c}^{2}); Entries", 2000, 1., 2.);
  fHistAntiSigmaMixedInvMassPtDiElectron =
      new TH2F("fHistAntiSigmaMixedInvMassPtDiElectron",
               "; #it{p}_{T} #bar{#Lambda}e^{+}e^{-} (GeV/#it{c}); "
               "M_{#bar{#Lambda}e^{+}e^{-}} (GeV/#it{c}^{2})",
               1000, 0, 20, 2000, 1., 2.);
  fHistAntiSigmaMixedInvMassEtaDiElectron =
      new TH2F("fHistAntiSigmaMixedInvMassEtaDiElectron",
               "; #eta; M_{#bar{#Lambda}e^{+}e^{-}} (GeV/#it{c}^{2})", 1000, 0,
               1, 2000, 1., 2.);

  fHistogramsAntiSigma->Add(fHistNAntiSigma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaPt);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMassCutPt);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMass);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassFinal);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassRec);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassPt);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassEta);
  fHistogramsAntiSigma->Add(fHistAntiSigmaRapidity);
  fHistogramsAntiSigma->Add(fHistAntiSigmaEtaPhi);
  fHistogramsAntiSigma->Add(fHistAntiSigmaNShared);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMassDiff);
  fHistogramsAntiSigma->Add(fHistAntiSigmaSharedInvMass);
  fHistogramsAntiSigma->Add(fHistAntiSigmaArmenterosBefore);
  fHistogramsAntiSigma->Add(fHistAntiSigmaArmenterosAfter);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedPt);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMass);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassPt);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassEta);
  fHistogramsAntiSigma->Add(fHistAntiSigmaPtTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassPtTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassEtaTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaArmenterosBeforeTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaArmenterosAfterTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedPtTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassPtTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassEtaTwoGamma);
  fHistogramsAntiSigma->Add(fHistAntiSigmaPtDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassPtDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaInvMassEtaDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaArmenterosBeforeDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaArmenterosAfterDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedPtDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassPtDiElectron);
  fHistogramsAntiSigma->Add(fHistAntiSigmaMixedInvMassEtaDiElectron);

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
  fHistograms->Add(fHistogramsAntiSigma);

  if (fIsMC) {
    if (fHistogramsSigmaMC != nullptr) {
      delete fHistogramsSigmaMC;
      fHistogramsSigmaMC = nullptr;
    }
    if (fHistogramsSigmaMC == nullptr) {
      fHistogramsSigmaMC = new TList();
      fHistogramsSigmaMC->SetOwner(kTRUE);
      fHistogramsSigmaMC->SetName("SigmaCut_MC");
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

    fHistMCTruthSigma0TwoElectronPtPos =
        new TH1F("fHistMCTruthSigma0TwoElectronPtPos",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0TwoElectronPPos =
        new TH1F("fHistMCTruthSigma0TwoElectronPPos",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0TwoElectronPtEtaPos =
        new TH2F("fHistMCTruthSigma0TwoElectronPtEtaPos",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistMCTruthSigma0TwoElectronPtNeg =
        new TH1F("fHistMCTruthSigma0TwoElectronPtNeg",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0TwoElectronPNeg =
        new TH1F("fHistMCTruthSigma0TwoElectronPNeg",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0TwoElectronPtEtaNeg =
        new TH2F("fHistMCTruthSigma0TwoElectronPtEtaNeg",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoElectronPtPos);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoElectronPPos);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoElectronPtEtaPos);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoElectronPtNeg);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoElectronPNeg);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoElectronPtEtaNeg);

    fHistMCTruthSigma0TwoGammaPt =
        new TH1F("fHistMCTruthSigma0TwoGammaPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0TwoGammaP =
        new TH1F("fHistMCTruthSigma0TwoGammaP",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0TwoGammaPtEta =
        new TH2F("fHistMCTruthSigma0TwoGammaPtEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistMCTruthSigma0TwoGammaPtEle =
        new TH1F("fHistMCTruthSigma0TwoGammaPtEle",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0TwoGammaPEle =
        new TH1F("fHistMCTruthSigma0TwoGammaPEle",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthSigma0TwoGammaPtEtaEle =
        new TH2F("fHistMCTruthSigma0TwoGammaPtEtaEle",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoGammaPt);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoGammaP);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoGammaPtEta);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoGammaPtEle);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoGammaPEle);
    fHistogramsSigmaMC->Add(fHistMCTruthSigma0TwoGammaPtEtaEle);

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
        new TH1F("fHistMCTruthSigma0PhotonPt", "; #it{p}_{T} [GeV/#it{c}]; Entries",
                 500, 0, 10);
    fHistMCTruthSigma0PhotonPtY =
        new TH2F("fHistMCTruthSigma0PhotonPtY", "; y; #it{p}_{T} [GeV/#it{c}]", 1000,
                 -10, 10, 500, 0, 10);
    fHistMCTruthSigma0PhotonPtEta =
        new TH2F("fHistMCTruthSigma0PhotonPtEta", "; #eta; #it{p}_{T} [GeV/#it{c}]",
                 500, -10, 10, 500, 0, 10);
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
      fHistogramsAntiSigmaMC->SetName("AntiSigmaCut_MC");
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

    fHistMCTruthAntiSigma0TwoElectronPtPos =
        new TH1F("fHistMCTruthAntiSigma0TwoElectronPtPos",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0TwoElectronPPos =
        new TH1F("fHistMCTruthAntiSigma0TwoElectronPPos",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0TwoElectronPtEtaPos =
        new TH2F("fHistMCTruthAntiSigma0TwoElectronPtEtaPos",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistMCTruthAntiSigma0TwoElectronPtNeg =
        new TH1F("fHistMCTruthAntiSigma0TwoElectronPtNeg",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0TwoElectronPNeg =
        new TH1F("fHistMCTruthAntiSigma0TwoElectronPNeg",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0TwoElectronPtEtaNeg =
        new TH2F("fHistMCTruthAntiSigma0TwoElectronPtEtaNeg",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoElectronPtPos);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoElectronPPos);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoElectronPtEtaPos);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoElectronPtNeg);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoElectronPNeg);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoElectronPtEtaNeg);

    fHistMCTruthAntiSigma0TwoGammaPt =
        new TH1F("fHistMCTruthAntiSigma0TwoGammaPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0TwoGammaP =
        new TH1F("fHistMCTruthAntiSigma0TwoGammaP",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0TwoGammaPtEta =
        new TH2F("fHistMCTruthAntiSigma0TwoGammaPtEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistMCTruthAntiSigma0TwoGammaPtEle =
        new TH1F("fHistMCTruthAntiSigma0TwoGammaPtEle",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0TwoGammaPEle =
        new TH1F("fHistMCTruthAntiSigma0TwoGammaPEle",
                 "; #it{p} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiSigma0TwoGammaPtEtaEle =
        new TH2F("fHistMCTruthAntiSigma0TwoGammaPtEtaEle",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoGammaPt);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoGammaP);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoGammaPtEta);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoGammaPtEle);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoGammaPEle);
    fHistogramsAntiSigmaMC->Add(fHistMCTruthAntiSigma0TwoGammaPtEtaEle);

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
        new TH2F("fHistMCTruthAntiSigma0PhotonPtY", "; y; #it{p}_{T} [GeV/#it{c}]",
                 1000, -10, 10, 500, 0, 10);
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
