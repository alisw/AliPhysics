#include "AliSigma0V0Cuts.h"
#include <iostream>
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "TMath.h"

ClassImp(AliSigma0V0Cuts)

    //____________________________________________________________________________________________________
    AliSigma0V0Cuts::AliSigma0V0Cuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fHistogramsBefore(nullptr),
      fHistogramsAfter(nullptr),
      fHistogramsPos(nullptr),
      fHistogramsNeg(nullptr),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fV0Container(),
      fIsLightweight(false),
      fCheckCutsMC(false),
      fV0cut(0),
      fPID(0),
      fPosPDG(0),
      fNegPDG(0),
      fIsMC(false),
      fPileUpRejectionMode(PileUpRejectionMode::BothDaughtersCombined),
      fPosPID(AliPID::kProton),
      fNegPID(AliPID::kPion),
      fV0OnFly(false),
      fK0Rejection(true),
      fUsePID(false),
      fUseArmenteros(false),
      fV0PtMin(0.f),
      fV0PtMax(999.f),
      fV0CosPAMin(0.f),
      fV0RadiusMax(999.f),
      fV0RadiusMin(0.f),
      fV0DecayVertexMax(999.f),
      fPIDnSigma(999.f),
      fEtaMax(999.f),
      fChi2Max(999.f),
      fTPCclusterMin(0.f),
      fTPCnCrossedRowsMin(0),
      fTPCratioFindable(0.f),
      fTPCfindableMin(0),
      fTPCnSharedMax(180),
      fDaughterDCAMax(999.f),
      fDaughterDCAPV(999.f),
      fK0RejectionLow(0.f),
      fK0RejectionUp(0.f),
      fArmenterosQtLow(0.f),
      fArmenterosQtUp(0.f),
      fArmenterosAlphaLow(0.f),
      fArmenterosAlphaUp(0.f),
      fLambdaSelectionLow(0.f),
      fLambdaSelectionUp(999.f),
      fPsiPairMax(999.f),
      fMCHighMultThreshold(5.f),
      fPIDResponse(nullptr),
      fHistCutBooking(nullptr),
      fHistCuts(nullptr),
      fHistNV0(nullptr),
      fHistLambdaMass(nullptr),
      fHistAntiLambdaMass(nullptr),
      fHistPhotonMass(nullptr),
      fHistPhotonMassRefit(nullptr),
      fHistK0Mass(nullptr),
      fHistV0Pt(nullptr),
      fHistV0Mass(nullptr),
      fHistV0MassPt(nullptr),
      fHistK0MassAfter(nullptr),
      fHistCosPA(nullptr),
      fHistEtaPhi(nullptr),
      fHistPsiPair(nullptr),
      fHistDecayVertexXBefore(nullptr),
      fHistDecayVertexYBefore(nullptr),
      fHistDecayVertexZBefore(nullptr),
      fHistDecayVertexXAfter(nullptr),
      fHistDecayVertexYAfter(nullptr),
      fHistDecayVertexZAfter(nullptr),
      fHistTransverseRadiusBefore(nullptr),
      fHistTransverseRadiusAfter(nullptr),
      fHistCosPABefore(nullptr),
      fHistCosPAAfter(nullptr),
      fHistDCADaughtersBefore(nullptr),
      fHistDCADaughtersAfter(nullptr),
      fHistDCA(nullptr),
      fHistDecayLength(nullptr),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistMCTruthV0PtY(nullptr),
      fHistMCTruthV0DaughterPtY(nullptr),
      fHistMCTruthV0DaughterPtYAccept(nullptr),
      fHistMCTruthPtYHighMult(nullptr),
      fHistMCTruthDaughterPtYHighMult(nullptr),
      fHistMCTruthDaughterPtYAcceptHighMult(nullptr),
      fHistMCV0Pt(nullptr),
      fHistV0Mother(nullptr),
      fHistV0MotherTrue(nullptr),
      fHistV0MassPtTrue(nullptr),
      fHistDecayVertexXTrue(nullptr),
      fHistDecayVertexYTrue(nullptr),
      fHistDecayVertexZTrue(nullptr),
      fHistTransverseRadiusTrue(nullptr),
      fHistCosPATrue(nullptr),
      fHistDCADaughtersTrue(nullptr),
      fHistArmenterosTrue(nullptr),
      fHistPsiPairTrue(nullptr),
      fHistSingleParticlePtTrue(),
      fHistSingleParticleDCAtoPVTrue(),
      fHistSingleParticleNclsTPCTrue(),
      fHistSingleParticleNclsTPCFindableTrue(),
      fHistSingleParticleNclsTPCRatioFindableTrue(),
      fHistSingleParticleNcrossedTPCTrue(),
      fHistSingleParticleNclsTPCSharedTrue(),
      fHistSingleParticleNclsITSSharedTrue(),
      fHistSingleParticleChi2True(),
      fHistSingleParticlePIDTrue(),
      fHistV0MassPtTrueSigma(nullptr),
      fHistDecayVertexXTrueSigma(nullptr),
      fHistDecayVertexYTrueSigma(nullptr),
      fHistDecayVertexZTrueSigma(nullptr),
      fHistTransverseRadiusTrueSigma(nullptr),
      fHistCosPATrueSigma(nullptr),
      fHistDCADaughtersTrueSigma(nullptr),
      fHistArmenterosTrueSigma(nullptr),
      fHistPsiPairTrueSigma(nullptr),
      fHistSingleParticlePtTrueSigma(),
      fHistSingleParticleDCAtoPVTrueSigma(),
      fHistSingleParticleNclsTPCTrueSigma(),
      fHistSingleParticleNclsTPCFindableTrueSigma(),
      fHistSingleParticleNclsTPCRatioFindableTrueSigma(),
      fHistSingleParticleNcrossedTPCTrueSigma(),
      fHistSingleParticleNclsTPCSharedTrueSigma(),
      fHistSingleParticleNclsITSSharedTrueSigma(),
      fHistSingleParticleChi2TrueSigma(),
      fHistSingleParticlePIDTrueSigma(),
      fHistV0MassPtBkg(nullptr),
      fHistDecayVertexXBkg(nullptr),
      fHistDecayVertexYBkg(nullptr),
      fHistDecayVertexZBkg(nullptr),
      fHistTransverseRadiusBkg(nullptr),
      fHistCosPABkg(nullptr),
      fHistDCADaughtersBkg(nullptr),
      fHistArmenterosBkg(nullptr),
      fHistPsiPairBkg(nullptr),
      fHistSingleParticlePtBkg(),
      fHistSingleParticleDCAtoPVBkg(),
      fHistSingleParticleNclsTPCBkg(),
      fHistSingleParticleNclsTPCFindableBkg(),
      fHistSingleParticleNclsTPCRatioFindableBkg(),
      fHistSingleParticleNcrossedTPCBkg(),
      fHistSingleParticleNclsTPCSharedBkg(),
      fHistSingleParticleNclsITSSharedBkg(),
      fHistSingleParticleChi2Bkg(),
      fHistSingleParticlePIDBkg(),
      fHistSingleParticleCuts(),
      fHistSingleParticlePt(),
      fHistSingleParticleEtaBefore(),
      fHistSingleParticleEtaAfter(),
      fHistSingleParticleChi2Before(),
      fHistSingleParticleChi2After(),
      fHistSingleParticleNclsTPCBefore(),
      fHistSingleParticleNclsTPCAfter(),
      fHistSingleParticleNclsTPCFindableBefore(),
      fHistSingleParticleNclsTPCFindableAfter(),
      fHistSingleParticleNclsTPCRatioFindableBefore(),
      fHistSingleParticleNclsTPCRatioFindableAfter(),
      fHistSingleParticleNcrossedTPCBefore(),
      fHistSingleParticleNcrossedTPCAfter(),
      fHistSingleParticleNclsTPCSharedBefore(),
      fHistSingleParticleNclsITSSharedBefore(),
      fHistSingleParticleNclsTPCSharedAfter(),
      fHistSingleParticleNclsITSSharedAfter(),
      fHistSingleParticleDCAtoPVBefore(),
      fHistSingleParticleDCAtoPVAfter(),
      fHistSingleParticlePileUp(),
      fHistSingleParticlePID() {}

//____________________________________________________________________________________________________
AliSigma0V0Cuts::AliSigma0V0Cuts(const AliSigma0V0Cuts &ref)
    : TObject(ref),
      fHistograms(nullptr),
      fHistogramsMC(nullptr),
      fHistogramsBefore(nullptr),
      fHistogramsAfter(nullptr),
      fHistogramsPos(nullptr),
      fHistogramsNeg(nullptr),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fDataBasePDG(),
      fV0Container(),
      fIsLightweight(false),
      fCheckCutsMC(false),
      fV0cut(0),
      fPID(0),
      fPosPDG(0),
      fNegPDG(0),
      fIsMC(false),
      fPileUpRejectionMode(PileUpRejectionMode::BothDaughtersCombined),
      fPosPID(AliPID::kProton),
      fNegPID(AliPID::kPion),
      fV0OnFly(false),
      fK0Rejection(true),
      fUsePID(false),
      fUseArmenteros(false),
      fV0PtMin(0.f),
      fV0PtMax(999.f),
      fV0CosPAMin(0.f),
      fV0RadiusMax(999.f),
      fV0RadiusMin(0.f),
      fV0DecayVertexMax(999.f),
      fPIDnSigma(999.f),
      fEtaMax(999.f),
      fChi2Max(999.f),
      fTPCclusterMin(0.f),
      fTPCnCrossedRowsMin(0),
      fTPCratioFindable(0.f),
      fTPCfindableMin(0),
      fTPCnSharedMax(180),
      fDaughterDCAMax(999.f),
      fDaughterDCAPV(999.f),
      fK0RejectionLow(0.f),
      fK0RejectionUp(0.f),
      fArmenterosQtLow(0.f),
      fArmenterosQtUp(0.f),
      fArmenterosAlphaLow(0.f),
      fArmenterosAlphaUp(0.f),
      fLambdaSelectionLow(0.f),
      fLambdaSelectionUp(999.f),
      fPsiPairMax(999.f),
      fMCHighMultThreshold(5.f),
      fPIDResponse(nullptr),
      fHistCutBooking(nullptr),
      fHistCuts(nullptr),
      fHistNV0(nullptr),
      fHistLambdaMass(nullptr),
      fHistAntiLambdaMass(nullptr),
      fHistPhotonMass(nullptr),
      fHistPhotonMassRefit(nullptr),
      fHistK0Mass(nullptr),
      fHistV0Pt(nullptr),
      fHistV0Mass(nullptr),
      fHistV0MassPt(nullptr),
      fHistK0MassAfter(nullptr),
      fHistCosPA(nullptr),
      fHistEtaPhi(nullptr),
      fHistPsiPair(nullptr),
      fHistDecayVertexXBefore(nullptr),
      fHistDecayVertexYBefore(nullptr),
      fHistDecayVertexZBefore(nullptr),
      fHistDecayVertexXAfter(nullptr),
      fHistDecayVertexYAfter(nullptr),
      fHistDecayVertexZAfter(nullptr),
      fHistTransverseRadiusBefore(nullptr),
      fHistTransverseRadiusAfter(nullptr),
      fHistCosPABefore(nullptr),
      fHistCosPAAfter(nullptr),
      fHistDCADaughtersBefore(nullptr),
      fHistDCADaughtersAfter(nullptr),
      fHistDCA(nullptr),
      fHistDecayLength(nullptr),
      fHistArmenterosBefore(nullptr),
      fHistArmenterosAfter(nullptr),
      fHistMCTruthV0PtY(nullptr),
      fHistMCTruthV0DaughterPtY(nullptr),
      fHistMCTruthV0DaughterPtYAccept(nullptr),
      fHistMCTruthPtYHighMult(nullptr),
      fHistMCTruthDaughterPtYHighMult(nullptr),
      fHistMCTruthDaughterPtYAcceptHighMult(nullptr),
      fHistMCV0Pt(nullptr),
      fHistV0Mother(nullptr),
      fHistV0MotherTrue(nullptr),
      fHistV0MassPtTrue(nullptr),
      fHistDecayVertexXTrue(nullptr),
      fHistDecayVertexYTrue(nullptr),
      fHistDecayVertexZTrue(nullptr),
      fHistTransverseRadiusTrue(nullptr),
      fHistCosPATrue(nullptr),
      fHistDCADaughtersTrue(nullptr),
      fHistArmenterosTrue(nullptr),
      fHistPsiPairTrue(nullptr),
      fHistSingleParticlePtTrue(),
      fHistSingleParticleDCAtoPVTrue(),
      fHistSingleParticleNclsTPCTrue(),
      fHistSingleParticleNclsTPCFindableTrue(),
      fHistSingleParticleNclsTPCRatioFindableTrue(),
      fHistSingleParticleNcrossedTPCTrue(),
      fHistSingleParticleNclsTPCSharedTrue(),
      fHistSingleParticleNclsITSSharedTrue(),
      fHistSingleParticleChi2True(),
      fHistSingleParticlePIDTrue(),
      fHistV0MassPtTrueSigma(nullptr),
      fHistDecayVertexXTrueSigma(nullptr),
      fHistDecayVertexYTrueSigma(nullptr),
      fHistDecayVertexZTrueSigma(nullptr),
      fHistTransverseRadiusTrueSigma(nullptr),
      fHistCosPATrueSigma(nullptr),
      fHistDCADaughtersTrueSigma(nullptr),
      fHistArmenterosTrueSigma(nullptr),
      fHistPsiPairTrueSigma(nullptr),
      fHistSingleParticlePtTrueSigma(),
      fHistSingleParticleDCAtoPVTrueSigma(),
      fHistSingleParticleNclsTPCTrueSigma(),
      fHistSingleParticleNclsTPCFindableTrueSigma(),
      fHistSingleParticleNclsTPCRatioFindableTrueSigma(),
      fHistSingleParticleNcrossedTPCTrueSigma(),
      fHistSingleParticleNclsTPCSharedTrueSigma(),
      fHistSingleParticleNclsITSSharedTrueSigma(),
      fHistSingleParticleChi2TrueSigma(),
      fHistSingleParticlePIDTrueSigma(),
      fHistV0MassPtBkg(nullptr),
      fHistDecayVertexXBkg(nullptr),
      fHistDecayVertexYBkg(nullptr),
      fHistDecayVertexZBkg(nullptr),
      fHistTransverseRadiusBkg(nullptr),
      fHistCosPABkg(nullptr),
      fHistDCADaughtersBkg(nullptr),
      fHistArmenterosBkg(nullptr),
      fHistPsiPairBkg(nullptr),
      fHistSingleParticlePtBkg(),
      fHistSingleParticleDCAtoPVBkg(),
      fHistSingleParticleNclsTPCBkg(),
      fHistSingleParticleNclsTPCFindableBkg(),
      fHistSingleParticleNclsTPCRatioFindableBkg(),
      fHistSingleParticleNcrossedTPCBkg(),
      fHistSingleParticleNclsTPCSharedBkg(),
      fHistSingleParticleNclsITSSharedBkg(),
      fHistSingleParticleChi2Bkg(),
      fHistSingleParticlePIDBkg(),
      fHistSingleParticleCuts(),
      fHistSingleParticlePt(),
      fHistSingleParticleEtaBefore(),
      fHistSingleParticleEtaAfter(),
      fHistSingleParticleChi2Before(),
      fHistSingleParticleChi2After(),
      fHistSingleParticleNclsTPCBefore(),
      fHistSingleParticleNclsTPCAfter(),
      fHistSingleParticleNclsTPCFindableBefore(),
      fHistSingleParticleNclsTPCFindableAfter(),
      fHistSingleParticleNclsTPCRatioFindableBefore(),
      fHistSingleParticleNclsTPCRatioFindableAfter(),
      fHistSingleParticleNcrossedTPCBefore(),
      fHistSingleParticleNcrossedTPCAfter(),
      fHistSingleParticleNclsTPCSharedBefore(),
      fHistSingleParticleNclsITSSharedBefore(),
      fHistSingleParticleNclsTPCSharedAfter(),
      fHistSingleParticleNclsITSSharedAfter(),
      fHistSingleParticleDCAtoPVBefore(),
      fHistSingleParticleDCAtoPVAfter(),
      fHistSingleParticlePileUp(),
      fHistSingleParticlePID() {}

//____________________________________________________________________________________________________
AliSigma0V0Cuts &AliSigma0V0Cuts::operator=(const AliSigma0V0Cuts &ref) {
  // Assignment operator
  if (this == &ref) return *this;

  return (*this);
}

//____________________________________________________________________________________________________
AliSigma0V0Cuts::~AliSigma0V0Cuts() {
  if (fPIDResponse) delete fPIDResponse;
}

//____________________________________________________________________________________________________
AliSigma0V0Cuts *AliSigma0V0Cuts::LambdaCuts() {
  AliSigma0V0Cuts *v0Cuts = new AliSigma0V0Cuts();
  v0Cuts->SetV0OnFlyStatus(false);
  v0Cuts->SetV0PtMin(0.3);
  v0Cuts->SetV0CosPAMin(0.99);
  v0Cuts->SetV0RadiusMax(100.f);
  v0Cuts->SetV0RadiusMin(0.2);
  v0Cuts->SetV0DecayVertexMax(100.f);
  v0Cuts->SetPIDnSigma(5.f);
  v0Cuts->SetTPCclusterMin(70.f);
  v0Cuts->SetEtaMax(0.9);
  v0Cuts->SetDaughterDCAMax(1.5);
  v0Cuts->SetDaughterDCAtoPV(0.05);
  v0Cuts->SetLambdaSelection(1.115683 - 0.006, 1.115683 + 0.006);
  v0Cuts->SetPileUpRejectionMode(OneDaughterCombined);
  v0Cuts->SetArmenterosCut(0.01, 0.12, 0.3, 0.95);
  return v0Cuts;
}

//____________________________________________________________________________________________________
AliSigma0V0Cuts *AliSigma0V0Cuts::PhotonCuts() {
  AliSigma0V0Cuts *v0Cuts = new AliSigma0V0Cuts();
  v0Cuts->SetV0OnFlyStatus(false);
  v0Cuts->SetV0PtMin(0.2);
  v0Cuts->SetV0PtMax(2.f);
  v0Cuts->SetV0CosPAMin(0.99);
  v0Cuts->SetV0DecayVertexMax(150.f);
  v0Cuts->SetV0RadiusMax(100.f);
  v0Cuts->SetV0RadiusMin(1.f);
  v0Cuts->SetTPCRatioFindable(0.6f);
  v0Cuts->SetEtaMax(0.8);
  v0Cuts->SetDaughterDCAMax(1.5);
  v0Cuts->SetDaughterDCAtoPV(0.05);
  v0Cuts->SetArmenterosCut(0, 0.06, -1, 1);
  v0Cuts->SetK0Rejection(0., 0.);
  v0Cuts->SetLambdaSelection(0., 0.);
  v0Cuts->SetPsiPairMax(0.2);
  v0Cuts->SetPileUpRejectionMode(None);
  v0Cuts->SetChi2Max(4);
  return v0Cuts;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::SelectV0(AliVEvent *inputEvent, AliMCEvent *mcEvent) {
  // Get the PID manager from the input event handler
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      static_cast<AliInputEventHandler *>(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fMCEvent = mcEvent;
  fInputEvent = static_cast<AliESDEvent *>(inputEvent);

  fV0Container.clear();

  if (fIsMC) {
    ProcessMC();
    if (fCheckCutsMC) {
      CheckCutsMC();
    }
  }

  for (int iV0 = 0; iV0 < fInputEvent->GetNumberOfV0s(); ++iV0) {
    AliESDv0 *v0 = fInputEvent->GetV0(iV0);

    // V0 Quality
    if (!V0QualityCuts(v0)) continue;

    // Daughter tracks
    AliESDtrack *pos = fInputEvent->GetTrack(v0->GetPindex());
    AliESDtrack *neg = fInputEvent->GetTrack(v0->GetNindex());
    if (!pos || !neg) continue;

    if (pos->Charge() < 0) {
      pos = neg;
      neg = fInputEvent->GetTrack(v0->GetPindex());
    }

    // Single Particle Quality
    if (!SingleParticleQualityCuts(pos) || !SingleParticleQualityCuts(neg))
      continue;
    if (!fIsLightweight)
      if (!fIsLightweight) fHistCuts->Fill(6);

    // Pile-up rejection
    if (!PileUpRejection(pos, neg)) continue;
    if (!fIsLightweight) fHistCuts->Fill(7);

    // Topological Selection
    if (!V0TopologicalSelection(v0)) continue;

    // PID
    if (!V0PID(v0, pos, neg)) continue;

    PlotMasses(v0);

    // V0 Selection
    if (TMath::Abs(fPID) == 3122 && !LambdaSelection(v0)) continue;
    if (TMath::Abs(fPID) == 22 && !PhotonSelection(v0)) continue;

    AliSigma0ParticleV0 v0Candidate(v0, pos, neg,
                                    fInputEvent->GetPrimaryVertex(), fPID,
                                    fInputEvent->GetMagneticField(), fMCEvent);

    if (fIsMC) {
      int label = v0Candidate.MatchToMC(fMCEvent, fPID, {{fPosPDG, fNegPDG}});

      AliMCParticle *mcParticle =
          static_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
      if (mcParticle){
      v0Candidate.ProcessMCInfo(mcParticle, fMCEvent);

      fHistMCV0Pt->Fill(v0->Pt());

      const int pdgMother = static_cast<AliMCParticle *>(
                                fMCEvent->GetTrack(mcParticle->GetMother()))
                                ->PdgCode();
      fHistV0Mother->Fill(v0->Pt(), TMath::Abs(pdgMother));
      }
    }

    v0Candidate.SetPDGMass(fDataBasePDG.GetParticle(fPID)->Mass());
    if (TMath::Abs(fPID) == 22) {
      const float massPhoton = ComputePhotonMassRefit(v0);
      v0Candidate.SetMass(massPhoton);
      v0Candidate.SetRecMass(massPhoton);
    }

    fV0Container.push_back(v0Candidate);
  }
  if (!fIsLightweight) fHistNV0->Fill(fV0Container.size());
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0QualityCuts(const AliESDv0 *v0) const {
  if (!fIsLightweight) fHistCuts->Fill(0);

  // v0 is neutral
  if (v0->Charge() != 0) return false;
  if (!fIsLightweight) fHistCuts->Fill(1);

  // on-fly status
  if (fV0OnFly) {
    // online v0
    if (!(v0->GetOnFlyStatus())) return false;
  } else {
    // offline v0
    if (v0->GetOnFlyStatus()) return false;
  }
  if (!fIsLightweight) fHistCuts->Fill(2);

  // pt min cut
  if (v0->Pt() < fV0PtMin) return false;
  if (!fIsLightweight) fHistCuts->Fill(3);

  // pt max cut
  if (v0->Pt() > fV0PtMax) return false;
  if (!fIsLightweight) fHistCuts->Fill(4);

  // check for daughter tracks
  AliESDtrack *pos = fInputEvent->GetTrack(v0->GetPindex());
  AliESDtrack *neg = fInputEvent->GetTrack(v0->GetNindex());
  if (!pos || !neg) return false;
  if (pos->Charge() == neg->Charge()) return false;
  if (!fIsLightweight) fHistCuts->Fill(5);

  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0PID(const AliESDv0 *v0, const AliESDtrack *pos,
                            const AliESDtrack *neg) const {
  float posProb = 999.f;
  float negProb = 999.f;

  const float armAlpha = v0->AlphaV0();
  const float armQt = v0->PtArmV0();
  if (!fIsLightweight) fHistArmenterosBefore->Fill(armAlpha, armQt);

  // When true, we use the full PID capabilities, when false just rough
  // Armenteros-Podolandksi cut
  if (fUsePID) {
    bool isParticle = (SingleParticlePID(pos, fPosPID, posProb) &&
                       SingleParticlePID(neg, fNegPID, negProb));

    if (!isParticle) return false;
  }
  if (fUseArmenteros) {
    // Armenteros cut
    if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) return false;
    float prefactorAlpha = (fPID > 0) ? 1.f : -1.f;  // for anti-particles
                                                     // negative to compensate
                                                     // for sign change of qt
    if (armAlpha * prefactorAlpha > fArmenterosAlphaUp ||
        armAlpha * prefactorAlpha < fArmenterosAlphaLow)
      return false;
  }

  if (!fIsLightweight) fHistCuts->Fill(15);
  if (!fIsLightweight) {
    fHistArmenterosAfter->Fill(armAlpha, armQt);
    PlotSingleParticlePID(pos, fPosPID);
    PlotSingleParticlePID(neg, fNegPID);
  }

  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::SingleParticlePID(const AliVTrack *track,
                                        AliPID::EParticleType particle,
                                        float &prob) const {
  AliPIDResponse::EDetPidStatus statusPosTPC =
      fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);

  // TPC signal must be available
  if (!(AliPIDResponse::kDetPidOk == statusPosTPC)) return false;

  prob = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, particle));

  // Proton PID cut
  if (prob > fPIDnSigma) return false;
  return true;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::PlotSingleParticlePID(
    const AliVTrack *track, AliPID::EParticleType particle) const {
  AliPIDResponse::EDetPidStatus statusPosTPC =
      fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);

  // TPC signal must be available
  if (!(AliPIDResponse::kDetPidOk == statusPosTPC)) return;
  const float nSigma = fPIDResponse->NumberOfSigmasTPC(track, particle);
  const int charge = track->Charge();
  const int histoPrefix = (charge > 0) ? 0 : 1;
  if (!fIsLightweight)
    fHistSingleParticlePID[histoPrefix]->Fill(track->P(), nSigma);
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::SingleParticleQualityCuts(AliESDtrack *track) const {
  const float pt = track->Pt();
  const int charge = track->Charge();
  const float eta = track->Eta();
  const short nClsTPC = track->GetTPCNcls();
  const float nCrossedRows = track->GetTPCClusterInfo(2, 1);
  const short nFindable = track->GetTPCNclsF();
  const short nClsSharedTPC = track->GetTPCnclsS();
  const float ratioFindable =
      nCrossedRows / (static_cast<float>(nFindable) + 1.e-3);
  const float chi2 =
      (nClsTPC > 5) ? (track->GetTPCchi2() / float(nClsTPC - 5)) : -1.;
  /// see AliAnalysisTaskESDfilter::Chi2perNDF
  int nClsSharedITS = 0;
  for (int i = 0; i < 6; ++i) {
    if (track->HasSharedPointOnITSLayer(i)) ++nClsSharedITS;
  }

  const float magField = fInputEvent->GetMagneticField();
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  const float dcaDaughterToPV =
      TMath::Abs(track->GetD(vertex->GetX(), vertex->GetY(), magField));

  const int histoPrefix = (charge > 0) ? 0 : 1;
  int qaHistoCounter = 0;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (!fIsLightweight) {
    fHistSingleParticleEtaBefore[histoPrefix]->Fill(pt, eta);
    fHistSingleParticleChi2Before[histoPrefix]->Fill(pt, chi2);
    fHistSingleParticleNclsTPCBefore[histoPrefix]->Fill(pt, nClsTPC);
    fHistSingleParticleNclsTPCFindableBefore[histoPrefix]->Fill(pt, nFindable);
    fHistSingleParticleNclsTPCSharedBefore[histoPrefix]->Fill(pt,
                                                              nClsSharedTPC);
    fHistSingleParticleNclsITSSharedBefore[histoPrefix]->Fill(pt,
                                                              nClsSharedITS);
    fHistSingleParticleNclsTPCRatioFindableBefore[histoPrefix]->Fill(
        pt, ratioFindable);
    fHistSingleParticleNcrossedTPCBefore[histoPrefix]->Fill(pt, nCrossedRows);
    fHistSingleParticleDCAtoPVBefore[histoPrefix]->Fill(pt, dcaDaughterToPV);
  }

  // Max eta cut
  if (TMath::Abs(eta) > fEtaMax) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // Max chi2 cut
  if (chi2 > fChi2Max) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC nCluster cut
  if (nClsTPC < fTPCclusterMin) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC Crossed rows
  if (nCrossedRows < fTPCnCrossedRowsMin) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC Ratio Findable
  if (ratioFindable < fTPCratioFindable) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC nCls Findable
  if (nFindable < fTPCfindableMin) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC nCls Shared
  if (nClsSharedTPC > fTPCnSharedMax) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // Minimal distance of daughters to primary vertex
  if (dcaDaughterToPV < fDaughterDCAPV) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (!(track->GetStatus() & AliESDtrack::kTPCrefit)) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (track->GetKinkIndex(0) > 0) return false;
  if (!fIsLightweight)
    fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (!fIsLightweight) {
    fHistSingleParticlePt[histoPrefix]->Fill(pt);
    fHistSingleParticleEtaAfter[histoPrefix]->Fill(pt, eta);
    fHistSingleParticleChi2After[histoPrefix]->Fill(pt, chi2);
    fHistSingleParticleNclsTPCAfter[histoPrefix]->Fill(pt, nClsTPC);
    fHistSingleParticleNclsTPCFindableAfter[histoPrefix]->Fill(pt, nFindable);
    fHistSingleParticleNclsTPCSharedAfter[histoPrefix]->Fill(pt, nClsSharedTPC);
    fHistSingleParticleNclsITSSharedAfter[histoPrefix]->Fill(pt, nClsSharedITS);
    fHistSingleParticleNclsTPCRatioFindableAfter[histoPrefix]->Fill(
        pt, ratioFindable);
    fHistSingleParticleNcrossedTPCAfter[histoPrefix]->Fill(pt, nCrossedRows);
    fHistSingleParticleDCAtoPVAfter[histoPrefix]->Fill(pt, dcaDaughterToPV);
  }
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::PileUpRejection(AliESDtrack *pos,
                                      AliESDtrack *neg) const {
  const float posPt = pos->Pt();
  const float negPt = neg->Pt();

  bool posTrackITS =
      (pos->HasPointOnITSLayer(0) || pos->HasPointOnITSLayer(1) ||
       pos->HasPointOnITSLayer(4) || pos->HasPointOnITSLayer(5));
  bool negTrackITS =
      (neg->HasPointOnITSLayer(0) || pos->HasPointOnITSLayer(1) ||
       neg->HasPointOnITSLayer(4) || pos->HasPointOnITSLayer(5));
  bool posTrackTOF = pos->GetTOFBunchCrossing() == 0;
  bool negTrackTOF = neg->GetTOFBunchCrossing() == 0;

  bool posTrackCombined = (posTrackITS || posTrackTOF);
  bool negTrackCombined = (negTrackITS || negTrackTOF);

  if (!fIsLightweight) {
    if (posTrackITS) fHistSingleParticlePileUp[0]->Fill(0.f, posPt);
    if (posTrackTOF) fHistSingleParticlePileUp[0]->Fill(1, posPt);
    if (posTrackCombined) fHistSingleParticlePileUp[0]->Fill(2, posPt);
    if (!posTrackITS && !posTrackTOF)
      fHistSingleParticlePileUp[0]->Fill(3, posPt);
    if (negTrackITS) fHistSingleParticlePileUp[1]->Fill(0.f, negPt);
    if (negTrackTOF) fHistSingleParticlePileUp[1]->Fill(1, negPt);
    if (negTrackCombined) fHistSingleParticlePileUp[1]->Fill(2, negPt);
    if (!negTrackITS && !negTrackTOF)
      fHistSingleParticlePileUp[1]->Fill(3, negPt);
  }

  switch (fPileUpRejectionMode) {
    case PileUpRejectionMode::BothDaughtersCombined: {
      if (posTrackCombined && negTrackCombined) {
        return true;
        break;
      } else {
        return false;
        break;
      }
    }
    case PileUpRejectionMode::OneDaughterCombined: {
      if (posTrackCombined || negTrackCombined) {
        return true;
        break;
      } else {
        return false;
        break;
      }
    }
    case PileUpRejectionMode::None: {
      return true;
      break;
    }
    case PileUpRejectionMode::BothDaughtersITSonly: {
      if (posTrackITS && negTrackITS) {
        return true;
        break;
      } else {
        return false;
        break;
      }
    }
    case PileUpRejectionMode::BothDaughtersTOFonly: {
      if (posTrackTOF && negTrackTOF) {
        return true;
        break;
      } else {
        return false;
        break;
      }
    }
    case PileUpRejectionMode::OneDaughterITSonly: {
      if (posTrackITS || negTrackITS) {
        return true;
        break;
      } else {
        return false;
        break;
      }
    }
    case PileUpRejectionMode::OneDaughterTOFonly: {
      if (posTrackTOF || negTrackTOF) {
        return true;
        break;
      } else {
        return false;
        break;
      }
    }
  }
  return false;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0TopologicalSelection(const AliESDv0 *v0) const {
  // Get the coordinates of the primary vertex
  const float pt = v0->Pt();

  Double_t xPV = fInputEvent->GetPrimaryVertex()->GetX();
  Double_t yPV = fInputEvent->GetPrimaryVertex()->GetY();
  Double_t zPV = fInputEvent->GetPrimaryVertex()->GetZ();
  Double_t decayVertexV0[3];

  v0->GetXYZ(decayVertexV0[0], decayVertexV0[1], decayVertexV0[2]);

  // Calculate vertex variables:
  const float point = v0->GetV0CosineOfPointingAngle(xPV, yPV, zPV);
  const float dcaV0Dau = v0->GetDcaV0Daughters();
  const float dcaPrim = v0->GetD(xPV, yPV, zPV);
  const float lenDecay = std::sqrt(decayVertexV0[0] * decayVertexV0[0] +
                                   decayVertexV0[1] * decayVertexV0[1] +
                                   decayVertexV0[2] * decayVertexV0[2]);
  const float transverseRadius = std::sqrt(decayVertexV0[0] * decayVertexV0[0] +
                                           decayVertexV0[1] * decayVertexV0[1]);

  if (!fIsLightweight) {
    fHistCosPABefore->Fill(pt, point);
    fHistDecayVertexXBefore->Fill(pt, decayVertexV0[0]);
    fHistDecayVertexYBefore->Fill(pt, decayVertexV0[1]);
    fHistDecayVertexZBefore->Fill(pt, decayVertexV0[2]);
    fHistTransverseRadiusBefore->Fill(pt, transverseRadius);
    fHistDCADaughtersBefore->Fill(pt, dcaV0Dau);
  }

  // Position of the decay vertex x, y & z
  if (TMath::Abs(decayVertexV0[0]) > fV0DecayVertexMax) return false;
  if (!fIsLightweight)
    if (!fIsLightweight) fHistCuts->Fill(8);
  if (TMath::Abs(decayVertexV0[1]) > fV0DecayVertexMax) return false;
  if (!fIsLightweight)
    if (!fIsLightweight) fHistCuts->Fill(9);
  if (TMath::Abs(decayVertexV0[2]) > fV0DecayVertexMax) return false;
  if (!fIsLightweight) fHistCuts->Fill(10);

  // Transverse decay radius min & max
  if (transverseRadius < fV0RadiusMin) return false;
  if (!fIsLightweight) fHistCuts->Fill(11);

  if (transverseRadius > fV0RadiusMax) return false;
  if (!fIsLightweight) fHistCuts->Fill(12);

  // DCA of the daughter tracks at the decay vertex
  if (dcaV0Dau > fDaughterDCAMax) return false;
  if (!fIsLightweight) fHistCuts->Fill(13);

  // Cos pointing angle
  if (point < fV0CosPAMin) return false;
  if (!fIsLightweight) fHistCuts->Fill(14);

  if (!fIsLightweight) {
    fHistCosPA->Fill(pt, point);
    fHistCosPAAfter->Fill(pt, point);
    fHistDecayVertexXAfter->Fill(pt, decayVertexV0[0]);
    fHistDecayVertexYAfter->Fill(pt, decayVertexV0[1]);
    fHistDecayVertexZAfter->Fill(pt, decayVertexV0[2]);
    fHistTransverseRadiusAfter->Fill(pt, transverseRadius);
    fHistDCADaughtersAfter->Fill(pt, dcaV0Dau);
    fHistDCA->Fill(pt, dcaPrim);
    fHistDecayLength->Fill(pt, lenDecay);
  }
  return true;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::PlotMasses(AliESDv0 *v0) const {
  v0->ChangeMassHypothesis(310);
  const float massK0 = v0->GetEffMass();
  v0->ChangeMassHypothesis(3122);
  const float massLambda = v0->GetEffMass();
  v0->ChangeMassHypothesis(-3122);
  const float massAntiLambda = v0->GetEffMass();
  const float massPhoton = ComputePhotonMass(v0);
  const float massPhotonRefit = ComputePhotonMassRefit(v0);

  if (!fIsLightweight) {
    fHistK0Mass->Fill(massK0);
    fHistLambdaMass->Fill(massLambda);
    fHistAntiLambdaMass->Fill(massAntiLambda);
    fHistPhotonMass->Fill(massPhoton);
    fHistPhotonMassRefit->Fill(massPhotonRefit);
  }
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::LambdaSelection(AliESDv0 *v0) const {
  v0->ChangeMassHypothesis(fPID);
  const float massLambda = v0->GetEffMass();
  v0->ChangeMassHypothesis(310);
  const float massK0 = v0->GetEffMass();

  // K0 rejection cut
  if (fK0Rejection && massK0 > fK0RejectionLow && massK0 < fK0RejectionUp) {
    return false;
  }
  if (!fIsLightweight) fHistCuts->Fill(16);

  fHistV0MassPt->Fill(v0->Pt(), massLambda);
  if (!fIsLightweight) {
    fHistK0MassAfter->Fill(massK0);
    fHistEtaPhi->Fill(v0->Eta(), v0->Phi());
  }

  if (!fIsLightweight) fHistV0Mass->Fill(massLambda);

  // Lambda selection cut
  if (massLambda < fLambdaSelectionLow || massLambda > fLambdaSelectionUp) {
    return false;
  }

  if (!fIsLightweight) fHistV0Pt->Fill(v0->Pt());
  if (!fIsLightweight) fHistCuts->Fill(17);
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::PhotonSelection(AliESDv0 *v0) const {
  v0->ChangeMassHypothesis(3122);
  const float massLambda = v0->GetEffMass();
  v0->ChangeMassHypothesis(310);
  const float massK0 = v0->GetEffMass();
  const float massPhoton = ComputePhotonMassRefit(v0);

  const float psiPair = ComputePsiPair(v0);
  if (!fIsLightweight) fHistPsiPair->Fill(v0->Pt(), psiPair);

  if (TMath::Abs(psiPair) > fPsiPairMax) return false;
  if (!fIsLightweight) fHistCuts->Fill(16);

  // K0 rejection cut
  if (fK0Rejection && massK0 > fK0RejectionLow && massK0 < fK0RejectionUp) {
    return false;
  }
  if (!fIsLightweight) fHistCuts->Fill(17);

  fHistV0MassPt->Fill(v0->Pt(), massPhoton);
  if (!fIsLightweight) {
    fHistK0MassAfter->Fill(massK0);
    fHistEtaPhi->Fill(v0->Eta(), v0->Phi());
  }

  // Lambda rejection cut
  if (massLambda > fLambdaSelectionLow && massLambda < fLambdaSelectionUp) {
    return false;
  }

  if (!fIsLightweight) {
    fHistV0Mass->Fill(massPhoton);
    fHistV0Pt->Fill(v0->Pt());
  }

  if (!fIsLightweight) fHistCuts->Fill(18);
  return true;
}

//____________________________________________________________________________________________________
float AliSigma0V0Cuts::ComputePhotonMass(const AliESDv0 *v0) const {
  static float eleMass = fDataBasePDG.GetParticle(kElectron)->Mass();
  float nmass = eleMass, pmass = eleMass;

  Double_t pxn, pyn, pzn;
  v0->GetNPxPyPz(pxn, pyn, pzn);
  Double_t pxp, pyp, pzp;
  v0->GetPPxPyPz(pxp, pyp, pzp);

  float en = std::sqrt(nmass * nmass + pxn * pxn + pyn * pyn + pzn * pzn);
  float ep = std::sqrt(pmass * pmass + pxp * pxp + pyp * pyp + pzp * pzp);
  float pxl = pxn + pxp, pyl = pyn + pyp, pzl = pzn + pzp;
  float pl = std::sqrt(pxl * pxl + pyl * pyl + pzl * pzl);

  return std::sqrt((en + ep) * (en + ep) - pl * pl);
}

//____________________________________________________________________________________________________
float AliSigma0V0Cuts::ComputePhotonMassRefit(const AliESDv0 *v0) const {
  const AliExternalTrackParam *parPos = v0->GetParamP();
  const AliExternalTrackParam *parNeg = v0->GetParamN();
  Double_t mass = -99.0, mass_width = -99.0, Pt = -99.0, Pt_width = -99.0;
  AliKFParticle neg(*(parNeg), 11);
  AliKFParticle pos(*(parPos), -11);
  AliKFParticle fCurrentMotherKFForMass(neg, pos);
  fCurrentMotherKFForMass.GetMass(mass, mass_width);
  fCurrentMotherKFForMass.GetPt(Pt, Pt_width);

  AliKFConversionPhoton *mother = new AliKFConversionPhoton();
  mother->ConstructGamma(neg, pos);

  const float invMassRefit = mother->M();
  delete mother;
  return invMassRefit;
}

/////________________________________________________________________________
float AliSigma0V0Cuts::ComputePsiPair(const AliESDv0 *v0) const {
  const float magneticField = fInputEvent->GetMagneticField();
  Double_t decayVertexV0[3];
  v0->GetXYZ(decayVertexV0[0], decayVertexV0[1], decayVertexV0[2]);

  Double_t m1[3] = {0, 0, 0};
  Double_t m2[3] = {0, 0, 0};
  v0->GetPPxPyPz(m1[0], m1[1], m1[2]);
  v0->GetNPxPyPz(m2[0], m2[1], m2[2]);

  const float deltat =
      std::atan(m2[2] / (std::sqrt(m2[0] * m2[0] + m2[1] * m2[1]) + 1.e-13)) -
      std::atan(m1[2] / (std::sqrt(m1[0] * m1[0] + m1[1] * m1[1]) +
                         1.e-13));  // difference of angles of the two
                                    // daughter tracks with z-axis

  // radius to which tracks shall be propagated

  // WHY 50???

  const float radiussum = std::sqrt(decayVertexV0[0] * decayVertexV0[0] +
                                    decayVertexV0[1] * decayVertexV0[1]) +
                          50;

  Double_t mom1Prop[3];
  Double_t mom2Prop[3];
  const AliExternalTrackParam *d1 = v0->GetParamP();
  const AliExternalTrackParam *d2 = v0->GetParamN();
  AliExternalTrackParam nt(*d1), pt(*d2);

  float psiPair = 4.;
  // propagate tracks to the outside
  if (nt.PropagateTo(radiussum, magneticField) == 0) {
    psiPair = -5.;
  }
  if (pt.PropagateTo(radiussum, magneticField) == 0) {
    psiPair = -5.;
  }
  // Get momentum vectors of tracks after propagation
  pt.GetPxPyPz(mom1Prop);
  nt.GetPxPyPz(mom2Prop);

  const float pEle =
      std::sqrt(mom2Prop[0] * mom2Prop[0] + mom2Prop[1] * mom2Prop[1] +
                mom2Prop[2] * mom2Prop[2]);
  const float pPos =
      std::sqrt(mom1Prop[0] * mom1Prop[0] + mom1Prop[1] * mom1Prop[1] +
                mom1Prop[2] * mom1Prop[2]);

  const float scalarproduct = mom1Prop[0] * mom2Prop[0] +
                              mom1Prop[1] * mom2Prop[1] +
                              mom1Prop[2] * mom2Prop[2];

  // Angle between propagated daughter tracks
  const float chipair = std::acos(scalarproduct / (pEle * pPos));

  psiPair = std::asin(deltat / chipair);
  return psiPair;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::ProcessMC() const {
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
    if (!mcParticle->IsPhysicalPrimary()) continue;
    if (mcParticle->PdgCode() != fPID) continue;
    fHistMCTruthV0PtY->Fill(mcParticle->Y(), mcParticle->Pt());
    if (lPercentile < fMCHighMultThreshold) {
      fHistMCTruthPtYHighMult->Fill(mcParticle->Y(), mcParticle->Pt());
    }

    if (!CheckDaughters(mcParticle)) continue;
    fHistMCTruthV0DaughterPtY->Fill(mcParticle->Y(), mcParticle->Pt());
    if (lPercentile < fMCHighMultThreshold) {
      fHistMCTruthDaughterPtYHighMult->Fill(mcParticle->Y(), mcParticle->Pt());
    }

    if (!CheckDaughtersInAcceptance(mcParticle)) continue;
    fHistMCTruthV0DaughterPtYAccept->Fill(mcParticle->Y(), mcParticle->Pt());
    if (lPercentile < fMCHighMultThreshold) {
      fHistMCTruthDaughterPtYAcceptHighMult->Fill(mcParticle->Y(),
                                                  mcParticle->Pt());
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::PhotonQA(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                               const TClonesArray *photons) {
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      static_cast<AliInputEventHandler *>(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  fHistNV0->Fill(photons->GetEntriesFast());

  for (int iGamma = 0; iGamma < photons->GetEntriesFast(); ++iGamma) {
    auto *PhotonCandidate =
        dynamic_cast<AliAODConversionPhoton *>(photons->At(iGamma));
    if (!PhotonCandidate) continue;

    const float pt = PhotonCandidate->GetPhotonPt();
    const float invMass = PhotonCandidate->GetPhotonMass();
    // commpute CPA
    double momV0[3] = {0, 0, 0};
    momV0[0] = PhotonCandidate->Px();
    momV0[1] = PhotonCandidate->Py();
    momV0[2] = PhotonCandidate->Pz();
    double conv[3]{0, 0, 0};
    conv[0] = PhotonCandidate->GetConversionX();
    conv[1] = PhotonCandidate->GetConversionY();
    conv[2] = PhotonCandidate->GetConversionZ();

    const AliVVertex *vertex = inputEvent->GetPrimaryVertex();
    Double_t xPV = vertex->GetX();
    Double_t yPV = vertex->GetY();
    Double_t zPV = vertex->GetZ();

    double PosV0[3] = {conv[0] - xPV, conv[1] - yPV, conv[2] - zPV};
    // Recalculated V0 Position vector

    double momV02 =
        momV0[0] * momV0[0] + momV0[1] * momV0[1] + momV0[2] * momV0[2];
    double PosV02 =
        PosV0[0] * PosV0[0] + PosV0[1] * PosV0[1] + PosV0[2] * PosV0[2];

    double cosinePointingAngle =
        (momV02 * PosV02 > 0.0)
            ? (PosV0[0] * momV0[0] + PosV0[1] * momV0[1] +
               PosV0[2] * momV0[2]) /
                  std::sqrt(momV02 * PosV02)
            : -999.f;

    int v0Index = PhotonCandidate->GetV0Index();
    AliESDv0 *v0 = static_cast<AliESDEvent *>(inputEvent)->GetV0(v0Index);
    if (!v0) continue;
    v0->ChangeMassHypothesis(3122);
    const float massLambda = v0->GetEffMass();
    v0->ChangeMassHypothesis(310);
    const float massK0 = v0->GetEffMass();
    v0->ChangeMassHypothesis(-3122);
    const float massAntiLambda = v0->GetEffMass();

    fHistK0Mass->Fill(massK0);
    fHistLambdaMass->Fill(massLambda);
    fHistAntiLambdaMass->Fill(massAntiLambda);

    // Calculate vertex variables:
    const float dcaPrim = v0->GetD(xPV, yPV, zPV);

    const float armAlpha = v0->AlphaV0();
    const float armQt = v0->PtArmV0();
    fHistArmenterosAfter->Fill(armAlpha, armQt);

    auto pos =
        (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel1());
    auto neg =
        (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel2());
    if (!pos || !neg) continue;

    if (fIsMC) {
      int label = PhotonCandidate->GetMCParticleLabel(mcEvent);
      if (label > 0) {
        AliMCParticle *mcParticle =
            static_cast<AliMCParticle *>(mcEvent->GetTrack(label));
        if (mcParticle && mcParticle->PdgCode() == fPID) {
          fHistMCV0Pt->Fill(v0->Pt());
          const int pdgMother = static_cast<AliMCParticle *>(
                                    mcEvent->GetTrack(mcParticle->GetMother()))
                                    ->PdgCode();
          fHistV0Mother->Fill(v0->Pt(), TMath::Abs(pdgMother));
        }
      }
    }

    fHistV0MassPt->Fill(pt, invMass);

    if (!fIsLightweight) {
      fHistPhotonMass->Fill(PhotonCandidate->M());
      fHistPhotonMassRefit->Fill(PhotonCandidate->GetInvMassPair());
      fHistV0Pt->Fill(pt);
      fHistV0Mass->Fill(invMass);
      fHistCosPA->Fill(pt, cosinePointingAngle);
      fHistEtaPhi->Fill(PhotonCandidate->GetPhotonEta(),
                        PhotonCandidate->GetPhotonPhi());
      fHistPsiPair->Fill(pt, PhotonCandidate->GetPsiPair());
      fHistDecayVertexXAfter->Fill(pt, conv[0]);
      fHistDecayVertexYAfter->Fill(pt, conv[1]);
      fHistDecayVertexZAfter->Fill(pt, conv[2]);
      fHistTransverseRadiusAfter->Fill(
          pt, std::sqrt(conv[0] * conv[0] + conv[1] * conv[1]));
      fHistCosPAAfter->Fill(pt, cosinePointingAngle);

      fHistDCADaughtersAfter->Fill(pt, v0->GetDcaV0Daughters());
      fHistDCA->Fill(pt, dcaPrim);
      fHistDecayLength->Fill(
          pt,
          std::sqrt(conv[0] * conv[0] + conv[1] * conv[1] + conv[2] * conv[2]));

      const float posPt = pos->Pt();
      const float negPt = neg->Pt();
      const float magField = inputEvent->GetMagneticField();
      const float dcaDaughterToPVPos =
          TMath::Abs(pos->GetD(vertex->GetX(), vertex->GetY(), magField));
      const short nClsTPCPos = pos->GetTPCNcls();
      const float nCrossedRowsPos = pos->GetTPCClusterInfo(2, 1);
      const short nFindablePos = pos->GetTPCNclsF();
      const float ratioFindablePos =
          nCrossedRowsPos / (static_cast<float>(nFindablePos) + 1.e-3);
      const int nClsSharedTPCPos = pos->GetTPCnclsS();
      const float chi2Pos =
          (nClsTPCPos > 5) ? (pos->GetTPCchi2() / float(nClsTPCPos - 5)) : -1.;

      const float dcaDaughterToPVNeg =
          TMath::Abs(neg->GetD(vertex->GetX(), vertex->GetY(), magField));
      const short nClsTPCNeg = neg->GetTPCNcls();
      const float nCrossedRowsNeg = neg->GetTPCClusterInfo(2, 1);
      const short nFindableNeg = neg->GetTPCNclsF();
      const float ratioFindableNeg =
          nCrossedRowsNeg / (static_cast<float>(nFindableNeg) + 1.e-3);
      const int nClsSharedTPCNeg = neg->GetTPCnclsS();
      const float chi2Neg =
          (nClsTPCNeg > 5) ? (neg->GetTPCchi2() / float(nClsTPCNeg - 5)) : -1.;

      int nClsSharedITSPos = 0;
      int nClsSharedITSNeg = 0;
      for (int i = 0; i < 6; ++i) {
        if (pos->HasSharedPointOnITSLayer(i)) ++nClsSharedITSPos;
        if (neg->HasSharedPointOnITSLayer(i)) ++nClsSharedITSNeg;
      }

      const float pidPos = fPIDResponse->NumberOfSigmasTPC(pos, fPosPID);
      const float pidNeg = fPIDResponse->NumberOfSigmasTPC(neg, fNegPID);

      bool posTrackITS =
          (pos->HasPointOnITSLayer(0) || pos->HasPointOnITSLayer(1) ||
           pos->HasPointOnITSLayer(4) || pos->HasPointOnITSLayer(5));
      bool negTrackITS =
          (neg->HasPointOnITSLayer(0) || pos->HasPointOnITSLayer(1) ||
           neg->HasPointOnITSLayer(4) || pos->HasPointOnITSLayer(5));
      bool posTrackTOF = pos->GetTOFBunchCrossing() == 0;
      bool negTrackTOF = neg->GetTOFBunchCrossing() == 0;

      bool posTrackCombined = (posTrackITS || posTrackTOF);
      bool negTrackCombined = (negTrackITS || negTrackTOF);

      if (!fIsLightweight) {
        if (posTrackITS) fHistSingleParticlePileUp[0]->Fill(0.f, posPt);
        if (posTrackTOF) fHistSingleParticlePileUp[0]->Fill(1, posPt);
        if (posTrackCombined) fHistSingleParticlePileUp[0]->Fill(2, posPt);
        if (!posTrackITS && !posTrackTOF)
          fHistSingleParticlePileUp[0]->Fill(3, posPt);
        if (negTrackITS) fHistSingleParticlePileUp[1]->Fill(0.f, negPt);
        if (negTrackTOF) fHistSingleParticlePileUp[1]->Fill(1, negPt);
        if (negTrackCombined) fHistSingleParticlePileUp[1]->Fill(2, negPt);
        if (!negTrackITS && !negTrackTOF)
          fHistSingleParticlePileUp[1]->Fill(3, negPt);
      }

      fHistSingleParticlePt[0]->Fill(posPt);
      fHistSingleParticleEtaAfter[0]->Fill(posPt, pos->Eta());
      fHistSingleParticleChi2After[0]->Fill(posPt, chi2Pos);
      fHistSingleParticleNclsTPCAfter[0]->Fill(posPt, nClsTPCPos);
      fHistSingleParticleNclsTPCFindableAfter[0]->Fill(posPt, nFindablePos);
      fHistSingleParticleNclsTPCRatioFindableAfter[0]->Fill(posPt,
                                                            ratioFindablePos);
      fHistSingleParticleNcrossedTPCAfter[0]->Fill(posPt, nCrossedRowsPos);
      fHistSingleParticleNclsTPCSharedAfter[0]->Fill(posPt, nClsSharedTPCPos);
      fHistSingleParticleNclsITSSharedAfter[0]->Fill(posPt, nClsSharedITSPos);
      fHistSingleParticleDCAtoPVAfter[0]->Fill(posPt, dcaDaughterToPVPos);
      fHistSingleParticlePID[0]->Fill(posPt, pidPos);

      fHistSingleParticlePt[1]->Fill(negPt);
      fHistSingleParticleEtaAfter[1]->Fill(negPt, neg->Eta());
      fHistSingleParticleChi2After[1]->Fill(negPt, chi2Neg);
      fHistSingleParticleNclsTPCAfter[1]->Fill(negPt, nClsTPCNeg);
      fHistSingleParticleNclsTPCFindableAfter[1]->Fill(negPt, nFindableNeg);
      fHistSingleParticleNclsTPCRatioFindableAfter[1]->Fill(negPt,
                                                            ratioFindableNeg);
      fHistSingleParticleNcrossedTPCAfter[1]->Fill(negPt, nCrossedRowsNeg);
      fHistSingleParticleNclsTPCSharedAfter[1]->Fill(negPt, nClsSharedTPCNeg);
      fHistSingleParticleNclsITSSharedAfter[1]->Fill(negPt, nClsSharedITSNeg);
      fHistSingleParticleDCAtoPVAfter[1]->Fill(negPt, dcaDaughterToPVNeg);
      fHistSingleParticlePID[1]->Fill(negPt, pidNeg);
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::CheckCutsMC() const {
  for (int iV0 = 0; iV0 < fInputEvent->GetNumberOfV0s(); ++iV0) {
    AliESDv0 *v0 = fInputEvent->GetV0(iV0);

    // On-fly status
    if (fV0OnFly) {
      // online v0
      if (!(v0->GetOnFlyStatus())) continue;
    } else {
      // offline v0
      if (v0->GetOnFlyStatus()) continue;
    }

    // Daughter tracks
    AliESDtrack *pos = fInputEvent->GetTrack(v0->GetPindex());
    AliESDtrack *neg = fInputEvent->GetTrack(v0->GetNindex());
    if (!pos || !neg) continue;

    if (pos->Charge() < 0) {
      pos = neg;
      neg = fInputEvent->GetTrack(v0->GetPindex());
    }
    AliSigma0ParticleV0 v0Candidate(v0, pos, neg,
                                    fInputEvent->GetPrimaryVertex(), fPID,
                                    fInputEvent->GetMagneticField(), fMCEvent);

    int label = v0Candidate.MatchToMC(fMCEvent, fPID, {{fPosPDG, fNegPDG}});

    const float pt = v0->Pt();
    const float posPt = pos->Pt();
    const float negPt = neg->Pt();
    const float magField = fInputEvent->GetMagneticField();
    const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
    const float dcaDaughterToPVPos =
        TMath::Abs(pos->GetD(vertex->GetX(), vertex->GetY(), magField));
    const short nClsTPCPos = pos->GetTPCNcls();
    const float nCrossedRowsPos = pos->GetTPCClusterInfo(2, 1);
    const short nFindablePos = pos->GetTPCNclsF();
    const float ratioFindablePos =
        nCrossedRowsPos / (static_cast<float>(nFindablePos) + 1.e-3);
    const int nClsSharedTPCPos = pos->GetTPCnclsS();
    const float chi2Pos =
        (nClsTPCPos > 5) ? (pos->GetTPCchi2() / float(nClsTPCPos - 5)) : -1.;

    const float dcaDaughterToPVNeg =
        TMath::Abs(neg->GetD(vertex->GetX(), vertex->GetY(), magField));
    const short nClsTPCNeg = neg->GetTPCNcls();
    const float nCrossedRowsNeg = neg->GetTPCClusterInfo(2, 1);
    const short nFindableNeg = neg->GetTPCNclsF();
    const float ratioFindableNeg =
        nCrossedRowsNeg / (static_cast<float>(nFindableNeg) + 1.e-3);
    const int nClsSharedTPCNeg = neg->GetTPCnclsS();
    const float chi2Neg =
        (nClsTPCNeg > 5) ? (neg->GetTPCchi2() / float(nClsTPCNeg - 5)) : -1.;

    int nClsSharedITSPos = 0;
    int nClsSharedITSNeg = 0;
    for (int i = 0; i < 6; ++i) {
      if (pos->HasSharedPointOnITSLayer(i)) ++nClsSharedITSPos;
      if (neg->HasSharedPointOnITSLayer(i)) ++nClsSharedITSNeg;
    }

    const float pidPos = fPIDResponse->NumberOfSigmasTPC(pos, fPosPID);
    const float pidNeg = fPIDResponse->NumberOfSigmasTPC(neg, fNegPID);

    Double_t xPV = vertex->GetX();
    Double_t yPV = vertex->GetY();
    Double_t zPV = vertex->GetZ();
    Double_t decayVertexV0[3];

    v0->GetXYZ(decayVertexV0[0], decayVertexV0[1], decayVertexV0[2]);

    // Calculate vertex variables:
    const float point = v0->GetV0CosineOfPointingAngle(xPV, yPV, zPV);
    const float dcaV0Dau = v0->GetDcaV0Daughters();
    const float transverseRadius =
        std::sqrt(decayVertexV0[0] * decayVertexV0[0] +
                  decayVertexV0[1] * decayVertexV0[1]);

    const float psiPair = ComputePsiPair(v0);

    float invMass = 0.f;
    if (fPID == 22) {
      invMass = ComputePhotonMassRefit(v0);
    } else {
      v0->ChangeMassHypothesis(fPID);
      invMass = v0->GetEffMass();
    }

    const float armAlpha = v0->AlphaV0();
    const float armQt = v0->PtArmV0();

    if (label > 0) {
      fHistV0MassPtTrue->Fill(pt, invMass);
      fHistDecayVertexXTrue->Fill(pt, decayVertexV0[0]);
      fHistDecayVertexYTrue->Fill(pt, decayVertexV0[1]);
      fHistDecayVertexZTrue->Fill(pt, decayVertexV0[2]);
      fHistTransverseRadiusTrue->Fill(pt, transverseRadius);
      fHistCosPATrue->Fill(pt, point);
      fHistDCADaughtersTrue->Fill(pt, dcaV0Dau);
      fHistArmenterosTrue->Fill(armAlpha, armQt);
      fHistPsiPairTrue->Fill(pt, psiPair);
      fHistSingleParticlePtTrue[0]->Fill(posPt);
      fHistSingleParticlePtTrue[1]->Fill(negPt);
      fHistSingleParticleDCAtoPVTrue[0]->Fill(posPt, dcaDaughterToPVPos);
      fHistSingleParticleDCAtoPVTrue[1]->Fill(negPt, dcaDaughterToPVNeg);

      fHistSingleParticleNclsTPCTrue[0]->Fill(posPt, nClsTPCPos);
      fHistSingleParticleNclsTPCTrue[1]->Fill(negPt, nClsTPCNeg);
      fHistSingleParticleNclsTPCFindableTrue[0]->Fill(posPt, nFindablePos);
      fHistSingleParticleNclsTPCFindableTrue[1]->Fill(negPt, nFindableNeg);
      fHistSingleParticleNclsTPCRatioFindableTrue[0]->Fill(posPt,
                                                           ratioFindablePos);
      fHistSingleParticleNclsTPCRatioFindableTrue[1]->Fill(negPt,
                                                           ratioFindableNeg);
      fHistSingleParticleNcrossedTPCTrue[0]->Fill(posPt, nCrossedRowsPos);
      fHistSingleParticleNcrossedTPCTrue[1]->Fill(negPt, nCrossedRowsNeg);
      fHistSingleParticleNclsTPCSharedTrue[0]->Fill(posPt, nClsSharedTPCPos);
      fHistSingleParticleNclsTPCSharedTrue[1]->Fill(negPt, nClsSharedTPCNeg);
      fHistSingleParticleNclsITSSharedTrue[0]->Fill(posPt, nClsSharedITSPos);
      fHistSingleParticleNclsITSSharedTrue[1]->Fill(negPt, nClsSharedITSNeg);
      fHistSingleParticleChi2True[0]->Fill(posPt, chi2Pos);
      fHistSingleParticleChi2True[1]->Fill(negPt, chi2Neg);
      fHistSingleParticlePIDTrue[0]->Fill(pos->P(), pidPos);
      fHistSingleParticlePIDTrue[1]->Fill(neg->P(), pidNeg);

      AliMCParticle *mcParticle =
          static_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
      if (!mcParticle) continue;
      const int pdgMother = static_cast<AliMCParticle *>(
                                fMCEvent->GetTrack(mcParticle->GetMother()))
                                ->PdgCode();
      fHistV0MotherTrue->Fill(pt, TMath::Abs(pdgMother));
      if (TMath::Abs(pdgMother) == 3212) {
        fHistV0MassPtTrueSigma->Fill(pt, invMass);
        fHistDecayVertexXTrueSigma->Fill(pt, decayVertexV0[0]);
        fHistDecayVertexYTrueSigma->Fill(pt, decayVertexV0[1]);
        fHistDecayVertexZTrueSigma->Fill(pt, decayVertexV0[2]);
        fHistTransverseRadiusTrueSigma->Fill(pt, transverseRadius);
        fHistCosPATrueSigma->Fill(pt, point);
        fHistDCADaughtersTrueSigma->Fill(pt, dcaV0Dau);
        fHistArmenterosTrueSigma->Fill(armAlpha, armQt);
        fHistPsiPairTrueSigma->Fill(pt, psiPair);
        fHistSingleParticlePtTrueSigma[0]->Fill(posPt);
        fHistSingleParticlePtTrueSigma[1]->Fill(negPt);
        fHistSingleParticleDCAtoPVTrueSigma[0]->Fill(posPt, dcaDaughterToPVPos);
        fHistSingleParticleDCAtoPVTrueSigma[1]->Fill(negPt, dcaDaughterToPVNeg);
        fHistSingleParticleNclsTPCTrueSigma[0]->Fill(posPt, nClsTPCPos);
        fHistSingleParticleNclsTPCTrueSigma[1]->Fill(negPt, nClsTPCNeg);
        fHistSingleParticleNclsTPCFindableTrueSigma[0]->Fill(posPt,
                                                             nFindablePos);
        fHistSingleParticleNclsTPCFindableTrueSigma[1]->Fill(negPt,
                                                             nFindableNeg);
        fHistSingleParticleNclsTPCRatioFindableTrueSigma[0]->Fill(
            posPt, ratioFindablePos);
        fHistSingleParticleNclsTPCRatioFindableTrueSigma[1]->Fill(
            negPt, ratioFindableNeg);
        fHistSingleParticleNcrossedTPCTrueSigma[0]->Fill(posPt,
                                                         nCrossedRowsPos);
        fHistSingleParticleNcrossedTPCTrueSigma[1]->Fill(negPt,
                                                         nCrossedRowsNeg);
        fHistSingleParticleNclsTPCSharedTrueSigma[0]->Fill(posPt,
                                                           nClsSharedTPCPos);
        fHistSingleParticleNclsTPCSharedTrueSigma[1]->Fill(negPt,
                                                           nClsSharedTPCNeg);
        fHistSingleParticleNclsITSSharedTrueSigma[0]->Fill(posPt,
                                                           nClsSharedITSPos);
        fHistSingleParticleNclsITSSharedTrueSigma[1]->Fill(negPt,
                                                           nClsSharedITSNeg);
        fHistSingleParticleChi2TrueSigma[0]->Fill(posPt, chi2Pos);
        fHistSingleParticleChi2TrueSigma[1]->Fill(negPt, chi2Neg);
        fHistSingleParticlePIDTrueSigma[0]->Fill(pos->P(), pidPos);
        fHistSingleParticlePIDTrueSigma[1]->Fill(neg->P(), pidNeg);
      }
    } else {
      fHistV0MassPtBkg->Fill(pt, invMass);
      fHistDecayVertexXBkg->Fill(pt, decayVertexV0[0]);
      fHistDecayVertexYBkg->Fill(pt, decayVertexV0[1]);
      fHistDecayVertexZBkg->Fill(pt, decayVertexV0[2]);
      fHistTransverseRadiusBkg->Fill(pt, transverseRadius);
      fHistCosPABkg->Fill(pt, point);
      fHistDCADaughtersBkg->Fill(pt, dcaV0Dau);
      fHistArmenterosBkg->Fill(armAlpha, armQt);
      fHistPsiPairBkg->Fill(pt, psiPair);
      fHistSingleParticlePtBkg[0]->Fill(posPt);
      fHistSingleParticlePtBkg[1]->Fill(negPt);
      fHistSingleParticleDCAtoPVBkg[0]->Fill(posPt, dcaDaughterToPVPos);
      fHistSingleParticleDCAtoPVBkg[1]->Fill(negPt, dcaDaughterToPVNeg);

      fHistSingleParticleNclsTPCBkg[0]->Fill(posPt, nClsTPCPos);
      fHistSingleParticleNclsTPCBkg[1]->Fill(negPt, nClsTPCNeg);
      fHistSingleParticleNclsTPCFindableBkg[0]->Fill(posPt, nFindablePos);
      fHistSingleParticleNclsTPCFindableBkg[1]->Fill(negPt, nFindableNeg);
      fHistSingleParticleNclsTPCRatioFindableBkg[0]->Fill(posPt,
                                                          ratioFindablePos);
      fHistSingleParticleNclsTPCRatioFindableBkg[1]->Fill(negPt,
                                                          ratioFindableNeg);
      fHistSingleParticleNcrossedTPCBkg[0]->Fill(posPt, nCrossedRowsPos);
      fHistSingleParticleNcrossedTPCBkg[1]->Fill(negPt, nCrossedRowsNeg);
      fHistSingleParticleNclsTPCSharedBkg[0]->Fill(posPt, nClsSharedTPCPos);
      fHistSingleParticleNclsTPCSharedBkg[1]->Fill(negPt, nClsSharedTPCNeg);
      fHistSingleParticleNclsITSSharedBkg[0]->Fill(posPt, nClsSharedITSPos);
      fHistSingleParticleNclsITSSharedBkg[1]->Fill(negPt, nClsSharedITSNeg);
      fHistSingleParticleChi2Bkg[0]->Fill(posPt, chi2Pos);
      fHistSingleParticleChi2Bkg[1]->Fill(negPt, chi2Neg);
      fHistSingleParticlePIDBkg[0]->Fill(pos->P(), pidPos);
      fHistSingleParticlePIDBkg[1]->Fill(neg->P(), pidNeg);
    }
  }
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::CheckDaughters(const AliMCParticle *particle) const {
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
    if (pdgCode == fPosPDG)
      posDaughter = tmpDaughter;
    else if (pdgCode == fNegPDG)
      negDaughter = tmpDaughter;
  }

  if (!posDaughter || !negDaughter) return false;

  if (particle->PdgCode() == 22 &&
      (posDaughter->Particle()->GetUniqueID() != 5 ||
       negDaughter->Particle()->GetUniqueID() != 5)) {
    return false;
  }
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::CheckDaughtersInAcceptance(
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
    if (pdgCode == fPosPDG)
      posDaughter = tmpDaughter;
    else if (pdgCode == fNegPDG)
      negDaughter = tmpDaughter;
  }

  if (!posDaughter || !negDaughter) return false;

  if (particle->PdgCode() == 22 &&
      (posDaughter->Particle()->GetUniqueID() != 5 ||
       negDaughter->Particle()->GetUniqueID() != 5)) {
    return false;
  }

  // Daughter Eta
  const float posEta = posDaughter->Eta();
  const float negEta = negDaughter->Eta();
  if (posEta > fEtaMax) return false;
  if (negEta > fEtaMax) return false;

  // Mother pT
  const float pT = particle->Pt();
  if (pT > fV0PtMax || pT < fV0PtMin) return false;

  // Mother decay vertex
  Double_t posVertex[3];
  posDaughter->XvYvZv(posVertex);
  const float rPos =
      std::sqrt(posVertex[0] * posVertex[0] + posVertex[1] * posVertex[1]);
  if (rPos > fV0RadiusMax || rPos < fV0RadiusMin) return false;
  if (TMath::Abs(posVertex[0]) > fV0DecayVertexMax) return false;
  if (TMath::Abs(posVertex[1]) > fV0DecayVertexMax) return false;
  if (TMath::Abs(posVertex[2]) > fV0DecayVertexMax) return false;

  Double_t negVertex[3];
  negDaughter->XvYvZv(negVertex);
  const float rNeg =
      std::sqrt(negVertex[0] * negVertex[0] + negVertex[1] * negVertex[1]);
  if (rNeg > fV0RadiusMax || rNeg < fV0RadiusMin) return false;
  if (TMath::Abs(negVertex[0]) > fV0DecayVertexMax) return false;
  if (TMath::Abs(negVertex[1]) > fV0DecayVertexMax) return false;
  if (TMath::Abs(negVertex[2]) > fV0DecayVertexMax) return false;

  return true;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::InitCutHistograms(TString appendix) {
  std::cout << "============================\n"
            << " V0 CUT CONFIGURATION " << appendix << "\n"
            << " On-fly       " << fV0OnFly << "\n"
            << " p_T min      " << fV0PtMin << "\n"
            << " p_T max      " << fV0PtMax << "\n"
            << " CPA min      " << fV0CosPAMin << "\n"
            << " r max        " << fV0RadiusMax << "\n"
            << " r min        " << fV0RadiusMin << "\n"
            << " decay vtx    " << fV0DecayVertexMax << "\n"
            << " PID nsigma   " << fPIDnSigma << "\n"
            << " nCls min     " << fTPCclusterMin << "\n"
            << " Crossed rows " << fTPCnCrossedRowsMin << "\n"
            << " Ratio find.  " << fTPCratioFindable << "\n"
            << " nFinable     " << fTPCfindableMin << "\n"
            << " nSharedMax   " << fTPCnSharedMax << "\n"
            << " eta max      " << fEtaMax << "\n"
            << " pile up rej  " << static_cast<int>(fPileUpRejectionMode)
            << "\n"
            << " daughter dca " << fDaughterDCAMax << "\n"
            << " daugher pv   " << fDaughterDCAPV << "\n"
            << " K0 rej       " << fK0RejectionLow << " " << fK0RejectionUp
            << "\n"
            << " Lambda sel   " << fLambdaSelectionLow << " "
            << fLambdaSelectionUp << "\n"
            << " Psi_pair max " << fPsiPairMax << "\n"
            << " PID pos      " << fPosPID << "\n"
            << " PDG pos      " << fPosPDG << "\n"
            << " PID pos      " << fNegPID << "\n"
            << " PDG pos      " << fNegPDG << "\n"
            << "============================\n";

  const float pi = TMath::Pi();

  TH1::AddDirectory(kFALSE);
  TString name;
  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    name = "V0_" + appendix;
    fHistograms->SetName(name);
  }

  if (!fIsLightweight) {
    if (fHistogramsPos != nullptr) {
      delete fHistogramsPos;
      fHistogramsPos = nullptr;
    }
    if (fHistogramsPos == nullptr) {
      fHistogramsPos = new TList();
      fHistogramsPos->SetOwner(kTRUE);
      fHistogramsPos->SetName("V0_PosDaughter");
    }

    if (fHistogramsNeg != nullptr) {
      delete fHistogramsNeg;
      fHistogramsNeg = nullptr;
    }
    if (fHistogramsNeg == nullptr) {
      fHistogramsNeg = new TList();
      fHistogramsNeg->SetOwner(kTRUE);
      fHistogramsNeg->SetName("V0_NegDaughter");
    }
  }

  fHistCutBooking = new TProfile("fHistCutBooking", ";;Cut value", 32, 0, 32);
  fHistCutBooking->GetXaxis()->SetBinLabel(1, "V0 on fly");
  fHistCutBooking->GetXaxis()->SetBinLabel(2, "#it{p}_{T} min");
  fHistCutBooking->GetXaxis()->SetBinLabel(3, "#it{p}_{T} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(4, "cos#alpha min");
  fHistCutBooking->GetXaxis()->SetBinLabel(5, "max r");
  fHistCutBooking->GetXaxis()->SetBinLabel(6, "min r");
  fHistCutBooking->GetXaxis()->SetBinLabel(7, "Decay vtx max");
  fHistCutBooking->GetXaxis()->SetBinLabel(8, "n#sigma PID");
  fHistCutBooking->GetXaxis()->SetBinLabel(9, "TPC nCls min");
  fHistCutBooking->GetXaxis()->SetBinLabel(10, "TPC nCrossed min");
  fHistCutBooking->GetXaxis()->SetBinLabel(11, "TPC Ratio Findable");
  fHistCutBooking->GetXaxis()->SetBinLabel(12, "TPC nCls Findable min");
  fHistCutBooking->GetXaxis()->SetBinLabel(13, "TPC nShared max");
  fHistCutBooking->GetXaxis()->SetBinLabel(14, "#eta max");
  fHistCutBooking->GetXaxis()->SetBinLabel(15, "#chi2 max");
  fHistCutBooking->GetXaxis()->SetBinLabel(16, "Daughter DCA max");
  fHistCutBooking->GetXaxis()->SetBinLabel(17, "Daughter DCA pvtx min");
  fHistCutBooking->GetXaxis()->SetBinLabel(18, "K^{0} rejection low");
  fHistCutBooking->GetXaxis()->SetBinLabel(19, "K^{0} rejection up");
  fHistCutBooking->GetXaxis()->SetBinLabel(20, "#Lambda selection low");
  fHistCutBooking->GetXaxis()->SetBinLabel(21, "#Lambda selection up");
  fHistCutBooking->GetXaxis()->SetBinLabel(22, "Pile-up rejection");
  fHistCutBooking->GetXaxis()->SetBinLabel(23, "Armenteros q_{T} low");
  fHistCutBooking->GetXaxis()->SetBinLabel(24, "Armenteros q_{T} up");
  fHistCutBooking->GetXaxis()->SetBinLabel(25, "Armenteros #alpha low");
  fHistCutBooking->GetXaxis()->SetBinLabel(26, "Armenteros #alpha up");
  fHistCutBooking->GetXaxis()->SetBinLabel(27, "#Psi_{pair} max");
  fHistCutBooking->GetXaxis()->SetBinLabel(28, "PID pos");
  fHistCutBooking->GetXaxis()->SetBinLabel(29, "PDG pos");
  fHistCutBooking->GetXaxis()->SetBinLabel(20, "PID neg");
  fHistCutBooking->GetXaxis()->SetBinLabel(31, "PDG neg");
  fHistCutBooking->GetXaxis()->SetBinLabel(32, "MC Mult for efficiency");
  fHistograms->Add(fHistCutBooking);

  fHistCutBooking->Fill(0.f, static_cast<double>(fV0OnFly));
  fHistCutBooking->Fill(1.f, fV0PtMin);
  fHistCutBooking->Fill(2.f, fV0PtMax);
  fHistCutBooking->Fill(3.f, fV0CosPAMin);
  fHistCutBooking->Fill(4.f, fV0RadiusMax);
  fHistCutBooking->Fill(5.f, fV0RadiusMin);
  fHistCutBooking->Fill(6.f, fV0DecayVertexMax);
  fHistCutBooking->Fill(7.f, fPIDnSigma);
  fHistCutBooking->Fill(8.f, fTPCclusterMin);
  fHistCutBooking->Fill(9.f, fTPCnCrossedRowsMin);
  fHistCutBooking->Fill(10.f, fTPCratioFindable);
  fHistCutBooking->Fill(11.f, fTPCfindableMin);
  fHistCutBooking->Fill(12.f, fTPCnSharedMax);
  fHistCutBooking->Fill(13.f, fEtaMax);
  fHistCutBooking->Fill(14.f, fChi2Max);
  fHistCutBooking->Fill(15.f, fDaughterDCAMax);
  fHistCutBooking->Fill(16.f, fDaughterDCAPV);
  fHistCutBooking->Fill(17.f, fK0RejectionLow);
  fHistCutBooking->Fill(18.f, fK0RejectionUp);
  fHistCutBooking->Fill(19.f, fLambdaSelectionLow);
  fHistCutBooking->Fill(20.f, fLambdaSelectionUp);
  fHistCutBooking->Fill(21.f, static_cast<double>(fPileUpRejectionMode));
  fHistCutBooking->Fill(22.f, fArmenterosQtLow);
  fHistCutBooking->Fill(23.f, fArmenterosQtUp);
  fHistCutBooking->Fill(24.f, fArmenterosAlphaLow);
  fHistCutBooking->Fill(25.f, fArmenterosAlphaUp);
  fHistCutBooking->Fill(26.f, fPsiPairMax);
  fHistCutBooking->Fill(27.f, fPosPID);
  fHistCutBooking->Fill(28.f, fPosPDG);
  fHistCutBooking->Fill(29.f, fNegPID);
  fHistCutBooking->Fill(30.f, fNegPDG);
  fHistCutBooking->Fill(31.f, fMCHighMultThreshold);

  if (TMath::Abs(fPID) == 3122) {
    fHistV0MassPt =
        new TH2F("InvMassPt",
                 "; #it{p}_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^{2})",
                 100, 0, 10, 100, 1.07, 1.17);
  } else if (TMath::Abs(fPID) == 22) {
    fHistV0MassPt =
        new TH2F("InvMassPt",
                 "; #it{p}_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^{2})",
                 100, 0, 10, 100, 0., .1);
  }
  fHistograms->Add(fHistV0MassPt);

  if (!fIsLightweight) {
    fHistCuts = new TH1F("fHistCuts", ";;Entries", 20, 0, 20);
    fHistCuts->GetXaxis()->SetBinLabel(1, "V0");
    fHistCuts->GetXaxis()->SetBinLabel(2, "Neutral");
    fHistCuts->GetXaxis()->SetBinLabel(3, "On-fly");
    fHistCuts->GetXaxis()->SetBinLabel(4, "#it{p}_{T, min}");
    fHistCuts->GetXaxis()->SetBinLabel(5, "#it{p}_{T, max}");
    fHistCuts->GetXaxis()->SetBinLabel(6, "Daughter tracks");
    fHistCuts->GetXaxis()->SetBinLabel(7, "SingleParticle QA");
    fHistCuts->GetXaxis()->SetBinLabel(8, "Pile-up");
    fHistCuts->GetXaxis()->SetBinLabel(9, "Decay vtx x");
    fHistCuts->GetXaxis()->SetBinLabel(10, "Decay vtx y");
    fHistCuts->GetXaxis()->SetBinLabel(11, "Decay vtx z");
    fHistCuts->GetXaxis()->SetBinLabel(12, "Radius min");
    fHistCuts->GetXaxis()->SetBinLabel(13, "Radius max");
    fHistCuts->GetXaxis()->SetBinLabel(14, "Daughter DCA");
    fHistCuts->GetXaxis()->SetBinLabel(15, "cos #alpha");
    fHistCuts->GetXaxis()->SetBinLabel(16, "PID");
    if (TMath::Abs(fPID) == 3122) {
      fHistCuts->GetXaxis()->SetBinLabel(17, "K^{0} rejection");
      fHistCuts->GetXaxis()->SetBinLabel(18, "#Lambda selection");
    } else if (TMath::Abs(fPID) == 22) {
      fHistCuts->GetXaxis()->SetBinLabel(17, "#Psi_{pair} selection");
      fHistCuts->GetXaxis()->SetBinLabel(18, "K^{0} rejection");
      fHistCuts->GetXaxis()->SetBinLabel(19, "#Lambda rejection");
    }
    fHistograms->Add(fHistCuts);

    fHistNV0 =
        new TH1F("fHistNV0", ";Number of V0 candidates; Entries", 15, 0, 15);
    fHistograms->Add(fHistNV0);

    fHistLambdaMass = new TH1F(
        "fHistLambdaMass",
        "; Invariant mass p#pi^{-} hypothesis (GeV/#it{c}^{2}); Entries", 125,
        1., 1.25);
    fHistograms->Add(fHistLambdaMass);

    fHistAntiLambdaMass = new TH1F(
        "fHistAntiLambdaMass",
        "; Invariant mass #bar{p}#pi hypothesis (GeV/#it{c}^{2}); Entries", 125,
        1., 1.25);
    fHistograms->Add(fHistAntiLambdaMass);

    fHistPhotonMass = new TH1F(
        "fHistPhotonMass",
        "; Invariant mass e^{+}e^{-} hypothesis (GeV/#it{c}^{2}); Entries", 125,
        0., 0.25);
    fHistograms->Add(fHistPhotonMass);

    fHistPhotonMassRefit = new TH1F("fHistPhotonMassRefit",
                                    "; Invariant mass e^{+}e^{-} refit "
                                    "hypothesis (GeV/#it{c}^{2}); Entries",
                                    125, 0., 0.05);
    fHistograms->Add(fHistPhotonMassRefit);

    fHistK0Mass =
        new TH1F("fHistK0Mass",
                 "; Invariant mass #pi#pi hypothesis (GeV/#it{c}^{2}); Entries",
                 125, 0.4, 0.65);
    fHistograms->Add(fHistK0Mass);

    fHistV0Pt =
        new TH1F("fHistV0Pt", "; #it{p}_{T} (GeV/#it{c}); Entries", 100, 0, 10);
    fHistograms->Add(fHistV0Pt);

    if (TMath::Abs(fPID) == 3122) {
      fHistV0Mass =
          new TH1F("fHistV0Mass", "; Invariant mass (GeV/#it{c}^{2}); Entries",
                   200, 1., 1.2);
    } else if (TMath::Abs(fPID) == 22) {
      fHistV0Mass =
          new TH1F("fHistV0Mass", "; Invariant mass (GeV/#it{c}^{2}); Entries",
                   200, 0., 0.2);
    }
    fHistograms->Add(fHistV0Mass);

    fHistK0MassAfter =
        new TH1F("fHistK0MassAfter",
                 "; Invariant mass #pi#pi hypothesis (GeV/#it{c}^{2}); Entries",
                 125, 0.4, 0.65);
    fHistograms->Add(fHistK0MassAfter);
    fHistCosPA =
        new TH2F("fHistCosPA", "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)", 100, 0,
                 10, 100, 0.95, 1);
    fHistograms->Add(fHistCosPA);

    fHistEtaPhi =
        new TH2F("fHistEtaPhi", "; #eta; #phi", 100, -1, 1, 100, 0, 2 * pi);
    fHistograms->Add(fHistEtaPhi);

    if (fPID == 22) {
      fHistPsiPair =
          new TH2F("fHistPsiPair", "; #it{p}_{T} (GeV/#it{c}); #Psi_{pair}",
                   100, 0, 10, 200, -pi / 2.f, pi / 2.f);
      fHistograms->Add(fHistPsiPair);
    }
  }

  if (!fIsLightweight) {
    fHistSingleParticleCuts[0] =
        new TH1F("fHistSingleParticleCuts_pos", ";;Entries", 11, 0, 11);
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(2, "#eta");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(3, "#chi^{2}");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(4, "nCls TPC");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(5, "nCrossed Rows TPC");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(6,
                                                        "Ratio Findable TPC");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(7, "nCls Findable");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(8, "nCls Shared max");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(9,
                                                        "Daughter DCA to PV");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(10, "TPC refit");
    fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(11, "Kink");
    fHistogramsPos->Add(fHistSingleParticleCuts[0]);

    fHistSingleParticleCuts[1] =
        new TH1F("fHistSingleParticleCuts_neg", ";;Entries", 11, 0, 11);
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(2, "#eta");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(3, "#chi^{2}");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(4, "nCls TPC");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(5, "nCrossed Rows TPC");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(6,
                                                        "Ratio Findable TPC");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(7, "nCls Findable");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(8, "nCls Shared max");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(9,
                                                        "Daughter DCA to PV");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(10, "TPC refit");
    fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(11, "Kink");
    fHistogramsNeg->Add(fHistSingleParticleCuts[1]);

    if (fHistogramsBefore != nullptr) {
      delete fHistogramsBefore;
      fHistogramsBefore = nullptr;
    }
    if (fHistogramsBefore == nullptr) {
      fHistogramsBefore = new TList();
      fHistogramsBefore->SetOwner(kTRUE);
      fHistogramsBefore->SetName("Before");
    }

    if (fHistogramsAfter != nullptr) {
      delete fHistogramsAfter;
      fHistogramsAfter = nullptr;
    }
    if (fHistogramsAfter == nullptr) {
      fHistogramsAfter = new TList();
      fHistogramsAfter->SetOwner(kTRUE);
      fHistogramsAfter->SetName("After");
    }

    fHistDecayVertexXBefore =
        new TH2F("fHistDecayVertexXBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex x (cm)", 50, 0, 10,
                 200, 0, 200);
    fHistogramsBefore->Add(fHistDecayVertexXBefore);

    fHistDecayVertexYBefore =
        new TH2F("fHistDecayVertexYBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex y (cm)", 50, 0, 10,
                 200, 0, 200);
    fHistogramsBefore->Add(fHistDecayVertexYBefore);

    fHistDecayVertexZBefore =
        new TH2F("fHistDecayVertexZBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex z (cm)", 50, 0, 10,
                 200, 0, 200);
    fHistogramsBefore->Add(fHistDecayVertexZBefore);

    fHistDecayVertexXAfter =
        new TH2F("fHistDecayVertexXAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex x (cm)", 50, 0, 10,
                 200, 0, 200);
    fHistogramsAfter->Add(fHistDecayVertexXAfter);

    fHistDecayVertexYAfter =
        new TH2F("fHistDecayVertexYAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex y (cm)", 50, 0, 10,
                 200, 0, 200);
    fHistogramsAfter->Add(fHistDecayVertexYAfter);

    fHistDecayVertexZAfter =
        new TH2F("fHistDecayVertexZAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex z (cm)", 50, 0, 10,
                 200, 0, 200);
    fHistogramsAfter->Add(fHistDecayVertexZAfter);

    fHistTransverseRadiusBefore =
        new TH2F("fHistTransverseRadiusBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Transverse radius (cm)", 50, 0, 10,
                 200, 0, 200);
    fHistogramsBefore->Add(fHistTransverseRadiusBefore);

    fHistTransverseRadiusAfter =
        new TH2F("fHistTransverseRadiusAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Transverse radius (cm)", 50, 0, 10,
                 200, 0, 200);
    fHistogramsAfter->Add(fHistTransverseRadiusAfter);

    fHistCosPABefore =
        new TH2F("fHistCosPABefore", "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)",
                 100, 0, 10, 100, 0.95, 1);
    fHistogramsBefore->Add(fHistCosPABefore);

    fHistCosPAAfter =
        new TH2F("fHistCosPAAfter", "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)",
                 50, 0, 10, 100, 0.95, 1);
    fHistogramsAfter->Add(fHistCosPAAfter);

    fHistDCADaughtersBefore =
        new TH2F("fHistDCADaughtersBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Daughter DCA at decay vertex (cm)",
                 50, 0, 10, 100, 0, 2);
    fHistogramsBefore->Add(fHistDCADaughtersBefore);

    fHistDCADaughtersAfter =
        new TH2F("fHistDCADaughtersAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Daughter DCA at decay vertex (cm)",
                 50, 0, 10, 100, 0, 2);
    fHistogramsAfter->Add(fHistDCADaughtersAfter);

    fHistDCA = new TH2F("fHistDCA", "; #it{p}_{T} (GeV/#it{c}); DCA to PV (cm)",
                        50, 0, 10, 100, 0, 10);
    fHistogramsAfter->Add(fHistDCA);

    fHistDecayLength = new TH2F("fHistDecayLength",
                                "; #it{p}_{T} (GeV/#it{c}); Decay length (cm)",
                                50, 0, 10, 200, 0, 200);
    fHistogramsAfter->Add(fHistDecayLength);

    fHistArmenterosBefore =
        new TH2F("fHistArmenterosBefore", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 200, -1, 1, 100, 0, 0.25);
    fHistogramsBefore->Add(fHistArmenterosBefore);

    fHistArmenterosAfter =
        new TH2F("fHistArmenterosAfter", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 200, -1, 1, 100, 0, 0.25);
    fHistogramsAfter->Add(fHistArmenterosAfter);

    fHistograms->Add(fHistogramsBefore);
    fHistograms->Add(fHistogramsAfter);

    fHistSingleParticlePt[0] = new TH1F("fHistSingleParticlePt_pos",
                                        ";#it{p}_{T}; Entries", 100, 0, 10);
    fHistSingleParticlePt[1] = new TH1F("fHistSingleParticlePt_neg",
                                        ";#it{p}_{T}; Entries", 100, 0, 10);
    fHistSingleParticleEtaBefore[0] =
        new TH2F("fHistSingleParticleEtaBefore_pos",
                 "; #it{p}_{T} (GeV/#it{c}); #eta", 50, 0, 10, 100, -2, 2);
    fHistSingleParticleEtaBefore[1] =
        new TH2F("fHistSingleParticleEtaBefore_neg",
                 "; #it{p}_{T} (GeV/#it{c}); #eta", 50, 0, 10, 100, -2, 2);
    fHistSingleParticleEtaAfter[0] =
        new TH2F("fHistSingleParticleEtaAfter_pos",
                 "; #it{p}_{T} (GeV/#it{c}); #eta", 50, 0, 10, 100, -2, 2);
    fHistSingleParticleEtaAfter[1] =
        new TH2F("fHistSingleParticleEtaAfter_neg",
                 "; #it{p}_{T} (GeV/#it{c}); #eta", 50, 0, 10, 100, -2, 2);
    fHistSingleParticleChi2Before[0] =
        new TH2F("fHistSingleParticleChi2Before_pos",
                 "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 50, 0, 10, 100, 0, 20);
    fHistSingleParticleChi2Before[1] =
        new TH2F("fHistSingleParticleChi2Before_neg",
                 "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 50, 0, 10, 100, 0, 20);
    fHistSingleParticleChi2After[0] =
        new TH2F("fHistSingleParticleChi2After_pos",
                 "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 50, 0, 10, 100, 0, 20);
    fHistSingleParticleChi2After[1] =
        new TH2F("fHistSingleParticleChi2After_neg",
                 "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 50, 0, 10, 100, 0, 20);
    fHistSingleParticleNclsTPCBefore[0] = new TH2F(
        "fHistSingleParticleNclsTPCBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCBefore[1] = new TH2F(
        "fHistSingleParticleNclsTPCBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCAfter[0] = new TH2F(
        "fHistSingleParticleNclsTPCAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCAfter[1] = new TH2F(
        "fHistSingleParticleNclsTPCAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCFindableBefore[0] =
        new TH2F("fHistSingleParticleNclsTPCFindableBefore_pos",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 50, 0, 10,
                 160, 0, 160);
    fHistSingleParticleNclsTPCFindableBefore[1] =
        new TH2F("fHistSingleParticleNclsTPCFindableBefore_neg",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 50, 0, 10,
                 160, 0, 160);
    fHistSingleParticleNclsTPCFindableAfter[0] =
        new TH2F("fHistSingleParticleNclsTPCFindableAfter_pos",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 50, 0, 10,
                 160, 0, 160);
    fHistSingleParticleNclsTPCFindableAfter[1] =
        new TH2F("fHistSingleParticleNclsTPCFindableAfter_neg",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 50, 0, 10,
                 160, 0, 160);
    fHistSingleParticleNclsTPCRatioFindableBefore[0] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 50, 0, 10, 200, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableBefore[1] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 50, 0, 10, 200, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableAfter[0] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 50, 0, 10, 200, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableAfter[1] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 50, 0, 10, 200, 0, 2);
    fHistSingleParticleNcrossedTPCBefore[0] = new TH2F(
        "fHistSingleParticleNcrossedTPCBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNcrossedTPCBefore[1] = new TH2F(
        "fHistSingleParticleNcrossedTPCBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNcrossedTPCAfter[0] = new TH2F(
        "fHistSingleParticleNcrossedTPCAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNcrossedTPCAfter[1] = new TH2F(
        "fHistSingleParticleNcrossedTPCAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCSharedBefore[0] = new TH2F(
        "fHistSingleParticleNclsTPCSharedBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCSharedBefore[1] = new TH2F(
        "fHistSingleParticleNclsTPCSharedBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsITSSharedBefore[0] = new TH2F(
        "fHistSingleParticleNclsITSSharedBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared", 50, 0, 10, 6, 0, 6);
    fHistSingleParticleNclsITSSharedBefore[1] = new TH2F(
        "fHistSingleParticleNclsITSSharedBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared", 50, 0, 10, 6, 0, 6);
    fHistSingleParticleNclsTPCSharedAfter[0] = new TH2F(
        "fHistSingleParticleNclsTPCSharedAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsTPCSharedAfter[1] = new TH2F(
        "fHistSingleParticleNclsTPCSharedAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared", 50, 0, 10, 160, 0, 160);
    fHistSingleParticleNclsITSSharedAfter[0] = new TH2F(
        "fHistSingleParticleNclsITSSharedAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared", 50, 0, 10, 6, 0, 6);
    fHistSingleParticleNclsITSSharedAfter[1] = new TH2F(
        "fHistSingleParticleNclsITSSharedAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared", 50, 0, 10, 6, 0, 6);
    fHistSingleParticleDCAtoPVBefore[0] = new TH2F(
        "fHistSingleParticleDCAtoPVBefore_pos",
        "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA)", 50, 0, 10, 100, 0, 10);
    fHistSingleParticleDCAtoPVBefore[1] = new TH2F(
        "fHistSingleParticleDCAtoPVBefore_neg",
        "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA)", 50, 0, 10, 100, 0, 10);
    fHistSingleParticleDCAtoPVAfter[0] = new TH2F(
        "fHistSingleParticleDCAtoPVAfter_pos",
        "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA)", 50, 0, 10, 100, 0, 10);
    fHistSingleParticleDCAtoPVAfter[1] = new TH2F(
        "fHistSingleParticleDCAtoPVAfter_neg",
        "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA)", 50, 0, 10, 100, 0, 10);
    fHistSingleParticlePileUp[0] =
        new TH2F("fHistSingleParticlePileUp_pos",
                 "; Pileup flag; #it{p}_{T} (GeV/#it{c})", 4, 0, 4, 50, 0, 10);
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(1, "ITS");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(2, "TOF");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(3, "Combined");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(4, "None");
    fHistSingleParticlePileUp[1] =
        new TH2F("fHistSingleParticlePileUp_neg",
                 "; Pileup flag; #it{p}_{T} (GeV/#it{c})", 4, 0, 4, 50, 0, 10);
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(1, "ITS");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(2, "TOF");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(3, "Combined");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(4, "None");

    fHistSingleParticlePID[0] = new TH2F(
        "fHistSingleParticlePID_pos", "; #it{p} (GeV/#it{c}); n_{#sigma} TPC",
        50, 0, 10, 250, -10, 10);
    fHistSingleParticlePID[1] = new TH2F(
        "fHistSingleParticlePID_neg", "; #it{p} (GeV/#it{c}); n_{#sigma} TPC",
        50, 0, 10, 250, -10, 10);

    fHistogramsPos->Add(fHistSingleParticlePt[0]);
    fHistogramsNeg->Add(fHistSingleParticlePt[1]);
    fHistogramsPos->Add(fHistSingleParticleEtaBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleEtaBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleChi2Before[0]);
    fHistogramsNeg->Add(fHistSingleParticleChi2Before[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCFindableBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCFindableBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCRatioFindableBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCRatioFindableBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNcrossedTPCBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNcrossedTPCBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCSharedBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCSharedBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsITSSharedBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsITSSharedBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleDCAtoPVBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleDCAtoPVBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleEtaAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleEtaAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleChi2After[0]);
    fHistogramsNeg->Add(fHistSingleParticleChi2After[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCFindableAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCFindableAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCRatioFindableAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCRatioFindableAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNcrossedTPCAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNcrossedTPCAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCSharedAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCSharedAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsITSSharedAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsITSSharedAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleDCAtoPVAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleDCAtoPVAfter[1]);
    fHistogramsPos->Add(fHistSingleParticlePileUp[0]);
    fHistogramsNeg->Add(fHistSingleParticlePileUp[1]);
    fHistogramsPos->Add(fHistSingleParticlePID[0]);
    fHistogramsNeg->Add(fHistSingleParticlePID[1]);
  }
  fHistograms->Add(fHistogramsPos);
  fHistograms->Add(fHistogramsNeg);

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

    fHistMCTruthV0PtY =
        new TH2F("fHistMCTruthV0PtY", "; y; #it{p}_{T} (GeV/#it{c})", 100, -10,
                 10, 100, 0, 10);
    fHistMCTruthV0DaughterPtY =
        new TH2F("fHistMCTruthV0DaughterPtY", "; y; #it{p}_{T} (GeV/#it{c})",
                 100, -10, 10, 100, 0, 10);
    fHistMCTruthV0DaughterPtYAccept =
        new TH2F("fHistMCTruthV0DaughterPtYAccept",
                 "; y; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);

    fHistMCTruthPtYHighMult =
        new TH2F("fHistMCTruthPtYHighMult", "; y; #it{p}_{T} (GeV/#it{c})", 100,
                 -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtYHighMult =
        new TH2F("fHistMCTruthDaughterPtYHighMult",
                 "; y; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);
    fHistMCTruthDaughterPtYAcceptHighMult =
        new TH2F("fHistMCTruthDaughterPtYAcceptHighMult",
                 "; y; #it{p}_{T} (GeV/#it{c})", 100, -10, 10, 100, 0, 10);

    fHistMCV0Pt = new TH1F("fHistMCV0Pt", "; #it{p}_{T} (GeV/#it{c}); Entries",
                           100, 0, 10);

    fHistV0Mother =
        new TH2F("fHistV0Mother", "; #it{p}_{T} (GeV/#it{c}); PDG code mother",
                 100, 0, 10, 4000, 0, 4000);

    fHistogramsMC->Add(fHistMCTruthV0PtY);
    fHistogramsMC->Add(fHistMCTruthV0DaughterPtY);
    fHistogramsMC->Add(fHistMCTruthV0DaughterPtYAccept);
    fHistogramsMC->Add(fHistMCTruthPtYHighMult);
    fHistogramsMC->Add(fHistMCTruthDaughterPtYHighMult);
    fHistogramsMC->Add(fHistMCTruthDaughterPtYAcceptHighMult);
    fHistogramsMC->Add(fHistMCV0Pt);
    fHistogramsMC->Add(fHistV0Mother);

    if (fCheckCutsMC) {
      // TRUE V0
      fHistV0MotherTrue = new TH2F("fHistV0MotherTrue",
                                   "; #it{p}_{T} (GeV/#it{c}); PDG code mother",
                                   100, 0, 10, 4000, 0, 4000);
      fHistogramsMC->Add(fHistV0MotherTrue);

      fHistV0MassPtTrue =
          new TH2F("fHistV0MassPtTrue",
                   "; #it{p}_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^{2})",
                   100, 0, 10, 1200, 0., 1.2);
      fHistogramsMC->Add(fHistV0MassPtTrue);

      fHistDecayVertexXTrue =
          new TH2F("fHistDecayVertexXTrue",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex x true (cm)", 100,
                   0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexXTrue);

      fHistDecayVertexYTrue =
          new TH2F("fHistDecayVertexYTrue",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex y true (cm)", 100,
                   0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexYTrue);

      fHistDecayVertexZTrue =
          new TH2F("fHistDecayVertexZTrue",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex z true (cm)", 100,
                   0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexZTrue);

      fHistTransverseRadiusTrue =
          new TH2F("fHistTransverseRadiusTrue",
                   "; #it{p}_{T} (GeV/#it{c}); Transverse radius true (cm)",
                   100, 0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistTransverseRadiusTrue);

      fHistCosPATrue = new TH2F("fHistCosPATrue",
                                "; #it{p}_{T} (GeV/#it{c}); cos(#alpha) true",
                                100, 0, 10, 500, 0.8, 1);
      fHistogramsMC->Add(fHistCosPATrue);

      fHistDCADaughtersTrue = new TH2F(
          "fHistDCADaughtersTrue",
          "; #it{p}_{T} (GeV/#it{c}); Daughter DCA at decay vertex true (cm)",
          100, 0, 10, 500, 0, 2);
      fHistogramsMC->Add(fHistDCADaughtersTrue);

      fHistArmenterosTrue =
          new TH2F("fHistArmenterosTrue", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                   500, -1, 1, 500, 0, 1);
      fHistogramsMC->Add(fHistArmenterosTrue);

      fHistPsiPairTrue =
          new TH2F("fHistPsiPairTrue", "; #it{p}_{T} (GeV/#it{c}); #Psi_{pair}",
                   100, 0, 10, 1000, -pi, pi);
      fHistogramsMC->Add(fHistPsiPairTrue);

      fHistSingleParticlePtTrue[0] =
          new TH1F("fHistSingleParticlePtTrue_pos",
                   ";#it{p}_{T} true ; Entries", 1000, 0, 10);
      fHistSingleParticlePtTrue[1] =
          new TH1F("fHistSingleParticlePtTrue_neg", ";#it{p}_{T} true; Entries",
                   1000, 0, 10);
      fHistogramsMC->Add(fHistSingleParticlePtTrue[0]);
      fHistogramsMC->Add(fHistSingleParticlePtTrue[1]);

      fHistSingleParticleDCAtoPVTrue[0] =
          new TH2F("fHistSingleParticleDCAtoPVTrue_pos",
                   "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA) true", 100, 0, 10,
                   500, 0, 10);
      fHistSingleParticleDCAtoPVTrue[1] =
          new TH2F("fHistSingleParticleDCAtoPVTrue_neg",
                   "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA) true", 100, 0, 10,
                   500, 0, 10);
      fHistogramsMC->Add(fHistSingleParticleDCAtoPVTrue[0]);
      fHistogramsMC->Add(fHistSingleParticleDCAtoPVTrue[1]);

      fHistSingleParticleNclsTPCTrue[0] = new TH2F(
          "fHistSingleParticleNclsTPCTrue_pos",
          "; #it{p}_{T} (GeV/#it{c}); # cls TPC true", 100, 0, 10, 160, 0, 160);
      fHistSingleParticleNclsTPCTrue[1] = new TH2F(
          "fHistSingleParticleNclsTPCTrue_neg",
          "; #it{p}_{T} (GeV/#it{c}); # cls TPC true", 100, 0, 10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCTrue[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCTrue[1]);

      fHistSingleParticleNclsTPCFindableTrue[0] =
          new TH2F("fHistSingleParticleNclsTPCFindableTrue_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable true", 100, 0,
                   10, 160, 0, 160);
      fHistSingleParticleNclsTPCFindableTrue[1] =
          new TH2F("fHistSingleParticleNclsTPCFindableTrue_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable true", 100, 0,
                   10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCFindableTrue[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCFindableTrue[1]);

      fHistSingleParticleNclsTPCRatioFindableTrue[0] =
          new TH2F("fHistSingleParticleNclsTPCRatioFindableTrue_pos",
                   "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable true", 100, 0,
                   10, 500, 0, 2);
      fHistSingleParticleNclsTPCRatioFindableTrue[1] =
          new TH2F("fHistSingleParticleNclsTPCRatioFindableTrue_neg",
                   "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable true", 100, 0,
                   10, 500, 0, 2);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCRatioFindableTrue[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCRatioFindableTrue[1]);

      fHistSingleParticleNcrossedTPCTrue[0] =
          new TH2F("fHistSingleParticleNcrossedTPCTrue_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed true", 100, 0,
                   10, 160, 0, 160);
      fHistSingleParticleNcrossedTPCTrue[1] =
          new TH2F("fHistSingleParticleNcrossedTPCTrue_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed true", 100, 0,
                   10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNcrossedTPCTrue[0]);
      fHistogramsMC->Add(fHistSingleParticleNcrossedTPCTrue[1]);

      fHistSingleParticleNclsTPCSharedTrue[0] =
          new TH2F("fHistSingleParticleNclsTPCSharedTrue_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared true", 100, 0,
                   10, 160, 0, 160);
      fHistSingleParticleNclsTPCSharedTrue[1] =
          new TH2F("fHistSingleParticleNclsTPCSharedTrue_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared true", 100, 0,
                   10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCSharedTrue[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCSharedTrue[1]);

      fHistSingleParticleNclsITSSharedTrue[0] =
          new TH2F("fHistSingleParticleNclsITSSharedTrue_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared true", 100, 0,
                   10, 6, 0, 6);
      fHistSingleParticleNclsITSSharedTrue[1] =
          new TH2F("fHistSingleParticleNclsITSSharedTrue_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared true", 100, 0,
                   10, 6, 0, 6);
      fHistogramsMC->Add(fHistSingleParticleNclsITSSharedTrue[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsITSSharedTrue[1]);

      fHistSingleParticleChi2True[0] = new TH2F(
          "fHistSingleParticleChi2True_pos",
          "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
      fHistSingleParticleChi2True[1] = new TH2F(
          "fHistSingleParticleChi2True_neg",
          "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
      fHistogramsMC->Add(fHistSingleParticleChi2True[0]);
      fHistogramsMC->Add(fHistSingleParticleChi2True[1]);

      fHistSingleParticlePIDTrue[0] = new TH2F(
          "fHistSingleParticlePIDTrue_pos",
          "; #it{p} (GeV/#it{c}); n_{#sigma} TPC", 100, 0, 10, 500, -10, 10);
      fHistSingleParticlePIDTrue[1] = new TH2F(
          "fHistSingleParticlePIDTrue_neg",
          "; #it{p} (GeV/#it{c}); n_{#sigma} TPC", 100, 0, 10, 500, -10, 10);
      fHistogramsMC->Add(fHistSingleParticlePIDTrue[0]);
      fHistogramsMC->Add(fHistSingleParticlePIDTrue[1]);

      // TRUE SIGMA DAUGHTER
      fHistV0MassPtTrueSigma =
          new TH2F("fHistV0MassPtTrueSigma",
                   "; #it{p}_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^{2})",
                   100, 0, 10, 1200, 0., 1.2);
      fHistogramsMC->Add(fHistV0MassPtTrueSigma);

      fHistDecayVertexXTrueSigma =
          new TH2F("fHistDecayVertexXTrueSigma",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex x true Sigma (cm)",
                   100, 0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexXTrueSigma);

      fHistDecayVertexYTrueSigma =
          new TH2F("fHistDecayVertexYTrueSigma",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex y true Sigma (cm)",
                   100, 0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexYTrueSigma);

      fHistDecayVertexZTrueSigma =
          new TH2F("fHistDecayVertexZTrueSigma",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex z true Sigma (cm)",
                   100, 0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexZTrueSigma);

      fHistTransverseRadiusTrueSigma = new TH2F(
          "fHistTransverseRadiusTrueSigma",
          "; #it{p}_{T} (GeV/#it{c}); Transverse radius true Sigma (cm)", 100,
          0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistTransverseRadiusTrueSigma);

      fHistCosPATrueSigma =
          new TH2F("fHistCosPATrueSigma",
                   "; #it{p}_{T} (GeV/#it{c}); cos(#alpha) true Sigma", 100, 0,
                   10, 500, 0.8, 1);
      fHistogramsMC->Add(fHistCosPATrueSigma);

      fHistDCADaughtersTrueSigma =
          new TH2F("fHistDCADaughtersTrueSigma",
                   "; #it{p}_{T} (GeV/#it{c}); Daughter DCA at decay vertex "
                   "true Sigma (cm)",
                   100, 0, 10, 500, 0, 2);
      fHistogramsMC->Add(fHistDCADaughtersTrueSigma);

      fHistArmenterosTrueSigma =
          new TH2F("fHistArmenterosTrueSigma",
                   " ; #alpha; #it{q}_{T} (GeV/#it{c})", 500, -1, 1, 500, 0, 1);
      fHistogramsMC->Add(fHistArmenterosTrueSigma);

      fHistPsiPairTrueSigma = new TH2F("fHistPsiPairTrueSigma",
                                       "; #it{p}_{T} (GeV/#it{c}); #Psi_{pair}",
                                       100, 0, 10, 1000, -pi, pi);
      fHistogramsMC->Add(fHistPsiPairTrueSigma);

      fHistSingleParticlePtTrueSigma[0] =
          new TH1F("fHistSingleParticlePtTrueSigma_pos",
                   ";#it{p}_{T} true Sigma; Entries", 100, 0, 10);
      fHistSingleParticlePtTrueSigma[1] =
          new TH1F("fHistSingleParticlePtTrueSigma_neg",
                   ";#it{p}_{T} true Sigma; Entries", 100, 0, 10);
      fHistogramsMC->Add(fHistSingleParticlePtTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticlePtTrueSigma[1]);

      fHistSingleParticleDCAtoPVTrueSigma[0] =
          new TH2F("fHistSingleParticleDCAtoPVTrueSigma_pos",
                   "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA) true Sigma", 100, 0,
                   10, 500, 0, 10);
      fHistSingleParticleDCAtoPVTrueSigma[1] =
          new TH2F("fHistSingleParticleDCAtoPVTrueSigma_neg",
                   "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA) true Sigma", 100, 0,
                   10, 500, 0, 10);
      fHistogramsMC->Add(fHistSingleParticleDCAtoPVTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticleDCAtoPVTrueSigma[1]);

      fHistSingleParticleNclsTPCTrueSigma[0] =
          new TH2F("fHistSingleParticleNclsTPCTrueSigma_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC true Sigma", 100, 0,
                   10, 160, 0, 160);
      fHistSingleParticleNclsTPCTrueSigma[1] =
          new TH2F("fHistSingleParticleNclsTPCTrueSigma_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC true Sigma", 100, 0,
                   10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCTrueSigma[1]);

      fHistSingleParticleNclsTPCFindableTrueSigma[0] =
          new TH2F("fHistSingleParticleNclsTPCFindableTrueSigma_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable true Sigma",
                   100, 0, 10, 160, 0, 160);
      fHistSingleParticleNclsTPCFindableTrueSigma[1] =
          new TH2F("fHistSingleParticleNclsTPCFindableTrueSigma_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable true Sigma",
                   100, 0, 10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCFindableTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCFindableTrueSigma[1]);

      fHistSingleParticleNclsTPCRatioFindableTrueSigma[0] =
          new TH2F("fHistSingleParticleNclsTPCRatioFindableTrueSigma_pos",
                   "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable true Sigma",
                   100, 0, 10, 500, 0, 2);
      fHistSingleParticleNclsTPCRatioFindableTrueSigma[1] =
          new TH2F("fHistSingleParticleNclsTPCRatioFindableTrueSigma_neg",
                   "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable true Sigma",
                   100, 0, 10, 500, 0, 2);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCRatioFindableTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCRatioFindableTrueSigma[1]);

      fHistSingleParticleNcrossedTPCTrueSigma[0] =
          new TH2F("fHistSingleParticleNcrossedTPCTrueSigma_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed true Sigma",
                   100, 0, 10, 160, 0, 160);
      fHistSingleParticleNcrossedTPCTrueSigma[1] =
          new TH2F("fHistSingleParticleNcrossedTPCTrueSigma_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed true Sigma",
                   100, 0, 10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNcrossedTPCTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticleNcrossedTPCTrueSigma[1]);

      fHistSingleParticleNclsTPCSharedTrueSigma[0] =
          new TH2F("fHistSingleParticleNclsTPCSharedTrueSigma_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared true Sigma",
                   100, 0, 10, 160, 0, 160);
      fHistSingleParticleNclsTPCSharedTrueSigma[1] =
          new TH2F("fHistSingleParticleNclsTPCSharedTrueSigma_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared true Sigma",
                   100, 0, 10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCSharedTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCSharedTrueSigma[1]);

      fHistSingleParticleNclsITSSharedTrueSigma[0] =
          new TH2F("fHistSingleParticleNclsITSSharedTrueSigma_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared true Sigma",
                   100, 0, 10, 6, 0, 6);
      fHistSingleParticleNclsITSSharedTrueSigma[1] =
          new TH2F("fHistSingleParticleNclsITSSharedTrueSigma_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared true Sigma",
                   100, 0, 10, 6, 0, 6);
      fHistogramsMC->Add(fHistSingleParticleNclsITSSharedTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsITSSharedTrueSigma[1]);

      fHistSingleParticleChi2TrueSigma[0] = new TH2F(
          "fHistSingleParticleChi2TrueSigma_pos",
          "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
      fHistSingleParticleChi2TrueSigma[1] = new TH2F(
          "fHistSingleParticleChi2TrueSigma_neg",
          "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
      fHistogramsMC->Add(fHistSingleParticleChi2TrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticleChi2TrueSigma[1]);

      fHistSingleParticlePIDTrueSigma[0] = new TH2F(
          "fHistSingleParticlePIDTrueSigma_pos",
          "; #it{p} (GeV/#it{c}); n_{#sigma} TPC", 100, 0, 10, 500, -10, 10);
      fHistSingleParticlePIDTrueSigma[1] = new TH2F(
          "fHistSingleParticlePIDTrueSigma_neg",
          "; #it{p} (GeV/#it{c}); n_{#sigma} TPC", 100, 0, 10, 500, -10, 10);
      fHistogramsMC->Add(fHistSingleParticlePIDTrueSigma[0]);
      fHistogramsMC->Add(fHistSingleParticlePIDTrueSigma[1]);

      // BACKGROUND
      fHistV0MassPtBkg =
          new TH2F("fHistV0MassPtBkg",
                   "; #it{p}_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^{2})",
                   100, 0, 10, 1200, 0., 1.2);
      fHistogramsMC->Add(fHistV0MassPtBkg);

      fHistDecayVertexXBkg =
          new TH2F("fHistDecayVertexXBkg",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex x Bkg (cm)", 100, 0,
                   10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexXBkg);

      fHistDecayVertexYBkg =
          new TH2F("fHistDecayVertexYBkg",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex y Bkg (cm)", 100, 0,
                   10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexYBkg);

      fHistDecayVertexZBkg =
          new TH2F("fHistDecayVertexZBkg",
                   "; #it{p}_{T} (GeV/#it{c}); Decay vertex z Bkg (cm)", 100, 0,
                   10, 500, 0, 200);
      fHistogramsMC->Add(fHistDecayVertexZBkg);

      fHistTransverseRadiusBkg =
          new TH2F("fHistTransverseRadiusBkg",
                   "; #it{p}_{T} (GeV/#it{c}); Transverse radius Bkg (cm)", 100,
                   0, 10, 500, 0, 200);
      fHistogramsMC->Add(fHistTransverseRadiusBkg);

      fHistCosPABkg = new TH2F("fHistCosPABkg",
                               "; #it{p}_{T} (GeV/#it{c}); cos(#alpha) Bkg",
                               100, 0, 10, 500, 0.8, 1);
      fHistogramsMC->Add(fHistCosPABkg);

      fHistDCADaughtersBkg = new TH2F(
          "fHistDCADaughtersBkg",
          "; #it{p}_{T} (GeV/#it{c}); Daughter DCA at decay vertex Bkg (cm)",
          100, 0, 10, 500, 0, 2);
      fHistogramsMC->Add(fHistDCADaughtersBkg);

      fHistArmenterosBkg =
          new TH2F("fHistArmenterosBkg", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                   500, -1, 1, 500, 0, 1);
      fHistogramsMC->Add(fHistArmenterosBkg);

      fHistPsiPairBkg =
          new TH2F("fHistPsiPairBkg", "; #it{p}_{T} (GeV/#it{c}); #Psi_{pair}",
                   1000, 0, 10, 1000, -pi, pi);
      fHistogramsMC->Add(fHistPsiPairBkg);

      fHistSingleParticlePtBkg[0] =
          new TH1F("fHistSingleParticlePtBkg_pos", ";#it{p}_{T} Bkg; Entries",
                   1000, 0, 10);
      fHistSingleParticlePtBkg[1] =
          new TH1F("fHistSingleParticlePtBkg_neg", ";#it{p}_{T} Bkg; Entries",
                   1000, 0, 10);
      fHistogramsMC->Add(fHistSingleParticlePtBkg[0]);
      fHistogramsMC->Add(fHistSingleParticlePtBkg[1]);

      fHistSingleParticleDCAtoPVBkg[0] =
          new TH2F("fHistSingleParticleDCAtoPVBkg_pos",
                   "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA) Bkg", 100, 0, 10,
                   500, 0, 10);
      fHistSingleParticleDCAtoPVBkg[1] =
          new TH2F("fHistSingleParticleDCAtoPVBkg_neg",
                   "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA) Bkg", 100, 0, 10,
                   500, 0, 10);
      fHistogramsMC->Add(fHistSingleParticleDCAtoPVBkg[0]);
      fHistogramsMC->Add(fHistSingleParticleDCAtoPVBkg[1]);

      fHistSingleParticleNclsTPCBkg[0] = new TH2F(
          "fHistSingleParticleNclsTPCBkg_pos",
          "; #it{p}_{T} (GeV/#it{c}); # cls TPC Bkg", 100, 0, 10, 160, 0, 160);
      fHistSingleParticleNclsTPCBkg[1] = new TH2F(
          "fHistSingleParticleNclsTPCBkg_neg",
          "; #it{p}_{T} (GeV/#it{c}); # cls TPC Bkg", 100, 0, 10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCBkg[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCBkg[1]);

      fHistSingleParticleNclsTPCFindableBkg[0] =
          new TH2F("fHistSingleParticleNclsTPCFindableBkg_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable Bkg", 100, 0,
                   10, 160, 0, 160);
      fHistSingleParticleNclsTPCFindableBkg[1] =
          new TH2F("fHistSingleParticleNclsTPCFindableBkg_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable Bkg", 100, 0,
                   10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCFindableBkg[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCFindableBkg[1]);

      fHistSingleParticleNclsTPCRatioFindableBkg[0] =
          new TH2F("fHistSingleParticleNclsTPCRatioFindableBkg_pos",
                   "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable Bkg", 100, 0,
                   10, 500, 0, 2);
      fHistSingleParticleNclsTPCRatioFindableBkg[1] =
          new TH2F("fHistSingleParticleNclsTPCRatioFindableBkg_neg",
                   "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable Bkg", 100, 0,
                   10, 500, 0, 2);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCRatioFindableBkg[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCRatioFindableBkg[1]);

      fHistSingleParticleNcrossedTPCBkg[0] =
          new TH2F("fHistSingleParticleNcrossedTPCBkg_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed Bkg", 100, 0,
                   10, 160, 0, 160);
      fHistSingleParticleNcrossedTPCBkg[1] =
          new TH2F("fHistSingleParticleNcrossedTPCBkg_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed Bkg", 100, 0,
                   10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNcrossedTPCBkg[0]);
      fHistogramsMC->Add(fHistSingleParticleNcrossedTPCBkg[1]);

      fHistSingleParticleNclsTPCSharedBkg[0] =
          new TH2F("fHistSingleParticleNclsTPCSharedBkg_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared Bkg", 100, 0,
                   10, 160, 0, 160);
      fHistSingleParticleNclsTPCSharedBkg[1] =
          new TH2F("fHistSingleParticleNclsTPCSharedBkg_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared Bkg", 100, 0,
                   10, 160, 0, 160);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCSharedBkg[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsTPCSharedBkg[1]);

      fHistSingleParticleNclsITSSharedBkg[0] =
          new TH2F("fHistSingleParticleNclsITSSharedBkg_pos",
                   "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared Bkg", 100, 0,
                   10, 6, 0, 6);
      fHistSingleParticleNclsITSSharedBkg[1] =
          new TH2F("fHistSingleParticleNclsITSSharedBkg_neg",
                   "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared Bkg", 100, 0,
                   10, 6, 0, 6);
      fHistogramsMC->Add(fHistSingleParticleNclsITSSharedBkg[0]);
      fHistogramsMC->Add(fHistSingleParticleNclsITSSharedBkg[1]);

      fHistSingleParticleChi2Bkg[0] = new TH2F(
          "fHistSingleParticleChi2Bkg_pos",
          "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
      fHistSingleParticleChi2Bkg[1] = new TH2F(
          "fHistSingleParticleChi2Bkg_neg",
          "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
      fHistogramsMC->Add(fHistSingleParticleChi2Bkg[0]);
      fHistogramsMC->Add(fHistSingleParticleChi2Bkg[1]);

      fHistSingleParticlePIDBkg[0] = new TH2F(
          "fHistSingleParticlePIDBkg_pos",
          "; #it{p} (GeV/#it{c}); n_{#sigma} TPC", 100, 0, 10, 500, -10, 10);
      fHistSingleParticlePIDBkg[1] = new TH2F(
          "fHistSingleParticlePIDBkg_neg",
          "; #it{p} (GeV/#it{c}); n_{#sigma} TPC", 100, 0, 10, 500, -10, 10);
      fHistogramsMC->Add(fHistSingleParticlePIDBkg[0]);
      fHistogramsMC->Add(fHistSingleParticlePIDBkg[1]);
    }

    fHistograms->Add(fHistogramsMC);
  }
}
