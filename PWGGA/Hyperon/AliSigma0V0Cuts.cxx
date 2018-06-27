#include "AliSigma0V0Cuts.h"
#include <iostream>
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

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
      fPDGDatabase(),
      fIsLightweight(false),
      fV0cut(0),
      fAntiV0cut(0),
      fPID(0),
      fIsMC(false),
      fPileUpRejectionMode(PileUpRejectionMode::BothDaughtersCombined),
      fPosPID(AliPID::kProton),
      fNegPID(AliPID::kPion),
      fV0OnFly(false),
      fK0Rejection(true),
      fUsePID(true),
      fV0PtMin(0.f),
      fV0PtMax(1E30),
      fV0CosPAMin(0.f),
      fV0RadiusMax(999.f),
      fV0RadiusMin(0.f),
      fV0DecayVertexMax(999.f),
      fPIDnSigma(999.f),
      fEtaMax(999.f),
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
      fPIDResponse(nullptr),
      fHistCutBooking(nullptr),
      fHistCuts(nullptr),
      fHistNV0(nullptr),
      fHistLambdaMass(nullptr),
      fHistLambdaPt(nullptr),
      fHistLambdaPtY(),
      fHistLambdaMassPt(nullptr),
      fHistLambdaMassK0Rej(nullptr),
      fHistK0Mass(nullptr),
      fHistK0MassAfter(nullptr),
      fHistCosPA(nullptr),
      fHistEtaPhi(nullptr),
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
      fHistMCTruthV0Pt(nullptr),
      fHistMCTruthV0PtY(nullptr),
      fHistMCTruthV0PtEta(nullptr),
      fHistMCTruthV0ProtonPionPt(nullptr),
      fHistMCTruthV0ProtonPionPtY(nullptr),
      fHistMCTruthV0ProtonPionPtEta(nullptr),
      fHistSingleParticleCuts(),
      fHistSingleParticlePt(),
      fHistSingleParticleEtaBefore(),
      fHistSingleParticleEtaAfter(),
      fHistSingleParticleNclsTPCBefore(),
      fHistSingleParticleNclsTPCAfter(),
      fHistSingleParticleNclsTPCFindableBefore(),
      fHistSingleParticleNclsTPCFindableAfter(),
      fHistSingleParticleNclsTPCRatioFindableBefore(),
      fHistSingleParticleNclsTPCRatioFindableAfter(),
      fHistSingleParticleNcrossedTPCBefore(),
      fHistSingleParticleNcrossedTPCAfter(),
      fHistSingleParticleNclsTPCShared(),
      fHistSingleParticleNclsITSShared(),
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
      fPDGDatabase(),
      fIsLightweight(false),
      fV0cut(0),
      fAntiV0cut(0),
      fPID(0),
      fIsMC(false),
      fPileUpRejectionMode(PileUpRejectionMode::BothDaughtersCombined),
      fPosPID(AliPID::kProton),
      fNegPID(AliPID::kPion),
      fV0OnFly(false),
      fK0Rejection(true),
      fUsePID(true),
      fV0PtMin(0.f),
      fV0PtMax(1E30),
      fV0CosPAMin(0.f),
      fV0RadiusMax(999.f),
      fV0RadiusMin(0.f),
      fV0DecayVertexMax(999.f),
      fPIDnSigma(999.f),
      fEtaMax(999.f),
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
      fPIDResponse(nullptr),
      fHistCutBooking(nullptr),
      fHistCuts(nullptr),
      fHistNV0(nullptr),
      fHistLambdaMass(nullptr),
      fHistLambdaPt(nullptr),
      fHistLambdaPtY(),
      fHistLambdaMassPt(nullptr),
      fHistLambdaMassK0Rej(nullptr),
      fHistK0Mass(nullptr),
      fHistK0MassAfter(nullptr),
      fHistCosPA(nullptr),
      fHistEtaPhi(nullptr),
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
      fHistMCTruthV0Pt(nullptr),
      fHistMCTruthV0PtY(nullptr),
      fHistMCTruthV0PtEta(nullptr),
      fHistMCTruthV0ProtonPionPt(nullptr),
      fHistMCTruthV0ProtonPionPtY(nullptr),
      fHistMCTruthV0ProtonPionPtEta(nullptr),
      fHistSingleParticleCuts(),
      fHistSingleParticlePt(),
      fHistSingleParticleEtaBefore(),
      fHistSingleParticleEtaAfter(),
      fHistSingleParticleNclsTPCBefore(),
      fHistSingleParticleNclsTPCAfter(),
      fHistSingleParticleNclsTPCFindableBefore(),
      fHistSingleParticleNclsTPCFindableAfter(),
      fHistSingleParticleNclsTPCRatioFindableBefore(),
      fHistSingleParticleNclsTPCRatioFindableAfter(),
      fHistSingleParticleNcrossedTPCBefore(),
      fHistSingleParticleNcrossedTPCAfter(),
      fHistSingleParticleNclsTPCShared(),
      fHistSingleParticleNclsITSShared(),
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
  v0Cuts->SetEtaMax(0.8);
  v0Cuts->SetDaughterDCAMax(1.5);
  v0Cuts->SetDaughterDCAtoPV(0.05);
  v0Cuts->SetK0Rejection(0.48, 0.515);
  v0Cuts->SetLambdaSelection(1.115683 - 0.004, 1.115683 + 0.004);
  v0Cuts->SetPileUpRejectionMode(BothDaughtersCombined);
  return v0Cuts;
}

//____________________________________________________________________________________________________
AliSigma0V0Cuts *AliSigma0V0Cuts::PhotonCuts() {
  AliSigma0V0Cuts *v0Cuts = new AliSigma0V0Cuts();
  v0Cuts->SetV0OnFlyStatus(false);
  v0Cuts->SetV0PtMin(0.3);
  v0Cuts->SetV0CosPAMin(0.99);
  v0Cuts->SetV0RadiusMax(200.f);
  v0Cuts->SetV0RadiusMin(0.);
  v0Cuts->SetPIDnSigma(5.f);
  v0Cuts->SetTPCRatioFindable(0.6f);
  v0Cuts->SetEtaMax(0.8);
  v0Cuts->SetDaughterDCAMax(1.5);
  v0Cuts->SetDaughterDCAtoPV(0.05);
  v0Cuts->SetPileUpRejectionMode(None);
  return v0Cuts;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::SelectLambda(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    std::vector<AliSigma0ParticleV0> &V0Container) {
  // Get the PID manager from the input event handler
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      static_cast<AliInputEventHandler *>(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fMCEvent = mcEvent;

  V0Container.clear();

  fInputEvent = static_cast<AliESDEvent *>(inputEvent);
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

    // PID
    if (!V0PID(v0, pos, neg, fPosPID, fNegPID)) continue;

    // Single Particle Quality
    if (!SingleParticleQualityCuts(pos) || !SingleParticleQualityCuts(neg))
      continue;
    fHistCuts->Fill(7);

    // Pile-up rejection
    if (!PileUpRejection(pos, neg)) continue;
    fHistCuts->Fill(8);

    // Topological Selection
    if (!V0TopologicalSelection(v0)) continue;

    // Lambda Selection
    if (!LambdaSelection(v0)) continue;

    AliSigma0ParticleV0 v0Candidate(v0, pos, neg,
                                    fInputEvent->GetPrimaryVertex(), fPID,
                                    fInputEvent->GetMagneticField(), fMCEvent);

    if (fIsMC) {
    }

    v0Candidate.SetPDGMass(fPDGDatabase.GetParticle(fPID)->Mass());
    V0Container.push_back(v0Candidate);
  }

  fHistNV0->Fill(V0Container.size());
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::SelectPhoton(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    std::vector<AliSigma0ParticleV0> &V0Container,
    AliPID::EParticleType particle1, AliPID::EParticleType particle2) {
  // Get the PID manager from the input event handler
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      static_cast<AliInputEventHandler *>(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fMCEvent = mcEvent;

  V0Container.clear();

  fInputEvent = static_cast<AliESDEvent *>(inputEvent);
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

    // PID
    /// DO BE DONE - DO PROPER PID
    if (!V0PID(v0, pos, neg, particle1, particle2)) continue;

    // Single Particle Quality
    if (!SingleParticleQualityCuts(pos) || !SingleParticleQualityCuts(neg))
      continue;
    fHistCuts->Fill(7);

    // Pile-up rejection
    if (!PileUpRejection(pos, neg)) continue;
    fHistCuts->Fill(8);

    // Topological Selection
    if (!V0TopologicalSelection(v0)) continue;

    // PHOTON SELECTION - TO BE DONE
    // REMOVE REMAINING LAMBDAS, CUT ON PSI_PAIR

    AliSigma0ParticleV0 v0Candidate(v0, pos, neg,
                                    fInputEvent->GetPrimaryVertex(), fPID,
                                    fInputEvent->GetMagneticField(), fMCEvent);

    if (fIsMC) {
    }

    v0Candidate.SetPDGMass(fPDGDatabase.GetParticle(fPID)->Mass());
    V0Container.push_back(v0Candidate);
  }
  fHistNV0->Fill(V0Container.size());
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0QualityCuts(const AliESDv0 *v0) {
  fHistCuts->Fill(0);

  // v0 is neutral
  if (v0->Charge() != 0) return false;
  fHistCuts->Fill(1);

  // on-fly status
  if (fV0OnFly) {
    // online v0
    if (!(v0->GetOnFlyStatus())) return false;
  } else {
    // offline v0
    if (v0->GetOnFlyStatus()) return false;
  }
  fHistCuts->Fill(2);

  // pt min cut
  if (v0->Pt() < fV0PtMin) return false;
  fHistCuts->Fill(3);

  // pt max cut
  if (v0->Pt() > fV0PtMax) return false;
  fHistCuts->Fill(4);

  // check for daughter tracks
  AliESDEvent *event = static_cast<AliESDEvent *>(fInputEvent);
  AliESDtrack *pos = event->GetTrack(v0->GetPindex());
  AliESDtrack *neg = event->GetTrack(v0->GetNindex());
  if (!pos || !neg) return false;
  if (pos->Charge() == neg->Charge()) return false;
  fHistCuts->Fill(5);

  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0PID(const AliESDv0 *v0, const AliESDtrack *pos,
                            const AliESDtrack *neg,
                            AliPID::EParticleType pidPos,
                            AliPID::EParticleType pidNeg) {
  float posProb = 999.f;
  float negProb = 999.f;

  const float armAlpha = v0->AlphaV0();
  const float armQt = v0->PtArmV0();
  if (!fIsLightweight) fHistArmenterosBefore->Fill(armAlpha, armQt);

  // When true, we use the full PID capabilities, when false just rough
  // Armenteros-Podolandksi cut
  if (fUsePID) {
    bool isParticle = (SingleParticlePID(pos, pidPos, posProb) &&
                       SingleParticlePID(neg, pidNeg, negProb));

    if (!isParticle) return false;
  } else {
    // Armenteros cut
    if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) return false;
    float prefactorAlpha = (fPID > 0) ? 1.f : -1.f;  // for anti-particles
                                                     // negative to compensate
                                                     // for sign change of qt
    if (armAlpha * prefactorAlpha > fArmenterosAlphaUp ||
        armAlpha * prefactorAlpha < fArmenterosAlphaLow)
      return false;
  }

  fHistCuts->Fill(6);
  if (!fIsLightweight) {
    fHistArmenterosAfter->Fill(armAlpha, armQt);
    PlotSingleParticlePID(pos, pidPos);
    PlotSingleParticlePID(neg, pidNeg);
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

  prob = std::abs(fPIDResponse->NumberOfSigmasTPC(track, particle));

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

  const float nSigma =
      std::abs(fPIDResponse->NumberOfSigmasTPC(track, particle));

  const int charge = track->Charge();
  const int histoPrefix = (charge > 0) ? 0 : 1;
  if (!fIsLightweight)
    fHistSingleParticlePID[histoPrefix]->Fill(track->P(), nSigma);
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::SingleParticleQualityCuts(AliESDtrack *track) {
  const float pt = track->Pt();
  const int charge = track->Charge();
  const float eta = track->Eta();
  const short nClsTPC = track->GetTPCNcls();
  const float nCrossedRows = track->GetTPCClusterInfo(2, 1);
  const short nFindable = track->GetTPCNclsF();
  const float ratioFindable = nCrossedRows / static_cast<float>(nFindable);

  const float magField = fInputEvent->GetMagneticField();
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  const float dcaDaughterToPV =
      std::abs(track->GetD(vertex->GetX(), vertex->GetY(), magField));

  const int histoPrefix = (charge > 0) ? 0 : 1;
  int qaHistoCounter = 0;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (!fIsLightweight) {
    fHistSingleParticleEtaBefore[histoPrefix]->Fill(eta);
    fHistSingleParticleNclsTPCBefore[histoPrefix]->Fill(nClsTPC);
    fHistSingleParticleNclsTPCFindableBefore[histoPrefix]->Fill(nFindable);
    fHistSingleParticleNclsTPCRatioFindableBefore[histoPrefix]->Fill(
        ratioFindable);
    fHistSingleParticleNcrossedTPCBefore[histoPrefix]->Fill(nCrossedRows);
    fHistSingleParticleDCAtoPVBefore[histoPrefix]->Fill(dcaDaughterToPV);
  }

  // Max eta cut
  if (std::abs(eta) > fEtaMax) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC nCluster cut
  if (nClsTPC < fTPCclusterMin) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC Crossed rows
  if (nCrossedRows < fTPCnCrossedRowsMin) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC Ratio Findable
  if (ratioFindable < fTPCratioFindable) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC nCls Findable
  if (nFindable < fTPCfindableMin) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // TPC nCls Shared
  const int nClsSharedTPC = track->GetTPCnclsS();
  if (nClsSharedTPC > fTPCnSharedMax) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // Minimal distance of daughters to primary vertex
  if (dcaDaughterToPV < fDaughterDCAPV) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // Keep track of shared clusters
  if (!fIsLightweight) {
    int nClsSharedITS = 0;
    for (int i = 0; i < 6; ++i) {
      if (track->HasSharedPointOnITSLayer(i)) ++nClsSharedITS;
    }
    const int histoPrefix = (charge > 0) ? 0 : 1;
    fHistSingleParticleNclsTPCShared[histoPrefix]->Fill(nClsSharedTPC);
    fHistSingleParticleNclsITSShared[histoPrefix]->Fill(nClsSharedITS);
  }

  if (!(track->GetStatus() & AliESDtrack::kTPCrefit)) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (track->GetKinkIndex(0) > 0) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (!fIsLightweight) {
    fHistSingleParticlePt[histoPrefix]->Fill(pt);
    fHistSingleParticleEtaAfter[histoPrefix]->Fill(eta);
    fHistSingleParticleNclsTPCAfter[histoPrefix]->Fill(nClsTPC);
    fHistSingleParticleNclsTPCFindableAfter[histoPrefix]->Fill(nFindable);
    fHistSingleParticleNclsTPCRatioFindableAfter[histoPrefix]->Fill(
        ratioFindable);
    fHistSingleParticleNcrossedTPCAfter[histoPrefix]->Fill(nCrossedRows);
    fHistSingleParticleDCAtoPVAfter[histoPrefix]->Fill(dcaDaughterToPV);
  }
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::PileUpRejection(AliESDtrack *pos, AliESDtrack *neg) {
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
bool AliSigma0V0Cuts::V0TopologicalSelection(const AliESDv0 *v0) {
  // Get the coordinates of the primary vertex
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
    fHistCosPABefore->Fill(point);
    fHistDecayVertexXBefore->Fill(decayVertexV0[0]);
    fHistDecayVertexYBefore->Fill(decayVertexV0[1]);
    fHistDecayVertexZBefore->Fill(decayVertexV0[2]);
    fHistTransverseRadiusBefore->Fill(transverseRadius);
    fHistDCADaughtersBefore->Fill(dcaV0Dau);
  }

  // Position of the decay vertex x, y & z
  if (std::abs(decayVertexV0[0]) > fV0DecayVertexMax) return false;
  fHistCuts->Fill(9);
  if (std::abs(decayVertexV0[1]) > fV0DecayVertexMax) return false;
  fHistCuts->Fill(10);
  if (std::abs(decayVertexV0[2]) > fV0DecayVertexMax) return false;
  fHistCuts->Fill(11);

  // Transverse decay radius min & max
  if (transverseRadius < fV0RadiusMin) return false;
  fHistCuts->Fill(12);

  if (transverseRadius > fV0RadiusMax) return false;
  fHistCuts->Fill(13);

  // DCA of the daughter tracks at the decay vertex
  if (dcaV0Dau > fDaughterDCAMax) return false;
  fHistCuts->Fill(14);

  // Cos pointing angle
  if (point < fV0CosPAMin) return false;
  fHistCuts->Fill(15);

  if (!fIsLightweight) {
    fHistCosPA->Fill(v0->Pt(), point);
    fHistCosPAAfter->Fill(point);
    fHistDecayVertexXAfter->Fill(decayVertexV0[0]);
    fHistDecayVertexYAfter->Fill(decayVertexV0[1]);
    fHistDecayVertexZAfter->Fill(decayVertexV0[2]);
    fHistTransverseRadiusAfter->Fill(transverseRadius);
    fHistDCADaughtersAfter->Fill(dcaV0Dau);
    fHistDCA->Fill(dcaPrim);
    fHistDecayLength->Fill(lenDecay);
  }
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::LambdaSelection(AliESDv0 *v0) {
  const float armAlpha = v0->AlphaV0();
  const float armQt = v0->PtArmV0();

  v0->ChangeMassHypothesis(fPID);
  const float massLambda = v0->GetEffMass();
  v0->ChangeMassHypothesis(310);
  const float massK0 = v0->GetEffMass();

  if (!fIsLightweight) {
    fHistLambdaMassK0Rej->Fill(massLambda);
    fHistK0Mass->Fill(massK0);
  }

  // K0 rejection cut
  if (fK0Rejection && massK0 > fK0RejectionLow && massK0 < fK0RejectionUp)
    return false;

  float rap = ComputeRapidity(v0->Pt(), v0->Pz(), massLambda);
  int rapBin = GetRapidityBin(rap);

  fHistLambdaMassPt->Fill(v0->Pt(), massLambda);
  if (!fIsLightweight) {
    fHistLambdaMass->Fill(massLambda);
    if (rapBin > -1) fHistLambdaPtY[rapBin]->Fill(v0->Pt(), massLambda);
    fHistK0MassAfter->Fill(massK0);
    fHistEtaPhi->Fill(v0->Eta(), v0->Phi());
  }

  // Lambda selection cut
  if (massLambda < fLambdaSelectionLow || massLambda > fLambdaSelectionUp)
    return false;

  if (!fIsLightweight) fHistLambdaPt->Fill(v0->Pt());
  fHistCuts->Fill(16);
  return true;
}

//____________________________________________________________________________________________________
float AliSigma0V0Cuts::ComputeRapidity(float pt, float pz, float m) const {
  // calculates rapidity keeping the sign in case E == pz

  float energy = std::sqrt(pt * pt + pz * pz + m * m);
  if (energy != std::fabs(pz))
    return 0.5 * std::log((energy + pz) / (energy - pz));
  return (pz >= 0) ? 1.e30 : -1.e30;
}

//____________________________________________________________________________________________________
int AliSigma0V0Cuts::GetRapidityBin(float rapidity) const {
  if (-1.f < rapidity && rapidity <= -0.9)
    return 0;
  else if (-0.9 < rapidity && rapidity <= -0.8)
    return 1;
  else if (-0.8 < rapidity && rapidity <= -0.7)
    return 2;
  else if (-0.7 < rapidity && rapidity <= -0.6)
    return 3;
  else if (-0.6 < rapidity && rapidity <= -0.5)
    return 4;
  else if (-0.5 < rapidity && rapidity <= -0.4)
    return 5;
  else if (-0.4 < rapidity && rapidity <= -0.3)
    return 6;
  else if (-0.3 < rapidity && rapidity <= -0.2)
    return 7;
  else if (-0.2 < rapidity && rapidity <= -0.1)
    return 8;
  else if (-0.1 < rapidity && rapidity <= 0.f)
    return 9;
  else if (0.f < rapidity && rapidity <= 0.1)
    return 10;
  else if (0.1 < rapidity && rapidity <= 0.2)
    return 11;
  else if (0.2 < rapidity && rapidity <= 0.3)
    return 12;
  else if (0.3 < rapidity && rapidity <= 0.4)
    return 13;
  else if (0.4 < rapidity && rapidity <= 0.5)
    return 14;
  else if (0.5 < rapidity && rapidity <= 0.6)
    return 15;
  else if (0.6 < rapidity && rapidity <= 0.7)
    return 16;
  else if (0.7 < rapidity && rapidity <= 0.8)
    return 17;
  else if (0.8 < rapidity && rapidity <= 0.9)
    return 18;
  else if (0.9 < rapidity && rapidity < 1.f)
    return 19;
  else
    return -1;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::InitCutHistograms(TString appendix) {
  std::cout << "============================\n"
            << " V0 CUT CONFIGURATION \n"
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
            << "============================\n";

  std::vector<float> rapBins = {{-1.,  -0.9, -0.8, -0.7, -0.6, -0.5, -0.4,
                                 -0.3, -0.2, -0.1, 0.f,  0.1,  0.2,  0.3,
                                 0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.f}};

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

  fHistCutBooking = new TProfile("fHistCutBooking", ";;Cut value", 25, 0, 25);
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
  fHistCutBooking->GetXaxis()->SetBinLabel(15, "Daughter DCA max");
  fHistCutBooking->GetXaxis()->SetBinLabel(16, "Daughter DCA pvtx min");
  fHistCutBooking->GetXaxis()->SetBinLabel(17, "K^{0} rejection low");
  fHistCutBooking->GetXaxis()->SetBinLabel(18, "K^{0} rejection up");
  fHistCutBooking->GetXaxis()->SetBinLabel(19, "#Lambda selection low");
  fHistCutBooking->GetXaxis()->SetBinLabel(20, "#Lambda selection up");
  fHistCutBooking->GetXaxis()->SetBinLabel(21, "Pile-up rejection");
  fHistCutBooking->GetXaxis()->SetBinLabel(22, "Armenteros q_{T} low");
  fHistCutBooking->GetXaxis()->SetBinLabel(23, "Armenteros q_{T} up");
  fHistCutBooking->GetXaxis()->SetBinLabel(24, "Armenteros #alpha low");
  fHistCutBooking->GetXaxis()->SetBinLabel(25, "Armenteros #alpha up");
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
  fHistCutBooking->Fill(14.f, fDaughterDCAMax);
  fHistCutBooking->Fill(15.f, fDaughterDCAPV);
  fHistCutBooking->Fill(16.f, fK0RejectionLow);
  fHistCutBooking->Fill(17.f, fK0RejectionUp);
  fHistCutBooking->Fill(18.f, fLambdaSelectionLow);
  fHistCutBooking->Fill(19.f, fLambdaSelectionUp);
  fHistCutBooking->Fill(20.f, static_cast<double>(fPileUpRejectionMode));
  fHistCutBooking->Fill(21.f, fArmenterosQtLow);
  fHistCutBooking->Fill(22.f, fArmenterosQtUp);
  fHistCutBooking->Fill(23.f, fArmenterosAlphaLow);
  fHistCutBooking->Fill(24.f, fArmenterosAlphaUp);

  fHistCuts = new TH1F("fHistCuts", ";;Entries", 18, 0, 18);
  fHistCuts->GetXaxis()->SetBinLabel(1, "V0");
  fHistCuts->GetXaxis()->SetBinLabel(2, "Neutral");
  fHistCuts->GetXaxis()->SetBinLabel(3, "On-fly");
  fHistCuts->GetXaxis()->SetBinLabel(4, "#it{p}_{T, min}");
  fHistCuts->GetXaxis()->SetBinLabel(5, "#it{p}_{T, max}");
  fHistCuts->GetXaxis()->SetBinLabel(6, "Daughter tracks");
  fHistCuts->GetXaxis()->SetBinLabel(7, "PID");
  fHistCuts->GetXaxis()->SetBinLabel(8, "SingleParticle QA");
  fHistCuts->GetXaxis()->SetBinLabel(9, "Pile-up");
  fHistCuts->GetXaxis()->SetBinLabel(10, "Decay vtx x");
  fHistCuts->GetXaxis()->SetBinLabel(11, "Decay vtx y");
  fHistCuts->GetXaxis()->SetBinLabel(12, "Decay vtx z");
  fHistCuts->GetXaxis()->SetBinLabel(13, "Radius min");
  fHistCuts->GetXaxis()->SetBinLabel(14, "Radius max");
  fHistCuts->GetXaxis()->SetBinLabel(15, "Daughter DCA");
  fHistCuts->GetXaxis()->SetBinLabel(16, "cos #alpha");
  fHistCuts->GetXaxis()->SetBinLabel(17, "#Lambda selection");
  fHistograms->Add(fHistCuts);

  fHistNV0 =
      new TH1F("fHistNV0", ";Number of V0 candidates; Entries", 25, 0, 25);
  fHistograms->Add(fHistNV0);

  fHistLambdaMassPt = new TH2F("fHistLambdaMassPt",
                               "; #it{p}_{T} (GeV/#it{c});Invariant mass "
                               "p#pi hypothesis (GeV/#it{c}^{2})",
                               1000, 0, 15, 400, 1., 1.2);
  fHistograms->Add(fHistLambdaMassPt);

  if (!fIsLightweight) {
    fHistLambdaMass =
        new TH1F("fHistLambdaMass",
                 "; Invariant mass p#pi hypothesis (GeV/#it{c}^{2}); Entries",
                 4000, 0., 2.);
    fHistLambdaPt = new TH1F("fHistLambdaPt",
                             "; #it{p}_{T} (GeV/#it{c}); Entries", 500, 0, 10);
    fHistLambdaMassK0Rej = new TH1F("fHistLambdaMassK0Rej",
                                    "; Invariant mass p#pi hypothesis w/o "
                                    "K^{0} rejection (GeV/#it{c}^{2}); Entries",
                                    1000, 0., 2.);

    fHistK0Mass =
        new TH1F("fHistK0Mass",
                 "; Invariant mass #pi#pi hypothesis (GeV/#it{c}^{2}); Entries",
                 1000, 0., 2.);

    fHistK0MassAfter =
        new TH1F("fHistK0MassAfter",
                 "; Invariant mass #pi#pi hypothesis (GeV/#it{c}^{2}); Entries",
                 1000, 0., 2.);
    fHistCosPA =
        new TH2F("fHistCosPA", "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)", 1000,
                 0, 10, 1000, 0.9, 1);
    fHistEtaPhi = new TH2F("fHistEtaPhi", "; #eta; #phi", 200, -1, 1, 200, 0,
                           2 * TMath::Pi());

    fHistograms->Add(fHistLambdaMass);
    fHistograms->Add(fHistLambdaPt);
    fHistograms->Add(fHistLambdaMassK0Rej);
    fHistograms->Add(fHistK0Mass);
    fHistograms->Add(fHistK0MassAfter);
    fHistograms->Add(fHistCosPA);
    fHistograms->Add(fHistEtaPhi);

    for (int i = 0; i < static_cast<int>(rapBins.size() - 1); ++i) {
      fHistLambdaPtY[i] =
          new TH2F(Form("fHistLambdaPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
                   Form("%.2f < y < %.2f ; #it{p}_{T} (GeV/#it{c}); Invariant "
                        "mass p#pi hypothesis (GeV/#it{c}^{2})",
                        rapBins[i], rapBins[i + 1]),
                   500, 0, 10, 400, 1., 1.2);
      fHistograms->Add(fHistLambdaPtY[i]);
    }
  }

  fHistSingleParticleCuts[0] =
      new TH1F("fHistSingleParticleCuts_pos", ";;Entries", 10, 0, 10);
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(2, "#eta");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(3, "nCls TPC");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(4, "nCrossed Rows TPC");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(5, "Ratio Findable TPC");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(6, "nCls Findable");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(7, "nCls Shared max");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(8, "Daughter DCA to PV");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(9, "TPC refit");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(10, "Kink");
  fHistogramsPos->Add(fHistSingleParticleCuts[0]);

  fHistSingleParticleCuts[1] =
      new TH1F("fHistSingleParticleCuts_neg", ";;Entries", 10, 0, 10);
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(2, "#eta");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(3, "nCls TPC");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(4, "nCrossed Rows TPC");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(5, "Ratio Findable TPC");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(6, "nCls Findable");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(7, "nCls Shared max");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(8, "Daughter DCA to PV");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(9, "TPC refit");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(10, "Kink");
  fHistogramsNeg->Add(fHistSingleParticleCuts[1]);

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

    fHistMCTruthV0Pt = new TH1F(
        "fHistMCTruthV0Pt", "; #it{p}_{T} (GeV/#it{c}); Entries", 500, 0, 10);
    fHistMCTruthV0PtY =
        new TH2F("fHistMCTruthV0PtY", "; y; #it{p}_{T} (GeV/#it{c})", 500, -10,
                 10, 500, 0, 10);
    fHistMCTruthV0PtEta =
        new TH2F("fHistMCTruthV0PtEta", "; #eta; #it{p}_{T} (GeV/#it{c})", 500,
                 -10, 10, 500, 0, 10);
    fHistMCTruthV0ProtonPionPt =
        new TH1F("fHistMCTruthV0ProtonPionPt",
                 "; #it{p}_{T} (GeV/#it{c}); Entries", 500, 0, 10);
    fHistMCTruthV0ProtonPionPtY =
        new TH2F("fHistMCTruthV0ProtonPionPtY", "; y; #it{p}_{T} (GeV/#it{c})",
                 500, -10, 10, 500, 0, 10);
    fHistMCTruthV0ProtonPionPtEta =
        new TH2F("fHistMCTruthV0ProtonPionPtEta",
                 "; #eta; #it{p}_{T} (GeV/#it{c})", 500, -10, 10, 500, 0, 10);
    fHistogramsMC->Add(fHistMCTruthV0Pt);
    fHistogramsMC->Add(fHistMCTruthV0PtY);
    fHistogramsMC->Add(fHistMCTruthV0PtEta);
    fHistogramsMC->Add(fHistMCTruthV0ProtonPionPt);
    fHistogramsMC->Add(fHistMCTruthV0ProtonPionPtY);
    fHistogramsMC->Add(fHistMCTruthV0ProtonPionPtEta);
    fHistograms->Add(fHistogramsMC);
  }

  if (!fIsLightweight) {
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
        new TH1F("fHistDecayVertexXBefore", "; Decay vertex x (cm); Entries",
                 1000, 0, 200);
    fHistDecayVertexYBefore =
        new TH1F("fHistDecayVertexYBefore", "; Decay vertex y (cm); Entries",
                 1000, 0, 200);
    fHistDecayVertexZBefore =
        new TH1F("fHistDecayVertexZBefore", "; Decay vertex z (cm); Entries",
                 1000, 0, 200);
    fHistDecayVertexXAfter =
        new TH1F("fHistDecayVertexXAfter", "; Decay vertex x (cm); Entries",
                 1000, 0, 200);
    fHistDecayVertexYAfter =
        new TH1F("fHistDecayVertexYAfter", "; Decay vertex y (cm); Entries",
                 1000, 0, 200);
    fHistDecayVertexZAfter =
        new TH1F("fHistDecayVertexZAfter", "; Decay vertex z (cm); Entries",
                 1000, 0, 200);
    fHistTransverseRadiusBefore =
        new TH1F("fHistTransverseRadiusBefore",
                 "; Transverse radius (cm); Entries", 1000, 0, 200);
    fHistTransverseRadiusAfter =
        new TH1F("fHistTransverseRadiusAfter",
                 "; Transverse radius (cm); Entries", 1000, 0, 200);
    fHistCosPABefore =
        new TH1F("fHistCosPABefore", "; cos(#alpha); Entries", 1000, 0.5, 1);
    fHistCosPAAfter =
        new TH1F("fHistCosPAAfter", "; cos(#alpha); Entries", 1000, 0.95, 1);
    fHistDCADaughtersBefore =
        new TH1F("fHistDCADaughtersBefore",
                 "; Daughter DCA at decay vertex (cm); Entries", 1000, 0, 10);
    fHistDCADaughtersAfter =
        new TH1F("fHistDCADaughtersAfter",
                 "; Daughter DCA at decay vertex (cm); Entries", 1000, 0, 10);
    fHistDCA = new TH1F("fHistDCA", "; DCA to PV (cm); Entries", 1000, 0, 10);
    fHistDecayLength = new TH1F("fHistDecayLength",
                                "; Decay length (cm); Entries", 1000, 0, 200);
    fHistArmenterosBefore = new TH2F("fHistArmenterosBefore",
                                     " ; #alpha; #it{q}_{T} p#pi (GeV/#it{c})",
                                     500, -1, 1, 500, 0, 1);
    fHistArmenterosAfter = new TH2F("fHistArmenterosAfter",
                                    " ; #alpha; #it{q}_{T} p#pi (GeV/#it{c})",
                                    500, -1, 1, 500, 0, 1);

    fHistogramsBefore->Add(fHistDecayVertexXBefore);
    fHistogramsBefore->Add(fHistDecayVertexYBefore);
    fHistogramsBefore->Add(fHistDecayVertexZBefore);
    fHistogramsBefore->Add(fHistTransverseRadiusBefore);
    fHistogramsBefore->Add(fHistCosPABefore);
    fHistogramsBefore->Add(fHistDCADaughtersBefore);
    fHistogramsAfter->Add(fHistDecayVertexXAfter);
    fHistogramsAfter->Add(fHistDecayVertexYAfter);
    fHistogramsAfter->Add(fHistDecayVertexZAfter);
    fHistogramsAfter->Add(fHistTransverseRadiusAfter);
    fHistogramsAfter->Add(fHistCosPAAfter);
    fHistogramsAfter->Add(fHistDCADaughtersAfter);
    fHistograms->Add(fHistDCA);
    fHistograms->Add(fHistDecayLength);
    fHistogramsBefore->Add(fHistArmenterosBefore);
    fHistogramsAfter->Add(fHistArmenterosAfter);

    fHistSingleParticlePt[0] = new TH1F("fHistSingleParticlePt_pos",
                                        ";#it{p}_{T}; Entries", 1000, 0, 10);
    fHistSingleParticlePt[1] = new TH1F("fHistSingleParticlePt_neg",
                                        ";#it{p}_{T}; Entries", 1000, 0, 10);
    fHistSingleParticleEtaBefore[0] =
        new TH1F("fHistSingleParticleEtaBefore_pos", "; #eta (before); Entries",
                 1000, -2, 2);
    fHistSingleParticleEtaBefore[1] =
        new TH1F("fHistSingleParticleEtaBefore_neg", "; #eta (before); Entries",
                 1000, -2, 2);
    fHistSingleParticleEtaAfter[0] =
        new TH1F("fHistSingleParticleEtaAfter_pos", "; #eta (after); Entries",
                 1000, -2, 2);
    fHistSingleParticleEtaAfter[1] =
        new TH1F("fHistSingleParticleEtaAfter_neg", "; #eta (after); Entries",
                 1000, -2, 2);
    fHistSingleParticleNclsTPCBefore[0] =
        new TH1F("fHistSingleParticleNclsTPCBefore_pos",
                 "; # cls TPC (before); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCBefore[1] =
        new TH1F("fHistSingleParticleNclsTPCBefore_neg",
                 "; # cls TPC (before); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCAfter[0] =
        new TH1F("fHistSingleParticleNclsTPCAfter_pos",
                 "; # cls TPC (after); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCAfter[1] =
        new TH1F("fHistSingleParticleNclsTPCAfter_neg",
                 "; # cls TPC (after); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCFindableBefore[0] =
        new TH1F("fHistSingleParticleNclsTPCFindableBefore_pos",
                 "; # cls TPC findable (before); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCFindableBefore[1] =
        new TH1F("fHistSingleParticleNclsTPCFindableBefore_neg",
                 "; # cls TPC findable (before); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCFindableAfter[0] =
        new TH1F("fHistSingleParticleNclsTPCFindableAfter_pos",
                 "; # cls TPC findable (after); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCFindableAfter[1] =
        new TH1F("fHistSingleParticleNclsTPCFindableAfter_neg",
                 "; # cls TPC findable (after); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCRatioFindableBefore[0] =
        new TH1F("fHistSingleParticleNclsTPCRatioFindableBefore_pos",
                 "; TPC ratio findable (before); Entries", 1000, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableBefore[1] =
        new TH1F("fHistSingleParticleNclsTPCRatioFindableBefore_neg",
                 ";  TPC ratio findable (before); Entries", 1000, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableAfter[0] =
        new TH1F("fHistSingleParticleNclsTPCRatioFindableAfter_pos",
                 "; TPC ratio findable (after); Entries", 1000, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableAfter[1] =
        new TH1F("fHistSingleParticleNclsTPCRatioFindableAfter_neg",
                 ";  TPC ratio findable (after); Entries", 1000, 0, 2);
    fHistSingleParticleNcrossedTPCBefore[0] =
        new TH1F("fHistSingleParticleNcrossedTPCBefore_pos",
                 "; # cls TPC crossed (before); Entries", 170, 0, 170);
    fHistSingleParticleNcrossedTPCBefore[1] =
        new TH1F("fHistSingleParticleNcrossedTPCBefore_neg",
                 "; # cls TPC crossed (before); Entries", 170, 0, 170);
    fHistSingleParticleNcrossedTPCAfter[0] =
        new TH1F("fHistSingleParticleNcrossedTPCAfter_pos",
                 "; # cls TPC crossed (after); Entries", 170, 0, 170);
    fHistSingleParticleNcrossedTPCAfter[1] =
        new TH1F("fHistSingleParticleNcrossedTPCAfter_neg",
                 "; # cls TPC crossed (after); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCShared[0] =
        new TH1F("fHistSingleParticleNclsTPCShared_pos",
                 "; # cls TPC shared; Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCShared[1] =
        new TH1F("fHistSingleParticleNclsTPCShared_neg",
                 "; # cls TPC shared; Entries", 170, 0, 170);
    fHistSingleParticleNclsITSShared[0] =
        new TH1F("fHistSingleParticleNclsITSShared_pos",
                 "; # cls ITS shared; Entries", 6, 0, 6);
    fHistSingleParticleNclsITSShared[1] =
        new TH1F("fHistSingleParticleNclsITSShared_neg",
                 "; # cls ITS shared; Entries", 6, 0, 6);
    fHistSingleParticleDCAtoPVBefore[0] =
        new TH1F("fHistSingleParticleDCAtoPVBefore_pos",
                 "; DCA (pos, PA) (before); Entries", 1000, 0, 10);
    fHistSingleParticleDCAtoPVBefore[1] =
        new TH1F("fHistSingleParticleDCAtoPVBefore_neg",
                 "; DCA (pos, PA) (before); Entries", 1000, 0, 10);
    fHistSingleParticleDCAtoPVAfter[0] =
        new TH1F("fHistSingleParticleDCAtoPVAfter_pos",
                 "; DCA (pos, PA) (after); Entries", 1000, 0, 10);
    fHistSingleParticleDCAtoPVAfter[1] =
        new TH1F("fHistSingleParticleDCAtoPVAfter_neg",
                 "; DCA (pos, PA) (after); Entries", 1000, 0, 10);
    fHistSingleParticlePileUp[0] =
        new TH2F("fHistSingleParticlePileUp_pos",
                 "; Pileup flag; #it{p}_{T} (GeV/#it{c})", 4, 0, 4, 100, 0, 10);
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(1, "ITS");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(2, "TOF");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(3, "Combined");
    fHistSingleParticlePileUp[0]->GetXaxis()->SetBinLabel(4, "None");
    fHistSingleParticlePileUp[1] =
        new TH2F("fHistSingleParticlePileUp_neg",
                 "; Pileup flag; #it{p}_{T} (GeV/#it{c})", 4, 0, 4, 100, 0, 10);
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(1, "ITS");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(2, "TOF");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(3, "Combined");
    fHistSingleParticlePileUp[1]->GetXaxis()->SetBinLabel(4, "None");

    fHistogramsPos->Add(fHistSingleParticlePt[0]);
    fHistogramsNeg->Add(fHistSingleParticlePt[1]);
    fHistogramsPos->Add(fHistSingleParticleEtaBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleEtaBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCFindableBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCFindableBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCRatioFindableBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCRatioFindableBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNcrossedTPCBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleNcrossedTPCBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCShared[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCShared[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsITSShared[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsITSShared[1]);
    fHistogramsPos->Add(fHistSingleParticleDCAtoPVBefore[0]);
    fHistogramsNeg->Add(fHistSingleParticleDCAtoPVBefore[1]);
    fHistogramsPos->Add(fHistSingleParticleEtaAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleEtaAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCFindableAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCFindableAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsTPCRatioFindableAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCRatioFindableAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleNcrossedTPCAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleNcrossedTPCAfter[1]);
    fHistogramsPos->Add(fHistSingleParticleDCAtoPVAfter[0]);
    fHistogramsNeg->Add(fHistSingleParticleDCAtoPVAfter[1]);
    fHistogramsPos->Add(fHistSingleParticlePileUp[0]);
    fHistogramsNeg->Add(fHistSingleParticlePileUp[1]);

    fHistograms->Add(fHistogramsBefore);
    fHistograms->Add(fHistogramsAfter);

    fHistSingleParticlePID[0] = new TH2F(
        "fHistSingleParticlePID_pos", "; #it{p} (GeV/#it{c}); n_{#sigma} TPC",
        1000, 0, 5, 1000, -10, 10);
    fHistSingleParticlePID[1] = new TH2F(
        "fHistSingleParticlePID_neg", "; #it{p} (GeV/#it{c}); n_{#sigma} TPC",
        1000, 0, 5, 1000, -10, 10);
    fHistogramsPos->Add(fHistSingleParticlePID[0]);
    fHistogramsNeg->Add(fHistSingleParticlePID[1]);
  }
  fHistograms->Add(fHistogramsPos);
  fHistograms->Add(fHistogramsNeg);
}
