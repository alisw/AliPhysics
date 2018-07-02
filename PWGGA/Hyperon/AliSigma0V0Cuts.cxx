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
      fDataBasePDG(),
      fIsLightweight(false),
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
      fUsePID(true),
      fV0PtMin(0.f),
      fV0PtMax(1E30),
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
      fPIDResponse(nullptr),
      fHistCutBooking(nullptr),
      fHistCuts(nullptr),
      fHistNV0(nullptr),
      fHistLambdaMass(nullptr),
      fHistAntiLambdaMass(nullptr),
      fHistPhotonMass(nullptr),
      fHistK0Mass(nullptr),
      fHistV0Pt(nullptr),
      fHistV0Mass(nullptr),
      fHistV0PtY(),
      fHistV0MassPt(nullptr),
      fHistLambdaMassK0Rej(nullptr),
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
      fHistMCTruthV0Pt(nullptr),
      fHistMCTruthV0PtY(nullptr),
      fHistMCTruthV0PtEta(nullptr),
      fHistMCTruthV0DaughterPt(nullptr),
      fHistMCTruthV0DaughterPtY(nullptr),
      fHistMCTruthV0DaughterPtEta(nullptr),
      fHistMCV0Pt(nullptr),
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
      fDataBasePDG(),
      fIsLightweight(false),
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
      fUsePID(true),
      fV0PtMin(0.f),
      fV0PtMax(1E30),
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
      fPIDResponse(nullptr),
      fHistCutBooking(nullptr),
      fHistCuts(nullptr),
      fHistNV0(nullptr),
      fHistLambdaMass(nullptr),
      fHistAntiLambdaMass(nullptr),
      fHistPhotonMass(nullptr),
      fHistK0Mass(nullptr),
      fHistV0Pt(nullptr),
      fHistV0Mass(nullptr),
      fHistV0PtY(),
      fHistV0MassPt(nullptr),
      fHistLambdaMassK0Rej(nullptr),
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
      fHistMCTruthV0Pt(nullptr),
      fHistMCTruthV0PtY(nullptr),
      fHistMCTruthV0PtEta(nullptr),
      fHistMCTruthV0DaughterPt(nullptr),
      fHistMCTruthV0DaughterPtY(nullptr),
      fHistMCTruthV0DaughterPtEta(nullptr),
      fHistMCV0Pt(nullptr),
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
  v0Cuts->SetV0PtMin(0.1);
  v0Cuts->SetV0CosPAMin(0.99);
  v0Cuts->SetV0RadiusMax(200.f);
  v0Cuts->SetV0RadiusMin(5.f);
  v0Cuts->SetPIDnSigma(5.f);
  v0Cuts->SetTPCRatioFindable(0.6f);
  v0Cuts->SetEtaMax(0.8);
  v0Cuts->SetDaughterDCAMax(1.5);
  v0Cuts->SetDaughterDCAtoPV(0.05);
  v0Cuts->SetArmenterosCut(0, 0.06, -1, 1);
  v0Cuts->SetK0Rejection(0., 0.);
  v0Cuts->SetLambdaSelection(0., 0.);
  v0Cuts->SetPsiPairMax(0.2);
  v0Cuts->SetPileUpRejectionMode(None);
  return v0Cuts;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::SelectV0(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                               std::vector<AliSigma0ParticleV0> &V0Container) {
  // Get the PID manager from the input event handler
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      static_cast<AliInputEventHandler *>(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fMCEvent = mcEvent;

  V0Container.clear();

  if (fIsMC) {
    ProcessMC();
  }

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

    // Single Particle Quality
    if (!SingleParticleQualityCuts(pos) || !SingleParticleQualityCuts(neg))
      continue;
    fHistCuts->Fill(6);

    // Pile-up rejection
    if (!PileUpRejection(pos, neg)) continue;
    fHistCuts->Fill(7);

    // Topological Selection
    if (!V0TopologicalSelection(v0)) continue;

    // PID
    if (!V0PID(v0, pos, neg)) continue;

    PlotMasses(v0);

    // V0 Selection
    if (std::abs(fPID) == 3122 && !LambdaSelection(v0)) continue;
    if (std::abs(fPID) == 22 && !PhotonSelection(v0)) continue;

    AliSigma0ParticleV0 v0Candidate(v0, pos, neg,
                                    fInputEvent->GetPrimaryVertex(), fPID,
                                    fInputEvent->GetMagneticField(), fMCEvent);

    if (fIsMC) {
      int label = v0Candidate.MatchToMC(fMCEvent, fPID, {{fPosPDG, fNegPDG}});
      if (label < 0) {
        // no mc info assigned to this track - don't use it
        v0Candidate.SetUse(false);
      }

      AliMCParticle *mcParticle =
          static_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
      if (!mcParticle) continue;
      v0Candidate.ProcessMCInfo(mcParticle, fMCEvent);

      fHistMCV0Pt->Fill(v0->Pt());
    }

    v0Candidate.SetPDGMass(fDataBasePDG.GetParticle(fPID)->Mass());
    if (std::abs(fPID) == 22) {
      const float massPhoton = ComputePhotonMass(v0);
      v0Candidate.SetMass(massPhoton);
      v0Candidate.SetRecMass(massPhoton);
    }

    V0Container.push_back(v0Candidate);
  }
  fHistNV0->Fill(V0Container.size());
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0QualityCuts(const AliESDv0 *v0) const {
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

  fHistCuts->Fill(15);
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
  const float ratioFindable =
      nCrossedRows / (static_cast<float>(nFindable) + 1.e-13);
  const float chi2 =
      (nClsTPC > 5) ? (track->GetTPCchi2() / float(nClsTPC - 5)) : -1.;
  /// see AliAnalysisTaskESDfilter::Chi2perNDF

  const float magField = fInputEvent->GetMagneticField();
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  const float dcaDaughterToPV =
      std::abs(track->GetD(vertex->GetX(), vertex->GetY(), magField));

  const int histoPrefix = (charge > 0) ? 0 : 1;
  int qaHistoCounter = 0;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (!fIsLightweight) {
    fHistSingleParticleEtaBefore[histoPrefix]->Fill(pt, eta);
    fHistSingleParticleChi2Before[histoPrefix]->Fill(pt, chi2);
    fHistSingleParticleNclsTPCBefore[histoPrefix]->Fill(pt, nClsTPC);
    fHistSingleParticleNclsTPCFindableBefore[histoPrefix]->Fill(pt, nFindable);
    fHistSingleParticleNclsTPCRatioFindableBefore[histoPrefix]->Fill(
        pt, ratioFindable);
    fHistSingleParticleNcrossedTPCBefore[histoPrefix]->Fill(pt, nCrossedRows);
    fHistSingleParticleDCAtoPVBefore[histoPrefix]->Fill(pt, dcaDaughterToPV);
  }

  // Max eta cut
  if (std::abs(eta) > fEtaMax) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  // Max chi2 cut
  if (chi2 > fChi2Max) return false;
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
    fHistSingleParticleNclsTPCShared[histoPrefix]->Fill(pt, nClsSharedTPC);
    fHistSingleParticleNclsITSShared[histoPrefix]->Fill(pt, nClsSharedITS);
  }

  if (!(track->GetStatus() & AliESDtrack::kTPCrefit)) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (track->GetKinkIndex(0) > 0) return false;
  fHistSingleParticleCuts[histoPrefix]->Fill(qaHistoCounter++);

  if (!fIsLightweight) {
    fHistSingleParticlePt[histoPrefix]->Fill(pt);
    fHistSingleParticleEtaAfter[histoPrefix]->Fill(pt, eta);
    fHistSingleParticleChi2After[histoPrefix]->Fill(pt, chi2);
    fHistSingleParticleNclsTPCAfter[histoPrefix]->Fill(pt, nClsTPC);
    fHistSingleParticleNclsTPCFindableAfter[histoPrefix]->Fill(pt, nFindable);
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
  if (std::abs(decayVertexV0[0]) > fV0DecayVertexMax) return false;
  fHistCuts->Fill(8);
  if (std::abs(decayVertexV0[1]) > fV0DecayVertexMax) return false;
  fHistCuts->Fill(9);
  if (std::abs(decayVertexV0[2]) > fV0DecayVertexMax) return false;
  fHistCuts->Fill(10);

  // Transverse decay radius min & max
  if (transverseRadius < fV0RadiusMin) return false;
  fHistCuts->Fill(11);

  if (transverseRadius > fV0RadiusMax) return false;
  fHistCuts->Fill(12);

  // DCA of the daughter tracks at the decay vertex
  if (dcaV0Dau > fDaughterDCAMax) return false;
  fHistCuts->Fill(13);

  // Cos pointing angle
  if (point < fV0CosPAMin) return false;
  fHistCuts->Fill(14);

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

  if (!fIsLightweight) {
    fHistK0Mass->Fill(massK0);
    fHistLambdaMass->Fill(massLambda);
    fHistAntiLambdaMass->Fill(massAntiLambda);
    fHistPhotonMass->Fill(massPhoton);
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
  fHistCuts->Fill(16);

  float rap = v0->RapLambda();
  int rapBin = GetRapidityBin(rap);

  fHistV0MassPt->Fill(v0->Pt(), massLambda);
  if (!fIsLightweight) {
    fHistLambdaMassK0Rej->Fill(massLambda);
    if (rapBin > -1) fHistV0PtY[rapBin]->Fill(v0->Pt(), massLambda);
    fHistK0MassAfter->Fill(massK0);
    fHistEtaPhi->Fill(v0->Eta(), v0->Phi());
  }

  if (!fIsLightweight) fHistV0Mass->Fill(massLambda);

  // Lambda selection cut
  if (massLambda < fLambdaSelectionLow || massLambda > fLambdaSelectionUp) {
    return false;
  }

  if (!fIsLightweight) fHistV0Pt->Fill(v0->Pt());
  fHistCuts->Fill(17);
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::PhotonSelection(AliESDv0 *v0) const {
  v0->ChangeMassHypothesis(3122);
  const float massLambda = v0->GetEffMass();
  v0->ChangeMassHypothesis(310);
  const float massK0 = v0->GetEffMass();
  const float massPhoton = ComputePhotonMass(v0);

  const float psiPair = ComputePsiPair(v0);
  if (!fIsLightweight) fHistPsiPair->Fill(v0->Pt(), psiPair);

  if (std::abs(psiPair) > fPsiPairMax) return false;
  fHistCuts->Fill(16);

  // K0 rejection cut
  if (fK0Rejection && massK0 > fK0RejectionLow && massK0 < fK0RejectionUp) {
    return false;
  }
  fHistCuts->Fill(17);

  float rap = ComputeRapidity(v0->Pt(), v0->Pz(), massPhoton);
  int rapBin = GetRapidityBin(rap);

  fHistV0MassPt->Fill(v0->Pt(), massPhoton);
  if (!fIsLightweight) {
    if (rapBin > -1) fHistV0PtY[rapBin]->Fill(v0->Pt(), massPhoton);
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

  fHistCuts->Fill(18);
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
  // Loop over the MC tracks
  for (int iPart = 1; iPart < (fMCEvent->GetNumberOfTracks()); iPart++) {
    AliMCParticle *mcParticle =
        static_cast<AliMCParticle *>(fMCEvent->GetTrack(iPart));
    if (!mcParticle) continue;
    if (!mcParticle->IsPhysicalPrimary()) continue;
    if (mcParticle->PdgCode() != fPID) continue;
    fHistMCTruthV0Pt->Fill(mcParticle->Pt());
    fHistMCTruthV0PtY->Fill(
        ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(), mcParticle->M()),
        mcParticle->Pt());
    fHistMCTruthV0PtEta->Fill(mcParticle->Eta(), mcParticle->Pt());

    if (!CheckDaughters(mcParticle)) continue;
    fHistMCTruthV0DaughterPt->Fill(mcParticle->Pt());
    fHistMCTruthV0DaughterPtY->Fill(
        ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(), mcParticle->M()),
        mcParticle->Pt());
    fHistMCTruthV0DaughterPtEta->Fill(mcParticle->Eta(), mcParticle->Pt());
  }
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::CheckDaughters(const AliMCParticle *particle) const {
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

  fHistCutBooking = new TProfile("fHistCutBooking", ";;Cut value", 31, 0, 31);
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
  fHistCutBooking->Fill(13.f, fChi2Max);
  fHistCutBooking->Fill(14.f, fDaughterDCAMax);
  fHistCutBooking->Fill(15.f, fDaughterDCAPV);
  fHistCutBooking->Fill(16.f, fK0RejectionLow);
  fHistCutBooking->Fill(17.f, fK0RejectionUp);
  fHistCutBooking->Fill(18.f, fLambdaSelectionLow);
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
  if (std::abs(fPID) == 3122) {
    fHistCuts->GetXaxis()->SetBinLabel(17, "K^{0} rejection");
    fHistCuts->GetXaxis()->SetBinLabel(18, "#Lambda selection");
  } else if (std::abs(fPID) == 22) {
    fHistCuts->GetXaxis()->SetBinLabel(17, "#Psi_{pair} selection");
    fHistCuts->GetXaxis()->SetBinLabel(18, "K^{0} rejection");
    fHistCuts->GetXaxis()->SetBinLabel(19, "#Lambda rejection");
  }
  fHistograms->Add(fHistCuts);

  fHistNV0 =
      new TH1F("fHistNV0", ";Number of V0 candidates; Entries", 25, 0, 25);
  fHistograms->Add(fHistNV0);

  fHistV0MassPt =
      new TH2F("fHistV0MassPt",
               "; #it{p}_{T} (GeV/#it{c});Invariant mass (GeV/#it{c}^{2})",
               1000, 0, 15, 2400, 0., 1.2);
  fHistograms->Add(fHistV0MassPt);

  if (!fIsLightweight) {
    fHistLambdaMass = new TH1F(
        "fHistLambdaMass",
        "; Invariant mass p#pi^{-} hypothesis (GeV/#it{c}^{2}); Entries", 2000,
        0., 2.);
    fHistograms->Add(fHistLambdaMass);

    fHistAntiLambdaMass = new TH1F(
        "fHistAntiLambdaMass",
        "; Invariant mass #bar{p}#pi hypothesis (GeV/#it{c}^{2}); Entries",
        2000, 0., 2.);
    fHistograms->Add(fHistAntiLambdaMass);

    fHistPhotonMass = new TH1F(
        "fHistPhotonMass",
        "; Invariant mass e^{+}e^{-} hypothesis (GeV/#it{c}^{2}); Entries",
        2000, 0., 2.);
    fHistograms->Add(fHistPhotonMass);

    fHistK0Mass =
        new TH1F("fHistK0Mass",
                 "; Invariant mass #pi#pi hypothesis (GeV/#it{c}^{2}); Entries",
                 2000, 0., 2.);
    fHistograms->Add(fHistK0Mass);

    fHistV0Pt =
        new TH1F("fHistV0Pt", "; #it{p}_{T} (GeV/#it{c}); Entries", 500, 0, 10);
    fHistograms->Add(fHistV0Pt);

    fHistV0Mass =
        new TH1F("fHistV0Mass", "; Invariant mass (GeV/#it{c}^{2}); Entries", 2400, 0., 1.2);
    fHistograms->Add(fHistV0Mass);

    fHistLambdaMassK0Rej = new TH1F("fHistLambdaMassK0Rej",
                                    "; Invariant mass p#pi hypothesis w/o "
                                    "K^{0} rejection (GeV/#it{c}^{2}); Entries",
                                    1000, 0., 2.);
    fHistograms->Add(fHistLambdaMassK0Rej);

    fHistK0MassAfter =
        new TH1F("fHistK0MassAfter",
                 "; Invariant mass #pi#pi hypothesis (GeV/#it{c}^{2}); Entries",
                 1000, 0., 2.);
    fHistograms->Add(fHistK0MassAfter);
    fHistCosPA =
        new TH2F("fHistCosPA", "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)", 1000,
                 0, 10, 1000, 0.9, 1);
    fHistograms->Add(fHistCosPA);

    fHistEtaPhi = new TH2F("fHistEtaPhi", "; #eta; #phi", 200, -1, 1, 200, 0,
                           2 * TMath::Pi());
    fHistograms->Add(fHistEtaPhi);

    if (fPID == 22) {
      fHistPsiPair =
          new TH2F("fHistPsiPair", "; #it{p}_{T} (GeV/#it{c}); #Psi_{pair}",
                   1000, 0, 10, 1000, TMath::Pi(), -TMath::Pi());
      fHistograms->Add(fHistPsiPair);
    }

    for (int i = 0; i < static_cast<int>(rapBins.size() - 1); ++i) {
      fHistV0PtY[i] =
          new TH2F(Form("fHistV0PtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
                   Form("%.2f < y < %.2f ; #it{p}_{T} (GeV/#it{c}); Invariant "
                        "mass p#pi hypothesis (GeV/#it{c}^{2})",
                        rapBins[i], rapBins[i + 1]),
                   500, 0, 10, 400, 1., 1.2);
      fHistograms->Add(fHistV0PtY[i]);
    }
  }

  fHistSingleParticleCuts[0] =
      new TH1F("fHistSingleParticleCuts_pos", ";;Entries", 11, 0, 11);
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(2, "#eta");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(3, "#chi^{2}");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(4, "nCls TPC");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(5, "nCrossed Rows TPC");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(6, "Ratio Findable TPC");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(7, "nCls Findable");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(8, "nCls Shared max");
  fHistSingleParticleCuts[0]->GetXaxis()->SetBinLabel(9, "Daughter DCA to PV");
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
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(6, "Ratio Findable TPC");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(7, "nCls Findable");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(8, "nCls Shared max");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(9, "Daughter DCA to PV");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(10, "TPC refit");
  fHistSingleParticleCuts[1]->GetXaxis()->SetBinLabel(11, "Kink");
  fHistogramsNeg->Add(fHistSingleParticleCuts[1]);

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
        new TH2F("fHistDecayVertexXBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex x (cm)", 100, 0, 10,
                 500, 0, 200);
    fHistogramsBefore->Add(fHistDecayVertexXBefore);

    fHistDecayVertexYBefore =
        new TH2F("fHistDecayVertexYBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex y (cm)", 100, 0, 10,
                 500, 0, 200);
    fHistogramsBefore->Add(fHistDecayVertexYBefore);

    fHistDecayVertexZBefore =
        new TH2F("fHistDecayVertexZBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex z (cm)", 100, 0, 10,
                 500, 0, 200);
    fHistogramsBefore->Add(fHistDecayVertexZBefore);

    fHistDecayVertexXAfter =
        new TH2F("fHistDecayVertexXAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex x (cm)", 100, 0, 10,
                 500, 0, 200);
    fHistogramsAfter->Add(fHistDecayVertexXAfter);

    fHistDecayVertexYAfter =
        new TH2F("fHistDecayVertexYAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex y (cm)", 100, 0, 10,
                 500, 0, 200);
    fHistogramsAfter->Add(fHistDecayVertexYAfter);

    fHistDecayVertexZAfter =
        new TH2F("fHistDecayVertexZAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Decay vertex z (cm)", 100, 0, 10,
                 500, 0, 200);
    fHistogramsAfter->Add(fHistDecayVertexZAfter);

    fHistTransverseRadiusBefore =
        new TH2F("fHistTransverseRadiusBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Transverse radius (cm)", 100, 0,
                 10, 500, 0, 200);
    fHistogramsBefore->Add(fHistTransverseRadiusBefore);

    fHistTransverseRadiusAfter =
        new TH2F("fHistTransverseRadiusAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Transverse radius (cm)", 100, 0,
                 10, 500, 0, 200);
    fHistogramsAfter->Add(fHistTransverseRadiusAfter);

    fHistCosPABefore =
        new TH2F("fHistCosPABefore", "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)",
                 100, 0, 10, 500, 0.8, 1);
    fHistogramsBefore->Add(fHistCosPABefore);

    fHistCosPAAfter =
        new TH2F("fHistCosPAAfter", "; #it{p}_{T} (GeV/#it{c}); cos(#alpha)",
                 100, 0, 10, 500, 0.8, 1);
    fHistogramsAfter->Add(fHistCosPAAfter);

    fHistDCADaughtersBefore =
        new TH2F("fHistDCADaughtersBefore",
                 "; #it{p}_{T} (GeV/#it{c}); Daughter DCA at decay vertex (cm)",
                 100, 0, 10, 500, 0, 2);
    fHistogramsBefore->Add(fHistDCADaughtersBefore);

    fHistDCADaughtersAfter =
        new TH2F("fHistDCADaughtersAfter",
                 "; #it{p}_{T} (GeV/#it{c}); Daughter DCA at decay vertex (cm)",
                 100, 0, 10, 500, 0, 2);
    fHistogramsAfter->Add(fHistDCADaughtersAfter);

    fHistDCA = new TH2F("fHistDCA", "; #it{p}_{T} (GeV/#it{c}); DCA to PV (cm)",
                        100, 0, 10, 500, 0, 10);
    fHistogramsAfter->Add(fHistDCA);

    fHistDecayLength = new TH2F("fHistDecayLength",
                                "; #it{p}_{T} (GeV/#it{c}); Decay length (cm)",
                                100, 0, 10, 500, 0, 200);
    fHistogramsAfter->Add(fHistDecayLength);

    fHistArmenterosBefore =
        new TH2F("fHistArmenterosBefore", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 500, -1, 1, 500, 0, 1);
    fHistogramsBefore->Add(fHistArmenterosBefore);

    fHistArmenterosAfter =
        new TH2F("fHistArmenterosAfter", " ; #alpha; #it{q}_{T} (GeV/#it{c})",
                 500, -1, 1, 500, 0, 1);
    fHistogramsAfter->Add(fHistArmenterosAfter);

    fHistograms->Add(fHistogramsBefore);
    fHistograms->Add(fHistogramsAfter);

    fHistSingleParticlePt[0] = new TH1F("fHistSingleParticlePt_pos",
                                        ";#it{p}_{T}; Entries", 1000, 0, 10);
    fHistSingleParticlePt[1] = new TH1F("fHistSingleParticlePt_neg",
                                        ";#it{p}_{T}; Entries", 1000, 0, 10);
    fHistSingleParticleEtaBefore[0] =
        new TH2F("fHistSingleParticleEtaBefore_pos",
                 "; #it{p}_{T} (GeV/#it{c}); #eta", 100, 0, 10, 500, -2, 2);
    fHistSingleParticleEtaBefore[1] =
        new TH2F("fHistSingleParticleEtaBefore_neg",
                 "; #it{p}_{T} (GeV/#it{c}); #eta", 100, 0, 10, 500, -2, 2);
    fHistSingleParticleEtaAfter[0] =
        new TH2F("fHistSingleParticleEtaAfter_pos",
                 "; #it{p}_{T} (GeV/#it{c}); #eta", 100, 0, 10, 500, -2, 2);
    fHistSingleParticleEtaAfter[1] =
        new TH2F("fHistSingleParticleEtaAfter_neg",
                 "; #it{p}_{T} (GeV/#it{c}); #eta", 100, 0, 10, 500, -2, 2);
    fHistSingleParticleChi2Before[0] =
        new TH2F("fHistSingleParticleChi2Before_pos",
                 "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
    fHistSingleParticleChi2Before[1] =
        new TH2F("fHistSingleParticleChi2Before_neg",
                 "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
    fHistSingleParticleChi2After[0] =
        new TH2F("fHistSingleParticleChi2After_pos",
                 "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
    fHistSingleParticleChi2After[1] =
        new TH2F("fHistSingleParticleChi2After_neg",
                 "; #it{p}_{T} (GeV/#it{c}); #chi^{2}", 100, 0, 10, 500, 0, 20);
    fHistSingleParticleNclsTPCBefore[0] = new TH2F(
        "fHistSingleParticleNclsTPCBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 100, 0, 10, 170, 0, 170);
    fHistSingleParticleNclsTPCBefore[1] = new TH2F(
        "fHistSingleParticleNclsTPCBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 100, 0, 10, 170, 0, 170);
    fHistSingleParticleNclsTPCAfter[0] = new TH2F(
        "fHistSingleParticleNclsTPCAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 100, 0, 10, 170, 0, 170);
    fHistSingleParticleNclsTPCAfter[1] = new TH2F(
        "fHistSingleParticleNclsTPCAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC ", 100, 0, 10, 170, 0, 170);
    fHistSingleParticleNclsTPCFindableBefore[0] =
        new TH2F("fHistSingleParticleNclsTPCFindableBefore_pos",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 100, 0, 10,
                 170, 0, 170);
    fHistSingleParticleNclsTPCFindableBefore[1] =
        new TH2F("fHistSingleParticleNclsTPCFindableBefore_neg",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 100, 0, 10,
                 170, 0, 170);
    fHistSingleParticleNclsTPCFindableAfter[0] =
        new TH2F("fHistSingleParticleNclsTPCFindableAfter_pos",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 100, 0, 10,
                 170, 0, 170);
    fHistSingleParticleNclsTPCFindableAfter[1] =
        new TH2F("fHistSingleParticleNclsTPCFindableAfter_neg",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC findable", 100, 0, 10,
                 170, 0, 170);
    fHistSingleParticleNclsTPCRatioFindableBefore[0] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableBefore_pos",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 100, 0, 10, 500, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableBefore[1] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableBefore_neg",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 100, 0, 10, 500, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableAfter[0] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableAfter_pos",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 100, 0, 10, 500, 0, 2);
    fHistSingleParticleNclsTPCRatioFindableAfter[1] = new TH2F(
        "fHistSingleParticleNclsTPCRatioFindableAfter_neg",
        "; #it{p}_{T} (GeV/#it{c}); TPC ratio findable", 100, 0, 10, 500, 0, 2);
    fHistSingleParticleNcrossedTPCBefore[0] =
        new TH2F("fHistSingleParticleNcrossedTPCBefore_pos",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 100, 0, 10,
                 170, 0, 170);
    fHistSingleParticleNcrossedTPCBefore[1] =
        new TH2F("fHistSingleParticleNcrossedTPCBefore_neg",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 100, 0, 10,
                 170, 0, 170);
    fHistSingleParticleNcrossedTPCAfter[0] =
        new TH2F("fHistSingleParticleNcrossedTPCAfter_pos",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 100, 0, 10,
                 170, 0, 170);
    fHistSingleParticleNcrossedTPCAfter[1] =
        new TH2F("fHistSingleParticleNcrossedTPCAfter_neg",
                 "; #it{p}_{T} (GeV/#it{c}); # cls TPC crossed", 100, 0, 10,
                 170, 0, 170);
    fHistSingleParticleNclsTPCShared[0] = new TH2F(
        "fHistSingleParticleNclsTPCShared_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared", 100, 0, 10, 170, 0, 170);
    fHistSingleParticleNclsTPCShared[1] = new TH2F(
        "fHistSingleParticleNclsTPCShared_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls TPC shared", 100, 0, 10, 170, 0, 170);
    fHistSingleParticleNclsITSShared[0] = new TH2F(
        "fHistSingleParticleNclsITSShared_pos",
        "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared", 100, 0, 10, 6, 0, 6);
    fHistSingleParticleNclsITSShared[1] = new TH2F(
        "fHistSingleParticleNclsITSShared_neg",
        "; #it{p}_{T} (GeV/#it{c}); # cls ITS shared", 100, 0, 10, 6, 0, 6);
    fHistSingleParticleDCAtoPVBefore[0] = new TH2F(
        "fHistSingleParticleDCAtoPVBefore_pos",
        "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA)", 100, 0, 10, 500, 0, 10);
    fHistSingleParticleDCAtoPVBefore[1] = new TH2F(
        "fHistSingleParticleDCAtoPVBefore_neg",
        "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA)", 100, 0, 10, 500, 0, 10);
    fHistSingleParticleDCAtoPVAfter[0] = new TH2F(
        "fHistSingleParticleDCAtoPVAfter_pos",
        "; #it{p}_{T} (GeV/#it{c});DCA (pos, PA)", 100, 0, 10, 500, 0, 10);
    fHistSingleParticleDCAtoPVAfter[1] = new TH2F(
        "fHistSingleParticleDCAtoPVAfter_neg",
        "; #it{p}_{T} (GeV/#it{c});DCA (neg, PA)", 100, 0, 10, 500, 0, 10);
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

    fHistSingleParticlePID[0] = new TH2F(
        "fHistSingleParticlePID_pos", "; #it{p} (GeV/#it{c}); n_{#sigma} TPC",
        1000, 0, 5, 1000, -10, 10);
    fHistSingleParticlePID[1] = new TH2F(
        "fHistSingleParticlePID_neg", "; #it{p} (GeV/#it{c}); n_{#sigma} TPC",
        1000, 0, 5, 1000, -10, 10);

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
    fHistogramsPos->Add(fHistSingleParticleNclsTPCShared[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsTPCShared[1]);
    fHistogramsPos->Add(fHistSingleParticleNclsITSShared[0]);
    fHistogramsNeg->Add(fHistSingleParticleNclsITSShared[1]);
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

    fHistMCTruthV0Pt = new TH1F(
        "fHistMCTruthV0Pt", "; #it{p}_{T} (GeV/#it{c}); Entries", 500, 0, 10);
    fHistMCTruthV0PtY =
        new TH2F("fHistMCTruthV0PtY", "; y; #it{p}_{T} (GeV/#it{c})", 500, -10,
                 10, 500, 0, 10);
    fHistMCTruthV0PtEta =
        new TH2F("fHistMCTruthV0PtEta", "; #eta; #it{p}_{T} (GeV/#it{c})", 500,
                 -10, 10, 500, 0, 10);
    fHistMCTruthV0DaughterPt =
        new TH1F("fHistMCTruthV0DaughterPt",
                 "; #it{p}_{T} (GeV/#it{c}); Entries", 500, 0, 10);
    fHistMCTruthV0DaughterPtY =
        new TH2F("fHistMCTruthV0DaughterPtY", "; y; #it{p}_{T} (GeV/#it{c})",
                 500, -10, 10, 500, 0, 10);
    fHistMCTruthV0DaughterPtEta =
        new TH2F("fHistMCTruthV0DaughterPtEta",
                 "; #eta; #it{p}_{T} (GeV/#it{c})", 500, -10, 10, 500, 0, 10);

    fHistMCV0Pt = new TH1F("fHistMCV0Pt", "; #it{p}_{T} (GeV/#it{c}); Entries",
                           500, 0, 10);

    fHistogramsMC->Add(fHistMCTruthV0Pt);
    fHistogramsMC->Add(fHistMCTruthV0PtY);
    fHistogramsMC->Add(fHistMCTruthV0PtEta);
    fHistogramsMC->Add(fHistMCTruthV0DaughterPt);
    fHistogramsMC->Add(fHistMCTruthV0DaughterPtY);
    fHistogramsMC->Add(fHistMCTruthV0DaughterPtEta);
    fHistogramsMC->Add(fHistMCV0Pt);
    fHistograms->Add(fHistogramsMC);
  }
}
