#include "AliSigma0V0Cuts.h"

#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include <iostream>

ClassImp(AliSigma0V0Cuts)

//____________________________________________________________________________________________________
AliSigma0V0Cuts::AliSigma0V0Cuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsV0(nullptr),
      fHistogramsV0MC(nullptr),
      fHistogramsV0Before(nullptr),
      fHistogramsV0After(nullptr),
      fHistogramsV0Pos(nullptr),
      fHistogramsV0Neg(nullptr),
      fHistogramsAntiV0(nullptr),
      fHistogramsAntiV0MC(nullptr),
      fHistogramsAntiV0Before(nullptr),
      fHistogramsAntiV0After(nullptr),
      fHistogramsAntiV0Pos(nullptr),
      fHistogramsAntiV0Neg(nullptr),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fSingleParticleCuts(nullptr),
      fV0Vector(),
      fAntiV0Vector(),
      fGlobalTrackReference(nullptr),
      fV0cut(0),
      fAntiV0cut(0),
      fPID(0),
      fIsMC(false),
	  fPileUpRejection(false),
      fIsExtendedQA(false),
      fV0OnFly(false),
      fArmenterosCut(false),
      fUsePID(false),
      fV0PtMin(0.f),
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
      fHistCuts(nullptr),
      fHistV0Cuts(nullptr),
      fHistNV0(nullptr),
      fHistV0LambdaMass(nullptr),
      fHistV0LambdaPt(nullptr),
      fHistV0LambdaPtY(),
      fHistV0LambdaMassPt(nullptr),
      fHistV0LambdaMassK0Rej(nullptr),
      fHistV0K0Mass(nullptr),
      fHistV0K0MassAfter(nullptr),
      fHistV0CosPA(nullptr),
      fHistV0EtaPhi(nullptr),
      fHistV0DecayVertexXBefore(nullptr),
      fHistV0DecayVertexYBefore(nullptr),
      fHistV0DecayVertexZBefore(nullptr),
      fHistV0DecayVertexXAfter(nullptr),
      fHistV0DecayVertexYAfter(nullptr),
      fHistV0DecayVertexZAfter(nullptr),
      fHistV0TransverseRadiusBefore(nullptr),
      fHistV0TransverseRadiusAfter(nullptr),
      fHistV0CosPABefore(nullptr),
      fHistV0CosPAAfter(nullptr),
      fHistV0DCADaughtersBefore(nullptr),
      fHistV0DCADaughtersAfter(nullptr),
      fHistV0DCA(nullptr),
      fHistV0DecayLength(nullptr),
      fHistV0ArmenterosBefore(nullptr),
      fHistV0ArmenterosAfter(nullptr),
      fHistMCTruthV0Pt(nullptr),
      fHistMCTruthV0PtY(nullptr),
      fHistMCTruthV0PtEta(nullptr),
      fHistMCTruthV0ProtonPionPt(nullptr),
      fHistMCTruthV0ProtonPionPtY(nullptr),
      fHistMCTruthV0ProtonPionPtEta(nullptr),
      fHistAntiV0Cuts(nullptr),
      fHistNAntiV0(nullptr),
      fHistAntiV0LambdaMass(nullptr),
      fHistAntiV0LambdaPt(nullptr),
      fHistAntiV0LambdaPtY(),
      fHistAntiV0LambdaMassPt(nullptr),
      fHistAntiV0LambdaMassK0Rej(nullptr),
      fHistAntiV0K0Mass(nullptr),
      fHistAntiV0K0MassAfter(nullptr),
      fHistAntiV0CosPA(nullptr),
      fHistAntiV0EtaPhi(nullptr),
      fHistAntiV0DecayVertexXBefore(nullptr),
      fHistAntiV0DecayVertexYBefore(nullptr),
      fHistAntiV0DecayVertexZBefore(nullptr),
      fHistAntiV0DecayVertexXAfter(nullptr),
      fHistAntiV0DecayVertexYAfter(nullptr),
      fHistAntiV0DecayVertexZAfter(nullptr),
      fHistAntiV0TransverseRadiusBefore(nullptr),
      fHistAntiV0TransverseRadiusAfter(nullptr),
      fHistAntiV0CosPABefore(nullptr),
      fHistAntiV0CosPAAfter(nullptr),
      fHistAntiV0DCADaughtersBefore(nullptr),
      fHistAntiV0DCADaughtersAfter(nullptr),
      fHistAntiV0DCA(nullptr),
      fHistAntiV0DecayLength(nullptr),
      fHistAntiV0ArmenterosBefore(nullptr),
      fHistAntiV0ArmenterosAfter(nullptr),
      fHistMCTruthAntiV0Pt(nullptr),
      fHistMCTruthAntiV0PtY(nullptr),
      fHistMCTruthAntiV0PtEta(nullptr),
      fHistMCTruthAntiV0ProtonPionPt(nullptr),
      fHistMCTruthAntiV0ProtonPionPtY(nullptr),
      fHistMCTruthAntiV0ProtonPionPtEta(nullptr),
      fHistV0SingleParticleCuts(),
      fHistV0SingleParticlePt(),
      fHistV0SingleParticleEtaBefore(),
      fHistV0SingleParticleEtaAfter(),
      fHistV0SingleParticleNclsTPCBefore(),
      fHistV0SingleParticleNclsTPCAfter(),
      fHistV0SingleParticleNclsTPCFindableBefore(),
      fHistV0SingleParticleNclsTPCFindableAfter(),
      fHistV0SingleParticleNclsTPCRatioFindableBefore(),
      fHistV0SingleParticleNclsTPCRatioFindableAfter(),
      fHistV0SingleParticleNcrossedTPCBefore(),
      fHistV0SingleParticleNcrossedTPCAfter(),
      fHistV0SingleParticleNclsTPCShared(),
      fHistV0SingleParticleNclsTPCSharedTiming(),
      fHistV0SingleParticleNclsITSShared(),
      fHistV0SingleParticleNclsITSSharedTiming(),
      fHistV0SingleParticleDCAtoPVBefore(),
      fHistV0SingleParticleDCAtoPVAfter(),
      fHistV0SingleParticlePID(),
      fHistAntiV0SingleParticleCuts(),
      fHistAntiV0SingleParticlePt(),
      fHistAntiV0SingleParticleEtaBefore(),
      fHistAntiV0SingleParticleEtaAfter(),
      fHistAntiV0SingleParticleNclsTPCBefore(),
      fHistAntiV0SingleParticleNclsTPCAfter(),
      fHistAntiV0SingleParticleNclsTPCFindableBefore(),
      fHistAntiV0SingleParticleNclsTPCFindableAfter(),
      fHistAntiV0SingleParticleNclsTPCRatioFindableBefore(),
      fHistAntiV0SingleParticleNclsTPCRatioFindableAfter(),
      fHistAntiV0SingleParticleNcrossedTPCBefore(),
      fHistAntiV0SingleParticleNcrossedTPCAfter(),
      fHistAntiV0SingleParticleNclsTPCShared(),
      fHistAntiV0SingleParticleNclsTPCSharedTiming(),
      fHistAntiV0SingleParticleNclsITSShared(),
      fHistAntiV0SingleParticleNclsITSSharedTiming(),
      fHistAntiV0SingleParticleDCAtoPVBefore(),
      fHistAntiV0SingleParticleDCAtoPVAfter(),
      fHistAntiV0SingleParticlePID() {}

//____________________________________________________________________________________________________
AliSigma0V0Cuts::AliSigma0V0Cuts(const AliSigma0V0Cuts &ref)
    : TObject(ref),
      fHistograms(nullptr),
      fHistogramsV0(nullptr),
      fHistogramsV0MC(nullptr),
      fHistogramsV0Before(nullptr),
      fHistogramsV0After(nullptr),
      fHistogramsV0Pos(nullptr),
      fHistogramsV0Neg(nullptr),
      fHistogramsAntiV0(nullptr),
      fHistogramsAntiV0MC(nullptr),
      fHistogramsAntiV0Before(nullptr),
      fHistogramsAntiV0After(nullptr),
      fHistogramsAntiV0Pos(nullptr),
      fHistogramsAntiV0Neg(nullptr),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fSingleParticleCuts(nullptr),
      fV0Vector(),
      fAntiV0Vector(),
      fGlobalTrackReference(nullptr),
      fV0cut(0),
      fAntiV0cut(0),
      fPID(0),
      fIsMC(false),
	  fPileUpRejection(false),
      fIsExtendedQA(false),
      fV0OnFly(false),
      fArmenterosCut(false),
      fUsePID(false),
      fV0PtMin(0.f),
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
      fHistCuts(nullptr),
      fHistV0Cuts(nullptr),
      fHistNV0(nullptr),
      fHistV0LambdaMass(nullptr),
      fHistV0LambdaPt(nullptr),
      fHistV0LambdaPtY(),
      fHistV0LambdaMassPt(nullptr),
      fHistV0LambdaMassK0Rej(nullptr),
      fHistV0K0Mass(nullptr),
      fHistV0K0MassAfter(nullptr),
      fHistV0CosPA(nullptr),
      fHistV0EtaPhi(nullptr),
      fHistV0DecayVertexXBefore(nullptr),
      fHistV0DecayVertexYBefore(nullptr),
      fHistV0DecayVertexZBefore(nullptr),
      fHistV0DecayVertexXAfter(nullptr),
      fHistV0DecayVertexYAfter(nullptr),
      fHistV0DecayVertexZAfter(nullptr),
      fHistV0TransverseRadiusBefore(nullptr),
      fHistV0TransverseRadiusAfter(nullptr),
      fHistV0CosPABefore(nullptr),
      fHistV0CosPAAfter(nullptr),
      fHistV0DCADaughtersBefore(nullptr),
      fHistV0DCADaughtersAfter(nullptr),
      fHistV0DCA(nullptr),
      fHistV0DecayLength(nullptr),
      fHistV0ArmenterosBefore(nullptr),
      fHistV0ArmenterosAfter(nullptr),
      fHistMCTruthV0Pt(nullptr),
      fHistMCTruthV0PtY(nullptr),
      fHistMCTruthV0PtEta(nullptr),
      fHistMCTruthV0ProtonPionPt(nullptr),
      fHistMCTruthV0ProtonPionPtY(nullptr),
      fHistMCTruthV0ProtonPionPtEta(nullptr),
      fHistAntiV0Cuts(nullptr),
      fHistNAntiV0(nullptr),
      fHistAntiV0LambdaMass(nullptr),
      fHistAntiV0LambdaPt(nullptr),
      fHistAntiV0LambdaPtY(),
      fHistAntiV0LambdaMassPt(nullptr),
      fHistAntiV0LambdaMassK0Rej(nullptr),
      fHistAntiV0K0Mass(nullptr),
      fHistAntiV0K0MassAfter(nullptr),
      fHistAntiV0CosPA(nullptr),
      fHistAntiV0EtaPhi(nullptr),
      fHistAntiV0DecayVertexXBefore(nullptr),
      fHistAntiV0DecayVertexYBefore(nullptr),
      fHistAntiV0DecayVertexZBefore(nullptr),
      fHistAntiV0DecayVertexXAfter(nullptr),
      fHistAntiV0DecayVertexYAfter(nullptr),
      fHistAntiV0DecayVertexZAfter(nullptr),
      fHistAntiV0TransverseRadiusBefore(nullptr),
      fHistAntiV0TransverseRadiusAfter(nullptr),
      fHistAntiV0CosPABefore(nullptr),
      fHistAntiV0CosPAAfter(nullptr),
      fHistAntiV0DCADaughtersBefore(nullptr),
      fHistAntiV0DCADaughtersAfter(nullptr),
      fHistAntiV0DCA(nullptr),
      fHistAntiV0DecayLength(nullptr),
      fHistAntiV0ArmenterosBefore(nullptr),
      fHistAntiV0ArmenterosAfter(nullptr),
      fHistMCTruthAntiV0Pt(nullptr),
      fHistMCTruthAntiV0PtY(nullptr),
      fHistMCTruthAntiV0PtEta(nullptr),
      fHistMCTruthAntiV0ProtonPionPt(nullptr),
      fHistMCTruthAntiV0ProtonPionPtY(nullptr),
      fHistMCTruthAntiV0ProtonPionPtEta(nullptr),
      fHistV0SingleParticleCuts(),
      fHistV0SingleParticlePt(),
      fHistV0SingleParticleEtaBefore(),
      fHistV0SingleParticleEtaAfter(),
      fHistV0SingleParticleNclsTPCBefore(),
      fHistV0SingleParticleNclsTPCAfter(),
      fHistV0SingleParticleNclsTPCFindableBefore(),
      fHistV0SingleParticleNclsTPCFindableAfter(),
      fHistV0SingleParticleNclsTPCRatioFindableBefore(),
      fHistV0SingleParticleNclsTPCRatioFindableAfter(),
      fHistV0SingleParticleNcrossedTPCBefore(),
      fHistV0SingleParticleNcrossedTPCAfter(),
      fHistV0SingleParticleNclsTPCShared(),
      fHistV0SingleParticleNclsTPCSharedTiming(),
      fHistV0SingleParticleNclsITSShared(),
      fHistV0SingleParticleNclsITSSharedTiming(),
      fHistV0SingleParticleDCAtoPVBefore(),
      fHistV0SingleParticleDCAtoPVAfter(),
      fHistV0SingleParticlePID(),
      fHistAntiV0SingleParticleCuts(),
      fHistAntiV0SingleParticlePt(),
      fHistAntiV0SingleParticleEtaBefore(),
      fHistAntiV0SingleParticleEtaAfter(),
      fHistAntiV0SingleParticleNclsTPCBefore(),
      fHistAntiV0SingleParticleNclsTPCAfter(),
      fHistAntiV0SingleParticleNclsTPCFindableBefore(),
      fHistAntiV0SingleParticleNclsTPCFindableAfter(),
      fHistAntiV0SingleParticleNclsTPCRatioFindableBefore(),
      fHistAntiV0SingleParticleNclsTPCRatioFindableAfter(),
      fHistAntiV0SingleParticleNcrossedTPCBefore(),
      fHistAntiV0SingleParticleNcrossedTPCAfter(),
      fHistAntiV0SingleParticleNclsTPCShared(),
      fHistAntiV0SingleParticleNclsTPCSharedTiming(),
      fHistAntiV0SingleParticleNclsITSShared(),
      fHistAntiV0SingleParticleNclsITSSharedTiming(),
      fHistAntiV0SingleParticleDCAtoPVBefore(),
      fHistAntiV0SingleParticleDCAtoPVAfter(),
      fHistAntiV0SingleParticlePID() {}

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
AliSigma0V0Cuts *AliSigma0V0Cuts::DefaultCuts() {
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
  return v0Cuts;
}

//____________________________________________________________________________________________________
AliSigma0V0Cuts *AliSigma0V0Cuts::Sigma0Cuts() {
  AliSigma0V0Cuts *v0Cuts = new AliSigma0V0Cuts();
  v0Cuts->SetV0OnFlyStatus(true);
  v0Cuts->SetV0PtMin(0.f);
  v0Cuts->SetDaughterDCAtoPV(0.06);
  v0Cuts->SetDaughterDCAMax(1.5);
  v0Cuts->SetV0CosPAMin(0.993f);
  v0Cuts->SetV0RadiusMin(0.f);
  v0Cuts->SetV0RadiusMax(180.f);
  v0Cuts->SetArmenterosCut(0.01, 0.17, 0.2, 0.9);
  v0Cuts->SetLambdaSelection(1.094, 1.136);
  return v0Cuts;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::SelectV0s(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    std::vector<AliSigma0ParticleV0> &V0Container,
    std::vector<AliSigma0ParticleV0> &AntiV0Container,
    AliPID::EParticleType particle1, AliPID::EParticleType particle2) {
  // Get the PID manager from the input event handler
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      static_cast<AliInputEventHandler *>(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fMCEvent = mcEvent;

  fV0Vector.clear();
  fAntiV0Vector.clear();

  fInputEvent = inputEvent;

  if (fInputEvent->IsA() == AliESDEvent::Class()) {
    ProcessESDs(particle1, particle2);
  }
  if (fInputEvent->IsA() == AliAODEvent::Class()) {
    ProcessAODs(particle1, particle2);
  }

  if (fIsMC) ProcessMC();

  fHistNV0->Fill(fV0Vector.size());
  fHistNAntiV0->Fill(fAntiV0Vector.size());

  V0Container.swap(fV0Vector);
  AntiV0Container.swap(fAntiV0Vector);
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::ProcessESDs(AliPID::EParticleType particle,
                                  AliPID::EParticleType antiParticle) {
  AliESDEvent *event = static_cast<AliESDEvent *>(fInputEvent);
  for (int iV0 = 0; iV0 < event->GetNumberOfV0s(); ++iV0) {
    AliESDv0 *v0 = event->GetV0(iV0);

    // V0 Quality
    if (!V0QualityCuts(v0)) continue;

    AliVTrack *pos = static_cast<AliVTrack *>(event->GetTrack(v0->GetPindex()));
    AliVTrack *neg = static_cast<AliVTrack *>(event->GetTrack(v0->GetNindex()));
    if (!pos || !neg) continue;

    if (pos->Charge() < 0) {
      pos = neg;
      neg = event->GetTrack(v0->GetPindex());
    }

    // PID
    const float armAlpha = v0->AlphaV0();
    const float armQt = v0->PtArmV0();
    fHistV0ArmenterosBefore->Fill(armAlpha, armQt);
    fHistAntiV0ArmenterosBefore->Fill(armAlpha, armQt);
    if (fUsePID && !V0PID(pos, neg, particle, antiParticle)) continue;
    if (!fUsePID) {
      if (v0->AlphaV0() > 0)
        fPID = 1;
      else
        fPID = -1;
      if (fPID == 1)
        fHistV0Cuts->Fill(7);
      else
        fHistAntiV0Cuts->Fill(7);
    }
    if (fPID == 1 && fIsExtendedQA) {
      PlotSingleParticlePID(pos, particle);
      PlotSingleParticlePID(neg, antiParticle);
    } else if (fPID == -1 && fIsExtendedQA) {
      PlotSingleParticlePID(pos, antiParticle);
      PlotSingleParticlePID(neg, particle);
    }

    // Single Particle Quality
    const float magField = fInputEvent->GetMagneticField();
    const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
    const float dcaPosToPV = std::fabs(static_cast<AliESDtrack *>(pos)->GetD(
        vertex->GetX(), vertex->GetY(), magField));
    const float dcaNegToPV = std::fabs(static_cast<AliESDtrack *>(neg)->GetD(
        vertex->GetX(), vertex->GetY(), magField));
    if (!SingleParticleQualityCuts(pos, dcaPosToPV) ||
        !SingleParticleQualityCuts(neg, dcaNegToPV))
      continue;
    if (fPID == 1)
      fHistV0Cuts->Fill(8);
    else
      fHistAntiV0Cuts->Fill(8);

    // Topological Selection
    if (!V0TopologicalSelection(v0)) continue;

    // Particle Selection
    if (fPID == 1)
      v0->ChangeMassHypothesis(3122);
    else
      v0->ChangeMassHypothesis(-3122);
    const float massLambda = v0->GetEffMass();
    v0->ChangeMassHypothesis(310);
    const float massK0 = v0->GetEffMass();
    if (!LambdaSelection(v0, massK0, massLambda)) continue;

    AliSigma0ParticleV0 *v0Candidate = new AliSigma0ParticleV0(
        *v0, *pos, *neg, fInputEvent->GetPrimaryVertex(), 3122 * fPID, fPID,
        fInputEvent->GetMagneticField(), fMCEvent);

    if (fIsMC) {
      // TClonesArray *mcarray = dynamic_cast<TClonesArray *>(
      //      fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      // if (!mcarray) continue;
      /// @todo Fix matching to MC!
      //      Int_t daugh[2] = {2212, 211};  // daughter IDs of Lambda
      //      Int_t label = v0->MatchToMC(3122,mcarray,2,daugh);
      //      if(label<0) {
      //        // no mc info assigned to this track - don't use it
      //        v0Candidate->SetUse(false);
      //      }

      //      AliMCParticle *mcParticle  =
      //      static_cast<AliMCParticle*>(fMCEvent->GetTrack(label));
      //      if(!mcParticle) continue;
      //      v0Candidate->ProcessMCInfo(mcParticle, fMCEvent);
    }

    /// @todo set mass
    v0Candidate->SetPDGMass(1.115683);

    if (fPID == 1) {
      fV0Vector.push_back(*v0Candidate);
    } else {
      fAntiV0Vector.push_back(*v0Candidate);
    }

    delete v0Candidate;
  }
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::ProcessAODs(AliPID::EParticleType particle,
                                  AliPID::EParticleType antiParticle) {
  // Get the pointer to the global track reference from the single particle cut
  // class
  fGlobalTrackReference = fSingleParticleCuts->GetGlobalTrackReference();

  for (auto v0obj : *(static_cast<AliAODEvent *>(fInputEvent)->GetV0s())) {
    AliAODv0 *v0 = static_cast<AliAODv0 *>(v0obj);

    // V0 Quality
    if (!V0QualityCuts(v0)) continue;

    AliVTrack *pos =
        static_cast<AliVTrack *>((*fGlobalTrackReference)[v0->GetPosID()]);
    AliVTrack *neg =
        static_cast<AliVTrack *>((*fGlobalTrackReference)[v0->GetNegID()]);
    if (!pos || !neg) continue;

    if (pos->Charge() < 0) {
      pos = neg;
      neg = (*fGlobalTrackReference)[v0->GetPosID()];
    }

    // PID
    const float armAlpha = v0->AlphaV0();
    const float armQt = v0->PtArmV0();
    fHistV0ArmenterosBefore->Fill(armAlpha, armQt);
    fHistAntiV0ArmenterosBefore->Fill(armAlpha, armQt);
    if (fUsePID && !V0PID(pos, neg, particle, antiParticle)) continue;
    if (!fUsePID) {
      if (armAlpha > 0)
        fPID = 1;
      else
        fPID = -1;
      if (fPID == 1)
        fHistV0Cuts->Fill(7);
      else
        fHistAntiV0Cuts->Fill(7);
    }
    if (fPID == 1 && fIsExtendedQA) {
      PlotSingleParticlePID(pos, particle);
      PlotSingleParticlePID(neg, antiParticle);
    } else if (fPID == -1 && fIsExtendedQA) {
      PlotSingleParticlePID(pos, antiParticle);
      PlotSingleParticlePID(neg, particle);
    }

    // Single Particle Quality
    const float dcaPosToPV = v0->DcaPosToPrimVertex();
    const float dcaNegToPV = v0->DcaNegToPrimVertex();
    if (!SingleParticleQualityCuts(pos, dcaPosToPV) ||
        !SingleParticleQualityCuts(neg, dcaNegToPV))
      continue;
    if (fPID == 1)
      fHistV0Cuts->Fill(8);
    else
      fHistAntiV0Cuts->Fill(8);

    // Topological Selection
    if (!V0TopologicalSelection(v0)) continue;

    // Particle Selection
    const float massK0 = v0->MassK0Short();
    const float massLambda =
        (fPID == 1) ? v0->MassLambda() : v0->MassAntiLambda();
    if (!LambdaSelection(v0, massK0, massLambda)) continue;

    AliSigma0ParticleV0 *v0Candidate = new AliSigma0ParticleV0(
        *v0, *pos, *neg, fInputEvent->GetPrimaryVertex(), 3122 * fPID, fPID,
        fInputEvent->GetMagneticField(), fMCEvent);

    if (fIsMC) {
      TClonesArray *mcarray = dynamic_cast<TClonesArray *>(
          fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!mcarray) continue;
      Int_t daugh[2] = {2212, 211};  // daughter IDs of Lambda
      Int_t label = v0->MatchToMC(3122, mcarray, 2, daugh);
      if (label < 0) {
        // no mc info assigned to this track - don't use it
        v0Candidate->SetUse(false);
      }

      AliMCParticle *mcParticle =
          static_cast<AliMCParticle *>(fMCEvent->GetTrack(label));
      if (!mcParticle) continue;
      v0Candidate->ProcessMCInfo(mcParticle, fMCEvent);
    }

    /// @todo set mass
    v0Candidate->SetPDGMass(1.115683);

    if (fPID == 1) {
      fV0Vector.push_back(*v0Candidate);
    } else {
      fAntiV0Vector.push_back(*v0Candidate);
    }

    delete v0Candidate;
  }
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0QualityCuts(const AliAODv0 *v0) {
  fHistV0Cuts->Fill(0);
  fHistAntiV0Cuts->Fill(0);

  // Two prongs
  if (v0->GetNProngs() > 2) return false;
  fHistV0Cuts->Fill(1);
  fHistAntiV0Cuts->Fill(1);

  // Two daughters
  if (v0->GetNDaughters() > 2) return false;
  fHistV0Cuts->Fill(2);
  fHistAntiV0Cuts->Fill(2);

  // v0 is neutral
  if (v0->GetCharge() != 0) return false;
  fHistV0Cuts->Fill(3);
  fHistAntiV0Cuts->Fill(3);

  // on-fly status
  if (fV0OnFly) {
    if (!(v0->GetOnFlyStatus())) return false;  // online v0
  } else {
    if (v0->GetOnFlyStatus()) return false;  // offline v0
  }
  fHistV0Cuts->Fill(4);
  fHistAntiV0Cuts->Fill(4);

  // pt min cut
  if (v0->Pt() < fV0PtMin) return false;
  fHistV0Cuts->Fill(5);
  fHistAntiV0Cuts->Fill(5);

  // check for daughter tracks
  AliAODTrack *pos = (*fGlobalTrackReference)[v0->GetPosID()];
  AliAODTrack *neg = (*fGlobalTrackReference)[v0->GetNegID()];
  if (!pos || !neg) return false;
  if (pos->Charge() == neg->Charge()) return false;
  fHistV0Cuts->Fill(6);
  fHistAntiV0Cuts->Fill(6);

  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0QualityCuts(const AliESDv0 *v0) {
  fHistV0Cuts->Fill(0);
  fHistAntiV0Cuts->Fill(0);

  // Two daughters
  if (v0->GetNDaughters() > 2) return false;
  fHistV0Cuts->Fill(2);
  fHistAntiV0Cuts->Fill(2);

  // v0 is neutral
  if (v0->Charge() != 0) return false;
  fHistV0Cuts->Fill(3);
  fHistAntiV0Cuts->Fill(3);

  // on-fly status
  if (fV0OnFly) {
    if (!(v0->GetOnFlyStatus())) return false;  // online v0
  } else {
    if (v0->GetOnFlyStatus()) return false;  // offline v0
  }
  fHistV0Cuts->Fill(4);
  fHistAntiV0Cuts->Fill(4);

  // pt min cut
  if (v0->Pt() < fV0PtMin) return false;
  fHistV0Cuts->Fill(5);
  fHistAntiV0Cuts->Fill(5);

  // check for daughter tracks
  AliESDEvent *event = static_cast<AliESDEvent *>(fInputEvent);
  AliESDtrack *pos = event->GetTrack(v0->GetPindex());
  AliESDtrack *neg = event->GetTrack(v0->GetNindex());
  if (!pos || !neg) return false;
  if (pos->Charge() == neg->Charge()) return false;
  fHistV0Cuts->Fill(6);
  fHistAntiV0Cuts->Fill(6);

  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0PID(const AliVTrack *pos, const AliVTrack *neg,
                            AliPID::EParticleType particle,
                            AliPID::EParticleType antiParticle) {
  float protonProb = 999.f;
  float piMinusProb = 999.f;
  float antiProtonProb = 999.f;
  float piPlusProb = 999.f;

  bool isParticle = (SingleParticlePID(pos, particle, protonProb) &&
                     SingleParticlePID(neg, antiParticle, piMinusProb));
  bool isAntiParticle = (SingleParticlePID(neg, particle, antiProtonProb) &&
                         SingleParticlePID(pos, antiParticle, piPlusProb));

  // weighting the probabilities in case both particle and anti-particle are
  // found - not done by oli
  float particleProb =
      std::sqrt(protonProb * protonProb + piMinusProb * piMinusProb);
  float antiParticleProb =
      std::sqrt(antiProtonProb * antiProtonProb + piPlusProb * piPlusProb);

  if (isParticle && !isAntiParticle)
    fPID = 1;
  else if (isAntiParticle && !isParticle)
    fPID = -1;
  else if (isParticle && isAntiParticle) {
    fPID = (particleProb > antiParticleProb) ? -1 : 1;
  } else {
    fPID = 0;
    return false;
  }

  if (fPID == 1)
    fHistV0Cuts->Fill(7);
  else
    fHistAntiV0Cuts->Fill(7);

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

  const float nSigma =
      std::fabs(fPIDResponse->NumberOfSigmasTPC(track, particle));

  prob = nSigma;

  // Proton PID cut
  if (nSigma > fPIDnSigma) return false;
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
      std::fabs(fPIDResponse->NumberOfSigmasTPC(track, particle));

  const int charge = track->Charge();
  if (fPID == 1) {
    if (charge > 0)
      fHistV0SingleParticlePID[0]->Fill(track->P(), nSigma);
    else
      fHistV0SingleParticlePID[1]->Fill(track->P(), nSigma);
  } else {
    if (charge > 0)
      fHistAntiV0SingleParticlePID[0]->Fill(track->P(), nSigma);
    else
      fHistAntiV0SingleParticlePID[1]->Fill(track->P(), nSigma);
  }
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::V0TopologicalSelection(const AliAODv0 *v0) {
  // Get the coordinates of the primary vertex
  Double_t xPV = fInputEvent->GetPrimaryVertex()->GetX();
  Double_t yPV = fInputEvent->GetPrimaryVertex()->GetY();
  Double_t zPV = fInputEvent->GetPrimaryVertex()->GetZ();
  Double_t PV[3] = {xPV, yPV, zPV};
  Double_t decayVertexV0[3];

  v0->GetXYZ(decayVertexV0);

  // Calculate vertex variables:
  const float point = v0->CosPointingAngle(PV);
  const float dcaV0Dau = v0->DcaV0Daughters();
  const float dcaPrim = v0->DcaV0ToPrimVertex();
  const float lenDecay = v0->DecayLengthV0(PV);
  const float transverseRadius = v0->DecayLengthXY(PV);

  if (fIsExtendedQA) {
    if (fPID == 1) {
      fHistV0CosPABefore->Fill(point);
      fHistV0DecayVertexXBefore->Fill(decayVertexV0[0]);
      fHistV0DecayVertexYBefore->Fill(decayVertexV0[1]);
      fHistV0DecayVertexZBefore->Fill(decayVertexV0[2]);
      fHistV0TransverseRadiusBefore->Fill(transverseRadius);
      fHistV0DCADaughtersBefore->Fill(dcaV0Dau);
    } else {
      fHistAntiV0CosPABefore->Fill(point);
      fHistAntiV0DecayVertexXBefore->Fill(decayVertexV0[0]);
      fHistAntiV0DecayVertexYBefore->Fill(decayVertexV0[1]);
      fHistAntiV0DecayVertexZBefore->Fill(decayVertexV0[2]);
      fHistAntiV0TransverseRadiusBefore->Fill(transverseRadius);
      fHistAntiV0DCADaughtersBefore->Fill(dcaV0Dau);
    }
  }

  // Position of the decay vertex x, y & z
  if (std::fabs(decayVertexV0[0]) > fV0DecayVertexMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(9);
  else
    fHistAntiV0Cuts->Fill(9);
  if (std::fabs(decayVertexV0[1]) > fV0DecayVertexMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(10);
  else
    fHistAntiV0Cuts->Fill(10);
  if (std::fabs(decayVertexV0[2]) > fV0DecayVertexMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(11);
  else
    fHistAntiV0Cuts->Fill(11);

  // Transverse decay radius min & max
  if (transverseRadius < fV0RadiusMin) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(12);
  else
    fHistAntiV0Cuts->Fill(12);
  if (transverseRadius > fV0RadiusMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(13);
  else
    fHistAntiV0Cuts->Fill(13);

  // DCA of the daughter tracks at the decay vertex
  if (dcaV0Dau > fDaughterDCAMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(14);
  else
    fHistAntiV0Cuts->Fill(14);

  // Cos pointing angle
  if (point < fV0CosPAMin) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(15);
  else
    fHistAntiV0Cuts->Fill(15);

  if (fIsExtendedQA) {
    if (fPID == 1) {
      fHistV0CosPA->Fill(v0->Pt(), point);
      fHistV0CosPAAfter->Fill(point);
      fHistV0DecayVertexXAfter->Fill(decayVertexV0[0]);
      fHistV0DecayVertexYAfter->Fill(decayVertexV0[1]);
      fHistV0DecayVertexZAfter->Fill(decayVertexV0[2]);
      fHistV0TransverseRadiusAfter->Fill(transverseRadius);
      fHistV0DCADaughtersAfter->Fill(dcaV0Dau);
      fHistV0DCA->Fill(dcaPrim);
      fHistV0DecayLength->Fill(lenDecay);
    } else {
      fHistAntiV0CosPA->Fill(v0->Pt(), point);
      fHistAntiV0CosPAAfter->Fill(point);
      fHistAntiV0DecayVertexXAfter->Fill(decayVertexV0[0]);
      fHistAntiV0DecayVertexYAfter->Fill(decayVertexV0[1]);
      fHistAntiV0DecayVertexZAfter->Fill(decayVertexV0[2]);
      fHistAntiV0TransverseRadiusAfter->Fill(transverseRadius);
      fHistAntiV0DCADaughtersAfter->Fill(dcaV0Dau);
      fHistAntiV0DCA->Fill(dcaPrim);
      fHistAntiV0DecayLength->Fill(lenDecay);
    }
  }
  return true;
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

  if (fIsExtendedQA) {
    if (fPID == 1) {
      fHistV0CosPABefore->Fill(point);
      fHistV0DecayVertexXBefore->Fill(decayVertexV0[0]);
      fHistV0DecayVertexYBefore->Fill(decayVertexV0[1]);
      fHistV0DecayVertexZBefore->Fill(decayVertexV0[2]);
      fHistV0TransverseRadiusBefore->Fill(transverseRadius);
      fHistV0DCADaughtersBefore->Fill(dcaV0Dau);
    } else {
      fHistAntiV0CosPABefore->Fill(point);
      fHistAntiV0DecayVertexXBefore->Fill(decayVertexV0[0]);
      fHistAntiV0DecayVertexYBefore->Fill(decayVertexV0[1]);
      fHistAntiV0DecayVertexZBefore->Fill(decayVertexV0[2]);
      fHistAntiV0TransverseRadiusBefore->Fill(transverseRadius);
      fHistAntiV0DCADaughtersBefore->Fill(dcaV0Dau);
    }
  }

  // Position of the decay vertex x, y & z
  if (std::fabs(decayVertexV0[0]) > fV0DecayVertexMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(9);
  else
    fHistAntiV0Cuts->Fill(9);
  if (std::fabs(decayVertexV0[1]) > fV0DecayVertexMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(10);
  else
    fHistAntiV0Cuts->Fill(10);
  if (std::fabs(decayVertexV0[2]) > fV0DecayVertexMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(11);
  else
    fHistAntiV0Cuts->Fill(11);

  // Transverse decay radius min & max
  if (transverseRadius < fV0RadiusMin) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(12);
  else
    fHistAntiV0Cuts->Fill(12);
  if (transverseRadius > fV0RadiusMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(13);
  else
    fHistAntiV0Cuts->Fill(13);

  // DCA of the daughter tracks at the decay vertex
  if (dcaV0Dau > fDaughterDCAMax) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(14);
  else
    fHistAntiV0Cuts->Fill(14);

  // Cos pointing angle
  if (point < fV0CosPAMin) return false;
  if (fPID == 1)
    fHistV0Cuts->Fill(15);
  else
    fHistAntiV0Cuts->Fill(15);

  if (fIsExtendedQA) {
    if (fPID == 1) {
      fHistV0CosPA->Fill(v0->Pt(), point);
      fHistV0CosPAAfter->Fill(point);
      fHistV0DecayVertexXAfter->Fill(decayVertexV0[0]);
      fHistV0DecayVertexYAfter->Fill(decayVertexV0[1]);
      fHistV0DecayVertexZAfter->Fill(decayVertexV0[2]);
      fHistV0TransverseRadiusAfter->Fill(transverseRadius);
      fHistV0DCADaughtersAfter->Fill(dcaV0Dau);
      fHistV0DCA->Fill(dcaPrim);
      fHistV0DecayLength->Fill(lenDecay);
    } else {
      fHistAntiV0CosPA->Fill(v0->Pt(), point);
      fHistAntiV0CosPAAfter->Fill(point);
      fHistAntiV0DecayVertexXAfter->Fill(decayVertexV0[0]);
      fHistAntiV0DecayVertexYAfter->Fill(decayVertexV0[1]);
      fHistAntiV0DecayVertexZAfter->Fill(decayVertexV0[2]);
      fHistAntiV0TransverseRadiusAfter->Fill(transverseRadius);
      fHistAntiV0DCADaughtersAfter->Fill(dcaV0Dau);
      fHistAntiV0DCA->Fill(dcaPrim);
      fHistAntiV0DecayLength->Fill(lenDecay);
    }
  }
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::SingleParticleQualityCuts(AliVTrack *track,
                                                const float dcaDaughterToPV) {
  const float pt = track->Pt();
  const int charge = track->Charge();
  const float eta = track->Eta();
  const short nClsTPC = track->GetTPCNcls();
  const float nCrossedRows = track->GetTPCClusterInfo(2, 1);
  const short nFindable = track->GetTPCNclsF();
  const float ratioFindable = nCrossedRows / static_cast<float>(nFindable);

  int qaHistoCounter = 0;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  if (fIsExtendedQA) {
    const int histoPrefix = (charge > 0) ? 0 : 1;
    if (fPID == 1) {
      fHistV0SingleParticleEtaBefore[histoPrefix]->Fill(eta);
      fHistV0SingleParticleNclsTPCBefore[histoPrefix]->Fill(nClsTPC);
      fHistV0SingleParticleNclsTPCFindableBefore[histoPrefix]->Fill(nFindable);
      fHistV0SingleParticleNclsTPCRatioFindableBefore[histoPrefix]->Fill(
          ratioFindable);
      fHistV0SingleParticleNcrossedTPCBefore[histoPrefix]->Fill(nCrossedRows);
      fHistV0SingleParticleDCAtoPVBefore[histoPrefix]->Fill(dcaDaughterToPV);
    } else if (fPID == -1) {
      fHistAntiV0SingleParticleEtaBefore[histoPrefix]->Fill(eta);
      fHistAntiV0SingleParticleNclsTPCBefore[histoPrefix]->Fill(nClsTPC);
      fHistAntiV0SingleParticleNclsTPCFindableBefore[histoPrefix]->Fill(
          nFindable);
      fHistAntiV0SingleParticleNclsTPCRatioFindableBefore[histoPrefix]->Fill(
          ratioFindable);
      fHistAntiV0SingleParticleNcrossedTPCBefore[histoPrefix]->Fill(
          nCrossedRows);
      fHistAntiV0SingleParticleDCAtoPVBefore[histoPrefix]->Fill(
          dcaDaughterToPV);
    }
  }

  // Max eta cut
  if (std::fabs(eta) > fEtaMax) return false;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  // TPC nCluster cut
  if (nClsTPC < fTPCclusterMin) return false;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  // TPC Crossed rows
  if (nCrossedRows < fTPCnCrossedRowsMin) return false;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  // TPC Ratio Findable
  if (ratioFindable < fTPCratioFindable) return false;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  // TPC nCls Findable
  if (nFindable < fTPCfindableMin) return false;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  // TPC nCls Shared
  const int nClsSharedTPC =
      (fInputEvent->IsA() == AliESDEvent::Class())
          ? (static_cast<AliESDtrack *>(track))->GetTPCnclsS()
          : (static_cast<AliAODTrack *>(track))->GetTPCnclsS();
  if (nClsSharedTPC > fTPCnSharedMax) return false;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  // Minimal distance of daughters to primary vertex
  if (dcaDaughterToPV < fDaughterDCAPV) return false;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  // Keep track of shared clusters
  int nClsSharedITS = 0;
  for (int i = 0; i < 6; ++i) {
    if (track->HasSharedPointOnITSLayer(i)) ++nClsSharedITS;
  }

  if (fIsExtendedQA) {
    const int histoPrefix = (charge > 0) ? 0 : 1;
    if (fPID == 1) {
      fHistV0SingleParticleNclsTPCShared[histoPrefix]->Fill(nClsSharedTPC);
      fHistV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
          0.f, nClsSharedTPC);
      fHistV0SingleParticleNclsITSShared[histoPrefix]->Fill(nClsSharedITS);
      fHistV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
          0.f, nClsSharedITS);
      if (track->HasPointOnITSLayer(0)) {
        fHistV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            1.f, nClsSharedTPC);
        fHistV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            1.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(1)) {
        fHistV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            2.f, nClsSharedTPC);
        fHistV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            2.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(2)) {
        fHistV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            3.f, nClsSharedTPC);
        fHistV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            3.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(3)) {
        fHistV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            4.f, nClsSharedTPC);
        fHistV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            4.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(4)) {
        fHistV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            5.f, nClsSharedTPC);
        fHistV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            5.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(5)) {
        fHistV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            6.f, nClsSharedTPC);
        fHistV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            6.f, nClsSharedITS);
      }
      if (track->GetTOFBunchCrossing() == 0) {
        fHistV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            7.f, nClsSharedTPC);
        fHistV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            7.f, nClsSharedITS);
      }
    } else if (fPID == -1) {
      fHistAntiV0SingleParticleNclsTPCShared[histoPrefix]->Fill(nClsSharedTPC);
      fHistAntiV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
          0.f, nClsSharedTPC);
      fHistAntiV0SingleParticleNclsITSShared[histoPrefix]->Fill(nClsSharedITS);
      fHistAntiV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
          0.f, nClsSharedITS);
      if (track->HasPointOnITSLayer(0)) {
        fHistAntiV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            1.f, nClsSharedTPC);
        fHistAntiV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            1.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(1)) {
        fHistAntiV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            2.f, nClsSharedTPC);
        fHistAntiV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            2.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(2)) {
        fHistAntiV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            3.f, nClsSharedTPC);
        fHistAntiV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            3.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(3)) {
        fHistAntiV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            4.f, nClsSharedTPC);
        fHistAntiV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            4.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(4)) {
        fHistAntiV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            5.f, nClsSharedTPC);
        fHistAntiV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            5.f, nClsSharedITS);
      }
      if (track->HasPointOnITSLayer(5)) {
        fHistAntiV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            6.f, nClsSharedTPC);
        fHistAntiV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            6.f, nClsSharedITS);
      }
      if (track->GetTOFBunchCrossing() == 0) {
        fHistAntiV0SingleParticleNclsTPCSharedTiming[histoPrefix]->Fill(
            7.f, nClsSharedTPC);
        fHistAntiV0SingleParticleNclsITSSharedTiming[histoPrefix]->Fill(
            7.f, nClsSharedITS);
      }
    }
  }

  // Pile-up cut
  if (fPileUpRejection && !track->HasPointOnITSLayer(0) &&
      !track->HasPointOnITSLayer(1) && !track->HasPointOnITSLayer(4) &&
      !track->HasPointOnITSLayer(5) && !(track->GetTOFBunchCrossing() == 0))
    return false;
  if (fPID == 1 && charge > 0)
    fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == 1 && charge < 0)
    fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge > 0)
    fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
  else if (fPID == -1 && charge < 0)
    fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);

  // For ESDs: TPC refit and Kinks
  if (fInputEvent->IsA() == AliESDEvent::Class()) {
    AliESDtrack *esdTrack = static_cast<AliESDtrack *>(track);
    if (!(esdTrack->GetStatus() & AliESDtrack::kTPCrefit)) return false;
    if (fPID == 1 && charge > 0)
      fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
    else if (fPID == 1 && charge < 0)
      fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
    else if (fPID == -1 && charge > 0)
      fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
    else if (fPID == -1 && charge < 0)
      fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
    if (esdTrack->GetKinkIndex(0) > 0) return false;
    if (fPID == 1 && charge > 0)
      fHistV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
    else if (fPID == 1 && charge < 0)
      fHistV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
    else if (fPID == -1 && charge > 0)
      fHistAntiV0SingleParticleCuts[0]->Fill(qaHistoCounter++);
    else if (fPID == -1 && charge < 0)
      fHistAntiV0SingleParticleCuts[1]->Fill(qaHistoCounter++);
  }

  if (fIsExtendedQA) {
    const int histoPrefix = (charge > 0) ? 0 : 1;
    if (fPID == 1) {
      fHistV0SingleParticlePt[histoPrefix]->Fill(pt);
      fHistV0SingleParticleEtaAfter[histoPrefix]->Fill(eta);
      fHistV0SingleParticleNclsTPCAfter[histoPrefix]->Fill(nClsTPC);
      fHistV0SingleParticleNclsTPCFindableAfter[histoPrefix]->Fill(nFindable);
      fHistV0SingleParticleNclsTPCRatioFindableAfter[histoPrefix]->Fill(
          ratioFindable);
      fHistV0SingleParticleNcrossedTPCAfter[histoPrefix]->Fill(nCrossedRows);
      fHistV0SingleParticleDCAtoPVAfter[histoPrefix]->Fill(dcaDaughterToPV);
    } else if (fPID == -1) {
      fHistAntiV0SingleParticlePt[histoPrefix]->Fill(pt);
      fHistAntiV0SingleParticleEtaAfter[histoPrefix]->Fill(eta);
      fHistAntiV0SingleParticleNclsTPCAfter[histoPrefix]->Fill(nClsTPC);
      fHistAntiV0SingleParticleNclsTPCFindableAfter[histoPrefix]->Fill(
          nFindable);
      fHistAntiV0SingleParticleNclsTPCRatioFindableAfter[histoPrefix]->Fill(
          ratioFindable);
      fHistAntiV0SingleParticleNcrossedTPCAfter[histoPrefix]->Fill(
          nCrossedRows);
      fHistAntiV0SingleParticleDCAtoPVAfter[histoPrefix]->Fill(dcaDaughterToPV);
    }
  }
  return true;
}

//____________________________________________________________________________________________________
template <typename T>
bool AliSigma0V0Cuts::LambdaSelection(const T *v0, float massK0,
                                      float massLambda) {
  const float armAlpha = v0->AlphaV0();
  const float armQt = v0->PtArmV0();

  if (fPID == 1) {
    fHistV0LambdaMassK0Rej->Fill(massLambda);
    fHistV0K0Mass->Fill(massK0);
  } else {
    fHistAntiV0LambdaMassK0Rej->Fill(massLambda);
    fHistAntiV0K0Mass->Fill(massK0);
  }

  // K0 rejection cut
  if (!fArmenterosCut && massK0 > fK0RejectionLow && massK0 < fK0RejectionUp)
    return false;

  // Armenteros cut
  if (fArmenterosCut) {
    if (armQt > fArmenterosQtUp || armQt < fArmenterosQtLow) return false;
    float prefactorAlpha = (fPID == 1) ? 1.f : -1.f;  // for anti-particles
                                                      // negative to compensate
                                                      // for sign change of qt
    if (armAlpha * prefactorAlpha > fArmenterosAlphaUp ||
        armAlpha * prefactorAlpha < fArmenterosAlphaLow)
      return false;
  }

  float rap = ComputeRapidity(v0->Pt(), v0->Pz(), massLambda);
  int rapBin = GetRapidityBin(rap);

  if (fPID == 1) {
    fHistV0ArmenterosAfter->Fill(armAlpha, armQt);
    fHistV0LambdaMass->Fill(massLambda);
    fHistV0LambdaMassPt->Fill(v0->Pt(), massLambda);
    if (rapBin > -1) fHistV0LambdaPtY[rapBin]->Fill(v0->Pt(), massLambda);
    fHistV0K0MassAfter->Fill(massK0);
    fHistV0EtaPhi->Fill(v0->Eta(), v0->Phi());
  } else {
    fHistAntiV0ArmenterosAfter->Fill(armAlpha, armQt);
    fHistAntiV0LambdaMass->Fill(massLambda);
    fHistAntiV0LambdaMassPt->Fill(v0->Pt(), massLambda);
    if (rapBin > -1) fHistAntiV0LambdaPtY[rapBin]->Fill(v0->Pt(), massLambda);
    fHistAntiV0K0MassAfter->Fill(massK0);
    fHistAntiV0EtaPhi->Fill(v0->Eta(), v0->Phi());
  }

  // Lambda selection cut
  if (massLambda < fLambdaSelectionLow || massLambda > fLambdaSelectionUp)
    return false;

  if (fPID == 1) {
    fHistV0LambdaPt->Fill(v0->Pt());
    fHistV0Cuts->Fill(16);
  } else {
    fHistAntiV0LambdaPt->Fill(v0->Pt());
    fHistAntiV0Cuts->Fill(16);
  }

  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::CheckIfRealV0(const AliAODv0 *v0, int PDGmother,
                                    int PDGdaugh1, int PDGdaugh2) const {
  // Check whether V0 is real or fake

  bool realV0 = false;
  TClonesArray *mcarray = static_cast<TClonesArray *>(
      fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcarray) return false;

  Int_t daugh[2] = {PDGdaugh1, PDGdaugh2};
  Int_t label = v0->MatchToMC(PDGmother, mcarray, 2, daugh);

  if (label >= 0) realV0 = true;

  return realV0;
}

//____________________________________________________________________________________________________
void AliSigma0V0Cuts::ProcessMC() const {
  // Loop over the MC tracks
  for (int iPart = 1; iPart < (fMCEvent->GetNumberOfTracks()); iPart++) {
    AliMCParticle *mcParticle =
        static_cast<AliMCParticle *>(fMCEvent->GetTrack(iPart));
    if (!mcParticle) continue;
    if (mcParticle->PdgCode() == 3122) {
      fHistMCTruthV0Pt->Fill(mcParticle->Pt());
      fHistMCTruthV0PtY->Fill(
          ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(), mcParticle->M()),
          mcParticle->Pt());
      fHistMCTruthV0PtEta->Fill(mcParticle->Eta(), mcParticle->Pt());

      if (IsLambdaProtonPion(mcParticle)) {
        fHistMCTruthV0ProtonPionPt->Fill(mcParticle->Pt());
        fHistMCTruthV0ProtonPionPtY->Fill(
            ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(),
                            mcParticle->M()),
            mcParticle->Pt());
        fHistMCTruthV0ProtonPionPtEta->Fill(mcParticle->Eta(),
                                            mcParticle->Pt());
      }
    }
    if (mcParticle->PdgCode() == -3122) {
      fHistMCTruthAntiV0Pt->Fill(mcParticle->Pt());
      fHistMCTruthAntiV0PtY->Fill(
          ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(), mcParticle->M()),
          mcParticle->Pt());
      fHistMCTruthAntiV0PtEta->Fill(mcParticle->Eta(), mcParticle->Pt());

      if (IsLambdaProtonPion(mcParticle)) {
        fHistMCTruthAntiV0ProtonPionPt->Fill(mcParticle->Pt());
        fHistMCTruthAntiV0ProtonPionPtY->Fill(
            ComputeRapidity(mcParticle->Pt(), mcParticle->Pz(),
                            mcParticle->M()),
            mcParticle->Pt());
        fHistMCTruthAntiV0ProtonPionPtEta->Fill(mcParticle->Eta(),
                                                mcParticle->Pt());
      }
    }
  }
}

//____________________________________________________________________________________________________
bool AliSigma0V0Cuts::IsLambdaProtonPion(AliMCParticle *particle) const {
  if (std::abs(particle->PdgCode()) != 3122) return false;

  // looking for proton and pion daughers
  AliMCParticle *ePion = nullptr;
  AliMCParticle *eProton = nullptr;
  if (particle->GetNDaughters() >= 2) {
    for (int daughterIndex = particle->GetFirstDaughter();
         daughterIndex <= particle->GetLastDaughter(); ++daughterIndex) {
      if (daughterIndex < 0) continue;
      AliMCParticle *tmpDaughter =
          static_cast<AliMCParticle *>(fMCEvent->GetTrack(daughterIndex));
      const int pdgCode = tmpDaughter->PdgCode();
      if (particle->PdgCode() > 0) {
        if (pdgCode == -211)
          ePion = tmpDaughter;
        else if (pdgCode == 2212)
          eProton = tmpDaughter;
      } else {
        if (pdgCode == 211)
          ePion = tmpDaughter;
        else if (pdgCode == -2212)
          eProton = tmpDaughter;
      }
    }
  }

  if (ePion == nullptr || eProton == nullptr) return false;
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
void AliSigma0V0Cuts::InitCutHistograms(const char *appendix) {
  std::cout << "============================\n"
            << " V0 CUT CONFIGURATION \n"
            << " Pile-up rej  " << fPileUpRejection << "\n"
            << " On-fly       " << fV0OnFly << "\n"
            << " p_T min      " << fV0PtMin << "\n"
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
            << " daughter dca " << fDaughterDCAMax << "\n"
            << " daugher pv   " << fDaughterDCAPV << "\n"
            << " K0 rej       " << fK0RejectionLow << " " << fK0RejectionUp
            << "\n"
            << " Lambda sel   " << fLambdaSelectionLow << " "
            << fLambdaSelectionUp << "\n"
            << "============================\n";

  TH1::AddDirectory(kFALSE);

  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    fHistograms->SetName(Form("V0Cut_%s_QA", appendix));
  }

  if (fHistogramsV0 != nullptr) {
    delete fHistogramsV0;
    fHistogramsV0 = nullptr;
  }
  if (fHistogramsV0 == nullptr) {
    fHistogramsV0 = new TList();
    fHistogramsV0->SetOwner(kTRUE);
    fHistogramsV0->SetName("V0_QA");
  }

  if (fHistogramsV0Pos != nullptr) {
    delete fHistogramsV0Pos;
    fHistogramsV0Pos = nullptr;
  }
  if (fHistogramsV0Pos == nullptr) {
    fHistogramsV0Pos = new TList();
    fHistogramsV0Pos->SetOwner(kTRUE);
    fHistogramsV0Pos->SetName("V0_SingleParticle_Pos_QA");
  }

  if (fHistogramsV0Neg != nullptr) {
    delete fHistogramsV0Neg;
    fHistogramsV0Neg = nullptr;
  }
  if (fHistogramsV0Neg == nullptr) {
    fHistogramsV0Neg = new TList();
    fHistogramsV0Neg->SetOwner(kTRUE);
    fHistogramsV0Neg->SetName("V0_SingleParticle_Neg_QA");
  }

  if (fHistogramsAntiV0 != nullptr) {
    delete fHistogramsAntiV0;
    fHistogramsAntiV0 = nullptr;
  }
  if (fHistogramsAntiV0 == nullptr) {
    fHistogramsAntiV0 = new TList();
    fHistogramsAntiV0->SetOwner(kTRUE);
    fHistogramsAntiV0->SetName("AntiV0_QA");
  }

  if (fHistogramsAntiV0Pos != nullptr) {
    delete fHistogramsAntiV0Pos;
    fHistogramsAntiV0Pos = nullptr;
  }
  if (fHistogramsAntiV0Pos == nullptr) {
    fHistogramsAntiV0Pos = new TList();
    fHistogramsAntiV0Pos->SetOwner(kTRUE);
    fHistogramsAntiV0Pos->SetName("AntiV0_SingleParticle_Pos_QA");
  }

  if (fHistogramsAntiV0Neg != nullptr) {
    delete fHistogramsAntiV0Neg;
    fHistogramsAntiV0Neg = nullptr;
  }
  if (fHistogramsAntiV0Neg == nullptr) {
    fHistogramsAntiV0Neg = new TList();
    fHistogramsAntiV0Neg->SetOwner(kTRUE);
    fHistogramsAntiV0Neg->SetName("AntiV0_SingleParticle_Neg_QA");
  }

  fHistCuts = new TProfile("fHistCuts", ";;Cut value", 24, 0, 24);
  fHistCuts->GetXaxis()->SetBinLabel(1, "V0 on fly");
  fHistCuts->GetXaxis()->SetBinLabel(2, "#it{p}_{T} min");
  fHistCuts->GetXaxis()->SetBinLabel(3, "cos#alpha min");
  fHistCuts->GetXaxis()->SetBinLabel(4, "max r");
  fHistCuts->GetXaxis()->SetBinLabel(5, "min r");
  fHistCuts->GetXaxis()->SetBinLabel(6, "Decay vtx max");
  fHistCuts->GetXaxis()->SetBinLabel(7, "n#sigma PID");
  fHistCuts->GetXaxis()->SetBinLabel(8, "TPC nCls min");
  fHistCuts->GetXaxis()->SetBinLabel(9, "TPC nCrossed min");
  fHistCuts->GetXaxis()->SetBinLabel(10, "TPC Ratio Findable");
  fHistCuts->GetXaxis()->SetBinLabel(11, "TPC nCls Findable min");
  fHistCuts->GetXaxis()->SetBinLabel(12, "TPC nShared max");
  fHistCuts->GetXaxis()->SetBinLabel(13, "#eta max");
  fHistCuts->GetXaxis()->SetBinLabel(14, "Daughter DCA max");
  fHistCuts->GetXaxis()->SetBinLabel(15, "Daughter DCA pvtx min");
  fHistCuts->GetXaxis()->SetBinLabel(16, "K^{0} rejection low");
  fHistCuts->GetXaxis()->SetBinLabel(17, "K^{0} rejection up");
  fHistCuts->GetXaxis()->SetBinLabel(18, "#Lambda selection low");
  fHistCuts->GetXaxis()->SetBinLabel(19, "#Lambda selection up");
  fHistCuts->GetXaxis()->SetBinLabel(20, "Pile-up rejection");
  fHistCuts->GetXaxis()->SetBinLabel(21, "Armenteros q_{T} low");
  fHistCuts->GetXaxis()->SetBinLabel(22, "Armenteros q_{T} up");
  fHistCuts->GetXaxis()->SetBinLabel(23, "Armenteros #alpha low");
  fHistCuts->GetXaxis()->SetBinLabel(24, "Armenteros #alpha up");
  fHistograms->Add(fHistCuts);

  fHistCuts->Fill(0.f, static_cast<double>(fV0OnFly));
  fHistCuts->Fill(1.f, fV0PtMin);
  fHistCuts->Fill(2.f, fV0CosPAMin);
  fHistCuts->Fill(3.f, fV0RadiusMax);
  fHistCuts->Fill(4.f, fV0RadiusMin);
  fHistCuts->Fill(5.f, fV0DecayVertexMax);
  fHistCuts->Fill(6.f, fPIDnSigma);
  fHistCuts->Fill(7.f, fTPCclusterMin);
  fHistCuts->Fill(8.f, fTPCnCrossedRowsMin);
  fHistCuts->Fill(9.f, fTPCratioFindable);
  fHistCuts->Fill(10.f, fTPCfindableMin);
  fHistCuts->Fill(11.f, fTPCnSharedMax);
  fHistCuts->Fill(12.f, fEtaMax);
  fHistCuts->Fill(13.f, fDaughterDCAMax);
  fHistCuts->Fill(14.f, fDaughterDCAPV);
  fHistCuts->Fill(15.f, fK0RejectionLow);
  fHistCuts->Fill(16.f, fK0RejectionUp);
  fHistCuts->Fill(17.f, fLambdaSelectionLow);
  fHistCuts->Fill(18.f, fLambdaSelectionUp);
  fHistCuts->Fill(19.f, static_cast<double>(fPileUpRejection));
  fHistCuts->Fill(20.f, fArmenterosQtLow);
  fHistCuts->Fill(21.f, fArmenterosQtUp);
  fHistCuts->Fill(22.f, fArmenterosAlphaLow);
  fHistCuts->Fill(23.f, fArmenterosAlphaUp);

  fHistV0Cuts = new TH1F("fHistV0Cuts", ";;Entries", 18, 0, 18);
  fHistV0Cuts->GetXaxis()->SetBinLabel(1, "V0");
  fHistV0Cuts->GetXaxis()->SetBinLabel(2, "Two prongs");
  fHistV0Cuts->GetXaxis()->SetBinLabel(3, "Two daughters");
  fHistV0Cuts->GetXaxis()->SetBinLabel(4, "Neutral");
  fHistV0Cuts->GetXaxis()->SetBinLabel(5, "On-fly");
  fHistV0Cuts->GetXaxis()->SetBinLabel(6, "#it{p}_{T, min}");
  fHistV0Cuts->GetXaxis()->SetBinLabel(7, "Daughter tracks");
  fHistV0Cuts->GetXaxis()->SetBinLabel(8, "PID");
  fHistV0Cuts->GetXaxis()->SetBinLabel(9, "SingleParticle QA");
  fHistV0Cuts->GetXaxis()->SetBinLabel(10, "Decay vtx x");
  fHistV0Cuts->GetXaxis()->SetBinLabel(11, "Decay vtx y");
  fHistV0Cuts->GetXaxis()->SetBinLabel(12, "Decay vtx z");
  fHistV0Cuts->GetXaxis()->SetBinLabel(13, "Radius min");
  fHistV0Cuts->GetXaxis()->SetBinLabel(14, "Radius max");
  fHistV0Cuts->GetXaxis()->SetBinLabel(15, "Daughter DCA");
  fHistV0Cuts->GetXaxis()->SetBinLabel(16, "cos #alpha");
  fHistV0Cuts->GetXaxis()->SetBinLabel(17, "#Lambda selection");
  fHistogramsV0->Add(fHistV0Cuts);

  fHistNV0 =
      new TH1F("fHistNV0", ";Number of V0 candidates; Entries", 25, 0, 25);
  fHistogramsV0->Add(fHistNV0);

  fHistV0LambdaMass =
      new TH1F("fHistV0LambdaMass",
               "; Invariant mass p#pi hypothesis [GeV/#it{c}^{2}]; Entries",
               4000, 0., 2.);
  fHistV0LambdaPt = new TH1F("fHistV0LambdaPt",
                             "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistV0LambdaMassPt = new TH2F("fHistV0LambdaMassPt",
                                 "; #it{p}_{T} [GeV/#it{c}];Invariant mass "
                                 "p#pi hypothesis [GeV/#it{c}^{2}]",
                                 1000, 0, 15, 400, 1., 1.2);
  fHistV0LambdaMassK0Rej = new TH1F("fHistV0LambdaMassK0Rej",
                                    "; Invariant mass p#pi hypothesis w/o "
                                    "K^{0} rejection [GeV/#it{c}^{2}]; Entries",
                                    1000, 0., 2.);
  fHistV0K0Mass =
      new TH1F("fHistV0K0Mass",
               "; Invariant mass #pi#pi hypothesis [GeV/#it{c}^{2}]; Entries",
               1000, 0., 2.);
  fHistV0K0MassAfter =
      new TH1F("fHistV0K0MassAfter",
               "; Invariant mass #pi#pi hypothesis [GeV/#it{c}^{2}]; Entries",
               1000, 0., 2.);
  fHistV0CosPA =
      new TH2F("fHistV0CosPA", "; #it{p}_{T} [GeV/#it{c}]; cos(#alpha)", 1000,
               0, 10, 1000, 0.9, 1);
  fHistV0EtaPhi = new TH2F("fHistV0EtaPhi", "; #eta; #phi", 200, -1, 1, 200, 0,
                           2 * TMath::Pi());
  fHistogramsV0->Add(fHistV0LambdaMass);
  fHistogramsV0->Add(fHistV0LambdaPt);
  fHistogramsV0->Add(fHistV0LambdaMassPt);
  fHistogramsV0->Add(fHistV0LambdaMassK0Rej);
  fHistogramsV0->Add(fHistV0K0Mass);
  fHistogramsV0->Add(fHistV0K0MassAfter);
  fHistogramsV0->Add(fHistV0CosPA);
  fHistogramsV0->Add(fHistV0EtaPhi);

  std::vector<float> rapBins = {{-1.,  -0.9, -0.8, -0.7, -0.6, -0.5, -0.4,
                                 -0.3, -0.2, -0.1, 0.f,  0.1,  0.2,  0.3,
                                 0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.f}};
  for (int i = 0; i < static_cast<int>(rapBins.size() - 1); ++i) {
    fHistV0LambdaPtY[i] =
        new TH2F(Form("fHistV0LambdaPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
                 Form("%.2f < y < %.2f ; #it{p}_{T} [GeV/#it{c}]; Invariant "
                      "mass p#pi hypothesis [GeV/#it{c}^{2}]",
                      rapBins[i], rapBins[i + 1]),
                 500, 0, 10, 400, 1., 1.2);
    fHistogramsV0->Add(fHistV0LambdaPtY[i]);
  }

  fHistV0SingleParticleCuts[0] =
      new TH1F("fHistV0SingleParticleCuts_pos", ";;Entries", 11, 0, 11);
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(2, "#eta");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(3, "nCls TPC");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(4, "nCrossed Rows TPC");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(5,
                                                        "Ratio Findable TPC");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(6, "nCls Findable");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(7, "nCls Shared max");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(8,
                                                        "Daughter DCA to PV");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(9, "Pile-up cut");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(10, "ESD: TPC refit");
  fHistV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(11, "ESD: Kink");
  fHistogramsV0Pos->Add(fHistV0SingleParticleCuts[0]);

  fHistV0SingleParticleCuts[1] =
      new TH1F("fHistV0SingleParticleCuts_neg", ";;Entries", 11, 0, 11);
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(1, "Daughter tracks");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(2, "#eta");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(3, "nCls TPC");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(4, "nCrossed Rows TPC");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(5,
                                                        "Ratio Findable TPC");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(6, "nCls Findable");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(7, "nCls Shared max");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(8,
                                                        "Daughter DCA to PV");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(9, "Pile-up cut");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(10, "ESD: TPC refit");
  fHistV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(11, "ESD: Kink");
  fHistogramsV0Neg->Add(fHistV0SingleParticleCuts[1]);

  fHistAntiV0Cuts = new TH1F("fHistAntiV0Cuts", ";;Entries", 18, 0, 18);
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(1, "V0");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(2, "Two prongs");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(3, "Two daughters");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(4, "Neutral");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(5, "On-fly");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(6, "#it{p}_{T, min}");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(7, "Daughter tracks");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(8, "PID");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(9, "SingleParticle QA");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(10, "Decay vtx x");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(11, "Decay vtx y");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(12, "Decay vtx z");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(13, "Radius min");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(14, "Radius max");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(15, "Daughter DCA");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(16, "cos #alpha");
  fHistAntiV0Cuts->GetXaxis()->SetBinLabel(17, "#Lambda selection");
  fHistogramsAntiV0->Add(fHistAntiV0Cuts);

  fHistNAntiV0 = new TH1F("fHistNAntiV0",
                          ";Number of Anti-V0 candidates; Entries", 25, 0, 25);
  fHistogramsAntiV0->Add(fHistNAntiV0);

  fHistAntiV0LambdaMass =
      new TH1F("fHistAntiV0LambdaMass",
               "; Invariant mass p#pi hypothesis [GeV/#it{c}^{2}]; Entries",
               4000, 0., 2.);
  fHistAntiV0LambdaPt = new TH1F(
      "fHistAntiV0LambdaPt", "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
  fHistAntiV0LambdaMassPt = new TH2F("fHistAntiV0LambdaMassPt",
                                     "; #it{p}_{T} [GeV/#it{c}];Invariant mass "
                                     "p#pi hypothesis [GeV/#it{c}^{2}]",
                                     1000, 0, 15, 400, 1., 1.2);
  fHistAntiV0LambdaMassK0Rej =
      new TH1F("fHistAntiV0LambdaMassK0Rej",
               "; Invariant mass p#pi hypothesis w/o K^{0} rejection "
               "[GeV/#it{c}^{2}]; Entries",
               1000, 0., 2.);
  fHistAntiV0K0Mass =
      new TH1F("fHistAntiV0K0Mass",
               "; Invariant mass #pi#pi hypothesis [GeV/#it{c}^{2}]; Entries",
               1000, 0., 2.);
  fHistAntiV0K0MassAfter =
      new TH1F("fHistAntiV0K0MassAfter",
               "; Invariant mass #pi#pi hypothesis [GeV/#it{c}^{2}]; Entries",
               1000, 0., 2.);
  fHistAntiV0CosPA =
      new TH2F("fHistAntiV0CosPA", "; #it{p}_{T} [GeV/#it{c}]; cos(#alpha)",
               1000, 0, 10, 1000, 0.9, 1);
  fHistAntiV0EtaPhi = new TH2F("fHistAntiV0EtaPhi", "; #eta; #phi", 200, -1, 1,
                               200, 0, 2 * TMath::Pi());
  fHistogramsAntiV0->Add(fHistAntiV0LambdaMass);
  fHistogramsAntiV0->Add(fHistAntiV0LambdaPt);
  fHistogramsAntiV0->Add(fHistAntiV0LambdaMassPt);
  fHistogramsAntiV0->Add(fHistAntiV0LambdaMassK0Rej);
  fHistogramsAntiV0->Add(fHistAntiV0K0Mass);
  fHistogramsAntiV0->Add(fHistAntiV0K0MassAfter);
  fHistogramsAntiV0->Add(fHistAntiV0CosPA);
  fHistogramsAntiV0->Add(fHistAntiV0EtaPhi);

  for (int i = 0; i < static_cast<int>(rapBins.size() - 1); ++i) {
    fHistAntiV0LambdaPtY[i] = new TH2F(
        Form("fHistAntiV0LambdaPtY_%.2f_%.2f", rapBins[i], rapBins[i + 1]),
        Form("%.2f < y < %.2f ; #it{p}_{T} [GeV/#it{c}]; Invariant mass p#pi "
             "hypothesis [GeV/#it{c}^{2}]",
             rapBins[i], rapBins[i + 1]),
        500, 0, 10, 400, 1., 1.2);
    fHistogramsAntiV0->Add(fHistAntiV0LambdaPtY[i]);
  }

  fHistAntiV0SingleParticleCuts[0] =
      new TH1F("fHistAntiV0SingleParticleCuts_pos", ";;Entries", 11, 0, 11);
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(1,
                                                            "Daughter tracks");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(2, "#eta");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(3, "nCls TPC");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(
      4, "nCrossed Rows TPC");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(
      5, "Ratio Findable TPC");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(6, "nCls Findable");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(7,
                                                            "nCls Shared max");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(
      8, "Daughter DCA to PV");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(9, "Pile-up cut");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(10,
                                                            "ESD: TPC refit");
  fHistAntiV0SingleParticleCuts[0]->GetXaxis()->SetBinLabel(11, "ESD: Kink");
  fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleCuts[0]);

  fHistAntiV0SingleParticleCuts[1] =
      new TH1F("fHistAntiV0SingleParticleCuts_pos", ";;Entries", 11, 0, 11);
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(1,
                                                            "Daughter tracks");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(2, "#eta");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(3, "nCls TPC");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(
      4, "nCrossed Rows TPC");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(
      5, "Ratio Findable TPC");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(6, "nCls Findable");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(7,
                                                            "nCls Shared max");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(
      8, "Daughter DCA to PV");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(9, "Pile-up cut");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(10,
                                                            "ESD: TPC refit");
  fHistAntiV0SingleParticleCuts[1]->GetXaxis()->SetBinLabel(11, "ESD: Kink");
  fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleCuts[1]);

  if (fIsMC) {
    if (fHistogramsV0MC != nullptr) {
      delete fHistogramsV0MC;
      fHistogramsV0MC = nullptr;
    }
    if (fHistogramsV0MC == nullptr) {
      fHistogramsV0MC = new TList();
      fHistogramsV0MC->SetOwner(kTRUE);
      fHistogramsV0MC->SetName("V0Cut_MC");
    }

    if (fHistogramsAntiV0MC != nullptr) {
      delete fHistogramsAntiV0MC;
      fHistogramsAntiV0MC = nullptr;
    }
    if (fHistogramsAntiV0MC == nullptr) {
      fHistogramsAntiV0MC = new TList();
      fHistogramsAntiV0MC->SetOwner(kTRUE);
      fHistogramsAntiV0MC->SetName("AntiV0Cut_MC");
    }

    fHistMCTruthV0Pt = new TH1F(
        "fHistMCTruthV0Pt", "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthV0PtY =
        new TH2F("fHistMCTruthV0PtY", "; y; #it{p}_{T} [GeV/#it{c}]", 500, -10,
                 10, 500, 0, 10);
    fHistMCTruthV0PtEta =
        new TH2F("fHistMCTruthV0PtEta", "; #eta; #it{p}_{T} [GeV/#it{c}]", 500,
                 -10, 10, 500, 0, 10);
    fHistMCTruthV0ProtonPionPt =
        new TH1F("fHistMCTruthV0ProtonPionPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthV0ProtonPionPtY =
        new TH2F("fHistMCTruthV0ProtonPionPtY", "; y; #it{p}_{T} [GeV/#it{c}]",
                 500, -10, 10, 500, 0, 10);
    fHistMCTruthV0ProtonPionPtEta =
        new TH2F("fHistMCTruthV0ProtonPionPtEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistogramsV0MC->Add(fHistMCTruthV0Pt);
    fHistogramsV0MC->Add(fHistMCTruthV0PtY);
    fHistogramsV0MC->Add(fHistMCTruthV0PtEta);
    fHistogramsV0MC->Add(fHistMCTruthV0ProtonPionPt);
    fHistogramsV0MC->Add(fHistMCTruthV0ProtonPionPtY);
    fHistogramsV0MC->Add(fHistMCTruthV0ProtonPionPtEta);

    fHistMCTruthAntiV0Pt =
        new TH1F("fHistMCTruthAntiV0Pt", "; #it{p}_{T} [GeV/#it{c}]; Entries",
                 500, 0, 10);
    fHistMCTruthAntiV0PtY =
        new TH2F("fHistMCTruthAntiV0PtY", "; y; #it{p}_{T} [GeV/#it{c}]", 1000,
                 -10, 10, 1000, 0, 10);
    fHistMCTruthAntiV0PtEta =
        new TH2F("fHistMCTruthAntiV0PtEta", "; #eta; #it{p}_{T} [GeV/#it{c}]",
                 500, -10, 10, 500, 0, 10);
    fHistMCTruthAntiV0ProtonPionPt =
        new TH1F("fHistMCTruthAntiV0ProtonPionPt",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCTruthAntiV0ProtonPionPtY =
        new TH2F("fHistMCTruthAntiV0ProtonPionPtY",
                 "; y; #it{p}_{T} [GeV/#it{c}]", 1000, -10, 10, 1000, 0, 10);
    fHistMCTruthAntiV0ProtonPionPtEta =
        new TH2F("fHistMCTruthAntiV0ProtonPionPtEta",
                 "; #eta; #it{p}_{T} [GeV/#it{c}]", 500, -10, 10, 500, 0, 10);
    fHistogramsAntiV0MC->Add(fHistMCTruthAntiV0Pt);
    fHistogramsAntiV0MC->Add(fHistMCTruthAntiV0PtY);
    fHistogramsAntiV0MC->Add(fHistMCTruthAntiV0PtEta);
    fHistogramsAntiV0MC->Add(fHistMCTruthAntiV0ProtonPionPt);
    fHistogramsAntiV0MC->Add(fHistMCTruthAntiV0ProtonPionPtY);
    fHistogramsAntiV0MC->Add(fHistMCTruthAntiV0ProtonPionPtEta);

    fHistogramsV0->Add(fHistogramsV0MC);
    fHistogramsAntiV0->Add(fHistogramsAntiV0MC);
  }

  if (fIsExtendedQA) {
    if (fHistogramsV0Before != nullptr) {
      delete fHistogramsV0Before;
      fHistogramsV0Before = nullptr;
    }
    if (fHistogramsV0Before == nullptr) {
      fHistogramsV0Before = new TList();
      fHistogramsV0Before->SetOwner(kTRUE);
      fHistogramsV0Before->SetName("Before");
    }

    if (fHistogramsV0After != nullptr) {
      delete fHistogramsV0After;
      fHistogramsV0After = nullptr;
    }
    if (fHistogramsV0After == nullptr) {
      fHistogramsV0After = new TList();
      fHistogramsV0After->SetOwner(kTRUE);
      fHistogramsV0After->SetName("After");
    }

    fHistV0DecayVertexXBefore =
        new TH1F("fHistV0DecayVertexXBefore", "; Decay vertex x [cm]; Entries",
                 1000, 0, 200);
    fHistV0DecayVertexYBefore =
        new TH1F("fHistV0DecayVertexYBefore", "; Decay vertex y [cm]; Entries",
                 1000, 0, 200);
    fHistV0DecayVertexZBefore =
        new TH1F("fHistV0DecayVertexZBefore", "; Decay vertex z [cm]; Entries",
                 1000, 0, 200);
    fHistV0DecayVertexXAfter =
        new TH1F("fHistV0DecayVertexXAfter", "; Decay vertex x [cm]; Entries",
                 1000, 0, 200);
    fHistV0DecayVertexYAfter =
        new TH1F("fHistV0DecayVertexYAfter", "; Decay vertex y [cm]; Entries",
                 1000, 0, 200);
    fHistV0DecayVertexZAfter =
        new TH1F("fHistV0DecayVertexZAfter", "; Decay vertex z [cm]; Entries",
                 1000, 0, 200);
    fHistV0TransverseRadiusBefore =
        new TH1F("fHistV0TransverseRadiusBefore",
                 "; Transverse radius [cm]; Entries", 1000, 0, 200);
    fHistV0TransverseRadiusAfter =
        new TH1F("fHistV0TransverseRadiusAfter",
                 "; Transverse radius [cm]; Entries", 1000, 0, 200);
    fHistV0CosPABefore =
        new TH1F("fHistV0CosPABefore", "; cos(#alpha); Entries", 1000, 0.5, 1);
    fHistV0CosPAAfter =
        new TH1F("fHistV0CosPAAfter", "; cos(#alpha); Entries", 1000, 0.95, 1);
    fHistV0DCADaughtersBefore =
        new TH1F("fHistV0DCADaughtersBefore",
                 "; Daughter DCA at decay vertex [cm]; Entries", 1000, 0, 10);
    fHistV0DCADaughtersAfter =
        new TH1F("fHistV0DCADaughtersAfter",
                 "; Daughter DCA at decay vertex [cm]; Entries", 1000, 0, 10);
    fHistV0DCA =
        new TH1F("fHistV0DCA", "; DCA to PV [cm]; Entries", 1000, 0, 10);
    fHistV0DecayLength = new TH1F("fHistV0DecayLength",
                                  "; Decay length [cm]; Entries", 1000, 0, 200);
    fHistV0ArmenterosBefore = new TH2F(
        "fHistV0ArmenterosBefore", " ; #alpha; #it{q}_{T} p#pi [GeV/#it{c}]",
        500, -1, 1, 500, 0, 1);
    fHistV0ArmenterosAfter = new TH2F("fHistV0ArmenterosAfter",
                                      " ; #alpha; #it{q}_{T} p#pi [GeV/#it{c}]",
                                      500, -1, 1, 500, 0, 1);

    fHistogramsV0Before->Add(fHistV0DecayVertexXBefore);
    fHistogramsV0Before->Add(fHistV0DecayVertexYBefore);
    fHistogramsV0Before->Add(fHistV0DecayVertexZBefore);
    fHistogramsV0Before->Add(fHistV0TransverseRadiusBefore);
    fHistogramsV0Before->Add(fHistV0CosPABefore);
    fHistogramsV0Before->Add(fHistV0DCADaughtersBefore);
    fHistogramsV0After->Add(fHistV0DecayVertexXAfter);
    fHistogramsV0After->Add(fHistV0DecayVertexYAfter);
    fHistogramsV0After->Add(fHistV0DecayVertexZAfter);
    fHistogramsV0After->Add(fHistV0TransverseRadiusAfter);
    fHistogramsV0After->Add(fHistV0CosPAAfter);
    fHistogramsV0After->Add(fHistV0DCADaughtersAfter);
    fHistogramsV0->Add(fHistV0DCA);
    fHistogramsV0->Add(fHistV0DecayLength);
    fHistogramsV0Before->Add(fHistV0ArmenterosBefore);
    fHistogramsV0After->Add(fHistV0ArmenterosAfter);

    fHistV0SingleParticlePt[0] = new TH1F("fHistV0SingleParticlePt_pos",
                                          ";#it{p}_{T}; Entries", 1000, 0, 10);
    fHistV0SingleParticlePt[1] = new TH1F("fHistV0SingleParticlePt_neg",
                                          ";#it{p}_{T}; Entries", 1000, 0, 10);
    fHistV0SingleParticleEtaBefore[0] =
        new TH1F("fHistV0SingleParticleEtaBefore_pos",
                 "; #eta (before); Entries", 1000, -2, 2);
    fHistV0SingleParticleEtaBefore[1] =
        new TH1F("fHistV0SingleParticleEtaBefore_neg",
                 "; #eta (before); Entries", 1000, -2, 2);
    fHistV0SingleParticleEtaAfter[0] =
        new TH1F("fHistV0SingleParticleEtaAfter_pos", "; #eta (after); Entries",
                 1000, -2, 2);
    fHistV0SingleParticleEtaAfter[1] =
        new TH1F("fHistV0SingleParticleEtaAfter_neg", "; #eta (after); Entries",
                 1000, -2, 2);
    fHistV0SingleParticleNclsTPCBefore[0] =
        new TH1F("fHistV0SingleParticleNclsTPCBefore_pos",
                 "; # cls TPC (before); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCBefore[1] =
        new TH1F("fHistV0SingleParticleNclsTPCBefore_neg",
                 "; # cls TPC (before); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCAfter[0] =
        new TH1F("fHistV0SingleParticleNclsTPCAfter_pos",
                 "; # cls TPC (after); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCAfter[1] =
        new TH1F("fHistV0SingleParticleNclsTPCAfter_neg",
                 "; # cls TPC (after); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCFindableBefore[0] =
        new TH1F("fHistV0SingleParticleNclsTPCFindableBefore_pos",
                 "; # cls TPC findable (before); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCFindableBefore[1] =
        new TH1F("fHistV0SingleParticleNclsTPCFindableBefore_neg",
                 "; # cls TPC findable (before); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCFindableAfter[0] =
        new TH1F("fHistV0SingleParticleNclsTPCFindableAfter_pos",
                 "; # cls TPC findable (after); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCFindableAfter[1] =
        new TH1F("fHistV0SingleParticleNclsTPCFindableAfter_neg",
                 "; # cls TPC findable (after); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCRatioFindableBefore[0] =
        new TH1F("fHistV0SingleParticleNclsTPCRatioFindableBefore_pos",
                 "; TPC ratio findable (before); Entries", 1000, 0, 2);
    fHistV0SingleParticleNclsTPCRatioFindableBefore[1] =
        new TH1F("fHistV0SingleParticleNclsTPCRatioFindableBefore_neg",
                 ";  TPC ratio findable (before); Entries", 1000, 0, 2);
    fHistV0SingleParticleNclsTPCRatioFindableAfter[0] =
        new TH1F("fHistV0SingleParticleNclsTPCRatioFindableAfter_pos",
                 "; TPC ratio findable (after); Entries", 1000, 0, 2);
    fHistV0SingleParticleNclsTPCRatioFindableAfter[1] =
        new TH1F("fHistV0SingleParticleNclsTPCRatioFindableAfter_neg",
                 ";  TPC ratio findable (after); Entries", 1000, 0, 2);
    fHistV0SingleParticleNcrossedTPCBefore[0] =
        new TH1F("fHistV0SingleParticleNcrossedTPCBefore_pos",
                 "; # cls TPC crossed (before); Entries", 170, 0, 170);
    fHistV0SingleParticleNcrossedTPCBefore[1] =
        new TH1F("fHistV0SingleParticleNcrossedTPCBefore_neg",
                 "; # cls TPC crossed (before); Entries", 170, 0, 170);
    fHistV0SingleParticleNcrossedTPCAfter[0] =
        new TH1F("fHistV0SingleParticleNcrossedTPCAfter_pos",
                 "; # cls TPC crossed (after); Entries", 170, 0, 170);
    fHistV0SingleParticleNcrossedTPCAfter[1] =
        new TH1F("fHistV0SingleParticleNcrossedTPCAfter_neg",
                 "; # cls TPC crossed (after); Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCShared[0] =
        new TH1F("fHistV0SingleParticleNclsTPCShared_pos",
                 "; # cls TPC shared; Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCShared[1] =
        new TH1F("fHistV0SingleParticleNclsTPCShared_neg",
                 "; # cls TPC shared; Entries", 170, 0, 170);
    fHistV0SingleParticleNclsTPCSharedTiming[0] =
        new TH2F("fHistV0SingleParticleNclsTPCSharedTiming_pos",
                 "; Timing info; # cls TPC shared", 8, 0, 8, 170, 0, 170);
    fHistV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(1,
                                                                         "all");
    fHistV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistV0SingleParticleNclsTPCSharedTiming[1] =
        new TH2F("fHistV0SingleParticleNclsTPCSharedTiming_neg",
                 "; Timing info; # cls TPC shared", 8, 0, 8, 170, 0, 170);
    fHistV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(1,
                                                                         "all");
    fHistV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistV0SingleParticleNclsITSShared[0] =
        new TH1F("fHistV0SingleParticleNclsITSShared_pos",
                 "; # cls ITS shared; Entries", 6, 0, 6);
    fHistV0SingleParticleNclsITSShared[1] =
        new TH1F("fHistV0SingleParticleNclsITSShared_neg",
                 "; # cls ITS shared; Entries", 6, 0, 6);
    fHistV0SingleParticleNclsITSSharedTiming[0] =
        new TH2F("fHistV0SingleParticleNclsITSSharedTiming_pos",
                 "; Timing info; # cls ITS shared", 8, 0, 8, 6, 0, 6);
    fHistV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(1,
                                                                         "all");
    fHistV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistV0SingleParticleNclsITSSharedTiming[1] =
        new TH2F("fHistV0SingleParticleNclsITSSharedTiming_neg",
                 "; Timing info; # cls ITS shared", 8, 0, 8, 6, 0, 6);
    fHistV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(1,
                                                                         "all");
    fHistV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistV0SingleParticleDCAtoPVBefore[0] =
        new TH1F("fHistV0SingleParticleDCAtoPVBefore_pos",
                 "; DCA (pos, PA) (before); Entries", 1000, 0, 10);
    fHistV0SingleParticleDCAtoPVBefore[1] =
        new TH1F("fHistV0SingleParticleDCAtoPVBefore_neg",
                 "; DCA (pos, PA) (before); Entries", 1000, 0, 10);
    fHistV0SingleParticleDCAtoPVAfter[0] =
        new TH1F("fHistV0SingleParticleDCAtoPVAfter_pos",
                 "; DCA (pos, PA) (after); Entries", 1000, 0, 10);
    fHistV0SingleParticleDCAtoPVAfter[1] =
        new TH1F("fHistV0SingleParticleDCAtoPVAfter_neg",
                 "; DCA (pos, PA) (after); Entries", 1000, 0, 10);

    fHistogramsV0Pos->Add(fHistV0SingleParticlePt[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticlePt[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleEtaBefore[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleEtaBefore[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsTPCBefore[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsTPCBefore[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsTPCFindableBefore[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsTPCFindableBefore[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsTPCRatioFindableBefore[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsTPCRatioFindableBefore[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNcrossedTPCBefore[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNcrossedTPCBefore[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsTPCShared[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsTPCShared[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsTPCSharedTiming[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsTPCSharedTiming[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsITSShared[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsITSShared[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsITSSharedTiming[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsITSSharedTiming[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleDCAtoPVBefore[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleDCAtoPVBefore[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleEtaAfter[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleEtaAfter[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsTPCAfter[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsTPCAfter[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsTPCFindableAfter[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsTPCFindableAfter[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNclsTPCRatioFindableAfter[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNclsTPCRatioFindableAfter[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleNcrossedTPCAfter[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleNcrossedTPCAfter[1]);
    fHistogramsV0Pos->Add(fHistV0SingleParticleDCAtoPVAfter[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticleDCAtoPVAfter[1]);

    fHistogramsV0->Add(fHistogramsV0Before);
    fHistogramsV0->Add(fHistogramsV0After);

    fHistV0SingleParticlePID[0] = new TH2F(
        "fHistV0SingleParticlePID_pos", "; #it{p} [GeV/#it{c}]; n_{#sigma} TPC",
        1000, 0, 5, 1000, 0, 10);
    fHistV0SingleParticlePID[1] = new TH2F(
        "fHistV0SingleParticlePID_neg", "; #it{p} [GeV/#it{c}]; n_{#sigma} TPC",
        1000, 0, 5, 1000, 0, 10);
    fHistogramsV0Pos->Add(fHistV0SingleParticlePID[0]);
    fHistogramsV0Neg->Add(fHistV0SingleParticlePID[1]);

    if (fHistogramsAntiV0Before != nullptr) {
      delete fHistogramsAntiV0Before;
      fHistogramsAntiV0Before = nullptr;
    }
    if (fHistogramsAntiV0Before == nullptr) {
      fHistogramsAntiV0Before = new TList();
      fHistogramsAntiV0Before->SetOwner(kTRUE);
      fHistogramsAntiV0Before->SetName("Before");
    }

    if (fHistogramsAntiV0After != nullptr) {
      delete fHistogramsAntiV0After;
      fHistogramsAntiV0After = nullptr;
    }
    if (fHistogramsAntiV0After == nullptr) {
      fHistogramsAntiV0After = new TList();
      fHistogramsAntiV0After->SetOwner(kTRUE);
      fHistogramsAntiV0After->SetName("After");
    }

    fHistAntiV0DecayVertexXBefore =
        new TH1F("fHistAntiV0DecayVertexXBefore",
                 "; Decay vertex x [cm]; Entries", 1000, 0, 200);
    fHistAntiV0DecayVertexYBefore =
        new TH1F("fHistAntiV0DecayVertexYBefore",
                 "; Decay vertex y [cm]; Entries", 1000, 0, 200);
    fHistAntiV0DecayVertexZBefore =
        new TH1F("fHistAntiV0DecayVertexZBefore",
                 "; Decay vertex z [cm]; Entries", 1000, 0, 200);
    fHistAntiV0DecayVertexXAfter =
        new TH1F("fHistAntiV0DecayVertexXAfter",
                 "; Decay vertex x [cm]; Entries", 1000, 0, 200);
    fHistAntiV0DecayVertexYAfter =
        new TH1F("fHistAntiV0DecayVertexYAfter",
                 "; Decay vertex y [cm]; Entries", 1000, 0, 200);
    fHistAntiV0DecayVertexZAfter =
        new TH1F("fHistAntiV0DecayVertexZAfter",
                 "; Decay vertex z [cm]; Entries", 1000, 0, 200);
    fHistAntiV0TransverseRadiusBefore =
        new TH1F("fHistAntiV0TransverseRadiusBefore",
                 "; Transverse radius [cm]; Entries", 1000, 0, 200);
    fHistAntiV0TransverseRadiusAfter =
        new TH1F("fHistAntiV0TransverseRadiusAfter",
                 "; Transverse radius [cm]; Entries", 1000, 0, 200);
    fHistAntiV0CosPABefore = new TH1F("fHistAntiV0CosPABefore",
                                      "; cos(#alpha); Entries", 1000, 0.5, 1);
    fHistAntiV0CosPAAfter = new TH1F("fHistAntiV0CosPAAfter",
                                     "; cos(#alpha); Entries", 1000, 0.95, 1);
    fHistAntiV0DCADaughtersBefore =
        new TH1F("fHistAntiV0DCADaughtersBefore",
                 "; Daughter DCA at decay vertex [cm]; Entries", 1000, 0, 10);
    fHistAntiV0DCADaughtersAfter =
        new TH1F("fHistAntiV0DCADaughtersAfter",
                 "; Daughter DCA at decay vertex [cm]; Entries", 1000, 0, 10);
    fHistAntiV0DCA =
        new TH1F("fHistAntiV0DCA", "; DCA to PV [cm]; Entries", 1000, 0, 10);
    fHistAntiV0DecayLength = new TH1F(
        "fHistAntiV0DecayLength", "; Decay length [cm]; Entries", 1000, 0, 200);
    fHistAntiV0ArmenterosBefore = new TH2F(
        "fHistAntiV0ArmenterosBefore",
        " ; #alpha; #it{q}_{T} p#pi [GeV/#it{c}]", 500, -1, 1, 500, 0, 1);
    fHistAntiV0ArmenterosAfter = new TH2F(
        "fHistAntiV0ArmenterosAfter", " ; #alpha; #it{q}_{T} p#pi [GeV/#it{c}]",
        500, -1, 1, 500, 0, 1);

    fHistogramsAntiV0Before->Add(fHistAntiV0DecayVertexXBefore);
    fHistogramsAntiV0Before->Add(fHistAntiV0DecayVertexYBefore);
    fHistogramsAntiV0Before->Add(fHistAntiV0DecayVertexZBefore);
    fHistogramsAntiV0Before->Add(fHistAntiV0TransverseRadiusBefore);
    fHistogramsAntiV0Before->Add(fHistAntiV0CosPABefore);
    fHistogramsAntiV0Before->Add(fHistAntiV0DCADaughtersBefore);
    fHistogramsAntiV0After->Add(fHistAntiV0DecayVertexXAfter);
    fHistogramsAntiV0After->Add(fHistAntiV0DecayVertexYAfter);
    fHistogramsAntiV0After->Add(fHistAntiV0DecayVertexZAfter);
    fHistogramsAntiV0After->Add(fHistAntiV0TransverseRadiusAfter);
    fHistogramsAntiV0After->Add(fHistAntiV0CosPAAfter);
    fHistogramsAntiV0After->Add(fHistAntiV0DCADaughtersAfter);
    fHistogramsAntiV0->Add(fHistAntiV0DCA);
    fHistogramsAntiV0->Add(fHistAntiV0DecayLength);
    fHistogramsAntiV0Before->Add(fHistAntiV0ArmenterosBefore);
    fHistogramsAntiV0After->Add(fHistAntiV0ArmenterosAfter);

    fHistAntiV0SingleParticlePt[0] = new TH1F(
        "fHistAntiV0SingleParticlePt_pos", ";#it{p}_{T}; Entries", 1000, 0, 10);
    fHistAntiV0SingleParticlePt[1] = new TH1F(
        "fHistAntiV0SingleParticlePt_neg", ";#it{p}_{T}; Entries", 1000, 0, 10);
    fHistAntiV0SingleParticleEtaBefore[0] =
        new TH1F("fHistAntiV0SingleParticleEtaBefore_pos",
                 "; #eta (before); Entries", 1000, -2, 2);
    fHistAntiV0SingleParticleEtaBefore[1] =
        new TH1F("fHistAntiV0SingleParticleEtaBefore_neg",
                 "; #eta (before); Entries", 1000, -2, 2);
    fHistAntiV0SingleParticleEtaAfter[0] =
        new TH1F("fHistAntiV0SingleParticleEtaAfter_pos",
                 "; #eta (after); Entries", 1000, -2, 2);
    fHistAntiV0SingleParticleEtaAfter[1] =
        new TH1F("fHistAntiV0SingleParticleEtaAfter_neg",
                 "; #eta (after); Entries", 1000, -2, 2);
    fHistAntiV0SingleParticleNclsTPCBefore[0] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCBefore_pos",
                 "; # cls TPC (before); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCBefore[1] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCBefore_neg",
                 "; # cls TPC (before); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCAfter[0] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCAfter_pos",
                 "; # cls TPC (after); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCAfter[1] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCAfter_neg",
                 "; # cls TPC (after); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCFindableBefore[0] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCFindableBefore_pos",
                 "; # cls TPC findable (before); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCFindableBefore[1] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCFindableBefore_neg",
                 "; # cls TPC findable (before); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCFindableAfter[0] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCFindableAfter_pos",
                 "; # cls TPC findable (after); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCFindableAfter[1] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCFindableAfter_neg",
                 "; # cls TPC findable (after); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCRatioFindableBefore[0] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCRatioFindableBefore_pos",
                 "; TPC ratio findable (before); Entries", 1000, 0, 2);
    fHistAntiV0SingleParticleNclsTPCRatioFindableBefore[1] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCRatioFindableBefore_neg",
                 ";  TPC ratio findable (before); Entries", 1000, 0, 2);
    fHistAntiV0SingleParticleNclsTPCRatioFindableAfter[0] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCRatioFindableAfter_pos",
                 "; TPC ratio findable (after); Entries", 1000, 0, 2);
    fHistAntiV0SingleParticleNclsTPCRatioFindableAfter[1] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCRatioFindableAfter_neg",
                 ";  TPC ratio findable (after); Entries", 1000, 0, 2);
    fHistAntiV0SingleParticleNcrossedTPCBefore[0] =
        new TH1F("fHistAntiV0SingleParticleNcrossedTPCBefore_pos",
                 "; # cls TPC crossed (before); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNcrossedTPCBefore[1] =
        new TH1F("fHistAntiV0SingleParticleNcrossedTPCBefore_neg",
                 "; # cls TPC crossed (before); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNcrossedTPCAfter[0] =
        new TH1F("fHistAntiV0SingleParticleNcrossedTPCAfter_pos",
                 "; # cls TPC crossed (after); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNcrossedTPCAfter[1] =
        new TH1F("fHistAntiV0SingleParticleNcrossedTPCAfter_neg",
                 "; # cls TPC crossed (after); Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCShared[0] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCShared_pos",
                 "; # cls TPC shared; Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCShared[1] =
        new TH1F("fHistAntiV0SingleParticleNclsTPCShared_neg",
                 "; # cls TPC shared; Entries", 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0] =
        new TH2F("fHistAntiV0SingleParticleNclsTPCSharedTiming_pos",
                 "; Timing info; # cls TPC shared", 8, 0, 8, 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        1, "all");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[0]->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1] =
        new TH2F("fHistAntiV0SingleParticleNclsTPCSharedTiming_neg",
                 "; Timing info; # cls TPC shared", 8, 0, 8, 170, 0, 170);
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        1, "all");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistAntiV0SingleParticleNclsTPCSharedTiming[1]->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistAntiV0SingleParticleNclsITSShared[0] =
        new TH1F("fHistAntiV0SingleParticleNclsITSShared_pos",
                 "; # cls ITS shared; Entries", 6, 0, 6);
    fHistAntiV0SingleParticleNclsITSShared[1] =
        new TH1F("fHistAntiV0SingleParticleNclsITSShared_neg",
                 "; # cls ITS shared; Entries", 6, 0, 6);
    fHistAntiV0SingleParticleNclsITSSharedTiming[0] =
        new TH2F("fHistAntiV0SingleParticleNclsITSSharedTiming_pos",
                 "; Timing info; # cls ITS shared", 8, 0, 8, 6, 0, 6);
    fHistAntiV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        1, "all");
    fHistAntiV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistAntiV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistAntiV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistAntiV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistAntiV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistAntiV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistAntiV0SingleParticleNclsITSSharedTiming[0]->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistAntiV0SingleParticleNclsITSSharedTiming[1] =
        new TH2F("fHistAntiV0SingleParticleNclsITSSharedTiming_neg",
                 "; Timing info; # cls ITS shared", 8, 0, 8, 6, 0, 6);
    fHistAntiV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        1, "all");
    fHistAntiV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistAntiV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistAntiV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistAntiV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistAntiV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistAntiV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistAntiV0SingleParticleNclsITSSharedTiming[1]->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistAntiV0SingleParticleDCAtoPVBefore[0] =
        new TH1F("fHistAntiV0SingleParticleDCAtoPVBefore_pos",
                 "; DCA (pos, PA) (before); Entries", 1000, 0, 10);
    fHistAntiV0SingleParticleDCAtoPVBefore[1] =
        new TH1F("fHistAntiV0SingleParticleDCAtoPVBefore_neg",
                 "; DCA (pos, PA) (before); Entries", 1000, 0, 10);
    fHistAntiV0SingleParticleDCAtoPVAfter[0] =
        new TH1F("fHistAntiV0SingleParticleDCAtoPVAfter_pos",
                 "; DCA (pos, PA) (after); Entries", 1000, 0, 10);
    fHistAntiV0SingleParticleDCAtoPVAfter[1] =
        new TH1F("fHistAntiV0SingleParticleDCAtoPVAfter_neg",
                 "; DCA (pos, PA) (after); Entries", 1000, 0, 10);

    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticlePt[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticlePt[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleEtaBefore[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleEtaBefore[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNclsTPCBefore[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNclsTPCBefore[1]);
    fHistogramsAntiV0Pos->Add(
        fHistAntiV0SingleParticleNclsTPCFindableBefore[0]);
    fHistogramsAntiV0Neg->Add(
        fHistAntiV0SingleParticleNclsTPCFindableBefore[1]);
    fHistogramsAntiV0Pos->Add(
        fHistAntiV0SingleParticleNclsTPCRatioFindableBefore[0]);
    fHistogramsAntiV0Neg->Add(
        fHistAntiV0SingleParticleNclsTPCRatioFindableBefore[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNcrossedTPCBefore[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNcrossedTPCBefore[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNclsTPCShared[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNclsTPCShared[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNclsTPCSharedTiming[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNclsTPCSharedTiming[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNclsITSShared[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNclsITSShared[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNclsITSSharedTiming[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNclsITSSharedTiming[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleDCAtoPVBefore[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleDCAtoPVBefore[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleEtaAfter[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleEtaAfter[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNclsTPCAfter[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNclsTPCAfter[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNclsTPCFindableAfter[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNclsTPCFindableAfter[1]);
    fHistogramsAntiV0Pos->Add(
        fHistAntiV0SingleParticleNclsTPCRatioFindableAfter[0]);
    fHistogramsAntiV0Neg->Add(
        fHistAntiV0SingleParticleNclsTPCRatioFindableAfter[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleNcrossedTPCAfter[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleNcrossedTPCAfter[1]);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticleDCAtoPVAfter[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticleDCAtoPVAfter[1]);

    fHistAntiV0SingleParticlePID[0] = new TH2F(
        "fHistAntiV0SingleParticlePID_pos",
        "; #it{p} [GeV/#it{c}]; n_{#sigma} TPC", 1000, 0, 5, 1000, 0, 10);
    fHistAntiV0SingleParticlePID[1] = new TH2F(
        "fHistAntiV0SingleParticlePID_neg",
        "; #it{p} [GeV/#it{c}]; n_{#sigma} TPC", 1000, 0, 5, 1000, 0, 10);
    fHistogramsAntiV0Pos->Add(fHistAntiV0SingleParticlePID[0]);
    fHistogramsAntiV0Neg->Add(fHistAntiV0SingleParticlePID[1]);

    fHistogramsAntiV0->Add(fHistogramsAntiV0Before);
    fHistogramsAntiV0->Add(fHistogramsAntiV0After);
  }

  fHistogramsV0->Add(fHistogramsV0Pos);
  fHistogramsV0->Add(fHistogramsV0Neg);
  fHistogramsAntiV0->Add(fHistogramsAntiV0Pos);
  fHistogramsAntiV0->Add(fHistogramsAntiV0Neg);

  fHistograms->Add(fHistogramsV0);
  fHistograms->Add(fHistogramsAntiV0);
}
