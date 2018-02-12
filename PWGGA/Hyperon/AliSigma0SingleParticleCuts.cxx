#include "AliSigma0SingleParticleCuts.h"

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"

#include <iostream>

ClassImp(AliSigma0SingleParticleCuts)

    //____________________________________________________________________________________________________
    AliSigma0SingleParticleCuts::AliSigma0SingleParticleCuts()
    : TObject(),
      fHistograms(nullptr),
      fHistogramsSingleParticle(nullptr),
      fHistogramsSingleParticleMC(nullptr),
      fHistogramsSingleParticleBefore(nullptr),
      fHistogramsSingleParticleAfter(nullptr),
      fHistogramsSingleAntiParticle(nullptr),
      fHistogramsSingleAntiParticleMC(nullptr),
      fHistogramsSingleAntiParticleBefore(nullptr),
      fHistogramsSingleAntiParticleAfter(nullptr),
      fIsMC(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fSingleParticleVector(),
      fSingleAntiParticleVector(),
      fGlobalTrackReference(),
      fTrackReferenceSize(1000),
      fParticleName(),
      fParticle(AliPID::kProton),
      fParticleVeto1(AliPID::kKaon),
      fParticleVeto2(AliPID::kPion),
      fParticleVeto3(AliPID::kElectron),
      fPDGCode(2212),
      fIsExtendedQA(false),
      fFilterbit(128),
      fTPCclusterMin(0.f),
      fTPCnCrossedRowsMin(0.f),
      fTPCfindableMin(0.f),
      fUseSharedMap(true),
      fDCArMax(999.f),
      fDCAzMax(999.f),
      fPtMin(0.f),
      fPtMax(999.f),
      fEtaMax(999.f),
      fPIDnSigmaStrict(0.f),
      fPIDMomentumSwitch(999.f),
      fPIDHypothesisRejection(true),
      fPIDAdditionalRejection(false),
      fPIDTOFif(false),
      fParticleAdditionalRejection(),
      fPIDAdditionalRejectionUpperBound(),
      fPIDAdditionalRejectionLowerBound(),
      fPIDAdditionalRejectionUpperpTBound(),
      fPIDAdditionalRejectionLowerpTBound(),
      fPIDResponse(nullptr),
      fHistCuts(nullptr),
      fHistSingleParticleCuts(nullptr),
      fHistSingleParticlePhi(nullptr),
      fHistSingleParticlePt(nullptr),
      fHistSingleParticleEta(nullptr),
      fHistSingleParticleEtaPhi(nullptr),
      fHistNSingleParticle(nullptr),
      fHistSingleParticleDCAxy(nullptr),
      fHistSingleParticleEtaBefore(nullptr),
      fHistSingleParticleEtaAfter(nullptr),
      fHistSingleParticlePtBefore(nullptr),
      fHistSingleParticlePtAfter(nullptr),
      fHistSingleParticleNclsTPCBefore(nullptr),
      fHistSingleParticleNclsTPCAfter(nullptr),
      fHistSingleParticleNclsTPCShared(nullptr),
      fHistSingleParticleNclsTPCSharedTiming(nullptr),
      fHistSingleParticleNclsITSShared(nullptr),
      fHistSingleParticleNclsITSSharedTiming(nullptr),
      fHistSingleParticleNcrossedTPCBefore(nullptr),
      fHistSingleParticleNcrossedTPCAfter(nullptr),
      fHistSingleParticleFindableTPCBefore(nullptr),
      fHistSingleParticleFindableTPCAfter(nullptr),
      fHistSingleParticleDCAxyBefore(nullptr),
      fHistSingleParticleDCAxyAfterDCAz(nullptr),
      fHistSingleParticleDCAxyAfter(nullptr),
      fHistSingleParticleDCAzBefore(nullptr),
      fHistSingleParticleDCAzAfter(nullptr),
      fHistSingleParticleNsigmaBefore(nullptr),
      fHistSingleParticleNsigmaAfter(nullptr),
      fHistSingleParticleNsigmaTPCBefore(nullptr),
      fHistSingleParticleNsigmaTPCAfter(nullptr),
      fHistSingleParticleNsigmaTOFBefore(nullptr),
      fHistSingleParticleNsigmaTOFAfter(nullptr),
      fHistSingleParticleTPCsignalBefore(nullptr),
      fHistSingleParticleTOFsignalBefore(nullptr),
      fHistSingleParticleTPCsignalAfter(nullptr),
      fHistSingleParticleTOFsignalAfter(nullptr),
      fHistSingleParticleNsigmaTPCRejBefore(),
      fHistSingleParticleNsigmaTPCRejAfter(),
      fHistMCRecSingleParticlePtTruth(nullptr),
      fHistMCRecSingleParticleMomentum(nullptr),
      fHistSingleAntiParticleCuts(nullptr),
      fHistSingleAntiParticlePhi(nullptr),
      fHistSingleAntiParticlePt(nullptr),
      fHistSingleAntiParticleEta(nullptr),
      fHistSingleAntiParticleEtaPhi(nullptr),
      fHistNSingleAntiParticle(nullptr),
      fHistSingleAntiParticleDCAxy(nullptr),
      fHistSingleAntiParticleEtaBefore(nullptr),
      fHistSingleAntiParticleEtaAfter(nullptr),
      fHistSingleAntiParticlePtBefore(nullptr),
      fHistSingleAntiParticlePtAfter(nullptr),
      fHistSingleAntiParticleNclsTPCBefore(nullptr),
      fHistSingleAntiParticleNclsTPCAfter(nullptr),
      fHistSingleAntiParticleNclsTPCShared(nullptr),
      fHistSingleAntiParticleNclsTPCSharedTiming(nullptr),
      fHistSingleAntiParticleNclsITSShared(nullptr),
      fHistSingleAntiParticleNclsITSSharedTiming(nullptr),
      fHistSingleAntiParticleNcrossedTPCBefore(nullptr),
      fHistSingleAntiParticleNcrossedTPCAfter(nullptr),
      fHistSingleAntiParticleFindableTPCBefore(nullptr),
      fHistSingleAntiParticleFindableTPCAfter(nullptr),
      fHistSingleAntiParticleDCAxyBefore(nullptr),
      fHistSingleAntiParticleDCAxyAfterDCAz(nullptr),
      fHistSingleAntiParticleDCAxyAfter(nullptr),
      fHistSingleAntiParticleDCAzBefore(nullptr),
      fHistSingleAntiParticleDCAzAfter(nullptr),
      fHistSingleAntiParticleNsigmaBefore(nullptr),
      fHistSingleAntiParticleNsigmaAfter(nullptr),
      fHistSingleAntiParticleNsigmaTPCBefore(nullptr),
      fHistSingleAntiParticleNsigmaTPCAfter(nullptr),
      fHistSingleAntiParticleNsigmaTOFBefore(nullptr),
      fHistSingleAntiParticleNsigmaTOFAfter(nullptr),
      fHistSingleAntiParticleTPCsignalBefore(nullptr),
      fHistSingleAntiParticleTOFsignalBefore(nullptr),
      fHistSingleAntiParticleTPCsignalAfter(nullptr),
      fHistSingleAntiParticleTOFsignalAfter(nullptr),
      fHistSingleAntiParticleNsigmaTPCRejBefore(),
      fHistSingleAntiParticleNsigmaTPCRejAfter(),
      fHistMCRecSingleAntiParticlePtTruth(nullptr),
      fHistMCRecSingleAntiParticleMomentum(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0SingleParticleCuts::AliSigma0SingleParticleCuts(
    const AliSigma0SingleParticleCuts &ref)
    : TObject(ref),
      fHistograms(nullptr),
      fHistogramsSingleParticle(nullptr),
      fHistogramsSingleParticleMC(nullptr),
      fHistogramsSingleParticleBefore(nullptr),
      fHistogramsSingleParticleAfter(nullptr),
      fHistogramsSingleAntiParticle(nullptr),
      fHistogramsSingleAntiParticleMC(nullptr),
      fHistogramsSingleAntiParticleBefore(nullptr),
      fHistogramsSingleAntiParticleAfter(nullptr),
      fIsMC(false),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fSingleParticleVector(),
      fSingleAntiParticleVector(),
      fGlobalTrackReference(),
      fTrackReferenceSize(1000),
      fParticleName(),
      fParticle(AliPID::kProton),
      fParticleVeto1(AliPID::kKaon),
      fParticleVeto2(AliPID::kPion),
      fParticleVeto3(AliPID::kElectron),
      fPDGCode(2212),
      fIsExtendedQA(false),
      fFilterbit(128),
      fTPCclusterMin(0),
      fTPCnCrossedRowsMin(0),
      fTPCfindableMin(0.f),
      fUseSharedMap(true),
      fDCArMax(999.f),
      fDCAzMax(999.f),
      fPtMin(0.f),
      fPtMax(999.f),
      fEtaMax(999.f),
      fPIDnSigmaStrict(0.f),
      fPIDMomentumSwitch(999.f),
      fPIDHypothesisRejection(true),
      fPIDAdditionalRejection(false),
      fPIDTOFif(false),
      fParticleAdditionalRejection(),
      fPIDAdditionalRejectionUpperBound(),
      fPIDAdditionalRejectionLowerBound(),
      fPIDAdditionalRejectionUpperpTBound(),
      fPIDAdditionalRejectionLowerpTBound(),
      fPIDResponse(nullptr),
      fHistCuts(nullptr),
      fHistSingleParticleCuts(nullptr),
      fHistSingleParticlePhi(nullptr),
      fHistSingleParticlePt(nullptr),
      fHistSingleParticleEta(nullptr),
      fHistSingleParticleEtaPhi(nullptr),
      fHistNSingleParticle(nullptr),
      fHistSingleParticleDCAxy(nullptr),
      fHistSingleParticleEtaBefore(nullptr),
      fHistSingleParticleEtaAfter(nullptr),
      fHistSingleParticlePtBefore(nullptr),
      fHistSingleParticlePtAfter(nullptr),
      fHistSingleParticleNclsTPCBefore(nullptr),
      fHistSingleParticleNclsTPCAfter(nullptr),
      fHistSingleParticleNclsTPCShared(nullptr),
      fHistSingleParticleNclsTPCSharedTiming(nullptr),
      fHistSingleParticleNclsITSShared(nullptr),
      fHistSingleParticleNclsITSSharedTiming(nullptr),
      fHistSingleParticleNcrossedTPCBefore(nullptr),
      fHistSingleParticleNcrossedTPCAfter(nullptr),
      fHistSingleParticleFindableTPCBefore(nullptr),
      fHistSingleParticleFindableTPCAfter(nullptr),
      fHistSingleParticleDCAxyBefore(nullptr),
      fHistSingleParticleDCAxyAfterDCAz(nullptr),
      fHistSingleParticleDCAxyAfter(nullptr),
      fHistSingleParticleDCAzBefore(nullptr),
      fHistSingleParticleDCAzAfter(nullptr),
      fHistSingleParticleNsigmaBefore(nullptr),
      fHistSingleParticleNsigmaAfter(nullptr),
      fHistSingleParticleNsigmaTPCBefore(nullptr),
      fHistSingleParticleNsigmaTPCAfter(nullptr),
      fHistSingleParticleNsigmaTOFBefore(nullptr),
      fHistSingleParticleNsigmaTOFAfter(nullptr),
      fHistSingleParticleTPCsignalBefore(nullptr),
      fHistSingleParticleTOFsignalBefore(nullptr),
      fHistSingleParticleTPCsignalAfter(nullptr),
      fHistSingleParticleTOFsignalAfter(nullptr),
      fHistSingleParticleNsigmaTPCRejBefore(),
      fHistSingleParticleNsigmaTPCRejAfter(),
      fHistMCRecSingleParticlePtTruth(nullptr),
      fHistMCRecSingleParticleMomentum(nullptr),
      fHistSingleAntiParticleCuts(nullptr),
      fHistSingleAntiParticlePhi(nullptr),
      fHistSingleAntiParticlePt(nullptr),
      fHistSingleAntiParticleEta(nullptr),
      fHistSingleAntiParticleEtaPhi(nullptr),
      fHistNSingleAntiParticle(nullptr),
      fHistSingleAntiParticleDCAxy(nullptr),
      fHistSingleAntiParticleEtaBefore(nullptr),
      fHistSingleAntiParticleEtaAfter(nullptr),
      fHistSingleAntiParticlePtBefore(nullptr),
      fHistSingleAntiParticlePtAfter(nullptr),
      fHistSingleAntiParticleNclsTPCBefore(nullptr),
      fHistSingleAntiParticleNclsTPCAfter(nullptr),
      fHistSingleAntiParticleNclsTPCShared(nullptr),
      fHistSingleAntiParticleNclsTPCSharedTiming(nullptr),
      fHistSingleAntiParticleNclsITSShared(nullptr),
      fHistSingleAntiParticleNclsITSSharedTiming(nullptr),
      fHistSingleAntiParticleNcrossedTPCBefore(nullptr),
      fHistSingleAntiParticleNcrossedTPCAfter(nullptr),
      fHistSingleAntiParticleFindableTPCBefore(nullptr),
      fHistSingleAntiParticleFindableTPCAfter(nullptr),
      fHistSingleAntiParticleDCAxyBefore(nullptr),
      fHistSingleAntiParticleDCAxyAfterDCAz(nullptr),
      fHistSingleAntiParticleDCAxyAfter(nullptr),
      fHistSingleAntiParticleDCAzBefore(nullptr),
      fHistSingleAntiParticleDCAzAfter(nullptr),
      fHistSingleAntiParticleNsigmaBefore(nullptr),
      fHistSingleAntiParticleNsigmaAfter(nullptr),
      fHistSingleAntiParticleNsigmaTPCBefore(nullptr),
      fHistSingleAntiParticleNsigmaTPCAfter(nullptr),
      fHistSingleAntiParticleNsigmaTOFBefore(nullptr),
      fHistSingleAntiParticleNsigmaTOFAfter(nullptr),
      fHistSingleAntiParticleTPCsignalBefore(nullptr),
      fHistSingleAntiParticleTOFsignalBefore(nullptr),
      fHistSingleAntiParticleTPCsignalAfter(nullptr),
      fHistSingleAntiParticleTOFsignalAfter(nullptr),
      fHistSingleAntiParticleNsigmaTPCRejBefore(),
      fHistSingleAntiParticleNsigmaTPCRejAfter(),
      fHistMCRecSingleAntiParticlePtTruth(nullptr),
      fHistMCRecSingleAntiParticleMomentum(nullptr) {}

//____________________________________________________________________________________________________
AliSigma0SingleParticleCuts &AliSigma0SingleParticleCuts::operator=(
    const AliSigma0SingleParticleCuts &ref) {
  // Assignment operator
  if (this == &ref) return *this;

  fSingleParticleVector = ref.fSingleParticleVector;
  fTPCclusterMin = ref.fTPCclusterMin;
  fTPCnCrossedRowsMin = ref.fTPCnCrossedRowsMin;
  fTPCfindableMin = ref.fTPCfindableMin;
  fUseSharedMap = ref.fUseSharedMap;
  fDCArMax = ref.fDCArMax;
  fDCAzMax = ref.fDCAzMax;
  fPtMin = ref.fPtMin;
  fPtMax = ref.fPtMax;
  fEtaMax = ref.fEtaMax;

  return (*this);
}

//____________________________________________________________________________________________________
AliSigma0SingleParticleCuts::~AliSigma0SingleParticleCuts() {
  if (fPIDResponse) delete fPIDResponse;
}

//____________________________________________________________________________________________________
AliSigma0SingleParticleCuts *AliSigma0SingleParticleCuts::DefaultCuts() {
  AliSigma0SingleParticleCuts *spCuts = new AliSigma0SingleParticleCuts();
  spCuts->SetParticle(AliPID::kProton);
  spCuts->SetParticleName("Proton");
  spCuts->SetParticleRejection(AliPID::kKaon, AliPID::kPion, AliPID::kElectron);
  spCuts->SetPDGCode(2212);
  spCuts->SetFilterBit(128);
  spCuts->SetTPCclusterMin(80);
  spCuts->SetTPCnCrossedRowsMin(70);
  spCuts->SetTPCFindableMin(0.83);
  spCuts->SetUseSharedMap(true);
  spCuts->SetDCArMax(0.1);
  spCuts->SetDCAzMax(0.2);
  spCuts->SetPtMin(0.5);
  spCuts->SetPtMax(4.05);
  spCuts->SetEtaMax(0.8);
  spCuts->SetPIDnSigmaStrict(3);
  spCuts->SetPIDHypothesisRejection(true);
  spCuts->SetPIDMomentumSwitch(0.75);
  spCuts->SetPIDTOFif(false);
  return spCuts;
}

//____________________________________________________________________________________________________
AliSigma0SingleParticleCuts *AliSigma0SingleParticleCuts::ElectronCuts() {
  AliSigma0SingleParticleCuts *spCuts = new AliSigma0SingleParticleCuts();
  spCuts->SetParticle(AliPID::kElectron);
  spCuts->SetParticleName("Electron");
  spCuts->SetParticleRejection(AliPID::kKaon, AliPID::kPion, AliPID::kProton);
  spCuts->SetAdditionalPIDRejection(AliPID::kPion, -1.f, 1.f, 0.f, 0.5f);
  spCuts->SetAdditionalPIDRejection(AliPID::kPion, -1E30f, 2.f, 0.5, 1E30f);
  spCuts->SetAdditionalPIDRejection(AliPID::kProton, -1.f, 1.f);
  spCuts->SetAdditionalPIDRejection(AliPID::kKaon, -1.f, 1.f);
  spCuts->SetPDGCode(11);
  spCuts->SetFilterBit(96);
  spCuts->SetTPCclusterMin(0);
  spCuts->SetTPCnCrossedRowsMin(0);
  spCuts->SetTPCFindableMin(0.35);
  spCuts->SetUseSharedMap(false);
  spCuts->SetDCArMax(1);
  spCuts->SetDCAzMax(1);
  spCuts->SetPtMin(0.025);
  spCuts->SetPtMax(1E30);
  spCuts->SetEtaMax(0.9);
  spCuts->SetPIDnSigmaStrict(3);
  spCuts->SetPIDHypothesisRejection(true);
  spCuts->SetPIDMomentumSwitch(0.);
  spCuts->SetPIDTOFif(true);
  return spCuts;
}

//____________________________________________________________________________________________________
void AliSigma0SingleParticleCuts::SelectSingleParticles(
    AliVEvent *inputEvent, AliMCEvent *mcEvent,
    std::vector<AliSigma0ParticleBase> &ParticleContainer,
    std::vector<AliSigma0ParticleBase> &AntiParticleContainer) {
  // Get the PID manager from the input event handler
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler =
      static_cast<AliInputEventHandler *>(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  fMCEvent = mcEvent;

  fSingleParticleVector.clear();
  fSingleAntiParticleVector.clear();

  fInputEvent = inputEvent;
  //  fMCEvent = static_cast<AliMCEvent*>(mcEvent);

  if (fInputEvent->IsA() == AliESDEvent::Class()) {
    ProcessESDs();
  }
  if (fInputEvent->IsA() == AliAODEvent::Class()) {
    ProcessAODs();
  }

  fHistNSingleParticle->Fill(fSingleParticleVector.size());
  fHistNSingleAntiParticle->Fill(fSingleAntiParticleVector.size());

  ParticleContainer.swap(fSingleParticleVector);
  AntiParticleContainer.swap(fSingleAntiParticleVector);
}

//____________________________________________________________________________________________________
void AliSigma0SingleParticleCuts::ProcessESDs() {
  AliESDEvent *esdEvent = static_cast<AliESDEvent *>(fInputEvent);
  for (int iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); ++iTrack) {
    AliESDtrack *track = static_cast<AliESDtrack *>(esdEvent->GetTrack(iTrack));
    if (!track) continue;

    const int charge = track->Charge();
    if (charge > 0)
      fHistSingleParticleCuts->Fill(0);
    else
      fHistSingleAntiParticleCuts->Fill(0);

    // PID
    if (!SingleParticlePID(track)) continue;
    if (charge > 0)
      fHistSingleParticleCuts->Fill(3);
    else
      fHistSingleAntiParticleCuts->Fill(3);

    // Quality cuts
    if (!SingleParticleQualityCuts(track)) continue;

    // Set DCA
    float DCAr = -999.f;
    float DCAz = -999.f;
    if (!SingleParticleDCASelection(track, DCAr, DCAz)) continue;

    // PARTICLE

    /// @todo set pid via cut
    AliSigma0ParticleBase *fCandidate = new AliSigma0ParticleBase(
        *(static_cast<AliVTrack *>(track)), track->Charge() * fPDGCode,
        fInputEvent->GetMagneticField(), fFilterbit);

    /// @todo set mass
    if (fPDGCode == 2212)
      fCandidate->SetMass(0.938272);
    else if (fPDGCode == 11)
      fCandidate->SetMass(0.00051099894);
    fCandidate->SetDCA(DCAr, DCAz);

    if (fIsMC) {
      AliMCParticle *mcParticle =
          static_cast<AliMCParticle *>(fMCEvent->GetTrack(track->GetLabel()));
      if (!mcParticle) continue;
      fCandidate->ProcessMCInfo(mcParticle, fMCEvent);
      ProcessMC(fCandidate, mcParticle);
    }

    if (charge > 0) {
      fHistSingleParticlePt->Fill(track->Pt());
      fHistSingleParticlePhi->Fill(track->Phi());
      fHistSingleParticleEta->Fill(track->Eta());
      fHistSingleParticleEtaPhi->Fill(track->Eta(), track->Phi());
    } else {
      fHistSingleAntiParticlePt->Fill(track->Pt());
      fHistSingleAntiParticlePhi->Fill(track->Phi());
      fHistSingleAntiParticleEta->Fill(track->Eta());
      fHistSingleAntiParticleEtaPhi->Fill(track->Eta(), track->Phi());
    }

    if (fCandidate) {
      if (track->Charge() > 0)
        fSingleParticleVector.push_back(*fCandidate);
      else
        fSingleAntiParticleVector.push_back(*fCandidate);
      delete fCandidate;
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0SingleParticleCuts::ProcessAODs() {
  // Store the reference to the global tracks
  StoreGlobalTrackReference();

  AliAODEvent *aodEvent = static_cast<AliAODEvent *>(fInputEvent);
  for (int iTrack = 0; iTrack < aodEvent->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(aodEvent->GetTrack(iTrack));
    if (!track) continue;

    const int charge = track->Charge();
    if (charge > 0)
      fHistSingleParticleCuts->Fill(0);
    else
      fHistSingleAntiParticleCuts->Fill(0);

    if (!track->TestFilterBit(fFilterbit)) continue;
    if (charge > 0)
      fHistSingleParticleCuts->Fill(1);
    else
      fHistSingleAntiParticleCuts->Fill(1);

    if (fFilterbit == 128 && !(fGlobalTrackReference[-track->GetID() - 1]))
      continue;
    if (charge > 0)
      fHistSingleParticleCuts->Fill(2);
    else
      fHistSingleAntiParticleCuts->Fill(2);

    // PID
    if (!SingleParticlePID(track)) continue;
    if (charge > 0)
      fHistSingleParticleCuts->Fill(3);
    else
      fHistSingleAntiParticleCuts->Fill(3);

    // Quality cuts
    if (!SingleParticleQualityCuts(track)) continue;

    // Set DCA
    float DCAr = -999.f;
    float DCAz = -999.f;
    if (!SingleParticleDCASelection(track, DCAr, DCAz)) continue;

    // PARTICLE

    /// @todo set pid via cut
    AliSigma0ParticleBase *fCandidate = new AliSigma0ParticleBase(
        *(static_cast<AliVTrack *>(track)), track->Charge() * fPDGCode,
        fInputEvent->GetMagneticField(), fFilterbit);

    /// @todo set mass
    if (fPDGCode == 2212)
      fCandidate->SetMass(0.938272);
    else if (fPDGCode == 11)
      fCandidate->SetMass(0.00051099894);
    fCandidate->SetDCA(DCAr, DCAz);

    if (fIsMC) {
      AliMCParticle *mcParticle =
          static_cast<AliMCParticle *>(fMCEvent->GetTrack(track->GetLabel()));
      if (!mcParticle) continue;
      fCandidate->ProcessMCInfo(mcParticle, fMCEvent);
      ProcessMC(fCandidate, mcParticle);
    }

    if (charge > 0) {
      fHistSingleParticlePt->Fill(track->Pt());
      fHistSingleParticlePhi->Fill(track->Phi());
      fHistSingleParticleEta->Fill(track->Eta());
      fHistSingleParticleEtaPhi->Fill(track->Eta(), track->Phi());
    } else {
      fHistSingleAntiParticlePt->Fill(track->Pt());
      fHistSingleAntiParticlePhi->Fill(track->Phi());
      fHistSingleAntiParticleEta->Fill(track->Eta());
      fHistSingleAntiParticleEtaPhi->Fill(track->Eta(), track->Phi());
    }

    if (fCandidate) {
      if (track->Charge() > 0)
        fSingleParticleVector.push_back(*fCandidate);
      else
        fSingleAntiParticleVector.push_back(*fCandidate);
      delete fCandidate;
    }
  }
}

//____________________________________________________________________________________________________
void AliSigma0SingleParticleCuts::StoreGlobalTrackReference() {
  // This method was inherited form H. Beck & O. Arnold analysis
  // Stores the pointer to the global track
  // Modified to work with vectors

  fGlobalTrackReference.clear();
  fGlobalTrackReference.resize(fTrackReferenceSize);
  AliAODEvent *aodEvent = static_cast<AliAODEvent *>(fInputEvent);
  for (int iTrack = 0; iTrack < aodEvent->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack *>(aodEvent->GetTrack(iTrack));
    if (!track) continue;

    // Check that the id is positive
    if (track->GetID() < 0) continue;

    // Check id is not too big for buffer
    if (track->GetID() >= static_cast<int>(fGlobalTrackReference.size()))
      fGlobalTrackReference.resize(track->GetID() + 1);

    // Warn if we overwrite a track
    auto *trackRef = fGlobalTrackReference[track->GetID()];
    if (trackRef) {
      // Seems like there are FilterMap 0 tracks
      // that have zero TPCNcls, don't store these!
      if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) continue;

      // Imagine the other way around, the zero map zero clusters track
      // is stored and the good one wants to be added. We ommit the warning
      // and just overwrite the 'bad' track
      if (fGlobalTrackReference[track->GetID()]->GetFilterMap() ||
          fGlobalTrackReference[track->GetID()]->GetTPCNcls()) {
        // If we come here, there's a problem
        std::cout << "Warning! global track info already there! ";
        std::cout << "TPCNcls track1 "
                  << (fGlobalTrackReference[track->GetID()])->GetTPCNcls()
                  << " track2 " << track->GetTPCNcls();
        std::cout << " FilterMap track1 "
                  << (fGlobalTrackReference[track->GetID()])->GetFilterMap()
                  << " track 2" << track->GetFilterMap() << "\n";
        fGlobalTrackReference[track->GetID()] = nullptr;
      }
    }  // Two tracks same id

    // Assign the pointer
    (fGlobalTrackReference.at(track->GetID())) = track;
  }
}

//____________________________________________________________________________________________________
bool AliSigma0SingleParticleCuts::SingleParticlePID(
    const AliVTrack *track) const {
  const AliVTrack *globaltrack =
      (fFilterbit == 128 && fInputEvent->IsA() == AliAODEvent::Class())
          ? fGlobalTrackReference[-track->GetID() - 1]
          : track;

  AliPIDResponse::EDetPidStatus statusPosTOF =
      fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, globaltrack);
  AliPIDResponse::EDetPidStatus statusPosTPC =
      fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, globaltrack);

  // TPC signal must be available
  if (!(AliPIDResponse::kDetPidOk == statusPosTPC)) return false;

  const float momentum = globaltrack->P();

  bool isPIDTOF = false;
  if (!fPIDTOFif && (momentum > fPIDMomentumSwitch)) {
    if (AliPIDResponse::kDetPidOk == statusPosTOF) {
      isPIDTOF = true;
    } else
      return false;
  }
  if (fPIDTOFif && AliPIDResponse::kDetPidOk == statusPosTOF) {
    isPIDTOF = true;
  }

  float nSigmaParticle = 999.f;
  const float nSigmaTPCParticle =
      fPIDResponse->NumberOfSigmasTPC(globaltrack, fParticle);
  const float nSigmaTOFParticle =
      fPIDResponse->NumberOfSigmasTOF(globaltrack, fParticle);
  nSigmaParticle = std::fabs(nSigmaTPCParticle);

  // TOF
  if (isPIDTOF)
    nSigmaParticle = std::sqrt(nSigmaTPCParticle * nSigmaTPCParticle +
                               nSigmaTOFParticle * nSigmaTOFParticle);

  if (fIsExtendedQA) {
    if (track->Charge() > 0) {
      fHistSingleParticleNsigmaBefore->Fill(momentum, nSigmaParticle);
      fHistSingleParticleNsigmaTPCBefore->Fill(momentum, nSigmaTPCParticle);
      fHistSingleParticleNsigmaTOFBefore->Fill(momentum, nSigmaTOFParticle);

      fHistSingleParticleTPCsignalBefore->Fill(momentum,
                                               globaltrack->GetTPCsignal());
      fHistSingleParticleTOFsignalBefore->Fill(momentum, GetBeta(globaltrack));
    } else {
      fHistSingleAntiParticleNsigmaBefore->Fill(momentum, nSigmaParticle);
      fHistSingleAntiParticleNsigmaTPCBefore->Fill(momentum, nSigmaTPCParticle);
      fHistSingleAntiParticleNsigmaTOFBefore->Fill(momentum, nSigmaTOFParticle);

      fHistSingleAntiParticleTPCsignalBefore->Fill(momentum,
                                                   globaltrack->GetTPCsignal());
      fHistSingleAntiParticleTOFsignalBefore->Fill(momentum,
                                                   GetBeta(globaltrack));
    }
  }

  // Proton PID cut
  if (nSigmaParticle > fPIDnSigmaStrict) return false;

  // check different hypotheses if TOF PID is available and  p > 0.75
  if (isPIDTOF && fPIDHypothesisRejection) {
    const float nSigmaTPCVeto1 =
        std::fabs(fPIDResponse->NumberOfSigmasTPC(globaltrack, fParticleVeto1));
    const float nSigmaTPCVeto2 =
        std::fabs(fPIDResponse->NumberOfSigmasTPC(globaltrack, fParticleVeto2));
    const float nSigmaTPCVeto3 =
        std::fabs(fPIDResponse->NumberOfSigmasTPC(globaltrack, fParticleVeto3));

    const float nSigmaTOFVeto1 =
        std::fabs(fPIDResponse->NumberOfSigmasTOF(globaltrack, fParticleVeto1));
    const float nSigmaTOFVeto2 =
        std::fabs(fPIDResponse->NumberOfSigmasTOF(globaltrack, fParticleVeto2));
    const float nSigmaTOFVeto3 =
        std::fabs(fPIDResponse->NumberOfSigmasTOF(globaltrack, fParticleVeto3));

    const float nSigmaVeto1Combined = std::sqrt(
        nSigmaTPCVeto1 * nSigmaTPCVeto1 + nSigmaTOFVeto1 * nSigmaTOFVeto1);
    const float nSigmaVeto2Combined = std::sqrt(
        nSigmaTPCVeto2 * nSigmaTPCVeto2 + nSigmaTOFVeto2 * nSigmaTOFVeto2);
    const float nSigmaVeto3Combined = std::sqrt(
        nSigmaTPCVeto3 * nSigmaTPCVeto3 + nSigmaTOFVeto3 * nSigmaTOFVeto3);

    if (nSigmaVeto1Combined < nSigmaParticle ||
        nSigmaVeto2Combined < nSigmaParticle ||
        nSigmaVeto3Combined < nSigmaParticle)
      return false;
  }

  if (!isPIDTOF && fPIDAdditionalRejection) {
    // write histos before cuts
    const float pT = globaltrack->Pt();
    for (int i = 0; i < static_cast<int>(fParticleAdditionalRejection.size());
         ++i) {
      if (pT < fPIDAdditionalRejectionLowerpTBound[i] ||
          pT > fPIDAdditionalRejectionUpperpTBound[i])
        continue;
      const float nSigmaTPCVeto = fPIDResponse->NumberOfSigmasTPC(
          globaltrack, fParticleAdditionalRejection[i]);
      if (track->Charge() > 0)
        fHistSingleParticleNsigmaTPCRejBefore[i]->Fill(momentum, nSigmaTPCVeto);
      else
        fHistSingleAntiParticleNsigmaTPCRejBefore[i]->Fill(momentum,
                                                           nSigmaTPCVeto);
    }

    for (int i = 0; i < static_cast<int>(fParticleAdditionalRejection.size());
         ++i) {
      const float pT = globaltrack->Pt();
      if (pT < fPIDAdditionalRejectionLowerpTBound[i] ||
          pT > fPIDAdditionalRejectionUpperpTBound[i])
        continue;
      const float nSigmaTPCVeto = fPIDResponse->NumberOfSigmasTPC(
          globaltrack, fParticleAdditionalRejection[i]);
      if (nSigmaTPCVeto > fPIDAdditionalRejectionLowerBound[i] &&
          nSigmaTPCVeto < fPIDAdditionalRejectionUpperBound[i])
        return false;
    }

    // write histos after cuts
    for (int i = 0; i < static_cast<int>(fParticleAdditionalRejection.size());
         ++i) {
      if (pT < fPIDAdditionalRejectionLowerpTBound[i] ||
          pT > fPIDAdditionalRejectionUpperpTBound[i])
        continue;
      const float nSigmaTPCVeto = fPIDResponse->NumberOfSigmasTPC(
          globaltrack, fParticleAdditionalRejection[i]);
      if (track->Charge() > 0)
        fHistSingleParticleNsigmaTPCRejAfter[i]->Fill(momentum, nSigmaTPCVeto);
      else
        fHistSingleAntiParticleNsigmaTPCRejAfter[i]->Fill(momentum,
                                                          nSigmaTPCVeto);
    }
  }

  if (fIsExtendedQA) {
    if (track->Charge() > 0) {
      fHistSingleParticleNsigmaAfter->Fill(momentum, nSigmaParticle);
      fHistSingleParticleNsigmaTPCAfter->Fill(momentum, nSigmaTPCParticle);
      fHistSingleParticleNsigmaTOFAfter->Fill(momentum, nSigmaTOFParticle);

      fHistSingleParticleTPCsignalAfter->Fill(momentum,
                                              globaltrack->GetTPCsignal());
      fHistSingleParticleTOFsignalAfter->Fill(momentum, GetBeta(globaltrack));
    } else {
      fHistSingleAntiParticleNsigmaAfter->Fill(momentum, nSigmaParticle);
      fHistSingleAntiParticleNsigmaTPCAfter->Fill(momentum, nSigmaTPCParticle);
      fHistSingleAntiParticleNsigmaTOFAfter->Fill(momentum, nSigmaTOFParticle);

      fHistSingleAntiParticleTPCsignalAfter->Fill(momentum,
                                                  globaltrack->GetTPCsignal());
      fHistSingleAntiParticleTOFsignalAfter->Fill(momentum,
                                                  GetBeta(globaltrack));
    }
  }
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0SingleParticleCuts::SingleParticleQualityCuts(
    AliVTrack *track) const {
  const int charge = track->Charge();
  const float eta = track->Eta();
  const float pt = track->Pt();
  const short nClsTPC = track->GetTPCNcls();
  const float nCrossedRows = track->GetTPCClusterInfo(2, 1);
  const short nFindable = track->GetTPCNclsF();
  const float ratioFindable = nCrossedRows / static_cast<float>(nFindable);

  if (charge > 0 && fIsExtendedQA) {
    fHistSingleParticleEtaBefore->Fill(eta);
    fHistSingleParticlePtBefore->Fill(pt);
    fHistSingleParticleNclsTPCBefore->Fill(nClsTPC);
    fHistSingleParticleNcrossedTPCBefore->Fill(nCrossedRows);
    fHistSingleParticleFindableTPCBefore->Fill(ratioFindable);
  } else if (charge < 0 && fIsExtendedQA) {
    fHistSingleAntiParticleEtaBefore->Fill(eta);
    fHistSingleAntiParticlePtBefore->Fill(pt);
    fHistSingleAntiParticleNclsTPCBefore->Fill(nClsTPC);
    fHistSingleAntiParticleNcrossedTPCBefore->Fill(nCrossedRows);
    fHistSingleAntiParticleFindableTPCBefore->Fill(ratioFindable);
  }

  // Eta
  if (std::fabs(eta) > fEtaMax) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(4);
  else
    fHistSingleAntiParticleCuts->Fill(4);

  // pT min & max
  if (pt < fPtMin) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(5);
  else
    fHistSingleAntiParticleCuts->Fill(5);
  if (pt > fPtMax) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(6);
  else
    fHistSingleAntiParticleCuts->Fill(6);

  // TPC shared clusters
  if (fInputEvent->IsA() == AliESDEvent::Class()) {
    if (!GoodTPCFitMapSharedMap(dynamic_cast<AliESDtrack *>(track)))
      return false;
  }
  if (fInputEvent->IsA() == AliAODEvent::Class()) {
    if (!GoodTPCFitMapSharedMap(static_cast<AliAODTrack *>(track)))
      return false;
  }
  if (charge > 0)
    fHistSingleParticleCuts->Fill(7);
  else
    fHistSingleAntiParticleCuts->Fill(7);

  // TPC nCls
  if (nClsTPC < fTPCclusterMin) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(8);
  else
    fHistSingleAntiParticleCuts->Fill(8);

  // TPC crossed rows
  if (nCrossedRows < fTPCnCrossedRowsMin) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(9);
  else
    fHistSingleAntiParticleCuts->Fill(9);

  // TPC findable clusters
  if (!track->GetTPCNclsF()) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(10);
  else
    fHistSingleAntiParticleCuts->Fill(10);

  // TPC ratio findable/found
  if (ratioFindable < fTPCfindableMin) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(11);
  else
    fHistSingleAntiParticleCuts->Fill(11);

  // For ESDs: TPC refit and Kinks
  if (fInputEvent->IsA() == AliESDEvent::Class()) {
    AliESDtrack *esdTrack = static_cast<AliESDtrack *>(track);
    if (!(esdTrack->GetStatus() & AliESDtrack::kTPCrefit)) return false;
    if (charge > 0)
      fHistSingleParticleCuts->Fill(12);
    else
      fHistSingleAntiParticleCuts->Fill(12);
    if (esdTrack->GetKinkIndex(0) > 0) return false;
    if (charge > 0)
      fHistSingleParticleCuts->Fill(13);
    else
      fHistSingleAntiParticleCuts->Fill(13);
  }

  const int nClsSharedTPC =
      (fInputEvent->IsA() == AliESDEvent::Class())
          ? (static_cast<AliESDtrack *>(track))->GetTPCnclsS()
          : (static_cast<AliAODTrack *>(track))->GetTPCnclsS();
  int nClsSharedITS = 0;
  for (int i = 0; i < 6; ++i) {
    if (track->HasSharedPointOnITSLayer(i)) ++nClsSharedITS;
  }

  if (charge > 0 && fIsExtendedQA) {
    fHistSingleParticleEtaAfter->Fill(eta);
    fHistSingleParticlePtAfter->Fill(pt);
    fHistSingleParticleNclsTPCAfter->Fill(nClsTPC);
    fHistSingleParticleNcrossedTPCAfter->Fill(nCrossedRows);
    fHistSingleParticleFindableTPCAfter->Fill(ratioFindable);
    fHistSingleParticleNclsTPCShared->Fill(nClsSharedTPC);
    fHistSingleParticleNclsTPCSharedTiming->Fill(0.f, nClsSharedTPC);
    fHistSingleParticleNclsITSShared->Fill(nClsSharedITS);
    fHistSingleParticleNclsITSSharedTiming->Fill(0.f, nClsSharedITS);
    if (track->HasPointOnITSLayer(0)) {
      fHistSingleParticleNclsTPCSharedTiming->Fill(1.f, nClsSharedTPC);
      fHistSingleParticleNclsITSSharedTiming->Fill(1.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(1)) {
      fHistSingleParticleNclsTPCSharedTiming->Fill(2.f, nClsSharedTPC);
      fHistSingleParticleNclsITSSharedTiming->Fill(2.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(2)) {
      fHistSingleParticleNclsTPCSharedTiming->Fill(3.f, nClsSharedTPC);
      fHistSingleParticleNclsITSSharedTiming->Fill(3.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(3)) {
      fHistSingleParticleNclsTPCSharedTiming->Fill(4.f, nClsSharedTPC);
      fHistSingleParticleNclsITSSharedTiming->Fill(4.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(4)) {
      fHistSingleParticleNclsTPCSharedTiming->Fill(5.f, nClsSharedTPC);
      fHistSingleParticleNclsITSSharedTiming->Fill(5.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(5)) {
      fHistSingleParticleNclsTPCSharedTiming->Fill(6.f, nClsSharedTPC);
      fHistSingleParticleNclsITSSharedTiming->Fill(6.f, nClsSharedITS);
    }
    if (track->GetTOFBunchCrossing() == 0) {
      fHistSingleParticleNclsTPCSharedTiming->Fill(7.f, nClsSharedTPC);
      fHistSingleParticleNclsITSSharedTiming->Fill(7.f, nClsSharedITS);
    }
  } else if (charge < 0 && fIsExtendedQA) {
    fHistSingleAntiParticleEtaAfter->Fill(eta);
    fHistSingleAntiParticlePtAfter->Fill(pt);
    fHistSingleAntiParticleNclsTPCAfter->Fill(nClsTPC);
    fHistSingleAntiParticleNcrossedTPCAfter->Fill(nCrossedRows);
    fHistSingleAntiParticleFindableTPCAfter->Fill(ratioFindable);
    fHistSingleAntiParticleNclsTPCShared->Fill(nClsSharedTPC);
    fHistSingleAntiParticleNclsTPCSharedTiming->Fill(0.f, nClsSharedTPC);
    fHistSingleAntiParticleNclsITSShared->Fill(nClsSharedITS);
    fHistSingleAntiParticleNclsITSSharedTiming->Fill(0.f, nClsSharedITS);
    if (track->HasPointOnITSLayer(0)) {
      fHistSingleAntiParticleNclsTPCSharedTiming->Fill(1.f, nClsSharedTPC);
      fHistSingleAntiParticleNclsITSSharedTiming->Fill(1.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(1)) {
      fHistSingleAntiParticleNclsTPCSharedTiming->Fill(2.f, nClsSharedTPC);
      fHistSingleAntiParticleNclsITSSharedTiming->Fill(2.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(2)) {
      fHistSingleAntiParticleNclsTPCSharedTiming->Fill(3.f, nClsSharedTPC);
      fHistSingleAntiParticleNclsITSSharedTiming->Fill(3.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(3)) {
      fHistSingleAntiParticleNclsTPCSharedTiming->Fill(4.f, nClsSharedTPC);
      fHistSingleAntiParticleNclsITSSharedTiming->Fill(4.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(4)) {
      fHistSingleAntiParticleNclsTPCSharedTiming->Fill(5.f, nClsSharedTPC);
      fHistSingleAntiParticleNclsITSSharedTiming->Fill(5.f, nClsSharedITS);
    }
    if (track->HasPointOnITSLayer(5)) {
      fHistSingleAntiParticleNclsTPCSharedTiming->Fill(6.f, nClsSharedTPC);
      fHistSingleAntiParticleNclsITSSharedTiming->Fill(6.f, nClsSharedITS);
    }
    if (track->GetTOFBunchCrossing() == 0) {
      fHistSingleAntiParticleNclsTPCSharedTiming->Fill(7.f, nClsSharedTPC);
      fHistSingleAntiParticleNclsITSSharedTiming->Fill(7.f, nClsSharedITS);
    }
  }

  return true;
}

//____________________________________________________________________________________________________
template <typename T>
bool AliSigma0SingleParticleCuts::GoodTPCFitMapSharedMap(const T *track) const {
  // This method was inherited form H. Beck analysis

  // Rejects tracks with shared clusters
  // This overload is used for primaries

  // Get the shared maps
  const TBits sharedMap = track->GetTPCSharedMap();
  if ((sharedMap.CountBits()) >= 1) {
    // Bad track, has too many shared clusters!
    return false;
  }
  return true;
}

//____________________________________________________________________________________________________
bool AliSigma0SingleParticleCuts::SingleParticleDCASelection(
    AliVTrack *track, float &DCAr, float &DCAz) const {
  if (fInputEvent->IsA() == AliESDEvent::Class()) {
    const AliESDtrack *esdTrack = static_cast<AliESDtrack *>(track);
    esdTrack->GetImpactParameters(DCAr, DCAz);
  }
  if (fInputEvent->IsA() == AliAODEvent::Class()) {
    // in case of Filterbit 128 (TPConly track) use the DCA of the global track
    const AliAODTrack *aodTrack = static_cast<AliAODTrack *>(track);
    const AliAODTrack *globalTrack =
        (fFilterbit == 128) ? fGlobalTrackReference[-track->GetID() - 1]
                            : aodTrack;
    if (!globalTrack) return false;

    // Create an external parameter from the AODtrack
    AliExternalTrackParam globalTrackParams;
    globalTrackParams.CopyFromVTrack(globalTrack);

    Double_t *DCAVals = new Double_t[2];
    Double_t covar[3] = {
        0.,
        0.,
        0.,
    };
    if (!globalTrackParams.PropagateToDCA(fInputEvent->GetPrimaryVertex(),
                                          fInputEvent->GetMagneticField(), 10.,
                                          DCAVals, covar))
      return false;

    DCAr = DCAVals[0];
    DCAz = DCAVals[1];
    delete[] DCAVals;
  }

  const int charge = track->Charge();
  const float pt = track->Pt();

  if (charge > 0 && fIsExtendedQA) {
    fHistSingleParticleDCAxyBefore->Fill(DCAr);
    fHistSingleParticleDCAzBefore->Fill(DCAz);
  } else if (charge < 0 && fIsExtendedQA) {
    fHistSingleAntiParticleDCAxyBefore->Fill(DCAr);
    fHistSingleAntiParticleDCAzBefore->Fill(DCAz);
  }

  // DCA z
  if (std::fabs(DCAz) > fDCAzMax) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(14);
  else
    fHistSingleAntiParticleCuts->Fill(14);
  if (charge > 0 && fIsExtendedQA)
    fHistSingleParticleDCAxyAfterDCAz->Fill(DCAr);
  else if (charge < 0 && fIsExtendedQA)
    fHistSingleAntiParticleDCAxyAfterDCAz->Fill(DCAr);

  if (charge > 0) {
    fHistSingleParticleDCAxy->Fill(pt, DCAr);
  } else {
    fHistSingleAntiParticleDCAxy->Fill(pt, DCAr);
  }

  // DCA r
  if (std::fabs(DCAr) > fDCArMax) return false;
  if (charge > 0)
    fHistSingleParticleCuts->Fill(15);
  else
    fHistSingleAntiParticleCuts->Fill(15);

  if (charge > 0 && fIsExtendedQA) {
    fHistSingleParticleDCAxyAfter->Fill(DCAr);
    fHistSingleParticleDCAzAfter->Fill(DCAz);
  } else if (charge < 0 && fIsExtendedQA) {
    fHistSingleAntiParticleDCAzAfter->Fill(DCAz);
    fHistSingleAntiParticleDCAxyAfter->Fill(DCAr);
  }

  return true;
}

//____________________________________________________________________________________________________
void AliSigma0SingleParticleCuts::ProcessMC(AliSigma0ParticleBase *baseParticle,
                                            AliMCParticle *mcParticle) const {
  if (baseParticle->GetCharge() > 0) {
    fHistMCRecSingleParticleMomentum->Fill(baseParticle->GetPtMC(),
                                           baseParticle->GetPt());
    if (std::abs(mcParticle->PdgCode()) == fPDGCode)
      fHistMCRecSingleParticlePtTruth->Fill(baseParticle->GetPt());
  } else {
    fHistMCRecSingleAntiParticleMomentum->Fill(baseParticle->GetPtMC(),
                                               baseParticle->GetPt());
    if (std::abs(mcParticle->PdgCode()) == fPDGCode)
      fHistMCRecSingleAntiParticlePtTruth->Fill(baseParticle->GetPt());
  }
}

//____________________________________________________________________________________________________
float AliSigma0SingleParticleCuts::GetBeta(const AliVTrack *track) const {
  float beta = -999;
  Double_t integratedTimes[9] = {-1.0, -1.0, -1.0, -1.0, -1.0,
                                 -1.0, -1.0, -1.0, -1.0};

  track->GetIntegratedTimes(integratedTimes);

  const float c = 2.99792457999999984e-02;
  float p = track->P();
  float l = integratedTimes[0] * c;

  float trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);

  float timeTOF = track->GetTOFsignal() - trackT0;
  if (timeTOF > 0) beta = l / timeTOF / c;
  return beta;
}

//____________________________________________________________________________________________________
void AliSigma0SingleParticleCuts::InitCutHistograms() {
  std::cout << "============================\n"
            << " SINGLE PARTICLE CUT CONFIGURATION \n"
            << " Filterbit     " << fFilterbit << "\n"
            << " nCls min      " << fTPCclusterMin << "\n"
            << " nCrossed rows " << fTPCnCrossedRowsMin << "\n"
            << " Findable min  " << fTPCfindableMin << "\n"
            << " Shared map    " << fUseSharedMap << "\n"
            << " DCA r max     " << fDCArMax << "\n"
            << " DCA z max     " << fDCAzMax << "\n"
            << " p_T min       " << fPtMin << "\n"
            << " p_T max       " << fPtMax << "\n"
            << " eta max       " << fEtaMax << "\n"
            << " PID nsigma    " << fPIDnSigmaStrict << "\n"
            << " PID rejection " << fPIDHypothesisRejection << "\n"
            << " PID p switch  " << fPIDMomentumSwitch << "\n"
            << " Add. reject.  " << fPIDAdditionalRejection << "\n";
  if (fPIDAdditionalRejection) {
    for (int i = 0; i < static_cast<int>(fParticleAdditionalRejection.size());
         ++i) {
      std::cout << " - " << fParticleAdditionalRejection[i] << " lower "
                << fPIDAdditionalRejectionLowerBound[i] << " upper "
                << fPIDAdditionalRejectionUpperBound[i] << " for "
                << fPIDAdditionalRejectionLowerpTBound[i] << " < pT < "
                << fPIDAdditionalRejectionUpperpTBound[i] << "\n";
    }
  }
  std::cout << "============================\n";

  TH1::AddDirectory(kFALSE);

  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    TString name = "SingleParticleCut_QA_" + fParticleName;
    fHistograms->SetName(name);
  }

  if (fHistogramsSingleParticle != nullptr) {
    delete fHistogramsSingleParticle;
    fHistogramsSingleParticle = nullptr;
  }
  if (fHistogramsSingleParticle == nullptr) {
    fHistogramsSingleParticle = new TList();
    fHistogramsSingleParticle->SetOwner(kTRUE);
    fHistogramsSingleParticle->SetName("Particle_QA");
  }

  if (fHistogramsSingleAntiParticle != nullptr) {
    delete fHistogramsSingleAntiParticle;
    fHistogramsSingleAntiParticle = nullptr;
  }
  if (fHistogramsSingleAntiParticle == nullptr) {
    fHistogramsSingleAntiParticle = new TList();
    fHistogramsSingleAntiParticle->SetOwner(kTRUE);
    fHistogramsSingleAntiParticle->SetName("AntiParticle_QA");
  }

  const float PI = std::atan(1.0) * 4.f;

  fHistCuts = new TProfile("fHistCuts", ";;Cut value", 14, 0, 14);
  fHistCuts->GetXaxis()->SetBinLabel(1, "Filter bit");
  fHistCuts->GetXaxis()->SetBinLabel(2, "TPC nCls min");
  fHistCuts->GetXaxis()->SetBinLabel(3, "TPC ncrossed rows min");
  fHistCuts->GetXaxis()->SetBinLabel(4, "TPC findable min");
  fHistCuts->GetXaxis()->SetBinLabel(5, "TPC use shared map");
  fHistCuts->GetXaxis()->SetBinLabel(6, "DCAr max");
  fHistCuts->GetXaxis()->SetBinLabel(7, "DCAz max");
  fHistCuts->GetXaxis()->SetBinLabel(8, "#it{p}_{T} min");
  fHistCuts->GetXaxis()->SetBinLabel(9, "#it{p}_{T} max");
  fHistCuts->GetXaxis()->SetBinLabel(10, "#eta max");
  fHistCuts->GetXaxis()->SetBinLabel(11, "PID n#sigma");
  fHistCuts->GetXaxis()->SetBinLabel(12, "PID #it{p} switch");
  fHistCuts->GetXaxis()->SetBinLabel(13, "PID rejection");
  fHistCuts->GetXaxis()->SetBinLabel(14, "Additional rejection");
  fHistCuts->GetXaxis()->SetBinLabel(15, "PDG code");
  fHistograms->Add(fHistCuts);

  fHistCuts->Fill(0.f, fFilterbit);
  fHistCuts->Fill(1.f, fTPCclusterMin);
  fHistCuts->Fill(2.f, fTPCnCrossedRowsMin);
  fHistCuts->Fill(3.f, fTPCfindableMin);
  fHistCuts->Fill(4.f, static_cast<double>(fUseSharedMap));
  fHistCuts->Fill(5.f, fDCArMax);
  fHistCuts->Fill(6.f, fDCAzMax);
  fHistCuts->Fill(7.f, fPtMin);
  fHistCuts->Fill(8.f, fPtMax);
  fHistCuts->Fill(9.f, fEtaMax);
  fHistCuts->Fill(10.f, fPIDnSigmaStrict);
  fHistCuts->Fill(11.f, fPIDMomentumSwitch);
  fHistCuts->Fill(12.f, static_cast<double>(fPIDHypothesisRejection));
  fHistCuts->Fill(13.f, static_cast<double>(fPIDAdditionalRejection));
  fHistCuts->Fill(14.f, static_cast<double>(fPDGCode));

  fHistSingleParticleCuts =
      new TH1F("fHistSingleParticleCuts", ";;Entries", 19, 0, 19);
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(1, "Tracks");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(
      2, Form("Filterbit %i", fFilterbit));
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(3, "Global ref");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(4, "PID");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(5, "Eta");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(6, "#it{p}_{T} min");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(7, "#it{p}_{T} max");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(8, "TPC shared map");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(9, "TPC cluster");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(10, "TPC rows");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(11, "?TPC findable");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(12, "TPC find. cut");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(13, "ESD: TPC refit");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(14, "ESD: Kink");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(15, "DCA z");
  fHistSingleParticleCuts->GetXaxis()->SetBinLabel(16, "DCA r");

  fHistSingleParticlePt =
      new TH1F("fHistSingleParticlePt", ";#it{p}_{T}; Entries", 500, 0, 10);
  fHistSingleParticlePhi = new TH1F(
      "fHistSingleParticlePhi", "; #phi [rad]; Entries", 100 * PI, 0, 2 * PI);
  fHistSingleParticleEta =
      new TH1F("fHistSingleParticleEta", "; #eta; Entries", 500, -2, 2);
  fHistNSingleParticle =
      new TH1F("fHistNSingleParticle", "; # particles; Entries", 50, 0, 50);
  fHistSingleParticleDCAxy = new TH2F(
      "fHistSingleParticleDCAxy", "; #it{p}_{T} [GeV/#it{c}]; DCA_{xy} [cm]",
      1000, 0, 10, 500, -5, 5);
  fHistSingleParticleEtaPhi =
      new TH2F("fHistSingleParticleEtaPhi", "; #eta; #phi", 200, -1, 1, 200, 0,
               2 * TMath::Pi());

  fHistogramsSingleParticle->Add(fHistSingleParticleCuts);
  fHistogramsSingleParticle->Add(fHistSingleParticlePt);
  fHistogramsSingleParticle->Add(fHistSingleParticlePhi);
  fHistogramsSingleParticle->Add(fHistSingleParticleEta);
  fHistogramsSingleParticle->Add(fHistNSingleParticle);
  fHistogramsSingleParticle->Add(fHistSingleParticleDCAxy);
  fHistogramsSingleParticle->Add(fHistSingleParticleEtaPhi);

  fHistSingleAntiParticleCuts =
      new TH1F("fHistSingleAntiParticleCuts", ";;Entries", 19, 0, 19);
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(1, "Tracks");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(2, "Filterbit 128");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(3, "Global ref");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(4, "PID");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(5, "Eta");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(6, "#it{p}_{T} min");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(7, "#it{p}_{T} max");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(8, "TPC shared map");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(9, "TPC cluster");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(10, "TPC rows");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(11, "?TPC findable");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(12, "TPC find. cut");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(13, "ESD: TPC refit");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(14, "ESD: Kink");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(15, "DCA z");
  fHistSingleAntiParticleCuts->GetXaxis()->SetBinLabel(16, "DCA r");

  fHistSingleAntiParticlePt =
      new TH1F("fHistSingleAntiParticlePt", ";#it{p}_{T}; Entries", 500, 0, 10);
  fHistSingleAntiParticlePhi =
      new TH1F("fHistSingleAntiParticlePhi", "; #phi [rad]; Entries", 100 * PI,
               0, 2 * PI);
  fHistSingleAntiParticleEta =
      new TH1F("fHistSingleAntiParticleEta", "; #eta; Entries", 500, -2, 2);
  fHistNSingleAntiParticle = new TH1F("fHistNSingleAntiParticle",
                                      "; # anti-particles; Entries", 50, 0, 50);
  fHistSingleAntiParticleDCAxy = new TH2F(
      "fHistSingleAntiParticleDCAxy",
      "; #it{p}_{T} [GeV/#it{c}]; DCA_{xy} [cm]", 1000, 0, 10, 500, -5, 5);
  fHistSingleAntiParticleEtaPhi =
      new TH2F("fHistSingleAntiParticleEtaPhi", "; #eta; #phi", 200, -1, 1, 200,
               0, 2 * TMath::Pi());

  fHistogramsSingleAntiParticle->Add(fHistSingleAntiParticleCuts);
  fHistogramsSingleAntiParticle->Add(fHistSingleAntiParticlePt);
  fHistogramsSingleAntiParticle->Add(fHistSingleAntiParticlePhi);
  fHistogramsSingleAntiParticle->Add(fHistSingleAntiParticleEta);
  fHistogramsSingleAntiParticle->Add(fHistNSingleAntiParticle);
  fHistogramsSingleAntiParticle->Add(fHistSingleAntiParticleDCAxy);
  fHistogramsSingleAntiParticle->Add(fHistSingleAntiParticleEtaPhi);

  if (fIsMC) {
    if (fHistogramsSingleParticleMC != nullptr) {
      delete fHistogramsSingleParticleMC;
      fHistogramsSingleParticleMC = nullptr;
    }
    if (fHistogramsSingleParticleMC == nullptr) {
      fHistogramsSingleParticleMC = new TList();
      fHistogramsSingleParticleMC->SetOwner(kTRUE);
      fHistogramsSingleParticleMC->SetName("Particle_MC");
    }

    if (fHistogramsSingleAntiParticleMC != nullptr) {
      delete fHistogramsSingleAntiParticleMC;
      fHistogramsSingleAntiParticleMC = nullptr;
    }
    if (fHistogramsSingleAntiParticleMC == nullptr) {
      fHistogramsSingleAntiParticleMC = new TList();
      fHistogramsSingleAntiParticleMC->SetOwner(kTRUE);
      fHistogramsSingleAntiParticleMC->SetName("AntiParticle_MC");
    }

    fHistMCRecSingleParticlePtTruth =
        new TH1F("fHistMCRecSingleParticlePtTruth",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCRecSingleParticleMomentum = new TH2F(
        "fHistMCRecSingleParticleMomentum",
        "; #it{p}_{T, truth} [GeV/#it{c}]; #it{p}_{T, rec} [GeV/#it{c}]", 2000,
        0, 2, 2000, 0, 2);
    fHistogramsSingleParticleMC->Add(fHistMCRecSingleParticlePtTruth);
    fHistogramsSingleParticleMC->Add(fHistMCRecSingleParticleMomentum);

    fHistMCRecSingleAntiParticlePtTruth =
        new TH1F("fHistMCRecSingleAntiParticlePtTruth",
                 "; #it{p}_{T} [GeV/#it{c}]; Entries", 500, 0, 10);
    fHistMCRecSingleAntiParticleMomentum = new TH2F(
        "fHistMCRecSingleAntiParticleMomentum",
        "; #it{p}_{T, truth} [GeV/#it{c}]; #it{p}_{T, rec} [GeV/#it{c}]", 2000,
        0, 2, 2000, 0, 2);
    fHistogramsSingleAntiParticleMC->Add(fHistMCRecSingleAntiParticlePtTruth);
    fHistogramsSingleAntiParticleMC->Add(fHistMCRecSingleAntiParticleMomentum);

    fHistogramsSingleParticle->Add(fHistogramsSingleParticleMC);
    fHistogramsSingleAntiParticle->Add(fHistogramsSingleAntiParticleMC);
  }

  if (fIsExtendedQA) {
    if (fHistogramsSingleParticleBefore != nullptr) {
      delete fHistogramsSingleParticleBefore;
      fHistogramsSingleParticleBefore = nullptr;
    }
    if (fHistogramsSingleParticleBefore == nullptr) {
      fHistogramsSingleParticleBefore = new TList();
      fHistogramsSingleParticleBefore->SetOwner(kTRUE);
      fHistogramsSingleParticleBefore->SetName("Before");
    }

    if (fHistogramsSingleParticleAfter != nullptr) {
      delete fHistogramsSingleParticleAfter;
      fHistogramsSingleParticleAfter = nullptr;
    }
    if (fHistogramsSingleParticleAfter == nullptr) {
      fHistogramsSingleParticleAfter = new TList();
      fHistogramsSingleParticleAfter->SetOwner(kTRUE);
      fHistogramsSingleParticleAfter->SetName("After");
    }

    fHistSingleParticleEtaBefore =
        new TH1F("fHistSingleParticleEtaBefore", "; #eta (before); Entries",
                 1000, -2, 2);
    fHistSingleParticleEtaAfter = new TH1F(
        "fHistSingleParticleEtaAfter", "; #eta (after); Entries", 1000, -2, 2);
    fHistSingleParticlePtBefore =
        new TH1F("fHistSingleParticlePtBefore", ";#it{p}_{T} (before); Entries",
                 1000, 0, 10);
    fHistSingleParticlePtAfter =
        new TH1F("fHistSingleParticlePtAfter", "#it{p}_{T} (after); Entries",
                 1000, 0, 10);
    fHistSingleParticleNclsTPCBefore =
        new TH1F("fHistSingleParticleNclsTPCBefore",
                 "; # cls TPC (before); Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCAfter =
        new TH1F("fHistSingleParticleNclsTPCAfter",
                 "; # cls TPC (after); Entries", 170, 0, 170);

    fHistSingleParticleNclsTPCShared =
        new TH1F("fHistSingleParticleNclsTPCShared",
                 "; # cls TPC shared; Entries", 170, 0, 170);
    fHistSingleParticleNclsTPCSharedTiming =
        new TH2F("fHistSingleParticleNclsTPCSharedTiming",
                 "; Timing info; # cls TPC shared", 8, 0, 8, 170, 0, 170);
    fHistSingleParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(1, "all");
    fHistSingleParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistSingleParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistSingleParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistSingleParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistSingleParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistSingleParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistSingleParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistSingleParticleNclsITSShared =
        new TH1F("fHistSingleParticleNclsITSShared",
                 "; # cls ITS shared; Entries", 6, 0, 6);
    fHistSingleParticleNclsITSSharedTiming =
        new TH2F("fHistSingleParticleNclsITSSharedTiming",
                 "; Timing info; # cls ITS shared", 8, 0, 8, 6, 0, 6);
    fHistSingleParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(1, "all");
    fHistSingleParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistSingleParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistSingleParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistSingleParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistSingleParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistSingleParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistSingleParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        8, "TOF timing");

    fHistSingleParticleNcrossedTPCBefore =
        new TH1F("fHistSingleParticleNcrossedTPCBefore",
                 ";# crossed rows TPC (before); Entries", 170, 0, 170);
    fHistSingleParticleNcrossedTPCAfter =
        new TH1F("fHistSingleParticleNcrossedTPCAfter",
                 ";# crossed rows TPC (after); Entries", 170, 0, 170);
    fHistSingleParticleFindableTPCBefore =
        new TH1F("fHistSingleParticleFindableTPCBefore",
                 ";Findable TPC (before); Entries", 1000, 0, 2);
    fHistSingleParticleFindableTPCAfter =
        new TH1F("fHistSingleParticleFindableTPCAfter",
                 "; Findable TPC (after); Entries) ", 1000, 0, 2);
    fHistSingleParticleDCAxyBefore =
        new TH1F("fHistSingleParticleDCAxyBefore",
                 "; DCA_{xy} (before); Entries", 1000, -10, 10);
    fHistSingleParticleDCAxyAfterDCAz =
        new TH1F("fHistSingleParticleDCAxyAfterDCAz",
                 "; DCA_{xy} (after DCA_{z} cut}); Entries", 1000, -10, 10);
    fHistSingleParticleDCAxyAfter =
        new TH1F("fHistSingleParticleDCAxyAfter", "; DCA_{xy} (after); Entries",
                 1000, -10, 10);
    fHistSingleParticleDCAzBefore =
        new TH1F("fHistSingleParticleDCAzBefore", ";DCA_{z} (before); Entries",
                 1000, -10, 10);
    fHistSingleParticleDCAzAfter =
        new TH1F("fHistSingleParticleDCAzAfter", ";DCA_{z} (after); Entries",
                 1000, -10, 10);
    fHistSingleParticleNsigmaBefore =
        new TH2F("fHistSingleParticleNsigmaBefore",
                 "; p [GeV/c]; n_{#sigma} before", 1000, 0, 5, 1000, 0, 10);
    fHistSingleParticleNsigmaAfter =
        new TH2F("fHistSingleParticleNsigmaAfter",
                 "; p [GeV/c]; n_{#sigma} after", 1000, 0, 5, 1000, 0, 10);
    fHistSingleParticleNsigmaTPCBefore =
        new TH2F("fHistSingleParticleNsigmaTPCBefore",
                 "; p [GeV/c]; n_{#sigma} TPC before", 1000, 0, 5, 1000, -5, 5);
    fHistSingleParticleNsigmaTPCAfter =
        new TH2F("fHistSingleParticleNsigmaTPCAfter",
                 "; p [GeV/c]; n_{#sigma} TPC after", 1000, 0, 5, 1000, -5, 5);
    fHistSingleParticleNsigmaTOFBefore =
        new TH2F("fHistSingleParticleNsigmaTOFBefore",
                 "; p [GeV/c]; n_{#sigma} TOF before", 1000, 0, 5, 1000, -5, 5);
    fHistSingleParticleNsigmaTOFAfter =
        new TH2F("fHistSingleParticleNsigmaTOFAfter",
                 "; p [GeV/c]; n_{#sigma} TOF after", 1000, 0, 5, 1000, -5, 5);

    fHistSingleParticleTPCsignalBefore = new TH2F(
        "fHistSingleParticleTPCsignalBefore",
        "; p [GeV/c]; TPC d#it{E}/d#it{x} before", 500, 0, 5, 500, 0, 200);
    fHistSingleParticleTOFsignalBefore =
        new TH2F("fHistSingleParticleTOFsignalBefore",
                 "; p [GeV/c]; TOF #beta before", 500, 0, 5, 500, 0, 1);
    fHistSingleParticleTPCsignalAfter = new TH2F(
        "fHistSingleParticleTPCsignalAfter",
        "; p [GeV/c];  TPC d#it{E}/d#it{x} after", 500, 0, 5, 500, 0, 200);
    fHistSingleParticleTOFsignalAfter =
        new TH2F("fHistSingleParticleTOFsignalAfter",
                 "; p [GeV/c]; TOF #beta after", 500, 0, 5, 500, 0, 1);

    const int nPIDbins = fParticleAdditionalRejection.size();
    TH2F *singleParticlePIDrejBefore[nPIDbins];
    TH2F *singleParticlePIDrejAfter[nPIDbins];
    for (int i = 0; i < static_cast<int>(fParticleAdditionalRejection.size());
         ++i) {
      singleParticlePIDrejBefore[i] =
          new TH2F(Form("fHistSingleParticlePIDRejBefore_%i_pt_%.1f_%1.1f",
                        fParticleAdditionalRejection[i],
                        fPIDAdditionalRejectionLowerpTBound[i],
                        fPIDAdditionalRejectionUpperpTBound[i]),
                   Form("; p [GeV/c]; n_{#sigma} %i before",
                        fParticleAdditionalRejection[i]),
                   1000, 0, 5, 1000, -5, 5);
      fHistSingleParticleNsigmaTPCRejBefore.push_back(
          singleParticlePIDrejBefore[i]);
      singleParticlePIDrejAfter[i] =
          new TH2F(Form("fHistSingleParticlePIDRejAfter_%i_pt_%.1f_%1.1f",
                        fParticleAdditionalRejection[i],
                        fPIDAdditionalRejectionLowerpTBound[i],
                        fPIDAdditionalRejectionUpperpTBound[i]),
                   Form("; p [GeV/c]; n_{#sigma} %i before",
                        fParticleAdditionalRejection[i]),
                   1000, 0, 5, 1000, -5, 5);
      fHistSingleParticleNsigmaTPCRejAfter.push_back(
          singleParticlePIDrejAfter[i]);
    }

    fHistogramsSingleParticleBefore->Add(fHistSingleParticleEtaBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleEtaAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticlePtBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticlePtAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleNclsTPCBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNclsTPCAfter);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNclsTPCShared);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNclsTPCSharedTiming);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNclsITSShared);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNclsITSSharedTiming);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleNcrossedTPCBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNcrossedTPCAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleFindableTPCBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleFindableTPCAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleDCAxyBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleDCAxyAfterDCAz);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleDCAxyAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleDCAzBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleDCAzAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleNsigmaBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNsigmaAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleNsigmaTPCBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNsigmaTPCAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleNsigmaTOFBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleNsigmaTOFAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleTPCsignalBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleTPCsignalAfter);
    fHistogramsSingleParticleBefore->Add(fHistSingleParticleTOFsignalBefore);
    fHistogramsSingleParticleAfter->Add(fHistSingleParticleTOFsignalAfter);

    fHistogramsSingleParticle->Add(fHistogramsSingleParticleBefore);
    fHistogramsSingleParticle->Add(fHistogramsSingleParticleAfter);

    for (int i = 0;
         i < static_cast<int>(fHistSingleParticleNsigmaTPCRejBefore.size());
         ++i) {
      fHistogramsSingleParticleBefore->Add(
          fHistSingleParticleNsigmaTPCRejBefore[i]);
      fHistogramsSingleParticleAfter->Add(
          fHistSingleParticleNsigmaTPCRejAfter[i]);
    }

    if (fHistogramsSingleAntiParticleBefore != nullptr) {
      delete fHistogramsSingleAntiParticleBefore;
      fHistogramsSingleAntiParticleBefore = nullptr;
    }
    if (fHistogramsSingleAntiParticleBefore == nullptr) {
      fHistogramsSingleAntiParticleBefore = new TList();
      fHistogramsSingleAntiParticleBefore->SetOwner(kTRUE);
      fHistogramsSingleAntiParticleBefore->SetName("Before");
    }

    if (fHistogramsSingleAntiParticleAfter != nullptr) {
      delete fHistogramsSingleAntiParticleAfter;
      fHistogramsSingleAntiParticleAfter = nullptr;
    }
    if (fHistogramsSingleAntiParticleAfter == nullptr) {
      fHistogramsSingleAntiParticleAfter = new TList();
      fHistogramsSingleAntiParticleAfter->SetOwner(kTRUE);
      fHistogramsSingleAntiParticleAfter->SetName("After");
    }

    fHistSingleAntiParticleEtaBefore =
        new TH1F("fHistSingleAntiParticleEtaBefore", "; #eta (before); Entries",
                 1000, -2, 2);
    fHistSingleAntiParticleEtaAfter =
        new TH1F("fHistSingleAntiParticleEtaAfter", "; #eta (after); Entries",
                 1000, -2, 2);
    fHistSingleAntiParticlePtBefore =
        new TH1F("fHistSingleAntiParticlePtBefore",
                 ";#it{p}_{T} (before); Entries", 1000, 0, 10);
    fHistSingleAntiParticlePtAfter =
        new TH1F("fHistSingleAntiParticlePtAfter",
                 "#it{p}_{T} (after); Entries", 1000, 0, 10);
    fHistSingleAntiParticleNclsTPCBefore =
        new TH1F("fHistSingleAntiParticleNclsTPCBefore",
                 "; # cls TPC (before); Entries", 170, 0, 170);
    fHistSingleAntiParticleNclsTPCAfter =
        new TH1F("fHistSingleAntiParticleNclsTPCAfter",
                 "; # cls TPC (after); Entries", 170, 0, 170);

    fHistSingleAntiParticleNclsTPCShared =
        new TH1F("fHistSingleAntiParticleNclsTPCShared",
                 "; # cls TPC shared; Entries", 170, 0, 170);
    fHistSingleAntiParticleNclsTPCSharedTiming =
        new TH2F("fHistSingleAntiParticleNclsTPCSharedTiming",
                 "; Timing info; # cls TPC shared", 8, 0, 8, 170, 0, 170);
    fHistSingleAntiParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(1,
                                                                        "all");
    fHistSingleAntiParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistSingleAntiParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistSingleAntiParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistSingleAntiParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistSingleAntiParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistSingleAntiParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistSingleAntiParticleNclsTPCSharedTiming->GetXaxis()->SetBinLabel(
        8, "TOF timing");
    fHistSingleAntiParticleNclsITSShared =
        new TH1F("fHistSingleAntiParticleNclsITSShared",
                 "; # cls ITS shared; Entries", 6, 0, 6);
    fHistSingleAntiParticleNclsITSSharedTiming =
        new TH2F("fHistSingleAntiParticleNclsITSSharedTiming",
                 "; Timing info; # cls ITS shared", 8, 0, 8, 6, 0, 6);
    fHistSingleAntiParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(1,
                                                                        "all");
    fHistSingleAntiParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        2, "Hit ITS 0");
    fHistSingleAntiParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        3, "Hit ITS 1");
    fHistSingleAntiParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        4, "Hit ITS 2");
    fHistSingleAntiParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        5, "Hit ITS 3");
    fHistSingleAntiParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        6, "Hit ITS 4");
    fHistSingleAntiParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        7, "Hit ITS 5");
    fHistSingleAntiParticleNclsITSSharedTiming->GetXaxis()->SetBinLabel(
        8, "TOF timing");

    fHistSingleAntiParticleNcrossedTPCBefore =
        new TH1F("fHistSingleAntiParticleNcrossedTPCBefore",
                 ";# crossed rows TPC (before); Entries", 170, 0, 170);
    fHistSingleAntiParticleNcrossedTPCAfter =
        new TH1F("fHistSingleAntiParticleNcrossedTPCAfter",
                 ";# crossed rows TPC (after); Entries", 170, 0, 170);
    fHistSingleAntiParticleFindableTPCBefore =
        new TH1F("fHistSingleAntiParticleFindableTPCBefore",
                 ";Findable TPC (before); Entries", 1000, 0, 2);
    fHistSingleAntiParticleFindableTPCAfter =
        new TH1F("fHistSingleAntiParticleFindableTPCAfter",
                 "; Findable TPC (after); Entries) ", 1000, 0, 2);
    fHistSingleAntiParticleDCAxyBefore =
        new TH1F("fHistSingleAntiParticleDCAxyBefore",
                 "; DCA_{xy} (before); Entries", 1000, -10, 10);
    fHistSingleAntiParticleDCAxyAfterDCAz =
        new TH1F("fHistSingleAntiParticleDCAxyAfterDCAz",
                 "; DCA_{xy} (after DCA_{z} cut}); Entries", 1000, -10, 10);
    fHistSingleAntiParticleDCAxyAfter =
        new TH1F("fHistSingleAntiParticleDCAxyAfter",
                 "; DCA_{xy} (after); Entries", 1000, -10, 10);
    fHistSingleAntiParticleDCAzBefore =
        new TH1F("fHistSingleAntiParticleDCAzBefore",
                 ";DCA_{z} (before); Entries", 1000, -10, 10);
    fHistSingleAntiParticleDCAzAfter =
        new TH1F("fHistSingleAntiParticleDCAzAfter",
                 ";DCA_{z} (after); Entries", 1000, -10, 10);
    fHistSingleAntiParticleNsigmaBefore =
        new TH2F("fHistSingleAntiParticleNsigmaBefore",
                 "; p [GeV/c]; n_{#sigma} before", 1000, 0, 5, 1000, 0, 10);
    fHistSingleAntiParticleNsigmaAfter =
        new TH2F("fHistSingleAntiParticleNsigmaAfter",
                 "; p [GeV/c]; n_{#sigma} after", 1000, 0, 5, 1000, 0, 10);
    fHistSingleAntiParticleNsigmaTPCBefore =
        new TH2F("fHistSingleAntiParticleNsigmaTPCBefore",
                 "; p [GeV/c]; n_{#sigma} TPC before", 1000, 0, 5, 1000, -5, 5);
    fHistSingleAntiParticleNsigmaTPCAfter =
        new TH2F("fHistSingleAntiParticleNsigmaTPCAfter",
                 "; p [GeV/c]; n_{#sigma} TPC after", 1000, 0, 5, 1000, -5, 5);
    fHistSingleAntiParticleNsigmaTOFBefore =
        new TH2F("fHistSingleAntiParticleNsigmaTOFBefore",
                 "; p [GeV/c]; n_{#sigma} TOF before", 1000, 0, 5, 1000, -5, 5);
    fHistSingleAntiParticleNsigmaTOFAfter =
        new TH2F("fHistSingleAntiParticleNsigmaTOFAfter",
                 "; p [GeV/c]; n_{#sigma} TOF after", 1000, 0, 5, 1000, -5, 5);

    fHistSingleAntiParticleTPCsignalBefore = new TH2F(
        "fHistSingleAntiParticleTPCsignalBefore",
        "; p [GeV/c]; TPC d#it{E}/d{x} before", 500, 0, 5, 500, 0, 200);
    fHistSingleAntiParticleTOFsignalBefore =
        new TH2F("fHistSingleAntiParticleTOFsignalBefore",
                 "; p [GeV/c]; TOF #beta before", 500, 0, 5, 500, 0, 1);
    fHistSingleAntiParticleTPCsignalAfter = new TH2F(
        "fHistSingleAntiParticleTPCsignalAfter",
        "; p [GeV/c];  TPC d#it{E}/d#it{x} after", 500, 0, 5, 500, 0, 200);
    fHistSingleAntiParticleTOFsignalAfter =
        new TH2F("fHistSingleAntiParticleTOFsignalAfter",
                 "; p [GeV/c]; TOF #beta after", 500, 0, 5, 500, 0, 1);

    TH2F *singleAntiParticlePIDrejBefore[nPIDbins];
    TH2F *singleAntiParticlePIDrejAfter[nPIDbins];
    for (int i = 0; i < static_cast<int>(fParticleAdditionalRejection.size());
         ++i) {
      singleAntiParticlePIDrejBefore[i] =
          new TH2F(Form("fHistSingleAntiParticlePIDRejBefore_%i_pt_%.1f_%1.1f",
                        fParticleAdditionalRejection[i],
                        fPIDAdditionalRejectionLowerpTBound[i],
                        fPIDAdditionalRejectionUpperpTBound[i]),
                   Form("; p [GeV/c]; n_{#sigma} %i before",
                        fParticleAdditionalRejection[i]),
                   1000, 0, 5, 1000, -5, 5);
      fHistSingleAntiParticleNsigmaTPCRejBefore.push_back(
          singleAntiParticlePIDrejBefore[i]);
      singleAntiParticlePIDrejAfter[i] =
          new TH2F(Form("fHistSingleAntiParticlePIDRejAfter_%i_pt_%.1f_%1.1f",
                        fParticleAdditionalRejection[i],
                        fPIDAdditionalRejectionLowerpTBound[i],
                        fPIDAdditionalRejectionUpperpTBound[i]),
                   Form("; p [GeV/c]; n_{#sigma} %i before",
                        fParticleAdditionalRejection[i]),
                   1000, 0, 5, 1000, -5, 5);
      fHistSingleAntiParticleNsigmaTPCRejAfter.push_back(
          singleAntiParticlePIDrejAfter[i]);
    }

    fHistogramsSingleAntiParticleBefore->Add(fHistSingleAntiParticleEtaBefore);
    fHistogramsSingleAntiParticleAfter->Add(fHistSingleAntiParticleEtaAfter);
    fHistogramsSingleAntiParticleBefore->Add(fHistSingleAntiParticlePtBefore);
    fHistogramsSingleAntiParticleAfter->Add(fHistSingleAntiParticlePtAfter);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleNclsTPCBefore);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleNclsTPCAfter);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleNclsTPCShared);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleNclsTPCSharedTiming);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleNclsITSShared);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleNclsITSSharedTiming);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleNcrossedTPCBefore);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleNcrossedTPCAfter);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleFindableTPCBefore);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleFindableTPCAfter);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleDCAxyBefore);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleDCAxyAfterDCAz);
    fHistogramsSingleAntiParticleAfter->Add(fHistSingleAntiParticleDCAxyAfter);
    fHistogramsSingleAntiParticleBefore->Add(fHistSingleAntiParticleDCAzBefore);
    fHistogramsSingleAntiParticleAfter->Add(fHistSingleAntiParticleDCAzAfter);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleNsigmaBefore);
    fHistogramsSingleAntiParticleAfter->Add(fHistSingleAntiParticleNsigmaAfter);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleNsigmaTPCBefore);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleNsigmaTPCAfter);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleNsigmaTOFBefore);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleNsigmaTOFAfter);

    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleTPCsignalBefore);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleTPCsignalAfter);
    fHistogramsSingleAntiParticleBefore->Add(
        fHistSingleAntiParticleTOFsignalBefore);
    fHistogramsSingleAntiParticleAfter->Add(
        fHistSingleAntiParticleTOFsignalAfter);

    fHistogramsSingleAntiParticle->Add(fHistogramsSingleAntiParticleBefore);
    fHistogramsSingleAntiParticle->Add(fHistogramsSingleAntiParticleAfter);

    for (int i = 0;
         i < static_cast<int>(fHistSingleAntiParticleNsigmaTPCRejBefore.size());
         ++i) {
      fHistogramsSingleAntiParticleBefore->Add(
          fHistSingleAntiParticleNsigmaTPCRejBefore[i]);
      fHistogramsSingleAntiParticleAfter->Add(
          fHistSingleAntiParticleNsigmaTPCRejAfter[i]);
    }
  }
  fHistograms->Add(fHistogramsSingleParticle);
  fHistograms->Add(fHistogramsSingleAntiParticle);
}

//____________________________________________________________________________________________________
void AliSigma0SingleParticleCuts::SetAdditionalPIDRejection(
    AliPID::EParticleType part, float lower, float upper, float pTlower,
    float pTupper) {
  fPIDAdditionalRejection = true;
  fParticleAdditionalRejection.emplace_back(part);
  fPIDAdditionalRejectionLowerBound.emplace_back(lower);
  fPIDAdditionalRejectionUpperBound.emplace_back(upper);
  fPIDAdditionalRejectionLowerpTBound.emplace_back(pTlower);
  fPIDAdditionalRejectionUpperpTBound.emplace_back(pTupper);
}
