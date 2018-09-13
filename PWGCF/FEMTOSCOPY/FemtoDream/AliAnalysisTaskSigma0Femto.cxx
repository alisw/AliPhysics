#include "AliAnalysisTaskSigma0Femto.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"

ClassImp(AliAnalysisTaskSigma0Femto)

    //____________________________________________________________________________________________________
    AliAnalysisTaskSigma0Femto::AliAnalysisTaskSigma0Femto()
    : AliAnalysisTaskSE("AliAnalysisTaskSigma0Femto"),
      fAliEventCuts(),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fV0Cuts(nullptr),
      fAntiV0Cuts(nullptr),
      fSigmaCuts(nullptr),
      fAntiSigmaCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(false),
      fIsHeavyIon(false),
      fIsLightweight(false),
      fPhotonLegPileUpCut(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fOutputContainer(nullptr),
      fQA(nullptr),
      fOutputFemto(nullptr),
      fHistCutQA(nullptr),
      fHistRunNumber(nullptr),
      fHistCutBooking(nullptr),
      fHistCentralityProfileBefore(nullptr),
      fHistCentralityProfileAfter(nullptr),
      fHistCentralityProfileCoarseAfter(nullptr),
      fHistTriggerBefore(nullptr),
      fHistTriggerAfter(nullptr),
      fHistPhotonPileUp(nullptr) {}

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Femto::AliAnalysisTaskSigma0Femto(const char *name)
    : AliAnalysisTaskSE(name),
      fAliEventCuts(),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fV0Cuts(nullptr),
      fAntiV0Cuts(nullptr),
      fSigmaCuts(nullptr),
      fAntiSigmaCuts(nullptr),
      fProtonTrack(nullptr),
      fTrackCutsPartProton(nullptr),
      fTrackCutsPartAntiProton(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fIsMC(false),
      fIsHeavyIon(false),
      fIsLightweight(false),
      fPhotonLegPileUpCut(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fOutputContainer(nullptr),
      fQA(nullptr),
      fOutputFemto(nullptr),
      fHistCutQA(nullptr),
      fHistRunNumber(nullptr),
      fHistCutBooking(nullptr),
      fHistCentralityProfileBefore(nullptr),
      fHistCentralityProfileAfter(nullptr),
      fHistCentralityProfileCoarseAfter(nullptr),
      fHistTriggerBefore(nullptr),
      fHistTriggerAfter(nullptr),
      fHistPhotonPileUp(nullptr) {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Femto::~AliAnalysisTaskSigma0Femto() {
  delete fPartColl;
  delete fPairCleaner;
  delete fProtonTrack;
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::UserExec(Option_t * /*option*/) {
  AliVEvent *fInputEvent = InputEvent();
  if (fIsMC) fMCEvent = MCEvent();

  // PREAMBLE - CHECK EVERYTHING IS THERE
  if (!fInputEvent) {
    AliError("No Input event");
    return;
  }

  if (fIsMC && !fMCEvent) {
    AliError("No MC event");
    return;
  }

  if (!fTrackCutsPartProton || !fTrackCutsPartAntiProton) {
    AliError("Proton Cuts missing");
    return;
  }

  if (!fV0Cuts || !fAntiV0Cuts) {
    AliError("V0 Cuts missing");
    return;
  }

  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (!fV0Reader) {
    AliError("No V0 reader");
    return;
  }

  if (!fSigmaCuts || !fAntiSigmaCuts) {
    AliError("Sigma0 Cuts missing");
    return;
  }

  // EVENT SELECTION
  if (!AcceptEvent(fInputEvent)) return;

  // PROTON SELECTION
  static std::vector<AliFemtoDreamBasePart> particles;
  particles.clear();
  static std::vector<AliFemtoDreamBasePart> antiParticles;
  antiParticles.clear();
  const AliESDEvent *evt = static_cast<AliESDEvent *>(fInputEvent);
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliESDtrack *track = static_cast<AliESDtrack *>(evt->GetTrack(iTrack));
    fProtonTrack->SetTrack(track, fMCEvent,
                           AliESDtrackCuts::GetReferenceMultiplicity(
                               evt, AliESDtrackCuts::kTrackletsITSTPC, 0.8, 0));
    if (fTrackCutsPartProton->isSelected(fProtonTrack)) {
      particles.push_back(*fProtonTrack);
    }
    if (fTrackCutsPartAntiProton->isSelected(fProtonTrack)) {
      antiParticles.push_back(*fProtonTrack);
    }
  }

  // LAMBDA SELECTION
  fV0Cuts->SelectV0(fInputEvent, fMCEvent);

  // LAMBDA SELECTION
  fAntiV0Cuts->SelectV0(fInputEvent, fMCEvent);

  // PHOTON SELECTION
  fGammaArray = fV0Reader->GetReconstructedGammas();  // Gammas from default Cut
  std::vector<AliSigma0ParticleV0> gammaConvContainer;
  CastToVector(gammaConvContainer, fInputEvent);

  // Sigma0 selection
  fSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, gammaConvContainer,
                                 fV0Cuts->GetV0s());

  // Sigma0 selection
  fAntiSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, gammaConvContainer,
                                     fAntiV0Cuts->GetV0s());

  // Convert the Sigma0 into Femto particles
  static std::vector<AliFemtoDreamBasePart> sigma0particles;
  sigma0particles.clear();
  static std::vector<AliFemtoDreamBasePart> antiSigma0particles;
  antiSigma0particles.clear();
  static std::vector<AliFemtoDreamBasePart> sigma0sidebandUp;
  sigma0sidebandUp.clear();
  static std::vector<AliFemtoDreamBasePart> antiSigma0sidebandUp;
  antiSigma0sidebandUp.clear();
  static std::vector<AliFemtoDreamBasePart> sigma0sidebandLow;
  sigma0sidebandLow.clear();
  static std::vector<AliFemtoDreamBasePart> antiSigma0sidebandLow;
  antiSigma0sidebandLow.clear();

  CastToVector(fSigmaCuts->GetSigma(), sigma0particles);
  CastToVector(fAntiSigmaCuts->GetSigma(), antiSigma0particles);
  CastToVector(fSigmaCuts->GetSidebandUp(), sigma0sidebandUp);
  CastToVector(fAntiSigmaCuts->GetSidebandUp(), antiSigma0sidebandUp);
  CastToVector(fSigmaCuts->GetSidebandDown(), sigma0sidebandLow);
  CastToVector(fAntiSigmaCuts->GetSidebandDown(), antiSigma0sidebandLow);

  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0particles, 0);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0particles, 1);
  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0sidebandUp, 2);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0sidebandUp, 3);
  fPairCleaner->CleanTrackAndDecay(&particles, &sigma0sidebandLow, 4);
  fPairCleaner->CleanTrackAndDecay(&antiParticles, &antiSigma0sidebandLow, 5);

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(particles);
  fPairCleaner->StoreParticle(antiParticles);
  fPairCleaner->StoreParticle(sigma0particles);
  fPairCleaner->StoreParticle(antiSigma0particles);
  fPairCleaner->StoreParticle(sigma0sidebandUp);
  fPairCleaner->StoreParticle(antiSigma0sidebandUp);
  fPairCleaner->StoreParticle(sigma0sidebandLow);
  fPairCleaner->StoreParticle(antiSigma0sidebandLow);

  AliMultSelection *MultSelection = 0x0;
  MultSelection =
      (AliMultSelection *)fInputEvent->FindListObject("MultSelection");
  if (MultSelection) {
    fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                        fInputEvent->GetPrimaryVertex()->GetZ(),
                        AliESDtrackCuts::GetReferenceMultiplicity(
                            evt, AliESDtrackCuts::kTrackletsITSTPC, 0.8, 0),
                        MultSelection->GetMultiplicityPercentile("V0M"));
  }

  // flush the data
  PostData(1, fOutputContainer);
  PostData(2, fOutputFemto);
}

//____________________________________________________________________________________________________
bool AliAnalysisTaskSigma0Femto::AcceptEvent(AliVEvent *event) {
  if (!fIsLightweight) {
    fHistRunNumber->Fill(0.f, event->GetRunNumber());
    FillTriggerHisto(fHistTriggerBefore);
  }
  fHistCutQA->Fill(0);

  // EVENT SELECTION
  if (!fAliEventCuts.AcceptEvent(event)) return false;
  if (!fIsLightweight) {
    FillTriggerHisto(fHistTriggerAfter);
  }
  fHistCutQA->Fill(1);

  Float_t lPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection *)event->FindListObject("MultSelection");
  if (!MultSelection) {
    // If you get this warning (and lPercentiles 300) please check that the
    // AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  } else {
    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }
  if (!fIsLightweight) fHistCentralityProfileBefore->Fill(lPercentile);

  // MULTIPLICITY SELECTION
  if (fTrigger == AliVEvent::kHighMultV0 && fV0PercentileMax < 100.f) {
    if (lPercentile > fV0PercentileMax) return false;
    if (!fIsLightweight) {
      fHistCentralityProfileAfter->Fill(lPercentile);
    }
    fHistCutQA->Fill(2);
  }

  bool isConversionEventSelected =
      ((AliConvEventCuts *)fV0Reader->GetEventCuts())
          ->EventIsSelected(event, static_cast<AliMCEvent *>(fMCEvent));
  if (!isConversionEventSelected) return false;

  if (!fIsLightweight) fHistCentralityProfileCoarseAfter->Fill(lPercentile);

  fHistCutQA->Fill(3);
  return true;
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::CastToVector(
    std::vector<AliSigma0ParticleV0> &container, const AliVEvent *inputEvent) {
  for (int iGamma = 0; iGamma < fGammaArray->GetEntriesFast(); ++iGamma) {
    auto *PhotonCandidate =
        dynamic_cast<AliAODConversionPhoton *>(fGammaArray->At(iGamma));
    if (!PhotonCandidate) continue;
    fHistPhotonPileUp->Fill(PhotonCandidate->Pt(), 0.5);

    // pile up check
    if (fPhotonLegPileUpCut) {
      auto pos =
          (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel1());
      auto neg =
          (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel2());
      if (!pos || !neg) continue;

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

      if (!posTrackCombined || !negTrackCombined) continue;
    }
    fHistPhotonPileUp->Fill(PhotonCandidate->Pt(), 1.5);

    AliSigma0ParticleV0 phot(PhotonCandidate, inputEvent);
    if (fIsMC) {
      const int label = phot.MatchToMC(fMCEvent, 22, {{11, -11}});
    }
    container.push_back(phot);
  }
}
//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::CastToVector(
    std::vector<AliSigma0ParticlePhotonMother> &sigmaContainer,
    std::vector<AliFemtoDreamBasePart> &particles) {
  for (const auto &sigma : sigmaContainer) {
    auto sigmaFemto = AliFemtoDreamBasePart(sigma);
    particles.push_back(sigmaFemto);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::UserCreateOutputObjects() {
  if (fOutputContainer != nullptr) {
    delete fOutputContainer;
    fOutputContainer = nullptr;
  }
  if (fOutputContainer == nullptr) {
    fOutputContainer = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  if (fOutputFemto != nullptr) {
    delete fOutputFemto;
    fOutputFemto = nullptr;
  }
  if (fOutputFemto == nullptr) {
    fOutputFemto = new TList();
    fOutputFemto->SetOwner(kTRUE);
  }

  fQA = new TList();
  fQA->SetName("EventCuts");
  fQA->SetOwner(true);

  if (fTrigger != AliVEvent::kINT7) {
    fAliEventCuts.SetManualMode();
    if (!fIsHeavyIon) fAliEventCuts.SetupRun2pp();
    fAliEventCuts.fTriggerMask = fTrigger;
  }

  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (!fV0Reader) {
    AliError("No V0 reader");
    return;
  }

  if (fV0Reader->GetEventCuts() &&
      fV0Reader->GetEventCuts()->GetCutHistograms()) {
    fOutputContainer->Add(fV0Reader->GetEventCuts()->GetCutHistograms());
  }
  if (fV0Reader->GetConversionCuts() &&
      fV0Reader->GetConversionCuts()->GetCutHistograms()) {
    fOutputContainer->Add(fV0Reader->GetConversionCuts()->GetCutHistograms());
  }
  if (fV0Reader->GetProduceV0FindingEfficiency() &&
      fV0Reader->GetV0FindingEfficiencyHistograms()) {
    fOutputContainer->Add(fV0Reader->GetV0FindingEfficiencyHistograms());
  }
  if (fV0Reader->GetProduceImpactParamHistograms()) {
    fOutputContainer->Add(fV0Reader->GetImpactParamHistograms());
  }

  fHistCutQA = new TH1F("fHistCutQA", ";;Entries", 5, 0, 5);
  fHistCutQA->GetXaxis()->SetBinLabel(1, "Event");
  fHistCutQA->GetXaxis()->SetBinLabel(2, "AliEventCuts");
  fHistCutQA->GetXaxis()->SetBinLabel(3, "Multiplicity selection");
  fHistCutQA->GetXaxis()->SetBinLabel(4, "AliConversionCuts");
  fQA->Add(fHistCutQA);

  fHistPhotonPileUp =
      new TH2F("fHistPhotonPileUp", ";#it{p}_{T} (GeV/#it{c}^{2}; PileUp", 100,
               0, 10, 2, 0, 2);
  fHistPhotonPileUp->GetYaxis()->SetBinLabel(1, "Before");
  fHistPhotonPileUp->GetYaxis()->SetBinLabel(2, "After");
  fQA->Add(fHistPhotonPileUp);

  if (!fIsLightweight) {
    fHistRunNumber = new TProfile("fHistRunNumber", ";;Run Number", 1, 0, 1);
    fQA->Add(fHistRunNumber);

    fHistCutBooking = new TProfile("fHistCutBooking", ";;Cut value", 1, 0, 1);
    fHistCutBooking->GetXaxis()->SetBinLabel(1, "V0 percentile");
    fQA->Add(fHistCutBooking);
    fHistCutBooking->Fill(0.f, fV0PercentileMax);

    fAliEventCuts.AddQAplotsToList(fQA);

    fHistCentralityProfileBefore =
        new TH1F("fHistCentralityProfileBefore", "; V0 percentile (%); Entries",
                 1000, 0, 5);
    fHistCentralityProfileAfter =
        new TH1F("fHistCentralityProfileAfter", "; V0 percentile (%); Entries",
                 1000, 0, 5);
    fHistCentralityProfileCoarseAfter =
        new TH1F("fHistCentralityProfileCoarseAfter",
                 "; V0 percentile (%); Entries", 100, 0, 100);
    fQA->Add(fHistCentralityProfileBefore);
    fQA->Add(fHistCentralityProfileAfter);
    fQA->Add(fHistCentralityProfileCoarseAfter);

    fHistTriggerBefore = new TH1F("fHistTriggerBefore", ";;Entries", 50, 0, 50);
    fHistTriggerBefore->GetXaxis()->LabelsOption("u");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(1, "kMB");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(2, "kINT1");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(3, "kINT7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(4, "kMUON");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(5, "kHighMult");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(6, "kHighMultSPD");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(7, "kEMC1");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(8, "kCINT5");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(9, "kINT5");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(10, "kCMUS5");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(11, "kMUSPB");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(12, "kINT7inMUON");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(13, "kMuonSingleHighPt7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(14, "kMUSH7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(15, "kMUSHPB");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(16, "kMuonLikeLowPt7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(17, "kMUL7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(18, "kMuonLikePB");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(19, "kMuonUnlikeLowPt7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(20, "kMUU7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(21, "kMuonUnlikePB");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(22, "kEMC7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(23, "kEMC8");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(24, "kMUS7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(25, "kMuonSingleLowPt7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(26, "kPHI1");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(27, "kPHI7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(28, "kPHI8");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(29, "kPHOSPb");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(30, "kEMCEJE");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(31, "kEMCEGA");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(32, "kHighMultV0");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(33, "kCentral");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(34, "kSemiCentral");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(35, "kDG");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(36, "kDG5");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(37, "kZED");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(38, "kSPI7");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(39, "kSPI");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(40, "kINT8");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(41, "kMuonSingleLowPt8");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(42, "kMuonSingleHighPt8");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(43, "kMuonLikeLowPt8");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(44, "kMuonUnlikeLowPt8");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(45, "kMuonUnlikeLowPt0");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(46, "kUserDefined");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(47, "kTRD");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(48, "kFastOnly");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(49, "kAny");
    fHistTriggerBefore->GetXaxis()->SetBinLabel(50, "kAnyINT");
    fQA->Add(fHistTriggerBefore);

    fHistTriggerAfter = new TH1F("fHistTriggerAfter", ";;Entries", 50, 0, 50);
    fHistTriggerAfter->GetXaxis()->LabelsOption("u");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(1, "kMB");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(2, "kINT1");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(3, "kINT7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(4, "kMUON");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(5, "kHighMult");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(6, "kHighMultSPD");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(7, "kEMC1");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(8, "kCINT5");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(9, "kINT5");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(10, "kCMUS5");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(11, "kMUSPB");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(12, "kINT7inMUON");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(13, "kMuonSingleHighPt7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(14, "kMUSH7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(15, "kMUSHPB");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(16, "kMuonLikeLowPt7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(17, "kMUL7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(18, "kMuonLikePB");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(19, "kMuonUnlikeLowPt7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(20, "kMUU7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(21, "kMuonUnlikePB");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(22, "kEMC7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(23, "kEMC8");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(24, "kMUS7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(25, "kMuonSingleLowPt7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(26, "kPHI1");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(27, "kPHI7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(28, "kPHI8");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(29, "kPHOSPb");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(30, "kEMCEJE");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(31, "kEMCEGA");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(32, "kHighMultV0");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(33, "kCentral");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(34, "kSemiCentral");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(35, "kDG");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(36, "kDG5");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(37, "kZED");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(38, "kSPI7");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(39, "kSPI");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(40, "kINT8");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(41, "kMuonSingleLowPt8");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(42, "kMuonSingleHighPt8");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(43, "kMuonLikeLowPt8");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(44, "kMuonUnlikeLowPt8");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(45, "kMuonUnlikeLowPt0");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(46, "kUserDefined");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(47, "kTRD");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(48, "kFastOnly");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(49, "kAny");
    fHistTriggerAfter->GetXaxis()->SetBinLabel(50, "kAnyINT");
    fQA->Add(fHistTriggerAfter);
  }

  fOutputContainer->Add(fQA);

  fProtonTrack = new AliFemtoDreamTrack();
  fProtonTrack->SetUseMCInfo(fIsMC);

  fTrackCutsPartProton->Init("Proton");
  // necessary for the non-min booking case
  fTrackCutsPartProton->SetName("Proton");

  if (fTrackCutsPartProton && fTrackCutsPartProton->GetQAHists()) {
    fOutputFemto->Add(fTrackCutsPartProton->GetQAHists());
    if (fIsMC && !fTrackCutsPartProton->GetMinimalBooking() &&
        fTrackCutsPartProton->GetMCQAHists()) {
      fTrackCutsPartProton->SetMCName("MC_Proton");
      fOutputFemto->Add(fTrackCutsPartProton->GetMCQAHists());
    }
  }

  fTrackCutsPartAntiProton->Init("Anti-proton");
  // necessary for the non-min booking case
  fTrackCutsPartAntiProton->SetName("Anti-proton");

  if (fTrackCutsPartAntiProton && fTrackCutsPartAntiProton->GetQAHists()) {
    fOutputFemto->Add(fTrackCutsPartAntiProton->GetQAHists());
    if (fIsMC && !fTrackCutsPartAntiProton->GetMinimalBooking() &&
        fTrackCutsPartAntiProton->GetMCQAHists()) {
      fTrackCutsPartAntiProton->SetMCName("MC_Anti-proton");
      fOutputFemto->Add(fTrackCutsPartAntiProton->GetMCQAHists());
    }
  }

  if (fV0Cuts) fV0Cuts->InitCutHistograms(TString("Lambda"));
  if (fAntiV0Cuts) fAntiV0Cuts->InitCutHistograms(TString("AntiLambda"));
  if (fSigmaCuts) fSigmaCuts->InitCutHistograms(TString("Sigma0"));
  if (fAntiSigmaCuts) fAntiSigmaCuts->InitCutHistograms(TString("AntiSigma0"));

  if (fV0Cuts && fV0Cuts->GetCutHistograms()) {
    fOutputContainer->Add(fV0Cuts->GetCutHistograms());
  }

  if (fAntiV0Cuts && fAntiV0Cuts->GetCutHistograms()) {
    fOutputContainer->Add(fAntiV0Cuts->GetCutHistograms());
  }

  if (fSigmaCuts && fSigmaCuts->GetCutHistograms()) {
    fOutputContainer->Add(fSigmaCuts->GetCutHistograms());
  }

  if (fAntiSigmaCuts && fAntiSigmaCuts->GetCutHistograms()) {
    fOutputContainer->Add(fAntiSigmaCuts->GetCutHistograms());
  }

  fPairCleaner =
      new AliFemtoDreamPairCleaner(6, 0, fConfig->GetMinimalBookingME());
  fPartColl =
      new AliFemtoDreamPartCollection(fConfig, fConfig->GetMinimalBookingME());

  if (!fConfig->GetMinimalBookingME() && fPairCleaner &&
      fPairCleaner->GetHistList() && !fIsMC) {
    fOutputFemto->Add(fPairCleaner->GetHistList());
  }

  if (fPartColl && fPartColl->GetHistList()) {
    fOutputFemto->Add(fPartColl->GetHistList());
  }
  if (!fConfig->GetMinimalBookingME() && fPartColl && fPartColl->GetQAList()) {
    fOutputFemto->Add(fPartColl->GetQAList());
  }

  PostData(1, fOutputContainer);
  PostData(2, fOutputFemto);
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Femto::FillTriggerHisto(TH1F *histo) {
  if (fInputHandler->IsEventSelected() & AliVEvent::kMB) histo->Fill(0);
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT1) histo->Fill(1);
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) histo->Fill(2);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMUON) histo->Fill(3);
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMult) histo->Fill(4);
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD)
    histo->Fill(5);
  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC1) histo->Fill(6);
  if (fInputHandler->IsEventSelected() & AliVEvent::kCINT5) histo->Fill(7);
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT5) histo->Fill(8);
  if (fInputHandler->IsEventSelected() & AliVEvent::kCMUS5) histo->Fill(9);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMUSPB) histo->Fill(10);
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT7inMUON)
    histo->Fill(11);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleHighPt7)
    histo->Fill(12);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMUSH7) histo->Fill(13);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMUSHPB) histo->Fill(14);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikeLowPt7)
    histo->Fill(15);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMUL7) histo->Fill(16);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikePB)
    histo->Fill(17);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt7)
    histo->Fill(18);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMUU7) histo->Fill(19);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikePB)
    histo->Fill(20);
  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC7) histo->Fill(21);
  if (fInputHandler->IsEventSelected() & AliVEvent::kEMC8) histo->Fill(22);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMUS7) histo->Fill(23);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleLowPt7)
    histo->Fill(24);
  if (fInputHandler->IsEventSelected() & AliVEvent::kPHI1) histo->Fill(25);
  if (fInputHandler->IsEventSelected() & AliVEvent::kPHI7) histo->Fill(26);
  if (fInputHandler->IsEventSelected() & AliVEvent::kPHI8) histo->Fill(27);
  if (fInputHandler->IsEventSelected() & AliVEvent::kPHOSPb) histo->Fill(28);
  if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE) histo->Fill(29);
  if (fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA) histo->Fill(30);
  if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)
    histo->Fill(31);
  if (fInputHandler->IsEventSelected() & AliVEvent::kCentral) histo->Fill(32);
  if (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)
    histo->Fill(33);
  if (fInputHandler->IsEventSelected() & AliVEvent::kDG) histo->Fill(34);
  if (fInputHandler->IsEventSelected() & AliVEvent::kDG5) histo->Fill(35);
  if (fInputHandler->IsEventSelected() & AliVEvent::kZED) histo->Fill(36);
  if (fInputHandler->IsEventSelected() & AliVEvent::kSPI7) histo->Fill(37);
  if (fInputHandler->IsEventSelected() & AliVEvent::kSPI) histo->Fill(38);
  if (fInputHandler->IsEventSelected() & AliVEvent::kINT8) histo->Fill(39);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleLowPt8)
    histo->Fill(40);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonSingleHighPt8)
    histo->Fill(41);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonLikeLowPt8)
    histo->Fill(42);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt8)
    histo->Fill(43);
  if (fInputHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt0)
    histo->Fill(44);
  if (fInputHandler->IsEventSelected() & AliVEvent::kUserDefined)
    histo->Fill(45);
  if (fInputHandler->IsEventSelected() & AliVEvent::kTRD) histo->Fill(46);
  if (fInputHandler->IsEventSelected() & AliVEvent::kFastOnly) histo->Fill(47);
  if (fInputHandler->IsEventSelected() & AliVEvent::kAny) histo->Fill(48);
  if (fInputHandler->IsEventSelected() & AliVEvent::kAnyINT) histo->Fill(49);
}
