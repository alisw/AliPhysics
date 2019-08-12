#include "AliAnalysisTaskSigma0Run2.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"

ClassImp(AliAnalysisTaskSigma0Run2)

    //____________________________________________________________________________________________________
    AliAnalysisTaskSigma0Run2::AliAnalysisTaskSigma0Run2()
    : AliAnalysisTaskSE("AliAnalysisTaskSigma0Run2"),
      fAliEventCuts(),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fV0Cuts(nullptr),
      fAntiV0Cuts(nullptr),
      fPhotonQA(nullptr),
      fSigmaCuts(nullptr),
      fAntiSigmaCuts(nullptr),
      fIsMC(false),
      fIsHeavyIon(false),
      fIsLightweight(false),
      fIsRun1(false),
      fPhotonLegPileUpCut(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fMultMode(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fOutputContainer(nullptr),
      fQA(nullptr),
      fHistCutQA(nullptr),
      fHistRunNumber(nullptr),
      fHistCutBooking(nullptr),
      fHistCentralityProfileBefore(nullptr),
      fHistCentralityProfileAfter(nullptr),
      fHistCentralityProfileCoarseAfter(nullptr),
      fHistMultiplicityRef08(nullptr),
      fHistTriggerBefore(nullptr),
      fHistTriggerAfter(nullptr),
      fHistMultiplicity(nullptr) {}

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Run2::AliAnalysisTaskSigma0Run2(const char *name)
    : AliAnalysisTaskSE(name),
      fAliEventCuts(),
      fInputEvent(nullptr),
      fMCEvent(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fV0Cuts(nullptr),
      fAntiV0Cuts(nullptr),
      fPhotonQA(nullptr),
      fSigmaCuts(nullptr),
      fAntiSigmaCuts(nullptr),
      fIsMC(false),
      fIsHeavyIon(false),
      fIsLightweight(false),
      fIsRun1(false),
      fPhotonLegPileUpCut(false),
      fV0PercentileMax(100.f),
      fTrigger(AliVEvent::kINT7),
      fMultMode(AliVEvent::kINT7),
      fGammaArray(nullptr),
      fOutputContainer(nullptr),
      fQA(nullptr),
      fHistCutQA(nullptr),
      fHistRunNumber(nullptr),
      fHistCutBooking(nullptr),
      fHistCentralityProfileBefore(nullptr),
      fHistCentralityProfileAfter(nullptr),
      fHistCentralityProfileCoarseAfter(nullptr),
      fHistMultiplicityRef08(nullptr),
      fHistTriggerBefore(nullptr),
      fHistTriggerAfter(nullptr),
      fHistMultiplicity(nullptr) {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________________________________________
AliAnalysisTaskSigma0Run2::~AliAnalysisTaskSigma0Run2() {}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Run2::UserExec(Option_t * /*option*/) {
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

  // LAMBDA SELECTION
  fV0Cuts->SelectV0(fInputEvent, fMCEvent);

  // LAMBDA SELECTION
  fAntiV0Cuts->SelectV0(fInputEvent, fMCEvent);

  // PHOTON SELECTION
  fGammaArray = fV0Reader->GetReconstructedGammas();  // Gammas from default Cut
  if (!fIsLightweight) {
    fPhotonQA->PhotonQA(fInputEvent, fMCEvent, fGammaArray);
  }
  std::vector<AliSigma0ParticleV0> gammaConvContainer;
  CastToVector(gammaConvContainer, fInputEvent);

  // Sigma0 selection
  fSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, gammaConvContainer,
                                 fV0Cuts->GetV0s());

  // Sigma0 selection
  fAntiSigmaCuts->SelectPhotonMother(fInputEvent, fMCEvent, gammaConvContainer,
                                     fAntiV0Cuts->GetV0s());

  // flush the data
  PostData(1, fOutputContainer);
}

//____________________________________________________________________________________________________
bool AliAnalysisTaskSigma0Run2::AcceptEvent(AliVEvent *event) {
  if (fIsRun1) {
    return AcceptEventRun1(event);
  } else {
    return AcceptEventRun2(event);
  }
}

//____________________________________________________________________________________________________
bool AliAnalysisTaskSigma0Run2::AcceptEventRun1(AliVEvent *event) {
  if (!fIsLightweight) {
    fHistRunNumber->Fill(0.f, event->GetRunNumber());
    FillTriggerHisto(fHistTriggerBefore);
  }
  fHistCutQA->Fill(0);

  auto esdEvent = static_cast<AliESDEvent *>(event);
  // Checks if we have a primary vertex  // Get primary vertices form ESD
  bool eventVtxExist = false;
  if (esdEvent->GetPrimaryVertexTracks()->GetNContributors() > 0) {
    eventVtxExist = true;
  } else if (esdEvent->GetPrimaryVertexSPD()->GetNContributors() > 0) {
    eventVtxExist = true;
  }
  if (!eventVtxExist) {
    return false;
  }

  const auto esdVertex5 = esdEvent->GetPrimaryVertex();

  if (!(std::abs(esdVertex5->GetZ()) < 10.)) {
    return false;
  }

  // Check for pileup and fill pileup histograms
  if (esdEvent->IsPileupFromSPD()) {
    return false;
  }
  if (!fIsLightweight) {
    FillTriggerHisto(fHistTriggerAfter);
  }
  fHistCutQA->Fill(2);

  bool isConversionEventSelected =
      ((AliConvEventCuts *)fV0Reader->GetEventCuts())
          ->EventIsSelected(event, static_cast<AliMCEvent *>(fMCEvent));
  if (!isConversionEventSelected) return false;

  fHistCutQA->Fill(4);
  return true;
}

//____________________________________________________________________________________________________
bool AliAnalysisTaskSigma0Run2::AcceptEventRun2(AliVEvent *event) {
  if (!fIsLightweight) {
    fHistRunNumber->Fill(0.f, event->GetRunNumber());
    FillTriggerHisto(fHistTriggerBefore);
  }
  fHistCutQA->Fill(0);

  // Remove overlap with HM trigger from MB
  if (!fIsMC && fTrigger == AliVEvent::kINT7 &&
      (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0))
    return false;
  fHistCutQA->Fill(1);

  // EVENT SELECTION
  if (!fAliEventCuts.AcceptEvent(event)) return false;
  if (!fIsLightweight) {
    FillTriggerHisto(fHistTriggerAfter);
  }
  fHistCutQA->Fill(2);

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
    fHistCutQA->Fill(3);
  }

  bool isConversionEventSelected =
      ((AliConvEventCuts *)fV0Reader->GetEventCuts())
          ->EventIsSelected(event, static_cast<AliMCEvent *>(fMCEvent));
  if (!isConversionEventSelected) return false;

  if (!fIsLightweight) {
    fHistCentralityProfileCoarseAfter->Fill(lPercentile);
    fHistMultiplicityRef08->Fill(AliESDtrackCuts::GetReferenceMultiplicity(
        static_cast<AliESDEvent *>(event), AliESDtrackCuts::kTrackletsITSTPC,
        0.8, 0));
  }

  fHistMultiplicity->Fill(
      AliSigma0PhotonMotherCuts::GetMultiplicityBin(lPercentile, fMultMode));

  fHistCutQA->Fill(4);
  return true;
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Run2::CastToVector(
    std::vector<AliSigma0ParticleV0> &container, const AliVEvent *inputEvent) {
  for (int iGamma = 0; iGamma < fGammaArray->GetEntriesFast(); ++iGamma) {
    auto *PhotonCandidate =
        dynamic_cast<AliAODConversionPhoton *>(fGammaArray->At(iGamma));
    if (!PhotonCandidate) continue;

    auto pos =
        (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel1());
    auto neg =
        (AliESDtrack *)inputEvent->GetTrack(PhotonCandidate->GetLabel2());
    if (!pos || !neg) continue;

    // pile up check
    if (fPhotonLegPileUpCut) {
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

    AliSigma0ParticleV0 phot(PhotonCandidate, pos, neg, inputEvent);
    if (fIsMC) {
      phot.MatchToMC(fMCEvent, 22, {{11, -11}});
    }
    container.push_back(phot);
  }
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Run2::UserCreateOutputObjects() {
  if (fOutputContainer != nullptr) {
    delete fOutputContainer;
    fOutputContainer = nullptr;
  }
  if (fOutputContainer == nullptr) {
    fOutputContainer = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fQA = new TList();
  fQA->SetName("EventCuts");
  fQA->SetOwner(true);

  if (fTrigger != AliVEvent::kINT7) {
    fAliEventCuts.OverrideAutomaticTriggerSelection(fTrigger);
  }

  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (!fV0Reader) {
    AliError("No V0 reader");
    return;
  }

  if (!fIsLightweight) {
    if (fV0Reader->GetEventCuts() &&
        fV0Reader->GetEventCuts()->GetCutHistograms()) {
      fOutputContainer->Add(fV0Reader->GetEventCuts()->GetCutHistograms());
    }
    if (fV0Reader->GetConversionCuts() &&
        fV0Reader->GetConversionCuts()->GetCutHistograms()) {
      fOutputContainer->Add(fV0Reader->GetConversionCuts()->GetCutHistograms());
    }
  }

  fHistCutQA = new TH1F("fHistCutQA", ";;Entries", 10, 0, 10);
  fHistCutQA->GetXaxis()->SetBinLabel(1, "Event");
  fHistCutQA->GetXaxis()->SetBinLabel(2, "Overlap with HM");
  fHistCutQA->GetXaxis()->SetBinLabel(3, "AliEventCuts");
  fHistCutQA->GetXaxis()->SetBinLabel(4, "Multiplicity selection");
  fHistCutQA->GetXaxis()->SetBinLabel(5, "AliConversionCuts");
  fQA->Add(fHistCutQA);

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

    fHistMultiplicityRef08 =
        new TH1F("fHistMultiplicityRef08", "; Multiplicity Ref08; Entries",
                 1000, 0, 1000);
    fQA->Add(fHistMultiplicityRef08);

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

  std::vector<float> multBinsLow, multBinsUp;
  if (fMultMode == AliVEvent::kINT7) {
    multBinsLow = {{0, 5., 15., 30., 50.}};
    multBinsUp = {{5., 15., 30., 50., 100.}};
  } else if (fMultMode == AliVEvent::kHighMultV0) {
    multBinsLow = {{0, 0.01, 0.05, 0.1, 1}};
    multBinsUp = {{0.01, 0.05, 0.1, 1, 100.}};
  }
  fHistMultiplicity =
      new TH1I("fHistMultiplicity", "; Multiplicity bin; Entries", 5, 0, 5);
  fHistMultiplicity->GetXaxis()->LabelsOption("u");
  for (int i = 0; i < static_cast<int>(multBinsLow.size()); i++) {
    fHistMultiplicity->GetXaxis()->SetBinLabel(
        i + 1, Form("V0M: %.2f - %.2f %%", multBinsLow[i], multBinsUp[i]));
  }
  fQA->Add(fHistMultiplicity);

  fOutputContainer->Add(fQA);

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

  if (!fIsLightweight) {
    fPhotonQA = new AliSigma0V0Cuts();
    fPhotonQA->SetIsMC(fIsMC);
    fPhotonQA->SetLightweight(false);
    fPhotonQA->SetPID(22);
    fPhotonQA->SetPosPID(AliPID::kElectron, 11);
    fPhotonQA->SetNegPID(AliPID::kElectron, -11);
    fPhotonQA->InitCutHistograms(TString("Photon"));
    fOutputContainer->Add(fPhotonQA->GetCutHistograms());
  }

  if (fSigmaCuts && fSigmaCuts->GetCutHistograms()) {
    fOutputContainer->Add(fSigmaCuts->GetCutHistograms());
  }

  if (fAntiSigmaCuts && fAntiSigmaCuts->GetCutHistograms()) {
    fOutputContainer->Add(fAntiSigmaCuts->GetCutHistograms());
  }

  PostData(1, fOutputContainer);
}

//____________________________________________________________________________________________________
void AliAnalysisTaskSigma0Run2::FillTriggerHisto(TH1F *histo) {
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
