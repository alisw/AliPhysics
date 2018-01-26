#include "AliSigma0EventCuts.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODTracklets.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"

ClassImp(AliSigma0EventCuts)

//____________________________________________________________________________________________________
AliSigma0EventCuts::AliSigma0EventCuts()
    : TObject(),
      fAliEventCuts(),
      fHistograms(nullptr),
      fQA(nullptr),
      fV0Reader(nullptr),
      fV0ReaderName("NoInit"),
      fUseAliEventCuts(true),
      fUseConversionCuts(false),
      fVertex(999.f),
      fNVertexContributors(0),
      fV0PercentileMax(100.f),
      fTrigger(0),
      fHistCuts(nullptr),
      fHistEventQA(nullptr),
      fHistTrigger(nullptr),
      fHistXVertex(nullptr),
      fHistYVertex(nullptr),
      fHistZVertex(nullptr),
      fHistXVertexAfter(nullptr),
      fHistYVertexAfter(nullptr),
      fHistZVertexAfter(nullptr),
      fHistVertexSeparation(nullptr),
      fHistSPDresolution(nullptr),
      fHistSPDTracklets(nullptr),
      fHistSPDTrackletsAfter(nullptr),
      fHistRunNumber(nullptr),
      fHistCentralityProfile(nullptr),
      fHistCentralityProfileAfter(nullptr),
      fHistCentralityProfileCoarseAfter(nullptr),
      fAnaUtils(new AliAnalysisUtils()),
      fInputHandler(nullptr) {
  //  fAnaUtils->SetMinPlpContribSPD(3);
  // configure analysis tools
  fAnaUtils->SetUseMVPlpSelection(false);
  fAnaUtils->SetMinPlpContribMV(
      5);  // probably not needed (see above) but better be safe
  fAnaUtils->SetMinPlpContribSPD(3);
  fAnaUtils->SetMaxVtxZ(fVertex);  // later configured by calling SetZVertexCut
}

//____________________________________________________________________________________________________
AliSigma0EventCuts::AliSigma0EventCuts(const AliSigma0EventCuts &ref)
    : TObject(ref),
      fAliEventCuts(),
      fHistograms(nullptr),
      fQA(nullptr),
      fV0Reader(ref.fV0Reader),
      fV0ReaderName(ref.fV0ReaderName.Data()),
      fUseAliEventCuts(ref.fUseAliEventCuts),
      fUseConversionCuts(ref.fUseConversionCuts),
      fVertex(ref.fVertex),
      fNVertexContributors(ref.fNVertexContributors),
      fV0PercentileMax(ref.fV0PercentileMax),
      fTrigger(ref.fTrigger),
      fHistCuts(nullptr),
      fHistEventQA(nullptr),
      fHistTrigger(nullptr),
      fHistXVertex(nullptr),
      fHistYVertex(nullptr),
      fHistZVertex(nullptr),
      fHistXVertexAfter(nullptr),
      fHistYVertexAfter(nullptr),
      fHistZVertexAfter(nullptr),
      fHistVertexSeparation(nullptr),
      fHistSPDresolution(nullptr),
      fHistSPDTracklets(nullptr),
      fHistSPDTrackletsAfter(nullptr),
      fHistRunNumber(nullptr),
      fHistCentralityProfile(nullptr),
      fHistCentralityProfileAfter(nullptr),
      fHistCentralityProfileCoarseAfter(nullptr),
      fAnaUtils(new AliAnalysisUtils()),
      fInputHandler(nullptr) {
  fAnaUtils->SetMinPlpContribSPD(3);
}

//____________________________________________________________________________________________________
AliSigma0EventCuts &AliSigma0EventCuts::operator=(
    const AliSigma0EventCuts &ref) {
  // Assignment operator
  if (this == &ref) return *this;

  fVertex = ref.fVertex;
  return (*this);
}

//____________________________________________________________________________________________________
AliSigma0EventCuts::~AliSigma0EventCuts() {
  if (fAnaUtils) delete fAnaUtils;
}

//____________________________________________________________________________________________________
bool AliSigma0EventCuts::EventIsSelected(AliVEvent *fInputEvent,
                                         AliVEvent *fMCEvent) {
  fInputHandler =
      (AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()
                                   ->GetInputEventHandler());
  if (fHistograms) fHistEventQA->Fill(0);
  if (fUseAliEventCuts)
    return EventIsSelectedAliEventCuts(fInputEvent, fMCEvent);
  else
    return EventIsSelectedCuts(fInputEvent, fMCEvent);
}

//____________________________________________________________________________________________________
bool AliSigma0EventCuts::EventIsSelectedAliEventCuts(AliVEvent *fInputEvent,
                                                     AliVEvent *fMCEvent) {
  fHistRunNumber->Fill(0.f, fInputEvent->GetRunNumber());

  FillTriggerHisto(fHistTrigger);
  bool isSelected = fAliEventCuts.AcceptEvent(fInputEvent);

  if (!fAliEventCuts.PassedCut(AliEventCuts::kTrigger)) return false;
  fHistEventQA->Fill(1);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kDAQincomplete)) return false;
  fHistEventQA->Fill(2);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kBfield)) return false;
  fHistEventQA->Fill(3);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kVertexSPD)) return false;
  fHistEventQA->Fill(4);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kVertexTracks)) return false;
  fHistEventQA->Fill(5);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kVertex)) return false;
  fHistEventQA->Fill(6);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kVertexPositionSPD)) return false;
  fHistEventQA->Fill(7);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kVertexPositionTracks))
    return false;
  fHistEventQA->Fill(8);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kVertexPosition)) return false;
  fHistEventQA->Fill(9);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kVertexQuality)) return false;
  fHistEventQA->Fill(10);
  fHistSPDTracklets->Fill(
      fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),
      fInputEvent->GetNumberOfITSClusters(1));
  if (!fAliEventCuts.PassedCut(AliEventCuts::kPileUp)) return false;
  fHistEventQA->Fill(11);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kMultiplicity)) return false;
  fHistEventQA->Fill(12);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kINELgt0)) return false;
  fHistEventQA->Fill(13);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kCorrelations)) return false;
  fHistEventQA->Fill(14);
  if (!fAliEventCuts.PassedCut(AliEventCuts::kAllCuts)) return false;
  fHistEventQA->Fill(15);
  fHistSPDTrackletsAfter->Fill(
      fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),
      fInputEvent->GetNumberOfITSClusters(1));

  Float_t lPercentile = 300;
  AliMultSelection *MultSelection = 0x0;
  MultSelection =
      (AliMultSelection *)fInputEvent->FindListObject("MultSelection");
  if (!MultSelection) {
    // If you get this warning (and lPercentiles 300) please check that the
    // AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  } else {
    lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }
  fHistCentralityProfile->Fill(lPercentile);
  if (lPercentile > fV0PercentileMax) return false;
  fHistCentralityProfileCoarseAfter->Fill(lPercentile);
  fHistCentralityProfileAfter->Fill(lPercentile);
  fHistEventQA->Fill(16);

  // Now check whether photons have been reconstructed for that event
  bool isConversionEventSelected = false;
  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (fV0Reader)
    isConversionEventSelected =
        ((AliConvEventCuts *)fV0Reader->GetEventCuts())
            ->EventIsSelected(fInputEvent, static_cast<AliMCEvent *>(fMCEvent));
  if (fUseConversionCuts && !isConversionEventSelected) return false;
  if (isConversionEventSelected) fHistEventQA->Fill(17);
  return isSelected;
}

//____________________________________________________________________________________________________
bool AliSigma0EventCuts::EventIsSelectedCuts(AliVEvent *fInputEvent,
                                             AliVEvent *fMCEvent) {
  // Trigger selection
  fHistRunNumber->Fill(0.f, fInputEvent->GetRunNumber());

  FillTriggerHisto(fHistTrigger);
  fHistEventQA->Fill(1);

  // Event has magnetic field
  if (std::fabs(fInputEvent->GetMagneticField()) < 0.001) return false;
  fHistEventQA->Fill(2);

  // Event has track and SPD vertex
  const AliVVertex *vertex = fInputEvent->GetPrimaryVertex();
  if (!vertex) return false;
  fHistEventQA->Fill(3);

  // Ncontrib
  if (vertex->GetNContributors() < fNVertexContributors) return false;
  fHistEventQA->Fill(4);

  const float zVertex = vertex->GetZ();
  if (fHistograms) {
    fHistXVertex->Fill(vertex->GetX());
    fHistYVertex->Fill(vertex->GetY());
    fHistZVertex->Fill(zVertex);
  }

  // z-vertex cut
  if (std::abs(zVertex) > fVertex) return false;
  if (fHistograms) fHistEventQA->Fill(5);

  if (fHistograms)
    fHistSPDTracklets->Fill(
        fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),
        fInputEvent->GetNumberOfITSClusters(1));
  if (fAnaUtils->IsPileUpEvent(fInputEvent)) return false;
  fHistEventQA->Fill(6);
  AliAODEvent *event = static_cast<AliAODEvent *>(fInputEvent);
  if (event->IsPileupFromSPD()) return false;
  fHistEventQA->Fill(7);
  if (fAnaUtils->IsSPDClusterVsTrackletBG(fInputEvent)) return false;
  fHistEventQA->Fill(8);

  if (fHistograms) {
    fHistXVertexAfter->Fill(vertex->GetX());
    fHistYVertexAfter->Fill(vertex->GetY());
    fHistZVertexAfter->Fill(zVertex);
    fHistSPDTrackletsAfter->Fill(
        fInputEvent->GetMultiplicity()->GetNumberOfTracklets(),
        fInputEvent->GetNumberOfITSClusters(1));
  }

  // Now check whether photons have been reconstructed for that event
  bool isConversionEventSelected = false;
  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (fV0Reader)
    isConversionEventSelected =
        ((AliConvEventCuts *)fV0Reader->GetEventCuts())
            ->EventIsSelected(fInputEvent, static_cast<AliMCEvent *>(fMCEvent));
  if (fUseConversionCuts && !isConversionEventSelected) return false;
  if (isConversionEventSelected) fHistEventQA->Fill(9);

  return true;
}

//____________________________________________________________________________________________________
void AliSigma0EventCuts::FillTriggerHisto(TH1F *histo) {
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

//____________________________________________________________________________________________________
void AliSigma0EventCuts::InitCutHistograms() {
  std::cout << "============================\n"
            << " EVENT CUT CONFIGURATION \n"
            << " AliEventCuts  " << fUseAliEventCuts << "\n"
            << " ConversionCut " << fUseConversionCuts << "\n"
            << " z-vertex      " << fVertex << "\n"
            << " Vertex contr  " << fNVertexContributors << "\n"
            << " V0 percentile " << fV0PercentileMax << "\n"
            << " Trigger       " << fTrigger << "\n"
            << "============================\n";

  TH1::AddDirectory(kFALSE);

  if (fHistograms != nullptr) {
    delete fHistograms;
    fHistograms = nullptr;
  }
  if (fHistograms == nullptr) {
    fHistograms = new TList();
    fHistograms->SetOwner(kTRUE);
    fHistograms->SetName("EventCut_QA");
  }

  fQA = new TList();
  fQA->SetName("AliEventCuts");
  fQA->SetOwner(true);

  if (fUseAliEventCuts) {
    fAliEventCuts.SetManualMode();
    fAliEventCuts.fTriggerMask = fTrigger;
    fAliEventCuts.AddQAplotsToList(fQA);
    fHistograms->Add(fQA);

    fHistEventQA = new TH1F("fHistEventQA", ";; Entries", 18, 0, 18);
    fHistEventQA->GetXaxis()->SetBinLabel(1, "Event");
    fHistEventQA->GetXaxis()->SetBinLabel(2, "Trigger selection");
    fHistEventQA->GetXaxis()->SetBinLabel(3, "Incomplete DAQ");
    fHistEventQA->GetXaxis()->SetBinLabel(4, "B-field");
    fHistEventQA->GetXaxis()->SetBinLabel(5, "SPD vertex");
    fHistEventQA->GetXaxis()->SetBinLabel(6, "Track vertex");
    fHistEventQA->GetXaxis()->SetBinLabel(7, "Vertex");
    fHistEventQA->GetXaxis()->SetBinLabel(8, "SPD vertex position");
    fHistEventQA->GetXaxis()->SetBinLabel(9, "Track vertex position");
    fHistEventQA->GetXaxis()->SetBinLabel(10, "Vertex position");
    fHistEventQA->GetXaxis()->SetBinLabel(11, "Vertex quality");
    fHistEventQA->GetXaxis()->SetBinLabel(12, "Pile-up");
    fHistEventQA->GetXaxis()->SetBinLabel(13, "Multiplicity");
    fHistEventQA->GetXaxis()->SetBinLabel(14, "INEL");
    fHistEventQA->GetXaxis()->SetBinLabel(15, "Correlations");
    fHistEventQA->GetXaxis()->SetBinLabel(16, "All");
    fHistEventQA->GetXaxis()->SetBinLabel(17, "V0 percentile");
    fHistEventQA->GetXaxis()->SetBinLabel(18, "AliConvEventCuts");
    fHistograms->Add(fHistEventQA);

    fHistCuts = new TProfile("fHistCuts", ";;Cut value", 2, 0, 2);
    fHistCuts->GetXaxis()->SetBinLabel(1, "V0 percentile");
    fHistCuts->GetXaxis()->SetBinLabel(2, "Conversion cuts");
    fHistograms->Add(fHistCuts);

    fHistCuts->Fill(0.f, fV0PercentileMax);
    fHistCuts->Fill(1.f, fUseConversionCuts);

    fHistCentralityProfile = new TH1F(
        "fHistCentralityProfile", "; V0 percentile [%]; Entries", 1000, 0, 5);
    fHistCentralityProfileAfter =
        new TH1F("fHistCentralityProfileAfter", "; V0 percentile [%]; Entries",
                 1000, 0, 5);
    fHistCentralityProfileCoarseAfter =
        new TH1F("fHistCentralityProfileCoarseAfter",
                 "; V0 percentile [%]; Entries", 100, 0, 100);
    fHistograms->Add(fHistCentralityProfile);
    fHistograms->Add(fHistCentralityProfileAfter);
    fHistograms->Add(fHistCentralityProfileCoarseAfter);
  } else {
    fHistCuts = new TProfile("fHistCuts", ";;Cut value", 4, 0, 4);
    fHistCuts->GetXaxis()->SetBinLabel(1, "N vertex contributors min");
    fHistCuts->GetXaxis()->SetBinLabel(2, "z-vertex");
    fHistCuts->GetXaxis()->SetBinLabel(3, "V0 percentile");
    fHistCuts->GetXaxis()->SetBinLabel(4, "Conversion cuts");
    fHistograms->Add(fHistCuts);

    fHistCuts->Fill(0.f, fNVertexContributors);
    fHistCuts->Fill(1.f, fVertex);
    fHistCuts->Fill(2.f, fV0PercentileMax);
    fHistCuts->Fill(3.f, fUseConversionCuts);

    fHistEventQA = new TH1F("fHistEventQA", ";; Entries", 10, 0, 10);
    fHistEventQA->GetXaxis()->SetBinLabel(1, "Event");
    fHistEventQA->GetXaxis()->SetBinLabel(2, "Trigger selection");
    fHistEventQA->GetXaxis()->SetBinLabel(3, "B-field");
    fHistEventQA->GetXaxis()->SetBinLabel(4, "Vertex");
    fHistEventQA->GetXaxis()->SetBinLabel(5, "Ncontrib vertex");
    fHistEventQA->GetXaxis()->SetBinLabel(6, "z-vertex");
    fHistEventQA->GetXaxis()->SetBinLabel(7, "AnaUtils::IsPileUpEvent");
    fHistEventQA->GetXaxis()->SetBinLabel(8, "AnaUtils::IsPileUpSPD");
    fHistEventQA->GetXaxis()->SetBinLabel(9, "AnaUtils::SPDclusters");
    fHistEventQA->GetXaxis()->SetBinLabel(10, "AliConvEventCuts");
    fHistograms->Add(fHistEventQA);

    fHistXVertex =
        new TH1F("fHistXVertex", ";x vertex [cm]; Entries", 1000, -15, 15);
    fHistYVertex =
        new TH1F("fHistYVertex", ";y vertex [cm]; Entries", 1000, -15, 15);
    fHistZVertex =
        new TH1F("fHistZVertex", ";z vertex [cm]; Entries", 1000, -15, 15);

    fHistograms->Add(fHistXVertex);
    fHistograms->Add(fHistYVertex);
    fHistograms->Add(fHistZVertex);

    fHistXVertexAfter =
        new TH1F("fHistXVertexAfter", ";x vertex [cm]; Entries", 1000, -15, 15);
    fHistYVertexAfter =
        new TH1F("fHistYVertexAfter", ";y vertex [cm]; Entries", 1000, -15, 15);
    fHistZVertexAfter =
        new TH1F("fHistZVertexAfter", ";z vertex [cm]; Entries", 1000, -15, 15);
    fHistograms->Add(fHistXVertexAfter);
    fHistograms->Add(fHistYVertexAfter);
    fHistograms->Add(fHistZVertexAfter);
  }

  // V0Reader QA
  fV0Reader =
      (AliV0ReaderV1 *)AliAnalysisManager::GetAnalysisManager()->GetTask(
          fV0ReaderName.Data());
  if (fV0Reader) {
    if ((AliConvEventCuts *)fV0Reader->GetEventCuts()) {
      if (((AliConvEventCuts *)fV0Reader->GetEventCuts())->GetCutHistograms()) {
        fHistograms->Add(((AliConvEventCuts *)fV0Reader->GetEventCuts())
                             ->GetCutHistograms());
      }
    }
  }

  fHistSPDTracklets =
      new TH2I("fHistSPDTracklets", "; SPD tracklets; SPD cluster", 200, 0, 200,
               200, 0, 200);
  fHistograms->Add(fHistSPDTracklets);

  fHistSPDTrackletsAfter =
      new TH2I("fHistSPDTrackletsAfter", "; SPD tracklets; SPD cluster", 200, 0,
               200, 200, 0, 200);
  fHistograms->Add(fHistSPDTrackletsAfter);

  fHistRunNumber = new TProfile("fHistRunNumber", ";;Run Number", 1, 0, 1);
  fHistograms->Add(fHistRunNumber);

  fHistTrigger = new TH1F("fHistTrigger", ";;Entries", 50, 0, 50);
  fHistTrigger->GetXaxis()->LabelsOption("u");
  fHistTrigger->GetXaxis()->SetBinLabel(1, "kMB");
  fHistTrigger->GetXaxis()->SetBinLabel(2, "kINT1");
  fHistTrigger->GetXaxis()->SetBinLabel(3, "kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(4, "kMUON");
  fHistTrigger->GetXaxis()->SetBinLabel(5, "kHighMult");
  fHistTrigger->GetXaxis()->SetBinLabel(6, "kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(7, "kEMC1");
  fHistTrigger->GetXaxis()->SetBinLabel(8, "kCINT5");
  fHistTrigger->GetXaxis()->SetBinLabel(9, "kINT5");
  fHistTrigger->GetXaxis()->SetBinLabel(10, "kCMUS5");
  fHistTrigger->GetXaxis()->SetBinLabel(11, "kMUSPB");
  fHistTrigger->GetXaxis()->SetBinLabel(12, "kINT7inMUON");
  fHistTrigger->GetXaxis()->SetBinLabel(13, "kMuonSingleHighPt7");
  fHistTrigger->GetXaxis()->SetBinLabel(14, "kMUSH7");
  fHistTrigger->GetXaxis()->SetBinLabel(15, "kMUSHPB");
  fHistTrigger->GetXaxis()->SetBinLabel(16, "kMuonLikeLowPt7");
  fHistTrigger->GetXaxis()->SetBinLabel(17, "kMUL7");
  fHistTrigger->GetXaxis()->SetBinLabel(18, "kMuonLikePB");
  fHistTrigger->GetXaxis()->SetBinLabel(19, "kMuonUnlikeLowPt7");
  fHistTrigger->GetXaxis()->SetBinLabel(20, "kMUU7");
  fHistTrigger->GetXaxis()->SetBinLabel(21, "kMuonUnlikePB");
  fHistTrigger->GetXaxis()->SetBinLabel(22, "kEMC7");
  fHistTrigger->GetXaxis()->SetBinLabel(23, "kEMC8");
  fHistTrigger->GetXaxis()->SetBinLabel(24, "kMUS7");
  fHistTrigger->GetXaxis()->SetBinLabel(25, "kMuonSingleLowPt7");
  fHistTrigger->GetXaxis()->SetBinLabel(26, "kPHI1");
  fHistTrigger->GetXaxis()->SetBinLabel(27, "kPHI7");
  fHistTrigger->GetXaxis()->SetBinLabel(28, "kPHI8");
  fHistTrigger->GetXaxis()->SetBinLabel(29, "kPHOSPb");
  fHistTrigger->GetXaxis()->SetBinLabel(30, "kEMCEJE");
  fHistTrigger->GetXaxis()->SetBinLabel(31, "kEMCEGA");
  fHistTrigger->GetXaxis()->SetBinLabel(32, "kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(33, "kCentral");
  fHistTrigger->GetXaxis()->SetBinLabel(34, "kSemiCentral");
  fHistTrigger->GetXaxis()->SetBinLabel(35, "kDG");
  fHistTrigger->GetXaxis()->SetBinLabel(36, "kDG5");
  fHistTrigger->GetXaxis()->SetBinLabel(37, "kZED");
  fHistTrigger->GetXaxis()->SetBinLabel(38, "kSPI7");
  fHistTrigger->GetXaxis()->SetBinLabel(39, "kSPI");
  fHistTrigger->GetXaxis()->SetBinLabel(40, "kINT8");
  fHistTrigger->GetXaxis()->SetBinLabel(41, "kMuonSingleLowPt8");
  fHistTrigger->GetXaxis()->SetBinLabel(42, "kMuonSingleHighPt8");
  fHistTrigger->GetXaxis()->SetBinLabel(43, "kMuonLikeLowPt8");
  fHistTrigger->GetXaxis()->SetBinLabel(44, "kMuonUnlikeLowPt8");
  fHistTrigger->GetXaxis()->SetBinLabel(45, "kMuonUnlikeLowPt0");
  fHistTrigger->GetXaxis()->SetBinLabel(46, "kUserDefined");
  fHistTrigger->GetXaxis()->SetBinLabel(47, "kTRD");
  fHistTrigger->GetXaxis()->SetBinLabel(48, "kFastOnly");
  fHistTrigger->GetXaxis()->SetBinLabel(49, "kAny");
  fHistTrigger->GetXaxis()->SetBinLabel(50, "kAnyINT");
  fHistograms->Add(fHistTrigger);
}
