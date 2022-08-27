#include "TCanvas.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TRefArray.h"
#include "TTree.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliEventPoolManager.h"
#include "AliEventPoolMuon.h"
#include "AliGenEventHeader.h"

#include "AliVEvent.h"
#include "AliVHeader.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliAODDimuon.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliAODv0.h"
#include "AliExternalTrackParam.h"
#include "AliPIDCombined.h"
#include "AliTOFHeader.h"

#include "AliAnalysisMuonUtility.h"
#include "AliMultSelection.h"

#include "AliEventCuts.h"
#include "AliPIDResponse.h"

#include "AliAnalysisTaskAODTrackPair.h"
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

#include "AliAnalysisTaskAODTrackPairUtils.h"

#include "iostream"
#include "memory"
// Authors: Satoshi Yano
// Reviewed:

using namespace std;

ClassImp(AliAnalysisTaskAODTrackPair)

    AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair()
    : AliAnalysisTaskSE(), fEvent(NULL), fPoolMuonTrackMgr(NULL), fUtils(NULL),
      fIsMC(false), fIsMidTrackAna(false), fIsV0TrackPairAna(false),
      fIsPrimTrackPairAna(false), fIsMixingAnalysis(false), fRunNumber(-99999),
      fTrackDepth(1000), fPoolSize(1), fReadyFraction(0.1),
      fTriggerMaskForSame(AliVEvent::kMuonUnlikeLowPt7 |
                          AliVEvent::kMuonLikeLowPt7),
      fTriggerMaskForMixing(AliVEvent::kMuonSingleLowPt7),
      onEvtMixingPoolVtxZ(true), onEvtMixingPoolCent(true),
      onEvtMixingPoolPsi(true), fIsCINT7(false), fIsCMSL7(false),
      fIsCMSH7(false), fIsCMUL7(false), fIsCMLL7(false),

      fOutputList(NULL), fEventCounter(NULL),

      fHistTrackPairPtBalance(NULL), fHistTrackPairLocalBoardPair(NULL),

      fHistTrackThetaAbs(NULL), fHistTrackTriggerMatch(NULL),
      fHistTrackPDCA(NULL), fHistTrackChiSquare(NULL),
      fHistTriggerChiSquare(NULL),

      fHistEventVtxZ(NULL), fHistEventCent(NULL), fHistEventMulti(NULL),
      fHistEventVtxCont(NULL),

      fTreeULSPair(NULL), fTreeLSppPair(NULL), fTreeLSmmPair(NULL),

      fSparseULSPairMassPt(NULL), fSparseLSppPairMassPt(NULL),
      fSparseLSmmPairMassPt(NULL),

      fTreeULSPair_TightCut(NULL), fTreeLSppPair_TightCut(NULL),
      fTreeLSmmPair_TightCut(NULL),

      fTreeULSPair_ProngV0(NULL), fTreeLSppPair_ProngV0(NULL),
      fTreeLSmmPair_ProngV0(NULL),

      fTreeMixULSPair(NULL), fTreeMixLSppPair(NULL), fTreeMixLSmmPair(NULL),

      fHistULSPairMassPt(NULL), fHistLSppPairMassPt(NULL),
      fHistLSmmPairMassPt(NULL),

      fHistULSPairMassPt_TightCut(NULL), fHistLSppPairMassPt_TightCut(NULL),
      fHistLSmmPairMassPt_TightCut(NULL),

      fHistULSPairMassPt_ProngV0(NULL), fHistLSppPairMassPt_ProngV0(NULL),
      fHistLSmmPairMassPt_ProngV0(NULL),

      fHistMixULSPairMassPt(NULL), fHistMixLSppPairMassPt(NULL),
      fHistMixLSmmPairMassPt(NULL),

      fSparseMixULSPairMassPt(NULL), fSparseMixLSppPairMassPt(NULL),
      fSparseMixLSmmPairMassPt(NULL),

      fHistMassK0s1K0s2(NULL),

      fHistTPCdEdxP(NULL), fHistBetaP(NULL), fHistTPCSigmaElectron(NULL),
      fHistTOFSigmaElectron(NULL), fHistTPCSigmaMuon(NULL),
      fHistTOFSigmaMuon(NULL), fHistTPCSigmaPion(NULL), fHistTOFSigmaPion(NULL),
      fHistTPCSigmaKaon(NULL), fHistTOFSigmaKaon(NULL),
      fHistTPCSigmaProton(NULL), fHistTOFSigmaProton(NULL),

      fHistSelTPCdEdxP(NULL), fHistSelBetaP(NULL),
      fHistSelTPCSigmaElectron(NULL), fHistSelTOFSigmaElectron(NULL),
      fHistSelTPCSigmaMuon(NULL), fHistSelTOFSigmaMuon(NULL),
      fHistSelTPCSigmaPion(NULL), fHistSelTOFSigmaPion(NULL),
      fHistSelTPCSigmaKaon(NULL), fHistSelTOFSigmaKaon(NULL),
      fHistSelTPCSigmaProton(NULL), fHistSelTOFSigmaProton(NULL),

      fHistTrackP(NULL), fHistTrackPt(NULL), fHistTrackEta(NULL),
      fHistTPCNClusts(NULL), fHistSPDNClusts(NULL),
      fHistTPCCrossRowsFindableRatio(NULL), fHistReducedChi2TPC(NULL),
      fHistReducedChi2ITS(NULL), fHistDCAz(NULL), fHistDCAxyPt(NULL),

      fHistArmenteros(NULL), fHistV0MassDecayLength(NULL),
      fHistV0MassPointingAngle(NULL), fHistV0MassV0DCA(NULL),
      fHistV0MassV0TrackDCA(NULL), fHistV0MassV0DecayRadius(NULL),
      fHistV0MassV0PropLifeTime(NULL), fHistSelArmenteros(NULL),
      fHistSelV0MassDecayLength(NULL), fHistSelV0MassPointingAngle(NULL),
      fHistSelV0MassV0DCA(NULL), fHistSelV0MassV0TrackDCA(NULL),
      fHistSelV0MassV0DecayRadius(NULL), fHistSelV0MassV0PropLifeTime(NULL),

      RecPairPt(0.), RecPairRap(0.), RecPairMass(0.),
      RecPairArmenterosArmPt(0.), RecPairArmenterosAlpha(0.), RecPairCent(0.),
      RecPairDS(0.) {}

AliAnalysisTaskAODTrackPair::AliAnalysisTaskAODTrackPair(const char *name)
    : AliAnalysisTaskSE(name), fEvent(NULL), fPoolMuonTrackMgr(NULL),
      fUtils(NULL), fIsMC(false), fIsMidTrackAna(false),
      fIsV0TrackPairAna(false), fIsPrimTrackPairAna(false),
      fIsMixingAnalysis(false), fRunNumber(-99999), fTrackDepth(1000),
      fPoolSize(1), fReadyFraction(0.1),
      fTriggerMaskForSame(AliVEvent::kMuonUnlikeLowPt7 |
                          AliVEvent::kMuonLikeLowPt7),
      fTriggerMaskForMixing(AliVEvent::kMuonSingleLowPt7),
      onEvtMixingPoolVtxZ(true), onEvtMixingPoolCent(true),
      onEvtMixingPoolPsi(true), fIsCINT7(false), fIsCMSL7(false),
      fIsCMSH7(false), fIsCMUL7(false), fIsCMLL7(false),

      fOutputList(NULL), fEventCounter(NULL),

      fHistTrackPairPtBalance(NULL), fHistTrackPairLocalBoardPair(NULL),

      fHistTrackThetaAbs(NULL), fHistTrackTriggerMatch(NULL),
      fHistTrackPDCA(NULL), fHistTrackChiSquare(NULL),
      fHistTriggerChiSquare(NULL),

      fHistEventVtxZ(NULL), fHistEventCent(NULL), fHistEventMulti(NULL),
      fHistEventVtxCont(NULL),

      fTreeULSPair(NULL), fTreeLSppPair(NULL), fTreeLSmmPair(NULL),

      fSparseULSPairMassPt(NULL), fSparseLSppPairMassPt(NULL),
      fSparseLSmmPairMassPt(NULL),

      fTreeULSPair_TightCut(NULL), fTreeLSppPair_TightCut(NULL),
      fTreeLSmmPair_TightCut(NULL),

      fTreeULSPair_ProngV0(NULL), fTreeLSppPair_ProngV0(NULL),
      fTreeLSmmPair_ProngV0(NULL),

      fTreeMixULSPair(NULL), fTreeMixLSppPair(NULL), fTreeMixLSmmPair(NULL),

      fHistULSPairMassPt(NULL), fHistLSppPairMassPt(NULL),
      fHistLSmmPairMassPt(NULL),

      fHistULSPairMassPt_TightCut(NULL), fHistLSppPairMassPt_TightCut(NULL),
      fHistLSmmPairMassPt_TightCut(NULL),

      fHistULSPairMassPt_ProngV0(NULL), fHistLSppPairMassPt_ProngV0(NULL),
      fHistLSmmPairMassPt_ProngV0(NULL),

      fHistMixULSPairMassPt(NULL), fHistMixLSppPairMassPt(NULL),
      fHistMixLSmmPairMassPt(NULL),

      fSparseMixULSPairMassPt(NULL), fSparseMixLSppPairMassPt(NULL),
      fSparseMixLSmmPairMassPt(NULL),

      fHistMassK0s1K0s2(NULL),

      fHistTPCdEdxP(NULL), fHistBetaP(NULL), fHistTPCSigmaElectron(NULL),
      fHistTOFSigmaElectron(NULL), fHistTPCSigmaMuon(NULL),
      fHistTOFSigmaMuon(NULL), fHistTPCSigmaPion(NULL), fHistTOFSigmaPion(NULL),
      fHistTPCSigmaKaon(NULL), fHistTOFSigmaKaon(NULL),
      fHistTPCSigmaProton(NULL), fHistTOFSigmaProton(NULL),

      fHistSelTPCdEdxP(NULL), fHistSelBetaP(NULL),
      fHistSelTPCSigmaElectron(NULL), fHistSelTOFSigmaElectron(NULL),
      fHistSelTPCSigmaMuon(NULL), fHistSelTOFSigmaMuon(NULL),
      fHistSelTPCSigmaPion(NULL), fHistSelTOFSigmaPion(NULL),
      fHistSelTPCSigmaKaon(NULL), fHistSelTOFSigmaKaon(NULL),
      fHistSelTPCSigmaProton(NULL), fHistSelTOFSigmaProton(NULL),

      fHistTrackP(NULL), fHistTrackPt(NULL), fHistTrackEta(NULL),
      fHistTPCNClusts(NULL), fHistSPDNClusts(NULL),
      fHistTPCCrossRowsFindableRatio(NULL), fHistReducedChi2TPC(NULL),
      fHistReducedChi2ITS(NULL), fHistDCAz(NULL), fHistDCAxyPt(NULL),

      fHistArmenteros(NULL), fHistV0MassDecayLength(NULL),
      fHistV0MassPointingAngle(NULL), fHistV0MassV0DCA(NULL),
      fHistV0MassV0TrackDCA(NULL), fHistV0MassV0DecayRadius(NULL),
      fHistV0MassV0PropLifeTime(NULL), fHistSelArmenteros(NULL),
      fHistSelV0MassDecayLength(NULL), fHistSelV0MassPointingAngle(NULL),
      fHistSelV0MassV0DCA(NULL), fHistSelV0MassV0TrackDCA(NULL),
      fHistSelV0MassV0DecayRadius(NULL), fHistSelV0MassV0PropLifeTime(NULL),

      RecPairPt(0.), RecPairRap(0.), RecPairMass(0.),
      RecPairArmenterosArmPt(0.), RecPairArmenterosAlpha(0.), RecPairCent(0.),
      RecPairDS(0.) {

  double fCentBins[] = {-1, 9, 15, 21, 26, 34, 42, 51, 61, 99999};
  double fVtxBins[] = {-50, -10.5, -6, -2, 0, 2, 6, 10.5, 50};
  double fPsiBins[] = {-10, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 10};

  int fNCentBins = sizeof(fCentBins) / sizeof(double) - 1;
  int fNVtxZBins = sizeof(fVtxBins) / sizeof(double) - 1;
  int fNPsiBins = sizeof(fPsiBins) / sizeof(double) - 1;

  fPoolMuonTrackMgr = new AliEventPoolManager(
      fPoolSize, fTrackDepth, fNCentBins, (double *)fCentBins, fNVtxZBins,
      (double *)fVtxBins, fNPsiBins, (double *)fPsiBins);
  fPoolMuonTrackMgr->SetTargetValues(fTrackDepth, (double)fReadyFraction,
                                     fPoolSize);

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskAODTrackPair::~AliAnalysisTaskAODTrackPair() {}
//________________________________________________________________________
void AliAnalysisTaskAODTrackPair::UserCreateOutputObjects() {
  // Create histograms
  // Called once
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  double bins_event_hist[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int binnum_event_hist = sizeof(bins_event_hist) / sizeof(double) - 1;

  std::string event_label[] = {"CMUL7",
                               "CMLL7",
                               "CMUL7orCMLL7",
                               "CMUL7andCMLL7",
                               "CMUL7withDS",
                               "CMLL7withDS",
                               "CMUL7orCMLL7withDS",
                               "CMUL7andCMLL7withDS"};
  fEventCounter = new TH2F("fEventCounter", "", 11, 0, 11, 200, 0, 200);
  for (unsigned int iname = 0;
       iname < sizeof(event_label) / sizeof(std::string); ++iname) {
    fEventCounter->GetXaxis()->SetBinLabel(iname + 1,
                                           event_label[iname].c_str());
  }

  fOutputList->Add(fEventCounter);

  float min_mass = 0.0;
  float max_mass = 3.0;
  float width_mass = 0.01;

  float min_pt = 0.0;
  float max_pt = 20.0;
  float width_pt = 0.5;

  if (!fIsMidTrackAna) {
    if (!fIsMixingAnalysis) {
      fTreeULSPair = new TTree("fTreeULSPair", "");
      fTreeULSPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeULSPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeULSPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeULSPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeULSPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeULSPair);
      fTreeLSppPair = new TTree("fTreeLSppPair", "");
      fTreeLSppPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeLSppPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeLSppPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeLSppPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeLSppPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeLSppPair);
      fTreeLSmmPair = new TTree("fTreeLSmmPair", "");
      fTreeLSmmPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeLSmmPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeLSmmPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeLSmmPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeLSmmPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeLSmmPair);
    } else {
      fTreeMixULSPair = new TTree("fTreeMixULSPair", "");
      fTreeMixULSPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeMixULSPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeMixULSPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeMixULSPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeMixULSPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeMixULSPair);
      fTreeMixLSppPair = new TTree("fTreeMixLSppPair", "");
      fTreeMixLSppPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeMixLSppPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeMixLSppPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeMixLSppPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeMixLSppPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeMixLSppPair);
      fTreeMixLSmmPair = new TTree("fTreeMixLSmmPair", "");
      fTreeMixLSmmPair->Branch("RecPairPt", &RecPairPt, "RecPairPt/F");
      fTreeMixLSmmPair->Branch("RecPairRap", &RecPairRap, "RecPairRap/F");
      fTreeMixLSmmPair->Branch("RecPairMass", &RecPairMass, "RecPairMass/F");
      fTreeMixLSmmPair->Branch("RecPairCent", &RecPairCent, "RecPairCent/F");
      fTreeMixLSmmPair->Branch("RecPairDS", &RecPairDS, "RecPairDS/F");
      fOutputList->Add(fTreeMixLSmmPair);
    }
    fHistTrackPairPtBalance =
        new TH2F("fHistTrackPairPtBalance", "", 50, 0, 5, 50, 0, 5);
    fHistTrackPairLocalBoardPair =
        new TH2F("fHistTrackPairLocalBoardPair", "", 240, 0, 240, 240, 0, 240);
    fOutputList->Add(fHistTrackPairPtBalance);
    fOutputList->Add(fHistTrackPairLocalBoardPair);

    fHistTrackThetaAbs =
        new TH2F("fHistTrackThetaAbs", "", 20, 0, 10, 60, 0, 15);
    fHistTrackTriggerMatch =
        new TH2F("fHistTrackTriggerMatch", "", 20, 0, 10, 5, 0, 5);
    fHistTrackPDCA = new TH2F("fHistTrackPDCA", "", 20, 0, 10, 200, 0, 20);
    fHistTrackChiSquare =
        new TH2F("fHistTrackChiSquare", "", 20, 0, 10, 100, 0, 10);
    fHistTriggerChiSquare =
        new TH2F("fHistTriggerChiSquare", "", 20, 0, 10, 100, 0, 10);
    fOutputList->Add(fHistTrackThetaAbs);
    fOutputList->Add(fHistTrackTriggerMatch);
    fOutputList->Add(fHistTrackPDCA);
    fOutputList->Add(fHistTrackChiSquare);
    fOutputList->Add(fHistTriggerChiSquare);
  }

  if (fIsMidTrackAna) {
    if (!fIsMixingAnalysis) {
      fHistULSPairMassPt =
          new TH2F("fHistULSPairMassPt", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistULSPairMassPt);
      fHistLSppPairMassPt =
          new TH2F("fHistLSppPairMassPt", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistLSppPairMassPt);
      fHistLSmmPairMassPt =
          new TH2F("fHistLSmmPairMassPt", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistLSmmPairMassPt);
      int bins[3] = {int((max_mass - min_mass) / width_mass),
                     int((max_pt - min_pt) / width_pt), 100};
      double min_bins[3] = {min_mass, min_pt, 0};
      double max_bins[3] = {max_mass, max_pt, 100};
      fSparseULSPairMassPt = new THnSparseF("fSparseULSPairMassPt", "", 3, bins,
                                            min_bins, max_bins);
      fSparseLSppPairMassPt = new THnSparseF("fSparseLSppPairMassPt", "", 3,
                                             bins, min_bins, max_bins);
      fSparseLSmmPairMassPt = new THnSparseF("fSparseLSmmPairMassPt", "", 3,
                                             bins, min_bins, max_bins);
      fOutputList->Add(fSparseULSPairMassPt);
      fOutputList->Add(fSparseLSppPairMassPt);
      fOutputList->Add(fSparseLSmmPairMassPt);
    }
    if (fIsMixingAnalysis) {
      fHistMixULSPairMassPt =
          new TH2F("fHistMixULSPairMassPt", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistMixULSPairMassPt);
      int bins[3] = {int((max_mass - min_mass) / width_mass),
                     int((max_pt - min_pt) / width_pt), 100};
      double min_bins[3] = {min_mass, min_pt, 0};
      double max_bins[3] = {max_mass, max_pt, 100};
      fSparseMixULSPairMassPt = new THnSparseF("fSparseMixULSPairMassPt", "", 3,
                                               bins, min_bins, max_bins);
      fOutputList->Add(fSparseMixULSPairMassPt);
    }

    if (fIsV0TrackPairAna) {
      min_mass = 0.0;
      max_mass = 3.0;
      width_mass = 0.001;

      fHistULSPairMassPt_ProngV0 =
          new TH2F("fHistULSPairMassPt_ProngV0", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistULSPairMassPt_ProngV0);

      fHistMassK0s1K0s2 =
          new TH2F("fHistMassK0s1K0s2", "", 200, 0.4, 0.6, 200, 0.4, 0.6);
      fOutputList->Add(fHistMassK0s1K0s2);
      fHistArmenteros =
          new TH2F("fHistArmenteros", "", 200, -3, 3, 400, 0, 0.4);
      fOutputList->Add(fHistArmenteros);
      fHistSelArmenteros =
          new TH2F("fHistSelArmenteros", "", 200, -3, 3, 400, 0, 0.4);
      fOutputList->Add(fHistSelArmenteros);

      fHistV0MassDecayLength = new TH2F("fHistV0MassDecayLength", "",
                                        int((max_mass - min_mass) / width_mass),
                                        min_mass, max_mass, 1000, 0., 200.);
      fHistV0MassPointingAngle =
          new TH2F("fHistV0MassPointingAngle", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   1000, 0.9, 1.0);
      fHistV0MassV0DCA = new TH2F("fHistV0MassV0DCA", "",
                                  int((max_mass - min_mass) / width_mass),
                                  min_mass, max_mass, 1000, 0., 200.);
      fHistV0MassV0TrackDCA = new TH2F("fHistV0MassV0TrackDCA", "",
                                       int((max_mass - min_mass) / width_mass),
                                       min_mass, max_mass, 1000, 0, 100);
      fHistV0MassV0DecayRadius =
          new TH2F("fHistV0MassV0DecayRadius", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   1000, 0., 200.);
      fHistV0MassV0PropLifeTime =
          new TH2F("fHistV0MassV0PropLifeTime", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   1000, 0, 200.);
      fOutputList->Add(fHistV0MassDecayLength);
      fOutputList->Add(fHistV0MassPointingAngle);
      fOutputList->Add(fHistV0MassV0DCA);
      fOutputList->Add(fHistV0MassV0TrackDCA);
      fOutputList->Add(fHistV0MassV0DecayRadius);
      fOutputList->Add(fHistV0MassV0PropLifeTime);

      fHistSelV0MassDecayLength =
          new TH2F("fHistSelV0MassDecayLength", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   1000, 0., 200.);
      fHistSelV0MassPointingAngle =
          new TH2F("fHistSelV0MassPointingAngle", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   1000, 0.9, 1.0);
      fHistSelV0MassV0DCA = new TH2F("fHistSelV0MassV0DCA", "",
                                     int((max_mass - min_mass) / width_mass),
                                     min_mass, max_mass, 1000, 0., 200.);
      fHistSelV0MassV0TrackDCA =
          new TH2F("fHistSelV0MassV0TrackDCA", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   1000, 0, 100);
      fHistSelV0MassV0DecayRadius =
          new TH2F("fHistSelV0MassV0DecayRadius", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   1000, 0., 200.);
      fHistSelV0MassV0PropLifeTime =
          new TH2F("fHistSelV0MassV0PropLifeTime", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   1000, 0, 200.);
      fOutputList->Add(fHistSelV0MassDecayLength);
      fOutputList->Add(fHistSelV0MassPointingAngle);
      fOutputList->Add(fHistSelV0MassV0DCA);
      fOutputList->Add(fHistSelV0MassV0TrackDCA);
      fOutputList->Add(fHistSelV0MassV0DecayRadius);
      fOutputList->Add(fHistSelV0MassV0PropLifeTime);
    }
    if (fIsPrimTrackPairAna) {
      fHistULSPairMassPt_TightCut =
          new TH2F("fHistULSPairMassPt_TightCut", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistULSPairMassPt_TightCut);
      fHistLSppPairMassPt_TightCut =
          new TH2F("fHistLSppPairMassPt_TightCut", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistLSppPairMassPt_TightCut);
      fHistLSmmPairMassPt_TightCut =
          new TH2F("fHistLSmmPairMassPt_TightCut", "",
                   int((max_mass - min_mass) / width_mass), min_mass, max_mass,
                   int((max_pt - min_pt) / width_pt), min_pt, max_pt);
      fOutputList->Add(fHistLSmmPairMassPt_TightCut);
    }

    fHistTrackP = new TH1F("fHistTrackP", "", 100, 0, 10);
    fHistTrackPt = new TH1F("fHistTrackPt", "", 100, 0, 10);
    fHistTrackEta = new TH1F("fHistTrackEta", "", 20, -1, 1);
    fHistTPCNClusts = new TH1F("fHistTPCNClusts", "", 160, 0, 160);
    fHistSPDNClusts = new TH1F("fHistSPDNClusts", "", 6, 0, 6);
    fHistTPCCrossRowsFindableRatio =
        new TH1F("fHistTPCCrossRowsFindableRatio", "", 200, 0, 2);
    fHistReducedChi2TPC = new TH1F("fHistReducedChi2TPC", "", 100, 0, 10);
    fHistReducedChi2ITS = new TH1F("fHistReducedChi2ITS", "", 250, 0, 50);
    fHistDCAz = new TH1F("fHistDCAz", "", 100, 0, 10);
    fHistDCAxyPt = new TH2F("fHistDCAxyPt", "", 200, -1, 1, 100, 0, 10);

    fOutputList->Add(fHistTrackP);
    fOutputList->Add(fHistTrackPt);
    fOutputList->Add(fHistTrackEta);
    fOutputList->Add(fHistTPCNClusts);
    fOutputList->Add(fHistSPDNClusts);
    fOutputList->Add(fHistTPCCrossRowsFindableRatio);
    fOutputList->Add(fHistReducedChi2TPC);
    fOutputList->Add(fHistReducedChi2ITS);
    fOutputList->Add(fHistDCAz);
    fOutputList->Add(fHistDCAxyPt);

    double min_dEdx = 0.;
    double max_dEdx = 300.;
    double width_dEdx = 0.1;

    double min_beta = 0.;
    double max_beta = 1.2;
    double width_beta = 0.01;

    double min_p = 0.;
    double max_p = 5.;
    double width_p = 0.1;

    double min_sigma = -10;
    double max_sigma = +10;
    double width_sigma = 0.1;

    fHistTPCdEdxP =
        new TH2F("fHistTPCdEdxP", "", (max_p - min_p) / width_p, min_p, max_p,
                 (max_dEdx - min_dEdx) / width_dEdx, min_dEdx, max_dEdx);
    fHistBetaP =
        new TH2F("fHistBetaP", "", (max_p - min_p) / width_p, min_p, max_p,
                 (max_beta - min_beta) / width_beta, min_beta, max_beta);
    fHistTPCSigmaElectron = new TH2F(
        "fHistTPCSigmaElectron", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaElectron = new TH2F(
        "fHistTOFSigmaElectron", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTPCSigmaMuon = new TH2F(
        "fHistTPCSigmaMuon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaMuon = new TH2F(
        "fHistTOFSigmaMuon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTPCSigmaPion = new TH2F(
        "fHistTPCSigmaPion", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaPion = new TH2F(
        "fHistTOFSigmaPion", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTPCSigmaKaon = new TH2F(
        "fHistTPCSigmaKaon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaKaon = new TH2F(
        "fHistTOFSigmaKaon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTPCSigmaProton = new TH2F(
        "fHistTPCSigmaProton", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistTOFSigmaProton = new TH2F(
        "fHistTOFSigmaProton", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);

    fHistSelTPCdEdxP =
        new TH2F("fHistSelTPCdEdxP", "", (max_p - min_p) / width_p, min_p,
                 max_p, (max_dEdx - min_dEdx) / width_dEdx, min_dEdx, max_dEdx);
    fHistSelBetaP =
        new TH2F("fHistSelBetaP", "", (max_p - min_p) / width_p, min_p, max_p,
                 (max_beta - min_beta) / width_beta, min_beta, max_beta);
    fHistSelTPCSigmaElectron = new TH2F(
        "fHistSelTPCSigmaElectron", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaElectron = new TH2F(
        "fHistSelTOFSigmaElectron", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTPCSigmaMuon = new TH2F(
        "fHistSelTPCSigmaMuon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaMuon = new TH2F(
        "fHistSelTOFSigmaMuon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTPCSigmaPion = new TH2F(
        "fHistSelTPCSigmaPion", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaPion = new TH2F(
        "fHistSelTOFSigmaPion", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTPCSigmaKaon = new TH2F(
        "fHistSelTPCSigmaKaon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaKaon = new TH2F(
        "fHistSelTOFSigmaKaon", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTPCSigmaProton = new TH2F(
        "fHistSelTPCSigmaProton", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);
    fHistSelTOFSigmaProton = new TH2F(
        "fHistSelTOFSigmaProton", "", (max_p - min_p) / width_p, min_p, max_p,
        (max_sigma - min_sigma) / width_sigma, min_sigma, max_sigma);

    fOutputList->Add(fHistTPCdEdxP);
    fOutputList->Add(fHistBetaP);
    fOutputList->Add(fHistTPCSigmaElectron);
    fOutputList->Add(fHistTOFSigmaElectron);
    fOutputList->Add(fHistTPCSigmaMuon);
    fOutputList->Add(fHistTOFSigmaMuon);
    fOutputList->Add(fHistTPCSigmaPion);
    fOutputList->Add(fHistTOFSigmaPion);
    fOutputList->Add(fHistTPCSigmaKaon);
    fOutputList->Add(fHistTOFSigmaKaon);
    fOutputList->Add(fHistTPCSigmaProton);
    fOutputList->Add(fHistTOFSigmaProton);

    fOutputList->Add(fHistSelTPCdEdxP);
    fOutputList->Add(fHistSelBetaP);
    fOutputList->Add(fHistSelTPCSigmaElectron);
    fOutputList->Add(fHistSelTOFSigmaElectron);
    fOutputList->Add(fHistSelTPCSigmaMuon);
    fOutputList->Add(fHistSelTOFSigmaMuon);
    fOutputList->Add(fHistSelTPCSigmaPion);
    fOutputList->Add(fHistSelTOFSigmaPion);
    fOutputList->Add(fHistSelTPCSigmaKaon);
    fOutputList->Add(fHistSelTOFSigmaKaon);
    fOutputList->Add(fHistSelTPCSigmaProton);
    fOutputList->Add(fHistSelTOFSigmaProton);
  }

  fHistEventVtxZ = new TH1F("fHistEventVtxZ", "", 60, -30, 30);
  fHistEventCent = new TH1F("fHistEventCent", "", 100, 0, 100);
  fHistEventMulti = new TH1F("fHistEventMulti", "", 200, 0, 200);
  fHistEventVtxCont = new TH1F("fHistEventVtxCont", "", 100, 0, 100);
  fOutputList->Add(fHistEventVtxZ);
  fOutputList->Add(fHistEventCent);
  fOutputList->Add(fHistEventMulti);
  fOutputList->Add(fHistEventVtxCont);

  PostData(1, fOutputList);
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPair::UserExec(Option_t *) {

  if (!Initialize()) {
    return;
  }
  if (!fUtils->isAcceptEvent()) {
    return;
  }

  EventQA();

  if (fIsMidTrackAna) {
    if (!fIsMixingAnalysis) {
      FwdMuonPairAnalysis();
    }
    if (fIsMixingAnalysis) {
      FwdMuonPairAnalysisEveMixing();
    }
  }
  if (fIsMidTrackAna) {
    if (!fIsMixingAnalysis) {
      if (fIsV0TrackPairAna) {
        MidV0Analysis(fUtils->getPairTargetPIDs(0),
                      fUtils->getPairTargetPIDs(1));
      }
      if (fIsPrimTrackPairAna) {
        MidPairAnalysis(fUtils->getPairTargetPIDs(0),
                        fUtils->getPairTargetPIDs(1));
      }
    }
    if (fIsMixingAnalysis) {
      if (fIsV0TrackPairAna) {
        MidV0AnalysisEventMixing(fUtils->getPairTargetPIDs(0),
                                 fUtils->getPairTargetPIDs(1));
      }
    }
  }
}

bool AliAnalysisTaskAODTrackPair::Initialize() {
  fEvent = dynamic_cast<AliAODEvent *>(InputEvent());

  if (!fUtils->setEvent(fEvent, fInputHandler)) {
    return false;
  }

  if (fRunNumber != fEvent->GetRunNumber()) {
    fRunNumber = fUtils->getRunnumber();
    if (!fIsMidTrackAna) {
      AliMuonTrackCuts *trackCut = fUtils->getMuonTrackCuts();
      trackCut->SetRun(fInputHandler);
    }
  }

  if (!fIsMidTrackAna) {
    fUtils->getTriggerInfo(fIsCINT7, fIsCMSL7, fIsCMSH7, fIsCMUL7, fIsCMLL7);
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::EventQA() {
  fHistEventVtxZ->Fill(fUtils->getVtxZ());
  fHistEventCent->Fill(fUtils->getCentClass());
  fHistEventMulti->Fill(fUtils->getNCorrSPDTrkInfo(1));
  fHistEventVtxCont->Fill(fUtils->getVtxCont());
  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonTrackQA(AliAODTrack *track) {
  fHistTrackEta->Fill(track->Pt(), track->Eta());
  fHistTrackThetaAbs->Fill(track->Pt(),
                           AliAnalysisMuonUtility::GetThetaAbsDeg(track));
  fHistTrackTriggerMatch->Fill(track->Pt(),
                               AliAnalysisMuonUtility::GetMatchTrigger(track));
  fHistTrackChiSquare->Fill(
      track->Pt(), AliAnalysisMuonUtility::GetChi2perNDFtracker(track));
  fHistTriggerChiSquare->Fill(
      track->Pt(), AliAnalysisMuonUtility::GetChi2MatchTrigger(track));
  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonPairQA(AliAODDimuon *dimuon) {

  AliAODTrack *track1 = dynamic_cast<AliAODTrack *>(dimuon->GetMu(0));
  AliAODTrack *track2 = dynamic_cast<AliAODTrack *>(dimuon->GetMu(1));

  int triggerLB1 = AliAnalysisMuonUtility::GetLoCircuit(track1);
  int triggerLB2 = AliAnalysisMuonUtility::GetLoCircuit(track2);

  float pt_max = track1->Pt();
  float pt_min = track2->Pt();

  if (track1->Pt() > track2->Pt()) {
    pt_max = track1->Pt();
    pt_min = track2->Pt();
  } else {
    pt_max = track2->Pt();
    pt_min = track1->Pt();
  }

  fHistTrackPairPtBalance->Fill(pt_max, pt_min);
  fHistTrackPairLocalBoardPair->Fill(triggerLB1, triggerLB2);

  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonPairAnalysisEveMixing() {

  if (!(fInputHandler->IsEventSelected() & fTriggerMaskForMixing))
    return false;

  TObjArray *fTrackArray = new TObjArray();
  fTrackArray->SetOwner();

  float poolCent = 0.;
  float poolVtxZ = 0.;
  float poolPsi = 0.;

  if (onEvtMixingPoolVtxZ) {
    poolVtxZ = fUtils->getVtxZ();
  }
  if (onEvtMixingPoolCent) {
    // poolCent=fUtils->getCentClass();
    poolCent = fUtils->getNCorrSPDTrkInfo(1);
  }
  if (onEvtMixingPoolPsi) {
    poolPsi = fUtils->getPsi();
  }

  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgr->GetEventPool(
      poolCent, poolVtxZ, poolPsi);

  Int_t nTrack = fEvent->GetNumberOfTracks();

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    AliAODTrack *track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    if (!fUtils->isAcceptFwdMuonTrack(track1))
      continue;

    if (pool->IsReady()) {

      for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {

        TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);

        for (Int_t iTrack2 = 0; iTrack2 < poolTracks->GetEntriesFast();
             ++iTrack2) {

          AliAODTrack *__track2__ = (AliAODTrack *)poolTracks->At(iTrack2);

          AliAODTrack *track2 = (AliAODTrack *)__track2__->Clone();

          if (!fUtils->isAcceptFwdMuonTrack(track2))
            continue;

          AliAODDimuon *dimuon = new AliAODDimuon();
          dimuon->SetMuons(track1, track2);

          if (!fUtils->isAcceptFwdDimuon(dimuon))
            continue;

          RecPairPt = dimuon->Pt();
          RecPairRap = fabs(dimuon->Y());
          RecPairMass = dimuon->M();
          // RecPairCent = fUtils->getCentClass();
          RecPairCent = fUtils->getNCorrSPDTrkInfo(1);
          RecPairDS = fUtils->getDS();

          string fFiredTrigName = string(fEvent->GetFiredTriggerClasses());

          if (dimuon->Charge() == 0) {
            fTreeMixULSPair->Fill();
          } else if (dimuon->Charge() > 0) {
            fTreeMixLSppPair->Fill();
          } else {
            fTreeMixLSmmPair->Fill();
          }

          delete track2;
          delete dimuon;

        } // end of loop track2
      }   // end of loop iMixEvt
    }     // poolPion->IsReady()

    fTrackArray->Add(track1);

  } // end of loop track1

  TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  if (fTrackArrayClone->GetEntriesFast() > 0) {
    pool->UpdatePool(fTrackArrayClone);
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::FwdMuonPairAnalysis() {

  if (!fIsMC && !(fInputHandler->IsEventSelected() & fTriggerMaskForSame)) {
    return false;
  }

  if (fIsCMUL7) {
    fEventCounter->Fill(0., fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(4., fUtils->getNCorrSPDTrkInfo(1),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMLL7) {
    fEventCounter->Fill(1., fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(5., fUtils->getNCorrSPDTrkInfo(1),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMUL7 || fIsCMLL7) {
    fEventCounter->Fill(2., fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(6., fUtils->getNCorrSPDTrkInfo(1),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMUL7 && fIsCMLL7) {
    fEventCounter->Fill(3., fUtils->getNCorrSPDTrkInfo(1));
    fEventCounter->Fill(7., fUtils->getNCorrSPDTrkInfo(1),
                        (double)1. / fUtils->getDS());
  }

  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack *track1;
  AliAODTrack *track2;

  AliAODDimuon *dimuon;

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    if (!fUtils->isAcceptFwdMuonTrack(track1))
      continue;

    FwdMuonTrackQA(track1);

    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTrack; ++iTrack2) {

      track2 = (AliAODTrack *)fEvent->GetTrack(iTrack2);

      if (!fUtils->isAcceptFwdMuonTrack(track2))
        continue;

      dimuon = new AliAODDimuon();
      dimuon->SetMuons(track1, track2);

      if (!fUtils->isAcceptFwdDimuon(dimuon))
        continue;

      FwdMuonPairQA(dimuon);

      RecPairPt = dimuon->Pt();
      RecPairRap = fabs(dimuon->Y());
      RecPairMass = dimuon->M();
      // RecPairCent = fUtils->getCentClass();
      RecPairCent = fUtils->getNCorrSPDTrkInfo(1);
      RecPairDS = fUtils->getDS();

      if (dimuon->Charge() == 0) {
        fTreeULSPair->Fill();
      } else if (dimuon->Charge() > 0) {
        fTreeLSppPair->Fill();
      } else {
        fTreeLSmmPair->Fill();
      }

      delete dimuon;

    } // end of loop track2
  }   // end of loop track1
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidTrackQualityChecker(AliAODTrack *track) {
  fHistTPCNClusts->Fill(track->GetTPCNcls());
  int nSPD = 0;
  if (track->HasPointOnITSLayer(0)) {
    ++nSPD;
  }
  if (track->HasPointOnITSLayer(1)) {
    ++nSPD;
  }
  fHistSPDNClusts->Fill(nSPD);
  fHistTPCCrossRowsFindableRatio->Fill(track->GetTPCNclsF() /
                                       track->GetTPCCrossedRows());
  fHistReducedChi2TPC->Fill(track->GetTPCchi2() / track->GetTPCNcls());

  if (track->GetITSNcls() > 0) {
    fHistReducedChi2ITS->Fill(track->GetITSchi2() / track->GetITSNcls());
  }

  float dca_xy = 9999;
  float dca_z = 9999;
  track->GetImpactParameters(dca_xy, dca_z);

  fHistDCAz->Fill(dca_z);
  fHistDCAxyPt->Fill(dca_xy, track->Pt());

  fHistTrackP->Fill(track->P());
  fHistTrackPt->Fill(track->Pt());
  fHistTrackEta->Fill(track->Eta());

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidV0Checker(AliAODv0 *v0, bool isSel) {

  double vtx[] = {fUtils->getVtxX(), fUtils->getVtxY(), fUtils->getVtxZ()};

  if (isSel) {
    fHistSelArmenteros->Fill(v0->Alpha(), v0->PtArmV0());
    fHistSelV0MassDecayLength->Fill(v0->MassK0Short(), v0->DecayLengthV0(vtx));
    fHistSelV0MassPointingAngle->Fill(v0->MassK0Short(),
                                      v0->CosPointingAngle(vtx));
    fHistSelV0MassV0DCA->Fill(v0->MassK0Short(), v0->DcaV0ToPrimVertex());
    fHistSelV0MassV0TrackDCA->Fill(v0->MassK0Short(), v0->DcaV0Daughters());
    fHistSelV0MassV0DecayRadius->Fill(v0->MassK0Short(), v0->RadiusV0());
    fHistSelV0MassV0PropLifeTime->Fill(v0->MassK0Short(),
                                       fUtils->fPdgK0sMass *
                                           v0->DecayLengthV0(vtx) / v0->P());
  } else {
    fHistArmenteros->Fill(v0->Alpha(), v0->PtArmV0());
    fHistV0MassDecayLength->Fill(v0->MassK0Short(), v0->DecayLengthV0(vtx));
    fHistV0MassPointingAngle->Fill(v0->MassK0Short(),
                                   v0->CosPointingAngle(vtx));
    fHistV0MassV0DCA->Fill(v0->MassK0Short(), v0->DcaV0ToPrimVertex());
    fHistV0MassV0TrackDCA->Fill(v0->MassK0Short(), v0->DcaV0Daughters());
    fHistV0MassV0DecayRadius->Fill(v0->MassK0Short(), v0->RadiusV0());
    fHistV0MassV0PropLifeTime->Fill(v0->MassK0Short(),
                                    fUtils->fPdgK0sMass *
                                        v0->DecayLengthV0(vtx) / v0->P());
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidTrackPIDChecker(AliAODTrack *track,
                                                     AliPID::EParticleType pid,
                                                     bool isSel) {

  float p = track->P();
  float sigTOF = track->GetTOFsignal();
  float length = track->GetIntegratedLength();
  float beta =
      (sigTOF > 0) ? (double)length / (2.99792457999999984e-02 * sigTOF) : -999;
  float dEdx = track->GetTPCsignal();

  if (isSel) {
    fHistSelTPCdEdxP->Fill(p, dEdx);
    fHistSelTPCSigmaElectron->Fill(
        p, fUtils->getTPCSigma(track, AliPID::kElectron));
    fHistSelTPCSigmaElectron->Fill(
        p, fUtils->getTPCSigma(track, AliPID::kElectron));
    fHistSelTPCSigmaMuon->Fill(p, fUtils->getTPCSigma(track, AliPID::kMuon));
    fHistSelTPCSigmaMuon->Fill(p, fUtils->getTPCSigma(track, AliPID::kMuon));
    fHistSelTPCSigmaPion->Fill(p, fUtils->getTPCSigma(track, AliPID::kPion));
    fHistSelTPCSigmaPion->Fill(p, fUtils->getTPCSigma(track, AliPID::kPion));
    fHistSelTPCSigmaKaon->Fill(p, fUtils->getTPCSigma(track, AliPID::kKaon));
    fHistSelTPCSigmaKaon->Fill(p, fUtils->getTPCSigma(track, AliPID::kKaon));
    fHistSelTPCSigmaProton->Fill(p,
                                 fUtils->getTPCSigma(track, AliPID::kProton));
    fHistSelTPCSigmaProton->Fill(p,
                                 fUtils->getTPCSigma(track, AliPID::kProton));
    if (beta > 0.) {
      fHistSelBetaP->Fill(p, beta);
      fHistSelTOFSigmaElectron->Fill(
          p, fUtils->getTOFSigma(track, AliPID::kElectron));
      fHistSelTOFSigmaElectron->Fill(
          p, fUtils->getTOFSigma(track, AliPID::kElectron));
      fHistSelTOFSigmaMuon->Fill(p, fUtils->getTOFSigma(track, AliPID::kMuon));
      fHistSelTOFSigmaMuon->Fill(p, fUtils->getTOFSigma(track, AliPID::kMuon));
      fHistSelTOFSigmaPion->Fill(p, fUtils->getTOFSigma(track, AliPID::kPion));
      fHistSelTOFSigmaPion->Fill(p, fUtils->getTOFSigma(track, AliPID::kPion));
      fHistSelTOFSigmaKaon->Fill(p, fUtils->getTOFSigma(track, AliPID::kKaon));
      fHistSelTOFSigmaKaon->Fill(p, fUtils->getTOFSigma(track, AliPID::kKaon));
      fHistSelTOFSigmaProton->Fill(p,
                                   fUtils->getTOFSigma(track, AliPID::kProton));
      fHistSelTOFSigmaProton->Fill(p,
                                   fUtils->getTOFSigma(track, AliPID::kProton));
    }
  } else {
    fHistTPCdEdxP->Fill(p, dEdx);
    fHistTPCSigmaElectron->Fill(p,
                                fUtils->getTPCSigma(track, AliPID::kElectron));
    fHistTPCSigmaElectron->Fill(p,
                                fUtils->getTPCSigma(track, AliPID::kElectron));
    fHistTPCSigmaMuon->Fill(p, fUtils->getTPCSigma(track, AliPID::kMuon));
    fHistTPCSigmaMuon->Fill(p, fUtils->getTPCSigma(track, AliPID::kMuon));
    fHistTPCSigmaPion->Fill(p, fUtils->getTPCSigma(track, AliPID::kPion));
    fHistTPCSigmaPion->Fill(p, fUtils->getTPCSigma(track, AliPID::kPion));
    fHistTPCSigmaKaon->Fill(p, fUtils->getTPCSigma(track, AliPID::kKaon));
    fHistTPCSigmaKaon->Fill(p, fUtils->getTPCSigma(track, AliPID::kKaon));
    fHistTPCSigmaProton->Fill(p, fUtils->getTPCSigma(track, AliPID::kProton));
    fHistTPCSigmaProton->Fill(p, fUtils->getTPCSigma(track, AliPID::kProton));
    if (beta > 0.) {
      fHistBetaP->Fill(p, beta);
      fHistTOFSigmaElectron->Fill(
          p, fUtils->getTOFSigma(track, AliPID::kElectron));
      fHistTOFSigmaElectron->Fill(
          p, fUtils->getTOFSigma(track, AliPID::kElectron));
      fHistTOFSigmaMuon->Fill(p, fUtils->getTOFSigma(track, AliPID::kMuon));
      fHistTOFSigmaMuon->Fill(p, fUtils->getTOFSigma(track, AliPID::kMuon));
      fHistTOFSigmaPion->Fill(p, fUtils->getTOFSigma(track, AliPID::kPion));
      fHistTOFSigmaPion->Fill(p, fUtils->getTOFSigma(track, AliPID::kPion));
      fHistTOFSigmaKaon->Fill(p, fUtils->getTOFSigma(track, AliPID::kKaon));
      fHistTOFSigmaKaon->Fill(p, fUtils->getTOFSigma(track, AliPID::kKaon));
      fHistTOFSigmaProton->Fill(p, fUtils->getTOFSigma(track, AliPID::kProton));
      fHistTOFSigmaProton->Fill(p, fUtils->getTOFSigma(track, AliPID::kProton));
    }
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidTrackQA(AliAODTrack *track) {

  float dca_xy = 9999;
  float dca_z = 9999;
  track->GetImpactParameters(dca_xy, dca_z);

  AliTOFHeader *tofHeader = (AliTOFHeader *)track->GetTOFHeader();
  /*
  fTrackPt = track->Pt();
  fTrackP = track->P();
  fTrackTheta = track->Theta();
  fTrackPhi = track->Phi();
  fTrackLength = track->GetIntegratedLength();
  fTrackBeta = beta;
  fTrackTrackChi2perNDF = track->Chi2perNDF();
  fTrackTrackITSNcls = track->GetITSNcls();
  fTrackTrackTPCNcls = track->GetTPCNcls();
  fTrackTrackTOFNcls = tofHeader->GetNumberOfTOFclusters();
  fTrackTrackTPCChi2 = track->GetTPCchi2();
  fTrackTrackITSChi2 = track->GetITSchi2();
  fTrackTPCCrossedRows = track->GetTPCCrossedRows();
  fTrackTPCFindableNcls = track->GetTPCNclsF();
  fTrackTOFBCTime = track->GetTOFBunchCrossing();
  fTrackTOFKinkIndex = track->GetKinkIndex(0);
  fTrackDCAxy = dca_xy;
  fTrackDCAz = dca_z;
  fTrackTPCsigmaMuon = fUtils->getTPCSigma(track,AliPID::kMuon);
  fTrackTOFsigmaMuon = fUtils->getTOFSigma(track,AliPID::kMuon);
  */
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidMuonPairQA(AliAODDimuon *dimuon) {
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidV0Analysis(AliPID::EParticleType pid1,
                                                AliPID::EParticleType pid2) {

  Int_t nV0 = fEvent->GetNumberOfV0s();

  AliAODv0 *v0_1;
  AliAODv0 *v0_2;

  TLorentzVector lv1, lv2, lv12;

  for (int iV0_1 = 0; iV0_1 < nV0; ++iV0_1) {

    v0_1 = (AliAODv0 *)fEvent->GetV0(iV0_1);

    if (!fUtils->isAcceptV0Kinematics(v0_1)) {
      continue;
    }

    RecPairPt = v0_1->Pt();
    RecPairMass = v0_1->MassK0Short();
    RecPairRap = v0_1->RapK0Short();
    RecPairArmenterosArmPt = v0_1->PtArmV0();
    RecPairArmenterosAlpha = v0_1->AlphaV0();

    AliAODTrack *pTrack = (AliAODTrack *)v0_1->GetDaughter(0);
    AliAODTrack *nTrack = (AliAODTrack *)v0_1->GetDaughter(1);

    MidV0Checker(v0_1, false);

    if (0.4 > RecPairMass || RecPairMass > 0.6) {
      continue;
    }
    if (!fUtils->isAcceptedK0s(v0_1, pid1, pid2, 0)) {
      continue;
    }
    if (!fUtils->isAcceptArmenterosK0s(v0_1)) {
      continue;
    }

    MidV0Checker(v0_1, true);

    MidTrackQualityChecker(pTrack);
    MidTrackQualityChecker(nTrack);

    MidTrackPIDChecker(pTrack, AliPID::kPion, false);
    MidTrackPIDChecker(nTrack, AliPID::kPion, false);

    fHistSelArmenteros->Fill(v0_1->Alpha(), v0_1->PtArmV0());
    fHistULSPairMassPt_ProngV0->Fill(RecPairMass, RecPairPt);

    if (!fUtils->isAcceptK0sCandidateMassRange(v0_1->MassK0Short())) {
      continue;
    }

    MidTrackPIDChecker(pTrack, AliPID::kPion, true);
    MidTrackPIDChecker(nTrack, AliPID::kPion, true);

    for (int iV0_2 = iV0_1 + 1; iV0_2 < nV0; ++iV0_2) {

      v0_2 = (AliAODv0 *)fEvent->GetV0(iV0_2);

      RecPairPt = v0_2->Pt();
      RecPairMass = v0_2->MassK0Short();
      RecPairRap = v0_2->RapK0Short();

      if (0.4 > RecPairMass || RecPairMass > 0.6) {
        continue;
      }

      if (!fUtils->isAcceptV0Kinematics(v0_2)) {
        continue;
      }
      if (!fUtils->isAcceptK0sCandidateMassRange(v0_2->MassK0Short())) {
        continue;
      }
      if (!fUtils->isAcceptedK0s(v0_2, pid1, pid2, 0)) {
        continue;
      }
      if (!fUtils->isAcceptArmenterosK0s(v0_2)) {
        continue;
      }

      fHistMassK0s1K0s2->Fill(v0_1->MassK0Short(), v0_2->MassK0Short());

      lv1.SetPtEtaPhiM(v0_1->Pt(), v0_1->Eta(), v0_1->Phi(),
                       TDatabasePDG::Instance()->GetParticle(310)->Mass());
      lv2.SetPtEtaPhiM(v0_2->Pt(), v0_2->Eta(), v0_2->Phi(),
                       TDatabasePDG::Instance()->GetParticle(310)->Mass());
      lv12 = lv1 + lv2;

      RecPairPt = lv12.Pt();
      RecPairMass = lv12.M();
      RecPairRap = lv12.Rapidity();

      // fHistULSPairMassPt->Fill(RecPairMass,RecPairPt);
      double fill[] = {RecPairMass, RecPairPt, fUtils->getCentClass(1)};
      fSparseULSPairMassPt->Fill(fill);
    }
  }
  return true;
}

bool AliAnalysisTaskAODTrackPair::MidV0AnalysisEventMixing(
    AliPID::EParticleType pid1, AliPID::EParticleType pid2) {

  TObjArray *fTrackArray = new TObjArray();
  fTrackArray->SetOwner();

  float poolCent = 0.;
  float poolVtxZ = 0.;
  float poolPsi = 0.;

  if (onEvtMixingPoolVtxZ) {
    poolVtxZ = fUtils->getVtxZ();
  }
  if (onEvtMixingPoolCent) {
    poolCent = fUtils->getNCorrSPDTrkInfo(1);
  }
  if (onEvtMixingPoolPsi) {
    poolPsi = fUtils->getPsi();
  }

  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgr->GetEventPool(
      poolCent, poolVtxZ, poolPsi);

  Int_t nV0 = fEvent->GetNumberOfV0s();

  AliAODv0 *v0_1;
  AliAODv0 *v0_2;

  TLorentzVector lv1, lv2, lv12;

  for (int iV0_1 = 0; iV0_1 < nV0; ++iV0_1) {

    v0_1 = (AliAODv0 *)fEvent->GetV0(iV0_1);

    RecPairPt = v0_1->Pt();
    RecPairMass = v0_1->MassK0Short();
    RecPairRap = v0_1->RapK0Short();

    AliAODTrack *pTrack = (AliAODTrack *)v0_1->GetDaughter(0);
    AliAODTrack *nTrack = (AliAODTrack *)v0_1->GetDaughter(1);

    MidV0Checker(v0_1, false);

    if (0.4 > RecPairMass || RecPairMass > 0.6) {
      continue;
    }
    if (!fUtils->isAcceptedK0s(v0_1, pid1, pid2, 0)) {
      continue;
    }
    if (!fUtils->isAcceptArmenterosK0s(v0_1)) {
      continue;
    }

    MidV0Checker(v0_1, true);

    MidTrackQualityChecker(pTrack);
    MidTrackQualityChecker(nTrack);

    MidTrackPIDChecker(pTrack, AliPID::kPion, false);
    MidTrackPIDChecker(nTrack, AliPID::kPion, false);

    fHistSelArmenteros->Fill(v0_1->Alpha(), v0_1->PtArmV0());
    fHistULSPairMassPt_ProngV0->Fill(RecPairMass, RecPairPt);

    if (!fUtils->isAcceptK0sCandidateMassRange(v0_1->MassK0Short())) {
      continue;
    }

    MidTrackPIDChecker(pTrack, AliPID::kPion, true);
    MidTrackPIDChecker(nTrack, AliPID::kPion, true);

    if (pool->IsReady()) {

      for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {

        TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);

        for (int iV0_2 = 0; iV0_2 < poolTracks->GetEntriesFast(); ++iV0_2) {

          v0_2 = (AliAODv0 *)poolTracks->At(iV0_2);

          lv1.SetPtEtaPhiM(v0_1->Pt(), v0_1->Eta(), v0_1->Phi(),
                           TDatabasePDG::Instance()->GetParticle(310)->Mass());
          lv2.SetPtEtaPhiM(v0_2->Pt(), v0_2->Eta(), v0_2->Phi(),
                           TDatabasePDG::Instance()->GetParticle(310)->Mass());
          lv12 = lv1 + lv2;

          RecPairPt = lv12.Pt();
          RecPairMass = lv12.M();
          RecPairRap = lv12.Rapidity();

          double fill[] = {RecPairMass, RecPairPt, fUtils->getCentClass(1)};
          fSparseULSPairMassPt->Fill(fill);
          // fHistMixULSPairMassPt->Fill(RecPairMass,RecPairPt);
        }
      }
    }

    fTrackArray->Add(v0_1);
  }

  TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  if (fTrackArrayClone->GetEntriesFast() > 0) {
    pool->UpdatePool(fTrackArrayClone);
  }

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidPairAnalysis(AliPID::EParticleType pid1,
                                                  AliPID::EParticleType pid2) {
  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack *track1;
  AliAODTrack *track2;

  TLorentzVector lv1, lv2, lv12;

  std::vector<TLorentzVector> tracks;
  std::vector<int> charges;
  std::vector<bool> skip_tracks_ids;

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    if (!fUtils->isAcceptMidPrimTrackQuality(track1)) {
      continue;
    }
    if (!fUtils->isAcceptTrackKinematics(track1)) {
      continue;
    }

    MidTrackPIDChecker(track1, pid1, false);

    if (!fUtils->isAcceptMidPid(track1, pid1)) {
      continue;
    }

    MidTrackPIDChecker(track1, pid1, true);
    MidTrackQualityChecker(track1);

    float mass1 = 0;

    if (pid1 == AliPID::kElectron) {
      mass1 = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    } else if (pid1 == AliPID::kPion) {
      mass1 = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    } else if (pid1 == AliPID::kKaon) {
      mass1 = TDatabasePDG::Instance()->GetParticle(321)->Mass();
    } else if (pid1 == AliPID::kProton) {
      mass1 = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    } else if (pid1 == AliPID::kMuon) {
      mass1 = TDatabasePDG::Instance()->GetParticle(13)->Mass();
    } else {
      continue;
    }

    TLorentzVector vec4;
    vec4.SetPtEtaPhiM(track1->Pt(), track1->Eta(), track1->Phi(), mass1);
    tracks.push_back(vec4);
    skip_tracks_ids.push_back(false);
    charges.push_back(track1->Charge());
  }

  nTrack = tracks.size();

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    lv1 = tracks[iTrack1];

    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTrack; ++iTrack2) {

      lv2 = tracks[iTrack2];

      if (charges[iTrack1] + charges[iTrack2] != 0) {
        continue;
      }

      lv12 = lv1 + lv2;
      /*
      if ( (pid1 == AliPID::kKaon && pid2 == AliPID::kKaon) && lv12.M() < 1.035
      ) { skip_tracks_ids[iTrack1] = true; skip_tracks_ids[iTrack2] = true; }
      else if ( (pid1 == AliPID::kPion && pid2 == AliPID::kPion) ) { if
      (0.48<lv12.M() && lv12.M()<0.51) { skip_tracks_ids[iTrack1] = true;
          skip_tracks_ids[iTrack2] = true;
        } else if (0.70<lv12.M() && lv12.M()<0.84) {
          skip_tracks_ids[iTrack1] = true;
          skip_tracks_ids[iTrack2] = true;
        }
      }
      */
    } // end of loop track2
  }   // end of loop track1

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    lv1 = tracks[iTrack1];

    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTrack; ++iTrack2) {

      lv2 = tracks[iTrack2];

      lv12 = lv1 + lv2;

      RecPairPt = lv12.Pt();
      RecPairMass = lv12.M();

      if (charges[iTrack1] + charges[iTrack2] == 0) {
        fHistULSPairMassPt->Fill(RecPairMass, RecPairPt);
      } else if (charges[iTrack1] + charges[iTrack2] > 0) {
        fHistLSppPairMassPt->Fill(RecPairMass, RecPairPt);
      } else {
        fHistLSmmPairMassPt->Fill(RecPairMass, RecPairPt);
      }

      if (skip_tracks_ids[iTrack1] && skip_tracks_ids[iTrack2]) {
        continue;
      }

      if (charges[iTrack1] + charges[iTrack2] == 0) {
        fHistULSPairMassPt_TightCut->Fill(RecPairMass, RecPairPt);
      } else if (charges[iTrack1] + charges[iTrack2] > 0) {
        fHistLSppPairMassPt_TightCut->Fill(RecPairMass, RecPairPt);
      } else {
        fHistLSmmPairMassPt_TightCut->Fill(RecPairMass, RecPairPt);
      }

    } // end of loop track2
  }   // end of loop track1

  return true;
}

bool AliAnalysisTaskAODTrackPair::MidMuonPairAnalysisEveMixing() {

  TObjArray *fTrackArray = new TObjArray();
  fTrackArray->SetOwner();

  float poolCent = 0.;
  float poolVtxZ = 0.;
  float poolPsi = 0.;

  if (onEvtMixingPoolVtxZ)
    poolVtxZ = fUtils->getVtxZ();
  if (onEvtMixingPoolCent)
    poolCent = fUtils->getCentClass();
  if (onEvtMixingPoolPsi)
    poolPsi = fUtils->getPsi();

  AliEventPool *pool = (AliEventPool *)fPoolMuonTrackMgr->GetEventPool(
      poolCent, poolVtxZ, poolPsi);

  Int_t nTrack = fEvent->GetNumberOfTracks();

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    AliAODTrack *track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    if (!fUtils->isAcceptMidMuonTrack(track1))
      continue;

    if (pool->IsReady()) {

      for (Int_t iMixEvt = 0; iMixEvt < pool->GetCurrentNEvents(); iMixEvt++) {

        TObjArray *poolTracks = (TObjArray *)pool->GetEvent(iMixEvt);

        for (Int_t iTrack2 = 0; iTrack2 < poolTracks->GetEntriesFast();
             ++iTrack2) {

          AliAODTrack *__track2__ = (AliAODTrack *)poolTracks->At(iTrack2);

          AliAODTrack *track2 = (AliAODTrack *)__track2__->Clone();

          if (!fUtils->isAcceptMidMuonTrack(track2))
            continue;

          AliAODDimuon *dimuon = new AliAODDimuon();
          dimuon->SetMuons(track1, track2);

          if (!fUtils->isAcceptMidDimuon(dimuon))
            continue;

          double fill[] = {dimuon->M(), fabs(dimuon->Y()), dimuon->Pt(),
                           fUtils->getCentClass()};

          RecPairPt = dimuon->Pt();
          RecPairRap = fabs(dimuon->Y());
          RecPairMass = dimuon->M();
          RecPairCent = fUtils->getCentClass();
          RecPairDS = fUtils->getDS();

          if (dimuon->Charge() == 0) {
            fTreeMixULSPair->Fill();
          } else if (dimuon->Charge() > 0) {
            fTreeMixLSppPair->Fill();
          } else {
            fTreeMixLSmmPair->Fill();
          }

          delete track2;
          delete dimuon;

        } // end of loop track2
      }   // end of loop iMixEvt
    }     // poolPion->IsReady()

    fTrackArray->Add(track1);

  } // end of loop track1

  TObjArray *fTrackArrayClone = (TObjArray *)fTrackArray->Clone();
  fTrackArrayClone->SetOwner();
  if (fTrackArrayClone->GetEntriesFast() > 0) {
    pool->UpdatePool(fTrackArrayClone);
  }

  return true;
}
