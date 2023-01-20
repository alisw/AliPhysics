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

#include "AliAnalysisTaskAODTrackPairMC.h"
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

#include "AliAnalysisTaskAODTrackPairUtils.h"

#include "iostream"
#include "memory"
// Authors: Satoshi Yano
// Reviewed:

using namespace std;
// AliAnalysisTaskAODTrackPairMCUtils.h
ClassImp(AliAnalysisTaskAODTrackPairMC)

    AliAnalysisTaskAODTrackPairMC::AliAnalysisTaskAODTrackPairMC()
    : AliAnalysisTaskSE(), fEvent(NULL), fRandom(NULL), fPoolMuonTrackMgr(NULL),
      fUtils(NULL), fMCTrackArray(NULL), fIsMC(false), fIsMixingAnalysis(false),
      fRunNumber(-99999), fTrackDepth(1000), fPoolSize(1), fReadyFraction(0.1),
      fTriggerMaskForSame(AliVEvent::kMuonUnlikeLowPt7 |
                          AliVEvent::kMuonLikeLowPt7),
      fTriggerMaskForMixing(AliVEvent::kMuonSingleLowPt7),
      onEvtMixingPoolVtxZ(true), onEvtMixingPoolCent(true),
      onEvtMixingPoolPsi(true), fIsCINT7(false), fIsCMSL7(false),
      fIsCMSH7(false), fIsCMUL7(false), fIsCMLL7(false),

      fOutputList(NULL), fEventCounter(NULL),

      fHistTrackEta(NULL), fHistTrackThetaAbs(NULL),
      fHistTrackTriggerMatch(NULL), fHistTrackPDCA(NULL),
      fHistTrackChiSquare(NULL), fHistTriggerChiSquare(NULL),

      fHistEventVtxZ(NULL), fHistEventCent(NULL), fHistEventMulti(NULL),
      fHistEventVtxCont(NULL),

      fTreeRecMuonP(NULL), fTreeRecMuonN(NULL), RecMuonPt(0.), RecMuonEta(0.),
      RecMuonRap(0.), RecMuonPhi(0.), RecMuonThetaAbs(0.),
      RecMuonTriggerMatch(0), RecMuonChiSquare(0.), RecMuonTriggerChiSquare(0.),
      RecMuonIsGoodTrack(0), RecMCMuonPt(0.), RecMCMuonEta(0.),
      RecMCMuonRap(0.), RecMCMuonPhi(0.),

      fTreeMCMuonP(NULL), fTreeMCMuonN(NULL), MCMuonPt(0.), MCMuonEta(0.),
      MCMuonRap(0.), MCMuonPhi(0.), MCMuonDetect(0.),

      fTreeULSDimuon(NULL), fTreeLSppDimuon(NULL), fTreeLSmmDimuon(NULL),

      fTreeMCULSDimuon(NULL), fTreeMCLSppDimuon(NULL), fTreeMCLSmmDimuon(NULL),

      fTreeMixULSDimuon(NULL), fTreeMixLSppDimuon(NULL),
      fTreeMixLSmmDimuon(NULL),

      RecDimuonPt(0.), RecDimuonRap(0.), RecDimuonMass(0.), RecDimuonCent(0.),
      RecDimuonDS(0.), RecMCDimuonPt(0.), RecMCDimuonRap(0.),
      RecMCDimuonMass(0.), RecMCMom2Body(0.), RecMCMomDalitz(0),
      RecMCMomPdgCode(0.), RecMCMomMass(0), RecMCMomEta(0), RecMCMomPt(0),

      MCDimuonPt(0.), MCDimuonRap(0.), MCDimuonMass(0.), MCMom2Body(0.),
      MCMomDalitz(0), MCMomPt(0), MCMomEta(0), MCMomPdgCode(0.),
      MCDimuonDetected(0),

      fIsMidTrackAna(0),

      fTreeTrack(NULL), fTreeTrackLight(NULL), fTrackPt(0.), fTrackP(0.),
      fTrackTheta(0.), fTrackEta(0.), fTrackPhi(0.), fTrackLength(0.),
      fTrackBeta(0.), fTrackdEdx(0.), fTrackTrackChi2perNDF(0.),
      fTrackTrigChi2perNDF(0.), fTrackTrackITSNcls(0.), fTrackTrackTPCNcls(0.),
      fTrackTrackTOFNcls(0.), fTrackTrackTPCChi2(0.), fTrackTrackITSChi2(0.),
      fTrackTPCCrossedRows(0.), fTrackTPCFindableNcls(0.), fTrackTOFBCTime(0.),
      fTrackTOFKinkIndex(0.), fTrackDCAxy(0.), fTrackDCAz(0.),
      fTrackTPCsigmaKaon(0.), fTrackTOFsigmaKaon(0.), fTrackTPCsigmaPion(0.),
      fTrackTOFsigmaPion(0.), fTrackTPCsigmaProton(0.),
      fTrackTOFsigmaProton(0.), fTrackTPCsigmaElectron(0.),
      fTrackTOFsigmaElectron(0.), fTrackTPCsigmaMuon(0.),
      fTrackTOFsigmaMuon(0.), fTrackThetaAbs(0.), fTrackGlobal(0),
      fTrackGlobalNoDCA(0), fTrackTPConly(0), fTrackGoodFwdQuarity(0),
      fTrackMCType(0), fTrackTrigMatch(0), fTrackPdgCode(0),
      fTrackMotherPdgCode(0) {}

AliAnalysisTaskAODTrackPairMC::AliAnalysisTaskAODTrackPairMC(const char *name)
    : AliAnalysisTaskSE(name), fEvent(NULL), fPoolMuonTrackMgr(NULL),
      fUtils(NULL), fMCTrackArray(NULL), fIsMC(false), fIsMixingAnalysis(false),
      fRunNumber(-99999), fTrackDepth(1000), fPoolSize(1), fReadyFraction(0.1),
      fTriggerMaskForSame(AliVEvent::kMuonUnlikeLowPt7 |
                          AliVEvent::kMuonLikeLowPt7),
      fTriggerMaskForMixing(AliVEvent::kMuonSingleLowPt7),
      onEvtMixingPoolVtxZ(true), onEvtMixingPoolCent(true),
      onEvtMixingPoolPsi(true), fIsCINT7(false), fIsCMSL7(false),
      fIsCMSH7(false), fIsCMUL7(false), fIsCMLL7(false),

      fOutputList(NULL), fEventCounter(NULL),

      fHistTrackEta(NULL), fHistTrackThetaAbs(NULL),
      fHistTrackTriggerMatch(NULL), fHistTrackPDCA(NULL),
      fHistTrackChiSquare(NULL), fHistTriggerChiSquare(NULL),

      fHistEventVtxZ(NULL), fHistEventCent(NULL), fHistEventMulti(NULL),
      fHistEventVtxCont(NULL),

      fTreeRecMuonP(NULL), fTreeRecMuonN(NULL), RecMuonPt(0.), RecMuonEta(0.),
      RecMuonRap(0.), RecMuonPhi(0.), RecMuonThetaAbs(0.),
      RecMuonTriggerMatch(0), RecMuonChiSquare(0.), RecMuonTriggerChiSquare(0.),
      RecMuonIsGoodTrack(0), RecMCMuonPt(0.), RecMCMuonEta(0.),
      RecMCMuonRap(0.), RecMCMuonPhi(0.),

      fTreeMCMuonP(NULL), fTreeMCMuonN(NULL), MCMuonPt(0.), MCMuonEta(0.),
      MCMuonRap(0.), MCMuonPhi(0.), MCMuonDetect(0.),

      fTreeULSDimuon(NULL), fTreeLSppDimuon(NULL), fTreeLSmmDimuon(NULL),
      fTreeMCULSDimuon(NULL), fTreeMCLSppDimuon(NULL), fTreeMCLSmmDimuon(NULL),
      fTreeMixULSDimuon(NULL), fTreeMixLSppDimuon(NULL),
      fTreeMixLSmmDimuon(NULL), RecDimuonPt(0.), RecDimuonRap(0.),
      RecDimuonMass(0.), RecDimuonCent(0.), RecDimuonDS(0.), RecMCDimuonPt(0.),
      RecMCDimuonRap(0.), RecMCDimuonMass(0.), RecMCMom2Body(0.),
      RecMCMomDalitz(0), RecMCMomPdgCode(0.), RecMCMomMass(0), RecMCMomEta(0),
      RecMCMomPt(0), MCDimuonPt(0.), MCDimuonRap(0.), MCDimuonMass(0.),
      MCMom2Body(0.), MCMomDalitz(0), MCMomPt(0), MCMomEta(0), MCMomPdgCode(0.),
      MCDimuonDetected(0),

      fIsMidTrackAna(0),

      fTreeTrack(NULL), fTreeTrackLight(NULL), fTrackPt(0.), fTrackP(0.),
      fTrackTheta(0.), fTrackEta(0.), fTrackPhi(0.), fTrackLength(0.),
      fTrackBeta(0.), fTrackdEdx(0.), fTrackTrackChi2perNDF(0.),
      fTrackTrigChi2perNDF(0.), fTrackTrackITSNcls(0.), fTrackTrackTPCNcls(0.),
      fTrackTrackTOFNcls(0.), fTrackTrackTPCChi2(0.), fTrackTrackITSChi2(0.),
      fTrackTPCCrossedRows(0.), fTrackTPCFindableNcls(0.), fTrackTOFBCTime(0.),
      fTrackTOFKinkIndex(0.), fTrackDCAxy(0.), fTrackDCAz(0.),
      fTrackTPCsigmaKaon(0.), fTrackTOFsigmaKaon(0.), fTrackTPCsigmaPion(0.),
      fTrackTOFsigmaPion(0.), fTrackTPCsigmaProton(0.),
      fTrackTOFsigmaProton(0.), fTrackTPCsigmaElectron(0.),
      fTrackTOFsigmaElectron(0.), fTrackTPCsigmaMuon(0.),
      fTrackTOFsigmaMuon(0.), fTrackThetaAbs(0.), fTrackGlobal(0),
      fTrackGlobalNoDCA(0), fTrackTPConly(0), fTrackGoodFwdQuarity(0),
      fTrackMCType(0), fTrackTrigMatch(0), fTrackPdgCode(0),
      fTrackMotherPdgCode(0) {
  double fCentBins[] = {-1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 101};
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

  fRandom = new TRandom1();
  // fRandom->SetSeed();

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskAODTrackPairMC::~AliAnalysisTaskAODTrackPairMC() {}
//________________________________________________________________________
void AliAnalysisTaskAODTrackPairMC::UserCreateOutputObjects() {
  // Create histograms
  // Called once
  fOutputList = new TList();
  fOutputList->SetOwner(true);

  double bins_cent_hist[] = {0.,  1.,  2.,  3.,  4.,  5.,
                             10., 20., 30., 50., 70., 100.};
  int binnum_cent_hist = sizeof(bins_cent_hist) / sizeof(double) - 1;

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
  fEventCounter = new TH2F("fEventCounter", "", binnum_event_hist,
                           bins_event_hist, binnum_cent_hist, bins_cent_hist);
  for (unsigned int iname = 0;
       iname < sizeof(event_label) / sizeof(std::string); ++iname) {
    fEventCounter->GetXaxis()->SetBinLabel(iname + 1,
                                           event_label[iname].c_str());
  }

  // fOutputList->Add(fEventCounter);

  if (!fIsMixingAnalysis) {
    fTreeULSDimuon = new TTree("fTreeULSDimuon", "");
    fTreeULSDimuon->Branch("RecDimuonPt", &RecDimuonPt, "RecDimuonPt/F");
    fTreeULSDimuon->Branch("RecDimuonRap", &RecDimuonRap, "RecDimuonRap/F");
    fTreeULSDimuon->Branch("RecDimuonMass", &RecDimuonMass, "RecDimuonMass/F");
    fTreeULSDimuon->Branch("RecDimuonCent", &RecDimuonCent, "RecDimuonCent/F");
    fTreeULSDimuon->Branch("RecDimuonDS", &RecDimuonDS, "RecDimuonDS/F");
    fTreeULSDimuon->Branch("RecMCDimuonPt", &RecMCDimuonPt, "RecMCDimuonPt/F");
    fTreeULSDimuon->Branch("RecMCDimuonRap", &RecMCDimuonRap,
                           "RecMCDimuonRap/F");
    fTreeULSDimuon->Branch("RecMCDimuonMass", &RecMCDimuonMass,
                           "RecMCDimuonMass/F");
    fTreeULSDimuon->Branch("RecMCMom2Body", &RecMCMom2Body, "RecMCMom2Body/I");
    fTreeULSDimuon->Branch("RecMCMomDalitz", &RecMCMomDalitz,
                           "RecMCMomDalitz/I");
    fTreeULSDimuon->Branch("RecMCMomPdgCode", &RecMCMomPdgCode,
                           "RecMCMomPdgCode/F");
    fTreeULSDimuon->Branch("RecMCMomMass", &RecMCMomMass, "RecMCMomMass/F");
    fTreeULSDimuon->Branch("RecMCMomEta", &RecMCMomEta, "RecMCMomEta/F");
    fTreeULSDimuon->Branch("RecMCMomPt", &RecMCMomPt, "RecMCMomPt/F");
    fTreeULSDimuon->Branch("RunNumber", &fRunNumber, "fRunNumber/I");
    // fOutputList->Add(fTreeULSDimuon);

    fTreeLSppDimuon = new TTree("fTreeLSppDimuon", "");
    fTreeLSppDimuon->Branch("RecDimuonPt", &RecDimuonPt, "RecDimuonPt/F");
    fTreeLSppDimuon->Branch("RecDimuonRap", &RecDimuonRap, "RecDimuonRap/F");
    fTreeLSppDimuon->Branch("RecDimuonMass", &RecDimuonMass, "RecDimuonMass/F");
    fTreeLSppDimuon->Branch("RecDimuonCent", &RecDimuonCent, "RecDimuonCent/F");
    fTreeLSppDimuon->Branch("RecDimuonDS", &RecDimuonDS, "RecDimuonDS/F");
    fTreeLSppDimuon->Branch("RecMCDimuonPt", &RecMCDimuonPt, "RecMCDimuonPt/F");
    fTreeLSppDimuon->Branch("RecMCDimuonRap", &RecMCDimuonRap,
                            "RecMCDimuonRap/F");
    fTreeLSppDimuon->Branch("RecMCDimuonMass", &RecMCDimuonMass,
                            "RecMCDimuonMass/F");
    fTreeLSppDimuon->Branch("RecMCMomPdgCode", &RecMCMomPdgCode,
                            "RecMCMomPdgCode/F");
    fTreeLSppDimuon->Branch("RecMCMomMass", &RecMCMomMass, "RecMCMomMass/F");
    fTreeLSppDimuon->Branch("RecMCMomEta", &RecMCMomEta, "RecMCMomEta/F");
    fTreeLSppDimuon->Branch("RecMCMomPt", &RecMCMomPt, "RecMCMomPt/F");
    fTreeLSppDimuon->Branch("RunNumber", &fRunNumber, "fRunNumber/I");
    // fOutputList->Add(fTreeLSppDimuon);

    fTreeLSmmDimuon = new TTree("fTreeLSmmDimuon", "");
    fTreeLSmmDimuon->Branch("RecDimuonPt", &RecDimuonPt, "RecDimuonPt/F");
    fTreeLSmmDimuon->Branch("RecDimuonRap", &RecDimuonRap, "RecDimuonRap/F");
    fTreeLSmmDimuon->Branch("RecDimuonMass", &RecDimuonMass, "RecDimuonMass/F");
    fTreeLSmmDimuon->Branch("RecDimuonCent", &RecDimuonCent, "RecDimuonCent/F");
    fTreeLSmmDimuon->Branch("RecDimuonDS", &RecDimuonDS, "RecDimuonDS/F");
    fTreeLSmmDimuon->Branch("RecMCDimuonPt", &RecMCDimuonPt, "RecMCDimuonPt/F");
    fTreeLSmmDimuon->Branch("RecMCDimuonRap", &RecMCDimuonRap,
                            "RecMCDimuonRap/F");
    fTreeLSmmDimuon->Branch("RecMCDimuonMass", &RecMCDimuonMass,
                            "RecMCDimuonMass/F");
    fTreeLSmmDimuon->Branch("RecMCMomPdgCode", &RecMCMomPdgCode,
                            "RecMCMomPdgCode/F");
    fTreeLSmmDimuon->Branch("RecMCMomMass", &RecMCMomMass, "RecMCMomMass/F");
    fTreeLSmmDimuon->Branch("RecMCMomEta", &RecMCMomEta, "RecMCMomEta/F");
    fTreeLSmmDimuon->Branch("RecMCMomPt", &RecMCMomPt, "RecMCMomPt/F");
    fTreeLSmmDimuon->Branch("RunNumber", &fRunNumber, "fRunNumber/I");
    // fOutputList->Add(fTreeLSmmDimuon);
  } else {
    fTreeMixULSDimuon = new TTree("fTreeMixULSDimuon", "");
    fTreeMixULSDimuon->Branch("RecDimuonPt", &RecDimuonPt, "RecDimuonPt/F");
    fTreeMixULSDimuon->Branch("RecDimuonRap", &RecDimuonRap, "RecDimuonRap/F");
    fTreeMixULSDimuon->Branch("RecDimuonMass", &RecDimuonMass,
                              "RecDimuonMass/F");
    fTreeMixULSDimuon->Branch("RecDimuonCent", &RecDimuonCent,
                              "RecDimuonCent/F");
    fTreeMixULSDimuon->Branch("RecDimuonDS", &RecDimuonDS, "RecDimuonDS/F");
    // fOutputList->Add(fTreeMixULSDimuon);

    fTreeMixLSppDimuon = new TTree("fTreeMixLSppDimuon", "");
    fTreeMixLSppDimuon->Branch("RecDimuonPt", &RecDimuonPt, "RecDimuonPt/F");
    fTreeMixLSppDimuon->Branch("RecDimuonRap", &RecDimuonRap, "RecDimuonRap/F");
    fTreeMixLSppDimuon->Branch("RecDimuonMass", &RecDimuonMass,
                               "RecDimuonMass/F");
    fTreeMixLSppDimuon->Branch("RecDimuonCent", &RecDimuonCent,
                               "RecDimuonCent/F");
    fTreeMixLSppDimuon->Branch("RecDimuonDS", &RecDimuonDS, "RecDimuonDS/F");
    // fOutputList->Add(fTreeMixLSppDimuon);

    fTreeMixLSmmDimuon = new TTree("fTreeMixLSmmDimuon", "");
    fTreeMixLSmmDimuon->Branch("RecDimuonPt", &RecDimuonPt, "RecDimuonPt/F");
    fTreeMixLSmmDimuon->Branch("RecDimuonRap", &RecDimuonRap, "RecDimuonRap/F");
    fTreeMixLSmmDimuon->Branch("RecDimuonMass", &RecDimuonMass,
                               "RecDimuonMass/F");
    fTreeMixLSmmDimuon->Branch("RecDimuonCent", &RecDimuonCent,
                               "RecDimuonCent/F");
    fTreeMixLSmmDimuon->Branch("RecDimuonDS", &RecDimuonDS, "RecDimuonDS/F");
    // fOutputList->Add(fTreeMixLSmmDimuon);
  }

  fTreeMCULSDimuon = new TTree("fTreeMCULSDimuon", "");
  fTreeMCULSDimuon->Branch("MCDimuonPt", &MCDimuonPt, "MCDimuonPt/F");
  fTreeMCULSDimuon->Branch("MCDimuonRap", &MCDimuonRap, "MCDimuonRap/F");
  fTreeMCULSDimuon->Branch("MCDimuonMass", &MCDimuonMass, "MCDimuonMass/F");
  fTreeMCULSDimuon->Branch("MCDimuonDetected", &MCDimuonDetected,
                           "MCDimuonDetected/I");
  fTreeMCULSDimuon->Branch("RecDimuonPt", &RecDimuonPt, "RecDimuonPt/F");
  fTreeMCULSDimuon->Branch("RecDimuonRap", &RecDimuonRap, "RecDimuonRap/F");
  fTreeMCULSDimuon->Branch("RecDimuonMass", &RecDimuonMass, "RecDimuonMass/F");
  fTreeMCULSDimuon->Branch("MCMom2Body", &MCMom2Body, "MCMom2Body/I");
  fTreeMCULSDimuon->Branch("MCMomDalitz", &MCMomDalitz, "MCMomDalitz/I");
  fTreeMCULSDimuon->Branch("MCMomPt", &MCMomPt, "MCMomPt/F");
  fTreeMCULSDimuon->Branch("MCMomEta", &MCMomEta, "MCMomEta/F");
  fTreeMCULSDimuon->Branch("MCMomPdgCode", &MCMomPdgCode, "MCMomPdgCode/F");
  fTreeMCULSDimuon->Branch("RunNumber", &fRunNumber, "fRunNumber/I");
  // fOutputList->Add(fTreeMCULSDimuon);

  fTreeMCLSppDimuon = new TTree("fTreeMCLSppDimuon", "");
  fTreeMCLSppDimuon->Branch("MCDimuonPt", &MCDimuonPt, "MCDimuonPt/F");
  fTreeMCLSppDimuon->Branch("MCDimuonRap", &MCDimuonRap, "MCDimuonRap/F");
  fTreeMCLSppDimuon->Branch("MCDimuonMass", &MCDimuonMass, "MCDimuonMass/F");
  fTreeMCLSppDimuon->Branch("MCMom2Body", &MCMom2Body, "MCMom2Body/I");
  fTreeMCLSppDimuon->Branch("MCMomDalitz", &MCMomDalitz, "MCMomDalitz/I");
  fTreeMCLSppDimuon->Branch("MCMomPt", &MCMomPt, "MCMomPt/F");
  fTreeMCLSppDimuon->Branch("MCMomEta", &MCMomEta, "MCMomEta/F");
  fTreeMCLSppDimuon->Branch("MCMomPdgCode", &MCMomPdgCode, "MCMomPdgCode/F");
  fTreeMCLSppDimuon->Branch("MCDimuonDetected", &MCDimuonDetected,
                            "MCDimuonDetected/I");
  fTreeMCLSppDimuon->Branch("RunNumber", &fRunNumber, "fRunNumber/I");
  // fOutputList->Add(fTreeMCLSppDimuon);

  fTreeMCLSmmDimuon = new TTree("fTreeMCLSmmDimuon", "");
  fTreeMCLSmmDimuon->Branch("MCDimuonPt", &MCDimuonPt, "MCDimuonPt/F");
  fTreeMCLSmmDimuon->Branch("MCDimuonRap", &MCDimuonRap, "MCDimuonRap/F");
  fTreeMCLSmmDimuon->Branch("MCDimuonMass", &MCDimuonMass, "MCDimuonMass/F");
  fTreeMCLSmmDimuon->Branch("MCMom2Body", &MCMom2Body, "MCMom2Body/I");
  fTreeMCLSmmDimuon->Branch("MCMomDalitz", &MCMomDalitz, "MCMomDalitz/I");
  fTreeMCLSmmDimuon->Branch("MCMomPt", &MCMomPt, "MCMomPt/F");
  fTreeMCLSmmDimuon->Branch("MCMomEta", &MCMomEta, "MCMomEta/F");
  fTreeMCLSmmDimuon->Branch("MCMomPdgCode", &MCMomPdgCode, "MCMomPdgCode/F");
  fTreeMCLSmmDimuon->Branch("MCDimuonDetected", &MCDimuonDetected,
                            "MCDimuonDetected/I");
  fTreeMCLSmmDimuon->Branch("RunNumber", &fRunNumber, "fRunNumber/I");
  // fOutputList->Add(fTreeMCLSmmDimuon);

  fTreeRecMuonP = new TTree("fTreeRecMuonP", "");
  fTreeRecMuonP->Branch("RecMuonPt", &RecMuonPt, "RecMuonPt/F");
  fTreeRecMuonP->Branch("RecMuonEta", &RecMuonEta, "RecMuonEta/F");
  fTreeRecMuonP->Branch("RecMuonRap", &RecMuonRap, "RecMuonRap/F");
  fTreeRecMuonP->Branch("RecMuonPhi", &RecMuonPhi, "RecMuonPhi/F");
  fTreeRecMuonP->Branch("RecMCMuonPt", &RecMCMuonPt, "RecMCMuonPt/F");
  fTreeRecMuonP->Branch("RecMCMuonEta", &RecMCMuonEta, "RecMCMuonEta/F");
  fTreeRecMuonP->Branch("RecMCMuonRap", &RecMCMuonRap, "RecMCMuonRap/F");
  fTreeRecMuonP->Branch("RecMCMuonPhi", &RecMCMuonPhi, "RecMCMuonPhi/F");
  fTreeRecMuonP->Branch("RecMuonThetaAbs", &RecMuonThetaAbs,
                        "RecMuonThetaAbs/F");
  fTreeRecMuonP->Branch("RecMuonTriggerMatch", &RecMuonTriggerMatch,
                        "RecMuonTriggerMatch/I");
  fTreeRecMuonP->Branch("RecMuonChiSquare", &RecMuonChiSquare,
                        "RecMuonChiSquare/F");
  fTreeRecMuonP->Branch("RecMuonTriggerChiSquare", &RecMuonTriggerChiSquare,
                        "RecMuonTriggerChiSquare/F");
  fTreeRecMuonP->Branch("RecMuonIsGoodTrack", &RecMuonIsGoodTrack,
                        "RecMuonIsGoodTrack/I");
  // fOutputList->Add(fTreeRecMuonP);

  fTreeRecMuonN = new TTree("fTreeRecMuonN", "");
  fTreeRecMuonN->Branch("RecMuonPt", &RecMuonPt, "RecMuonPt/F");
  fTreeRecMuonN->Branch("RecMuonEta", &RecMuonEta, "RecMuonEta/F");
  fTreeRecMuonN->Branch("RecMuonRap", &RecMuonRap, "RecMuonRap/F");
  fTreeRecMuonN->Branch("RecMuonPhi", &RecMuonPhi, "RecMuonPhi/F");
  fTreeRecMuonN->Branch("RecMCMuonPt", &RecMCMuonPt, "RecMCMuonPt/F");
  fTreeRecMuonN->Branch("RecMCMuonEta", &RecMCMuonEta, "RecMCMuonEta/F");
  fTreeRecMuonN->Branch("RecMCMuonRap", &RecMCMuonRap, "RecMCMuonRap/F");
  fTreeRecMuonN->Branch("RecMCMuonPhi", &RecMCMuonPhi, "RecMCMuonPhi/F");
  fTreeRecMuonN->Branch("RecMuonThetaAbs", &RecMuonThetaAbs,
                        "RecMuonThetaAbs/F");
  fTreeRecMuonN->Branch("RecMuonTriggerMatch", &RecMuonTriggerMatch,
                        "RecMuonTriggerMatch/I");
  fTreeRecMuonN->Branch("RecMuonChiSquare", &RecMuonChiSquare,
                        "RecMuonChiSquare/F");
  fTreeRecMuonN->Branch("RecMuonTriggerChiSquare", &RecMuonTriggerChiSquare,
                        "RecMuonTriggerChiSquare/F");
  fTreeRecMuonN->Branch("RecMuonIsGoodTrack", &RecMuonIsGoodTrack,
                        "RecMuonIsGoodTrack/I");
  // fOutputList->Add(fTreeRecMuonN);

  fTreeMCMuonP = new TTree("fTreeMCMuonP", "");
  fTreeMCMuonP->Branch("MCMuonPt", &MCMuonPt, "MCMuonPt/F");
  fTreeMCMuonP->Branch("MCMuonEta", &MCMuonEta, "MCMuonEta/F");
  fTreeMCMuonP->Branch("MCMuonRap", &MCMuonRap, "MCMuonRap/F");
  fTreeMCMuonP->Branch("MCMuonPhi", &MCMuonPhi, "MCMuonPhi/F");
  fTreeMCMuonP->Branch("MCMuonDetect", &MCMuonDetect, "MCMuonDetect/I");
  // fOutputList->Add(fTreeMCMuonP);

  fTreeMCMuonN = new TTree("fTreeMCMuonN", "");
  fTreeMCMuonN->Branch("MCMuonPt", &MCMuonPt, "MCMuonPt/F");
  fTreeMCMuonN->Branch("MCMuonEta", &MCMuonEta, "MCMuonEta/F");
  fTreeMCMuonN->Branch("MCMuonRap", &MCMuonRap, "MCMuonRap/F");
  fTreeMCMuonN->Branch("MCMuonPhi", &MCMuonPhi, "MCMuonPhi/F");
  fTreeMCMuonN->Branch("MCMuonDetect", &MCMuonDetect, "MCMuonDetect/I");
  // fOutputList->Add(fTreeMCMuonN);

  fHistEventVtxZ = new TH1F("fHistEventVtxZ", "", 60, -30, 30);
  fHistEventCent = new TH1F("fHistEventCent", "", 100, 0, 100);
  fHistEventMulti = new TH1F("fHistEventMulti", "", 200, 0, 200);
  fHistEventVtxCont = new TH1F("fHistEventVtxCont", "", 100, 0, 100);
  fOutputList->Add(fHistEventVtxZ);
  fOutputList->Add(fHistEventCent);
  fOutputList->Add(fHistEventMulti);
  fOutputList->Add(fHistEventVtxCont);

  fHistTrackEta = new TH2F("fHistTrackEta", "", 20, 0, 10, 25, -4.5, -2.0);
  fHistTrackThetaAbs = new TH2F("fHistTrackThetaAbs", "", 20, 0, 10, 60, 0, 15);
  fHistTrackTriggerMatch =
      new TH2F("fHistTrackTriggerMatch", "", 20, 0, 10, 5, 0, 5);
  fHistTrackPDCA = new TH2F("fHistTrackPDCA", "", 20, 0, 10, 200, 0, 20);
  fHistTrackChiSquare =
      new TH2F("fHistTrackChiSquare", "", 20, 0, 10, 100, 0, 10);
  fHistTriggerChiSquare =
      new TH2F("fHistTriggerChiSquare", "", 20, 0, 10, 100, 0, 10);
  /*
  fOutputList->Add(fHistTrackEta);
  fOutputList->Add(fHistTrackThetaAbs);
  fOutputList->Add(fHistTrackTriggerMatch);
  fOutputList->Add(fHistTrackPDCA);
  fOutputList->Add(fHistTrackChiSquare);
  fOutputList->Add(fHistTriggerChiSquare);
  */

  fTreeTrack = new TTree("fTreeTrack", "Tree for machine leraning for PID");
  fTreeTrackLight =
      new TTree("fTreeTrackLight", "Tree for machine leraning for PID");

  if (fIsMidTrackAna) {
    fTreeTrack->Branch("fTrackPt", &fTrackPt, "fTrackPt/F");
    fTreeTrack->Branch("fTrackP", &fTrackP, "fTrackP/F");
    fTreeTrack->Branch("fTrackEta", &fTrackEta, "fTrackEta/F");
    fTreeTrack->Branch("fTrackPhi", &fTrackPhi, "fTrackPhi/F");
    fTreeTrack->Branch("fTrackLength", &fTrackLength, "fTrackLength/F");
    fTreeTrack->Branch("fTrackBeta", &fTrackBeta, "fTrackBeta/F");
    fTreeTrack->Branch("fTrackdEdx", &fTrackdEdx, "fTrackdEdx/F");
    fTreeTrack->Branch("fTrackTrackChi2perNDF", &fTrackTrackChi2perNDF,
                       "fTrackTrackChi2perNDF/F");
    fTreeTrack->Branch("fTrackTrackITSNcls", &fTrackTrackITSNcls,
                       "fTrackTrackITSNcls/F");
    fTreeTrack->Branch("fTrackTrackTPCNcls", &fTrackTrackTPCNcls,
                       "fTrackTrackTPCNcls/F");
    fTreeTrack->Branch("fTrackTrackTOFNcls", &fTrackTrackTOFNcls,
                       "fTrackTrackTOFNcls/F");
    fTreeTrack->Branch("fTrackTrackTPCChi2", &fTrackTrackTPCChi2,
                       "fTrackTrackTPCChi2/F");
    fTreeTrack->Branch("fTrackTrackITSChi2", &fTrackTrackITSChi2,
                       "fTrackTrackITSChi2/F");
    fTreeTrack->Branch("fTrackTPCCrossedRows", &fTrackTPCCrossedRows,
                       "fTrackTPCCrossedRows/F");
    fTreeTrack->Branch("fTrackTPCFindableNcls", &fTrackTPCFindableNcls,
                       "fTrackTPCFindableNcls/F");
    fTreeTrack->Branch("fTrackTOFBCTime", &fTrackTOFBCTime,
                       "fTrackTOFBCTime/F");
    fTreeTrack->Branch("fTrackDCAxy", &fTrackDCAxy, "fTrackDCAxy/F");
    fTreeTrack->Branch("fTrackDCAz", &fTrackDCAz, "fTrackDCAz/F");
    fTreeTrack->Branch("fTrackTPCsigmaPion", &fTrackTPCsigmaPion,
                       "fTrackTPCsigmaPion/F");
    fTreeTrack->Branch("fTrackTOFsigmaPion", &fTrackTOFsigmaPion,
                       "fTrackTOFsigmaPion/F");
    fTreeTrack->Branch("fTrackTPCsigmaKaon", &fTrackTPCsigmaKaon,
                       "fTrackTPCsigmaKaon/F");
    fTreeTrack->Branch("fTrackTOFsigmaKaon", &fTrackTOFsigmaKaon,
                       "fTrackTOFsigmaKaon/F");
    fTreeTrack->Branch("fTrackTPCsigmaProton", &fTrackTPCsigmaProton,
                       "fTrackTPCsigmaProton/F");
    fTreeTrack->Branch("fTrackTOFsigmaProton", &fTrackTOFsigmaProton,
                       "fTrackTOFsigmaProton/F");
    fTreeTrack->Branch("fTrackTPCsigmaElectron", &fTrackTPCsigmaElectron,
                       "fTrackTPCsigmaElectron/F");
    fTreeTrack->Branch("fTrackTOFsigmaElectron", &fTrackTOFsigmaElectron,
                       "fTrackTOFsigmaElectron/F");
    fTreeTrack->Branch("fTrackTPCsigmaMuon", &fTrackTPCsigmaMuon,
                       "fTrackTPCsigmaMuon/F");
    fTreeTrack->Branch("fTrackTOFsigmaMuon", &fTrackTOFsigmaMuon,
                       "fTrackTOFsigmaMuon/F");
    fTreeTrack->Branch("fTrackGlobal", &fTrackGlobal, "fTrackGlobal/O");
    fTreeTrack->Branch("fTrackGlobalNoDCA", &fTrackGlobalNoDCA,
                       "fTrackGlobalNoDCA/O");
    fTreeTrack->Branch("fTrackTPConly", &fTrackTPConly, "fTrackTPConly/O");
    fTreeTrack->Branch("fTrackPdgCode", &fTrackPdgCode, "fTrackPdgCode/I");
    fTreeTrack->Branch("fTrackMotherPdgCode", &fTrackMotherPdgCode,
                       "fTrackMotherPdgCode/I");
    fTreeTrack->Branch("fTrackMCType", &fTrackMCType, "fTrackMCType/I");
    fOutputList->Add(fTreeTrack);

    fTreeTrackLight->Branch("fTrackPt", &fTrackPt, "fTrackPt/F");
    fTreeTrackLight->Branch("fTrackEta", &fTrackEta, "fTrackEta/F");
    fTreeTrackLight->Branch("fTrackPhi", &fTrackPhi, "fTrackPhi/F");
    fTreeTrackLight->Branch("fTrackBeta", &fTrackBeta, "fTrackBeta/F");
    fTreeTrackLight->Branch("fTrackdEdx", &fTrackdEdx, "fTrackdEdx/F");
    fTreeTrackLight->Branch("fTrackDCAxy", &fTrackDCAxy, "fTrackDCAxy/F");
    fTreeTrackLight->Branch("fTrackDCAz", &fTrackDCAz, "fTrackDCAz/F");
    fTreeTrackLight->Branch("fTrackTPCsigmaPion", &fTrackTPCsigmaPion,
                            "fTrackTPCsigmaPion/F");
    fTreeTrackLight->Branch("fTrackTOFsigmaPion", &fTrackTOFsigmaPion,
                            "fTrackTOFsigmaPion/F");
    fTreeTrackLight->Branch("fTrackTPCsigmaKaon", &fTrackTPCsigmaKaon,
                            "fTrackTPCsigmaKaon/F");
    fTreeTrackLight->Branch("fTrackTOFsigmaKaon", &fTrackTOFsigmaKaon,
                            "fTrackTOFsigmaKaon/F");
    fTreeTrackLight->Branch("fTrackTPCsigmaProton", &fTrackTPCsigmaProton,
                            "fTrackTPCsigmaProton/F");
    fTreeTrackLight->Branch("fTrackTOFsigmaProton", &fTrackTOFsigmaProton,
                            "fTrackTOFsigmaProton/F");
    fTreeTrackLight->Branch("fTrackTPCsigmaElectron", &fTrackTPCsigmaElectron,
                            "fTrackTPCsigmaElectron/F");
    fTreeTrackLight->Branch("fTrackTOFsigmaElectron", &fTrackTOFsigmaElectron,
                            "fTrackTOFsigmaElectron/F");
    fTreeTrackLight->Branch("fTrackTPCsigmaMuon", &fTrackTPCsigmaMuon,
                            "fTrackTPCsigmaMuon/F");
    fTreeTrackLight->Branch("fTrackTOFsigmaMuon", &fTrackTOFsigmaMuon,
                            "fTrackTOFsigmaMuon/F");
    fTreeTrackLight->Branch("fTrackGlobal", &fTrackGlobal, "fTrackGlobal/O");
    fTreeTrackLight->Branch("fTrackGlobalNoDCA", &fTrackGlobalNoDCA,
                            "fTrackGlobalNoDCA/O");
    fTreeTrackLight->Branch("fTrackTPConly", &fTrackTPConly, "fTrackTPConly/O");
    fTreeTrackLight->Branch("fTrackMCType", &fTrackMCType, "fTrackMCType/I");
    fTreeTrackLight->Branch("fTrackPdgCode", &fTrackPdgCode, "fTrackPdgCode/I");
    fTreeTrackLight->Branch("fTrackMotherPdgCode", &fTrackMotherPdgCode,
                            "fTrackMotherPdgCode/I");
    fOutputList->Add(fTreeTrackLight);
  } else {
    fTreeTrack->Branch("fTrackPt", &fTrackPt, "fTrackPt/F");
    fTreeTrack->Branch("fTrackP", &fTrackP, "fTrackP/F");
    fTreeTrack->Branch("fTrackEta", &fTrackEta, "fTrackEta/F");
    fTreeTrack->Branch("fTrackPhi", &fTrackPhi, "fTrackPhi/F");
    fTreeTrack->Branch("fTrackTrackChi2perNDF", &fTrackTrackChi2perNDF,
                       "fTrackTrackChi2perNDF/F");
    fTreeTrack->Branch("fTrackTrigChi2perNDF", &fTrackTrigChi2perNDF,
                       "fTrackTrigChi2perNDF/F");
    fTreeTrack->Branch("fTrackThetaAbs", &fTrackThetaAbs, "fTrackThetaAbs/F");
    fTreeTrack->Branch("fTrackDCAxy", &fTrackDCAxy, "fTrackDCAxy/F");
    fTreeTrack->Branch("fTrackDCAz", &fTrackDCAz, "fTrackDCAz/F");
    fTreeTrack->Branch("fTrackTrigMatch", &fTrackTrigMatch,
                       "fTrackTrigMatch/I");
    fTreeTrack->Branch("fTrackGoodFwdQuarity", &fTrackGoodFwdQuarity,
                       "fTrackGoodFwdQuarity/O");
    fTreeTrack->Branch("fTrackMCType", &fTrackMCType, "fTrackMCType/I");
    fTreeTrack->Branch("fTrackPdgCode", &fTrackPdgCode, "fTrackPdgCode/I");
    fTreeTrack->Branch("fTrackMotherPdgCode", &fTrackMotherPdgCode,
                       "fTrackMotherPdgCode/I");
    fOutputList->Add(fTreeTrack);
  }

  PostData(1, fOutputList);
}

//________________________________________________________________________

void AliAnalysisTaskAODTrackPairMC::UserExec(Option_t *) {

  if (!fIsMC) {
    return;
  }

  if (!Initialize()) {
    return;
  }

  if (!fUtils->isAcceptEvent()) {
    return;
  }

  // processMC();

  if (fIsMidTrackAna) {
    MidTrackQA();
  } else {
    FwdTrackQA();
  }

  if (!fIsMixingAnalysis) {
    // FwdMuonPairAnalysis();
  }
}

bool AliAnalysisTaskAODTrackPairMC::Initialize() {
  fEvent = dynamic_cast<AliAODEvent *>(InputEvent());

  if (!fUtils->setEvent(fEvent, fInputHandler)) {
    return false;
  }

  if (fRunNumber != fEvent->GetRunNumber()) {
    fRunNumber = fUtils->getRunnumber();
    AliMuonTrackCuts *trackCut = fUtils->getMuonTrackCuts();
    trackCut->SetRun(fInputHandler);
  }

  if (fIsMC) {

    fMCTrackArray = dynamic_cast<TClonesArray *>(
        fEvent->FindListObject(AliAODMCParticle::StdBranchName()));

    if (!fMCTrackArray) {
      return false;
    }

    fUtils->setMCArray(fMCTrackArray);

    AliAODMCHeader *mcHeader = (AliAODMCHeader *)fEvent->GetList()->FindObject(
        AliAODMCHeader::StdBranchName());
    if (!mcHeader)
      return false;

    TList *headerList = mcHeader->GetCocktailHeaders();

    if (!headerList)
      return false;

    for (Int_t i = 0; i < headerList->GetEntries(); i++) {
      AliGenEventHeader *eventHeader2 = (AliGenEventHeader *)headerList->At(i);
      TString name = eventHeader2->GetName();

      if (name.Contains("EtaDalitz")) {
        fUtils->setDalitzProd(true);
      } else if (name.Contains("EtaDirect")) {
        fUtils->set2BodyProd(true);
      } else if (name.Contains("OmegaDalitz")) {
        fUtils->setDalitzProd(true);
      } else if (name.Contains("OmegaDirect")) {
        fUtils->set2BodyProd(true);
      } else if (name.Contains("EtaPrimeDalitz")) {
        fUtils->setDalitzProd(true);
      } else if (name.Contains("PhiDirect")) {
        fUtils->set2BodyProd(true);
      } else if (name.Contains("RhoDirect")) {
        fUtils->set2BodyProd(true);
      } else {
        fUtils->setDalitzProd(false);
        fUtils->set2BodyProd(false);
      }
    }
  }

  fUtils->getTriggerInfo(fIsCINT7, fIsCMSL7, fIsCMSH7, fIsCMUL7, fIsCMLL7);

  return true;
}

bool AliAnalysisTaskAODTrackPairMC::EventQA() {
  fHistEventVtxZ->Fill(fUtils->getVtxZ());
  fHistEventCent->Fill(fUtils->getCentClass());
  fHistEventMulti->Fill(fUtils->getNSPDTrkInfo(4));
  fHistEventVtxCont->Fill(fUtils->getVtxCont());
  return true;
}

bool AliAnalysisTaskAODTrackPairMC::FwdMuonTrackQA(AliAODTrack *track) {
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

bool AliAnalysisTaskAODTrackPairMC::MidTrackQA() {

  AliAODMCParticle *particle;
  AliAODTrack *track;

  int nTrack = fEvent->GetNumberOfTracks();

  for (int iTrack = 0; iTrack < nTrack; ++iTrack) {

    track = (AliAODTrack *)fEvent->GetTrack(iTrack);

    if (track->TestFilterBit(AliAODTrack::kTrkGlobal) == false &&
        track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA) == false &&
        track->TestFilterBit(AliAODTrack::kTrkTPCOnly) == false) {
      continue;
    }

    if (fabs(track->Eta()) > 0.8 || track->Pt() < 0.15) {
      continue;
    }

    float sigTOF = track->GetTOFsignal();

    if (sigTOF > 99998.5) {
      continue;
    }

    float length = track->GetIntegratedLength();
    float beta =
        (sigTOF > 0) ? (double)length / (2.99792457999999984e-02 * sigTOF) : 0;

    float dca_xy = 9999;
    float dca_z = 9999;
    track->GetImpactParameters(dca_xy, dca_z);

    if (fIsMC && fMCTrackArray && track->GetLabel() > 0) {
      particle = (AliAODMCParticle *)fMCTrackArray->At(track->GetLabel());
      if (particle) {
        fTrackPdgCode = particle->GetPdgCode();
      } else {
        fTrackPdgCode = 0;
      }
    } else {
      fTrackPdgCode = 0;
    }

    if (fabs(fTrackPdgCode) == 211) {
      if (fRandom->Rndm() > 0.001) {
        continue;
      }
    } else if (fabs(fTrackPdgCode) == 321) {
      if (fRandom->Rndm() > 0.01) {
        continue;
      }
    } else if (fabs(fTrackPdgCode) == 2212) {
      if (fRandom->Rndm() > 0.05) {
        continue;
      }
    }

    AliTOFHeader *tofHeader = (AliTOFHeader *)track->GetTOFHeader();
    fTrackPt = track->Pt();
    fTrackP = track->P();
    fTrackEta = track->Eta();
    fTrackPhi = track->Phi();
    fTrackLength = track->GetIntegratedLength();
    fTrackBeta = beta;
    fTrackdEdx = track->GetTPCsignal();
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
    fTrackTPCsigmaPion = fUtils->getTPCSigma(track, AliPID::kPion);
    fTrackTOFsigmaPion = fUtils->getTOFSigma(track, AliPID::kPion);
    fTrackTPCsigmaKaon = fUtils->getTPCSigma(track, AliPID::kKaon);
    fTrackTOFsigmaKaon = fUtils->getTOFSigma(track, AliPID::kKaon);
    fTrackTPCsigmaProton = fUtils->getTPCSigma(track, AliPID::kProton);
    fTrackTOFsigmaProton = fUtils->getTOFSigma(track, AliPID::kProton);
    fTrackTPCsigmaElectron = fUtils->getTPCSigma(track, AliPID::kElectron);
    fTrackTOFsigmaElectron = fUtils->getTOFSigma(track, AliPID::kElectron);
    fTrackTPCsigmaMuon = fUtils->getTPCSigma(track, AliPID::kMuon);
    fTrackTOFsigmaMuon = fUtils->getTOFSigma(track, AliPID::kMuon);
    fTrackGlobal = track->TestFilterBit(AliAODTrack::kTrkGlobal) ? 1 : 0;
    fTrackGlobalNoDCA =
        track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA) ? 1 : 0;
    fTrackTPConly = track->TestFilterBit(AliAODTrack::kTrkTPCOnly) ? 1 : 0;

    if (fIsMC && fMCTrackArray && track->GetLabel() > 0) {
      particle = (AliAODMCParticle *)fMCTrackArray->At(track->GetLabel());
      if (particle) {
        fTrackPdgCode = particle->GetPdgCode();
        if (particle->IsPrimary()) {
          fTrackMCType = 1;
        } else if (particle->IsSecondaryFromWeakDecay()) {
          fTrackMCType = 2;
        } else if (particle->IsSecondaryFromMaterial()) {
          fTrackMCType = 3;
        } else {
          fTrackMCType = 0;
        }
        if (particle->GetMother() > 0) {
          particle =
              (AliAODMCParticle *)fMCTrackArray->At(particle->GetMother());
          if (particle) {
            fTrackMotherPdgCode = particle->GetPdgCode();
          } else {
            fTrackMotherPdgCode = 0;
          }
        } else {
          fTrackMotherPdgCode = 0;
        }
      } else {
        fTrackMCType = 0;
        fTrackPdgCode = 0;
        fTrackMotherPdgCode = 0;
      }
    } else {
      fTrackMCType = 0;
      fTrackPdgCode = 0;
      fTrackMotherPdgCode = 0;
    }

    fTreeTrack->Fill();
    fTreeTrackLight->Fill();
  }

  return true;
}

bool AliAnalysisTaskAODTrackPairMC::FwdTrackQA() {

  AliAODMCParticle *particle;
  AliAODTrack *track;

  int nTrack = fEvent->GetNumberOfTracks();

  for (int iTrack = 0; iTrack < nTrack; ++iTrack) {

    track = (AliAODTrack *)fEvent->GetTrack(iTrack);

    if (!AliAnalysisMuonUtility::IsMuonTrack(track)) {
      continue;
    }

    if (fIsMC && fMCTrackArray && track->GetLabel() > 0) {
      particle = (AliAODMCParticle *)fMCTrackArray->At(track->GetLabel());
      if (particle) {
        fTrackPdgCode = particle->GetPdgCode();
      } else {
        fTrackPdgCode = 0;
      }
    } else {
      fTrackPdgCode = 0;
    }

    fTrackPt = track->Pt();
    fTrackP = track->P();
    fTrackEta = track->Eta();
    fTrackPhi = track->Phi();
    fTrackThetaAbs = AliAnalysisMuonUtility::GetThetaAbsDeg(track);
    fTrackTrigMatch = AliAnalysisMuonUtility::GetMatchTrigger(track);
    fTrackTrackChi2perNDF = AliAnalysisMuonUtility::GetChi2perNDFtracker(track);
    fTrackTrigChi2perNDF = AliAnalysisMuonUtility::GetChi2MatchTrigger(track);
    fTrackGoodFwdQuarity = fUtils->isAcceptFwdMuonTrack(track) ? 1 : 0;
    fTrackDCAxy =
        pow(AliAnalysisMuonUtility::GetXatVertex(track), 2) +
                    pow(AliAnalysisMuonUtility::GetYatVertex(track), 2) >
                0
            ? sqrt(pow(AliAnalysisMuonUtility::GetXatVertex(track), 2) +
                   pow(AliAnalysisMuonUtility::GetYatVertex(track), 2))
            : 0;
    fTrackDCAz = AliAnalysisMuonUtility::GetZatVertex(track);

    if (fIsMC && fMCTrackArray && track->GetLabel() > 0) {
      particle = (AliAODMCParticle *)fMCTrackArray->At(track->GetLabel());
      if (particle) {
        fTrackPdgCode = particle->GetPdgCode();
        int mom = particle->GetMother();
        if (mom > 0) {
          particle = (AliAODMCParticle *)fMCTrackArray->At(mom);
          if (particle) {
            fTrackMotherPdgCode = particle->GetPdgCode();
          } else {
            fTrackMotherPdgCode = 0;
          }
        } else {
          fTrackMotherPdgCode = 0;
        }
      } else {
        fTrackPdgCode = 0;
        fTrackMotherPdgCode = 0;
      }
    } else {
      fTrackPdgCode = 0;
      fTrackMotherPdgCode = 0;
    }

    fTreeTrack->Fill();
    fTreeTrackLight->Fill();
  }

  return true;
}

bool AliAnalysisTaskAODTrackPairMC::FillingRecMuonTree(AliAODTrack *track) {

  int label = track->GetLabel();

  if (label < 0) {
    return false;
  }

  AliAODMCParticle *particle1 = (AliAODMCParticle *)fMCTrackArray->At(label);

  if (!particle1) {
    return false;
  }

  RecMCMuonPt = particle1->Pt();
  RecMCMuonEta = particle1->Eta();
  RecMCMuonRap = particle1->Y();
  RecMCMuonPhi = particle1->Phi();
  RecMuonPt = track->Pt();
  RecMuonEta = track->Eta();
  RecMuonRap = track->Y();
  RecMuonPhi = track->Phi();
  RecMuonThetaAbs = AliAnalysisMuonUtility::GetThetaAbsDeg(track);
  RecMuonTriggerMatch = AliAnalysisMuonUtility::GetMatchTrigger(track);
  RecMuonChiSquare = AliAnalysisMuonUtility::GetChi2perNDFtracker(track);
  RecMuonTriggerChiSquare = AliAnalysisMuonUtility::GetChi2MatchTrigger(track);

  if (fUtils->isAcceptFwdMuonTrack(track)) {
    RecMuonIsGoodTrack = 1;
  } else {
    RecMuonIsGoodTrack = 0;
  }

  if (particle1->Charge() > 0) {
    fTreeRecMuonP->Fill();
  } else {
    fTreeRecMuonN->Fill();
  }

  return true;
}

bool AliAnalysisTaskAODTrackPairMC::FwdMuonPairAnalysisEveMixing() {

  if (!(fInputHandler->IsEventSelected() & fTriggerMaskForMixing))
    return false;

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

          double fill[] = {dimuon->M(), fabs(dimuon->Y()), dimuon->Pt(),
                           fUtils->getCentClass()};

          RecDimuonPt = dimuon->Pt();
          RecDimuonRap = fabs(dimuon->Y());
          RecDimuonMass = dimuon->M();
          RecDimuonCent = fUtils->getCentClass();
          RecDimuonDS = fUtils->getDS();

          string fFiredTrigName = string(fEvent->GetFiredTriggerClasses());

          if (dimuon->Charge() == 0) {
            // fSparseULSDimuon->Fill(fill,(double)1./fUtils->getDS());
            fTreeMixULSDimuon->Fill();
          } else if (dimuon->Charge() > 0) {
            // fSparseLSppDimuon->Fill(fill,(double)1./fUtils->getDS());
            fTreeMixLSppDimuon->Fill();
          } else {
            // fSparseLSmmDimuon->Fill(fill,(double)1./fUtils->getDS());
            fTreeMixLSmmDimuon->Fill();
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

bool AliAnalysisTaskAODTrackPairMC::FwdMuonPairAnalysis() {

  if (!fIsMC && !(fInputHandler->IsEventSelected() & fTriggerMaskForSame))
    return false;

  if (fIsCMUL7) {
    fEventCounter->Fill(0., fUtils->getCentClass());
    fEventCounter->Fill(4., fUtils->getCentClass(),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMLL7) {
    fEventCounter->Fill(1., fUtils->getCentClass());
    fEventCounter->Fill(5., fUtils->getCentClass(),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMUL7 || fIsCMLL7) {
    fEventCounter->Fill(2., fUtils->getCentClass());
    fEventCounter->Fill(6., fUtils->getCentClass(),
                        (double)1. / fUtils->getDS());
  }
  if (fIsCMUL7 && fIsCMLL7) {
    fEventCounter->Fill(3., fUtils->getCentClass());
    fEventCounter->Fill(7., fUtils->getCentClass(),
                        (double)1. / fUtils->getDS());
  }

  Int_t nTrack = fEvent->GetNumberOfTracks();

  AliAODTrack *track1;
  AliAODTrack *track2;

  AliAODDimuon *dimuon;
  AliAODMCParticle *particle1;
  AliAODMCParticle *particle2;
  AliAODMCParticle *particle12;

  TLorentzVector muon1, muon2, muon12;

  for (Int_t iTrack1 = 0; iTrack1 < nTrack; ++iTrack1) {

    track1 = (AliAODTrack *)fEvent->GetTrack(iTrack1);

    FillingRecMuonTree(track1);

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

      if (track1->GetLabel() < 0 || track2->GetLabel() < 0)
        continue;

      RecDimuonPt = dimuon->Pt();
      RecDimuonRap = fabs(dimuon->Y());
      RecDimuonMass = dimuon->M();
      RecDimuonCent = fUtils->getCentClass();
      RecDimuonDS = fUtils->getDS();

      if (fUtils->isSameMotherPair(track1, track2)) {

        int mom_pdg = fUtils->getMotherPdgCode(track1);
        int mom_label = fUtils->getMotherLabel(track1);

        if (mom_label < 0) {
          continue;
        }

        particle12 = (AliAODMCParticle *)fMCTrackArray->At(mom_label);

        if (!particle12) {
          continue;
        }

        RecMCMomDalitz = 0;
        RecMCMom2Body = 0;
        RecMCMomPdgCode = mom_pdg;

        if (mom_pdg == fUtils->fPdgCodeEta) {
          if (fUtils->isDalitzProd()) {
            RecMCMomDalitz = 1;
          } else {
            RecMCMom2Body = 1;
          }
        } else if (mom_pdg == fUtils->fPdgCodeRho) {
          if (fUtils->is2BodyProd()) {
            RecMCMom2Body = 1;
          }
        } else if (mom_pdg == fUtils->fPdgCodeOmega) {
          if (fUtils->isDalitzProd()) {
            RecMCMomDalitz = 1;
          } else {
            RecMCMom2Body = 1;
          }
        } else if (mom_pdg == fUtils->fPdgCodePhi) {
          if (fUtils->is2BodyProd()) {
            RecMCMom2Body = 1;
          }
        } else if (mom_pdg == fUtils->fPdgCodeEtaPrime) {
          if (fUtils->isDalitzProd()) {
            RecMCMomDalitz = 1;
          }
        }

        particle1 = (AliAODMCParticle *)fMCTrackArray->At(track1->GetLabel());
        particle2 = (AliAODMCParticle *)fMCTrackArray->At(track2->GetLabel());

        muon1.SetPtEtaPhiM(particle1->Pt(), particle1->Eta(), particle1->Phi(),
                           particle1->M());
        muon2.SetPtEtaPhiM(particle2->Pt(), particle2->Eta(), particle2->Phi(),
                           particle2->M());
        muon12 = muon1 + muon2;

        RecMCDimuonPt = muon12.Pt();
        RecMCDimuonRap = fabs(muon12.Rapidity());
        RecMCDimuonMass = muon12.M();

        RecMCMomPt = particle12->Pt();
        RecMCMomEta = particle12->Eta();

        if (mom_pdg == fUtils->fPdgCodeRho) {
          RecMCMomMass = RecMCDimuonMass;
        } else {
          RecMCMomMass = particle12->M();
        }

        if (dimuon->Charge() == 0) {
          fTreeULSDimuon->Fill();
        } else if (dimuon->Charge() > 0) {
          fTreeLSppDimuon->Fill();
        } else {
          fTreeLSmmDimuon->Fill();
        }

      } else {

        particle1 = (AliAODMCParticle *)fMCTrackArray->At(track1->GetLabel());
        particle2 = (AliAODMCParticle *)fMCTrackArray->At(track2->GetLabel());

        if (!particle1 || !particle2) {
          continue;
        }
        if (!isPrimaryMuonTrack(particle1) || !isPrimaryMuonTrack(particle2)) {
          continue;
        }

        if (!fUtils->isHeavyFlavorOrigin(particle1) ||
            !fUtils->isHeavyFlavorOrigin(particle2)) {
          // continue;
        }

        TLorentzVector muon1, muon2, muon12;
        muon1.SetPtEtaPhiM(particle1->Pt(), particle1->Eta(), particle1->Phi(),
                           particle1->M());
        muon2.SetPtEtaPhiM(particle2->Pt(), particle2->Eta(), particle2->Phi(),
                           particle2->M());
        muon12 = muon1 + muon2;

        RecMCDimuonPt = muon12.Pt();
        RecMCDimuonRap = fabs(muon12.Rapidity());
        RecMCDimuonMass = muon12.M();

        if (dimuon->Charge() == 0) {
          fTreeULSDimuon->Fill();
        } else if (dimuon->Charge() > 0) {
          fTreeLSppDimuon->Fill();
        } else {
          fTreeLSmmDimuon->Fill();
        }
      }

      delete dimuon;

    } // end of loop track2
  }   // end of loop track1
  return true;
}

bool AliAnalysisTaskAODTrackPairMC::isPrimaryMuonTrack(
    AliAODMCParticle *particle1) {

  if (!particle1)
    return false;

  if (!TDatabasePDG::Instance()->GetParticle(particle1->GetPdgCode()))
    return false;
  if (!fUtils->isPrimary(particle1))
    return false;
  if (fabs(particle1->GetPdgCode()) != 13)
    return false;

  return true;
}

bool AliAnalysisTaskAODTrackPairMC::processMC() {

  AliAODMCParticle *particle1;
  AliAODMCParticle *particle2;
  AliAODMCParticle *particle12;

  TLorentzVector muon1, muon2, muon12;

  bool detect[2] = {};
  bool accept[2] = {};
  int detect_track_label[2] = {};

  for (Int_t iTrack1 = 0; iTrack1 < fMCTrackArray->GetEntries(); ++iTrack1) {

    detect[0] = false;
    accept[0] = false;

    particle1 = (AliAODMCParticle *)fMCTrackArray->At(iTrack1);

    if (!isPrimaryMuonTrack(particle1)) {
      continue;
    }

    for (Int_t iMCHTrack = 0; iMCHTrack < fEvent->GetNumberOfTracks();
         ++iMCHTrack) {
      AliAODTrack *track = (AliAODTrack *)fEvent->GetTrack(iMCHTrack);
      int label = track->GetLabel();
      if (iTrack1 == label) {
        detect[0] = true;
        detect_track_label[0] = iMCHTrack;
      }
    }

    if (171.0 * TMath::Pi() / 180. < particle1->Theta() &&
        particle1->Theta() < 178.0 * TMath::Pi() / 180.) {
      accept[0] = true;
    }

    MCMuonPt = particle1->Pt();
    MCMuonEta = particle1->Eta();
    MCMuonRap = particle1->Y();
    MCMuonPhi = particle1->Phi();
    MCMuonDetect = detect[0];

    if (particle1->Charge() > 0) {
      fTreeMCMuonP->Fill();
    } else {
      fTreeMCMuonN->Fill();
    }

    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < fMCTrackArray->GetEntries();
         ++iTrack2) {

      detect[1] = false;
      accept[1] = false;

      particle2 = (AliAODMCParticle *)fMCTrackArray->At(iTrack2);

      if (!isPrimaryMuonTrack(particle2)) {
        continue;
      }

      for (Int_t iMCHTrack = 0; iMCHTrack < fEvent->GetNumberOfTracks();
           ++iMCHTrack) {
        AliAODTrack *track = (AliAODTrack *)fEvent->GetTrack(iMCHTrack);
        int label = track->GetLabel();
        if (iTrack2 == label) {
          detect[1] = true;
          detect_track_label[1] = iMCHTrack;
        }
      }
      if (171.0 * TMath::Pi() / 180. < particle2->Theta() &&
          particle2->Theta() < 178.0 * TMath::Pi() / 180.) {
        accept[1] = true;
      }

      int mom_pdg1 = fUtils->getMotherPdgCode(particle1);
      int mom_pdg2 = fUtils->getMotherPdgCode(particle2);
      int mom_label1 = fUtils->getMotherLabel(particle1);
      int mom_label2 = fUtils->getMotherLabel(particle2);

      if (fUtils->isSameMotherPair(particle1, particle2)) {

        if (mom_label1 < 0) {
          continue;
        }

        particle12 = (AliAODMCParticle *)fMCTrackArray->At(mom_label1);

        MCMomPt = particle12->Pt();
        MCMomEta = fabs(particle12->Eta());

        if (!particle12) {
          continue;
        }

        MCMomDalitz = 0;
        MCMom2Body = 0;
        MCMomPdgCode = mom_pdg1;

        if (accept[0] == true && accept[1] == true && detect[0] == true &&
            detect[1] == true) {
          MCDimuonDetected = 2;
        } else if (accept[0] == true && accept[1] == true) {
          MCDimuonDetected = 1;
        } else {
          MCDimuonDetected = 0;
        }

        if (mom_pdg1 == fUtils->fPdgCodeEta) {
          if (fUtils->isDalitzProd()) {
            MCMomDalitz = 1;
          } else {
            MCMom2Body = 1;
          }
        } else if (mom_pdg1 == fUtils->fPdgCodeRho) {
          if (fUtils->is2BodyProd()) {
            MCMom2Body = 1;
          }
        } else if (mom_pdg1 == fUtils->fPdgCodeOmega) {
          if (fUtils->isDalitzProd()) {
            MCMomDalitz = 1;
          } else {
            MCMom2Body = 1;
          }
        } else if (mom_pdg1 == fUtils->fPdgCodePhi) {
          if (fUtils->is2BodyProd()) {
            MCMom2Body = 1;
          }
        } else if (mom_pdg1 == fUtils->fPdgCodeEtaPrime) {
          if (fUtils->isDalitzProd()) {
            MCMomDalitz = 1;
          }
        }

        muon1.SetPtEtaPhiM(particle1->Pt(), particle1->Eta(), particle1->Phi(),
                           particle1->M());
        muon2.SetPtEtaPhiM(particle2->Pt(), particle2->Eta(), particle2->Phi(),
                           particle2->M());
        muon12 = muon1 + muon2;

        MCDimuonPt = muon12.Pt();
        MCDimuonRap = fabs(muon12.Rapidity());
        MCDimuonMass = muon12.M();

        RecDimuonPt = 0;
        RecDimuonRap = 0;
        RecDimuonMass = 0;

        if (MCDimuonDetected &&
            detect_track_label[0] != detect_track_label[1]) {
          AliAODTrack *track1 =
              (AliAODTrack *)fEvent->GetTrack(detect_track_label[0]);
          AliAODTrack *track2 =
              (AliAODTrack *)fEvent->GetTrack(detect_track_label[1]);
          if (track1 && track2) {
            AliAODDimuon *dimuon = new AliAODDimuon();
            dimuon->SetMuons(track1, track2);
            if (dimuon->Charge() == 0) {
              RecDimuonPt = dimuon->Pt();
              RecDimuonRap = fabs(dimuon->Y());
              RecDimuonMass = dimuon->M();
            }
            delete dimuon;
          }
        }

        if (particle1->Charge() + particle2->Charge() == 0) {
          fTreeMCULSDimuon->Fill();
        } else if (particle1->Charge() + particle2->Charge() > 0) {
          fTreeMCLSppDimuon->Fill();
        } else {
          fTreeMCLSmmDimuon->Fill();
        }

      } else {

        if (!particle1 || !particle2) {
          continue;
        }

        if (!fUtils->isHeavyFlavorOrigin(particle1) ||
            !fUtils->isHeavyFlavorOrigin(particle2)) {
          continue;
        }

        muon1.SetPtEtaPhiM(particle1->Pt(), particle1->Eta(), particle1->Phi(),
                           particle1->M());
        muon2.SetPtEtaPhiM(particle2->Pt(), particle2->Eta(), particle2->Phi(),
                           particle2->M());
        muon12 = muon1 + muon2;

        MCDimuonPt = muon12.Pt();
        MCDimuonRap = fabs(muon12.Rapidity());
        MCDimuonMass = muon12.M();

        if (particle1->Charge() + particle2->Charge() == 0) {
          fTreeMCULSDimuon->Fill();
        } else if (particle1->Charge() + particle2->Charge() > 0) {
          fTreeMCLSppDimuon->Fill();
        } else {
          fTreeMCLSmmDimuon->Fill();
        }
      }
    }
  }

  return true;
}
