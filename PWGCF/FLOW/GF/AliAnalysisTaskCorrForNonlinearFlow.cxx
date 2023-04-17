/* -------------------------------------------
 * Maintainer: Mingrui Zhao
 */
#include "AliAnalysisTaskCorrForNonlinearFlow.h"
#include "AliGFWCuts.h"
#include "AliGFWNFCuts.h"
#include "AliGFWWeights.h"
#include "CorrelationCalculator.h"

#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TComplex.h>
#include <TBits.h>
// AliRoot includes
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TComplex.h>

// AliRoot includes
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODITSsaTrackCuts.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"

// STL includes
#include <iostream>
using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskCorrForNonlinearFlow)
// ---------------------------------------------------------------------------------
AliAnalysisTaskCorrForNonlinearFlow::AliAnalysisTaskCorrForNonlinearFlow():
AliAnalysisTaskSE(),
  fEventCuts(),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fAOD(0),
  fitssatrackcuts(0),
  fEtaCut(0.8),
  fVtxCut(10.0),
  fVtxCutDefault(10.0),
  fMinPt(0.2),
  fMaxPt(3.0),
  fSample(1),
  fTrigger(0),
  fAliTrigger(0),
  fNUE(0),
  fNUA(0),
  fIsMC(0),
  fUseTPCTruth(0),
  fUseFMDTruth(0),
  fNtrksName("Mult"),
//....
  fPeriod("LHC15o"),
  fCurrSystFlag(0),
  fSpringMode(false),
  fLowMultiplicityMode(false),
  fAddTPCPileupCuts(false),
  fESDvsTPConlyLinearCut(15000.),
  fUseCorrectedNTracks(true),
  fUseFlippedEta(false),
  fUseNarrowBin(false),
  fExtremeEfficiency(0),
  fTPCchi2perCluster(4.0),
  fUseAdditionalDCACut(false),
  fUseDefaultWeight(false),
  fEtaGap3Sub(0.4),

  fListOfObjects(0),
  fListOfProfile(0),

  fMultTOFLowCut(0),
  fMultTOFHighCut(0),
  fMultCentLowCut(0),

  fTrackEfficiency(0),
  hTrackEfficiency(0),
  hTrackEfficiencyRun(0),

  fFlowRunByRunWeights(false),
  fFlowPeriodWeights(false),
  fFlowUse3Dweights(false),
  fFlowWeightsList(nullptr),
  fFlowPtWeightsList(nullptr),
  fFlowFeeddownList(nullptr),
  fFlowPtWeightsFile(nullptr),

  fPhiWeight(0),
  fPhiWeightFile(0),
  fPhiWeightPlus(0),
  fPhiWeightMinus(0),
  fWeightsSystematics(0),
  fPtWeightsSystematics(0),
  hPhiWeight(0),
  hPhiWeightRun(0),
  hPhiWeight1D(0),

  hEventCount(0),
  hMult(0),
  fVtxAfterCuts(0),
  fCentralityDis(0),
  fV0CentralityDis(0),
  hMultV0vsNtrksAfterCuts(0),
  hMultSPDvsNtrksAfterCuts(0),
  hNtrksVSmultPercentile(0),
  fCentralityV0MCL1(0),
  fCentralityV0MCL0(0),
  fCentralityCL0CL1(0),
  fMultvsCentr(0),
  fMult128vsCentr(0),
  fMultTPCvsTOF(0),
  fMultTPCvsESD(0),

  hSPDClsVsTrk(0),
  hV0C012vsTkl(0),
  hV0C012vsV0C3(0),
  hV0MOnVsOf(0),
  hSPDOnVsOf(0),

  fPhiDis1D(0),
  fPhiDis(0),
  fEtaTriDis(0),
  fEtaTriDisBefore(0),
  fPtTriDis(0),
  fPtTriDisBefore(0),
  fEtaAssDis(0),
  fEtaAssDisBefore(0),
  fPtAssDis(0),
  fPtAssDisBefore(0),
  hDCAxyBefore(0),
  hDCAzBefore(0),
  hITSclustersBefore(0),
  hChi2Before(0),
  hDCAxy(0),
  hDCAz(0),
  hITSclusters(0),
  hChi2(0),
  fPoolMaxNEvents(2000),
  fPoolMinNTracks(50000),
  fMinEventsToMix(5),
  fBootstrapStat(true),
  fUsePhiStarCut(kTRUE),
  fUseFMDcut(kTRUE),
  fFMDcutapar0(1.64755),
  fFMDcutapar1(119.602),
  fFMDcutcpar0(2.73426),
  fFMDcutcpar1(150.31),
  fFMDAacceptanceCutLower(1.8),
  fFMDAacceptanceCutUpper(4.8),
  fFMDCacceptanceCutLower(-3.2),
  fFMDCacceptanceCutUpper(-1.8),
  nSamples(10),
  sampleLow(-100),
  sampleHigh(100),
  fFastMode(false),
  fAlternativePhiBinning(false),
  rand(2333)
{
}
//______________________________________________________________________________
AliAnalysisTaskCorrForNonlinearFlow::AliAnalysisTaskCorrForNonlinearFlow(const char *name, int _fNUA, int _fNUE, TString _fPeriod):
  AliAnalysisTaskSE(name),
  fEventCuts(),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fAOD(0),
  fitssatrackcuts(0),
  fEtaCut(0.8),
  fVtxCut(10.0),
  fVtxCutDefault(10.0),
  fMinPt(0.2),
  fMaxPt(3.0),
  fSample(1),
  fTrigger(0),
  fAliTrigger(0),
  fNUE(_fNUE),
  fNUA(_fNUA),
  fIsMC(0),
  fUseTPCTruth(0),
  fUseFMDTruth(0),
  fNtrksName("Mult"),
  //....
  fPeriod(_fPeriod),
  fCurrSystFlag(0),
  fSpringMode(false),
  fLowMultiplicityMode(false),
  fAddTPCPileupCuts(false),
  fESDvsTPConlyLinearCut(15000.),
  fUseCorrectedNTracks(true),
  fUseFlippedEta(false),
  fUseNarrowBin(false),
  fExtremeEfficiency(0),
  fTPCchi2perCluster(4.0),
  fUseAdditionalDCACut(false),
  fUseDefaultWeight(false),
  fEtaGap3Sub(0.4),

  fListOfObjects(0),
  fListOfProfile(0),

  fMultTOFLowCut(0),
  fMultTOFHighCut(0),
  fMultCentLowCut(0),

  fTrackEfficiency(0),
  hTrackEfficiency(0),
  hTrackEfficiencyRun(0),

  fFlowRunByRunWeights(false),
  fFlowPeriodWeights(false),
  fFlowUse3Dweights(false),
  fFlowWeightsList(nullptr),
  fFlowPtWeightsList(nullptr),
  fFlowFeeddownList(nullptr),
  fFlowPtWeightsFile(nullptr),

  fPhiWeight(0),
  fPhiWeightFile(0),
  fPhiWeightPlus(0),
  fPhiWeightMinus(0),
  fWeightsSystematics(0),
  fPtWeightsSystematics(0),
  hPhiWeight(0),
  hPhiWeightRun(0),
  hPhiWeight1D(0),
  hEventCount(0),
  hMult(0),
  fVtxAfterCuts(0),
  fCentralityDis(0),
  fV0CentralityDis(0),
  hMultV0vsNtrksAfterCuts(0),
  hMultSPDvsNtrksAfterCuts(0),
  hNtrksVSmultPercentile(0),
  fCentralityV0MCL1(0),
  fCentralityV0MCL0(0),
  fCentralityCL0CL1(0),
  fMultvsCentr(0),
  fMult128vsCentr(0),
  fMultTPCvsTOF(0),
  fMultTPCvsESD(0),

  hSPDClsVsTrk(0),
  hV0C012vsTkl(0),
  hV0C012vsV0C3(0),
  hV0MOnVsOf(0),
  hSPDOnVsOf(0),

  fPhiDis1D(0),
  fPhiDis(0),
  fEtaTriDis(0),
  fEtaTriDisBefore(0),
  fPtTriDis(0),
  fPtTriDisBefore(0),
  fEtaAssDis(0),
  fEtaAssDisBefore(0),
  fPtAssDis(0),
  fPtAssDisBefore(0),
  hDCAxyBefore(0),
  hDCAzBefore(0),
  hITSclustersBefore(0),
  hChi2Before(0),
  hDCAxy(0),
  hDCAz(0),
  hITSclusters(0),
  hChi2(0),
  fPoolMaxNEvents(2000),
fPoolMinNTracks(50000),
fMinEventsToMix(5),
fBootstrapStat(true),
  fUsePhiStarCut(kTRUE),
  fUseFMDcut(kTRUE),
  fFMDcutapar0(1.64755),
  fFMDcutapar1(119.602),
  fFMDcutcpar0(2.73426),
  fFMDcutcpar1(150.31),
  fFMDAacceptanceCutLower(1.8),
  fFMDAacceptanceCutUpper(4.8),
  fFMDCacceptanceCutLower(-3.2),
  fFMDCacceptanceCutUpper(-1.8),
  nSamples(10),
  sampleLow(-100),
  sampleHigh(100),
  fFastMode(false),
  fAlternativePhiBinning(false),
  rand(2333)
{

  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  int inputslot = 1;
  if (fNUA) {
    DefineInput(inputslot, TList::Class());
    inputslot++;
  }
  if (fNUE) {
    DefineInput(inputslot, TList::Class());
    inputslot++;
    if (fPeriod.EqualTo("LHC16qt")) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    }
    if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") ||
        fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")
        ) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    }
  }
}

//______________________________________________________________________________
AliAnalysisTaskCorrForNonlinearFlow::AliAnalysisTaskCorrForNonlinearFlow(const char *name):
  AliAnalysisTaskSE(name),
  fEventCuts(),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fAOD(0),
  fitssatrackcuts(0),
  fEtaCut(0.8),
  fVtxCut(10.0),
  fVtxCutDefault(10.0),
  fMinPt(0.2),
  fMaxPt(3.0),
  fSample(1),
  fTrigger(0),
  fAliTrigger(0),
  fNUE(0),
  fNUA(0),
  fNtrksName("Mult"),
  //....
  fPeriod("LHC15o"),
  fCurrSystFlag(0),
  fSpringMode(false),
  fLowMultiplicityMode(false),
  fAddTPCPileupCuts(false),
  fESDvsTPConlyLinearCut(15000.),
  fUseCorrectedNTracks(true),
  fUseFlippedEta(false),
  fUseNarrowBin(false),
  fExtremeEfficiency(0),
  fTPCchi2perCluster(4.0),
  fUseAdditionalDCACut(false),
  fUseDefaultWeight(false),
  fEtaGap3Sub(0.4),

  fListOfObjects(0),
  fListOfProfile(0),

  fMultTOFLowCut(0),
  fMultTOFHighCut(0),
  fMultCentLowCut(0),

  fTrackEfficiency(0),
  hTrackEfficiency(0),
  hTrackEfficiencyRun(0),

  fFlowRunByRunWeights(false),
  fFlowPeriodWeights(false),
  fFlowUse3Dweights(false),
  fFlowWeightsList(nullptr),
  fFlowPtWeightsList(nullptr),
  fFlowFeeddownList(nullptr),
  fFlowPtWeightsFile(nullptr),

  fPhiWeight(0),
  fPhiWeightFile(0),
  fPhiWeightPlus(0),
  fPhiWeightMinus(0),

  fWeightsSystematics(0),
  fPtWeightsSystematics(0),

  hPhiWeight(0),
  hPhiWeightRun(0),
  hPhiWeight1D(0),

  hEventCount(0),
  hMult(0),
  fVtxAfterCuts(0),
  fCentralityDis(0),
  fV0CentralityDis(0),
  hMultV0vsNtrksAfterCuts(0),
  hMultSPDvsNtrksAfterCuts(0),
  hNtrksVSmultPercentile(0),
  fCentralityV0MCL1(0),
  fCentralityV0MCL0(0),
  fCentralityCL0CL1(0),
  fMultvsCentr(0),
  fMult128vsCentr(0),
  fMultTPCvsTOF(0),
  fMultTPCvsESD(0),

  hSPDClsVsTrk(0),
  hV0C012vsTkl(0),
  hV0C012vsV0C3(0),
  hV0MOnVsOf(0),
  hSPDOnVsOf(0),

  fPhiDis1D(0),
  fPhiDis(0),
  fEtaTriDis(0),
  fEtaTriDisBefore(0),
  fPtTriDis(0),
  fPtTriDisBefore(0),
  fEtaAssDis(0),
  fEtaAssDisBefore(0),
  fPtAssDis(0),
  fPtAssDisBefore(0),
  hDCAxyBefore(0),
  hDCAzBefore(0),
  hITSclustersBefore(0),
  hChi2Before(0),
  hDCAxy(0),
  hDCAz(0),
  hITSclusters(0),
  hChi2(0),
  fPoolMaxNEvents(2000),
fPoolMinNTracks(50000),
  fMinEventsToMix(5),
fBootstrapStat(true),
  fUsePhiStarCut(kTRUE),
  fUseFMDcut(kTRUE),
  fFMDcutapar0(1.64755),
  fFMDcutapar1(119.602),
  fFMDcutcpar0(2.73426),
  fFMDcutcpar1(150.31),
  fFMDAacceptanceCutLower(1.8),
  fFMDAacceptanceCutUpper(4.8),
  fFMDCacceptanceCutLower(-3.2),
  fFMDCacceptanceCutUpper(-1.8),
  nSamples(10),
  sampleLow(-100),
  sampleHigh(100),
  fFastMode(false),
  fAlternativePhiBinning(false),
  rand(2333)
{

  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  int outputslot = 2;

  // int inputslot = 1;
  DefineInput(1, TList::Class());
  DefineInput(2, TList::Class());
}


// ---------------------------------------------------------------------------------
AliAnalysisTaskCorrForNonlinearFlow::~AliAnalysisTaskCorrForNonlinearFlow() {
  // Destructor
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskCorrForNonlinearFlow::UserCreateOutputObjects() {
  // Create output objects
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  // Setting for AliEventCuts:
  fEventCuts.AddQAplotsToList(fListOfObjects);

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n") ) {
    // Only for LHC15o pass1
    fGFWSelection15o = new AliGFWNFCuts();
    fGFWSelection15o->PrintSetup();
  } else {
    fGFWSelection = new AliGFWCuts();
    fGFWSelection->PrintSetup();
  }

  hEventCount = new TH1D("hEventCount", "; centrality;;", 1, 0, 1);
  fListOfObjects->Add(hEventCount);

  // BinMethod 1:Bootstrap, 2:Pt, 4:Nch, 8:Vtx
  std::vector<Double_t>   fCentBins; 
  std::vector<Double_t>   fCentBinsForMixing;
  fCentBinsForMixing.assign({0,5});
  if (fBinMethod & 4) fCentBins.assign({0,5,10,15,20,25,30,40,50,60,70,80,90,100,110,120,130,140,150}); 
  else fCentBins.assign({0, 150});
  for (int i = 0; i < fCentBins.size(); i++) fCentBins[i] += 0.5;
  for (int i = 0; i < fCentBinsForMixing.size(); i++) fCentBinsForMixing[i] += 0.5;

  hMult = new TH1F("hMult", ";number of tracks; entries", fCentBins.size()-1, fCentBins.data());
  hMult->Sumw2();
  fListOfObjects->Add(hMult);

  // fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
  // fVtxAfterCuts->Sumw2();
  // fListOfObjects->Add(fVtxAfterCuts);

  // fCentralityDis = new TH1F("fCentralityDis", "centrality distribution; centrality; Counts", 100, 0, 100);
  // fListOfObjects->Add(fCentralityDis);

    
  Int_t nSteps = 1;
  Double_t binning_deta_tpctpc[37] = {-1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,  0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6, 1.7, 1.8};
  Double_t binning_deta_tpcfmd[43]={-6.,-5.8, -5.6, -5.4, -5.2, -5.0, -4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4, -3.2, -3., -2.8, -2.6, -2.4, -2.2, -2., -1.8, -1.6, -1.4, -1.2, -1., -0.8, 1., 1.2, 1.4, 1.6, 1.8, 2. , 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4.};
  std::vector<Double_t>   fPtBinsTrigCharged;
  if (fBinMethod & 2) fPtBinsTrigCharged.assign({0.2, 1.0, 3.0, 4.0});
  else fPtBinsTrigCharged.assign({0.2, 3.0});
  std::vector<Double_t>   fzVtxBins;
  if (fBinMethod & 8) fzVtxBins.assign({-10.0, -8.0, -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0, 8.0, 10.0});
  else fzVtxBins.assign({-10.0, 10.0});
  std::vector<Double_t>   fzVtxBinsForMixing;
  fzVtxBinsForMixing.assign({-10.0, -8.0, -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0, 8.0, 10.0});

  Int_t sizeEta = 0;
  if (anaType.EqualTo("TPCTPC")) sizeEta = 36;
  else if (anaType.EqualTo("TPCFMD")) sizeEta = 42;
  else if (anaType.EqualTo("FMDFMD")) sizeEta = 24;

  const Int_t sizeCent = fCentBins.size() - 1;
  Int_t sizeOfSamples = (fBinMethod & 1) ? nSamples : 1;
  Int_t sizeOfVtxZbins = fzVtxBins.size() - 1; // (Int_t) fNOfSamples; 

  Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
  int nPhiBins = 80;
  if (fAlternativePhiBinning) nPhiBins = 72;
  const Int_t iBinning[] = {sizeEta,nPhiBins,sizeOfVtxZbins,sizeOfSamples,sizePtTrig,sizeCent};

  fEtaAssDisBefore = new TH1D("hEtaAssDisBefore", "eta distribution", 1000, -10, 10);
  fListOfObjects->Add(fEtaAssDisBefore);
  fPtAssDisBefore = new TH1D("hPtAssDisBefore", "pt distribution", 100, 0, 5);
  fListOfObjects->Add(fPtAssDisBefore);
  fPhiAssDisBefore = new TH1D("hPhiAssDisBefore", "phi distribution", nPhiBins, 0, TMath::Pi() / 2 * 4);
  fListOfObjects->Add(fPhiAssDisBefore);
  fEtaTriDisBefore = new TH1D("hEtaTriDisBefore", "eta distribution", 1000, -10, 10);
  fListOfObjects->Add(fEtaTriDisBefore);
  fPtTriDisBefore = new TH1D("hPtTriDisBefore", "pt distribution", 100, 0, 5);
  fListOfObjects->Add(fPtTriDisBefore);
  fPhiTriDisBefore = new TH1D("hPhiTriDisBefore", "phi distribution", nPhiBins, 0, TMath::Pi() / 2 * 4);
  fListOfObjects->Add(fPhiTriDisBefore);


  fEtaAssDis = new TH1D("hEtaAssDis", "eta distribution", 1000, -10, 10);
  fListOfObjects->Add(fEtaAssDis);
  fPtAssDis = new TH1D("hPtAssDis", "pt distribution", 100, 0, 5);
  fListOfObjects->Add(fPtAssDis);
  fPhiAssDis = new TH1D("hPhiAssDis", "phi distribution", nPhiBins, 0, TMath::Pi() / 2 * 4);
  fListOfObjects->Add(fPhiAssDis);
  fEtaTriDis = new TH1D("hEtaTriDis", "eta distribution", 1000, -10, 10);
  fListOfObjects->Add(fEtaTriDis);
  fPtTriDis = new TH1D("hPtTriDis", "pt distribution", 100, 0, 5);
  fListOfObjects->Add(fPtTriDis);
  fPhiTriDis = new TH1D("hPhiTriDis", "phi distribution", nPhiBins, 0, TMath::Pi() / 2 * 4);
  fListOfObjects->Add(fPhiTriDis);

  hFMDAvsV0 = new TH2D("hFMDAvsV0", "FMDA V0A correlation", 100, 0, 100, 100, 0, 100);
  hFMDCvsV0 = new TH2D("hFMDCvsV0", "FMDC V0C correlation", 100, 0, 100, 100, 0, 100);

  Int_t inSlotCounter=1;
  if(fNUA) {
    if (fPeriod.EqualTo("LHC15oKatarina") ) {
      fFlowWeightsList = (TList*) GetInputData(inSlotCounter);
    } else {
      fFlowWeightsList = (TList*) GetInputData(inSlotCounter);
    }
    inSlotCounter++;
  };
  if(fNUE) {
    if (fPeriod.EqualTo("LHC15oKatarina") ) {
      fFlowPtWeightsList = (TList*) GetInputData(inSlotCounter);
    } else if (fPeriod.EqualTo("LHC16qt") ||
               fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") ||
               fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")
               ) {
      fFlowPtWeightsList = (TList*) GetInputData(inSlotCounter);
      inSlotCounter++;
      fFlowFeeddownList = (TList*) GetInputData(inSlotCounter);
    } else {
      fFlowPtWeightsList = (TList*) GetInputData(inSlotCounter);
    }
    inSlotCounter++;
  };

  fListOfProfile = new TList();
  fListOfProfile->SetOwner();

  fhTracksTrigCent = new TH2D("fhTrigTracks", "fhTrigTracks; cent; PVz", sizeCent, fCentBins.data(), fzVtxBins.size()-1, fzVtxBins.data());
  fListOfProfile->Add(fhTracksTrigCent);

  // mixing
  int fNCentBinsForMixing = fCentBinsForMixing.size() - 1;
  int fNzVtxBinsForMixing = fzVtxBinsForMixing.size() - 1;
  fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBinsForMixing, fCentBinsForMixing.data(), fNzVtxBinsForMixing, fzVtxBinsForMixing.data());

  if (!fPoolMgr) {
    AliError("Event Pool manager not created!");
    return;
  }
  fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);

  fhChargedSE = new AliTHn("fhChargedSE", "fhChargedSE", nSteps, 6, iBinning);
  if (anaType.EqualTo("TPCTPC")) {
    fhChargedSE->SetBinLimits(0, binning_deta_tpctpc);
  } else if (anaType.EqualTo("TPCFMD")) {
    fhChargedSE->SetBinLimits(0, binning_deta_tpcfmd);
  } else if (anaType.EqualTo("FMDFMD")) {
    fhChargedSE->SetBinLimits(0, 3.4, 8.2);
  }
  fhChargedSE->SetBinLimits(1, -TMath::Pi() / 2, TMath::Pi() / 2 * 3);
  fhChargedSE->SetBinLimits(2, -10,10);
  fhChargedSE->SetBinLimits(3,  -0.5,sizeOfSamples-0.5);
  fhChargedSE->SetBinLimits(4, fPtBinsTrigCharged.data());
  fhChargedSE->SetBinLimits(5, fCentBins.data());
  fhChargedSE->SetVarTitle(0, "#Delta#eta");
  fhChargedSE->SetVarTitle(1, "#Delta#phi");
  fhChargedSE->SetVarTitle(2, "PVz");
  fhChargedSE->SetVarTitle(3, "sample");
  fhChargedSE->SetVarTitle(4, "p_{T} [GeV/c] (trig)");
  fhChargedSE->SetVarTitle(5, "Mult/Cent");
  fListOfProfile->Add(fhChargedSE);
  //
  fhChargedME = new AliTHn("fhChargedME", "fhChargedME", nSteps, 6, iBinning);
  if (anaType.EqualTo("TPCTPC")) {
    fhChargedME->SetBinLimits(0, binning_deta_tpctpc);
  } else if (anaType.EqualTo("TPCFMD")) {
    fhChargedME->SetBinLimits(0, binning_deta_tpcfmd);
  } else if (anaType.EqualTo("FMDFMD")) {
    fhChargedME->SetBinLimits(0, 3.4, 8.2);
  }

  fhChargedME->SetBinLimits(1, -TMath::Pi() / 2, TMath::Pi() / 2 * 3);
  fhChargedME->SetBinLimits(2, -10.,10.);
  fhChargedME->SetBinLimits(3,  -0.5,sizeOfSamples-0.5);
  fhChargedME->SetBinLimits(4, fPtBinsTrigCharged.data());
  fhChargedME->SetBinLimits(5, fCentBins.data());
  fhChargedME->SetVarTitle(0, "#Delta#eta");
  fhChargedME->SetVarTitle(1, "#Delta#phi");
  fhChargedME->SetVarTitle(2, "PVz");
  fhChargedME->SetVarTitle(3, "sample");
  fhChargedME->SetVarTitle(4, "p_{T} [GeV/c] (trig)");
  fhChargedME->SetVarTitle(5, "Mult/Cent");
  fListOfProfile->Add(fhChargedME);

  PostData(1, fListOfObjects);
  PostData(2, fListOfProfile);
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskCorrForNonlinearFlow::NotifyRun() {
  if (fAddTPCPileupCuts) {
    Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    fEventCuts.fESDvsTPConlyLinearCut[0] = fESDvsTPConlyLinearCut;
  }
}

// ---------------------------------------------------------------------------------
void AliAnalysisTaskCorrForNonlinearFlow::UserExec(Option_t *) {
  // Mingrui: apply the bootstrap later
  int sizeOfSamples = 1;
  if (fBootstrapStat) sizeOfSamples = nSamples;

  // Check if it can pass the trigger
  //..apply physics selection
  UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isTrigselected = false;
  if (fTrigger == 0) {
    isTrigselected = fSelectMask&AliVEvent::kINT7;
    fAliTrigger = AliVEvent::kINT7;
  } else if (fTrigger == 1) {
    isTrigselected = fSelectMask&AliVEvent::kHighMultV0;
    fAliTrigger = AliVEvent::kHighMultV0;
  }
  if(isTrigselected == false) return;

  //..check if I have AOD
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    return;
  }

  // Check if it passed the standard AOD selection
  if (!AcceptAOD(fAOD) ) {
    PostData(1, fListOfObjects);
    PostData(2, fListOfProfile);
    return;
  }
  hEventCount->Fill("after fEventCuts", 1.);

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n") ) { // Only for LHC15o pass1
    fGFWSelection15o->ResetCuts();
  } else {
    fGFWSelection->ResetCuts();
  }
  //..filling Vz distribution
  AliVVertex *vtx = fAOD->GetPrimaryVertex();
  float fVtxZ = vtx->GetZ();
  fPVz = fVtxZ;

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n") ) { // Only for LHC15o pass1
    if (!fGFWSelection15o->AcceptVertex(fAOD)) {
      PostData(1, fListOfObjects);
      PostData(2, fListOfProfile);
      return;
    }
  } else {
    if (!fGFWSelection->AcceptVertex(fAOD)) {
      PostData(1, fListOfObjects);
      PostData(2, fListOfProfile);
      return;
    }
  }

  if (!fIsMC) {
    AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    if(!multSelection) { return; }
    fCentrality = multSelection->GetMultiplicityPercentile("V0M");
    if (fCentrality < fCentMin || fCentrality > fCentMax) {
      PostData(1, fListOfObjects);
      PostData(2, fListOfProfile);
      return;
    }
  }

  // checking the run number for aplying weights & loading TList with weights
  //
  if (lastRunNumber != fAOD->GetRunNumber()) {
    lastRunNumber = fAOD->GetRunNumber();
    if (fPeriod.EqualTo("LHC15oKatarina")) {
      if (fNUA && !LoadWeightsKatarina()) {
        AliFatal("Trying to Load Systematics but weights not loaded!");
        return;
      }
      if (fNUE && !LoadPtWeightsKatarina()) {
        AliFatal("PtWeights not loaded!");
        return;
      }

    } else {
      if (fNUA && !LoadWeightsSystematics()) {
        AliFatal("Trying to Load Systematics but weights not loaded!");
        return;
      }
      if (fNUE && !LoadPtWeights()) {
        AliFatal("PtWeights not loaded!");
        return;
      }
    }

  }


  // Here we calcuate the multiplicity distribution

  if (fNtrksName.EqualTo("Mult")) {
    NTracksCalculation(fInputEvent);

    // Put a Ntrks cut to 100 for PbPb and XeXe
    if (fPeriod.EqualTo("LHC15o") ||
        fPeriod.EqualTo("LHC15o_pass2") ||
        fPeriod.EqualTo("LHC18qr_pass3") ||
        fPeriod.EqualTo("LHC17n")) {
      if (NtrksCounter > 100) {
        PostData(1, fListOfObjects);
        PostData(2, fListOfProfile);
        return;
      }
    }
  } else {
    NtrksCounter = fCentrality;
  }

  

  fbSign = (InputEvent()->GetMagneticField() > 0) ? 1 : -1;



  PrepareTPCFMDTracks();
  if (fIsMC) {
    PrepareTPCFMDTracksMCTruth();
  }


  if (fTracksTrigCharged) {
    FillCorrelations();
    FillCorrelationsMixed();
  }

  fTracksTrigCharged->Clear();
  delete fTracksTrigCharged;

  fTracksAss->Clear();
  delete fTracksAss;

  PostData(1, fListOfObjects);
  PostData(2, fListOfProfile);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskCorrForNonlinearFlow::Terminate(Option_t *)
{
  if (fPoolMgr) {
    // delete fPoolMgr;
  }
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForNonlinearFlow::RangePhiFMD(Double_t DPhi) {
  DPhi = TMath::ATan2(TMath::Sin(DPhi), TMath::Cos(DPhi));
  if (DPhi < (-0.5*TMath::Pi()-0.0001))    DPhi += 2 * TMath::Pi();
  return DPhi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForNonlinearFlow::RangePhi(Double_t DPhi) {
  if (DPhi < -TMath::Pi() / 2)   DPhi += 2 * TMath::Pi();
  if (DPhi > 3 * TMath::Pi() / 2) DPhi -= 2*TMath::Pi();
  return DPhi;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskCorrForNonlinearFlow::GetDPhiStar(Double_t phi1, Double_t pt1, Double_t charge1, Double_t phi2, Double_t pt2, Double_t charge2, Double_t radius){
  // calculates delta phi *
  Double_t dPhiStar = phi1 - phi2 - charge1 * fbSign * TMath::ASin(0.075 * radius / pt1) + charge2 * fbSign * TMath::ASin(0.075 * radius / pt2);

  if (dPhiStar > TMath::Pi()) dPhiStar = 2.0*TMath::Pi() - dPhiStar;
  if (dPhiStar < -TMath::Pi()) dPhiStar = -2.0*TMath::Pi() - dPhiStar;

  return dPhiStar;
}

//________________________________________________________________________
void AliAnalysisTaskCorrForNonlinearFlow::FillCorrelations() {
  if (!fTracksTrigCharged || !fTracksAss) {
    AliError("Necessary inputs missing, terminating!"); return;
  }

  /* don't check this Mingrui
     if(!fhChargedSE) { 
     AliError(Form("Output AliTHn missing for ch , terminating!")); return; 
     return;
     }
  */

  if (anaType.EqualTo("TPCTPC")) {

    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = bootstrap_value;

    // Start the two loop
    for (Int_t iTrig = 0; iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++) {
      AliVParticle* trackTrig = dynamic_cast<AliVParticle*>(fTracksTrigCharged->At(iTrig));
      
      AliAODTrack* trackTrigAOD = nullptr;
      if (!fIsMC) trackTrigAOD = dynamic_cast<AliAODTrack*>(fTracksTrigCharged->At(iTrig));

      Double_t ptTrig = trackTrig->Pt();
      Double_t phiTrig = trackTrig->Phi();
      Double_t etaTrig = trackTrig->Eta();
      Double_t chargeTrig = trackTrig->Charge();
      Double_t trigEff = 1.0; // Efficiency
      if(fNUA == 1) trigEff /= GetFlowWeightSystematics(trackTrig, fPVz, kRefs);
      if(fNUE == 1) trigEff /= GetPtWeight(ptTrig, etaTrig, fPVz, fInputEvent->GetRunNumber());

      for (Int_t iAss = 0; iAss < fTracksAss->GetEntriesFast(); iAss++) {
        AliVParticle* trackAss = dynamic_cast<AliVParticle*>(fTracksAss->At(iAss));
        AliAODTrack* trackAssAOD = nullptr;
        if (!fIsMC) trackAssAOD = dynamic_cast<AliAODTrack*>(fTracksAss->At(iAss));
        Double_t ptAss = trackAss->Pt();
        Double_t phiAss = trackAss->Phi();
        Double_t etaAss = trackAss->Eta();
        Double_t chargeAss = trackAss->Charge();
        Double_t assEff = 1.0; // Efficiency
        if(fNUA == 1) assEff /= GetFlowWeightSystematics(trackAss, fPVz, kRefs);
        if(fNUE == 1) assEff /= GetPtWeight(ptAss, etaAss, fPVz, fInputEvent->GetRunNumber());

        if (!fIsMC && trackTrigAOD->GetID() == trackAssAOD->GetID()) {
          continue;
        }

        //..check if the tracks are the same
        // Mingrui: I don't see Zuzana uses this
        // if (trackTrig == trackAss) continue;
        // Here these's a complicated way to check whether the pair pass or not
        binscont[0] = etaTrig - etaAss;
        binscont[1] = RangePhi(phiTrig - phiAss);
        binscont[4] = ptTrig;
        binscont[5] = NtrksCounter;

        if (fUsePhiStarCut) {
          double fMergingCut = 0.02;

          if(TMath::Abs(binscont[0]) < fMergingCut){
            Double_t dPhiStarLow = GetDPhiStar(phiTrig, ptTrig, chargeTrig, phiAss, ptAss, chargeAss, 0.8);
            Double_t dPhiStarHigh = GetDPhiStar(phiTrig, ptTrig, chargeTrig, phiAss, ptAss, chargeAss, 2.5);

            const Double_t kLimit = 3.0*fMergingCut;

            if(TMath::Abs(dPhiStarLow) < kLimit || TMath::Abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0 ) {
              Bool_t bIsBelow = kFALSE;
              for(Double_t rad(0.8); rad < 2.51; rad+=0.01){
                Double_t dPhiStar = GetDPhiStar(phiTrig, ptTrig, chargeTrig, phiAss, ptAss, chargeAss, rad);
                if(TMath::Abs(dPhiStar) < fMergingCut) {
                  bIsBelow = kTRUE;
                  break;
                }
              } // end loop radius
              if(bIsBelow) continue;
            }
          } 
        }
        fhChargedSE->Fill(binscont, 0, 1.0/(trigEff*assEff));
      } 
    } // end loop particle pairs
  } // endif TPC-TPC
  else if (anaType.EqualTo("TPCFMD")) {

    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = bootstrap_value;

    // Start the two loop for TPC-FMD
    for (Int_t iTrig = 0; iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++) {
      AliVParticle* trackTrig = dynamic_cast<AliVParticle*>(fTracksTrigCharged->At(iTrig));

      Double_t ptTrig = trackTrig->Pt();
      Double_t phiTrig = trackTrig->Phi();
      Double_t etaTrig = trackTrig->Eta();
      Double_t trigEff = 1.0; // Efficiency
      if(fNUA == 1) trigEff /= GetFlowWeightSystematics(trackTrig, fPVz, kRefs);
      if(fNUE == 1) trigEff /= GetPtWeight(ptTrig, etaTrig, fPVz, fInputEvent->GetRunNumber());

      for (Int_t iAss = 0; iAss < fTracksAss->GetEntriesFast(); iAss++) {

        AliPartSimpleForCorr* trackAss = dynamic_cast<AliPartSimpleForCorr*>(fTracksAss->At(iAss));
        Double_t phiAss = trackAss->Phi();
        Double_t etaAss = trackAss->Eta();

        Double_t assEff = 1.0;
        Double_t assMult = trackAss->Multiplicity(); // Efficiency

        //..check if the tracks are the same
        // Mingrui: I don't see Zuzana uses this
        // if (trackTrig == trackAss) continue;
        // Here these's a complicated way to check whether the pair pass or not
        binscont[0] = etaTrig - etaAss;
        binscont[1] = RangePhi(phiTrig - phiAss);
        binscont[4] = ptTrig;
        binscont[5] = NtrksCounter;

        fhChargedSE->Fill(binscont, 0, assMult/(trigEff*assEff));
      }
    } // end two loops
  } // endif TPC-FMD
  else if (anaType.EqualTo("FMDFMD")) {

    Double_t binscont[6];
    binscont[2] = fPVz;
    binscont[3] = bootstrap_value;

    // Start the two loop for FMD-FMD
    for (Int_t iTrig = 0; iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++) {

      AliPartSimpleForCorr* trackTrig = (AliPartSimpleForCorr*)fTracksTrigCharged->At(iTrig);
      if(!trackTrig) continue;

      Double_t phiTrig = trackTrig->Phi();
      Double_t etaTrig = trackTrig->Eta();
      Double_t trigMult = trackTrig->Multiplicity();

      for (Int_t iAss = 0; iAss < fTracksAss->GetEntriesFast(); iAss++) {

        AliPartSimpleForCorr* trackAss = dynamic_cast<AliPartSimpleForCorr*>(fTracksAss->At(iAss));
        Double_t phiAss = trackAss->Phi();
        Double_t etaAss = trackAss->Eta();

        Double_t assEff = 1.0;
        Double_t assMult = trackAss->Multiplicity(); // Efficiency

        //..check if the tracks are the same
        // Mingrui: I don't see Zuzana uses this
        // if (trackTrig == trackAss) continue;
        // Here these's a complicated way to check whether the pair pass or not
        binscont[0] = etaTrig - etaAss;
        binscont[1] = RangePhiFMD(phiTrig - phiAss);
        binscont[4] = 1.0;
        binscont[5] = NtrksCounter;

        fhChargedSE->Fill(binscont, 0, assMult*trigMult);
      }
    } // end two loops
  } // end FMD-FMD
}

//________________________________________________________________________
void AliAnalysisTaskCorrForNonlinearFlow::FillCorrelationsMixed() {
  if (!fTracksTrigCharged  || !fTracksAss) {
    AliError("Necessary inputs missing, terminating!"); return;
  }

  Double_t binscont[6];
  binscont[2] = fPVz;
  binscont[3] = bootstrap_value;

  AliEventPool* pool = fPoolMgr->GetEventPool(NtrksCounter, fPVz);
  if (!pool) {
    return;
  }

  if(pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks ||  pool->GetCurrentNEvents() > fMinEventsToMix) {
    int nMix = pool->GetCurrentNEvents();


    if (anaType.EqualTo("TPCTPC")) {
      for (Int_t iTrig = 0; iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++) {
        AliVParticle* trackTrig = dynamic_cast<AliVParticle*>(fTracksTrigCharged->At(iTrig));

        Double_t ptTrig = trackTrig->Pt();
        Double_t phiTrig = trackTrig->Phi();
        Double_t etaTrig = trackTrig->Eta();
        Double_t chargeTrig = trackTrig->Charge();

        Double_t trigEff = 1;

        if(fNUA == 1) trigEff /= GetFlowWeightSystematics(trackTrig, fPVz, kRefs);
        if(fNUE == 1) trigEff /= GetPtWeight(ptTrig, etaTrig, fPVz, fInputEvent->GetRunNumber());

        for (Int_t iMix = 0; iMix < nMix; iMix++) {
          TObjArray* mixTracks = pool->GetEvent(iMix);
          for (Int_t iAss = 0; iAss < mixTracks->GetEntriesFast(); iAss++) {

            AliVParticle* trackAss = dynamic_cast<AliVParticle*>(mixTracks->At(iAss));

            Double_t ptAss = trackAss->Pt();
            Double_t phiAss = trackAss->Phi();
            Double_t etaAss = trackAss->Eta();
            Double_t chargeAss = trackAss->Charge();

            Double_t assEff = 1;
            if(fNUA == 1) assEff /= GetFlowWeightSystematics(trackAss, fPVz, kRefs);
            if(fNUE == 1) assEff /= GetPtWeight(ptAss, etaAss, fPVz, fInputEvent->GetRunNumber());

            // We should not use this to reject self correlation in mixed event
            /*
              if (trackTrig->GetID() == trackAss->GetID()) {
              continue;
              }
            */

            //..check if the tracks are the same
            //
            // if (trackTrig == trackAss) continue;
            // Here these's a complicated way to check whether the pair pass or not
            binscont[0] = etaTrig - etaAss;
            binscont[1] = RangePhi(phiTrig - phiAss);
            binscont[4] = ptTrig;
            binscont[5] = NtrksCounter;

            if (fUsePhiStarCut) {
              double fMergingCut = 0.02;

              if(TMath::Abs(binscont[0]) < fMergingCut){
                Double_t dPhiStarLow = GetDPhiStar(phiTrig, ptTrig, chargeTrig, phiAss, ptAss, chargeAss, 0.8);
                Double_t dPhiStarHigh = GetDPhiStar(phiTrig, ptTrig, chargeTrig, phiAss, ptAss, chargeAss, 2.5);

                const Double_t kLimit = 3.0*fMergingCut;

                if(TMath::Abs(dPhiStarLow) < kLimit || TMath::Abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0 ) {
                  Bool_t bIsBelow = kFALSE;
                  for(Double_t rad(0.8); rad < 2.51; rad+=0.01){
                    Double_t dPhiStar = GetDPhiStar(phiTrig, ptTrig, chargeTrig, phiAss, ptAss, chargeAss, rad);
                    if(TMath::Abs(dPhiStar) < fMergingCut) {
                      bIsBelow = kTRUE;
                      break;
                    }
                  } // end loop radius
                  if(bIsBelow) continue;
                }
              }
            }
            fhChargedME->Fill(binscont, 0, 1.0/(trigEff*assEff)/nMix);
          } // end loop Ass
        } // end loop Mix
      } // end loop Trig
    } // end TPC-TPC
    else if (anaType.EqualTo("TPCFMD")) {

      for (Int_t iTrig = 0; iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++) {
        AliVParticle* trackTrig = dynamic_cast<AliVParticle*>(fTracksTrigCharged->At(iTrig));

        Double_t ptTrig = trackTrig->Pt();
        Double_t phiTrig = trackTrig->Phi();
        Double_t etaTrig = trackTrig->Eta();

        for (Int_t iMix = 0; iMix < nMix; iMix++) {
          TObjArray* mixTracks = pool->GetEvent(iMix);
          for (Int_t iAss = 0; iAss < mixTracks->GetEntriesFast(); iAss++) {


            AliPartSimpleForCorr* trackAss = dynamic_cast<AliPartSimpleForCorr*>(mixTracks->At(iAss));
            Double_t phiAss = trackAss->Phi();
            Double_t etaAss = trackAss->Eta();

            Double_t assEff = 1.0;
            Double_t assMult = trackAss->Multiplicity(); // Efficiency

            //..check if the tracks are the same
            //
            // if (trackTrig == trackAss) continue;
            // Here these's a complicated way to check whether the pair pass or not
            binscont[0] = etaTrig - etaAss;
            binscont[1] = RangePhi(phiTrig - phiAss);
            binscont[4] = ptTrig;
            binscont[5] = NtrksCounter;

            fhChargedME->Fill(binscont, 0, assMult/nMix);
          } // end loop Ass
        } // end loop Mix
      } // end loop Trig

    } else if (anaType.EqualTo("FMDFMD")) {
      Double_t binscont[6];
      binscont[2] = fPVz;
      binscont[3] = bootstrap_value;

      // Start the two loop for FMD-FMD
      for (Int_t iTrig = 0; iTrig < fTracksTrigCharged->GetEntriesFast(); iTrig++) {

        AliPartSimpleForCorr* trackTrig = (AliPartSimpleForCorr*)fTracksTrigCharged->At(iTrig);
        if(!trackTrig) continue;

        Double_t phiTrig = trackTrig->Phi();
        Double_t etaTrig = trackTrig->Eta();
        Double_t trigMult = trackTrig->Multiplicity();


        for (Int_t iMix = 0; iMix < nMix; iMix++) {
          TObjArray* mixTracks = pool->GetEvent(iMix);
          for (Int_t iAss = 0; iAss < mixTracks->GetEntriesFast(); iAss++) {


            AliPartSimpleForCorr* trackAss = dynamic_cast<AliPartSimpleForCorr*>(mixTracks->At(iAss));
            Double_t phiAss = trackAss->Phi();
            Double_t etaAss = trackAss->Eta();

            Double_t assEff = 1.0;
            Double_t assMult = trackAss->Multiplicity(); // Efficiency

            //..check if the tracks are the same
            //
            // if (trackTrig == trackAss) continue;
            // Here these's a complicated way to check whether the pair pass or not
            binscont[0] = etaTrig - etaAss;
            binscont[1] = RangePhiFMD(phiTrig - phiAss);
            binscont[4] = 1.0;
            binscont[5] = NtrksCounter;

            fhChargedME->Fill(binscont, 0, assMult*trigMult/nMix);
          } // end loop Ass
        } // end loop Mix
      } // end loop Trig
    }
    if (fFastMode) {
      pool->Clear();
    }
  } else {
    if (fFastMode) {
      TObjArray* cloneArray = (TObjArray*)fTracksAss->Clone();
      cloneArray->SetOwner(kTRUE);
      pool->UpdatePool(cloneArray);
    }
  }
  if (!fFastMode) {
    TObjArray* cloneArray = (TObjArray*)fTracksAss->Clone();
    cloneArray->SetOwner(kTRUE);
    pool->UpdatePool(cloneArray);
  }

  return;
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::PrepareTPCFMDTracks() {


  Int_t nTracks = fInputEvent->GetNumberOfTracks();
  fTracksTrigCharged = new TObjArray;
  fTracksAss = new TObjArray;
  //..for DCA
  double pos[3], vz, vx, vy;
  vz = fInputEvent->GetPrimaryVertex()->GetZ();
  vx = fInputEvent->GetPrimaryVertex()->GetX();
  vy = fInputEvent->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};

  // TPC-TPC
  int fNofTracksAss = 0;
  int fNofTracksTrig = 0;

  if (anaType.EqualTo("TPCTPC") && (!fIsMC || !fUseTPCTruth)) {
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      AliAODTrack* track = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(iTrack));

      // Require track to be existing and pass the track selection
      if (!track) continue;

      // Fill the QA plot before cuts
      double fb = (fCurrSystFlag == 1) ? 768 : 96;
      if (track->TestFilterBit(fb)) {
        fPtTriDisBefore->Fill(track->Pt());
        fPtAssDisBefore->Fill(track->Pt());
        fEtaTriDisBefore->Fill(track->Eta());
        fEtaAssDisBefore->Fill(track->Eta());
        fPhiTriDisBefore->Fill(track->Phi());
        fPhiAssDisBefore->Fill(track->Phi());
      }

      track->GetXYZ(pos);
      if (!AcceptAODTrack(track, pos,vtxp)) continue;

      Double_t pt = track->Pt();
      // Only if we consider the TPC-TPC correlation, we use this
      if (pt > fPtMinAss && pt < fPtMaxAss) {
        // Mingrui Polarity ??
        fTracksAss->Add(track);
        fNofTracksAss++; // number of associate tracks in the event
        // Fill the QA plot after cuts
        double weight = 1;
        if(fNUA == 1) weight = GetFlowWeightSystematics(track, fPVz, kRefs);
        fPtAssDis->Fill(track->Pt(),weight);
        fEtaAssDis->Fill(track->Eta(),weight);
        fPhiAssDis->Fill(track->Phi(),weight);
      }

      if (pt > fPtMinTrig && pt < fPtMaxTrig) {
        fTracksTrigCharged->Add(track);
        fNofTracksTrig++; // number of trigger tracks in the event
        fhTracksTrigCent->Fill(NtrksCounter, fPVz);
        double weight = 1;
        if(fNUA == 1) weight= GetFlowWeightSystematics(track, fPVz, kRefs);
        fPtTriDis->Fill(track->Pt(),weight);
        fEtaTriDis->Fill(track->Eta(),weight);
        fPhiTriDis->Fill(track->Phi(),weight);
      }
    }
  } else if (anaType.EqualTo("TPCFMD")) {
    if (!fIsMC || !fUseTPCTruth) {
      for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
        AliAODTrack* track = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(iTrack));

        // Require track to be existing and pass the track selection
        if (!track) continue;

        fPtTriDisBefore->Fill(track->Pt());
        fEtaTriDisBefore->Fill(track->Eta());
        fPhiTriDisBefore->Fill(track->Phi());

        track->GetXYZ(pos);
        if (!AcceptAODTrack(track, pos,vtxp)) continue;

        Double_t pt = track->Pt();

        if (pt > fPtMinTrig && pt < fPtMaxTrig) {
          fTracksTrigCharged->Add(track);
          fPtTriDis->Fill(track->Pt());
          fEtaTriDis->Fill(track->Eta());
          fPhiTriDis->Fill(track->Phi());
          fNofTracksTrig++; // number of trigger tracks in the event
          fhTracksTrigCent->Fill(NtrksCounter, fPVz);
        }
      }
    }

    if (!fIsMC || !fUseFMDTruth) {
      AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
      const TH2D& d2Ndetadphi = aodForward->GetHistogram();
      int nEta = d2Ndetadphi.GetXaxis()->GetNbins();
      int nPhi = d2Ndetadphi.GetYaxis()->GetNbins();


      for (int iEta = 1; iEta <= nEta; iEta++) {
        for (int iPhi = 1; iPhi <= nPhi; iPhi++) {


          double eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta); 
          double phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi); 

          double etaAccepted = false;
          if ( fFMDAacceptanceCutLower < eta && eta < fFMDAacceptanceCutUpper) etaAccepted = true;
          if ( fFMDCacceptanceCutLower < eta && eta < fFMDCacceptanceCutUpper) etaAccepted = true;
          if (!etaAccepted) continue;

          double mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
          if (mostProbableN > 0) {
            fTracksAss->Add(new AliPartSimpleForCorr(eta, phi, mostProbableN)); 
            fEtaAssDis->Fill(eta, mostProbableN);
            fPhiAssDis->Fill(phi, mostProbableN);
          }
        }
      } 
    }
  } else if (anaType.EqualTo("FMDFMD") && (!fIsMC || !fUseFMDTruth)) {
    AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
    const TH2D& d2Ndetadphi = aodForward->GetHistogram();
    int nEta = d2Ndetadphi.GetXaxis()->GetNbins();
    int nPhi = d2Ndetadphi.GetYaxis()->GetNbins();

    for (int iEta = 1; iEta <= nEta; iEta++) {
      for (int iPhi = 1; iPhi <= nPhi; iPhi++) {
        double eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta);
        double phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi);

        double etaAccepted = false;
        if ( fFMDAacceptanceCutLower < eta && eta < fFMDAacceptanceCutUpper) etaAccepted = true;
        if ( fFMDCacceptanceCutLower < eta && eta < fFMDCacceptanceCutUpper) etaAccepted = true;
        if (!etaAccepted) continue;

        double mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
        if (mostProbableN > 0) {
          if (eta > 0) {
            fTracksTrigCharged->Add(new AliPartSimpleForCorr(eta, phi, mostProbableN));
            // Set the pt to 1.0
            fhTracksTrigCent->Fill(NtrksCounter, fPVz, mostProbableN);
            fEtaTriDis->Fill(eta, mostProbableN);
            fPhiTriDis->Fill(phi, mostProbableN);
          } else {
            fTracksAss->Add(new AliPartSimpleForCorr(eta, phi, mostProbableN));
            fEtaAssDis->Fill(eta, mostProbableN);
            fPhiAssDis->Fill(phi, mostProbableN);
          }
        }
      } // end loop phi
    } // end loop eta
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::PrepareTPCFMDTracksMCTruth() {

  AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!mcEvent) return kFALSE;

  Int_t nTracks = mcEvent->GetNumberOfTracks();

  fTracksTrigCharged = new TObjArray;
  fTracksAss = new TObjArray;

  // TPC-TPC
  int fNofTracksAss = 0;
  int fNofTracksTrig = 0;

  if (anaType.EqualTo("TPCTPC") && fUseTPCTruth) {
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      AliAODMCParticle* track = (AliAODMCParticle*)mcEvent->GetTrack(iTrack);

      // Require track to be existing and pass the track selection
      if (!track) continue;

      // Fill the QA plot before cuts
      fPtTriDisBefore->Fill(track->Pt());
      fPtAssDisBefore->Fill(track->Pt());
      fEtaTriDisBefore->Fill(track->Eta());
      fEtaAssDisBefore->Fill(track->Eta());
      fPhiTriDisBefore->Fill(track->Phi());
      fPhiAssDisBefore->Fill(track->Phi());

      if (!AcceptMCTruthTrack(track)) continue;
      // Make sure it is in the TPC region
      if (track->Eta() < -0.8 || track->Eta() > 0.8) continue;

      Double_t pt = track->Pt();
      // Only if we consider the TPC-TPC correlation, we use this
      if (pt > fPtMinAss && pt < fPtMaxAss) {
        // Mingrui Polarity ??
        fTracksAss->Add(track);
        fNofTracksAss++; // number of associate tracks in the event
        // Fill the QA plot after cuts
        fPtAssDis->Fill(track->Pt());
        fEtaAssDis->Fill(track->Eta());
        fPhiAssDis->Fill(track->Phi());
      }

      if (pt > fPtMinTrig && pt < fPtMaxTrig) {
        fTracksTrigCharged->Add(track);
        fNofTracksTrig++; // number of trigger tracks in the event
        fhTracksTrigCent->Fill(NtrksCounter, fPVz);
        fPtTriDis->Fill(track->Pt());
        fEtaTriDis->Fill(track->Eta());
        fPhiTriDis->Fill(track->Phi());
      }
    }
  } else if (anaType.EqualTo("TPCFMD")) {
    if (fUseTPCTruth) {
      for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

        AliAODMCParticle* track = (AliAODMCParticle*)mcEvent->GetTrack(iTrack);
        if (!AcceptMCTruthTrack(track)) continue;
        // Make sure it is in the TPC region
        if (track->Eta() < -0.8 || track->Eta() > 0.8) continue;

        // Require track to be existing and pass the track selection
        if (!track) continue;

        fPtTriDisBefore->Fill(track->Pt());
        fEtaTriDisBefore->Fill(track->Eta());
        fPhiTriDisBefore->Fill(track->Phi());

        Double_t pt = track->Pt();

        if (pt > fPtMinTrig && pt < fPtMaxTrig) {
          fTracksTrigCharged->Add(track);
          fPtTriDis->Fill(track->Pt());
          fEtaTriDis->Fill(track->Eta());
          fPhiTriDis->Fill(track->Phi());
          fNofTracksTrig++; // number of trigger tracks in the event
          fhTracksTrigCent->Fill(NtrksCounter, fPVz);
        }
      }
    }
    if (fUseFMDTruth) {
      for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

        AliAODMCParticle* track = (AliAODMCParticle*)mcEvent->GetTrack(iTrack);
        // Require track to be existing and pass the track selection
        if (!track) continue;

        if (!AcceptMCTruthTrack(track)) continue;

        double eta = track->Eta();
        double phi = track->Phi();
        double etaAccepted = false;
        if ( fFMDAacceptanceCutLower < eta && eta < fFMDAacceptanceCutUpper) etaAccepted = true;
        if ( fFMDCacceptanceCutLower < eta && eta < fFMDCacceptanceCutUpper) etaAccepted = true;
        if (!etaAccepted) continue;

        fTracksAss->Add(new AliPartSimpleForCorr(eta, phi, 1)); 
        fEtaAssDis->Fill(eta, 1);
        fPhiAssDis->Fill(phi, 1);
      }
    }
  } else if (anaType.EqualTo("FMDFMD") && fUseFMDTruth) {
    for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

      AliAODMCParticle* track = (AliAODMCParticle*)mcEvent->GetTrack(iTrack);
      // Require track to be existing and pass the track selection
      if (!track) continue;

      if (!AcceptMCTruthTrack(track)) continue;
      // Make sure it is in the TPC region

      double eta = track->Eta();
      double phi = track->Phi();
      double etaAccepted = false;
      if ( fFMDAacceptanceCutLower < eta && eta < fFMDAacceptanceCutUpper) etaAccepted = true;
      if ( fFMDCacceptanceCutLower < eta && eta < fFMDCacceptanceCutUpper) etaAccepted = true;
      if (!etaAccepted) continue;

      if (eta > 0) {
        fTracksTrigCharged->Add(new AliPartSimpleForCorr(eta, phi, 1)); 
        // Set the pt to 1.0
        fhTracksTrigCent->Fill(NtrksCounter, fPVz, 1);
        fEtaTriDis->Fill(eta, 1);
        fPhiTriDis->Fill(phi, 1);
      } else {
        fTracksAss->Add(new AliPartSimpleForCorr(eta, phi, 1)); 
        fEtaAssDis->Fill(eta, 1);
        fPhiAssDis->Fill(phi, 1);
      }
    }
  }
  return kTRUE;
}


Bool_t AliAnalysisTaskCorrForNonlinearFlow::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp) {
  // Pt cut
  // if(mtr->Pt() < fMinPt) return kFALSE;
  // if(mtr->Pt() > fMaxPt) return kFALSE;

  // DCA cut
  if(ltrackXYZ && vtxp) {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = ltrackXYZ[2]-vtxp[2];
  } else return kFALSE; //DCA cut is a must for now

  // Additional cut for TPCchi2perCluster
  if (mtr->GetTPCchi2perCluster()>fTPCchi2perCluster) return kFALSE;

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n") ) { // Only for LHC15o pass1 and LHC17n
    return fGFWSelection15o->AcceptTrack(mtr,ltrackXYZ,0,kFALSE);
  } else {
    return fGFWSelection->AcceptTrack(mtr,ltrackXYZ,0,kFALSE);
  }
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::AcceptAOD(AliAODEvent *inEv) {
  // LHC15i, LHC15l, LHC16, LHC17, LHC18: means: pp sample
  if (fPeriod.EqualTo("LHC15i") ||
      fPeriod.EqualTo("LHC15l") ||
      fPeriod.EqualTo("LHC16Preview") ||
      fPeriod.EqualTo("LHC17Preview") ||
      fPeriod.EqualTo("LHC18Preview") || 
      fPeriod.EqualTo("LHC16") ||
      fPeriod.EqualTo("LHC17") ||
      fPeriod.EqualTo("LHC18") || 
      fPeriod.EqualTo("LHC16ZM") ||
      fPeriod.EqualTo("LHC17ZM") ||
      fPeriod.EqualTo("LHC18ZM") ) {
    if (fTrigger == 1) {
      fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
    }
  }


  if(!fEventCuts.AcceptEvent(inEv)) return false;

  // Primary vertex
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;

  // SPD Vertex
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertexSPD());
  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;

  if (fPeriod.EqualTo("LHC15o") ||
      fPeriod.EqualTo("LHC15o_pass2") ||
      fPeriod.EqualTo("LHC18qr_pass3") ||
      fPeriod.EqualTo("LHC16qt") ||
      fPeriod.EqualTo("LHC17n") ||
      fPeriod.EqualTo("LHC15oKatarina")) {
    // return false;
  } else {
    // if(fAOD->IsPileupFromSPDInMultBins() ) { return false; }

    // AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    // if (!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return false; }

    // if(!multSelection->GetThisEventIsNotPileup() || !multSelection->GetThisEventIsNotPileupInMultBins() || !multSelection->GetThisEventHasNoInconsistentVertices() || !multSelection->GetThisEventPassesTrackletVsCluster()) { return false; }

    // Int_t nTracksPrim = fAOD->GetPrimaryVertex()->GetNContributors();
    // if(nTracksPrim < 0.5) { return false; }
  }

  // Vertex Z
  const Double_t aodVtxZ = vtx->GetZ();
  if(TMath::Abs(aodVtxZ) > 10) return kFALSE;


  // FMD cut should be checked **before** the Preparation of FMD tracks
  if(fUseFMDcut){

    double nFMD_fwd_hits = 0;
    double nFMD_bwd_hits = 0;
    AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
    const TH2D& d2Ndetadphi = aodForward->GetHistogram();
    int nEta = d2Ndetadphi.GetXaxis()->GetNbins();
    int nPhi = d2Ndetadphi.GetYaxis()->GetNbins();

    for (int iEta = 1; iEta <= nEta; iEta++) {
      for (int iPhi = 1; iPhi <= nPhi; iPhi++) {
        double eta = d2Ndetadphi.GetXaxis()->GetBinCenter(iEta); 
        double phi = d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi); 

        double etaAccepted = false;
        if ( fFMDAacceptanceCutLower < eta && eta < fFMDAacceptanceCutUpper) etaAccepted = true;
        if ( fFMDCacceptanceCutLower < eta && eta < fFMDCacceptanceCutUpper) etaAccepted = true;
        if (!etaAccepted) continue;

        double mostProbableN = d2Ndetadphi.GetBinContent(iEta, iPhi);
        if (mostProbableN > 0) {
          if (eta > 0) nFMD_fwd_hits += mostProbableN;
          else nFMD_bwd_hits += mostProbableN;
        }
      }
    }
    if(nFMD_fwd_hits==0 || nFMD_bwd_hits==0) {
      return kFALSE;
    }

    AliAODVZERO *fvzero = fAOD->GetVZEROData();
    if(!fvzero) { AliError("Problem with VZEROData, terminating!"); return kFALSE; }
    Float_t nV0A_hits = fvzero->GetMTotV0A();
    Float_t nV0C_hits = fvzero->GetMTotV0C();

    // cout << nV0A_hits << endl;
    // cout << nV0C_hits << endl;
    // cout << nFMD_fwd_hits << endl;
    // cout << nFMD_bwd_hits << endl;

    if((nV0A_hits<(fFMDcutapar0*nFMD_fwd_hits-fFMDcutapar1)) || (nV0C_hits<(fFMDcutcpar0*nFMD_bwd_hits-fFMDcutcpar1))){
      return kFALSE;
    }
    hFMDAvsV0->Fill(nFMD_fwd_hits, nV0A_hits);
    hFMDCvsV0->Fill(nFMD_bwd_hits, nV0C_hits);
  }


  bootstrap_value = (((int)(aodVtxZ * 233)) % 10 + 10) % 10;
  if (bootstrap_value < sampleLow || bootstrap_value >= sampleHigh) return false;

  return kTRUE;
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::AcceptMCTruthTrack(AliAODMCParticle *mtrk) {
  // Pt cut
  if(mtrk->Pt() < fMinPt) return kFALSE;
  if(mtrk->Pt() > fMaxPt) return kFALSE;

  // if(TMath::Abs(mtrk->Eta()) > fEtaCut) return kFALSE;

  if (!(mtrk->IsPhysicalPrimary())) return kFALSE;
  if (mtrk->Charge() == 0) return kFALSE;
  return kTRUE;
}



Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadWeightsSystematics() {
  hWeight2D = (TH2D*)fFlowWeightsList->FindObject("hNUA");
  return kTRUE;
}

Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadPtWeights() {
  int EvFlag = 0, TrFlag = 0;

  // If the period is **NOT** pPb LHC16qt
  if ( !fPeriod.EqualTo("LHC16qt") &&
       !fPeriod.EqualTo("LHC16") && !fPeriod.EqualTo("LHC17") && !fPeriod.EqualTo("LHC18") &&
       !fPeriod.EqualTo("LHC16Preview") && !fPeriod.EqualTo("LHC17Preview") && !fPeriod.EqualTo("LHC18Preview")
       ) {
    if(fCurrSystFlag == 0) fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("EffRescaled_Cent0"));
    else fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("EffRescaled_Cent0_SystFlag%i_", fCurrSystFlag));
    if(!fPtWeightsSystematics)
      {
        printf("PtWeights could not be found in list!\n");
        return kFALSE;
      }
  } 
  // If it is the pPb LHC16qt 
  else {
    if (fCurrSystFlag == 0) EvFlag = 0, TrFlag = 0;
    if (fCurrSystFlag == 1) EvFlag = 0, TrFlag = 1;
    if (fCurrSystFlag == 2) EvFlag = 0, TrFlag = 3;
    if (fCurrSystFlag == 3) EvFlag = 0, TrFlag = 0; // Abandoned
    if (fCurrSystFlag == 4) EvFlag = 0, TrFlag = 2;
    if (fCurrSystFlag == 5) EvFlag = 0, TrFlag = 3;
    if (fCurrSystFlag == 6) EvFlag = 0, TrFlag = 8;
    if (fCurrSystFlag == 7) EvFlag = 0, TrFlag = 9;
    if (fCurrSystFlag == 8) EvFlag = 0, TrFlag = 10;

    if (fCurrSystFlag == 17) EvFlag = 1, TrFlag = 0;
    if (fCurrSystFlag == 18) EvFlag = 2, TrFlag = 0;
    if (fCurrSystFlag == 19) EvFlag = 3, TrFlag = 0;

    if (fPeriod.EqualTo("LHC16qt")) {
      fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("LHC17f2b_ch_Eta_0020_Ev%d_Tr%d", EvFlag, TrFlag));

      cout << "Trying to load" << Form("LHC17f2b_ch_Eta_0020_Ev%d_Tr%d", EvFlag, TrFlag) << endl;
      fPtWeightsFeeddown    = (TH1D*)fFlowFeeddownList->FindObject(Form("LHC17f2b_ch_Eta_0020_Ev%d_Tr%d", EvFlag, TrFlag));
      if(!fPtWeightsSystematics)
        {
          printf("pPb: PtWeights could not be found in list!\n");
          return kFALSE;
        }
    } else {
      std::string period = ReturnPPperiodMC(fAOD->GetRunNumber());
      fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("LHC%s_ch_Ev%d_Tr%d", period.c_str(), EvFlag, TrFlag));
      fPtWeightsFeeddown    = (TH1D*)fFlowFeeddownList->FindObject(Form("LHC%s_ch_Ev%d_Tr%d", period.c_str(), EvFlag, TrFlag));
      if(!fPtWeightsSystematics)
        {
          printf("pp: PtWeights could not be found in list!\n");
          return kFALSE;
        }
    }
  }
  return kTRUE;

}

double AliAnalysisTaskCorrForNonlinearFlow::GetWeightKatarina(double phi, double eta, double vz) {
  double weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
                                               hPhiWeightRun->GetYaxis()->FindBin(eta),
                                               hPhiWeightRun->GetZaxis()->FindBin(vz));
  return weight;
}

// Load Katarina's weights
Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadWeightsKatarina() {
  hPhiWeightRun = (TH3F*)fFlowWeightsList->FindObject(Form("fPhiWeight_%0.lf", (double)(fAOD->GetRunNumber())));
  if (!hPhiWeightRun) {
    printf("Weights could not be found in list!\n");
    return kFALSE;
  }
  return kTRUE;
}

double AliAnalysisTaskCorrForNonlinearFlow::GetPtWeightKatarina(double pt, double eta, double vz)
{
  double weight = 1;
  double binPt = hTrackEfficiencyRun->GetXaxis()->FindBin(pt);
  double binEta = hTrackEfficiencyRun->GetYaxis()->FindBin(eta);
  double binVz = hTrackEfficiencyRun->GetZaxis()->FindBin(vz);
  //..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
  double eff = hTrackEfficiencyRun->GetBinContent(binPt, binEta, binVz);
  double error = hTrackEfficiencyRun->GetBinError(binPt, binEta, binVz);

  if((eff < 0.03) || ((error/eff) > 0.1)) weight = 1;
  else {
    TRandom3 r(0);
    double efficiency = 0;
    efficiency = r.Gaus(eff, error);
    weight = 1./efficiency; //..taking into account errors
    //weight = 1./eff;
  }
  return weight;
}

// Load Katarina's pt weights
Bool_t AliAnalysisTaskCorrForNonlinearFlow::LoadPtWeightsKatarina() {
  hTrackEfficiencyRun = (TH3F*)fFlowPtWeightsList->FindObject(Form("eff_LHC15o_HIJING_%.0lf", (double)(fAOD->GetRunNumber())));
  if (!hTrackEfficiencyRun) {
    printf("Pt Weights could not be found in list!\n");
    return kFALSE;
  }
  return kTRUE;
}

Double_t AliAnalysisTaskCorrForNonlinearFlow::GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species) {

  double dPhi = track->Phi();
  double dEta = track->Eta();
  double dVz = fVtxZ;
  double dWeight = 1.0;
  int nEta = hWeight2D->GetXaxis()->FindBin(dEta);
  int nPhi = hWeight2D->GetYaxis()->FindBin(dPhi);
  dWeight = hWeight2D->GetBinContent(nEta,nPhi);
  // cout << nEta << " " << nPhi << " " << dWeight << endl;
  return dWeight;
}


Double_t AliAnalysisTaskCorrForNonlinearFlow::GetFlowWeight(const AliVParticle* track, double fVtxZ, const PartSpecies species) {
  // if not applying for reconstructed
  // if(!fFlowWeightsApplyForReco && HasMass(species)) { return 1.0; }

  Double_t dWeight = 1.0;
  if(fFlowUse3Dweights) {
    Int_t iBin = fh3Weights[species]->FindFixBin(track->Phi(),track->Eta(),fVtxZ);
    dWeight = fh3Weights[species]->GetBinContent(iBin);
  } else {
    Int_t iBin = fh2Weights[species]->FindFixBin(track->Phi(),track->Eta());
    dWeight = fh2Weights[species]->GetBinContent(iBin);
  }

  if(dWeight <= 0.0) { dWeight = 1.0; }
  return dWeight;
}

const char* AliAnalysisTaskCorrForNonlinearFlow::GetSpeciesName(const PartSpecies species) const {
  const char* name;

  switch(species) {
  case kRefs: name = "Refs"; break;
  case kCharged: name = "Charged"; break;
  case kPion: name = "Pion"; break;
  case kKaon: name = "Kaon"; break;
  case kProton: name = "Proton"; break;
  case kCharUnidentified: name = "UnidentifiedCharged"; break;
  case kK0s: name = "K0s"; break;
  case kLambda: name = "Lambda"; break;
  case kPhi: name = "Phi"; break;
  default: name = "Unknown";
  }

  return name;
}

void AliAnalysisTaskCorrForNonlinearFlow::NTracksCalculation(AliVEvent* aod) {
  const int nAODTracks = aod->GetNumberOfTracks();
  NtrksCounter = 0;
  NTracksCorrected = 0;
  NTracksUncorrected = 0;

  //..for DCA
  double pos[3], vz, vx, vy;
  vz = aod->GetPrimaryVertex()->GetZ();
  vx = aod->GetPrimaryVertex()->GetX();
  vy = aod->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};
  float fVtxZ = vz;
  double runNumber = fInputEvent->GetRunNumber();

  //..LOOP OVER TRACKS........
  //........................................
  for(Int_t nt = 0; nt < nAODTracks; nt++)
    {
      AliAODTrack *aodTrk = (AliAODTrack*) fInputEvent->GetTrack(nt);

      if (!aodTrk) {
        continue;
      }

      aodTrk->GetXYZ(pos);
      if (!AcceptAODTrack(aodTrk, pos, vtxp)) continue;
      if(TMath::Abs(aodTrk->Eta()) > fEtaCut) continue;

      double weightPt = 1;
      if (fPeriod.EqualTo("LHC15oKatarina") ) {
        if(fNUE == 1) weightPt = GetPtWeightKatarina(aodTrk->Pt(), aodTrk->Eta(), fVtxZ);
      } else {
        // This is to make sure the extreme efficiency is only applied when calculating the Qs.
        double extremeEfficiency = fExtremeEfficiency;
        fExtremeEfficiency = 0;
        if(fNUE == 1) weightPt = GetPtWeight(aodTrk->Pt(), aodTrk->Eta(), fVtxZ, runNumber);
        fExtremeEfficiency = extremeEfficiency;
      }

      NTracksUncorrected += 1;
      NTracksCorrected += weightPt;
    } // end loop of all track
  if (!fUseCorrectedNTracks) {
    NtrksCounter = NTracksUncorrected;
  } else {
    NtrksCounter = NTracksCorrected; 
  }
  // hTracksCorrection2d->Fill(NTracksUncorrected, NTracksCorrected);
  // hnCorrectedTracks->Fill(NtrksCounter, NTracksCorrected);
}

//____________________________________________________________________
double AliAnalysisTaskCorrForNonlinearFlow::GetPtWeight(double pt, double eta, float vz, double runNumber)
{
  double binPt = fPtWeightsSystematics->GetXaxis()->FindBin(pt);
  double eff = fPtWeightsSystematics->GetBinContent(binPt);
  double error = fPtWeightsSystematics->GetBinError(binPt);
  double weight = 1;
  //..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
  if((eff < 0.03) || ((error/eff) > 0.1)) error = 0.00001;
  if((eff < 0.03)) return 1;

  TRandom3 r(0);
  double efficiency = 0;
  efficiency = r.Gaus(eff, error);
  weight = 1./efficiency; //..taking into account errors
  //weight = 1./eff;

  if (fPeriod.EqualTo("LHC16qt") ||
      fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") ||
      fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview") 
      ) {
    double binPt = fPtWeightsFeeddown->GetXaxis()->FindBin(pt);
    double feeddown = fPtWeightsFeeddown->GetBinContent(binPt);
    weight /= feeddown;
  }
  if (fExtremeEfficiency == 1) {
    // Soft: Lower region: higher efficiency, lower weight
    if (pt < 1.5) return weight * 0.98; 
    if (pt > 1.5) return weight * 1.02; 
  } else if (fExtremeEfficiency == 2) {
    if (pt < 1.5) return weight * 1.02;
    if (pt > 1.5) return weight * 0.98;
  } else if (fExtremeEfficiency == 3) {
    if (pt < 1.5) return weight * 0.96; 
    if (pt > 1.5) return weight * 1.04; 
  } else if (fExtremeEfficiency == 4) {
    if (pt < 1.5) return weight * 1.04;
    if (pt > 1.5) return weight * 0.96;
  }
  return weight;

}
//
//____________________________________________________________________
double AliAnalysisTaskCorrForNonlinearFlow::GetWeight(double phi, double eta, double pt, int fRun, bool fPlus, double vz, double runNumber) {
  TList* weights_list = dynamic_cast<TList*>(fPhiWeight);
  // cout << "weights_list" << weights_list << endl;
  // weights_list->ls();

  TList* averaged_list = dynamic_cast<TList*>(weights_list->FindObject("averaged"));
  // cout << "averaged_list" << averaged_list << endl;
  TH2D* hPhiWeightRun = dynamic_cast<TH2D*>(averaged_list->FindObject("Charged"));
  // cout << "hist_list" << hPhiWeightRun << endl;

  double weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
                                               hPhiWeightRun->GetYaxis()->FindBin(eta));
  // , hPhiWeightRun->GetZaxis()->FindBin(vz));
  return weight;
}

const char* AliAnalysisTaskCorrForNonlinearFlow::ReturnPPperiodMC(const Int_t runNumber) const
{
  if(runNumber >= 252235 && runNumber <= 264347){ // LHC16
    if(runNumber >= 252235 && runNumber <= 252375) return "17f6";
    if(runNumber >= 253437 && runNumber <= 253591) return "17f9";
    if(runNumber >= 254128 && runNumber <= 254332) return "17d17";
    if(runNumber >= 254604 && runNumber <= 255467) return "17f5";
    if(runNumber >= 255539 && runNumber <= 255618) return "17d3";
    if(runNumber >= 256219 && runNumber <= 256418) return "17e5";
    if(runNumber >= 256941 && runNumber <= 258537) return "18f1";
    if(runNumber >= 258962 && runNumber <= 259888) return "18d8";
    if(runNumber >= 262424 && runNumber <= 264035) return "17d16";
    if(runNumber >= 264076 && runNumber <= 264347) return "17d18";
  }

  if(runNumber >= 270581 && runNumber <= 282704){ // LHC17
    if(runNumber >= 270581 && runNumber <= 270667) return "18d3";
    if(runNumber >= 270822 && runNumber <= 270830) return "17h1";
    if(runNumber >= 270854 && runNumber <= 270865) return "18d3";
    if(runNumber >= 271870 && runNumber <= 273103) return "18c12";
    if(runNumber >= 273591 && runNumber <= 274442) return "17k4";
    if(runNumber >= 274593 && runNumber <= 274671) return "17h11";
    if(runNumber >= 274690 && runNumber <= 276508) return "18c13";
    if(runNumber >= 276551 && runNumber <= 278216) return "18a8";
    if(runNumber >= 278914 && runNumber <= 280140) return "17l5";
    if(runNumber >= 280282 && runNumber <= 281961) return "18a9";
    if(runNumber >= 282528 && runNumber <= 282704) return "18a1";
  }

  if(runNumber >= 285009 && runNumber <= 294925){ // LHC18
    if(runNumber >= 285009 && runNumber <= 285396) return "18g4";
    if(runNumber >= 285978 && runNumber <= 286350) return "18g5";
    if(runNumber >= 286380 && runNumber <= 286937) return "18g6";
    if(runNumber >= 287000 && runNumber <= 287658) return "18h2";
    if(runNumber >= 288619 && runNumber <= 289201) return "18h4"; //g,h,i,j,k
    if(runNumber >= 289240 && runNumber <= 289971) return "18j1";
    if(runNumber >= 290323 && runNumber <= 292839) return "18j4";
    if(runNumber >= 293357 && runNumber <= 293359) return "18k1";
    if(runNumber >= 293475 && runNumber <= 293898) return "18k2";
    if(runNumber >= 294009 && runNumber <= 294925) return "18k3";
  }

  AliWarning("PP period identifier was called and based on the run number did not pick up the correct efficiency. Setting up efficiencies from LHC18j4.");
  return "18j4";
}
