/* -------------------------------------------
 * Maintainer: Mingrui Zhao
 */
#include "AliAnalysisTaskVnPtCorr.h"
#include "AliGFWMCuts.h"
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
ClassImp(AliAnalysisTaskVnPtCorr)
  //___________________________________________________________________________
  AliAnalysisTaskVnPtCorr::AliAnalysisTaskVnPtCorr():
    AliAnalysisTaskSE(),
    fEventCuts(),
    fGFWSelection(NULL),
    fGFWSelection15o(NULL),
    fAOD(0),
    fEtaCut(0.8),
    fVtxCut(10.0),
    fMinPt(0.2),
    fMaxPt(3.0),
    fTrigger(0),
    fAliTrigger(0),
    fNUE(0),
    fNUA(0),
    fIsMC(0),
    fNtrksName("Mult"),
    //....
    fPeriod("LHC15o"),
    fCurrSystFlag(0),
    fSpringMode(false),
    fAddTPCPileupCuts(false),
    fESDvsTPConlyLinearCut(15000.),
    fUseOutOfBunchPileupCut(false),
    fUseCorrectedNTracks(false),
    binning_factor(1),
    fCentralityCut(101.1),
    fUseNarrowBin(false),
    fExtremeEfficiency(0),
    fTPCchi2perCluster(4.0),
    fUseAdditionalDCACut(false),
    fUseDefaultWeight(false),
    bUseLikeSign(false),
    iSign(0),
    fExtendV0MAcceptance(false),
    fV0MRatioCut(0),
    fEtaGap3Sub1(0.4),
    fEtaGap3Sub2(0.4),
    fOnTheFly(false),

    fListOfObjects(0),
    fListOfProfile(0),


    hTrackEfficiencyRun(0),

    fFlowWeightsList(nullptr),
    fFlowPtWeightsList(nullptr),
    fFlowFeeddownList(nullptr),

    fPhiWeight(0),
    fWeightsSystematics(0),
    fPtWeightsSystematics(0),
    hPhiWeightRun(0),

    hEventCount(0),
    hMult(0),
    fVtxAfterCuts(0),
    fCentralityDis(0),
    fV0CentralityDis(0),
    fV0MMultiplicity(0),
    fV0MRatio(0),

    fPhiDis1D(0),
    fPhiDis(0),
    fEtaDis(0),
    fEtaBefore(0),
    fPtDis(0),
    fPtBefore(0),
    hDCAxyBefore(0),
    hDCAzBefore(0),
    hITSclustersBefore(0),
    hChi2Before(0),
    hnTPCClu(0),
    hDCAxy(0),
    hDCAz(0),
    hITSclusters(0),
    hChi2(0),
    hTracksCorrection2d(0),
    hnCorrectedTracks(0),
    multProfile(),
    correlator(),
    rand(32213)
{
  for (int i = 0; i < 30; i++) fListOfProfiles[i] = NULL;
}
//______________________________________________________________________________
AliAnalysisTaskVnPtCorr::AliAnalysisTaskVnPtCorr(const char *name, int _fNUA, int _fNUE, TString _fPeriod):
  AliAnalysisTaskSE(name),
  fEventCuts(),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fAOD(0),
  fEtaCut(0.8),
  fVtxCut(10.0),
  fMinPt(0.2),
  fMaxPt(3.0),
  fTrigger(0),
  fAliTrigger(0),
  fNUE(_fNUE),
  fNUA(_fNUA),
  fIsMC(0),
  fNtrksName("Mult"),
  //....
  fPeriod(_fPeriod),
  fCurrSystFlag(0),
  fSpringMode(false),
  fAddTPCPileupCuts(false),
  fESDvsTPConlyLinearCut(15000.),
  fUseOutOfBunchPileupCut(false),
  fUseCorrectedNTracks(false),
  binning_factor(1),
  fCentralityCut(101.1),
  fUseNarrowBin(false),
  fExtremeEfficiency(0),
  fTPCchi2perCluster(4.0),
  fUseAdditionalDCACut(false),
  fUseDefaultWeight(false),
  bUseLikeSign(false),
  iSign(0),
  fExtendV0MAcceptance(false),
  fV0MRatioCut(0),
  fEtaGap3Sub1(0.4),
  fEtaGap3Sub2(0.4),
  fOnTheFly(false),

  fListOfObjects(0),
  fListOfProfile(0),


  hTrackEfficiencyRun(0),

  fFlowWeightsList(nullptr),
  fFlowPtWeightsList(nullptr),
  fFlowFeeddownList(nullptr),

  fPhiWeight(0),
  fWeightsSystematics(0),
  fPtWeightsSystematics(0),
  hPhiWeightRun(0),
  hEventCount(0),
  hMult(0),
  fVtxAfterCuts(0),
  fCentralityDis(0),
  fV0CentralityDis(0),
  fV0MMultiplicity(0),
  fV0MRatio(0),

  fPhiDis1D(0),
  fPhiDis(0),
  fEtaDis(0),
  fEtaBefore(0),
  fPtDis(0),
  fPtBefore(0),
  hDCAxyBefore(0),
  hDCAzBefore(0),
  hITSclustersBefore(0),
  hChi2Before(0),
  hnTPCClu(0),
  hDCAxy(0),
  hDCAz(0),
  hITSclusters(0),
  hChi2(0),
  hTracksCorrection2d(0),
  hnCorrectedTracks(0),
  multProfile(),
  correlator(),
  rand(32213) {

    for (int i = 0; i < 30; i++) fListOfProfiles[i] = NULL;

    // Output slot #1 writes into a TList
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    int outputslot = 2;
    for (int i = 0; i < 30; i++) {
      outputslot++;
      DefineOutput(outputslot, TList::Class());
    }
    // DefineOutput(2, TList::Class());
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
AliAnalysisTaskVnPtCorr::AliAnalysisTaskVnPtCorr(const char *name):
  AliAnalysisTaskSE(name),
  fEventCuts(),
  fGFWSelection(NULL),
  fGFWSelection15o(NULL),
  fAOD(0),
  fEtaCut(0.8),
  fVtxCut(10.0),
  fMinPt(0.2),
  fMaxPt(3.0),
  fTrigger(0),
  fAliTrigger(0),
  fNUE(0),
  fNUA(0),
  fNtrksName("Mult"),
  //....
  fPeriod("LHC15o"),
  fCurrSystFlag(0),
  fSpringMode(false),
  fAddTPCPileupCuts(false),
  fESDvsTPConlyLinearCut(15000.),
  fUseOutOfBunchPileupCut(false),
  fUseCorrectedNTracks(false),
  binning_factor(1),
  fCentralityCut(101.1),
  fUseNarrowBin(false),
  fExtremeEfficiency(0),
  fTPCchi2perCluster(4.0),
  fUseAdditionalDCACut(false),
  fUseDefaultWeight(false),
  bUseLikeSign(false),
  iSign(0),
  fExtendV0MAcceptance(false),
  fV0MRatioCut(0),
  fEtaGap3Sub1(0.4),
  fEtaGap3Sub2(0.4),
  fOnTheFly(false),

  fListOfObjects(0),
  fListOfProfile(0),

  hTrackEfficiencyRun(0),

  fFlowWeightsList(nullptr),
  fFlowPtWeightsList(nullptr),
  fFlowFeeddownList(nullptr),

  fPhiWeight(0),

  fWeightsSystematics(0),
  fPtWeightsSystematics(0),

  hPhiWeightRun(0),

  hEventCount(0),
  hMult(0),
  fVtxAfterCuts(0),
  fCentralityDis(0),
  fV0CentralityDis(0),
  fV0MMultiplicity(0),
  fV0MRatio(0),


  fPhiDis1D(0),
  fPhiDis(0),
  fEtaDis(0),
  fEtaBefore(0),
  fPtDis(0),
  fPtBefore(0),
  hDCAxyBefore(0),
  hDCAzBefore(0),
  hITSclustersBefore(0),
  hChi2Before(0),
  hnTPCClu(0),
  hDCAxy(0),
  hDCAz(0),
  hITSclusters(0),
  hChi2(0),
  hTracksCorrection2d(0),
  hnCorrectedTracks(0),
  multProfile(),
  correlator(),
  rand(32213) {

    for (int i = 0; i < 30; i++) fListOfProfiles[i] = NULL;

    // Output slot #1 writes into a TList
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    int outputslot = 2;
    for (int i = 0; i < 30; i++) {
      outputslot++;
      DefineOutput(outputslot, TList::Class());
    }
    // DefineOutput(2, TList::Class());
    // int inputslot = 1;
    DefineInput(1, TList::Class());
    DefineInput(2, TList::Class());
  }

//_____________________________________________________________________________
AliAnalysisTaskVnPtCorr::~AliAnalysisTaskVnPtCorr()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects) delete fListOfObjects;
  if (fListOfProfile) delete fListOfProfile;
  for (int i = 0; i < 30; i++) {
    if (fListOfProfiles[i]) delete fListOfProfiles[i];
  }

  if (fGFWSelection) delete fGFWSelection;
  if (fGFWSelection15o) delete fGFWSelection15o;
}

//______________________________________________________________________________
void AliAnalysisTaskVnPtCorr::UserCreateOutputObjects()
{

  //OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  //..Settings for AliEventCuts:
  //..This adds QA plots to the output
  fEventCuts.AddQAplotsToList(fListOfObjects);
  //..kINT7 is set in the class as default, if I want to have kHigHMultV0 in pp, I have to switch to manual mode

  if (fSpringMode) {
    fEventCuts.SetManualMode();
    fEventCuts.fRequireTrackVertex = false; // !!
    fEventCuts.fMinVtz = -10.f;
    fEventCuts.fMaxVtz = 10.f;
    fEventCuts.fMaxResolutionSPDvertex = 0.25f;
    // Distance between track and SPD vertex < 0.2 cm
    fEventCuts.fPileUpCutMV = true;
  }

  if (fExtendV0MAcceptance) {
    fEventCuts.OverrideCentralityFramework(1);
    fEventCuts.SetCentralityEstimators("V0M", "CL0");
    fEventCuts.SetCentralityRange(0.f,101.f);
  }

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1 and LHC17n
    fGFWSelection15o = new AliGFWNFCuts();
    fGFWSelection15o->PrintSetup();
  } else {
    fGFWSelection = new AliGFWMCuts();
    fGFWSelection->PrintSetup();
  }

  if (fOnTheFly) {
    nn = 1000;
    for (int i = 0; i <= 1000; i++) {
      xbins[i] = 30.0/nn*i;
      xbins_uncorr[i] = 30.0/nn*i;
    }
  } else {
    if (fNtrksName == "Mult") {
      if (!fUseNarrowBin) {
        nn = 200 + 56;
        // 56 = (3000-200)/50
        for (int i = 0; i <= 200; i++) {
          xbins[i] = (i + 0.5)*binning_factor;
          xbins_uncorr[i] = (i + 0.5);
        }
        for (int i = 1; i <= 56; i++) {
          xbins[200+i] = (50*i + 200 + 0.5)*binning_factor;
          xbins_uncorr[200+i] = (50*i + 200 + 0.5);
        }
      } else {
        nn = 3000;
        for (int i = 0; i <= 3000; i++) {
          xbins[i] = i;
          xbins_uncorr[i] = i;
        }
      }
    } else {
      nn = 100;
      for (int i = 0; i <= 100; i++) {
        xbins[i] = i;
        xbins_uncorr[i] = i;
      }
    }
  }

  hEventCount = new TH1D("hEventCount", "; centrality;;", 5, 0.5, 5.5);
  fListOfObjects->Add(hEventCount);

  hMult = new TH1F("hMult", ";number of tracks; entries", nn, xbins);
  hMult->Sumw2();
  fListOfObjects->Add(hMult);

  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fVtxAfterCuts->Sumw2();
  fListOfObjects->Add(fVtxAfterCuts);

  fCentralityDis = new TH1F("fCentralityDis", "centrality distribution; centrality; Counts", 101, 0, 101);
  fListOfObjects->Add(fCentralityDis);

  fV0CentralityDis = new TH1F("fV0CentralityDis", "centrality V0/<V0> distribution; centrality; Counts", 101, 0, 101);
  fListOfObjects->Add(fV0CentralityDis);

  fV0CentralityDisNarrow = new TH1F("fV0CentralityDisNarrow", "centrality V0/<V0> distribution; centrality; Counts", 1000, 0, 10);
  fListOfObjects->Add(fV0CentralityDisNarrow);

  fV0MMultiplicity = new TH1F("fV0MMultiplicity", "V0 Multiplicity distribution", 1000, 0, 1000);
  fListOfObjects->Add(fV0MMultiplicity);

  fV0MRatio = new TH1F("fV0MRatio", "V0M / <V0M> distribution", 100, 0, 10);
  fListOfObjects->Add(fV0MRatio);


  fPhiDis1DBefore = new TH1D("hPhiDisBefore", "phi distribution before the weight correction", 60, 0, 2*3.1415926);
  fListOfObjects->Add(fPhiDis1DBefore);
  fPhiDis1D  = new TH1D("hPhiDis", "phi distribution after the weight correction", 60, 0, 2*3.1415926);
  fListOfObjects->Add(fPhiDis1D);
  fEtaDis = new TH1D("hEtaDis", "eta distribution", 100, -2, 2);
  fListOfObjects->Add(fEtaDis);
  fPtDis = new TH1D("hPtDis", "pt distribution", 100, 0, 5);
  fListOfObjects->Add(fPtDis);
  hTracksCorrection2d = new TH2D("hTracksCorrection2d", "Correlation table for number of tracks table", nn, xbins, nn, xbins_uncorr);
  fListOfObjects->Add(hTracksCorrection2d);
  hnCorrectedTracks = new TProfile("hnCorrectedTracks", "Number of corrected tracks in a ntracks bin", nn, xbins);
  fListOfObjects->Add(hnCorrectedTracks);

  hDCAxy = new TH2D("hDCAxy", "DCAxy distribution", 100, 0, 0.2, 600, 0, 3);
  fListOfObjects->Add(hDCAxy);
  hDCAz  = new TH1D("hDCAz",  "DCAz distribution", 100, 0, 5);
  fListOfObjects->Add(hDCAz);
  hDCAxyBefore = new TH2D("hDCAxyBefore", "DCAxy distribution", 100, 0, 0.2, 100, 0, 3);
  fListOfObjects->Add(hDCAxyBefore);
  hDCAzBefore  = new TH1D("hDCAzBefore",  "DCAz distribution", 100, 0, 5);
  fListOfObjects->Add(hDCAzBefore);
  hChi2  = new TH1D("hChi2", "TPC chi2 per cluster", 100, 0, 5);
  fListOfObjects->Add(hChi2);
  hChi2Before  = new TH1D("hChi2Before", "TPC chi2 per cluster", 100, 0, 5);
  fListOfObjects->Add(hChi2Before);
  hnTPCClu  = new TH1D("hnTPCClu",  "Number of TPC clusters", 100, 40, 140);
  fListOfObjects->Add(hnTPCClu);


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
      cout << "Got feeddown list" << endl;
    } else {
      fFlowPtWeightsList = (TList*) GetInputData(inSlotCounter);
    }
    inSlotCounter++;
  };

  fListOfProfile = new TList();
  fListOfProfile->SetOwner();
  for (int i = 0; i < 30; i++) {
    fListOfProfiles[i] = new TList();
    fListOfProfiles[i]->SetOwner();
  }

  // Create Q Distribution

  // Physics profiles
  //	NL response
  InitProfile(multProfile, "", fListOfProfile);
  for (int i = 0; i < 30; i++) InitProfile(multProfile_bin[i], Form("_%d", i), fListOfProfiles[i]);

  // Post output data.
  PostData(1, fListOfObjects);
  int outputslot = 2;
  PostData(2, fListOfProfile);
  for (int i = 0; i < 30; i++) {
    outputslot++;
    PostData(outputslot, fListOfProfiles[i]);
  }
}

//_________________________________________________________________
void AliAnalysisTaskVnPtCorr::NotifyRun() {
  if (fAddTPCPileupCuts) {
    Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
    fEventCuts.fUseVariablesCorrelationCuts = true;
    fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    fEventCuts.fESDvsTPConlyLinearCut[0] = fESDvsTPConlyLinearCut;
  }
}

//______________________________________________________________________________
void AliAnalysisTaskVnPtCorr::UserExec(Option_t *)
{
  // bootstrap_value = rand.Integer(30);

  // Check if it can pass the trigger
  //..apply physics selection
  if (!fOnTheFly) {
    UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isTrigselected = false;
    if (fTrigger == 0) {
      isTrigselected = fSelectMask&(AliVEvent::kINT7+AliVEvent::kMB);
      fAliTrigger = AliVEvent::kINT7+AliVEvent::kMB;
    } else if (fTrigger == 1) {
      isTrigselected = fSelectMask&AliVEvent::kHighMultV0;
      fAliTrigger = AliVEvent::kHighMultV0;
    }
    if(isTrigselected == false) return;
  }

  if (fOnTheFly) { fMCEvent = getMCEvent(); }
  else {
    //..check if I have AOD
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) {
      Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      return;
    }
  }

  if (fOnTheFly) {
    hEventCount->Fill("after fEventCuts", 1.);
    // bootstrap_value = rand.Integer(30);
  } else {
    //..standard event plots (cent. percentiles, mult-vs-percentile)
    const auto pms(static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection")));
    double dCentrality = 0;
    if (pms) {
      dCentrality = pms->GetMultiplicityPercentile("V0M");
    }
    float centrV0 = dCentrality;
    float cent = dCentrality;
    float centSPD = 0;

    fCentralityDis->Fill(cent);
    fV0CentralityDis->Fill(centrV0);
    fV0CentralityDisNarrow->Fill(centrV0);

    // Initialize the estimator
    AliMultEstimator *v0Est = 0;
    if (pms) {
      v0Est = pms->GetEstimator("V0M");
      fV0MMultiplicity->Fill(v0Est->GetValue());
      fV0MRatio->Fill(v0Est->GetValue()/v0Est->GetMean());
    }

    // Check if it passed the standard AOD selection
    if (!AcceptAOD(fAOD) ) {
      PostData(1,fListOfObjects);
      int outputslot = 2;
      PostData(2, fListOfProfile);
      for (int i = 0; i < 30; i++) {
        outputslot++;
        PostData(outputslot, fListOfProfiles[i]);
      }
      return;
    }
    hEventCount->Fill("after fEventCuts", 1.);

    if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
      fGFWSelection15o->ResetCuts();
    } else {
      fGFWSelection->ResetCuts();
    }
    //..filling Vz distribution
    AliVVertex *vtx = fAOD->GetPrimaryVertex();
    float fVtxZ = vtx->GetZ();

    if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
      if (!fGFWSelection15o->AcceptVertex(fAOD)) {
        PostData(1,fListOfObjects);
        int outputslot = 2;
        PostData(2, fListOfProfile);
        for (int i = 0; i < 30; i++) {
          outputslot++;
          PostData(outputslot, fListOfProfiles[i]);
        }
        return;
      }
    } else {
      if (!fGFWSelection->AcceptVertex(fAOD)) {
        PostData(1,fListOfObjects);
        int outputslot = 2;
        PostData(2, fListOfProfile);
        for (int i = 0; i < 30; i++) {
          outputslot++;
          PostData(outputslot, fListOfProfiles[i]);
        }
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


    NTracksCalculation(fInputEvent);

    // Setup AliGFWMCuts for a specific systematics
    if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
      fGFWSelection15o->SetupCuts(fCurrSystFlag);
      if (!fGFWSelection15o->AcceptVertex(fAOD)) {
        PostData(1,fListOfObjects);
        int outputslot = 2;
        PostData(2, fListOfProfile);
        for (int i = 0; i < 30; i++) {
          outputslot++;
          PostData(outputslot, fListOfProfiles[i]);
        }
        return;
      }
    } else {
      fGFWSelection->SetupCuts(fCurrSystFlag);
      if (!fGFWSelection->AcceptVertex(fAOD)) {
        PostData(1,fListOfObjects);
        int outputslot = 2;
        PostData(2, fListOfProfile);
        for (int i = 0; i < 30; i++) {
          outputslot++;
          PostData(outputslot, fListOfProfiles[i]);
        }
        return;
      }
    }
    // Check the VtxZ distribution
    fVtxAfterCuts->Fill(fVtxZ);

    hMult->Fill(NtrksCounter);


    if (cent > fCentralityCut) {
      PostData(1,fListOfObjects);
      int outputslot = 2;
      PostData(2, fListOfProfile);
      for (int i = 0; i < 30; i++) {
        outputslot++;
        PostData(outputslot, fListOfProfiles[i]);
      }
      return;
    }

    //..all charged particles
    if (!fIsMC) {
      AnalyzeAOD(fInputEvent, centrV0, cent, centSPD, fVtxZ, false);
    } else {
      AnalyzeMCTruth(fInputEvent, centrV0, cent, centSPD, fVtxZ, false);
    }
  }

  if (fOnTheFly) {
    AnalyzeMCOnTheFly(fMCEvent);
  }

  // Post output data.
  PostData(1, fListOfObjects);
  int outputslot = 2;
  PostData(2, fListOfProfile);
  for (int i = 0; i < 30; i++) {
    outputslot++;
    PostData(outputslot, fListOfProfiles[i]);
  }

}

void AliAnalysisTaskVnPtCorr::NTracksCalculation(AliVEvent* aod) {
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

    /*
    //..get phi-weight for NUA correction
    double weight = 1;
    if (fPeriod.EqualTo("LHC15oKatarina") ) {
    if(fNUA == 1) weight = GetWeightKatarina(aodTrk->Phi(), aodTrk->Eta(), fVtxZ);
    } else {
    if(fNUA == 1) weight = GetFlowWeightSystematics(aodTrk, fVtxZ, kRefs);
    }
    */
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
  hTracksCorrection2d->Fill(NTracksUncorrected, NTracksCorrected);
  hnCorrectedTracks->Fill(NtrksCounter, NTracksCorrected);
}

//________________________________________________________________________
void AliAnalysisTaskVnPtCorr::AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus)
{

  const int nAODTracks = aod->GetNumberOfTracks();

  // Init the number of tracks
  NtrksAfter = 0;
  NtrksAfterGap0M = 0;
  NtrksAfterGap0P = 0;
  NtrksAfterGap2M = 0;
  NtrksAfterGap2P = 0;
  NtrksAfterGap4M = 0;
  NtrksAfterGap4P = 0;
  NtrksAfterGap6M = 0;
  NtrksAfterGap6P = 0;
  NtrksAfterGap8M = 0;
  NtrksAfterGap8P = 0;
  NtrksAfterGap10M = 0;
  NtrksAfterGap10P = 0;
  NtrksAfterGap14M = 0;
  NtrksAfterGap14P = 0;
  NtrksAfter3subL = 0;
  NtrksAfter3subM = 0;
  NtrksAfter3subR = 0;

  sumPtw = 0;
  sumPtw2 = 0;
  sumPt2w2 = 0;
  sumWeight = 0;
  sumWeight2 = 0;
  eventWeight  = 0;
  eventWeight2 = 0;

  //..for DCA
  double pos[3], vz, vx, vy;
  vz = aod->GetPrimaryVertex()->GetZ();
  vx = aod->GetPrimaryVertex()->GetX();
  vy = aod->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};

  double QcosGap8M[20][20] = {0};
  double QsinGap8M[20][20] = {0};
  double QcosGap8P[20][20] = {0};
  double QsinGap8P[20][20] = {0};
  double QcosSubLeft[20][20] = {0};
  double QsinSubLeft[20][20] = {0};
  double QcosSubMiddle[20][20] = {0};
  double QsinSubMiddle[20][20] = {0};
  double QcosSubRight[20][20] = {0};
  double QsinSubRight[20][20] = {0};


  double runNumber = fInputEvent->GetRunNumber();

  //..LOOP OVER TRACKS........
  //........................................
  for(Int_t nt = 0; nt < nAODTracks; nt++) {

    AliAODTrack *aodTrk = (AliAODTrack*) fInputEvent->GetTrack(nt);

    if (!aodTrk) {
      continue;
    }

    aodTrk->GetXYZ(pos);
    double dcaX = pos[0] - vtxp[0];
    double dcaY = pos[1] - vtxp[1];
    double dcaZ = abs(pos[2] - vtxp[2]);
    double dcaXY = TMath::Sqrt(dcaX*dcaX+dcaY*dcaY);

    double fb = (fCurrSystFlag == 1) ? 768 : 96;
    if (aodTrk->TestFilterBit(fb)) {
      hDCAxyBefore->Fill(dcaXY, aodTrk->Pt());
      hDCAzBefore->Fill(dcaZ);
      hChi2Before->Fill(aodTrk->GetTPCchi2perCluster());
    }

    if (!AcceptAODTrack(aodTrk, pos, vtxp)) continue;
    if (fUseAdditionalDCACut) {
      if (dcaXY > 1) continue;
      if (dcaZ > 1) continue;
    }
    if(bUseLikeSign)
    {
      if(!(aodTrk->Charge() == iSign)) continue;
    }

    hDCAxy->Fill(dcaXY, aodTrk->Pt());
    hDCAz->Fill(dcaZ);
    hChi2->Fill(aodTrk->GetTPCchi2perCluster());
    hnTPCClu->Fill(aodTrk->GetTPCNclsF());
    NtrksAfter += 1;

    //..get phi-weight for NUA correction
    double weight = 1;
    if (fPeriod.EqualTo("LHC15oKatarina") ) {
      if(fNUA == 1) weight = GetWeightKatarina(aodTrk->Phi(), aodTrk->Eta(), fVtxZ);
    } else {
      if(fNUA == 1) weight = GetFlowWeightSystematics(aodTrk, fVtxZ, kRefs);
    }
    double weightPt = 1;
    if (fPeriod.EqualTo("LHC15oKatarina") ) {
      if(fNUE == 1) weightPt = GetPtWeightKatarina(aodTrk->Pt(), aodTrk->Eta(), fVtxZ);
    } else {
      if(fNUE == 1) weightPt = GetPtWeight(aodTrk->Pt(), aodTrk->Eta(), fVtxZ, runNumber);
    }

    fPhiDis1DBefore->Fill(aodTrk->Phi());
    fPtDis->Fill(aodTrk->Pt());
    fEtaDis->Fill(aodTrk->Eta());
    fPhiDis1D->Fill(aodTrk->Phi(), weight*weightPt);



    //..Gap > 0.8
    if(aodTrk->Eta() < -0.4) {
      NtrksAfterGap8M++;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          QsinGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
        }
      }
    }
    if(aodTrk->Eta() > 0.4) {
      NtrksAfterGap8P++;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          QsinGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
        }
      }
    }

    //..3-subevent method
    if(aodTrk->Eta() < -fEtaGap3Sub1) {//..left part
      NtrksAfter3subL += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
        }
      }
    }
    if(aodTrk->Eta() >= -fEtaGap3Sub2 && aodTrk->Eta() <= fEtaGap3Sub2) {//..middle part
                                                                         // eventWeight += weightPt;
      sumPtw+=weightPt*aodTrk->Pt();
      sumPtw2+=weightPt*weightPt*aodTrk->Pt();
      sumPt2w2 += weightPt*weightPt*aodTrk->Pt()*aodTrk->Pt();
      sumWeight += weightPt;
      sumWeight2 += weightPt*weightPt;

      NtrksAfter3subM += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
        }
      }
    }
    if(aodTrk->Eta() > fEtaGap3Sub1) {//..right part
      NtrksAfter3subR += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
        }
      }
    }
  } // end loop of all track

  //............................
  //..GENERIC FRAMEWORK RP
  //............................

  //..calculate Q-vector for each harmonics n and power p
  correlator.FillQVector(correlator.Qvector8M, QcosGap8M, QsinGap8M);
  correlator.FillQVector(correlator.Qvector8P, QcosGap8P, QsinGap8P);
  correlator.FillQVector(correlator.QvectorSubLeft, QcosSubLeft, QsinSubLeft);
  correlator.FillQVector(correlator.QvectorSubRight, QcosSubRight, QsinSubRight);
  correlator.FillQVector(correlator.QvectorSubMiddle, QcosSubMiddle, QsinSubMiddle);

  if (fNtrksName == "Mult") {
    CalculateProfile(multProfile, NtrksCounter, NTracksUncorrected);
    CalculateProfile(multProfile_bin[bootstrap_value], NtrksCounter, NTracksUncorrected);
  }
}

//________________________________________________________________________
void AliAnalysisTaskVnPtCorr::AnalyzeMCTruth(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus)
{

  TClonesArray* farray = (TClonesArray*)aod->FindListObject("mcparticles");
  const int nAODTracks = farray->GetEntries();
  // AliAODMCParticle *trk = (AliAODMCParticle*) farray->At(TMath::Abs(label));

  // Init the number of tracks
  NtrksAfter = 0;
  NtrksAfterGap0M = 0;
  NtrksAfterGap0P = 0;
  NtrksAfterGap2M = 0;
  NtrksAfterGap2P = 0;
  NtrksAfterGap4M = 0;
  NtrksAfterGap4P = 0;
  NtrksAfterGap6M = 0;
  NtrksAfterGap6P = 0;
  NtrksAfterGap8M = 0;
  NtrksAfterGap8P = 0;
  NtrksAfterGap10M = 0;
  NtrksAfterGap10P = 0;
  NtrksAfterGap14M = 0;
  NtrksAfterGap14P = 0;
  NtrksAfter3subL = 0;
  NtrksAfter3subM = 0;
  NtrksAfter3subR = 0;

  sumPtw = 0;
  sumPtw2 = 0;
  sumPt2w2 = 0;
  sumWeight = 0;
  sumWeight2 = 0;
  eventWeight  = 0;
  eventWeight2 = 0;


  //..for DCA
  // double pos[3], vz, vx, vy;
  // vz = aod->GetPrimaryVertex()->GetZ();
  // vx = aod->GetPrimaryVertex()->GetX();
  // vy = aod->GetPrimaryVertex()->GetY();
  // double vtxp[3] = {vx, vy, vz};
  // Assume that DCA cuts not needed here

  double QcosGap8M[20][20] = {0};
  double QsinGap8M[20][20] = {0};
  double QcosGap8P[20][20] = {0};
  double QsinGap8P[20][20] = {0};
  double QcosSubLeft[20][20] = {0};
  double QsinSubLeft[20][20] = {0};
  double QcosSubMiddle[20][20] = {0};
  double QsinSubMiddle[20][20] = {0};
  double QcosSubRight[20][20] = {0};
  double QsinSubRight[20][20] = {0};



  // double runNumber = fInputEvent->GetRunNumber();
  // Weight is not needed

  //..LOOP OVER TRACKS........
  //........................................
  for(Int_t nt = 0; nt < nAODTracks; nt++) {

    AliAODMCParticle *track = (AliAODMCParticle*) farray->At(TMath::Abs(nt));

    if (!track) {
      continue;
    }

    if (fUseOutOfBunchPileupCut && AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(nt, fMCEvent)) continue;
    // track->GetXYZ(pos);
    if (!AcceptMCTruthTrack(track)) continue;


    NtrksAfter += 1;

    //..get phi-weight for NUA correction
    double weight = 1;
    double weightPt = 1;

    //..calculate Q-vectors
    //..Gap > 0.8
    if(track->Eta() < -0.4) {
      NtrksAfterGap8M++;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
    if(track->Eta() > 0.4) {
      NtrksAfterGap8P++;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }

    //..3-subevent method
    if(track->Eta() < -fEtaGap3Sub1) {//..left part
      NtrksAfter3subL += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
    if(track->Eta() >= -fEtaGap3Sub2 && track->Eta() <= fEtaGap3Sub2) {//..middle part
      sumPtw+=weightPt*track->Pt();
      sumPtw2+=weightPt*weightPt*track->Pt();
      sumPt2w2 += weightPt*weightPt*track->Pt()*track->Pt();
      sumWeight += weightPt;
      sumWeight2 += weightPt*weightPt;

      NtrksAfter3subM += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
    if(track->Eta() > fEtaGap3Sub1) {//..right part
      NtrksAfter3subR += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
  } // end loop of all track

  //............................
  //..GENERIC FRAMEWORK RP
  //............................
  correlator.FillQVector(correlator.Qvector8M, QcosGap8M, QsinGap8M);
  correlator.FillQVector(correlator.Qvector8P, QcosGap8P, QsinGap8P);
  correlator.FillQVector(correlator.QvectorSubLeft, QcosSubLeft, QsinSubLeft);
  correlator.FillQVector(correlator.QvectorSubRight, QcosSubRight, QsinSubRight);
  correlator.FillQVector(correlator.QvectorSubMiddle, QcosSubMiddle, QsinSubMiddle);

  if (fNtrksName == "Mult") {
    CalculateProfile(multProfile, NtrksCounter, NTracksUncorrected);
    CalculateProfile(multProfile_bin[bootstrap_value], NtrksCounter, NTracksUncorrected);
  } else {
  }

}


//________________________________________________________________________
void AliAnalysisTaskVnPtCorr::AnalyzeMCOnTheFly(AliMCEvent* aod) {

  NtrksCounter = fImpactParameterMC;
  bootstrap_value = (((int)(NtrksCounter * 233)) % 30 + 30) % 30;

  // Init the number of tracks
  NtrksAfter = 0;
  NtrksAfterGap0M = 0;
  NtrksAfterGap0P = 0;
  NtrksAfterGap2M = 0;
  NtrksAfterGap2P = 0;
  NtrksAfterGap4M = 0;
  NtrksAfterGap4P = 0;
  NtrksAfterGap6M = 0;
  NtrksAfterGap6P = 0;
  NtrksAfterGap8M = 0;
  NtrksAfterGap8P = 0;
  NtrksAfterGap10M = 0;
  NtrksAfterGap10P = 0;
  NtrksAfterGap14M = 0;
  NtrksAfterGap14P = 0;
  NtrksAfter3subL = 0;
  NtrksAfter3subM = 0;
  NtrksAfter3subR = 0;


  //..for DCA
  // double pos[3], vz, vx, vy;
  // vz = aod->GetPrimaryVertex()->GetZ();
  // vx = aod->GetPrimaryVertex()->GetX();
  // vy = aod->GetPrimaryVertex()->GetY();
  // double vtxp[3] = {vx, vy, vz};
  // Assume that DCA cuts not needed here

  double QcosGap8M[20][20] = {0};
  double QsinGap8M[20][20] = {0};
  double QcosGap8P[20][20] = {0};
  double QsinGap8P[20][20] = {0};
  double QcosSubLeft[20][20] = {0};
  double QsinSubLeft[20][20] = {0};
  double QcosSubMiddle[20][20] = {0};
  double QsinSubMiddle[20][20] = {0};
  double QcosSubRight[20][20] = {0};
  double QsinSubRight[20][20] = {0};



  // double runNumber = fInputEvent->GetRunNumber();
  // Weight is not needed

  int nAODTracks = aod->GetNumberOfPrimaries();
  //..LOOP OVER TRACKS........
  //........................................
  for(Int_t nt = 0; nt < nAODTracks; nt++) {

    AliMCParticle *track = (AliMCParticle*)(aod->GetTrack(nt));

    if (!track) {
      continue;
    }

    // track->GetXYZ(pos);
    // if (!AcceptMCTruthTrack(track)) continue;
    if(track->Pt() < fMinPt) continue;
    if(track->Pt() > fMaxPt) continue;
    if(TMath::Abs(track->Eta()) > fEtaCut) continue;
    if (!(track->IsPhysicalPrimary())) continue;
    // if (!(track->IsPhysicalPrimary())) return kFALSE;
    if (track->Charge() == 0) continue;

    NtrksAfter += 1;

    //..get phi-weight for NUA correction
    double weight = 1;
    double weightPt = 1;

    //..Gap > 0.8
    if(track->Eta() < -0.4) {
      NtrksAfterGap8M++;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
    if(track->Eta() > 0.4) {
      NtrksAfterGap8P++;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }

    //..3-subevent method
    if(track->Eta() < -fEtaGap3Sub1) {//..left part
      NtrksAfter3subL += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
    if(track->Eta() >= -fEtaGap3Sub2 && track->Eta() <= fEtaGap3Sub2) {//..middle part
      sumPtw+=weightPt*track->Pt();
      sumPtw2+=weightPt*weightPt*track->Pt();
      sumPt2w2 += weightPt*weightPt*track->Pt()*track->Pt();
      sumWeight += weightPt;
      sumWeight2 += weightPt*weightPt;

      NtrksAfter3subM += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
    if(track->Eta() > fEtaGap3Sub1) {//..right part
      NtrksAfter3subR += 1;
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
  } // end loop of all track

  //............................
  //..GENERIC FRAMEWORK RP
  //............................

  //..calculate Q-vector for each harmonics n and power p

  correlator.FillQVector(correlator.Qvector8M, QcosGap8M, QsinGap8M);
  correlator.FillQVector(correlator.Qvector8P, QcosGap8P, QsinGap8P);
  correlator.FillQVector(correlator.QvectorSubLeft, QcosSubLeft, QsinSubLeft);
  correlator.FillQVector(correlator.QvectorSubRight, QcosSubRight, QsinSubRight);
  correlator.FillQVector(correlator.QvectorSubMiddle, QcosSubMiddle, QsinSubMiddle);

  if (fNtrksName == "Mult") {
    CalculateProfile(multProfile, NtrksCounter, NTracksUncorrected);
    CalculateProfile(multProfile_bin[bootstrap_value], NtrksCounter, NTracksUncorrected);
  } else {
  }
}


//____________________________________________________________________
//	END OF MAIN PROGRAM
//____________________________________________________________________

double AliAnalysisTaskVnPtCorr::GetPtWeight(double pt, double eta, float vz, double runNumber)
{
  double binPt = 0;
  double eff = 1;
  double error = 1;
  if (fPeriod.EqualTo("LHC16qt")) {
    binPt = fEtaPtWeightsSystematics[GetEtaPtFlag(eta)]->GetXaxis()->FindBin(pt);
    eff = fEtaPtWeightsSystematics[GetEtaPtFlag(eta)]->GetBinContent(binPt);
    error = fEtaPtWeightsSystematics[GetEtaPtFlag(eta)]->GetBinError(binPt);
  } else {
    binPt = fPtWeightsSystematics->GetXaxis()->FindBin(pt);
    eff = fPtWeightsSystematics->GetBinContent(binPt);
    error = fPtWeightsSystematics->GetBinError(binPt);
  }
  double weight = 1;
  //..take into account error on efficiency: randomly get number from gaussian distribution of eff. where width = error
  if((eff < 0.03) || ((error/eff) > 0.1)) error = 0.00001;
  if((eff < 0.03)) return 1;

  TRandom3 r(0);
  double efficiency = 0;
  efficiency = r.Gaus(eff, error);
  weight = 1./efficiency; //..taking into account errors
                          //weight = 1./eff;

  if (fPeriod.EqualTo("LHC16qt")) {
    double binPt = fEtaPtWeightsFeeddown[GetEtaPtFlag(eta)]->GetXaxis()->FindBin(pt);
    double feeddown = fEtaPtWeightsFeeddown[GetEtaPtFlag(eta)]->GetBinContent(binPt);
    weight /= feeddown;

  } else if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") ||
      fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")) {
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


Bool_t AliAnalysisTaskVnPtCorr::LoadWeightsSystematics() {

  // If the period is not pPb LHC16qt
  if (! (fPeriod.EqualTo("LHC16qt") || fPeriod.EqualTo("LHC16qt_Closure")) ) {
    // Only if it is the new LHC16,17,18, We need the period NUA
    if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18") ||
        fPeriod.EqualTo("LHC16_simp") || fPeriod.EqualTo("LHC17_simp") || fPeriod.EqualTo("LHC18_simp") ||
        fPeriod.EqualTo("LHC16_Closure") || fPeriod.EqualTo("LHC17_Closure") || fPeriod.EqualTo("LHC18_Closure")
       ) {
      std::string ppperiod = ReturnPPperiod(fAOD->GetRunNumber());
      // Old code: change to new one is because DCAxy < 10 is almost no cut
      // if(fCurrSystFlag == 0) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%s_SystFlag2_", ppperiod.c_str()));
      if(fCurrSystFlag == 0) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%s", ppperiod.c_str()));
      else fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%s_SystFlag%i_", ppperiod.c_str(), fCurrSystFlag));
      if(!fWeightsSystematics)
      {
        printf("Weights could not be found in list!\n");
        return kFALSE;
      }
      fWeightsSystematics->CreateNUA();
    } else if (fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")) {

      // Old code: change to new one becuase DCAxy < 10 is almost no cut
      // if(fCurrSystFlag == 0 || fUseDefaultWeight) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i",fAOD->GetRunNumber()));
      if(fCurrSystFlag == 0 || fUseDefaultWeight) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_SystFlag2_",fAOD->GetRunNumber()));
      else if (fCurrSystFlag >= 17)
        fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_SystFlag%i_",fAOD->GetRunNumber(), fCurrSystFlag-7));
      else fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_SystFlag%i_",fAOD->GetRunNumber(), fCurrSystFlag));
      // This special control is because the track flag is 1-16, but in NUA file it is 1-9
      if(!fWeightsSystematics)
      {
        printf("Weights could not be found in list!\n");
        return kFALSE;
      }
      fWeightsSystematics->CreateNUA();
    } else {
      // if(fCurrSystFlag == 0 || fUseDefaultWeight) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i",fAOD->GetRunNumber()));
      if(fCurrSystFlag == 0 || fUseDefaultWeight) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_SystFlag2_",fAOD->GetRunNumber()));
      else fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_SystFlag%i_",fAOD->GetRunNumber(), fCurrSystFlag));
      if(!fWeightsSystematics)
      {
        printf("Weights could not be found in list!\n");
        return kFALSE;
      }
      fWeightsSystematics->CreateNUA();
    }

  }
  // If it is the pPb LHC16qt
  else {
    int EvFlag = 0, TrFlag = 0;
    if (fCurrSystFlag == 0) EvFlag = 0, TrFlag = 2; // 0->2, because in the NUA file, DCAxy=10 is no cut
    if (fCurrSystFlag == 1) EvFlag = 0, TrFlag = 1;
    if (fCurrSystFlag == 2) EvFlag = 0, TrFlag = 5;
    if (fCurrSystFlag == 3) EvFlag = 0, TrFlag = 0; // Abandoned
    if (fCurrSystFlag == 4) EvFlag = 0, TrFlag = 2;
    if (fCurrSystFlag == 5) EvFlag = 0, TrFlag = 3;
    if (fCurrSystFlag == 6) EvFlag = 0, TrFlag = 8;
    if (fCurrSystFlag == 7) EvFlag = 0, TrFlag = 9;
    if (fCurrSystFlag == 8) EvFlag = 0, TrFlag = 10;

    if (fCurrSystFlag == 17) EvFlag = 1, TrFlag = 0;
    if (fCurrSystFlag == 18) EvFlag = 2, TrFlag = 0;
    if (fCurrSystFlag == 19) EvFlag = 3, TrFlag = 0;

    if (fUseDefaultWeight) EvFlag = 0, TrFlag = 0;

    fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i_Ev%d_Tr%d",fAOD->GetRunNumber(),EvFlag,TrFlag));
    if(!fWeightsSystematics)
    {
      printf("Weights could not be found in list!\n");
      return kFALSE;
    }
    fWeightsSystematics->CreateNUA();
  }
  return kTRUE;
}

Bool_t AliAnalysisTaskVnPtCorr::LoadPtWeights() {

  if (fPeriod.EqualTo("LHC15o_simp") || fPeriod.EqualTo("LHC18q_simp") || fPeriod.EqualTo("LHC18r_simp") ||
      fPeriod.EqualTo("LHC16qt_simp") ||
      fPeriod.EqualTo("LHC16_simp") || fPeriod.EqualTo("LHC17_simp") || fPeriod.EqualTo("LHC18_simp") ||
      fPeriod.Contains("Closure")
     ) {
    if(fCurrSystFlag == 0) fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("Default"));
    else fPtWeightsSystematics = (TH1D*)fFlowPtWeightsList->FindObject(Form("Sys%i", fCurrSystFlag));

    if(!fPtWeightsSystematics)
    {
      printf("PtWeights could not be found in list!\n");
      return kFALSE;
    }

    return kTRUE;
  }

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
  // If it is the pPb LHC16qt or pp
  else {
    if (fCurrSystFlag == 0) EvFlag = 0, TrFlag = 0;
    if (fCurrSystFlag == 1) EvFlag = 0, TrFlag = 1;
    if (fCurrSystFlag == 2) EvFlag = 0, TrFlag = 5;
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

      TString etaReg[8] = {"0020", "0200", "0204", "0402", "0406", "0604", "0608", "0806"};

      for (int flag = 0; flag < 8; flag++) {
        fEtaPtWeightsSystematics[flag] = (TH1D*)fFlowPtWeightsList->FindObject(Form("LHC17f2b_ch_Eta_%s_Ev%d_Tr%d", etaReg[flag].Data(), EvFlag, TrFlag));
        fEtaPtWeightsFeeddown[flag]    = (TH1D*)fFlowFeeddownList->FindObject(Form("LHC17f2b_ch_Eta_%s_Ev%d_Tr%d", etaReg[flag].Data(), EvFlag, TrFlag));
      }
      /* Too lazy to add the check
         if(!fPtWeightsSystematics)
         {
         printf("pPb: PtWeights could not be found in list!\n");
         return kFALSE;
         }
         */
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

double AliAnalysisTaskVnPtCorr::GetWeightKatarina(double phi, double eta, double vz) {
  double weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
      hPhiWeightRun->GetYaxis()->FindBin(eta),
      hPhiWeightRun->GetZaxis()->FindBin(vz));
  return weight;
}

// Load Katarina's weights
Bool_t AliAnalysisTaskVnPtCorr::LoadWeightsKatarina() {
  hPhiWeightRun = (TH3F*)fFlowWeightsList->FindObject(Form("fPhiWeight_%0.lf", (double)(fAOD->GetRunNumber())));
  if (!hPhiWeightRun) {
    printf("Weights could not be found in list!\n");
    return kFALSE;
  }
  return kTRUE;
}

double AliAnalysisTaskVnPtCorr::GetPtWeightKatarina(double pt, double eta, double vz)
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
Bool_t AliAnalysisTaskVnPtCorr::LoadPtWeightsKatarina() {
  hTrackEfficiencyRun = (TH3F*)fFlowPtWeightsList->FindObject(Form("eff_LHC15o_HIJING_%.0lf", (double)(fAOD->GetRunNumber())));
  if (!hTrackEfficiencyRun) {
    printf("Pt Weights could not be found in list!\n");
    return kFALSE;
  }
  return kTRUE;
}

Double_t AliAnalysisTaskVnPtCorr::GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species) {

  double dPhi = track->Phi();
  double dEta = track->Eta();
  double dVz = fVtxZ;
  double dWeight = 1.0;
  dWeight = fWeightsSystematics->GetNUA(dPhi, dEta, dVz);
  return dWeight;
}


void AliAnalysisTaskVnPtCorr::InitProfile(PhysicsProfileVnPt& multProfile, TString label, TList* listOfProfile) {

  multProfile.fMeanPt = new TProfile(Form("fPt%s", label.Data()), "Mean Pt", nn, xbins);
  multProfile.fMeanPt->Sumw2();
  listOfProfile->Add(multProfile.fMeanPt);

  multProfile.fc22w = new TProfile(Form("fChc22w%s", label.Data()), "v2^2 with event weight", nn, xbins);
  multProfile.fc22w->Sumw2();
  listOfProfile->Add(multProfile.fc22w);

  multProfile.fPcc = new TProfile(Form("fV2Pt%s", label.Data()), "v2-Pt correlation", nn, xbins);
  multProfile.fPcc->Sumw2();
  listOfProfile->Add(multProfile.fPcc);

  multProfile.fc22nw = new TProfile(Form("fChc22nw%s", label.Data()), "v2^2 without event weight", nn, xbins);
  multProfile.fc22nw->Sumw2();
  listOfProfile->Add(multProfile.fc22nw);

  multProfile.fc24nw = new TProfile(Form("fChc24nw%s", label.Data()), "v2^4 with event weight", nn, xbins);
  multProfile.fc24nw->Sumw2();
  listOfProfile->Add(multProfile.fc24nw);

  multProfile.fPtVariancea = new TProfile(Form("fdPt2a%s", label.Data()), "Pt variance", nn, xbins);
  multProfile.fPtVariancea->Sumw2();
  listOfProfile->Add(multProfile.fPtVariancea);

  multProfile.fPtVarianceb = new TProfile(Form("fdPt2b%s", label.Data()), "Pt variance", nn, xbins);
  multProfile.fPtVarianceb->Sumw2();
  listOfProfile->Add(multProfile.fPtVarianceb);

  multProfile.fPtVariancec = new TProfile(Form("fdPt2c%s", label.Data()), "Pt variance", nn, xbins);
  multProfile.fPtVariancec->Sumw2();
  listOfProfile->Add(multProfile.fPtVariancec);

  multProfile.fMeanPt2D = new TProfile2D(Form("fPt2D%s", label.Data()), "Mean Pt", nn, xbins, nn, xbins_uncorr);
  multProfile.fMeanPt2D->Sumw2();
  listOfProfile->Add(multProfile.fMeanPt2D);

  multProfile.fc22w2D = new TProfile2D(Form("fChc22w2D%s", label.Data()), "v2^2 with event weight", nn, xbins, nn, xbins_uncorr);
  multProfile.fc22w2D->Sumw2();
  listOfProfile->Add(multProfile.fc22w2D);

  multProfile.fPcc2D = new TProfile2D(Form("fV2Pt2D%s", label.Data()), "v2-Pt correlation", nn, xbins, nn, xbins_uncorr);
  multProfile.fPcc2D->Sumw2();
  listOfProfile->Add(multProfile.fPcc2D);

  multProfile.fc22nw2D = new TProfile2D(Form("fChc22nw2D%s", label.Data()), "v2^2 without event weight", nn, xbins, nn, xbins_uncorr);
  multProfile.fc22nw2D->Sumw2();
  listOfProfile->Add(multProfile.fc22nw2D);

  multProfile.fc24nw2D = new TProfile2D(Form("fChc24nw2D%s", label.Data()), "v2^4 with event weight", nn, xbins, nn, xbins_uncorr);
  multProfile.fc24nw2D->Sumw2();
  listOfProfile->Add(multProfile.fc24nw2D);

  multProfile.fPtVariancea2D = new TProfile2D(Form("fdPt2a2D%s", label.Data()), "Pt variance", nn, xbins, nn, xbins_uncorr);
  multProfile.fPtVariancea2D->Sumw2();
  listOfProfile->Add(multProfile.fPtVariancea2D);

  multProfile.fPtVarianceb2D = new TProfile2D(Form("fdPt2b2D%s", label.Data()), "Pt variance", nn, xbins, nn, xbins_uncorr);
  multProfile.fPtVarianceb2D->Sumw2();
  listOfProfile->Add(multProfile.fPtVarianceb2D);

  multProfile.fPtVariancec2D = new TProfile2D(Form("fdPt2c2D%s", label.Data()), "Pt variance", nn, xbins, nn, xbins_uncorr);
  multProfile.fPtVariancec2D->Sumw2();
  listOfProfile->Add(multProfile.fPtVariancec2D);

}

Bool_t AliAnalysisTaskVnPtCorr::AcceptAOD(AliAODEvent *inEv) {
  // LHC15i, LHC15l, LHC16, LHC17, LHC18: means: pp sample
  if (fPeriod.EqualTo("LHC15i") ||
      fPeriod.EqualTo("LHC15l") ||
      fPeriod.EqualTo("LHC16Preview") ||
      fPeriod.EqualTo("LHC17Preview") ||
      fPeriod.EqualTo("LHC18Preview") ||
      fPeriod.EqualTo("LHC16") ||
      fPeriod.EqualTo("LHC17") ||
      fPeriod.EqualTo("LHC18") ||
      fPeriod.EqualTo("LHC16_Closure") ||
      fPeriod.EqualTo("LHC17_Closure") ||
      fPeriod.EqualTo("LHC18_Closure") ||
      fPeriod.EqualTo("LHC16_simp") ||
      fPeriod.EqualTo("LHC17_simp") ||
      fPeriod.EqualTo("LHC18_simp") ||
      fPeriod.EqualTo("LHC16ZM") ||
      fPeriod.EqualTo("LHC17ZM") ||
      fPeriod.EqualTo("LHC18ZM") ) {
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);

    const auto pms(static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection")));
    const auto v0Est = pms->GetEstimator("V0M");
    if (v0Est->GetValue()/v0Est->GetMean() < fV0MRatioCut) return kFALSE;
  }

  if (fPeriod.EqualTo("LHC15o_pass2") || fPeriod.EqualTo("LHC15o_pass2_Closure")) {
    int currentRun = fAOD->GetRunNumber();
    if (currentRun == 245729 ||
        currentRun == 245731 ||
        currentRun == 245752 ||
        currentRun == 245759 ||
        currentRun == 245766 ||
        currentRun == 245775 ||
        currentRun == 245785 ||
        currentRun == 245793) {
      return kFALSE;
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
      fPeriod.EqualTo("LHC15o_pass2_Closure") ||
      fPeriod.EqualTo("LHC18qr_pass3_Closure") ||
      fPeriod.EqualTo("LHC15o_simp") ||
      fPeriod.EqualTo("LHC18q_simp") ||
      fPeriod.EqualTo("LHC18r_simp") ||
      fPeriod.EqualTo("LHC16qt") ||
      fPeriod.EqualTo("LHC16qt_Closure") ||
      fPeriod.EqualTo("LHC16qt_simp") ||
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

  bootstrap_value = (((int)(aodVtxZ * 233)) % 30 + 30) % 30;
  return kTRUE;
}

Bool_t AliAnalysisTaskVnPtCorr::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp) {
  // Pt cut
  if(mtr->Pt() < fMinPt) return kFALSE;
  if(mtr->Pt() > fMaxPt) return kFALSE;

  // DCA cut
  if(ltrackXYZ && vtxp) {
    mtr->GetXYZ(ltrackXYZ);
    ltrackXYZ[0] = ltrackXYZ[0]-vtxp[0];
    ltrackXYZ[1] = ltrackXYZ[1]-vtxp[1];
    ltrackXYZ[2] = abs(ltrackXYZ[2]-vtxp[2]);
  } else return kFALSE; //DCA cut is a must for now

  // Additional cut for TPCchi2perCluster
  if (mtr->GetTPCchi2perCluster()>fTPCchi2perCluster) return kFALSE;

  // Disable check DCAxy because we want to use the cut in FB96
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
    return fGFWSelection15o->AcceptTrack(mtr,ltrackXYZ,0, fCurrSystFlag ? kFALSE : kTRUE);
  } else {
    return fGFWSelection->AcceptTrack(mtr,ltrackXYZ,0, fCurrSystFlag ? kFALSE : kTRUE);
  }
}

Bool_t AliAnalysisTaskVnPtCorr::AcceptMCTruthTrack(AliAODMCParticle *mtrk) {
  // Pt cut
  if(mtrk->Pt() < fMinPt) return kFALSE;
  if(mtrk->Pt() > fMaxPt) return kFALSE;

  if(TMath::Abs(mtrk->Eta()) > fEtaCut) return kFALSE;

  if (!(mtrk->IsPhysicalPrimary())) return kFALSE;
  if (mtrk->Charge() == 0) return kFALSE;
  return kTRUE;
}


void AliAnalysisTaskVnPtCorr::CalculateProfile(PhysicsProfileVnPt& profile, double Ntrks, double Ntrks_uncorr) {
  //..calculate the PCC
  double Dn2GapLR = correlator.TwoGap8(0, 0).Re();
  double Dn4GapLR = correlator.FourGap8(0, 0, 0, 0).Re();
  if(NtrksAfter3subL > 1 && NtrksAfter3subM > 1 && NtrksAfter3subR > 1
      && Dn2GapLR != 0 && Dn4GapLR != 0)
  {
    eventWeight = sumWeight;
    eventWeight2 = sumWeight*sumWeight - sumWeight2;
    TComplex v22_3subLR = correlator.TwoGap8(2, -2);
    double v22Re_3subLR = v22_3subLR.Re()/Dn2GapLR;

    double meanPt = sumPtw/sumWeight;
    profile.fc22w->Fill(Ntrks, v22Re_3subLR, Dn2GapLR*eventWeight);
    profile.fPcc->Fill(Ntrks,  v22Re_3subLR*meanPt, Dn2GapLR*eventWeight);
    profile.fMeanPt->Fill(Ntrks, meanPt, eventWeight);

    profile.fc22w2D->Fill(Ntrks, Ntrks_uncorr, v22Re_3subLR, Dn2GapLR*eventWeight);
    profile.fPcc2D->Fill(Ntrks, Ntrks_uncorr, v22Re_3subLR*meanPt, Dn2GapLR*eventWeight);
    profile.fMeanPt2D->Fill(Ntrks, Ntrks_uncorr, meanPt, eventWeight);

    // Variance of Pt
    profile.fPtVariancea->Fill(Ntrks, (sumPtw*sumPtw-sumPt2w2)   / eventWeight2, eventWeight2);
    profile.fPtVarianceb->Fill(Ntrks, (sumWeight*sumPtw-sumPtw2) / eventWeight2, eventWeight2);
    profile.fPtVariancec->Fill(Ntrks, 1, eventWeight2);

    profile.fPtVariancea2D->Fill(Ntrks, Ntrks_uncorr, (sumPtw*sumPtw-sumPt2w2)   / eventWeight2, eventWeight2);
    profile.fPtVarianceb2D->Fill(Ntrks, Ntrks_uncorr, (sumWeight*sumPtw-sumPtw2) / eventWeight2, eventWeight2);
    profile.fPtVariancec2D->Fill(Ntrks, Ntrks_uncorr, 1, eventWeight2);

    // Variance of vn^2
    profile.fc22nw->Fill(Ntrks, v22Re_3subLR, Dn2GapLR);
    profile.fc22nw2D->Fill(Ntrks, Ntrks_uncorr, v22Re_3subLR, Dn2GapLR);

    TComplex v24_3subLR = correlator.FourGap8(2, 2, -2, -2);
    double v24Re_3subLR = v24_3subLR.Re()/Dn4GapLR;
    profile.fc24nw->Fill(Ntrks, v24Re_3subLR, Dn4GapLR);
    profile.fc24nw2D->Fill(Ntrks, Ntrks_uncorr, v24Re_3subLR, Dn4GapLR);
  }
}

const char* AliAnalysisTaskVnPtCorr::ReturnPPperiodMC(const Int_t runNumber) const
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

const char* AliAnalysisTaskVnPtCorr::ReturnPPperiod(const Int_t runNumber) const
{
  Bool_t isHM = kFALSE;
  if(fAliTrigger == AliVEvent::kHighMultV0) isHM = kTRUE;

  if(runNumber >= 252235 && runNumber <= 264347){ // LHC16
    if(!isHM && runNumber >= 252235 && runNumber <= 252375) return "LHC16de"; //d
    if(!isHM && runNumber >= 253437 && runNumber <= 253591) return "LHC16de"; //e
    if(runNumber >= 254128 && runNumber <= 254332) return "LHC16ghi"; //g
    if(runNumber >= 254604 && runNumber <= 255467) return "LHC16ghi"; //h
    if(runNumber >= 255539 && runNumber <= 255618) return "LHC16ghi"; //i
    if(runNumber >= 256219 && runNumber <= 256418) return "LHC16j";
    if(runNumber >= 256941 && runNumber <= 258537) return "LHC16k";
    if(runNumber >= 258962 && runNumber <= 259888) return "LHC16l";
    if(runNumber >= 262424 && runNumber <= 264035) return "LHC16o";
    if(runNumber >= 264076 && runNumber <= 264347) return "LHC16p";
  }

  if(runNumber >= 270581 && runNumber <= 282704){ // LHC17
    if(!isHM && runNumber >= 270581 && runNumber <= 270667) return "LHC17ce";
    if(runNumber >= 270822 && runNumber <= 270830){
      if(isHM) return "averaged";
      else return "LHC17ce";
    }
    if(runNumber >= 270854 && runNumber <= 270865){
      if(isHM) return "averaged";
      else return "LHC17f";
    }
    if(runNumber >= 271870 && runNumber <= 273103) return "LHC17h";
    if(runNumber >= 273591 && runNumber <= 274442) return "LHC17i";
    if(!isHM && runNumber >= 274593 && runNumber <= 274671) return "LHC17j";
    if(runNumber >= 274690 && runNumber <= 276508) return "LHC17k";
    if(runNumber >= 276551 && runNumber <= 278216) return "LHC17l";
    if(runNumber >= 278914 && runNumber <= 280140) return "LHC17m";
    if(runNumber >= 280282 && runNumber <= 281961) return "LHC17o";
    if(runNumber >= 282528 && runNumber <= 282704) return "LHC17r";
  }

  if(runNumber >= 285009 && runNumber <= 294925){ // LHC18
    if(runNumber >= 285009 && runNumber <= 285396){
      if(isHM) return "LHC18bd";
      else return "LHC18b";
    }
    if(runNumber >= 285978 && runNumber <= 286350){
      if(isHM) return "LHC18bd";
      else return "LHC18d";
    }
    if(runNumber >= 286380 && runNumber <= 286937) return "LHC18e";
    if(runNumber >= 287000 && runNumber <= 287658) return "LHC18f";
    if(runNumber >= 288804 && runNumber <= 288806){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(runNumber == 288943){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(runNumber >= 289165 && runNumber <= 289201){
      if(isHM) return "LHC18hjk";
      else return "LHC18ghijk";
    }
    if(!isHM && runNumber >= 288619 && runNumber <= 288750) return "LHC18ghijk"; //g, no HM event, only MB
    if(!isHM && runNumber >= 288861 && runNumber <= 288909) return "LHC18ghijk"; //i, no HM event, only MB
    if(runNumber >= 289240 && runNumber <= 289971) return "LHC18l";
    if(runNumber >= 290323 && runNumber <= 292839){
      if(isHM) return "LHC18m";
      else return "LHC18mn";
    }
    if(!isHM && runNumber >= 293357 && runNumber <= 293359) return "LHC18mn"; //n, no HM event, only MB
    if(runNumber >= 293475 && runNumber <= 293898) return "LHC18o";
    if(runNumber >= 294009 && runNumber <= 294925) return "LHC18p";
  }

  AliWarning("Unknown period! Returning averaged weights");
  return "averaged";
}


//_____________________________________________________________________________
void AliAnalysisTaskVnPtCorr::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate()");
}

ClassImp(PhysicsProfileVnPt);
PhysicsProfileVnPt::PhysicsProfileVnPt() :
  fMeanPt(nullptr),
  fc22w(nullptr),
  fPcc(nullptr),
  fc22nw(nullptr),
  fc24nw(nullptr),
  fPtVariancea(nullptr),
  fPtVarianceb(nullptr),
  fPtVariancec(nullptr),
  fMeanPt2D(nullptr),
  fc22w2D(nullptr),
  fPcc2D(nullptr),
  fc22nw2D(nullptr),
  fc24nw2D(nullptr),
  fPtVariancea2D(nullptr),
  fPtVarianceb2D(nullptr),
  fPtVariancec2D(nullptr) {}

  PhysicsProfileVnPt::PhysicsProfileVnPt(const PhysicsProfileVnPt& profile) :
    fMeanPt(nullptr),
    fc22w(nullptr),
    fPcc(nullptr),
    fc22nw(nullptr),
    fc24nw(nullptr),
    fPtVariancea(nullptr),
    fPtVarianceb(nullptr),
    fPtVariancec(nullptr),
    fMeanPt2D(nullptr),
    fc22w2D(nullptr),
    fPcc2D(nullptr),
    fc22nw2D(nullptr),
    fc24nw2D(nullptr),
    fPtVariancea2D(nullptr),
    fPtVarianceb2D(nullptr),
    fPtVariancec2D(nullptr) {}

    int AliAnalysisTaskVnPtCorr::GetEtaPtFlag(double dEta) {
      if(dEta > 0.0){
        if(dEta > 0.6) return 6;
        if(dEta > 0.4) return 4;
        if(dEta > 0.2) return 2;
        return 0;
      }
      else{
        if(dEta < -0.6) return 7;
        if(dEta < -0.4) return 5;
        if(dEta < -0.2) return 3;
        return 1;
      }
    }

AliMCEvent *AliAnalysisTaskVnPtCorr::getMCEvent() {
  AliMCEvent* ev = dynamic_cast<AliMCEvent*>(MCEvent());
  if(!ev) { AliFatal("MC event not found!"); return 0; }
  AliGenEventHeader *header = dynamic_cast<AliGenEventHeader*>(ev->GenEventHeader());
  if(!header) { AliFatal("MC event not generated!"); return 0; }
  AliCollisionGeometry* headerH;
  TString genName;
  TList *ltgen = (TList*)ev->GetCocktailList();
  if (ltgen) {
    for(auto&& listObject: *ltgen){
      genName = Form("%s",listObject->GetName());
      if (genName.Contains("Hijing")) {
        headerH = dynamic_cast<AliCollisionGeometry*>(listObject);
        break;
      }
    }
  }
  else headerH = dynamic_cast<AliCollisionGeometry*>(ev->GenEventHeader());
  if(headerH){
    fImpactParameterMC = headerH->ImpactParameter();
  }
  return ev;
}
