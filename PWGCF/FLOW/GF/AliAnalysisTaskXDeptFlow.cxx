/* -------------------------------------------
 * Maintainer: Mingrui Zhao
 Some naive try
 */
#include "AliAnalysisTaskNonlinearFlow.h"
#include "AliAnalysisTaskXDeptFlow.h"
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
ClassImp(AliAnalysisTaskXDeptFlow)
//___________________________________________________________________________
AliAnalysisTaskXDeptFlow::AliAnalysisTaskXDeptFlow():
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
    fUseCorrectedNTracks(false),
	  fCentralityCut(100),
    fUseNarrowBin(false),
    fExtremeEfficiency(0),
    fTPCchi2perCluster(4.0),
    fUseAdditionalDCACut(false),
    fUseDefaultWeight(false),
    fEtaGap3Sub(0.4),

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
  for (int i = 0; i < 10; i++) QDis[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap0P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap0M[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap10P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap10M[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap14P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap14M[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subL[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subM[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subR[i] = NULL;
}
//______________________________________________________________________________
AliAnalysisTaskXDeptFlow::AliAnalysisTaskXDeptFlow(const char *name, int _fNUA, int _fNUE, TString _fPeriod):
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
  fUseCorrectedNTracks(false),
	fCentralityCut(100),
  fUseNarrowBin(false),
  fExtremeEfficiency(0),
  fTPCchi2perCluster(4.0),
  fUseAdditionalDCACut(false),
  fUseDefaultWeight(false),
  fEtaGap3Sub(0.4),

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
  for (int i = 0; i < 10; i++) QDis[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap0P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap0M[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap10P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap10M[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap14P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap14M[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subL[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subM[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subR[i] = NULL;

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
AliAnalysisTaskXDeptFlow::AliAnalysisTaskXDeptFlow(const char *name):
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
  fUseCorrectedNTracks(false),
	fCentralityCut(100),
  fUseNarrowBin(false),
  fExtremeEfficiency(0),
  fTPCchi2perCluster(4.0),
  fUseAdditionalDCACut(false),
  fUseDefaultWeight(false),
  fEtaGap3Sub(0.4),

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
  for (int i = 0; i < 10; i++) QDis[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap0P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap0M[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap10P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap10M[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap14P[i] = NULL;
  for (int i = 0; i < 10; i++) QDisGap14M[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subL[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subM[i] = NULL;
  for (int i = 0; i < 10; i++) QDis3subR[i] = NULL;

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
AliAnalysisTaskXDeptFlow::~AliAnalysisTaskXDeptFlow()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects) delete fListOfObjects;
  if (fListOfProfile) delete fListOfProfile;
  for (int i = 0; i < 30; i++) {
    if (fListOfProfiles[i]) delete fListOfProfiles[i];
  }
  for (int i = 0; i < 10; i++) if (QDis[i]) delete QDis[i];
  for (int i = 0; i < 10; i++) if (QDisGap0P[i]) delete QDisGap0P[i];
  for (int i = 0; i < 10; i++) if (QDisGap0M[i]) delete QDisGap0M[i];
  for (int i = 0; i < 10; i++) if (QDisGap10P[i]) delete QDisGap10P[i];
  for (int i = 0; i < 10; i++) if (QDisGap10M[i]) delete QDisGap10M[i];
  for (int i = 0; i < 10; i++) if (QDisGap14P[i]) delete QDisGap14P[i];
  for (int i = 0; i < 10; i++) if (QDisGap14M[i]) delete QDisGap14M[i];
  for (int i = 0; i < 10; i++) if (QDis3subL[i]) delete QDis3subL[i];
  for (int i = 0; i < 10; i++) if (QDis3subM[i]) delete QDis3subM[i];
  for (int i = 0; i < 10; i++) if (QDis3subR[i]) delete QDis3subR[i];

  if (fGFWSelection) delete fGFWSelection;
  if (fGFWSelection15o) delete fGFWSelection15o;
}

//______________________________________________________________________________
void AliAnalysisTaskXDeptFlow::UserCreateOutputObjects()
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

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1 and LHC17n
    fGFWSelection15o = new AliGFWNFCuts();
    fGFWSelection15o->PrintSetup();
  } else {
    fGFWSelection = new AliGFWCuts();
    fGFWSelection->PrintSetup();
  }

  if (fNtrksName == "Mult") {
    nn = 2000;
    for (int i = 0; i <= 2000; i++) {
      xbins[i] = i/2000.0;
    }
  }

  // Calculate the used observables:
  // unsigned fgFlowHarmonics;        //! calculate v2, v3, v4, v5
  // unsigned fgFlowHarmonicsHigher;  //! calculate v6, v7, v8 ..
  // unsigned fgFlowHarmonicsMult;    //! calculate v2{4}  v2{6}, v2{8}
  // unsigned fgNonlinearFlow;        //! calculate v_4,22, v_5,32
  // unsigned fgSymmetricCumulants;   //! calculate SC(3,2), SC(4,2)
  //

  fgTwoParticleCorrelation       = fgFlowHarmonics | fgSymmetricCumulants | fgNonlinearFlow;
  fgTwoParticleCorrelationHigher = fgFlowHarmonicsHigher;
  fgThreeParticleCorrelation     = fgNonlinearFlow;
  fgFourParticleCorrelation      = fgNonlinearFlow | fgSymmetricCumulants | fgFlowHarmonicsMult;
  fgSixParticleCorrelation       = fgFlowHarmonicsMult;
  fgEightParticleCorrelation     = fgFlowHarmonicsMult;

  fuTwoParticleCorrelationStandard = fgTwoParticleCorrelation & knStandard;
  fuTwoParticleCorrelation0Gap     = fgTwoParticleCorrelation & kn0Gap;
  fuTwoParticleCorrelationLargeGap = fgTwoParticleCorrelation & knLargeGap;
  fuTwoParticleCorrelationThreeSub = fgTwoParticleCorrelation & knThreeSub;
  fuTwoParticleCorrelationGapScan  = fgTwoParticleCorrelation & knGapScan;

  fuTwoParticleCorrelationHigherStandard = fgTwoParticleCorrelationHigher & knStandard;
  fuTwoParticleCorrelationHigher0Gap     = fgTwoParticleCorrelationHigher & kn0Gap;
  fuTwoParticleCorrelationHigherLargeGap = fgTwoParticleCorrelationHigher & knLargeGap;
  fuTwoParticleCorrelationHigherThreeSub = fgTwoParticleCorrelationHigher & knThreeSub;
  fuTwoParticleCorrelationHigherGapScan  = fgTwoParticleCorrelationHigher & knGapScan;

  fuThreeParticleCorrelationStandard = fgThreeParticleCorrelation & knStandard;
  fuThreeParticleCorrelation0Gap     = fgThreeParticleCorrelation & kn0Gap;
  fuThreeParticleCorrelationLargeGap = fgThreeParticleCorrelation & knLargeGap;
  fuThreeParticleCorrelationThreeSub = fgThreeParticleCorrelation & knThreeSub;
  fuThreeParticleCorrelationGapScan  = fgThreeParticleCorrelation & knGapScan;

  fuFourParticleCorrelationStandard = fgFourParticleCorrelation & knStandard;
  fuFourParticleCorrelation0Gap     = fgFourParticleCorrelation & kn0Gap;
  fuFourParticleCorrelationLargeGap = fgFourParticleCorrelation & knLargeGap;
  fuFourParticleCorrelationThreeSub = fgFourParticleCorrelation & knThreeSub;
  fuFourParticleCorrelationGapScan  = fgFourParticleCorrelation & knGapScan;

  fuSixParticleCorrelationStandard  = fgSixParticleCorrelation & knStandard;
  fuSixParticleCorrelation0Gap      = fgSixParticleCorrelation & kn0Gap;

  fuEightParticleCorrelationStandard = fgEightParticleCorrelation & knStandard;
  fuEightParticleCorrelation0Gap     = fgEightParticleCorrelation & kn0Gap;

  if (fuTwoParticleCorrelationGapScan) {
    fuTwoParticleCorrelationStandard = true;
    fuTwoParticleCorrelation0Gap     = true;
    fuTwoParticleCorrelationLargeGap = true;
  }
  if (fuTwoParticleCorrelationHigherGapScan) {
    fuTwoParticleCorrelationHigherStandard = true;
    fuTwoParticleCorrelationHigher0Gap     = true;
    fuTwoParticleCorrelationHigherLargeGap = true;
  }

  if (fuThreeParticleCorrelationGapScan) {
    fuThreeParticleCorrelationStandard = true;
    fuThreeParticleCorrelation0Gap     = true;
    fuThreeParticleCorrelationLargeGap = true;
  }
  if (fuFourParticleCorrelationGapScan) {
    fuFourParticleCorrelationStandard = true;
    fuFourParticleCorrelation0Gap     = true;
    fuFourParticleCorrelationLargeGap = true;
  }

  if (fuTwoParticleCorrelationStandard || fuTwoParticleCorrelationHigherStandard
      || fuThreeParticleCorrelationStandard || fuFourParticleCorrelationStandard) fuQStandard = true;
  if (fuTwoParticleCorrelation0Gap || fuTwoParticleCorrelationHigher0Gap
      || fuThreeParticleCorrelation0Gap || fuFourParticleCorrelation0Gap) fuQ0Gap = true;
  if (fuTwoParticleCorrelationLargeGap || fuTwoParticleCorrelationHigherLargeGap
      || fuThreeParticleCorrelationLargeGap || fuFourParticleCorrelationLargeGap) fuQLargeGap = true;
  if (fuTwoParticleCorrelationThreeSub || fuTwoParticleCorrelationHigherThreeSub
      || fuThreeParticleCorrelationThreeSub || fuFourParticleCorrelationThreeSub) fuQThreeSub = true;
  if (fuTwoParticleCorrelationGapScan || fuTwoParticleCorrelationHigherGapScan
      || fuThreeParticleCorrelationGapScan || fuFourParticleCorrelationGapScan) fuQGapScan = true;

  hEventCount = new TH1D("hEventCount", "; centrality;;", 5, 0.5, 5.5);
  fListOfObjects->Add(hEventCount);

  hMult = new TH1F("hMult", ";number of tracks; entries", nn, xbins);
  hMult->Sumw2();
  fListOfObjects->Add(hMult);

  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
  fVtxAfterCuts->Sumw2();
  fListOfObjects->Add(fVtxAfterCuts);

  fCentralityDis = new TH1F("fCentralityDis", "centrality distribution; centrality; Counts", 100, 0, 100);
  fListOfObjects->Add(fCentralityDis);

  fV0CentralityDis = new TH1F("fV0CentralityDis", "centrality V0/<V0> distribution; centrality; Counts", 100, 0, 100);
  fListOfObjects->Add(fV0CentralityDis);

  fV0CentralityDisNarrow = new TH1F("fV0CentralityDisNarrow", "centrality V0/<V0> distribution; centrality; Counts", 1000, 0, 10);
  fListOfObjects->Add(fV0CentralityDisNarrow);

  fPhiDis1DBefore = new TH1D("hPhiDisBefore", "phi distribution before the weight correction", 60, 0, 2*3.1415926);
  fListOfObjects->Add(fPhiDis1DBefore);
  fPhiDis1D  = new TH1D("hPhiDis", "phi distribution after the weight correction", 60, 0, 2*3.1415926);
  fListOfObjects->Add(fPhiDis1D);
  fEtaDis = new TH1D("hEtaDis", "eta distribution", 100, -2, 2);
  fListOfObjects->Add(fEtaDis);
  fPtDis = new TH1D("hPtDis", "pt distribution", 100, 0, 5);
  fListOfObjects->Add(fPtDis);
  hTracksCorrection2d = new TH2D("hTracksCorrection2d", "Correlation table for number of tracks table", nn, xbins, nn, xbins);
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
  for (int h = 0; h < 6; h++) {
    QDis[h] = new TH2D(Form("Q%dDis", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDisGap0P[h] = new TH2D(Form("Q%dDisGap0P", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDisGap0M[h] = new TH2D(Form("Q%dDisGap0M", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDisGap10P[h] = new TH2D(Form("Q%dDisGap10P", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDisGap10M[h] = new TH2D(Form("Q%dDisGap10M", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDisGap14P[h] = new TH2D(Form("Q%dDisGap14P", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDisGap14M[h] = new TH2D(Form("Q%dDisGap14M", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDis3subL[h] = new TH2D(Form("Q%dDis3subL", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDis3subM[h] = new TH2D(Form("Q%dDis3subM", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    QDis3subR[h] = new TH2D(Form("Q%dDis3subR", h+2), "Q distribution", 100, -1, 1, 100, -1, 1);
    fListOfObjects->Add(QDis[h]);
    fListOfObjects->Add(QDisGap0P[h]);
    fListOfObjects->Add(QDisGap0M[h]);
    fListOfObjects->Add(QDisGap10P[h]);
    fListOfObjects->Add(QDisGap10M[h]);
    fListOfObjects->Add(QDisGap14P[h]);
    fListOfObjects->Add(QDisGap14M[h]);
    fListOfObjects->Add(QDis3subL[h]);
    fListOfObjects->Add(QDis3subM[h]);
    fListOfObjects->Add(QDis3subR[h]);
  }

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse) { AliError("AliPIDResponse not found!"); return; }

  fPIDCombined = new AliPIDCombined();
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask


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
void AliAnalysisTaskXDeptFlow::NotifyRun() {
    if (fAddTPCPileupCuts) {
      Bool_t dummy = fEventCuts.AcceptEvent(InputEvent());
	  fEventCuts.fUseVariablesCorrelationCuts = true;
      fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
      fEventCuts.fESDvsTPConlyLinearCut[0] = fESDvsTPConlyLinearCut;
    }
}

//______________________________________________________________________________
void AliAnalysisTaskXDeptFlow::UserExec(Option_t *)
{
  // bootstrap_value = rand.Integer(30);

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

  // Setup AliGFWCuts for a specific systematics
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

  //..standard event plots (cent. percentiles, mult-vs-percentile)
  const auto pms(static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection")));
  const auto dCentrality(pms->GetMultiplicityPercentile("V0M"));
  float centrV0 = dCentrality;
  float cent = dCentrality;
  float centSPD = 0;

  fCentralityDis->Fill(cent);
  fV0CentralityDis->Fill(centrV0);
  fV0CentralityDisNarrow->Fill(centrV0);

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

  // Post output data.
  PostData(1, fListOfObjects);
  int outputslot = 2;
  PostData(2, fListOfProfile);
  for (int i = 0; i < 30; i++) {
    outputslot++;
    PostData(outputslot, fListOfProfiles[i]);
  }

}

void AliAnalysisTaskXDeptFlow::NTracksCalculation(AliVEvent* aod) {
  const int nAODTracks = aod->GetNumberOfTracks();
  NtrksCounter = 0;
  NTracksCorrected = 0.01;
  NTracksUncorrected = 0.01;
  NKaon = 0;

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
    if (IdentifyTrack(aodTrk) == 2) NKaon += 1;
    NTracksCorrected += weightPt;
  } // end loop of all track
  if (!fUseCorrectedNTracks) {
    NtrksCounter = NKaon / NTracksUncorrected;
  } else {
    NtrksCounter = NKaon / NTracksCorrected; 
  }
  hTracksCorrection2d->Fill(NTracksUncorrected, NTracksCorrected);
  hnCorrectedTracks->Fill(NtrksCounter, NTracksCorrected);
}

//________________________________________________________________________
void AliAnalysisTaskXDeptFlow::AnalyzeAOD(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus)
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


  //..for DCA
  double pos[3], vz, vx, vy;
  vz = aod->GetPrimaryVertex()->GetZ();
  vx = aod->GetPrimaryVertex()->GetX();
  vy = aod->GetPrimaryVertex()->GetY();
  double vtxp[3] = {vx, vy, vz};

  double Qcos[20][20] = {0};
  double Qsin[20][20] = {0};
  double QcosGap0M[20][20] = {0};
  double QsinGap0M[20][20] = {0};
  double QcosGap0P[20][20] = {0};
  double QsinGap0P[20][20] = {0};
  double QcosGap2M[20][20] = {0};
  double QsinGap2M[20][20] = {0};
  double QcosGap2P[20][20] = {0};
  double QsinGap2P[20][20] = {0};
  double QcosGap4M[20][20] = {0};
  double QsinGap4M[20][20] = {0};
  double QcosGap4P[20][20] = {0};
  double QsinGap4P[20][20] = {0};
  double QcosGap6M[20][20] = {0};
  double QsinGap6M[20][20] = {0};
  double QcosGap6P[20][20] = {0};
  double QsinGap6P[20][20] = {0};
  double QcosGap8M[20][20] = {0};
  double QsinGap8M[20][20] = {0};
  double QcosGap8P[20][20] = {0};
  double QsinGap8P[20][20] = {0};
  double QcosGap10M[20][20] = {0};
  double QsinGap10M[20][20] = {0};
  double QcosGap10P[20][20] = {0};
  double QsinGap10P[20][20] = {0};
  double QcosGap14M[20][20] = {0};
  double QsinGap14M[20][20] = {0};
  double QcosGap14P[20][20] = {0};
  double QsinGap14P[20][20] = {0};
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

    //..calculate Q-vectors
    //..no eta gap
    // Calculate the values upto v7
    if (fuQStandard) {
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          Qcos[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
          Qsin[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
        }
      }
    }
    //..Gap > 0.0
    if (fuQ0Gap) {
      if(aodTrk->Eta() < 0) {
        NtrksAfterGap0M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
      if(aodTrk->Eta() > 0) {
        NtrksAfterGap0P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
    }

    if (fuQGapScan) {
      //..Gap > 0.2
      if(aodTrk->Eta() < -0.1) {
        NtrksAfterGap2M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
      if(aodTrk->Eta() > 0.1) {
        NtrksAfterGap2P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }

      //..Gap > 0.4
      if(aodTrk->Eta() < -0.2) {
        NtrksAfterGap4M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
      if(aodTrk->Eta() > 0.2) {
        NtrksAfterGap4P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }

      //..Gap > 0.6
      if(aodTrk->Eta() < -0.3) {
        NtrksAfterGap6M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap6M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap6M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
      if(aodTrk->Eta() > 0.3) {
        NtrksAfterGap6P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap6P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap6P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }

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
    }

    if (fuQLargeGap) {
      //..Gap > 1.0
      if(aodTrk->Eta() < -0.5) {
        NtrksAfterGap10M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap10M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap10M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
      if(aodTrk->Eta() > 0.5) {
        NtrksAfterGap10P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap10P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap10P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }

      //..Gap > 1.4
      if(aodTrk->Eta() < -0.7) {
        NtrksAfterGap14M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
      if(aodTrk->Eta() > 0.7) {
        NtrksAfterGap14P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
    }

    if (fuQThreeSub) {
      //..3-subevent method
      if(aodTrk->Eta() < -fEtaGap3Sub) {//..left part
        NtrksAfter3subL += 1;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
      if(aodTrk->Eta() >= -fEtaGap3Sub && aodTrk->Eta() <= fEtaGap3Sub) {//..middle part
        NtrksAfter3subM += 1;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
      if(aodTrk->Eta() > fEtaGap3Sub) {//..right part
        NtrksAfter3subR += 1;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*aodTrk->Phi());
            QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*aodTrk->Phi());
          }
        }
      }
    }
  } // end loop of all track

  //............................
  //..GENERIC FRAMEWORK RP
  //............................

  //..calculate Q-vector for each harmonics n and power p
  if (fuQStandard) correlator.FillQVector(correlator.Qvector, Qcos, Qsin);
  if (fuQ0Gap) {
    correlator.FillQVector(correlator.Qvector0M, QcosGap0M, QsinGap0M);
    correlator.FillQVector(correlator.Qvector0P, QcosGap0P, QsinGap0P);
  }
  if (fuQGapScan) {
    correlator.FillQVector(correlator.Qvector2M, QcosGap2M, QsinGap2M);
    correlator.FillQVector(correlator.Qvector2P, QcosGap2P, QsinGap2P);
    correlator.FillQVector(correlator.Qvector4M, QcosGap4M, QsinGap4M);
    correlator.FillQVector(correlator.Qvector4P, QcosGap4P, QsinGap4P);
    correlator.FillQVector(correlator.Qvector6M, QcosGap6M, QsinGap6M);
    correlator.FillQVector(correlator.Qvector6P, QcosGap6P, QsinGap6P);
    correlator.FillQVector(correlator.Qvector8M, QcosGap8M, QsinGap8M);
    correlator.FillQVector(correlator.Qvector8P, QcosGap8P, QsinGap8P);
  }
  if (fuQLargeGap) {
    correlator.FillQVector(correlator.Qvector10M, QcosGap10M, QsinGap10M);
    correlator.FillQVector(correlator.Qvector10P, QcosGap10P, QsinGap10P);
    correlator.FillQVector(correlator.Qvector14M, QcosGap14M, QsinGap14M);
    correlator.FillQVector(correlator.Qvector14P, QcosGap14P, QsinGap14P);
  }
  if (fuQThreeSub) {
    correlator.FillQVector(correlator.QvectorSubLeft, QcosSubLeft, QsinSubLeft);
    correlator.FillQVector(correlator.QvectorSubRight, QcosSubRight, QsinSubRight);
    correlator.FillQVector(correlator.QvectorSubMiddle, QcosSubMiddle, QsinSubMiddle);
  }

  for (int h = 0; h < 6; h++) {
    QDis[h]->Fill(correlator.Q(h+2,1).Re(), correlator.Q(h+2,1).Im());
    QDisGap0P[h]->Fill(correlator.QGap0P(h+2,1).Re(), correlator.QGap0P(h+2,1).Im());
    QDisGap0M[h]->Fill(correlator.QGap0M(h+2,1).Re(), correlator.QGap0M(h+2,1).Im());
    QDisGap10P[h]->Fill(correlator.QGap10P(h+2,1).Re(), correlator.QGap10P(h+2,1).Im());
    QDisGap10M[h]->Fill(correlator.QGap10M(h+2,1).Re(), correlator.QGap10M(h+2,1).Im());
    QDisGap14P[h]->Fill(correlator.QGap14P(h+2,1).Re(), correlator.QGap14P(h+2,1).Im());
    QDisGap14M[h]->Fill(correlator.QGap14M(h+2,1).Re(), correlator.QGap14M(h+2,1).Im());
    QDis3subL[h]->Fill(correlator.QsubLeft(h+2,1).Re(), correlator.QsubLeft(h+2,1).Im());
    QDis3subM[h]->Fill(correlator.QsubMiddle(h+2,1).Re(), correlator.QsubMiddle(h+2,1).Im());
    QDis3subR[h]->Fill(correlator.QsubRight(h+2,1).Re(), correlator.QsubRight(h+2,1).Im());
  }


  if (fNtrksName == "Mult") {
    CalculateProfile(multProfile, NtrksCounter);
    CalculateProfile(multProfile_bin[bootstrap_value], NtrksCounter);
  } else {
    CalculateProfile(multProfile, cent);
    CalculateProfile(multProfile_bin[bootstrap_value], cent);
  }

}

//________________________________________________________________________
void AliAnalysisTaskXDeptFlow::AnalyzeMCTruth(AliVEvent* aod, float centrV0, float cent, float centSPD, float fVtxZ, bool fPlus)
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


  //..for DCA
  // double pos[3], vz, vx, vy;
  // vz = aod->GetPrimaryVertex()->GetZ();
  // vx = aod->GetPrimaryVertex()->GetX();
  // vy = aod->GetPrimaryVertex()->GetY();
  // double vtxp[3] = {vx, vy, vz};
  // Assume that DCA cuts not needed here

  double Qcos[20][20] = {0};
  double Qsin[20][20] = {0};
  double QcosGap0M[20][20] = {0};
  double QsinGap0M[20][20] = {0};
  double QcosGap0P[20][20] = {0};
  double QsinGap0P[20][20] = {0};
  double QcosGap2M[20][20] = {0};
  double QsinGap2M[20][20] = {0};
  double QcosGap2P[20][20] = {0};
  double QsinGap2P[20][20] = {0};
  double QcosGap4M[20][20] = {0};
  double QsinGap4M[20][20] = {0};
  double QcosGap4P[20][20] = {0};
  double QsinGap4P[20][20] = {0};
  double QcosGap6M[20][20] = {0};
  double QsinGap6M[20][20] = {0};
  double QcosGap6P[20][20] = {0};
  double QsinGap6P[20][20] = {0};
  double QcosGap8M[20][20] = {0};
  double QsinGap8M[20][20] = {0};
  double QcosGap8P[20][20] = {0};
  double QsinGap8P[20][20] = {0};
  double QcosGap10M[20][20] = {0};
  double QsinGap10M[20][20] = {0};
  double QcosGap10P[20][20] = {0};
  double QsinGap10P[20][20] = {0};
  double QcosGap14M[20][20] = {0};
  double QsinGap14M[20][20] = {0};
  double QcosGap14P[20][20] = {0};
  double QsinGap14P[20][20] = {0};
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

    // track->GetXYZ(pos);
    if (!AcceptMCTruthTrack(track)) continue;

    NtrksAfter += 1;

    //..get phi-weight for NUA correction
    double weight = 1;
    double weightPt = 1;

    //..calculate Q-vectors
    //..no eta gap
    // Calculate the values upto v7
    if (fuQStandard) {
      for(int iharm=0; iharm<8; iharm++) {
        for(int ipow=0; ipow<6; ipow++) {
          Qcos[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
          Qsin[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
        }
      }
    }
    //..Gap > 0.0
    if (fuQ0Gap) {
      if(track->Eta() < 0) {
        NtrksAfterGap0M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
      if(track->Eta() > 0) {
        NtrksAfterGap0P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
    }

    if (fuQGapScan) {
      //..Gap > 0.2
      if(track->Eta() < -0.1) {
        NtrksAfterGap2M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
      if(track->Eta() > 0.1) {
        NtrksAfterGap2P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }

      //..Gap > 0.4
      if(track->Eta() < -0.2) {
        NtrksAfterGap4M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
      if(track->Eta() > 0.2) {
        NtrksAfterGap4P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }

      //..Gap > 0.6
      if(track->Eta() < -0.3) {
        NtrksAfterGap6M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap6M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap6M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
      if(track->Eta() > 0.3) {
        NtrksAfterGap6P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap6P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap6P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }

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
    }

    if (fuQLargeGap) {
      //..Gap > 1.0
      if(track->Eta() < -0.5) {
        NtrksAfterGap10M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap10M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap10M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
      if(track->Eta() > 0.5) {
        NtrksAfterGap10P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap10P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap10P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }

      //..Gap > 1.4
      if(track->Eta() < -0.7) {
        NtrksAfterGap14M++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
      if(track->Eta() > 0.7) {
        NtrksAfterGap14P++;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
    }

    if (fuQThreeSub) {
      //..3-subevent method
      if(track->Eta() < -fEtaGap3Sub) {//..left part
        NtrksAfter3subL += 1;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
      if(track->Eta() >= -fEtaGap3Sub && track->Eta() <= fEtaGap3Sub) {//..middle part
        NtrksAfter3subM += 1;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
      if(track->Eta() > fEtaGap3Sub) {//..right part
        NtrksAfter3subR += 1;
        for(int iharm=0; iharm<8; iharm++) {
          for(int ipow=0; ipow<6; ipow++) {
            QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*track->Phi());
            QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*track->Phi());
          }
        }
      }
    }
  } // end loop of all track

  //............................
  //..GENERIC FRAMEWORK RP
  //............................

  //..calculate Q-vector for each harmonics n and power p
  if (fuQStandard) correlator.FillQVector(correlator.Qvector, Qcos, Qsin);
  if (fuQ0Gap) {
    correlator.FillQVector(correlator.Qvector0M, QcosGap0M, QsinGap0M);
    correlator.FillQVector(correlator.Qvector0P, QcosGap0P, QsinGap0P);
  }
  if (fuQGapScan) {
    correlator.FillQVector(correlator.Qvector2M, QcosGap2M, QsinGap2M);
    correlator.FillQVector(correlator.Qvector2P, QcosGap2P, QsinGap2P);
    correlator.FillQVector(correlator.Qvector4M, QcosGap4M, QsinGap4M);
    correlator.FillQVector(correlator.Qvector4P, QcosGap4P, QsinGap4P);
    correlator.FillQVector(correlator.Qvector6M, QcosGap6M, QsinGap6M);
    correlator.FillQVector(correlator.Qvector6P, QcosGap6P, QsinGap6P);
    correlator.FillQVector(correlator.Qvector8M, QcosGap8M, QsinGap8M);
    correlator.FillQVector(correlator.Qvector8P, QcosGap8P, QsinGap8P);
  }
  if (fuQLargeGap) {
    correlator.FillQVector(correlator.Qvector10M, QcosGap10M, QsinGap10M);
    correlator.FillQVector(correlator.Qvector10P, QcosGap10P, QsinGap10P);
    correlator.FillQVector(correlator.Qvector14M, QcosGap14M, QsinGap14M);
    correlator.FillQVector(correlator.Qvector14P, QcosGap14P, QsinGap14P);
  }
  if (fuQThreeSub) {
    correlator.FillQVector(correlator.QvectorSubLeft, QcosSubLeft, QsinSubLeft);
    correlator.FillQVector(correlator.QvectorSubRight, QcosSubRight, QsinSubRight);
    correlator.FillQVector(correlator.QvectorSubMiddle, QcosSubMiddle, QsinSubMiddle);
  }

  if (fNtrksName == "Mult") {
    CalculateProfile(multProfile, NtrksCounter);
    CalculateProfile(multProfile_bin[bootstrap_value], NtrksCounter);
  } else {
    CalculateProfile(multProfile, cent);
    CalculateProfile(multProfile_bin[bootstrap_value], cent);
  }

}

//____________________________________________________________________
//	END OF MAIN PROGRAM
//____________________________________________________________________

double AliAnalysisTaskXDeptFlow::GetPtWeight(double pt, double eta, float vz, double runNumber)
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


Bool_t AliAnalysisTaskXDeptFlow::LoadWeightsSystematics() {

  // If the period is not pPb LHC16qt
  if (!fPeriod.EqualTo("LHC16qt")) {
    // Only if it is the new LHC16,17,18, We need the period NUA
    if (fPeriod.EqualTo("LHC16") || fPeriod.EqualTo("LHC17") || fPeriod.EqualTo("LHC18")) {
      std::string ppperiod = ReturnPPperiod(fAOD->GetRunNumber());
      if(fCurrSystFlag == 0) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%s", ppperiod.c_str()));
      if(!fWeightsSystematics)
      {
        printf("Weights could not be found in list!\n");
        return kFALSE;
      }
      fWeightsSystematics->CreateNUA();
    } else if (fPeriod.EqualTo("LHC16Preview") || fPeriod.EqualTo("LHC17Preview") || fPeriod.EqualTo("LHC18Preview")) {

      if(fCurrSystFlag == 0 || fUseDefaultWeight) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i",fAOD->GetRunNumber()));
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
      if(fCurrSystFlag == 0 || fUseDefaultWeight) fWeightsSystematics = (AliGFWWeights*)fFlowWeightsList->FindObject(Form("w%i",fAOD->GetRunNumber()));
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

Bool_t AliAnalysisTaskXDeptFlow::LoadPtWeights() {
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

double AliAnalysisTaskXDeptFlow::GetWeightKatarina(double phi, double eta, double vz) {
  double weight = hPhiWeightRun->GetBinContent(hPhiWeightRun->GetXaxis()->FindBin(phi),
      hPhiWeightRun->GetYaxis()->FindBin(eta),
      hPhiWeightRun->GetZaxis()->FindBin(vz));
  return weight;
}

// Load Katarina's weights
Bool_t AliAnalysisTaskXDeptFlow::LoadWeightsKatarina() {
  hPhiWeightRun = (TH3F*)fFlowWeightsList->FindObject(Form("fPhiWeight_%0.lf", (double)(fAOD->GetRunNumber())));
  if (!hPhiWeightRun) {
    printf("Weights could not be found in list!\n");
    return kFALSE;
  }
  return kTRUE;
}

double AliAnalysisTaskXDeptFlow::GetPtWeightKatarina(double pt, double eta, double vz)
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
Bool_t AliAnalysisTaskXDeptFlow::LoadPtWeightsKatarina() {
  hTrackEfficiencyRun = (TH3F*)fFlowPtWeightsList->FindObject(Form("eff_LHC15o_HIJING_%.0lf", (double)(fAOD->GetRunNumber())));
  if (!hTrackEfficiencyRun) {
    printf("Pt Weights could not be found in list!\n");
    return kFALSE;
  }
  return kTRUE;
}

Double_t AliAnalysisTaskXDeptFlow::GetFlowWeightSystematics(const AliVParticle* track, double fVtxZ, const PartSpecies species) {

  double dPhi = track->Phi();
  double dEta = track->Eta();
  double dVz = fVtxZ;
  double dWeight = 1.0;
  dWeight = fWeightsSystematics->GetNUA(dPhi, dEta, dVz);
  return dWeight;
}


void AliAnalysisTaskXDeptFlow::InitProfile(PhysicsProfile& multProfile, TString label, TList* listOfProfile) {

  // h = 0 -> 5 : v2, v3, v4, v5, v6, v7
  for(int h=0; h<6; h++)
  {
    if (fuTwoParticleCorrelationStandard) {
      multProfile.fChcn2[h] = new TProfile(Form("fChc%d{2}%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
      multProfile.fChcn2[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn2[h]);
    }

    if (fuTwoParticleCorrelation0Gap) {
      multProfile.fChcn2_Gap0[h] = new TProfile(Form("fChc%d{2}_Gap0%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
      multProfile.fChcn2_Gap0[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn2_Gap0[h]);
    }

    if (fuTwoParticleCorrelationLargeGap) {
      multProfile.fChcn2_Gap10[h] = new TProfile(Form("fChc%d{2}_Gap10%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
      multProfile.fChcn2_Gap10[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn2_Gap10[h]);
    }

    if (fuTwoParticleCorrelationLargeGap && h == 0) {
      multProfile.fChcn2_Gap14[h] = new TProfile(Form("fChc%d{2}_Gap14%s", h+2, label.Data()), "<<2>> Re; # of tracks", nn, xbins);
      multProfile.fChcn2_Gap14[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn2_Gap14[h]);
    }

    if (fuTwoParticleCorrelationThreeSub) {
      multProfile.fChcn2_3subLM[h] = new TProfile(Form("fChc%d{2}_3subLM%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, left+middle; # of tracks", nn, xbins);
      multProfile.fChcn2_3subLM[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn2_3subLM[h]);

      multProfile.fChcn2_3subRM[h] = new TProfile(Form("fChc%d{2}_3subRM%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, right+middle; # of tracks", nn, xbins);
      multProfile.fChcn2_3subRM[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn2_3subRM[h]);

      multProfile.fChcn2_3subLR[h] = new TProfile(Form("fChc%d{2}_3subLR%s", h+2, label.Data()), "<<2>> Re for 3-subevent method, left+right; # of tracks", nn, xbins);
      multProfile.fChcn2_3subLR[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn2_3subLR[h]);
    }

    if (fuFourParticleCorrelationStandard && h == 0) {
      multProfile.fChcn4[h] = new TProfile(Form("fChc%d{4}%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
      multProfile.fChcn4[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4[h]);
    }

    if (fuFourParticleCorrelation0Gap && h == 0) {
      multProfile.fChcn4_Gap0[h] = new TProfile(Form("fChc%d{4}_Gap0%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
      multProfile.fChcn4_Gap0[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_Gap0[h]);
    }

    if (fuFourParticleCorrelationLargeGap && h == 0) {
      multProfile.fChcn4_Gap10[h] = new TProfile(Form("fChc%d{4}_Gap10%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
      multProfile.fChcn4_Gap10[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_Gap10[h]);
    }

    if (fuFourParticleCorrelationGapScan && h == 0) {
      multProfile.fChcn4_Gap2[h] = new TProfile(Form("fChc%d{4}_Gap2%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
      multProfile.fChcn4_Gap2[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_Gap2[h]);

      multProfile.fChcn4_Gap4[h] = new TProfile(Form("fChc%d{4}_Gap4%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
      multProfile.fChcn4_Gap4[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_Gap4[h]);

      multProfile.fChcn4_Gap6[h] = new TProfile(Form("fChc%d{4}_Gap6%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
      multProfile.fChcn4_Gap6[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_Gap6[h]);

      multProfile.fChcn4_Gap8[h] = new TProfile(Form("fChc%d{4}_Gap8%s", h+2, label.Data()), "<<4>> Re; # of tracks", nn, xbins);
      multProfile.fChcn4_Gap8[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_Gap8[h]);
    }


    if (fuFourParticleCorrelationThreeSub && h == 0) {
      multProfile.fChcn4_3subLLMR[h] = new TProfile(Form("fChc%d{4}_3subLLMR%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
      multProfile.fChcn4_3subLLMR[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_3subLLMR[h]);

      multProfile.fChcn4_3subRRML[h] = new TProfile(Form("fChc%d{4}_3subRRML%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
      multProfile.fChcn4_3subRRML[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_3subRRML[h]);

      multProfile.fChcn4_3subMMLR[h] = new TProfile(Form("fChc%d{4}_3subMMLR%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
      multProfile.fChcn4_3subMMLR[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_3subMMLR[h]);

      multProfile.fChcn4_3subGap2[h] = new TProfile(Form("fChc%d{4}_3subGap2%s", h+2, label.Data()), "<<4>> 3-subevent method; # of tracks", nn, xbins);
      multProfile.fChcn4_3subGap2[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn4_3subGap2[h]);
    }


    if (fuSixParticleCorrelationStandard && h == 0) {
      multProfile.fChcn6[h] = new TProfile(Form("fChc%d{6}%s", h+2, label.Data()), "<<6>> Re; # of tracks", nn, xbins);
      multProfile.fChcn6[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn6[h]);
    }

    if (fuSixParticleCorrelation0Gap && h == 0) {
      multProfile.fChcn6_Gap0[h] = new TProfile(Form("fChc%d{6}_Gap0%s", h+2, label.Data()), "<<6>> Re; # of tracks", nn, xbins);
      multProfile.fChcn6_Gap0[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn6_Gap0[h]);
    }

    if (fuEightParticleCorrelationStandard && h == 0) {
      multProfile.fChcn8[h] = new TProfile(Form("fChc%d{8}%s", h+2, label.Data()), "<<8>> Re; # of tracks", nn, xbins);
      multProfile.fChcn8[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn8[h]);
    }

    if (fuEightParticleCorrelation0Gap && h == 0) {
      multProfile.fChcn8_Gap0[h] = new TProfile(Form("fChc%d{8}_Gap0%s", h+2, label.Data()), "<<8>> Re; # of tracks", nn, xbins);
      multProfile.fChcn8_Gap0[h]->Sumw2();
      listOfProfile->Add(multProfile.fChcn8_Gap0[h]);
    }

  } // harmonics

  if (fuThreeParticleCorrelationStandard) {
    multProfile.fChc422 = new TProfile(Form("fChc422%s", label.Data()), "", nn, xbins);
    multProfile.fChc422->Sumw2();
    listOfProfile->Add(multProfile.fChc422);

    multProfile.fChc532 = new TProfile(Form("fChc532%s", label.Data()), "", nn, xbins);
    multProfile.fChc532->Sumw2();
    listOfProfile->Add(multProfile.fChc532);
  }

  if (fuThreeParticleCorrelation0Gap) {
    multProfile.fChc422_Gap0A = new TProfile(Form("fChc422_Gap0A%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap0A->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap0A);

    multProfile.fChc422_Gap0B = new TProfile(Form("fChc422_Gap0B%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap0B->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap0B);

    multProfile.fChc532_Gap0A = new TProfile(Form("fChc532_Gap0A%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap0A->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap0A);

    multProfile.fChc532_Gap0B = new TProfile(Form("fChc532_Gap0B%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap0B->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap0B);
  }

  if (fuThreeParticleCorrelationGapScan) {
    multProfile.fChc422_Gap2A = new TProfile(Form("fChc422_Gap2A%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap2A->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap2A);

    multProfile.fChc422_Gap2B = new TProfile(Form("fChc422_Gap2B%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap2B->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap2B);

    multProfile.fChc532_Gap2A = new TProfile(Form("fChc532_Gap2A%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap2A->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap2A);

    multProfile.fChc532_Gap2B = new TProfile(Form("fChc532_Gap2B%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap2B->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap2B);

    multProfile.fChc422_Gap4A = new TProfile(Form("fChc422_Gap4A%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap4A->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap4A);

    multProfile.fChc422_Gap4B = new TProfile(Form("fChc422_Gap4B%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap4B->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap4B);

    multProfile.fChc532_Gap4A = new TProfile(Form("fChc532_Gap4A%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap4A->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap4A);

    multProfile.fChc532_Gap4B = new TProfile(Form("fChc532_Gap4B%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap4B->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap4B);

    multProfile.fChc422_Gap6A = new TProfile(Form("fChc422_Gap6A%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap6A->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap6A);

    multProfile.fChc422_Gap6B = new TProfile(Form("fChc422_Gap6B%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap6B->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap6B);

    multProfile.fChc532_Gap6A = new TProfile(Form("fChc532_Gap6A%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap6A->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap6A);

    multProfile.fChc532_Gap6B = new TProfile(Form("fChc532_Gap6B%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap6B->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap6B);

    multProfile.fChc422_Gap8A = new TProfile(Form("fChc422_Gap8A%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap8A->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap8A);

    multProfile.fChc422_Gap8B = new TProfile(Form("fChc422_Gap8B%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap8B->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap8B);

    multProfile.fChc532_Gap8A = new TProfile(Form("fChc532_Gap8A%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap8A->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap8A);

    multProfile.fChc532_Gap8B = new TProfile(Form("fChc532_Gap8B%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap8B->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap8B);
  }

  if (fuThreeParticleCorrelationLargeGap) {
    multProfile.fChc422_Gap10A = new TProfile(Form("fChc422_Gap10A%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap10A->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap10A);

    multProfile.fChc422_Gap10B = new TProfile(Form("fChc422_Gap10B%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_Gap10B->Sumw2();
    listOfProfile->Add(multProfile.fChc422_Gap10B);

    multProfile.fChc532_Gap10A = new TProfile(Form("fChc532_Gap10A%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap10A->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap10A);

    multProfile.fChc532_Gap10B = new TProfile(Form("fChc532_Gap10B%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_Gap10B->Sumw2();
    listOfProfile->Add(multProfile.fChc532_Gap10B);
  }

  if (fuThreeParticleCorrelationThreeSub) {
    multProfile.fChc422_3subL = new TProfile(Form("fChc422_3subL%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_3subL->Sumw2();
    listOfProfile->Add(multProfile.fChc422_3subL);

    multProfile.fChc422_3subM = new TProfile(Form("fChc422_3subM%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_3subM->Sumw2();
    listOfProfile->Add(multProfile.fChc422_3subM);

    multProfile.fChc422_3subR = new TProfile(Form("fChc422_3subR%s", label.Data()), "", nn, xbins);
    multProfile.fChc422_3subR->Sumw2();
    listOfProfile->Add(multProfile.fChc422_3subR);

    multProfile.fChc532_3subLA = new TProfile(Form("fChc532_3subLA%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_3subLA->Sumw2();
    listOfProfile->Add(multProfile.fChc532_3subLA);

    multProfile.fChc532_3subLB = new TProfile(Form("fChc532_3subLB%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_3subLB->Sumw2();
    listOfProfile->Add(multProfile.fChc532_3subLB);

    multProfile.fChc532_3subMA = new TProfile(Form("fChc532_3subMA%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_3subMA->Sumw2();
    listOfProfile->Add(multProfile.fChc532_3subMA);

    multProfile.fChc532_3subMB = new TProfile(Form("fChc532_3subMB%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_3subMB->Sumw2();
    listOfProfile->Add(multProfile.fChc532_3subMB);

    multProfile.fChc532_3subRA = new TProfile(Form("fChc532_3subRA%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_3subRA->Sumw2();
    listOfProfile->Add(multProfile.fChc532_3subRA);

    multProfile.fChc532_3subRB = new TProfile(Form("fChc532_3subRB%s", label.Data()), "", nn, xbins);
    multProfile.fChc532_3subRB->Sumw2();
    listOfProfile->Add(multProfile.fChc532_3subRB);
  }

  // SC(n,m): SC(3,2)
  if (fuFourParticleCorrelationStandard) {
    multProfile.fChsc3232 = new TProfile(Form("fChsc3232%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232);
  }

  if (fuFourParticleCorrelation0Gap) {
    multProfile.fChsc3232_Gap0 = new TProfile(Form("fChsc3232_Gap0%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_Gap0->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_Gap0);

    multProfile.fChsc3232_Gap2 = new TProfile(Form("fChsc3232_Gap2%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_Gap2->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_Gap2);

    multProfile.fChsc3232_Gap4 = new TProfile(Form("fChsc3232_Gap4%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_Gap4->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_Gap4);

    multProfile.fChsc3232_Gap6 = new TProfile(Form("fChsc3232_Gap6%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_Gap6->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_Gap6);

    multProfile.fChsc3232_Gap8 = new TProfile(Form("fChsc3232_Gap8%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_Gap8->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_Gap8);
  }

  if (fuFourParticleCorrelationLargeGap) {
    multProfile.fChsc3232_Gap10 = new TProfile(Form("fChsc3232_Gap10%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_Gap10->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_Gap10);
  }

  if (fuFourParticleCorrelationThreeSub) {
    multProfile.fChsc3232_3subMMLRA = new TProfile(Form("fChsc3232_3subMMLRA%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_3subMMLRA->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_3subMMLRA);

    multProfile.fChsc3232_3subMMLRB = new TProfile(Form("fChsc3232_3subMMLRB%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_3subMMLRB->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_3subMMLRB);

    multProfile.fChsc3232_3subLLMRA = new TProfile(Form("fChsc3232_3subLLMRA%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_3subLLMRA->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_3subLLMRA);

    multProfile.fChsc3232_3subLLMRB = new TProfile(Form("fChsc3232_3subLLMRB%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_3subLLMRB->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_3subLLMRB);

    multProfile.fChsc3232_3subRRMLA = new TProfile(Form("fChsc3232_3subRRMLA%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_3subRRMLA->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_3subRRMLA);

    multProfile.fChsc3232_3subRRMLB = new TProfile(Form("fChsc3232_3subRRMLB%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc3232_3subRRMLB->Sumw2();
    listOfProfile->Add(multProfile.fChsc3232_3subRRMLB);
  }

  // SC(n,m): SC(4,2)
  if (fuFourParticleCorrelationStandard) {
    multProfile.fChsc4242 = new TProfile(Form("fChsc4242%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242);
  }

  if (fuFourParticleCorrelation0Gap) {
    multProfile.fChsc4242_Gap0 = new TProfile(Form("fChsc4242_Gap0%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_Gap0->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_Gap0);
  }

  if (fuFourParticleCorrelationGapScan) {
    multProfile.fChsc4242_Gap2 = new TProfile(Form("fChsc4242_Gap2%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_Gap2->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_Gap2);

    multProfile.fChsc4242_Gap4 = new TProfile(Form("fChsc4242_Gap4%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_Gap4->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_Gap4);

    multProfile.fChsc4242_Gap6 = new TProfile(Form("fChsc4242_Gap6%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_Gap6->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_Gap6);

    multProfile.fChsc4242_Gap8 = new TProfile(Form("fChsc4242_Gap8%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_Gap8->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_Gap8);
  }

  if (fuFourParticleCorrelationLargeGap) {
    multProfile.fChsc4242_Gap10 = new TProfile(Form("fChsc4242_Gap10%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_Gap10->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_Gap10);
  }

  if (fuFourParticleCorrelationThreeSub) {
    multProfile.fChsc4242_3subMMLRA = new TProfile(Form("fChsc4242_3subMMLRA%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_3subMMLRA->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_3subMMLRA);

    multProfile.fChsc4242_3subMMLRB = new TProfile(Form("fChsc4242_3subMMLRB%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_3subMMLRB->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_3subMMLRB);

    multProfile.fChsc4242_3subLLMRA = new TProfile(Form("fChsc4242_3subLLMRA%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_3subLLMRA->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_3subLLMRA);

    multProfile.fChsc4242_3subLLMRB = new TProfile(Form("fChsc4242_3subLLMRB%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_3subLLMRB->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_3subLLMRB);

    multProfile.fChsc4242_3subRRMLA = new TProfile(Form("fChsc4242_3subRRMLA%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_3subRRMLA->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_3subRRMLA);

    multProfile.fChsc4242_3subRRMLB = new TProfile(Form("fChsc4242_3subRRMLB%s", label.Data()), "# of tracks", nn, xbins);
    multProfile.fChsc4242_3subRRMLB->Sumw2();
    listOfProfile->Add(multProfile.fChsc4242_3subRRMLB);
  }
}

Bool_t AliAnalysisTaskXDeptFlow::AcceptAOD(AliAODEvent *inEv) {
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
    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
  }

  if (fPeriod.EqualTo("LHC15o_pass2")) {
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

  bootstrap_value = (((int)(aodVtxZ * 233)) % 30 + 30) % 30;
  return kTRUE;
}

Bool_t AliAnalysisTaskXDeptFlow::AcceptAODTrack(AliAODTrack *mtr, Double_t *ltrackXYZ, Double_t *vtxp) {
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

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC17n")) { // Only for LHC15o pass1
    return fGFWSelection15o->AcceptTrack(mtr,ltrackXYZ,0,kFALSE);
  } else {
    return fGFWSelection->AcceptTrack(mtr,ltrackXYZ,0,kFALSE);
  }
}

Bool_t AliAnalysisTaskXDeptFlow::AcceptMCTruthTrack(AliAODMCParticle *mtrk) {
  // Pt cut
  if(mtrk->Pt() < fMinPt) return kFALSE;
  if(mtrk->Pt() > fMaxPt) return kFALSE;

  if(TMath::Abs(mtrk->Eta()) > fEtaCut) return kFALSE;

  if (!(mtrk->IsPhysicalPrimary())) return kFALSE;
  if (mtrk->Charge() == 0) return kFALSE;
  return kTRUE;
}


void AliAnalysisTaskXDeptFlow::CalculateProfile(PhysicsProfile& profile, double Ntrks) {
  //..calculate 2-particle correlations
  //..................................
  double Dn2 = 0, Dn2Gap0 = 0, Dn2Gap10 = 0, Dn2Gap14 = 0, Dn2_3subLM = 0, Dn2_3subRM = 0, Dn2_3subLR = 0;
  if (fuTwoParticleCorrelationStandard || fuTwoParticleCorrelationHigherStandard) {
    Dn2 = correlator.Two(0, 0).Re();
  }
  if (fuTwoParticleCorrelation0Gap || fuTwoParticleCorrelationHigher0Gap) {
    Dn2Gap0 = correlator.TwoGap0(0, 0).Re();
  }
  if (fuTwoParticleCorrelationLargeGap || fuTwoParticleCorrelationHigherLargeGap) {
    Dn2Gap10 = correlator.TwoGap10(0, 0).Re();
    Dn2Gap14 = correlator.TwoGap14(0, 0).Re();
  }
  if (fuTwoParticleCorrelationThreeSub || fuTwoParticleCorrelationThreeSub) {
    Dn2_3subLM = correlator.Two_3SubLM(0, 0).Re();
    Dn2_3subRM = correlator.Two_3SubRM(0, 0).Re();
    Dn2_3subLR = correlator.Two_3SubLR(0, 0).Re();
  }

  if (fuTwoParticleCorrelationStandard) {
    if(NtrksAfter > 1 && Dn2 != 0) {
      //..v2{2} = <cos2(phi1 - phi2)>
      TComplex v22 = correlator.Two(2, -2);
      double v22Re = v22.Re()/Dn2;
      profile.fChcn2[0]->Fill(Ntrks, v22Re, Dn2);

      //..v3{2} = <cos3(phi1 - phi2)>
      TComplex v32 = correlator.Two(3, -3);
      double v32Re = v32.Re()/Dn2;
      profile.fChcn2[1]->Fill(Ntrks, v32Re, Dn2);

      //..v4{2} = <cos4(phi1 - phi2)>
      TComplex v42 = correlator.Two(4, -4);
      double v42Re = v42.Re()/Dn2;
      profile.fChcn2[2]->Fill(Ntrks, v42Re, Dn2);

      //..v5{2} = <cos5(phi1 - phi2)>
      TComplex v52 = correlator.Two(5, -5);
      double v52Re = v52.Re()/Dn2;
      profile.fChcn2[3]->Fill(Ntrks, v52Re, Dn2);

      //..v6{2} = <cos6(phi1 - phi2)>
      TComplex v62 = correlator.Two(6, -6);
      double v62Re = v62.Re()/Dn2;
      profile.fChcn2[4]->Fill(Ntrks, v62Re, Dn2);
    }
  }

  if (fuTwoParticleCorrelation0Gap) {
    if(NtrksAfterGap0M > 0 && NtrksAfterGap0P > 0 && Dn2Gap0 != 0)
    {
      //..v2{2} with eta Gap > 1.0
      TComplex v22Gap0 = correlator.TwoGap0(2, -2);
      double v22ReGap0 = v22Gap0.Re()/Dn2Gap0;
      profile.fChcn2_Gap0[0]->Fill(Ntrks, v22ReGap0, Dn2Gap0);

      //..v3{2} with eta Gap > 1.0
      TComplex v32Gap0 = correlator.TwoGap0(3, -3);
      double v32ReGap0 = v32Gap0.Re()/Dn2Gap0;
      profile.fChcn2_Gap0[1]->Fill(Ntrks, v32ReGap0, Dn2Gap0);

      //..v4{2} with eta Gap > 1.0
      TComplex v42Gap0 = correlator.TwoGap0(4, -4);
      double v42ReGap0 = v42Gap0.Re()/Dn2Gap0;
      profile.fChcn2_Gap0[2]->Fill(Ntrks, v42ReGap0, Dn2Gap0);

      //..v5{2} with eta Gap > 1.0
      TComplex v52Gap0 = correlator.TwoGap0(5, -5);
      double v52ReGap0 = v52Gap0.Re()/Dn2Gap0;
      profile.fChcn2_Gap0[3]->Fill(Ntrks, v52ReGap0, Dn2Gap0);

      //..v6{2} with eta Gap > 1.0
      TComplex v62Gap0 = correlator.TwoGap0(6, -6);
      double v62ReGap0 = v62Gap0.Re()/Dn2Gap0;
      profile.fChcn2_Gap0[4]->Fill(Ntrks, v62ReGap0, Dn2Gap0);
    }
  }

  if (fuTwoParticleCorrelationLargeGap) {
    if(NtrksAfterGap10M > 0 && NtrksAfterGap10P > 0 && Dn2Gap10 != 0)
    {
      //..v2{2} with eta Gap > 1.0
      TComplex v22Gap10 = correlator.TwoGap10(2, -2);
      double v22ReGap10 = v22Gap10.Re()/Dn2Gap10;
      profile.fChcn2_Gap10[0]->Fill(Ntrks, v22ReGap10, Dn2Gap10);

      //..v3{2} with eta Gap > 1.0
      TComplex v32Gap10 = correlator.TwoGap10(3, -3);
      double v32ReGap10 = v32Gap10.Re()/Dn2Gap10;
      profile.fChcn2_Gap10[1]->Fill(Ntrks, v32ReGap10, Dn2Gap10);

      //..v4{2} with eta Gap > 1.0
      TComplex v42Gap10 = correlator.TwoGap10(4, -4);
      double v42ReGap10 = v42Gap10.Re()/Dn2Gap10;
      profile.fChcn2_Gap10[2]->Fill(Ntrks, v42ReGap10, Dn2Gap10);

      //..v5{2} with eta Gap > 1.0
      TComplex v52Gap10 = correlator.TwoGap10(5, -5);
      double v52ReGap10 = v52Gap10.Re()/Dn2Gap10;
      profile.fChcn2_Gap10[3]->Fill(Ntrks, v52ReGap10, Dn2Gap10);

      //..v6{2} with eta Gap > 1.0
      TComplex v62Gap10 = correlator.TwoGap10(6, -6);
      double v62ReGap10 = v62Gap10.Re()/Dn2Gap10;
      profile.fChcn2_Gap10[4]->Fill(Ntrks, v62ReGap10, Dn2Gap10);
    }

    if(NtrksAfterGap14M > 0 && NtrksAfterGap14P > 0 && Dn2Gap14 != 0)
    {
      //..v2{2} with eta Gap > 1.4
      TComplex v22Gap14 = correlator.TwoGap14(2, -2);
      double v22ReGap14 = v22Gap14.Re()/Dn2Gap14;
      profile.fChcn2_Gap14[0]->Fill(Ntrks, v22ReGap14, Dn2Gap14);
    }
  }

  //..for 3-subevent method, Gap0
  if (fuTwoParticleCorrelationThreeSub) {
    if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && Dn2_3subLM != 0)
    {//..left+middle
      TComplex v22_3subLM = correlator.Two_3SubLM(2, -2);
      double v22Re_3subLM = v22_3subLM.Re()/Dn2_3subLM;
      profile.fChcn2_3subLM[0]->Fill(Ntrks, v22Re_3subLM, Dn2_3subLM);

      TComplex v32_3subLM = correlator.Two_3SubLM(3, -3);
      double v32Re_3subLM = v32_3subLM.Re()/Dn2_3subLM;
      profile.fChcn2_3subLM[1]->Fill(Ntrks, v32Re_3subLM, Dn2_3subLM);

      TComplex v42_3subLM = correlator.Two_3SubLM(4, -4);
      double v42Re_3subLM = v42_3subLM.Re()/Dn2_3subLM;
      profile.fChcn2_3subLM[2]->Fill(Ntrks, v42Re_3subLM, Dn2_3subLM);
    }

    if(NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn2_3subRM != 0)
    {//..right+middle
      TComplex v22_3subRM = correlator.Two_3SubRM(2, -2);
      double v22Re_3subRM = v22_3subRM.Re()/Dn2_3subRM;
      profile.fChcn2_3subRM[0]->Fill(Ntrks, v22Re_3subRM, Dn2_3subRM);

      TComplex v32_3subRM = correlator.Two_3SubRM(3, -3);
      double v32Re_3subRM = v32_3subRM.Re()/Dn2_3subRM;
      profile.fChcn2_3subRM[1]->Fill(Ntrks, v32Re_3subRM, Dn2_3subRM);

      TComplex v42_3subRM = correlator.Two_3SubRM(4, -4);
      double v42Re_3subRM = v42_3subRM.Re()/Dn2_3subRM;
      profile.fChcn2_3subRM[2]->Fill(Ntrks, v42Re_3subRM, Dn2_3subRM);
    }

    if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && Dn2_3subLR != 0)
    {//..left+right
      TComplex v22_3subLR = correlator.Two_3SubLR(2, -2);
      double v22Re_3subLR = v22_3subLR.Re()/Dn2_3subLR;
      profile.fChcn2_3subLR[0]->Fill(Ntrks, v22Re_3subLR, Dn2_3subLR);

      TComplex v32_3subLR = correlator.Two_3SubLR(3, -3);
      double v32Re_3subLR = v32_3subLR.Re()/Dn2_3subLR;
      profile.fChcn2_3subLR[1]->Fill(Ntrks, v32Re_3subLR, Dn2_3subLR);

      TComplex v42_3subLR = correlator.Two_3SubLR(4, -4);
      double v42Re_3subLR = v42_3subLR.Re()/Dn2_3subLR;
      profile.fChcn2_3subLR[2]->Fill(Ntrks, v42Re_3subLR, Dn2_3subLR);
    }
  }

  //..calculate 3-particle correlations
  //................................
  double Dn3 = 0, Dn3Gap0A = 0, Dn3Gap0B = 0, Dn3Gap2A = 0, Dn3Gap2B = 0, Dn3Gap4A = 0, Dn3Gap4B = 0, Dn3Gap6A = 0,
         Dn3Gap6B = 0, Dn3Gap8A = 0, Dn3Gap8B = 0, Dn3Gap10A = 0, Dn3Gap10B = 0, Dn3_3sub = 0;
  if (fuThreeParticleCorrelationStandard) {
    Dn3 = correlator.Three(0, 0, 0).Re();
  }
  if (fuThreeParticleCorrelation0Gap) {
    Dn3Gap0A = correlator.ThreeGap0A(0, 0, 0).Re();
    Dn3Gap0B = correlator.ThreeGap0B(0, 0, 0).Re();
  }
  if (fuThreeParticleCorrelationGapScan) {
    Dn3Gap2A = correlator.ThreeGap2A(0, 0, 0).Re();
    Dn3Gap2B = correlator.ThreeGap2B(0, 0, 0).Re();
    Dn3Gap4A = correlator.ThreeGap4A(0, 0, 0).Re();
    Dn3Gap4B = correlator.ThreeGap4B(0, 0, 0).Re();
    Dn3Gap6A = correlator.ThreeGap6A(0, 0, 0).Re();
    Dn3Gap6B = correlator.ThreeGap6B(0, 0, 0).Re();
    Dn3Gap8A = correlator.ThreeGap8A(0, 0, 0).Re();
    Dn3Gap8B = correlator.ThreeGap8B(0, 0, 0).Re();
  }
  if (fuThreeParticleCorrelationLargeGap) {
    Dn3Gap10A = correlator.ThreeGap10A(0, 0, 0).Re();
    Dn3Gap10B = correlator.ThreeGap10B(0, 0, 0).Re();
  }
  if (fuThreeParticleCorrelationThreeSub) {
    Dn3_3sub = correlator.Three_3Sub(0,0,0).Re();
  }

  if (fuThreeParticleCorrelationStandard) {
    if(NtrksAfter > 2 && Dn3 != 0)
    {
      //..v4{psi2}
      TComplex v422 = correlator.Three(4, -2, -2);
      double v422Re = v422.Re()/Dn3;
      profile.fChc422->Fill(Ntrks, v422Re, Dn3);

      //..v5{psi32}
      TComplex v532 = correlator.Three(5, -3, -2);
      double v532Re = v532.Re()/Dn3;
      profile.fChc532->Fill(Ntrks, v532Re, Dn3);
    }
  }

  // A-type
  if (fuThreeParticleCorrelation0Gap) {
    if(NtrksAfterGap0M > 0 && NtrksAfterGap0P > 1 && Dn3Gap0A != 0) {

      TComplex v422Gap0A = correlator.ThreeGap0A(4, -2, -2);
      double v422Gap0ARe = v422Gap0A.Re()/Dn3Gap0A;
      profile.fChc422_Gap0A->Fill(Ntrks, v422Gap0ARe, Dn3Gap0A);

      TComplex v532Gap0A = correlator.ThreeGap0A(5, -3, -2);
      double v532Gap0ARe = v532Gap0A.Re()/Dn3Gap0A;
      profile.fChc532_Gap0A->Fill(Ntrks, v532Gap0ARe, Dn3Gap0A);
    }

    // B-type
    if(NtrksAfterGap0P > 0 && NtrksAfterGap0M > 1 && Dn3Gap0B != 0)
    {

      TComplex v422Gap0B = correlator.ThreeGap0B(4, -2, -2);
      double v422Gap0BRe = v422Gap0B.Re()/Dn3Gap0B;
      profile.fChc422_Gap0B->Fill(Ntrks, v422Gap0BRe, Dn3Gap0B);

      TComplex v532Gap0B = correlator.ThreeGap0B(5, -3, -2);
      double v532Gap0BRe = v532Gap0B.Re()/Dn3Gap0B;
      profile.fChc532_Gap0B->Fill(Ntrks, v532Gap0BRe, Dn3Gap0B);
    }
  }

  if (fuThreeParticleCorrelationGapScan) {
    // A-type
    if(NtrksAfterGap2M > 0 && NtrksAfterGap2P > 1 && Dn3Gap2A != 0)
    {

      TComplex v422Gap2A = correlator.ThreeGap2A(4, -2, -2);
      double v422Gap2ARe = v422Gap2A.Re()/Dn3Gap2A;
      profile.fChc422_Gap2A->Fill(Ntrks, v422Gap2ARe, Dn3Gap2A);

      TComplex v532Gap2A = correlator.ThreeGap2A(5, -3, -2);
      double v532Gap2ARe = v532Gap2A.Re()/Dn3Gap2A;
      profile.fChc532_Gap2A->Fill(Ntrks, v532Gap2ARe, Dn3Gap2A);
    }

    // B-type
    if(NtrksAfterGap2P > 0 && NtrksAfterGap2M > 1 && Dn3Gap2B != 0)
    {

      TComplex v422Gap2B = correlator.ThreeGap2B(4, -2, -2);
      double v422Gap2BRe = v422Gap2B.Re()/Dn3Gap2B;
      profile.fChc422_Gap2B->Fill(Ntrks, v422Gap2BRe, Dn3Gap2B);

      TComplex v532Gap2B = correlator.ThreeGap2B(5, -3, -2);
      double v532Gap2BRe = v532Gap2B.Re()/Dn3Gap2B;
      profile.fChc532_Gap2B->Fill(Ntrks, v532Gap2BRe, Dn3Gap2B);
    }

    // A-type
    if(NtrksAfterGap4M > 0 && NtrksAfterGap4P > 1 && Dn3Gap4A != 0)
    {

      TComplex v422Gap4A = correlator.ThreeGap4A(4, -2, -2);
      double v422Gap4ARe = v422Gap4A.Re()/Dn3Gap4A;
      profile.fChc422_Gap4A->Fill(Ntrks, v422Gap4ARe, Dn3Gap4A);

      TComplex v532Gap4A = correlator.ThreeGap4A(5, -3, -2);
      double v532Gap4ARe = v532Gap4A.Re()/Dn3Gap4A;
      profile.fChc532_Gap4A->Fill(Ntrks, v532Gap4ARe, Dn3Gap4A);
    }

    // B-type
    if(NtrksAfterGap4P > 0 && NtrksAfterGap4M > 1 && Dn3Gap4B != 0)
    {

      TComplex v422Gap4B = correlator.ThreeGap4B(4, -2, -2);
      double v422Gap4BRe = v422Gap4B.Re()/Dn3Gap4B;
      profile.fChc422_Gap4B->Fill(Ntrks, v422Gap4BRe, Dn3Gap4B);

      TComplex v532Gap4B = correlator.ThreeGap4B(5, -3, -2);
      double v532Gap4BRe = v532Gap4B.Re()/Dn3Gap4B;
      profile.fChc532_Gap4B->Fill(Ntrks, v532Gap4BRe, Dn3Gap4B);
    }

    // A-type
    if(NtrksAfterGap6M > 0 && NtrksAfterGap6P > 1 && Dn3Gap6A != 0)
    {

      TComplex v422Gap6A = correlator.ThreeGap6A(4, -2, -2);
      double v422Gap6ARe = v422Gap6A.Re()/Dn3Gap6A;
      profile.fChc422_Gap6A->Fill(Ntrks, v422Gap6ARe, Dn3Gap6A);

      TComplex v532Gap6A = correlator.ThreeGap6A(5, -3, -2);
      double v532Gap6ARe = v532Gap6A.Re()/Dn3Gap6A;
      profile.fChc532_Gap6A->Fill(Ntrks, v532Gap6ARe, Dn3Gap6A);
    }

    // B-type
    if(NtrksAfterGap6P > 0 && NtrksAfterGap6M > 1 && Dn3Gap6B != 0)
    {

      TComplex v422Gap6B = correlator.ThreeGap6B(4, -2, -2);
      double v422Gap6BRe = v422Gap6B.Re()/Dn3Gap6B;
      profile.fChc422_Gap6B->Fill(Ntrks, v422Gap6BRe, Dn3Gap6B);

      TComplex v532Gap6B = correlator.ThreeGap6B(5, -3, -2);
      double v532Gap6BRe = v532Gap6B.Re()/Dn3Gap6B;
      profile.fChc532_Gap6B->Fill(Ntrks, v532Gap6BRe, Dn3Gap6B);
    }

    // A-type
    if(NtrksAfterGap8M > 0 && NtrksAfterGap8P > 1 && Dn3Gap8A != 0)
    {

      TComplex v422Gap8A = correlator.ThreeGap8A(4, -2, -2);
      double v422Gap8ARe = v422Gap8A.Re()/Dn3Gap8A;
      profile.fChc422_Gap8A->Fill(Ntrks, v422Gap8ARe, Dn3Gap8A);

      TComplex v532Gap8A = correlator.ThreeGap8A(5, -3, -2);
      double v532Gap8ARe = v532Gap8A.Re()/Dn3Gap8A;
      profile.fChc532_Gap8A->Fill(Ntrks, v532Gap8ARe, Dn3Gap8A);
    }

    // B-type
    if(NtrksAfterGap8P > 0 && NtrksAfterGap8M > 1 && Dn3Gap8B != 0)
    {

      TComplex v422Gap8B = correlator.ThreeGap8B(4, -2, -2);
      double v422Gap8BRe = v422Gap8B.Re()/Dn3Gap8B;
      profile.fChc422_Gap8B->Fill(Ntrks, v422Gap8BRe, Dn3Gap8B);

      TComplex v532Gap8B = correlator.ThreeGap8B(5, -3, -2);
      double v532Gap8BRe = v532Gap8B.Re()/Dn3Gap8B;
      profile.fChc532_Gap8B->Fill(Ntrks, v532Gap8BRe, Dn3Gap8B);
    }
  }

  if (fuThreeParticleCorrelationLargeGap) {
    // A-type
    if(NtrksAfterGap10M > 0 && NtrksAfterGap10P > 1 && Dn3Gap10A != 0)
    {

      TComplex v422Gap10A = correlator.ThreeGap10A(4, -2, -2);
      double v422Gap10ARe = v422Gap10A.Re()/Dn3Gap10A;
      profile.fChc422_Gap10A->Fill(Ntrks, v422Gap10ARe, Dn3Gap10A);

      TComplex v532Gap10A = correlator.ThreeGap10A(5, -3, -2);
      double v532Gap10ARe = v532Gap10A.Re()/Dn3Gap10A;
      profile.fChc532_Gap10A->Fill(Ntrks, v532Gap10ARe, Dn3Gap10A);
    }

    // B-type
    if(NtrksAfterGap10P > 0 && NtrksAfterGap10M > 1 && Dn3Gap10B != 0)
    {

      TComplex v422Gap10B = correlator.ThreeGap10B(4, -2, -2);
      double v422Gap10BRe = v422Gap10B.Re()/Dn3Gap10B;
      profile.fChc422_Gap10B->Fill(Ntrks, v422Gap10BRe, Dn3Gap10B);

      TComplex v532Gap10B = correlator.ThreeGap10B(5, -3, -2);
      double v532Gap10BRe = v532Gap10B.Re()/Dn3Gap10B;
      profile.fChc532_Gap10B->Fill(Ntrks, v532Gap10BRe, Dn3Gap10B);
    }
  }

  if (fuThreeParticleCorrelationThreeSub) {
    //..3-subevent method
    if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn3_3sub != 0)
    {
      // v422
      TComplex v422_3subL = correlator.Three_3Sub(4, -2, -2);
      double v422_3subLRe = v422_3subL.Re()/Dn3_3sub;
      profile.fChc422_3subL->Fill(Ntrks, v422_3subLRe, Dn3_3sub);

      TComplex v422_3subM = correlator.Three_3Sub(-2, 4, -2);
      double v422_3subMRe = v422_3subM.Re()/Dn3_3sub;
      profile.fChc422_3subM->Fill(Ntrks, v422_3subMRe, Dn3_3sub);

      TComplex v422_3subR = correlator.Three_3Sub(4, -2, -2);
      double v422_3subRRe = v422_3subR.Re()/Dn3_3sub;
      profile.fChc422_3subR->Fill(Ntrks, v422_3subRRe, Dn3_3sub);

      // v532
      TComplex v532_3subLA = correlator.Three_3Sub(5, -3, -2);
      double v532_3subLARe = v532_3subLA.Re()/Dn3_3sub;
      profile.fChc532_3subLA->Fill(Ntrks, v532_3subLARe, Dn3_3sub);

      TComplex v532_3subLB = correlator.Three_3Sub(5, -2, -3);
      double v532_3subLBRe = v532_3subLB.Re()/Dn3_3sub;
      profile.fChc532_3subLB->Fill(Ntrks, v532_3subLBRe, Dn3_3sub);

      TComplex v532_3subMA = correlator.Three_3Sub(-3, 5, -2);
      double v532_3subMARe = v532_3subMA.Re()/Dn3_3sub;
      profile.fChc532_3subMA->Fill(Ntrks, v532_3subMARe, Dn3_3sub);

      TComplex v532_3subMB = correlator.Three_3Sub(-2, 5, -3);
      double v532_3subMBRe = v532_3subMB.Re()/Dn3_3sub;
      profile.fChc532_3subMB->Fill(Ntrks, v532_3subMBRe, Dn3_3sub);

      TComplex v532_3subRA = correlator.Three_3Sub(-2, -3, 5);
      double v532_3subRARe = v532_3subRA.Re()/Dn3_3sub;
      profile.fChc532_3subRA->Fill(Ntrks, v532_3subRARe, Dn3_3sub);

      TComplex v532_3subRB = correlator.Three_3Sub(-3, -2, 5);
      double v532_3subRBRe = v532_3subRB.Re()/Dn3_3sub;
      profile.fChc532_3subRB->Fill(Ntrks, v532_3subRBRe, Dn3_3sub);
    }
  }


  //..calculate 4-particle correlations
  //................................


  double Dn4, Dn4Gap0, Dn4Gap2, Dn4Gap4, Dn4Gap6, Dn4Gap8, Dn4Gap10, Dn4_3subMMLR, Dn4_3subLLMR, Dn4_3subRRML;
  Dn4 = correlator.Four(0, 0, 0, 0).Re();
  Dn4Gap0 = correlator.FourGap0(0, 0, 0, 0).Re();
  Dn4Gap2 = correlator.FourGap2(0, 0, 0, 0).Re();
  Dn4Gap4 = correlator.FourGap4(0, 0, 0, 0).Re();
  Dn4Gap6 = correlator.FourGap6(0, 0, 0, 0).Re();
  Dn4Gap8 = correlator.FourGap8(0, 0, 0, 0).Re();
  Dn4Gap10 = correlator.FourGap10(0, 0, 0, 0).Re();
  Dn4_3subMMLR = correlator.Four_3SubMMLR(0, 0, 0, 0).Re();
  Dn4_3subLLMR = correlator.Four_3SubLLMR(0, 0, 0, 0).Re();
  Dn4_3subRRML = correlator.Four_3SubRRML(0, 0, 0, 0).Re();

  if (fuFourParticleCorrelationStandard) {
    if(NtrksAfter > 3 && Dn4 != 0)
    {

      TComplex v24 = correlator.Four(2, 2, -2, -2);
      double v24Re = v24.Re()/Dn4;
      profile.fChcn4[0]->Fill(Ntrks, v24Re, Dn4);
      // fcn4Ntrks1bin[0][fBin]->Fill(NtrksAfter, v24Re, Dn4);

      // TComplex v34 = correlator.Four(3, 3, -3, -3);
      // double v34Re = v34.Re()/Dn4;
      // profile.fChcn4[1]->Fill(Ntrks, v34Re, Dn4);
      // fcn4Ntrks1bin[1][fBin]->Fill(NtrksAfter, v34Re, Dn4);

      // TComplex v44 = correlator.Four(4, 4, -4, -4);
      // double v44Re = v44.Re()/Dn4;
      // profile.fChcn4[2]->Fill(Ntrks, v44Re, Dn4);
      // fcn4Ntrks1bin[2][fBin]->Fill(NtrksAfter, v44Re, Dn4);

      //..SC(3,2,-3,-2)
      TComplex sc3232 = correlator.Four(3, 2, -3, -2);
      double sc3232Re = sc3232.Re()/Dn4;
      profile.fChsc3232->Fill(Ntrks, sc3232Re, Dn4);

      //..SC(4,2,-4,-2)
      TComplex sc4242 = correlator.Four(4, 2, -4, -2);
      double sc4242Re = sc4242.Re()/Dn4;
      profile.fChsc4242->Fill(Ntrks, sc4242Re, Dn4);
    }
  }

  if (fuFourParticleCorrelation0Gap) {
    if(NtrksAfterGap0M > 1 && NtrksAfterGap0P > 1 && Dn4Gap0 !=0)
    {
      TComplex v24Gap0 = correlator.FourGap0(2, 2, -2, -2);
      double v24Gap0Re = v24Gap0.Re()/Dn4Gap0;
      profile.fChcn4_Gap0[0]->Fill(Ntrks, v24Gap0Re, Dn4Gap0);

      // TComplex v34Gap0 = correlator.FourGap0(3, 3, -3, -3);
      // double v34Gap0Re = v34Gap0.Re()/Dn4Gap0;
      // profile.fChcn4_Gap0[1]->Fill(Ntrks, v34Gap0Re, Dn4Gap0);

      // TComplex v44Gap0 = correlator.FourGap0(4, 4, -4, -4);
      // double v44Gap0Re = v44Gap0.Re()/Dn4Gap0;
      // profile.fChcn4_Gap0[2]->Fill(Ntrks, v44Gap0Re, Dn4Gap0);

      TComplex sc3232Gap0 = correlator.FourGap0(3, 2, -3, -2);
      double sc3232Gap0Re = sc3232Gap0.Re()/Dn4Gap0;
      profile.fChsc3232_Gap0->Fill(Ntrks, sc3232Gap0Re, Dn4Gap0);

      TComplex sc4242Gap0 = correlator.FourGap0(4, 2, -4, -2);
      double sc4242Gap0Re = sc4242Gap0.Re()/Dn4Gap0;
      profile.fChsc4242_Gap0->Fill(Ntrks, sc4242Gap0Re, Dn4Gap0);
    }
  }

  if (fuFourParticleCorrelationGapScan) {
    if(NtrksAfterGap2M > 1 && NtrksAfterGap2P > 1 && Dn4Gap2 !=0)
    {
      TComplex v24Gap2 = correlator.FourGap2(2, 2, -2, -2);
      double v24Gap2Re = v24Gap2.Re()/Dn4Gap2;
      profile.fChcn4_Gap2[0]->Fill(Ntrks, v24Gap2Re, Dn4Gap2);

      // TComplex v34Gap2 = correlator.FourGap2(3, 3, -3, -3);
      // double v34Gap2Re = v34Gap2.Re()/Dn4Gap2;
      // profile.fChcn4_Gap2[1]->Fill(Ntrks, v34Gap2Re, Dn4Gap2);

      // TComplex v44Gap2 = correlator.FourGap2(4, 4, -4, -4);
      // double v44Gap2Re = v44Gap2.Re()/Dn4Gap2;
      // profile.fChcn4_Gap2[2]->Fill(Ntrks, v44Gap2Re, Dn4Gap2);

      TComplex sc3232Gap2 = correlator.FourGap2(3, 2, -3, -2);
      double sc3232Gap2Re = sc3232Gap2.Re()/Dn4Gap2;
      profile.fChsc3232_Gap2->Fill(Ntrks, sc3232Gap2Re, Dn4Gap2);

      TComplex sc4242Gap2 = correlator.FourGap2(4, 2, -4, -2);
      double sc4242Gap2Re = sc4242Gap2.Re()/Dn4Gap2;
      profile.fChsc4242_Gap2->Fill(Ntrks, sc4242Gap2Re, Dn4Gap2);
    }

    if(NtrksAfterGap4M > 1 && NtrksAfterGap4P > 1 && Dn4Gap4 !=0)
    {
      TComplex v24Gap4 = correlator.FourGap4(2, 2, -2, -2);
      double v24Gap4Re = v24Gap4.Re()/Dn4Gap4;
      profile.fChcn4_Gap4[0]->Fill(Ntrks, v24Gap4Re, Dn4Gap4);

      // TComplex v34Gap4 = correlator.FourGap4(3, 3, -3, -3);
      // double v34Gap4Re = v34Gap4.Re()/Dn4Gap4;
      // profile.fChcn4_Gap4[1]->Fill(Ntrks, v34Gap4Re, Dn4Gap4);

      // TComplex v44Gap4 = correlator.FourGap4(4, 4, -4, -4);
      // double v44Gap4Re = v44Gap4.Re()/Dn4Gap4;
      // profile.fChcn4_Gap4[2]->Fill(Ntrks, v44Gap4Re, Dn4Gap4);

      TComplex sc3232Gap4 = correlator.FourGap4(3, 2, -3, -2);
      double sc3232Gap4Re = sc3232Gap4.Re()/Dn4Gap4;
      profile.fChsc3232_Gap4->Fill(Ntrks, sc3232Gap4Re, Dn4Gap4);

      TComplex sc4242Gap4 = correlator.FourGap4(4, 2, -4, -2);
      double sc4242Gap4Re = sc4242Gap4.Re()/Dn4Gap4;
      profile.fChsc4242_Gap4->Fill(Ntrks, sc4242Gap4Re, Dn4Gap4);
    }

    if(NtrksAfterGap6M > 1 && NtrksAfterGap6P > 1 && Dn4Gap6 !=0)
    {
      TComplex v24Gap6 = correlator.FourGap6(2, 2, -2, -2);
      double v24Gap6Re = v24Gap6.Re()/Dn4Gap6;
      profile.fChcn4_Gap6[0]->Fill(Ntrks, v24Gap6Re, Dn4Gap6);

      // TComplex v34Gap6 = correlator.FourGap6(3, 3, -3, -3);
      // double v34Gap6Re = v34Gap6.Re()/Dn4Gap6;
      // profile.fChcn4_Gap6[1]->Fill(Ntrks, v34Gap6Re, Dn4Gap6);

      // TComplex v44Gap6 = correlator.FourGap6(4, 4, -4, -4);
      // double v44Gap6Re = v44Gap6.Re()/Dn4Gap6;
      // profile.fChcn4_Gap6[2]->Fill(Ntrks, v44Gap6Re, Dn4Gap6);

      TComplex sc3232Gap6 = correlator.FourGap6(3, 2, -3, -2);
      double sc3232Gap6Re = sc3232Gap6.Re()/Dn4Gap6;
      profile.fChsc3232_Gap6->Fill(Ntrks, sc3232Gap6Re, Dn4Gap6);

      TComplex sc4242Gap6 = correlator.FourGap6(4, 2, -4, -2);
      double sc4242Gap6Re = sc4242Gap6.Re()/Dn4Gap6;
      profile.fChsc4242_Gap6->Fill(Ntrks, sc4242Gap6Re, Dn4Gap6);
    }

    if(NtrksAfterGap8M > 1 && NtrksAfterGap8P > 1 && Dn4Gap8 !=0)
    {
      TComplex v24Gap8 = correlator.FourGap8(2, 2, -2, -2);
      double v24Gap8Re = v24Gap8.Re()/Dn4Gap8;
      profile.fChcn4_Gap8[0]->Fill(Ntrks, v24Gap8Re, Dn4Gap8);

      // TComplex v34Gap8 = correlator.FourGap8(3, 3, -3, -3);
      // double v34Gap8Re = v34Gap8.Re()/Dn4Gap8;
      // profile.fChcn4_Gap8[1]->Fill(Ntrks, v34Gap8Re, Dn4Gap8);

      // TComplex v44Gap8 = correlator.FourGap8(4, 4, -4, -4);
      // double v44Gap8Re = v44Gap8.Re()/Dn4Gap8;
      // profile.fChcn4_Gap8[2]->Fill(Ntrks, v44Gap8Re, Dn4Gap8);

      TComplex sc3232Gap8 = correlator.FourGap8(3, 2, -3, -2);
      double sc3232Gap8Re = sc3232Gap8.Re()/Dn4Gap8;
      profile.fChsc3232_Gap8->Fill(Ntrks, sc3232Gap8Re, Dn4Gap8);

      TComplex sc4242Gap8 = correlator.FourGap8(4, 2, -4, -2);
      double sc4242Gap8Re = sc4242Gap8.Re()/Dn4Gap8;
      profile.fChsc4242_Gap8->Fill(Ntrks, sc4242Gap8Re, Dn4Gap8);
    }
  }

  if (fuFourParticleCorrelationLargeGap) {
    if(NtrksAfterGap10M > 1 && NtrksAfterGap10P > 1 && Dn4Gap10 !=0)
    {
      TComplex v24Gap10 = correlator.FourGap10(2, 2, -2, -2);
      double v24Gap10Re = v24Gap10.Re()/Dn4Gap10;
      profile.fChcn4_Gap10[0]->Fill(Ntrks, v24Gap10Re, Dn4Gap10);

      // TComplex v34Gap10 = correlator.FourGap10(3, 3, -3, -3);
      // double v34Gap10Re = v34Gap10.Re()/Dn4Gap10;
      // profile.fChcn4_Gap10[1]->Fill(Ntrks, v34Gap10Re, Dn4Gap10);

      // TComplex v44Gap10 = correlator.FourGap10(4, 4, -4, -4);
      // double v44Gap10Re = v44Gap10.Re()/Dn4Gap10;
      // profile.fChcn4_Gap10[2]->Fill(Ntrks, v44Gap10Re, Dn4Gap10);

      TComplex sc3232Gap10 = correlator.FourGap10(3, 2, -3, -2);
      double sc3232Gap10Re = sc3232Gap10.Re()/Dn4Gap10;
      profile.fChsc3232_Gap10->Fill(Ntrks, sc3232Gap10Re, Dn4Gap10);

      TComplex sc4242Gap10 = correlator.FourGap10(4, 2, -4, -2);
      double sc4242Gap10Re = sc4242Gap10.Re()/Dn4Gap10;
      profile.fChsc4242_Gap10->Fill(Ntrks, sc4242Gap10Re, Dn4Gap10);
    }
  }

  if (fuFourParticleCorrelationThreeSub) {
    //..3-subevent method
    if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && NtrksAfter3subM > 1 && Dn4_3subMMLR != 0)
    {
      TComplex v24_3sub = correlator.Four_3SubMMLR(2, 2, -2, -2);
      double v24_3subRe = v24_3sub.Re()/Dn4_3subMMLR;
      profile.fChcn4_3subMMLR[0]->Fill(Ntrks, v24_3subRe, Dn4_3subMMLR);

      // TComplex v34_3sub = correlator.Four_3SubMMLR(3, 3, -3, -3);
      // double v34_3subRe = v34_3sub.Re()/Dn4_3subMMLR;
      // profile.fChcn4_3subMMLR[1]->Fill(Ntrks, v34_3subRe, Dn4_3subMMLR);

      // TComplex v44_3sub = correlator.Four_3SubMMLR(4, 4, -4, -4);
      // double v44_3subRe = v44_3sub.Re()/Dn4_3subMMLR;
      // profile.fChcn4_3subMMLR[2]->Fill(Ntrks, v44_3subRe, Dn4_3subMMLR);

      TComplex sc3232_3subA = correlator.Four_3SubMMLR(3, 2, -3, -2);
      double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subMMLR;
      profile.fChsc3232_3subMMLRA->Fill(Ntrks, sc3232_3subARe, Dn4_3subMMLR);

      TComplex sc3232_3subB = correlator.Four_3SubMMLR(3, 2, -2, -3);
      double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subMMLR;
      profile.fChsc3232_3subMMLRB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subMMLR);

      TComplex sc4242_3subA = correlator.Four_3SubMMLR(4, 2, -4, -2);
      double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subMMLR;
      profile.fChsc4242_3subMMLRA->Fill(Ntrks, sc4242_3subARe, Dn4_3subMMLR);

      TComplex sc4242_3subB = correlator.Four_3SubMMLR(4, 2, -2, -4);
      double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subMMLR;
      profile.fChsc4242_3subMMLRB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subMMLR);
    }

    //..3-subevent method
    if(NtrksAfter3subL > 1 && NtrksAfter3subR > 0 && NtrksAfter3subM > 0 && Dn4_3subLLMR != 0)
    {
      TComplex v24_3sub = correlator.Four_3SubLLMR(2, 2, -2, -2);
      double v24_3subRe = v24_3sub.Re()/Dn4_3subLLMR;
      profile.fChcn4_3subLLMR[0]->Fill(Ntrks, v24_3subRe, Dn4_3subLLMR);

      // TComplex v34_3sub = correlator.Four_3SubLLMR(3, 3, -3, -3);
      // double v34_3subRe = v34_3sub.Re()/Dn4_3subLLMR;
      // profile.fChcn4_3subLLMR[1]->Fill(Ntrks, v34_3subRe, Dn4_3subLLMR);

      // TComplex v44_3sub = correlator.Four_3SubLLMR(4, 4, -4, -4);
      // double v44_3subRe = v44_3sub.Re()/Dn4_3subLLMR;
      // profile.fChcn4_3subLLMR[2]->Fill(Ntrks, v44_3subRe, Dn4_3subLLMR);

      TComplex sc3232_3subA = correlator.Four_3SubLLMR(3, 2, -3, -2);
      double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subLLMR;
      profile.fChsc3232_3subLLMRA->Fill(Ntrks, sc3232_3subARe, Dn4_3subLLMR);

      TComplex sc3232_3subB = correlator.Four_3SubLLMR(3, 2, -2, -3);
      double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subLLMR;
      profile.fChsc3232_3subLLMRB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subLLMR);

      TComplex sc4242_3subA = correlator.Four_3SubLLMR(4, 2, -4, -2);
      double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subLLMR;
      profile.fChsc4242_3subLLMRA->Fill(Ntrks, sc4242_3subARe, Dn4_3subLLMR);

      TComplex sc4242_3subB = correlator.Four_3SubLLMR(4, 2, -2, -4);
      double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subLLMR;
      profile.fChsc4242_3subLLMRB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subLLMR);
    }

    //..3-subevent method
    if(NtrksAfter3subL > 0 && NtrksAfter3subR > 1 && NtrksAfter3subM > 0 && Dn4_3subRRML != 0)
    {
      TComplex v24_3sub = correlator.Four_3SubRRML(2, 2, -2, -2);
      double v24_3subRe = v24_3sub.Re()/Dn4_3subRRML;
      profile.fChcn4_3subRRML[0]->Fill(Ntrks, v24_3subRe, Dn4_3subRRML);

      // TComplex v34_3sub = correlator.Four_3SubRRML(3, 3, -3, -3);
      // double v34_3subRe = v34_3sub.Re()/Dn4_3subRRML;
      // profile.fChcn4_3subRRML[1]->Fill(Ntrks, v34_3subRe, Dn4_3subRRML);

      // TComplex v44_3sub = correlator.Four_3SubRRML(4, 4, -4, -4);
      // double v44_3subRe = v44_3sub.Re()/Dn4_3subRRML;
      // profile.fChcn4_3subRRML[2]->Fill(Ntrks, v44_3subRe, Dn4_3subRRML);

      TComplex sc3232_3subA = correlator.Four_3SubRRML(3, 2, -3, -2);
      double sc3232_3subARe = sc3232_3subA.Re()/Dn4_3subRRML;
      profile.fChsc3232_3subRRMLA->Fill(Ntrks, sc3232_3subARe, Dn4_3subRRML);

      TComplex sc3232_3subB = correlator.Four_3SubRRML(3, 2, -2, -3);
      double sc3232_3subBRe = sc3232_3subB.Re()/Dn4_3subRRML;
      profile.fChsc3232_3subRRMLB->Fill(Ntrks, sc3232_3subBRe, Dn4_3subRRML);

      TComplex sc4242_3subA = correlator.Four_3SubRRML(4, 2, -4, -2);
      double sc4242_3subARe = sc4242_3subA.Re()/Dn4_3subRRML;
      profile.fChsc4242_3subRRMLA->Fill(Ntrks, sc4242_3subARe, Dn4_3subRRML);

      TComplex sc4242_3subB = correlator.Four_3SubRRML(4, 2, -2, -4);
      double sc4242_3subBRe = sc4242_3subB.Re()/Dn4_3subRRML;
      profile.fChsc4242_3subRRMLB->Fill(Ntrks, sc4242_3subBRe, Dn4_3subRRML);
    }
  }

  //..calculate 6-particle correlations and 8-particle corrleations
  //................................


  double Dn6, Dn6Gap0;
  double Dn8, Dn8Gap0;
  Dn6 = correlator.Six(0, 0, 0, 0, 0, 0).Re();
  Dn6Gap0 = correlator.SixGap0(0, 0, 0, 0, 0, 0).Re();
  Dn8 = correlator.Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
  Dn8Gap0 = correlator.EightGap0(0, 0, 0, 0, 0, 0, 0, 0).Re();

  if (fuSixParticleCorrelationStandard) {
    if(NtrksAfter > 5 && Dn6 != 0)
    {

      TComplex v26 = correlator.Six(2, 2, 2, -2, -2, -2);
      double v26Re = v26.Re()/Dn6;
      profile.fChcn6[0]->Fill(Ntrks, v26Re, Dn6);
    }
  }
  if (fuSixParticleCorrelation0Gap) {
    if (NtrksAfterGap0M > 2 && NtrksAfterGap0P > 2 && Dn4Gap0 !=0)
    {
        TComplex v26Gap0 = correlator.SixGap0(2, 2, 2, -2, -2, -2);
        double v26Gap0Re = v26Gap0.Re()/Dn6Gap0;
        profile.fChcn6_Gap0[0]->Fill(Ntrks, v26Gap0Re, Dn6Gap0);
    }
  }
  if (fuEightParticleCorrelationStandard) {
    if(NtrksAfter > 7 && Dn8 != 0)
    {

      TComplex v28 = correlator.Eight(2, 2, 2, 2, -2, -2, -2, -2);
      double v28Re = v28.Re()/Dn8;
      profile.fChcn8[0]->Fill(Ntrks, v28Re, Dn8);
    }
  }
  if (fuEightParticleCorrelation0Gap) {
    if (NtrksAfterGap0M > 3 && NtrksAfterGap0P > 3 && Dn8Gap0 !=0)
    {
        TComplex v28Gap0 = correlator.EightGap0(2, 2, 2, 2, -2, -2, -2, -2);
        double v28Gap0Re = v28Gap0.Re()/Dn8Gap0;
        profile.fChcn8_Gap0[0]->Fill(Ntrks, v28Gap0Re, Dn8Gap0);
    }
  }

}

const char* AliAnalysisTaskXDeptFlow::ReturnPPperiodMC(const Int_t runNumber) const
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

const char* AliAnalysisTaskXDeptFlow::ReturnPPperiod(const Int_t runNumber) const
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
void AliAnalysisTaskXDeptFlow::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate()");
}

int AliAnalysisTaskXDeptFlow::GetEtaPtFlag(double dEta) {
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

//_____________________________________________________________________________
Int_t AliAnalysisTaskXDeptFlow::IdentifyTrack(const AliAODTrack* track) const// called inside prepareTPCTracks() once after IsTrackSelected() to identify Pi, Ka Pr
{
  // checking detector statuses
  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  if(!bIsTPCok) { return -1; }

  Double_t l_Probs[AliPID::kSPECIES];
  double fPIDbayesPion(0.95), fPIDbayesKaon(0.85), fPIDbayesProton(0.85);

  Double_t l_MaxProb[] = {fPIDbayesPion,fPIDbayesKaon,fPIDbayesProton};
  Bool_t l_TOFUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, l_Probs) & AliPIDResponse::kDetTOF;
  Int_t pidInd = 0;
  for(Int_t i(0); i < AliPID::kSPECIES; i++) pidInd=(l_Probs[i]>l_Probs[pidInd])?i:pidInd;
  Int_t retInd = pidInd-AliPID::kPion+1; //realigning
  if(retInd<1 || retInd>3) return -1;
  if(l_Probs[pidInd] < l_MaxProb[retInd-1]) return -1;
  //check nsigma cuts
  if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)pidInd))>3) return -1;
  if(bIsTOFok && l_TOFUsed) if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)pidInd))>3) return -1;

  return retInd;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskXDeptFlow::HasTrackPIDTPC(const AliAODTrack* track) const //called inside IdentifyTrack(), IsK0s() and IsLambda() for PID purposes
{
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  return (pidStatusTPC == AliPIDResponse::kDetPidOk);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskXDeptFlow::HasTrackPIDTOF(const AliAODTrack* track) const //called inside IdentifyTrack() only
{
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  return ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME));
}
