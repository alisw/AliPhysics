///
/// \file V0LamAnalysis/AliAnalysisV0Lam.cxx
///

#include <iostream>
#include <math.h>
#include "TChain.h"
#include "TFile.h"
#include "TCollection.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODv0.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisV0Lam.h"
#include "AliAODMCHeader.h"
#include "AliGenHijingEventHeader.h"

// Analysis task for studying lambda-lambda femtoscopic correlations
// Author: Jai Salzwedel, jai.salzwedel@cern.ch
using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisV0Lam)
/// \endcond

//________________________________________________________________________
AliAnalysisV0Lam::AliAnalysisV0Lam():
  AliAnalysisTaskSE(),
  nEventsToMix(5),
  fAOD(nullptr),
  fOutputList(nullptr),
  fpidAOD(nullptr),
  fCutProcessor(nullptr),
  fMaxV0Mult(700),
  fEtaDaughter(0.8),
  fMassWindowLam(0.00568),
  fTOFLow(0.8),

  fSigmaCutTOFPion(4.0),
  fSigmaCutTPCPion(3.0),
  fSigmaCutTOFProton(4.0),
  fSigmaCutTPCProton(3.0),

  fPDGLambda(1.115683),
  fPDGProton(.938272),
  fPDGPion(.1395702),

  fEventCount(0),
  fSysStudyType(kNoStudy),

  fNumberOfTopologicalCutValues(1),
  fNumberOfCfVariableCutValues(1),
  fVariableCutType(0),
  fNominalTopCutIndex(0),

  fTestNoTTC(kFALSE),
  fNumberVariableAvgSepCuts(3),
  fIsUsingVariableAvgSepCut(kFALSE),
  fIsMCEvent(kFALSE),
  fFlattenCent(kTRUE),
  fEC(NULL),
  fEvt(NULL),

  fTotalLambda(0),
  fTotalAntiLambda(0),
  fV0Candidates(0),

  fTPCVsPPosLam(NULL), fTPCVsPNegLam(NULL),
  fTPCVsPPosALam(NULL), fTPCVsPNegALam(NULL),
  fMultDistRough(NULL),
  fMultDistLambda(NULL), fMultDistAntiLambda(NULL),
  fMultCentLambda(NULL), fMultCentAntiLambda(NULL),
  fCentrality(NULL),
  fRemainingFromBeginningToV0Finder(NULL),
  fRemainingFromBeginningToRecon(NULL),
  fRemainingFromV0FinderToRecon(NULL),
  fMCTruthOfOriginalParticles(NULL),
  fMCTruthOfV0FinderParticles(NULL),
  fMCTruthOfReconstructedParticles(NULL),
  fMCFakeParticleIdentity(NULL),
  fMCOtherV0Identity(NULL),
// fDataCompeted(NULL), fDataCulled(NULL),
// fRemainingV0s(NULL), fRemainingFrac(NULL),
  fResMatrixLLSameAll(NULL),
  fResMatrixAASameAll(NULL),
  fResMatrixLASameAll(NULL),
  fResMatrixLLMixedAll(NULL),
  fResMatrixAAMixedAll(NULL),
  fResMatrixLAMixedAll(NULL),
  fResMatrixLLSamePure(NULL),
  fResMatrixAASamePure(NULL),
  fResMatrixLASamePure(NULL),
  fResMatrixLLMixedPure(NULL),
  fResMatrixAAMixedPure(NULL),
  fResMatrixLAMixedPure(NULL),
  fMCTruthPtLam(NULL), fMCTruthPtALam(NULL),
  fMCTruthPhiLam(NULL), fMCTruthPhiALam(NULL),
  fMCTruthEtaLam(NULL), fMCTruthEtaALam(NULL),
  fSignalLamLam(NULL), fBkgLamLam(NULL),
  fSignalALamALam(NULL), fBkgALamALam(NULL),
  fSignalLamALam(NULL), fBkgLamALam(NULL),
  fSignalKtVsKstarLamLam(NULL), fBkgKtVsKstarLamLam(NULL),
  fSignalKtVsKstarALamALam(NULL), fBkgKtVsKstarALamALam(NULL),
  fSignalKtVsKstarLamALam(NULL), fBkgKtVsKstarLamALam(NULL),
  fHistSignalProperDecayLengthDiffLamLam(NULL),
  fHistSignalProperDecayLengthDiffALamALam(NULL),
  fHistSignalProperDecayLengthDiffLamALam(NULL),
  fHistBkgProperDecayLengthDiffLamLam(NULL),
  fHistBkgProperDecayLengthDiffALamALam(NULL),
  fHistBkgProperDecayLengthDiffLamALam(NULL),
  fSignalLamLamProtSep(NULL), fSignalLamLamPiMinusSep(NULL),
  fSignalALamALamAntiProtSep(NULL), fSignalALamALamPiPlusSep(NULL),
  fSignalLamALamProtPiPlusSep(NULL), fSignalLamALamAntiProtPiMinusSep(NULL),
  fBkgLamLamProtSep(NULL), fBkgLamLamPiMinusSep(NULL),
  fBkgALamALamAntiProtSep(NULL), fBkgALamALamPiPlusSep(NULL),
  fBkgLamALamProtPiPlusSep(NULL), fBkgLamALamAntiProtPiMinusSep(NULL),
  fSignalLamLamPlusMinusSep(NULL), fSignalALamALamPlusMinusSep(NULL),
  fSignalLamALamProtSep(NULL), fSignalLamALamPionSep(NULL),
  fBkgLamLamPlusMinusSep(NULL), fBkgALamALamPlusMinusSep(NULL),
  fBkgLamALamProtSep(NULL), fBkgLamALamPionSep(NULL)
{
}
//________________________________________________________________________

AliAnalysisV0Lam::AliAnalysisV0Lam(const char *name, SysStudy sysStudyType, Int_t varCutType, Bool_t flattenCent, Int_t nMixingEvents, Bool_t testNoTwoTrackCuts):
  AliAnalysisTaskSE(name),
  nEventsToMix(nMixingEvents),
  fAOD(nullptr),
  fOutputList(nullptr),
  fpidAOD(nullptr),
  fCutProcessor(nullptr),
  fMaxV0Mult(700),
  fEtaDaughter(0.8),
  fMassWindowLam(0.00568),
  fTOFLow(0.8),

  fSigmaCutTOFPion(4.0),
  fSigmaCutTPCPion(3.0),
  fSigmaCutTOFProton(4.0),
  fSigmaCutTPCProton(3.0),

  fPDGLambda(1.115683),
  fPDGProton(.938272),
  fPDGPion(.1395702),

  fEventCount(0),
  fSysStudyType(sysStudyType),
  fNumberOfTopologicalCutValues(1),
  fNumberOfCfVariableCutValues(1),

  fVariableCutType(varCutType),
  fNominalTopCutIndex(0),

  fTestNoTTC(testNoTwoTrackCuts),
  fNumberVariableAvgSepCuts(3),
  fIsUsingVariableAvgSepCut(kFALSE),
  fIsMCEvent(kFALSE),
  fFlattenCent(flattenCent),

  fEC(NULL),
  fEvt(NULL),

  fTotalLambda(0),
  fTotalAntiLambda(0),
  fV0Candidates(0),

  fTPCVsPPosLam(NULL), fTPCVsPNegLam(NULL),
  fTPCVsPPosALam(NULL), fTPCVsPNegALam(NULL),
  fMultDistRough(NULL),
  fMultDistLambda(NULL), fMultDistAntiLambda(NULL),
  fMultCentLambda(NULL), fMultCentAntiLambda(NULL),
  fCentrality(NULL),
  fRemainingFromBeginningToV0Finder(NULL),
  fRemainingFromBeginningToRecon(NULL),
  fRemainingFromV0FinderToRecon(NULL),
  fMCTruthOfOriginalParticles(NULL),
  fMCTruthOfV0FinderParticles(NULL),
  fMCTruthOfReconstructedParticles(NULL),
  fMCFakeParticleIdentity(NULL),
  fMCOtherV0Identity(NULL),
  // fDataCompeted(NULL), fDataCulled(NULL),
  // fRemainingV0s(NULL), fRemainingFrac(NULL),
  fResMatrixLLSameAll(NULL),
  fResMatrixAASameAll(NULL),
  fResMatrixLASameAll(NULL),
  fResMatrixLLMixedAll(NULL),
  fResMatrixAAMixedAll(NULL),
  fResMatrixLAMixedAll(NULL),
  fResMatrixLLSamePure(NULL),
  fResMatrixAASamePure(NULL),
  fResMatrixLASamePure(NULL),
  fResMatrixLLMixedPure(NULL),
  fResMatrixAAMixedPure(NULL),
  fResMatrixLAMixedPure(NULL),
  fMCTruthPtLam(NULL), fMCTruthPtALam(NULL),
  fMCTruthPhiLam(NULL), fMCTruthPhiALam(NULL),
  fMCTruthEtaLam(NULL), fMCTruthEtaALam(NULL),
  fSignalLamLam(NULL), fBkgLamLam(NULL),
  fSignalALamALam(NULL), fBkgALamALam(NULL),
  fSignalLamALam(NULL), fBkgLamALam(NULL),
  fSignalKtVsKstarLamLam(NULL), fBkgKtVsKstarLamLam(NULL),
  fSignalKtVsKstarALamALam(NULL), fBkgKtVsKstarALamALam(NULL),
  fSignalKtVsKstarLamALam(NULL), fBkgKtVsKstarLamALam(NULL),
  fHistSignalProperDecayLengthDiffLamLam(NULL),
  fHistSignalProperDecayLengthDiffALamALam(NULL),
  fHistSignalProperDecayLengthDiffLamALam(NULL),
  fHistBkgProperDecayLengthDiffLamLam(NULL),
  fHistBkgProperDecayLengthDiffALamALam(NULL),
  fHistBkgProperDecayLengthDiffLamALam(NULL),
  fSignalLamLamProtSep(NULL), fSignalLamLamPiMinusSep(NULL),
  fSignalALamALamAntiProtSep(NULL), fSignalALamALamPiPlusSep(NULL),
  fSignalLamALamProtPiPlusSep(NULL), fSignalLamALamAntiProtPiMinusSep(NULL),
  fBkgLamLamProtSep(NULL), fBkgLamLamPiMinusSep(NULL),
  fBkgALamALamAntiProtSep(NULL), fBkgALamALamPiPlusSep(NULL),
  fBkgLamALamProtPiPlusSep(NULL), fBkgLamALamAntiProtPiMinusSep(NULL),
  fSignalLamLamPlusMinusSep(NULL), fSignalALamALamPlusMinusSep(NULL),
  fSignalLamALamProtSep(NULL), fSignalLamALamPionSep(NULL),
  fBkgLamLamPlusMinusSep(NULL), fBkgALamALamPlusMinusSep(NULL),
  fBkgLamALamProtSep(NULL), fBkgLamALamPionSep(NULL)
{
  // Define output slots here
  // Output slot #1
  DefineOutput(1, TList::Class());
  if (kTopologicalStudy == fSysStudyType) {
    fNominalTopCutIndex = 1;
  }
}
//________________________________________________________________________
AliAnalysisV0Lam::~AliAnalysisV0Lam()
{
  // Destructor

  for(unsigned short i=0; i<zVertexBins; i++)
  {
    for(unsigned short j=0; j<nCentBins; j++)
    {
      delete fEC[i][j];
    }
    delete[] fEC[i];
  }
  delete[] fEC;
  delete fCutProcessor;
  if(fOutputList) delete fOutputList; //This cleans up all output hists
}
//________________________________________________________________________
void AliAnalysisV0Lam::MyInit()
{
  // Setup V0 cut processor
  AliAnalysisV0LamCutProcessing::CutType_t variedTopologicalCut = AliAnalysisV0LamCutProcessing::kNoCut;
  if(kTopologicalStudy == fSysStudyType) {
    variedTopologicalCut = (AliAnalysisV0LamCutProcessing::CutType_t)fVariableCutType;
  }
  fCutProcessor = new AliAnalysisV0LamCutProcessing(fOutputList, variedTopologicalCut);
  if (kTopologicalStudy == fSysStudyType) {
    fNumberOfTopologicalCutValues = fCutProcessor->GetNumberOfVariableCutValues();
    fNumberOfCfVariableCutValues = fNumberOfTopologicalCutValues;
  } else if (kTwoTrackStudy == fSysStudyType) {
    fNumberOfTopologicalCutValues = 1;
    fNumberOfCfVariableCutValues = 3;
  } else {
    fNumberOfTopologicalCutValues = 1;
    fNumberOfCfVariableCutValues = 1;
  }
  // cout<<"Number of variable cf cut values: "<<fNumberOfCfVariableCutValues<<endl;

  //setup event collection for event mixing
  fEC = new AliAnalysisV0LamEventCollection **[zVertexBins];
  for(unsigned short i=0; i<zVertexBins; i++)
  {
    fEC[i] = new AliAnalysisV0LamEventCollection *[nCentBins];
    for(unsigned short j=0; j<nCentBins; j++)
    {
      fEC[i][j] = new AliAnalysisV0LamEventCollection(nEventsToMix+1, fMaxV0Mult);
    }
  }
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fpidAOD = aodH->GetAODpidUtil();
}
//________________________________________________________________________
void AliAnalysisV0Lam::UserCreateOutputObjects()
{
  // Create output histograms.
  // Histograms are added to fOutputList
  // When fOutputList is deleted, it automatically cleans up all
  // associated histograms
  // cout<<"Create histograms"<<endl;
  fOutputList = new TList();
  fOutputList->SetOwner();
  MyInit();// Initialize my settings

  fMultDistRough = new TH1F("fMultDistRough","Multiplicity Distribution",301,-.5,3001-.5);
  fMultDistRough->GetXaxis()->SetTitle("Event Multiplicity (pions)");
  fMultDistRough->GetYaxis()->SetTitle("# of events");
  fOutputList->Add(fMultDistRough);

  fCentrality = new TH1F("fCentrality", "Centrality Percentage of Event;Centrality %;# of Events", 100, 0., 100.);
  fOutputList->Add(fCentrality);

  //TPC signal vs track momentum
  fTPCVsPPosLam = new TH2F("fTPCVsPPosLam","TPC dE/dx Lambda Protons;p (GeV/c);TPC Signal",50,0.,2.,250,0.,500.);
  fOutputList->Add(fTPCVsPPosLam);
    //TPC signal vs track momentum
  fTPCVsPNegLam = new TH2F("fTPCVsPNegLam","TPC dE/dx Lambda Pions;p (GeV/c);TPC Signal",50,0.,2.,250,0.,500.);
  fOutputList->Add(fTPCVsPNegLam);
    //TPC signal vs track momentum
  fTPCVsPPosALam = new TH2F("fTPCVsPPosALam","TPC dE/dx Antilambda Pions;p (GeV/c);TPC Signal",50,0.,2.,250,0.,500.);
  fOutputList->Add(fTPCVsPPosALam);
    //TPC signal vs track momentum
  fTPCVsPNegALam = new TH2F("fTPCVsPNegALam","TPC dE/dx Antilambda Protons;p (GeV/c);TPC Signal",50,0.,2.,250,0.,500.);
  fOutputList->Add(fTPCVsPNegALam);

  //V0 Shared daughter culling statistics
  // fDataCompeted = new TH1F("fDataCompeted","Thunderdome Particles Competing;Particles Competing;# of Events",26, -0.5, 25.5);
  // fOutputList->Add(fDataCompeted);

  // fDataCulled = new TH1F("fDataCulled","Thunderdome Particles Removed;Particles Removed;# of Events",26, -0.5, 25.5);
  // fOutputList->Add(fDataCulled);

  // // fRemainingV0s = new TH1F("fRemainingV0s","",26, -0.5, 25.5);
  // fOutputList->Add(fRemainingV0s);

  // fRemainingFrac = new TH1F("fRemainingFrac","",101, -.005, 1.005);
  // fOutputList->Add(fRemainingFrac);

  fRemainingFromBeginningToV0Finder = new TH2F("fRemainingFromBeginningToV0Finder","Fraction Remaining At V0Finder Stage", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, 21, -0.025, 1.025);
  fOutputList->Add(fRemainingFromBeginningToV0Finder);

  fMCTruthOfOriginalParticles = new TH1F("fMCTruthOfOriginalParticles","MC Truth of Original Particles in Event", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1);
  fOutputList->Add(fMCTruthOfOriginalParticles);

  fMCTruthOfV0FinderParticles = new TH1F("fMCTruthOfV0FinderParticles","MC Truth of V0Finder Particles in Event", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1);
  fOutputList->Add(fMCTruthOfV0FinderParticles);

  SetBinsOnOriginHists(fRemainingFromBeginningToV0Finder);
  SetBinsOnOriginHists(fMCTruthOfOriginalParticles);
  SetBinsOnOriginHists(fMCTruthOfV0FinderParticles);
  //The first dimension is the index of the variable cut value
  fRemainingFromBeginningToRecon = new TH3F("fRemainingFromBeginningToRecon", "Fraction Remaining After Reconstruction", fNumberOfTopologicalCutValues, -0.5, fNumberOfTopologicalCutValues -0.5,  AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, 21, -0.025, 1.025);
  fRemainingFromV0FinderToRecon = new TH3F("fRemainingFromV0FinderToRecon", "Fraction From V0Finder That Remain After Reconstruction Stage", fNumberOfTopologicalCutValues, -0.5, fNumberOfTopologicalCutValues -0.5,  AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, 21, -0.025, 1.025);
  fMCTruthOfReconstructedParticles = new TH2F("fMCTruthOfReconstructedParticles", "MC Truth of Reconstructed Particles in Event", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, fNumberOfTopologicalCutValues, -0.5, fNumberOfTopologicalCutValues -0.5);
  SetBinsOnOriginHists(fRemainingFromBeginningToRecon);
  SetBinsOnOriginHists(fRemainingFromV0FinderToRecon);
  SetBinsOnOriginHists(fMCTruthOfReconstructedParticles);

  // Particle multiplicities
  fMultDistLambda = new TH2F("fMultDistLambda", "Lambda multiplicity;Cut Bin;Lambda Found;# of Events", fNumberOfTopologicalCutValues, -0.5, fNumberOfTopologicalCutValues -0.5,  21, -0.5, 21-0.5);
  fMultDistAntiLambda = new TH2F("fMultDistAntiLambda", "AntiLambda multiplicity;Cut Bin;Antilambda Found;# of Events", fNumberOfTopologicalCutValues, -0.5, fNumberOfTopologicalCutValues -0.5,  21, -0.5, 21-0.5);
  fMultCentLambda = new TH3F("fMultCentLambda", "Lambda multiplicity vs centrality;Cut Bin; Centrality Bin;Lambdas Found", fNumberOfTopologicalCutValues, -0.5, fNumberOfTopologicalCutValues -0.5,  nCentBins, .5, nCentBins+.5, 21, -0.5, 21-0.5);
  fMultCentAntiLambda =  new TH3F("fMultCentAntiLambda", "AntiLambda multiplicity vs centrality;Cut Bin; Centrality Bin;Antilambdas Found", fNumberOfTopologicalCutValues, -0.5, fNumberOfTopologicalCutValues -0.5,  nCentBins, .5, nCentBins+.5, 21, -0.5, 21-0.5);

  // fMultDistLambda->GetXaxis()->SetTitle("Var Bin");
  // fMultDistLambda->GetYaxis()->SetTitle("Event Multiplicity (Lambdas)");
  // fMultDistAntiLambda->GetXaxis()->SetTitle("Var Bin");
  // fMultDistAntiLambda->GetXaxis()->SetTitle("Event Multiplicity (AntiLambdas)");
  // fMultCentLambda->GetXaxis()->SetTitle("Var Bin");
  // fMultCentLambda->GetYaxis()->SetTitle("Centrality");
  // fMultCentLambda->GetZaxis()->SetTitle("Event Multiplicity (Lambdas)");
  // fMultCentAntiLambda->GetXaxis()->SetTitle("Var Bin");
  // fMultCentAntiLambda->GetYaxis()->SetTitle("Centrality");
  // fMultCentAntiLambda->GetZaxis()->SetTitle("Event Multiplicity (AntiLambdas)");
  fOutputList->Add(fRemainingFromBeginningToRecon);
  fOutputList->Add(fRemainingFromV0FinderToRecon);
  fOutputList->Add(fMCTruthOfReconstructedParticles);
  fOutputList->Add(fMultDistLambda);
  fOutputList->Add(fMultDistAntiLambda);
  fOutputList->Add(fMultCentLambda);
  fOutputList->Add(fMultCentAntiLambda);


  fMCFakeParticleIdentity = new TH1F("fMCFakeParticleIdentity", "Breakdown of fake particles", 3, 0, 3);
  fOutputList->Add(fMCFakeParticleIdentity);
  fMCFakeParticleIdentity->GetXaxis()->SetBinLabel(1,"Fake Lambda");
  fMCFakeParticleIdentity->GetXaxis()->SetBinLabel(2,"Fake AntiLambda");
  fMCFakeParticleIdentity->GetXaxis()->SetBinLabel(3,"Fake K0Short");

  fMCOtherV0Identity = new TH1F("fMCOtherV0Identity", "Breakdown of otherV0s particles", 3, 0, 3);
  fOutputList->Add(fMCOtherV0Identity);
  fMCOtherV0Identity->GetXaxis()->SetBinLabel(1,"Fake Lambda");
  fMCOtherV0Identity->GetXaxis()->SetBinLabel(2,"Fake AntiLambda");
  fMCOtherV0Identity->GetXaxis()->SetBinLabel(3,"Fake K0Short");
  int kTBins = 20;
  double maxKtBin = 4.;
  int kStarBins = 800;
  double maxKStar = 2.;
  //Pair kT Tracking: kT bins, centrality bins, k* bins
  fSignalKtVsKstarLamLam = new TH3F("fSignalKtVsKstarLamLam", "LamLam Pair Kt Same Event;k_T;CentBin;k*", kTBins, 0., maxKtBin, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalKtVsKstarLamLam);
  fSignalKtVsKstarALamALam = new TH3F("fSignalKtVsKstarALamALam", "ALamALam Pair Kt Same Event;k_T;CentBin;k*", kTBins, 0., maxKtBin, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalKtVsKstarALamALam);
  fSignalKtVsKstarLamALam = new TH3F("fSignalKtVsKstarLamALam", "LamALam Pair Kt Same Event;k_T;CentBin;k*", kTBins, 0., maxKtBin, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalKtVsKstarLamALam);

  fBkgKtVsKstarLamLam = new TH3F("fBkgKtVsKstarLamLam", "LamLam Pair Kt Mixed Event;k_T;CentBin;k*", kTBins, 0., maxKtBin, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgKtVsKstarLamLam);
  fBkgKtVsKstarALamALam = new TH3F("fBkgKtVsKstarALamALam", "ALamALam Pair Kt Mixed Event;k_T;CentBin;k*", kTBins, 0., maxKtBin, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgKtVsKstarALamALam);
  fBkgKtVsKstarLamALam = new TH3F("fBkgKtVsKstarLamALam", "LamALam Pair Kt Mixed Event;k_T;CentBin;k*", kTBins, 0., maxKtBin, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgKtVsKstarLamALam);





  /////////Signal Distributions///////////////////
  //First bin is variable cut value, second bin is centrality, third bin is Kstar
  fSignalLamLam = new TH3F("fSignalLamLam","Same Event Pair Distribution;VarBin;CentBin;k*", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamLam);

  fBkgLamLam = new TH3F("fBkgLamLam","Mixed Event Pair Distribution;VarBin;CentBin;k*", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamLam);

  fSignalALamALam = new TH3F("fSignalALamALam","Same Event Pair Distribution;VarBin;CentBin;k*", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalALamALam);

  fBkgALamALam = new TH3F("fBkgALamALam","Mixed Event Pair Distribution;VarBin;CentBin;k*", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgALamALam);

  fSignalLamALam = new TH3F("fSignalLamALam","Same Event Pair Distribution;VarBin;CentBin;k*", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALam);

  fBkgLamALam = new TH3F("fBkgLamALam","Mixed Event Pair Distribution;VarBin;CentBin;k*", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALam);


  // Distributions of difference in proper decaylength
  TObjArray *dlDiffDir = new TObjArray();
  dlDiffDir->SetName("ProperDecayLengthDiff");
  fOutputList->Add(dlDiffDir);

  Int_t dlDiffBins = 500;
  Double_t dlDiffMaxValue = 50.;
  fHistSignalProperDecayLengthDiffLamLam = new TH1F ("fHistSignalProperDecayLengthDiffLamLam", "Difference in proper decay length of #Lambda#Lambda pairs;ProperDecayLengthDiff", dlDiffBins, 0., dlDiffMaxValue);
  dlDiffDir->Add(fHistSignalProperDecayLengthDiffLamLam);
  fHistSignalProperDecayLengthDiffALamALam = new TH1F ("fHistSignalProperDecayLengthDiffALamALam", "Difference in proper decay length of #bar{#Lambda}#bar{#Lambda} pairs;ProperDecayLengthDiff", dlDiffBins, 0., dlDiffMaxValue);
  dlDiffDir->Add(fHistSignalProperDecayLengthDiffALamALam);
  fHistSignalProperDecayLengthDiffLamALam = new TH1F ("fHistSignalProperDecayLengthDiffLamALam", "Difference in proper decay length of #Lambda#bar{#Lambda} pairs;ProperDecayLengthDiff", dlDiffBins, 0., dlDiffMaxValue);
  dlDiffDir->Add(fHistSignalProperDecayLengthDiffLamALam);

  fHistBkgProperDecayLengthDiffLamLam = new TH1F ("fHistBkgProperDecayLengthDiffLamLam", "Difference in proper decay length of #Lambda#Lambda pairs;ProperDecayLengthDiff", dlDiffBins, 0., dlDiffMaxValue);
  dlDiffDir->Add(fHistBkgProperDecayLengthDiffLamLam);
  fHistBkgProperDecayLengthDiffALamALam = new TH1F ("fHistBkgProperDecayLengthDiffALamALam", "Difference in proper decay length of #bar{#Lambda}#bar{#Lambda} pairs;ProperDecayLengthDiff", dlDiffBins, 0., dlDiffMaxValue);
  dlDiffDir->Add(fHistBkgProperDecayLengthDiffALamALam);
  fHistBkgProperDecayLengthDiffLamALam = new TH1F ("fHistBkgProperDecayLengthDiffLamALam", "Difference in proper decay length of #Lambda#bar{#Lambda} pairs;ProperDecayLengthDiff", dlDiffBins, 0., dlDiffMaxValue);
  dlDiffDir->Add(fHistBkgProperDecayLengthDiffLamALam);
  // Daughter pair separation distributions
  // Binned according to (average separation, pT1, pT2)
  TObjArray *sepDirNew = new TObjArray();
  sepDirNew->SetName("AvgSepNew");
  fOutputList->Add(sepDirNew);


  Int_t avgSepBins = 400;
  Double_t avgSepMaxValue = 40.;
  Int_t ptBins = 30;
  Double_t ptMax = 3.;
  fSignalLamLamProtSep = new TH3F ("fSignalLamLamProtSep","Proton pair sep for Lam-Lam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalLamLamProtSep);

  fSignalLamLamPiMinusSep = new TH3F ("fSignalLamLamPiMinusSep","PiMinus pair sep for Lam-Lam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalLamLamPiMinusSep);

  fSignalALamALamAntiProtSep = new TH3F ("fSignalALamALamAntiProtSep","AntiProton pair sep for ALam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalALamALamAntiProtSep);

  fSignalALamALamPiPlusSep = new TH3F ("fSignalALamALamPiPlusSep","PiPlus pair sep for ALam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalALamALamPiPlusSep);

  fSignalLamALamAntiProtPiMinusSep = new TH3F ("fSignalLamALamAntiProtPiMinusSep","Neg particle pair sep for Lam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalLamALamAntiProtPiMinusSep);

  fSignalLamALamProtPiPlusSep = new TH3F ("fSignalLamALamProtPiPlusSep","Pos pair sep for Lam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalLamALamProtPiPlusSep);

  fBkgLamLamProtSep = new TH3F ("fBkgLamLamProtSep","Proton pair sep for Lam-Lam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgLamLamProtSep);

  fBkgLamLamPiMinusSep = new TH3F ("fBkgLamLamPiMinusSep","PiMinus pair sep for Lam-Lam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgLamLamPiMinusSep);

  fBkgALamALamAntiProtSep = new TH3F ("fBkgALamALamAntiProtSep","AntiProton pair sep for ALam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgALamALamAntiProtSep);

  fBkgALamALamPiPlusSep = new TH3F ("fBkgALamALamPiPlusSep","PiPlus pair sep for ALam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgALamALamPiPlusSep);

  fBkgLamALamAntiProtPiMinusSep = new TH3F ("fBkgLamALamAntiProtPiMinusSep","Neg particle pair sep for Lam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgLamALamAntiProtPiMinusSep);

  fBkgLamALamProtPiPlusSep = new TH3F ("fBkgLamALamProtPiPlusSep","Pos pair sep for Lam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgLamALamProtPiPlusSep);



  //opposite charged pair separation
  fSignalLamLamPlusMinusSep = new TH3F ("fSignalLamLamPlusMinusSep","Proton Pion pair sep for Lam-Lam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalLamLamPlusMinusSep);

  fSignalALamALamPlusMinusSep = new TH3F ("fSignalALamALamPlusMinusSep","Proton Pion pair sep for ALam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalALamALamPlusMinusSep);

  fSignalLamALamProtSep = new TH3F ("fSignalLamALamProtSep","Proton pair sep for Lam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalLamALamProtSep);

  fSignalLamALamPionSep = new TH3F ("fSignalLamALamPionSep","Pion pair sep for Lam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fSignalLamALamPionSep);


  fBkgLamLamPlusMinusSep = new TH3F ("fBkgLamLamPlusMinusSep","Proton Pion pair sep for Lam-Lam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgLamLamPlusMinusSep);

  fBkgALamALamPlusMinusSep = new TH3F ("fBkgALamALamPlusMinusSep","Proton Pion pair sep for ALam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgALamALamPlusMinusSep);

  fBkgLamALamProtSep = new TH3F ("fBkgLamALamProtSep","Proton pair sep for Lam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgLamALamProtSep);

  fBkgLamALamPionSep = new TH3F ("fBkgLamALamPionSep","Pion pair sep for Lam-ALam;AvgSep;pT1;pT2", avgSepBins, 0., avgSepMaxValue, ptBins, 0., ptMax, ptBins, 0., ptMax);
  sepDirNew->Add(fBkgLamALamPionSep);


  // Directory for MC Truth info
  TObjArray *arrMC = new TObjArray();
  arrMC->SetName("MCTruthInfo");
  fOutputList->Add(arrMC);

  fMCTruthPtLam = new TH1F("fMCTruthPtLam","p_{T} Lambda;p_{T};counts", 500, 0., 10.);
  arrMC->Add(fMCTruthPtLam);
  fMCTruthPtALam = new TH1F("fMCTruthPtALam","p_{T} AntiLambda;p_{T};counts", 500, 0., 10.);
  arrMC->Add(fMCTruthPtALam);

  fMCTruthPhiLam = new TH1F("fMCTruthPhiLam","Phi Lambda;Phi;counts", 200, -1.*TMath::Pi(), TMath::Pi());
  arrMC->Add(fMCTruthPhiLam);
  fMCTruthPhiALam = new TH1F("fMCTruthPhiALam","Phi AntiLambda;Phi;counts",  200, -1.*TMath::Pi(), TMath::Pi());
  arrMC->Add(fMCTruthPhiLam);

  fMCTruthEtaLam = new TH1F("fMCTruthEtaLam","Eta Lambda;Eta;counts", 200, -1., 1.);
  arrMC->Add(fMCTruthEtaLam);
  fMCTruthEtaALam = new TH1F("fMCTruthEtaALam","Eta AntiLambda;Eta;counts", 200, -1., 1.);
  arrMC->Add(fMCTruthEtaALam);

  // Directory for momentum resolution correction histograms
  TObjArray *resArr = new TObjArray();
  resArr->SetName("ResolutionMatrices");
  fOutputList->Add(resArr);
  fResMatrixLLSameAll = new TH2F("fResMatrixLLSameAll","ResolutionMatrix LLSameAll;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixAASameAll = new TH2F("fResMatrixAASameAll","ResolutionMatrix AASameAll;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixLASameAll = new TH2F("fResMatrixLASameAll","ResolutionMatrix LASameAll;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixLLMixedAll = new TH2F("fResMatrixLLMixedAll","ResolutionMatrix LLMixedAll;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixAAMixedAll = new TH2F("fResMatrixAAMixedAll","ResolutionMatrix AAMixedAll;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixLAMixedAll = new TH2F("fResMatrixLAMixedAll","ResolutionMatrix LAMixedAll;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixLLSamePure = new TH2F("fResMatrixLLSamePure","ResolutionMatrix LLSamePure;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixAASamePure = new TH2F("fResMatrixAASamePure","ResolutionMatrix AASamePure;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixLASamePure = new TH2F("fResMatrixLASamePure","ResolutionMatrix LASamePure;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixLLMixedPure = new TH2F("fResMatrixLLMixedPure","ResolutionMatrix LLMixedPure;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixAAMixedPure = new TH2F("fResMatrixAAMixedPure","ResolutionMatrix AAMixedPure;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fResMatrixLAMixedPure = new TH2F("fResMatrixLAMixedPure","ResolutionMatrix LAMixedPure;k^{*}_{true};k^{*}_{rec}", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  resArr->Add(fResMatrixLLSameAll);
  resArr->Add(fResMatrixAASameAll);
  resArr->Add(fResMatrixLASameAll);
  resArr->Add(fResMatrixLLMixedAll);
  resArr->Add(fResMatrixAAMixedAll);
  resArr->Add(fResMatrixLAMixedAll);
  resArr->Add(fResMatrixLLSamePure);
  resArr->Add(fResMatrixAASamePure);
  resArr->Add(fResMatrixLASamePure);
  resArr->Add(fResMatrixLLMixedPure);
  resArr->Add(fResMatrixAAMixedPure);
  resArr->Add(fResMatrixLAMixedPure);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisV0Lam::Exec(Option_t *)
{
  // Main loop
  // Called for each event
  // cout<<"Exec"<<endl;
  //Make sure we are using the correct event triggers
  if(!IsCorrectEventTrigger()) return;
  //  cout<<"===========  Event # "<<fEventCount+1<<"  ==========="<<endl;
  fEventCount++;
  fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!fAOD) {Printf("ERROR: fAOD not available"); return;}

  //Get Monte Carlo data if available
  TClonesArray *mcArray = 0x0;
  mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
  //Check for MC event headers and get # of original Hijing particles
  //This is used later for rejecting injected signals
  Int_t numberOfLastHijingLabel = 0;
  AliAODMCHeader *mcHeader = (AliAODMCHeader*)fAOD->FindListObject(AliAODMCHeader::StdBranchName());
  if(mcHeader){
    // Get the iterator on the list of cocktail headers
    TIter next(mcHeader->GetCocktailHeaders());
    // Loop over the cocktail headers
    while (const TObject *obj=next()){
      // Check whether it's a Hijing header
      const AliGenHijingEventHeader* hijingHeader = dynamic_cast<const AliGenHijingEventHeader*>(obj);
      if(hijingHeader) {
	numberOfLastHijingLabel=hijingHeader->NProduced()-1;
	// We're done!
	//printf("Good! Last Hijing particle label %d\n",numberOfLastHijingLabel);
      } // End of found the hijing header
    } // End of loop over cocktail headers
  } // End of MC header exists
  else if(mcArray) cerr<<"Could not find mcHeader!"<<endl;
  if(mcArray) fIsMCEvent = kTRUE;
  else fIsMCEvent = kFALSE;

  //Centrality selection
  AliCentrality *centrality = fAOD->GetCentrality();
  float centralityPercentile = centrality->GetCentralityPercentile("V0M");


  // flatten centrality dist.
  if(fFlattenCent && centralityPercentile < 9){
    if(RejectEventCentFlat(fAOD->GetMagneticField(),centralityPercentile)) {
      // cout<<"\t\t\tSkipping Event: Centrality flattening\n\n";
      return;}
  }

  int centralityBin=0;
  //Printf("Centrality percent = %f", centralityPercentile);
  fCentrality->Fill(centralityPercentile);
  AliAODVZERO *aodV0 = fAOD->GetVZEROData();
  Float_t multV0A=aodV0->GetMTotV0A();
  Float_t multV0C=aodV0->GetMTotV0C();
  if(centralityPercentile < 0) {
    //Printf("No centrality info");
    return;}
  else if(centralityPercentile == 0 && (multV0A + multV0C < 19500)) {
    //Printf("No centrality info");
    return;}
  // Centrality info is good.  Now find the correct 5% centrality bin
  else if(centralityPercentile <= 5.) centralityBin=0;
  else if(centralityPercentile <= 10.) centralityBin=1;
  else if(centralityPercentile <= 15.) centralityBin=2;
  else if(centralityPercentile <= 20.) centralityBin=3;
  else if(centralityPercentile <= 25.) centralityBin=4;
  else if(centralityPercentile <= 30.) centralityBin=5;
  else if(centralityPercentile <= 35.) centralityBin=6;
  else if(centralityPercentile <= 40.) centralityBin=7;
  else if(centralityPercentile <= 45.) centralityBin=8;
  else if(centralityPercentile <= 50.) centralityBin=9;
  else if(centralityPercentile <= 55.) centralityBin=10;
  else if(centralityPercentile <= 60.) centralityBin=11;
  else if(centralityPercentile <= 65.) centralityBin=12;
  else if(centralityPercentile <= 70.) centralityBin=13;
  else if(centralityPercentile <= 75.) centralityBin=14;
  else if(centralityPercentile <= 80.) centralityBin=15;
  else if(centralityPercentile <= 85.) centralityBin=16;
  else if(centralityPercentile <= 90.) centralityBin=17;
  else if(centralityPercentile <= 95.) centralityBin=18;
  else if(centralityPercentile <= 100.) centralityBin=19;
  else {/*Printf("Skipping Peripheral Event");*/ return;}

  //Vertexing
  AliAODVertex *primaryVertexAOD = fAOD->GetPrimaryVertex();
  TVector3 vertex(primaryVertexAOD->GetX(),
		  primaryVertexAOD->GetY(),
		  primaryVertexAOD->GetZ());

  double zStep=20./double(zVertexBins), zStart=-10.;
  if(vertex.x()<10e-5 && vertex.y()<10e-5 &&  vertex.z()<10e-5) {
    // cout<<"\t\t\tSkipping Event: Vertex at 0,0,0\n\n";
    return;
  }
  if(fabs(vertex.z()) > fabs(zStart)) {
    // cout<<"\t\t\tSkipping Event: Outside ZVertex range\n\n";
    return; // Z-Vertex Cut
  }

  if(!primaryVertexAOD || primaryVertexAOD->GetNContributors() < 1){
    // cout<<"\t\t\tSkipping Event: No vertex\n\n";
   return;
  }

  int zBin=0;
  for(int i=0; i<zVertexBins; i++)
  {
    if((vertex.z() > zStart+i*zStep) && (vertex.z() < zStart+(i+1)*zStep))
    {
      zBin=i;
      break;
    }
  }
  double bfield = fAOD->GetMagneticField();
  //Printf("Rough multiplicity = %d", fAOD->GetNumberOfTracks());
  fMultDistRough->Fill(fAOD->GetNumberOfTracks());
  /////////////////////////////////////////////////////////
  //Add Event to buffer - this is for event mixing
  fEC[zBin][centralityBin]->FifoShift();
  fEvt = fEC[zBin][centralityBin]->fEvt;
  fEvt->fPrimaryVertex = vertex;
  fCutProcessor->SetCentralityBin(centralityBin+1);
//////////////////////////////////////////////////////////////////
  //v0 tester
////////////////////////////////////////////////////////////////

  int v0Count = 0;
  vector<int> lambdaCount(fNumberOfTopologicalCutValues,0);
  vector<int> antiLambdaCount(fNumberOfTopologicalCutValues,0);
  TH1F *mcTruthOriginHist = NULL;
  if(fIsMCEvent) mcTruthOriginHist = CreateLambdaOriginHist(mcArray,numberOfLastHijingLabel); //Find MC truths of all the MC particles (before detector effects)
  TH1F *v0OriginHist = new TH1F("v0OriginHist", "Lambda Origins", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1);
  TH2F *v0PassedCutsOriginHist = new TH2F("v0PassedCutsOriginHist", "Lambda Origins", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, fNumberOfTopologicalCutValues, -0.5, fNumberOfTopologicalCutValues -0.5);
  for(int i = 0; i < fAOD->GetNumberOfV0s(); i++)
  {
    // cout<<"Testing v0s"<<endl;
    //Loop over all the V0 candidates to look for (anti)Lambdas
    bool hasPiPlusDaughter     = kFALSE;
    bool hasPiMinusDaughter    = kFALSE;
    bool hasProtonDaughter     = kFALSE;
    bool hasAntiProtonDaughter = kFALSE;
    AliAODv0* v0 = fAOD->GetV0(i);
    if(!v0) continue;
    //Make sure the V0 satisifies a few basic criteria
    if(v0->GetNDaughters() > 2)                          continue;
    if(v0->GetNProngs() > 2)                             continue;
    if(v0->GetCharge() != 0)                             continue;
    if(v0->ChargeProng(0) == v0->ChargeProng(1))         continue;
    if(v0->CosPointingAngle(primaryVertexAOD) < 0.998) continue;
    //Now look at daughter tracks
    AliAODTrack* daughterTrackPos = (AliAODTrack*)v0->GetDaughterLabel(0);
    AliAODTrack* daughterTrackNeg = (AliAODTrack*)v0->GetDaughterLabel(1);
    if(!daughterTrackPos) continue; //Daughter tracks must exist
    if(!daughterTrackNeg) continue;

    // The V0 has passed the most basic track cuts.
    fV0Candidates++;
    if(v0->GetOnFlyStatus())				 continue;

    daughterTrackPos->SetAODEvent(fAOD); //Need to set this for PID purposes
    daughterTrackNeg->SetAODEvent(fAOD);
    //Need to manually apply ITS refit cut for hybrid tracks in AOD115
    if((daughterTrackPos->GetStatus() & AliVTrack::kTPCrefit)==0) continue;
    if((daughterTrackNeg->GetStatus() & AliVTrack::kTPCrefit)==0) continue;
    AliPIDResponse::EDetPidStatus statusPosTPC = fpidAOD->CheckPIDStatus(AliPIDResponse::kTPC,daughterTrackPos);
    AliPIDResponse::EDetPidStatus statusNegTPC = fpidAOD->CheckPIDStatus(AliPIDResponse::kTPC,daughterTrackNeg);
    if(AliPIDResponse::kDetPidOk != statusPosTPC) continue;
    if(AliPIDResponse::kDetPidOk != statusNegTPC) continue;
    if(daughterTrackPos->GetTPCNcls() < 80) continue;
    if(daughterTrackNeg->GetTPCNcls() < 80) continue;
    //Need to manually apply shared cluster cut for hybrid tracks in AOD115 and 124
    Double_t fracSharedClustersPos = Double_t(daughterTrackPos->GetTPCnclsS()) / Double_t(daughterTrackPos->GetTPCncls());
    Double_t fracSharedClustersNeg = Double_t(daughterTrackNeg->GetTPCnclsS()) / Double_t(daughterTrackNeg->GetTPCncls());
    if(fracSharedClustersPos  > 0.4) continue;
    if(fracSharedClustersNeg  > 0.4) continue;
    if(daughterTrackPos->Pt() < .16) continue;
    if(daughterTrackNeg->Pt() < .16) continue;
    if(fabs(daughterTrackPos->Eta()) > fEtaDaughter) continue;
    if(fabs(daughterTrackNeg->Eta()) > fEtaDaughter) continue;

    //Now we'll get particle origin and momentum truth for MC particles
    AliReconstructedV0::MCV0Origin_t mcV0Origin = AliReconstructedV0::kUnassigned;
    TVector3 v0MomentumTruth(0., 0., 0.);
    if(fIsMCEvent){
      //first reject injected particles
      if(IsInjectedParticle(v0,mcArray,numberOfLastHijingLabel)) continue;
      mcV0Origin = DetermineV0Origin(v0, mcArray);
      v0OriginHist->Fill(mcV0Origin);
      GetMCParticleMomentumTruth(v0MomentumTruth, v0, mcArray);
    }
    //Now perform daughter track PID using TPC
    if (daughterTrackPos->Pt() > 0.5) {// min-pt cut to fix proton PID issues in LHC10h data
      if(fabs(fpidAOD->NumberOfSigmasTPC(daughterTrackPos,AliPID::kProton))
	 < fSigmaCutTPCProton)   hasProtonDaughter = kTRUE;
    }
    if(fabs(fpidAOD->NumberOfSigmasTPC(daughterTrackPos,AliPID::kPion)) < fSigmaCutTPCPion)  hasPiPlusDaughter = kTRUE;
    if (daughterTrackNeg->Pt() > 0.5){
      if(fabs(fpidAOD->NumberOfSigmasTPC(daughterTrackNeg,AliPID::kProton))
	 < fSigmaCutTPCProton) hasAntiProtonDaughter = kTRUE;
    }
    if(fabs(fpidAOD->NumberOfSigmasTPC(daughterTrackNeg,AliPID::kPion)) < fSigmaCutTPCPion) hasPiMinusDaughter = kTRUE;
    //Use TOF PID info if available.  This overrides TPC PID results.
    if(daughterTrackPos->P() > fTOFLow)
    { // positive daughter PID
      AliPIDResponse::EDetPidStatus statusPosTOF = fpidAOD->CheckPIDStatus(AliPIDResponse::kTOF,daughterTrackPos);
      if (AliPIDResponse::kDetPidOk == statusPosTOF) { // TOF signal is available for PID
	Float_t probMis = fpidAOD->GetTOFMismatchProbability(daughterTrackPos);
	if (probMis < 0.01) { // avoid TOF-TPC mismatch
	  //PiPlus
	  if(fabs(fpidAOD->NumberOfSigmasTOF(daughterTrackPos,AliPID::kPion)) < fSigmaCutTOFPion) hasPiPlusDaughter = kTRUE;
	  else hasPiPlusDaughter = kFALSE;
	  //Proton
	  if(fabs(fpidAOD->NumberOfSigmasTOF(daughterTrackPos,AliPID::kProton)) < fSigmaCutTOFProton) hasProtonDaughter = kTRUE;
	  else hasProtonDaughter = kFALSE;
	}
      }
    }
    if(daughterTrackNeg->P() > fTOFLow)
    { //negative daughter PID
      AliPIDResponse::EDetPidStatus statusNegTOF = fpidAOD->CheckPIDStatus(AliPIDResponse::kTOF,daughterTrackNeg);
      if (AliPIDResponse::kDetPidOk == statusNegTOF) { // TOF signal is available for PID
	Float_t probMis = fpidAOD->GetTOFMismatchProbability(daughterTrackNeg);
	if (probMis < 0.01) {
	  //PiMinus
	  if(fabs(fpidAOD->NumberOfSigmasTOF(daughterTrackNeg,AliPID::kPion)) < fSigmaCutTOFPion) hasPiMinusDaughter = kTRUE;
	  else hasPiMinusDaughter = kFALSE;
	  //AntiProton
	  if(fabs(fpidAOD->NumberOfSigmasTOF(daughterTrackNeg,AliPID::kProton)) < fSigmaCutTOFProton) hasAntiProtonDaughter = kTRUE;
	  else hasAntiProtonDaughter = kFALSE;
	}
      }
    }
    //If V0 doesn't have the right daughter combinations,
    //move on to the next candidate
    if(!((hasProtonDaughter && hasPiMinusDaughter) || (hasAntiProtonDaughter && hasPiPlusDaughter))) continue;

    //Save V0 information
    AliReconstructedV0 &thisV0 = fEvt->fReconstructedV0[v0Count];

    thisV0.v0Momentum.SetXYZ(v0->Px(), v0->Py(), v0->Pz());
    thisV0.v0MomentumTruth = v0MomentumTruth;
    thisV0.v0Pt     = v0->Pt();
    thisV0.v0Eta    = v0->Eta();
    thisV0.v0Phi    = v0->Phi();
    thisV0.massLam  = v0->MassLambda();
    thisV0.massALam = v0->MassAntiLambda();
    thisV0.massLamDifference  = fabs(v0->MassLambda() - fPDGLambda);
    thisV0.massALamDifference = fabs(v0->MassAntiLambda() - fPDGLambda);
    thisV0.massK0  = v0->MassK0Short();
    thisV0.lorentzGammaLam = v0->ELambda()/fPDGLambda;
    thisV0.v0DCA   = v0->DcaV0ToPrimVertex();
    thisV0.decayLength = v0->DecayLength(primaryVertexAOD);
    // thisV0.decayVertexPosition[0] = v0->DecayVertexV0X();
    // thisV0.decayVertexPosition[1] = v0->DecayVertexV0Y();
    // thisV0.decayVertexPosition[2] = v0->DecayVertexV0Z();
    thisV0.cosPointing = v0->CosPointingAngle(primaryVertexAOD);
    thisV0.hasProtonDaughter     = hasProtonDaughter;
    thisV0.hasAntiProtonDaughter = hasAntiProtonDaughter;
    thisV0.hasPiPlusDaughter     = hasPiPlusDaughter;
    thisV0.hasPiMinusDaughter    = hasPiMinusDaughter;
    thisV0.mcOriginType = mcV0Origin;
    //Save Daughter information
    thisV0.daughter1ID = v0->GetNegID();
    thisV0.daughter2ID = v0->GetPosID();
    thisV0.daughterPosMomentum = TVector3(v0->MomPosX(), v0->MomPosY(), v0->MomPosZ());
    thisV0.daughterPosProtonE = v0->EPosProton();
    thisV0.daughterPosPionE = v0->EPosPion();
    thisV0.daughterNegMomentum = TVector3(v0->MomNegX(), v0->MomNegY(), v0->MomNegZ());
    thisV0.daughterNegProtonE = v0->ENegProton();
    thisV0.daughterNegPionE = v0->ENegPion();
    thisV0.daughtersDCA = v0->DcaV0Daughters();
    thisV0.daughterPosDCAPrimaryVertex = v0->DcaPosToPrimVertex();
    thisV0.daughterNegDCAPrimaryVertex = v0->DcaNegToPrimVertex();
    thisV0.daughterPosGlobalPositions = GetGlobalPositionAtGlobalRadiiThroughTPC(daughterTrackPos, bfield);// used for merging cuts later
    thisV0.daughterNegGlobalPositions = GetGlobalPositionAtGlobalRadiiThroughTPC(daughterTrackNeg, bfield);
    daughterTrackPos->GetXYZ(thisV0.daughterPosPositionDCA);
    daughterTrackNeg->GetXYZ(thisV0.daughterNegPositionDCA);
    daughterTrackPos->GetPxPyPz(thisV0.daughterPosMomentumDCA);
    daughterTrackNeg->GetPxPyPz(thisV0.daughterNegMomentumDCA);
    daughterTrackPos->GetCovarianceXYZPxPyPz(thisV0.daughterPosCovariance);
    daughterTrackNeg->GetCovarianceXYZPxPyPz(thisV0.daughterNegCovariance);

    vector<TVector3> posCorrectedVector;
    vector<TVector3> negCorrectedVector;
    for(int iRad = 0; iRad <9; iRad++) {
      //find the track positions relative to the primary vertex position at different radii.
      TVector3 posCorrected = thisV0.daughterPosGlobalPositions[iRad] - vertex;
      TVector3 negCorrected = thisV0.daughterNegGlobalPositions[iRad] - vertex;
      posCorrectedVector.push_back(posCorrected);
      negCorrectedVector.push_back(negCorrected);
    }
    thisV0.daughterPosCorrectedGlobalPositions = posCorrectedVector;
    thisV0.daughterNegCorrectedGlobalPositions = negCorrectedVector;

    //Now analyze and histogram the V0
    fCutProcessor->CheckIfV0PassesCuts(& thisV0);
    fCutProcessor->DoV0Histogramming(& thisV0);
    FillReconstructedV0MCOrigin(& thisV0, v0PassedCutsOriginHist);
    AddV0ToMultiplicityCounts(& thisV0, lambdaCount, antiLambdaCount);
    FillTPCSignalHists(& thisV0, daughterTrackPos->P(), daughterTrackPos->GetTPCsignal(), daughterTrackNeg->P(), daughterTrackNeg->GetTPCsignal());
    if(fIsMCEvent) CheckForFakeV0s(& thisV0, fMCFakeParticleIdentity, fMCOtherV0Identity, mcV0Origin);

    if(fIsMCEvent) {
      // Fill some MC information for lambdas, fill for antilambdas
      if(thisV0.isLamCenter[fNominalTopCutIndex]) {
      	fMCTruthPtLam->Fill(thisV0.v0MomentumTruth.Pt());
      	fMCTruthPhiLam->Fill(thisV0.v0MomentumTruth.Phi());
    	fMCTruthEtaLam->Fill(thisV0.v0MomentumTruth.Eta());
      }
      if(thisV0.isALamCenter[fNominalTopCutIndex]) {
      	fMCTruthPtALam->Fill(thisV0.v0MomentumTruth.Phi());
    	fMCTruthPhiALam->Fill(thisV0.v0MomentumTruth.Phi());
    	fMCTruthEtaALam->Fill(thisV0.v0MomentumTruth.Eta());
      }
    }

    //Increment V0 count and check that we don't exceed size of V0 array
    v0Count++;
    if(fMaxV0Mult <= v0Count){
      cerr<<"V0 Count has exceeded"<<fMaxV0Mult<<"!"<<endl;
      break;
    }
  } //End of V0 loop
  //cout<<"Finished with V0 storage.  V0 candidate count is "<<v0Count<<endl;

  fEvt->fNumberCandidateV0 = v0Count;
  if(fIsMCEvent) BinOriginInformationForMCParticles(mcTruthOriginHist, v0OriginHist, v0PassedCutsOriginHist);

  //The following histograms don't get used again, so clean them up
  if(mcTruthOriginHist){
    delete mcTruthOriginHist;
    mcTruthOriginHist = NULL;
  }
  if(v0OriginHist){
    delete v0OriginHist;
    v0OriginHist = NULL;
  }
  if(v0PassedCutsOriginHist){
    delete v0PassedCutsOriginHist;
    v0PassedCutsOriginHist = NULL;
  }


  DoV0JudgmentCuts(fEvt, v0Count);
  HistogramEventMultiplicities(lambdaCount, antiLambdaCount, centralityBin);
  fTotalLambda += lambdaCount[fNominalTopCutIndex];
  fTotalAntiLambda += antiLambdaCount[fNominalTopCutIndex];
  //Printf("Reconstruction Finished. Starting pair studies.");

  //Now look at pairs for correlation function binning
  DoPairStudies(fEvt, centralityBin);
  //cout<<"Pair studies completed.  Event finished"<<endl;

  // Post output data.
  PostData(1, fOutputList);
}
//________________________________________________________________________
void AliAnalysisV0Lam::Terminate(Option_t *)
{
  // Called once at the end of the query
  // cout<<"Total Lambdas found:\t"<<fTotalLambda<<"."<<endl
  //     <<"Total AntiLambdas found:\t"<<fTotalAntiLambda<<"."<<endl
  //     <<"Total V0 Candidates found:\t"<<fV0Candidates<<"."<<endl
  //     <<"Done"<<endl;
}






//________________________________________________________________________
vector<TVector3> AliAnalysisV0Lam::GetGlobalPositionAtGlobalRadiiThroughTPC(const AliAODTrack *track, const Float_t bfield)
{
  // Gets the global position of the track at nine different radii in the TPC
  // track is the track you want to propagate
  // bfield is the magnetic field of your event
  // cout<<"Getting global position"<<endl;


  // The radii at which we get the global positions
  // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
  Float_t Rwanted[9]={85.,105.,125.,145.,165.,185.,205.,225.,245.};

  // Calculate positions with with GetXYZatR
  vector<TVector3> positions;
  // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam myEtp;
  myEtp.CopyFromVTrack(track);
  for(Int_t iRad = 0; iRad < 9; iRad++) {
    Double_t thisPosition[3] = {0., 0., 0.};
    myEtp.GetXYZatR(Rwanted[iRad], bfield, thisPosition, NULL);
    TVector3 posVec(thisPosition);
    positions.push_back(posVec);
  }

  return positions;
}



//________________________________________________________________________
Double_t AliAnalysisV0Lam::GetAverageSeparation(const vector<TVector3> &globalPositions1st, const vector<TVector3> &globalPositions2nd)
{
  // cout<<"Calculating avg sep"<<endl;
  //Compute the separation of two daughter tracks, averaged over 9 different positions
  Double_t radii[9]={85.,105.,125.,145.,165.,185.,205.,225.,245.};



  double avgSeparation = 0.;
  double pointsUsed = 0.;
  for(int iRad = 0; iRad < 9; iRad++) {
    Bool_t pointIsBad = kFALSE;
    if(fabs(globalPositions1st[iRad].Perp() - radii[iRad]) > 0.5) {
      pointIsBad = kTRUE;
      continue;
    }
    if(fabs(globalPositions2nd[iRad].Perp() - radii[iRad]) > 0.5) {
      pointIsBad = kTRUE;
      continue;
    }

    if(pointIsBad) continue;
    TVector3 diff = globalPositions1st[iRad] - globalPositions2nd[iRad];
    avgSeparation += diff.Mag();
    pointsUsed++;
  }

  if(pointsUsed < 1) return 0.;

  return avgSeparation/pointsUsed;
}


//________________________________________________________________________
void AliAnalysisV0Lam::DoV0JudgmentCuts(const AliAnalysisV0LamEvent * const event, const int totalV0s)
{
  // Looks at all V0s in a given event, and selectively removes V0s such
  // that each daughter track is claimed by no more than one V0.  This is
  // done by making judgment cuts.  The judgment cut compares a
  // characteristic (e.g. cosine of pointing angle) of two V0s that share a
  // daughter.  The V0 closer to the ideal value (e.g. cos(pointing) = 1) is
  // kept, while the other V0 is removed.

  // Occasionally several V0s will share a single daughter, or several V0s
  // will share several daughters.  Because of this, the possibility
  // exists for a V0 to be removed, and subsequently the V0 that it
  // competed with is also removed.  In that case the original V0 should be
  // "restored". This DoV0JudgmentCuts method includes an iterative process
  // which first removes V0s that fail the judgment cuts, and then
  // subsequently restores V0s which no longer compete with any other V0s.
  // This process of removing and restoring V0s continues until the
  // process stabilizes or 20 iterations occurs.

  //"Removed" V0s have a boolean flag "isDeemedUnworthy" which is set to
  // true.  Those V0s do not get used in correlation function pairs.

  // The selectionCriterion is used to set which of V0 DCA, daughter DCA to
  // each other, V0 cosine of pointing angle, or V0 mass is used as the
  // judgment cut.
  const int selectionCriterion = 0;
  // Start by looping over variable reconstruction cuts.  There will be
  // different lists of V0s for each reconstruction cut value, which may
  // lead to different sets of V0s competing over daughters.
  for (int cutIndex = 0; cutIndex < fNumberOfTopologicalCutValues; cutIndex++){
    bool converged;
    int iterations = 0;
    do { //Loop until the judgment cuts converge or 20 iterations pass
      converged = kTRUE;
      iterations++;
      for (int currentV0Number = 0; currentV0Number <totalV0s; currentV0Number++) { // Loop over each V0 in event
	// Mark V0s as bad if they are worse than other V0s with shared
	// daughters
	if(!(event->fReconstructedV0[currentV0Number].isLamCenter[cutIndex] || event->fReconstructedV0[currentV0Number].isALamCenter[cutIndex])) continue;
	// Don't bother if the current v0 isn't a center V0 (a center V0 has an m_inv inside the accepted mass window)
	if(event->fReconstructedV0[currentV0Number].isDeemedUnworthy[cutIndex]) continue;
	for (int comparisonV0Number = 0; comparisonV0Number <totalV0s; comparisonV0Number++)
	{ // Loop over all other V0s in event.
	  if(comparisonV0Number == currentV0Number) continue;
	  if(!(event->fReconstructedV0[comparisonV0Number].isLamCenter[cutIndex] || event->fReconstructedV0[comparisonV0Number].isALamCenter[cutIndex])) continue; //Don't bother if the comparison V0 isn't a center v0
	  if(!((event->fReconstructedV0[currentV0Number].daughter1ID == event->fReconstructedV0[comparisonV0Number].daughter1ID) || (event->fReconstructedV0[currentV0Number].daughter2ID == event->fReconstructedV0[comparisonV0Number].daughter2ID))) continue; //Don't bother if they don't share daughters
	  if(event->fReconstructedV0[comparisonV0Number].isDeemedUnworthy[cutIndex]) continue;
	  // If we reach this point in the loop, then these V0s compete
	  // over a daughter.  Compare them and determine which V0 needs
	  // to be removed.
	  int worseV0 = DetermineWhichV0IsWorse(event, currentV0Number, comparisonV0Number, selectionCriterion, cutIndex);
	  if (worseV0 != -1) event->fReconstructedV0[worseV0].isDeemedUnworthy[cutIndex] = kTRUE;
	  // A V0 has been removed, so process has not converged.
	  converged = kFALSE;
	  if(event->fReconstructedV0[currentV0Number].isDeemedUnworthy[cutIndex]) break;
	} //End loop over comparison V0s
      }// End V0 removal
      // Now Restore V0s if they no long compete over a daughter OR if they
      // compete and are judged to be good.
      for (int currentV0Number = 0; currentV0Number <totalV0s; currentV0Number++) { // Loop over each V0 in event
	if(!event->fReconstructedV0[currentV0Number].isDeemedUnworthy[cutIndex]) continue; //Only look at V0s that have been removed
	bool stillCompeting = kFALSE; // Assume that they no longer compete
	for (int comparisonV0Number = 0; comparisonV0Number < totalV0s; comparisonV0Number++)
	{ // Loop over all other V0s
	  if(comparisonV0Number == currentV0Number) continue;
	  if(!(event->fReconstructedV0[comparisonV0Number].isLamCenter[cutIndex] || event->fReconstructedV0[comparisonV0Number].isALamCenter[cutIndex])) continue; //Don't bother if the V0 is a center V0
	  if(event->fReconstructedV0[comparisonV0Number].isDeemedUnworthy[cutIndex]) continue; //Only compare with V0s that HAVE NOT been removed
	  if(!((event->fReconstructedV0[currentV0Number].daughter1ID == event->fReconstructedV0[comparisonV0Number].daughter1ID) || (event->fReconstructedV0[currentV0Number].daughter2ID == event->fReconstructedV0[comparisonV0Number].daughter2ID))) continue; //Don't bother if they don't share daughters
	  // If we reach this point in the loop, then these V0s compete
	  // over a daughter.  Compare them and determine which V0 needs
	  // to be removed.
	  int worseV0 = DetermineWhichV0IsWorse(event, currentV0Number, comparisonV0Number, selectionCriterion, cutIndex);
	  if(worseV0 == -1)
	  {
	    //Something has gone wrong.
	    cerr<<"Could not determine which V0 is worse"<<endl;
	  }
	  else if (worseV0 == currentV0Number)
	  {
	    // The V0 is still competing with another V0, and it has failed
	    // the judgment cut, so it stays removed.
	    stillCompeting = kTRUE;
	    break; //No need to keep comparing with other V0s
	  }
	  else {
	    //The comparison V0 is worse.  However, do nothing here.
	    //The comparison V0 will be removed in the removal loop
	    //if it still competes at that time.
	  }
	} //end comparison V0 loop
	if(!stillCompeting)
	{
	  //The V0 no longer fails any judgment cuts.  Restore it.
	  event->fReconstructedV0[currentV0Number].isDeemedUnworthy[cutIndex] = kFALSE;
	  converged = kFALSE;
	}
      }//End V0 restoration.
    } while ((converged == kFALSE) && (iterations < 20));
  } //End loop over variable reconstruction cuts.
  return;
}

//________________________________________________________________________
int AliAnalysisV0Lam::DetermineWhichV0IsWorse(const AliAnalysisV0LamEvent * const event, const int V01, const int V02, const int Criterion, const int cutIndex)
{
  // Performs a judgment cut on two V0 by comparing characteristics of those V0
  // and looking to see which of those V0 is further from the ideal value.
  // Cut is only performed on two V0 that claim the same daughter track.
  int worseV0 = -1;
  if (Criterion == 0)//compare using DCA to primary vertex
  {
    if (event->fReconstructedV0[V01].v0DCA <
	event->fReconstructedV0[V02].v0DCA) worseV0=V02;
    else worseV0=V01;
  }
  else if (Criterion == 1)//compare using DCA of daughters
  {
    if(event->fReconstructedV0[V01].daughtersDCA <
       event->fReconstructedV0[V02].daughtersDCA) worseV0 = V02;
    else worseV0 = V01;
  }
  else if (Criterion == 2)//compare using cos(pointing) of V0s
  {
    if(event->fReconstructedV0[V01].cosPointing >
       event->fReconstructedV0[V02].cosPointing) worseV0 = V02;
    else worseV0 = V01;
  }

  else if (Criterion == 3)//compare using Minv
  {
    double deltaM1=500.;
    double deltaM2=500.;
    if(event->fReconstructedV0[V01].isLamCenter[cutIndex]){
      deltaM1 = event->fReconstructedV0[V01].massLamDifference;
    }
    else if(event->fReconstructedV0[V01].isALamCenter[cutIndex]){
      deltaM1 = event->fReconstructedV0[V01].massALamDifference;
    }
    if(event->fReconstructedV0[V02].isLamCenter[cutIndex]){
      deltaM2 = event->fReconstructedV0[V02].massLamDifference;
    }
    else if(event->fReconstructedV0[V02].isALamCenter[cutIndex]){
      deltaM2 = event->fReconstructedV0[V02].massALamDifference;
    }
    if(deltaM1 <= deltaM2) worseV0 = V02;
    else worseV0 = V01;
  }
  else cerr<<"Invalid judgment cut criterion selected: "<<Criterion<<endl;
  return worseV0;
}

//________________________________________________________________________
TH1F *AliAnalysisV0Lam::CreateLambdaOriginHist(TClonesArray *mcArray, Int_t numberOfLastHijingLabel)
{
  //Create a histogram of the MC truth origin of each (anti)Lambda in the event
  //This allows us to see how many primary and secondary lambda there are in the
  //event.  We'll count again after all reconstruction is done to get an idea
  //of our reconstruction efficiency.
  TH1F *mcTruthOriginHist = new TH1F("mcTruthOriginHist", "Lambda Origins", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1);
  for (int i=0; i < mcArray->GetEntriesFast(); i++){
    AliAODMCParticle *mcParticle = (AliAODMCParticle*)mcArray->At(i);
    if(mcParticle->GetNDaughters() != 2) continue;
    //Reject injected particles.  Injected particles have a label greater
    //than numberOfLastHijingLabel.  Many secondary particles also have a label
    //greater than numberOfLastHijingLabel.  So first check if a particle has a
    //parent.  If it does, check if that parent is "original".  If not
    //original, reject that particle (e.g. this will reject lambdas with
    //injected cascade parents).  If a particle has no parent and it has
    //a label greater than numberOfLastHijingLabel, reject it (because it is
    //injected)
    if(mcParticle->GetMother() > -1){ // This MCParticle has a mother
      AliAODMCParticle *mcMother = (AliAODMCParticle*)mcArray->At(mcParticle->GetMother());
      if(!mcMother) continue;
      //Reject this MCParticle if its mother is injected
      if(mcMother->GetLabel() > numberOfLastHijingLabel) continue;
    }
    //If this MCParticle has no mother but it has a label > LastHijingLabel,
    //reject it
    else if(mcParticle->GetLabel() > numberOfLastHijingLabel) continue;
    AliAODMCParticle *mcDaughter1 = (AliAODMCParticle*)mcArray->At(mcParticle->GetDaughterLabel(0));
    AliAODMCParticle *mcDaughter2 = (AliAODMCParticle*)mcArray->At(mcParticle->GetDaughterLabel(1));
    //We won't count any MC Particles that have daughters outside the acceptance
    //region
    if(fabs(mcDaughter1->Eta()) > fEtaDaughter) continue;
    if(fabs(mcDaughter2->Eta()) > fEtaDaughter) continue;
    if(mcDaughter1->Pt() < 0.16) continue;
    if(mcDaughter2->Pt() < 0.16) continue;
    //Finally, get the PDG code of the MCParticle (or of its parent in the case
    //of secondary particles)
    AliReconstructedV0::MCV0Origin_t mcParticleOrigin = DeterminePdgCodeOfMcParticle(mcParticle,mcArray);
    mcTruthOriginHist->Fill(mcParticleOrigin);
  }
  return mcTruthOriginHist;
}

//________________________________________________________________________
void AliAnalysisV0Lam::FillReconstructedV0MCOrigin(const AliReconstructedV0 * v0, TH2F *histPassedCutsOrigin)
{
  //Make a histogram showing the MCTruth particle type of reconstructed V0s
  //(or the type of the mother particle if the V0 is secondary).
  for(int i = 0; i < fNumberOfTopologicalCutValues; i++){
    if(v0->isLamCenter[i]
       && (AliReconstructedV0::kFake == v0->mcOriginType))
    {
      histPassedCutsOrigin->Fill(AliReconstructedV0::kFakeLambda,i);
    }
    else if(v0->isLamCenter[i]){
      histPassedCutsOrigin->Fill(v0->mcOriginType,i);
    }
    if(v0->isALamCenter[i]
       && (AliReconstructedV0::kFake == v0->mcOriginType))
    {
      histPassedCutsOrigin->Fill(AliReconstructedV0::kFakeAntiLambda,i);
    }
    else if(v0->isALamCenter[i]){
      histPassedCutsOrigin->Fill(v0->mcOriginType,i);
    }
  }
}

//________________________________________________________________________
AliReconstructedV0::MCV0Origin_t AliAnalysisV0Lam::DetermineV0Origin(AliAODv0 *v0, TClonesArray *mcArray)
{
  // Determines the particle type (identity of it or of its parent particle)
  // from MC truth information
  AliReconstructedV0::MCV0Origin_t mcV0Origin = AliReconstructedV0::kUnassigned;
  //Get the MCParticle index for the V0
  int v0Id = GetV0MCParticleID(v0,mcArray);
  if (v0Id > 0){ // A real MC particle exists for this V0
    AliAODMCParticle* mcV0 = (AliAODMCParticle*)mcArray->At(v0Id);
    // Get the PDG code for this particle (or get the PDG code of its
    // mother in the case of secondary particles).
    // The PDG code gets converted into an MCV0Origin_t object
    mcV0Origin = DeterminePdgCodeOfMcParticle(mcV0,mcArray);
  }
  else{ // No MC truth exists for this V0.  It is fake.
    mcV0Origin = AliReconstructedV0::kFake;
  }
  return mcV0Origin;
}

//________________________________________________________________________
void AliAnalysisV0Lam::GetMCParticleMomentumTruth(TVector3 &pTruth, AliAODv0 *v0, TClonesArray *mcArray)
{
  // Get the MC truth of the 3D momentum of a V0.
  int v0Id = GetV0MCParticleID(v0,mcArray);
  if (v0Id > 0){
    AliAODMCParticle* mcV0 = (AliAODMCParticle*)mcArray->At(v0Id);
    pTruth.SetXYZ(mcV0->Px(), mcV0->Py(), mcV0->Pz());
  }
}

//________________________________________________________________________
int AliAnalysisV0Lam::GetV0MCParticleID(AliAODv0 *v0, TClonesArray *mcArray)
{
  // Returns the MCParticle index of a V0. Do this by finding the
  // corresponding MCParticles of the daughter tracks. If those daughter
  // MCParticles have the same mother, return the MCParticle index of that
  // mother.  If they don't have the same mother (or if both tracks are
  // primary) the V0 is a fake.  In that case, return -1.
  AliAODTrack* daughterTrackPos = (AliAODTrack*)v0->GetDaughterLabel(0);
  AliAODTrack* daughterTrackNeg = (AliAODTrack*)v0->GetDaughterLabel(1);
  daughterTrackPos->SetAODEvent(fAOD);
  daughterTrackNeg->SetAODEvent(fAOD);
  AliAODMCParticle* mcParticlePos = (AliAODMCParticle*)mcArray->At(abs(daughterTrackPos->GetLabel()));
  AliAODMCParticle* mcParticleNeg = (AliAODMCParticle*)mcArray->At(abs(daughterTrackNeg->GetLabel()));
  if(!(mcParticlePos) || !(mcParticleNeg)){
    //if either of these does not exist, V0 is fake.
    return -1;
  }
  //mcparticle->GetMother() will return a "-1" if the particle doesn't have a true mother (i.e. it's a fake track or primary)
  int motherOfPosID = mcParticlePos->GetMother();
  int motherOfNegID = mcParticleNeg->GetMother();
  if ((motherOfPosID > 0) && (motherOfPosID == motherOfNegID)){
    // Both daughter tracks refer to the same mother.  Return the MCParticle
    // index of that mother.
    return motherOfPosID;
  }
  else return -1; //Mother does not exist, or they refer to different
  // mothers. So this V0 is a fake
}

//________________________________________________________________________
AliReconstructedV0::MCV0Origin_t AliAnalysisV0Lam::DeterminePdgCodeOfMcParticle(AliAODMCParticle *mcParticle, TClonesArray *mcArray)
{
  // Get the PDG code for this particle (or get the PDG code of its
  // mother in the case of secondary particles)
  // The PDG code gets converted into an MCV0Origin_t object
  AliReconstructedV0::MCV0Origin_t mcParticleOrigin = AliReconstructedV0::kUnassigned;
  int v0PDG = mcParticle->GetPdgCode();
  //find if it has a parent and note the parent's pdg code
  int motherOfV0ID = mcParticle->GetMother();

  if(3122 == v0PDG){ //V0 is a Lambda
    if (motherOfV0ID <= 0) mcParticleOrigin = AliReconstructedV0::kPrimaryLambda;
    else { //V0 has a mother
      AliAODMCParticle* mcMotherOfV0 = (AliAODMCParticle*)mcArray->At(motherOfV0ID);
      int motherOfV0PDG = mcMotherOfV0->GetPdgCode();
      if(3212 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kPrimarySigmaZero;
      else if(3322 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kPrimaryCascadeZero;
      else if(3312 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kPrimaryCascadeMinus;
      else if(3334 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kPrimaryOmega;
      else if(3224 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kExcitedSigma;
      else if(3214 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kExcitedSigma;
      else if(3114 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kExcitedSigma;
      else {
	mcParticleOrigin = AliReconstructedV0::kOtherOriginLambda;
      }
    }
  }
  else if(-3122 == v0PDG){ // V0 is an Antilambda
    if (motherOfV0ID <= 0) mcParticleOrigin = AliReconstructedV0::kPrimaryAntiLambda;
    else { // V0 has a mother
      AliAODMCParticle* mcMotherOfV0 = (AliAODMCParticle*)mcArray->At(motherOfV0ID);
      int motherOfV0PDG = mcMotherOfV0->GetPdgCode();
      if(-3212 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kPrimaryAntiSigmaZero;
      else if(-3322 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kPrimaryAntiCascadeZero;
      else if(-3312 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kPrimaryAntiCascadePlus;
      else if(-3334 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kPrimaryAntiOmega;
      else if(-3224 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kExcitedAntiSigma;
      else if(-3214 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kExcitedAntiSigma;
      else if(-3114 == motherOfV0PDG) mcParticleOrigin = AliReconstructedV0::kExcitedAntiSigma;
      else {
	mcParticleOrigin = AliReconstructedV0::kOtherOriginAntiLambda;
      }
    }
  }
  else if(310 == mcParticle->GetPdgCode()) mcParticleOrigin = AliReconstructedV0::kKZeroShort;
  return mcParticleOrigin;
}

//________________________________________________________________________
void AliAnalysisV0Lam::BinOriginInformationForMCParticles(TH1F *mcOriginalV0Hist, TH1F *mcV0FinderHist, TH2F *mcV0PassedCutsHist)
{
  //Bin number of particles of each V0 type (e.g. primary lambda, lambda from Xi decay, etc.)
  //Also bin fraction of particles of each type still remaining at different stages
  SetBinsOnOriginHists(mcOriginalV0Hist); //Put labels on these histograms
  SetBinsOnOriginHists(mcV0FinderHist);
  SetBinsOnOriginHists(mcV0PassedCutsHist);
  TH1F *originalToV0Ratio = (TH1F*)mcV0FinderHist->Clone("originalToV0Ratio");
  originalToV0Ratio->Divide(mcOriginalV0Hist);
  for(int i = 0; i < AliReconstructedV0::kOriginTypeMax+1; i++){
    if(mcOriginalV0Hist->GetBinContent(i+1) >= 1){
      //Bin fraction of particles remaining at the V0 Finder stage
      //Only fill bin if original V0 distribution had content in that bin.
      //This avoids divide by zero problems
      fRemainingFromBeginningToV0Finder->Fill(i,originalToV0Ratio->GetBinContent(i+1));
    }
  }
  //Add yield results into output histograms
  fMCTruthOfOriginalParticles->Add(mcOriginalV0Hist);
  fMCTruthOfV0FinderParticles->Add(mcV0FinderHist);
  fMCTruthOfReconstructedParticles->Add(mcV0PassedCutsHist);
  delete originalToV0Ratio;
  for(int i = 0; i < fNumberOfTopologicalCutValues; i++){
    //Need to use a loop here because different variable cuts lead
    //to different distributions of reconstructed particles
    TH1F *originalToReconstructedRatio = (TH1F*)mcV0PassedCutsHist->ProjectionX("originalToReconstructedRatio",i+1,i+1);
    TH1F *v0FinderToReconstructedRatio = (TH1F*)mcV0PassedCutsHist->ProjectionX("v0FinderToReconstructedRatio",i+1,i+1);
    originalToReconstructedRatio->Divide(mcOriginalV0Hist);
    v0FinderToReconstructedRatio->Divide(mcV0FinderHist);
    for(int j = 0; j < AliReconstructedV0::kOriginTypeMax+1; j++){
      // Bin fraction of particles remaining after all reconstruction cuts
      // have been applied.  Only report a fraction remaining as zero if
      // the original event had one or more particles of that type
      if(mcOriginalV0Hist->GetBinContent(j+1) >= 1){
	fRemainingFromBeginningToRecon->Fill(i,j,originalToReconstructedRatio->GetBinContent(j+1));
      }
      if(mcV0FinderHist->GetBinContent(j+1) >= 1){
	fRemainingFromV0FinderToRecon->Fill(i,j,v0FinderToReconstructedRatio->GetBinContent(j+1));
      }
    }
    delete originalToReconstructedRatio;
    delete v0FinderToReconstructedRatio;
  }
}

//________________________________________________________________________
void AliAnalysisV0Lam::SetBinsOnOriginHists(TH1 *mcHist)
{
  mcHist->GetXaxis()->SetBinLabel(1,"OtherV0");
  mcHist->GetXaxis()->SetBinLabel(2,"Fake");
  mcHist->GetXaxis()->SetBinLabel(3,"Fake #Lambda");
  mcHist->GetXaxis()->SetBinLabel(4,"#Lambda");
  mcHist->GetXaxis()->SetBinLabel(5,"#Sigma0");
  mcHist->GetXaxis()->SetBinLabel(6,"#Sigma*");
  mcHist->GetXaxis()->SetBinLabel(7,"#Xi0");
  mcHist->GetXaxis()->SetBinLabel(8,"#Xi-");
  mcHist->GetXaxis()->SetBinLabel(9,"#Omega");
  mcHist->GetXaxis()->SetBinLabel(10,"Other #Lambda");
  mcHist->GetXaxis()->SetBinLabel(11,"Fake #bar{#Lambda}");
  mcHist->GetXaxis()->SetBinLabel(12,"#bar{#Lambda}");
  mcHist->GetXaxis()->SetBinLabel(13,"#bar{#Sigma}0");
  mcHist->GetXaxis()->SetBinLabel(14,"#bar{#Sigma}*");
  mcHist->GetXaxis()->SetBinLabel(15,"#bar{#Xi}0");
  mcHist->GetXaxis()->SetBinLabel(16,"#bar{#Xi}+");
  mcHist->GetXaxis()->SetBinLabel(17,"#bar{#Omega}");
  mcHist->GetXaxis()->SetBinLabel(18,"Other #bar{#Lambda}");
  mcHist->GetXaxis()->SetBinLabel(19,"K0s");
  return;
}

//________________________________________________________________________
void AliAnalysisV0Lam::SetBinsOnOriginHists(TH3 *mcHist)
{
  mcHist->GetYaxis()->SetBinLabel(1,"OtherV0");
  mcHist->GetYaxis()->SetBinLabel(2,"Fake");
  mcHist->GetYaxis()->SetBinLabel(3,"Fake #Lambda");
  mcHist->GetYaxis()->SetBinLabel(4,"#Lambda");
  mcHist->GetYaxis()->SetBinLabel(5,"#Sigma0");
  mcHist->GetYaxis()->SetBinLabel(6,"#Sigma*");
  mcHist->GetYaxis()->SetBinLabel(7,"#Xi0");
  mcHist->GetYaxis()->SetBinLabel(8,"#Xi-");
  mcHist->GetYaxis()->SetBinLabel(9,"#Omega");
  mcHist->GetYaxis()->SetBinLabel(10,"Other #Lambda");
  mcHist->GetYaxis()->SetBinLabel(11,"Fake #bar{#Lambda}");
  mcHist->GetYaxis()->SetBinLabel(12,"#bar{#Lambda}");
  mcHist->GetYaxis()->SetBinLabel(13,"#bar{#Sigma}0");
  mcHist->GetYaxis()->SetBinLabel(14,"#bar{#Sigma}*");
  mcHist->GetYaxis()->SetBinLabel(15,"#bar{#Xi}0");
  mcHist->GetYaxis()->SetBinLabel(16,"#bar{#Xi}+");
  mcHist->GetYaxis()->SetBinLabel(17,"#bar{#Omega}");
  mcHist->GetYaxis()->SetBinLabel(18,"Other #bar{#Lambda}");
  mcHist->GetYaxis()->SetBinLabel(19,"K0s");
  return;
}

//________________________________________________________________________
bool AliAnalysisV0Lam::IsInjectedParticle(AliAODv0 *v0, TClonesArray *mcArray, Int_t numberOfLastHijingLabel)
{
  //Check if a Monte Carlo particle comes from the base MC event, or if it
  //was injected.  Primary particles are injected if they have an
  //AliAODMCParticle::GetLabel() greater than AliGenHijingEventHeader::NProduced()-1
  bool isInjected = false;
  Int_t v0ID = GetV0MCParticleID(v0,mcArray);
  if(v0ID > -1){ //if the v0 comes from an actual MC particle V0
    AliAODMCParticle *mcParticle = (AliAODMCParticle*)mcArray->At(v0ID);
    if(mcParticle->GetMother() > -1){ //if it has a mother
      AliAODMCParticle *mcMother = (AliAODMCParticle*)mcArray->At(mcParticle->GetMother());
      if(!mcMother) return true; // if this doesn't exist, there was an error
      if(mcMother->GetLabel() > numberOfLastHijingLabel) isInjected = true;
    }
    else if(mcParticle->GetLabel() > numberOfLastHijingLabel) isInjected = true;
  }
  return isInjected;
}

//________________________________________________________________________s
bool AliAnalysisV0Lam::IsCorrectEventTrigger()
{
  //Pick out Central, SemiCentral, and MinBias events.  False if not using one of those event triggers.
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral));
  return isSelected;
}

//________________________________________________________________________
void AliAnalysisV0Lam::AddV0ToMultiplicityCounts(AliReconstructedV0 *v0, vector<int> & lambdaCount, vector<int> & antiLambdaCount)
{
  //If a V0 is a Lambda or an Antilambda, add one to the respective V0
  //yields for this event.  Depending on the variable cut value, the V0
  //may or may not get categorized as a (anti)Lambda.
  //This information is used for histogramming event multiplicities.
  for(int i = 0; i < fNumberOfTopologicalCutValues; i++){
    if(v0->isLamCenter[i]) lambdaCount[i]++;
    if(v0->isALamCenter[i]) antiLambdaCount[i]++;
  }
}

//________________________________________________________________________
void AliAnalysisV0Lam::HistogramEventMultiplicities(vector<int> & lambdaCount, vector<int> & antiLambdaCount, int centralityBin)
{
  //Add the event yields to the output yield histograms
  for(int i = 0; i < fNumberOfTopologicalCutValues; i++){
    //Centrality integrated histograms
    fMultDistLambda->Fill(i, lambdaCount[i]);
    fMultDistAntiLambda->Fill(i, antiLambdaCount[i]);
    //Centrality differential histograms
    fMultCentLambda->Fill(i, centralityBin+1, lambdaCount[i]);
    fMultCentAntiLambda->Fill(i, centralityBin+1, antiLambdaCount[i]);
  }
}

//________________________________________________________________________
void AliAnalysisV0Lam::FillTPCSignalHists(const AliReconstructedV0 *v0, const double posDaughterP, const double posDaughterTPCSignal, const double negDaughterP, const double negDaughterTPCSignal)
{
  //Histogram P vs TPCsignal for V0 daughters
  if(v0->isLamCenter[fNominalTopCutIndex])
  {
    fTPCVsPPosLam->Fill(posDaughterP,posDaughterTPCSignal);
    fTPCVsPNegLam->Fill(negDaughterP,negDaughterTPCSignal);
  }
  if(v0->isALamCenter[fNominalTopCutIndex])
  {
    fTPCVsPPosALam->Fill(posDaughterP,posDaughterTPCSignal);
    fTPCVsPNegALam->Fill(negDaughterP,negDaughterTPCSignal);
  }
}

//________________________________________________________________________
void AliAnalysisV0Lam::CheckForFakeV0s(const AliReconstructedV0 *v0, TH1F *mcFakeParticleIdentity, TH1F *mcOtherV0Identity, const AliReconstructedV0::MCV0Origin_t mcV0Origin)
{
  // Used in MC studies to determine how many reconstructed Lambdas and
  // Antilambdas are actually fake.  For simplicity, this is only done for
  // the default value of the variable reconstruction cut.
  if(v0->isLamCenter[fNominalTopCutIndex]
     || v0->isALamCenter[fNominalTopCutIndex])
  {
    if(AliReconstructedV0::kFake == mcV0Origin){
      //These V0s are fake
      if(v0->isLamCenter[fNominalTopCutIndex]) mcFakeParticleIdentity->Fill(0);
      if(v0->isALamCenter[fNominalTopCutIndex]) mcFakeParticleIdentity->Fill(1);
    }
    if(AliReconstructedV0::kUnassigned == mcV0Origin){
      // These V0s correspond to actual MC particles, but they aren't
      // Lambdas, Antilambdas, or K0s.
      if(v0->isLamCenter[fNominalTopCutIndex]) mcOtherV0Identity->Fill(0);
      if(v0->isALamCenter[fNominalTopCutIndex]) mcOtherV0Identity->Fill(1);
    }
  }
}

//________________________________________________________________________
double AliAnalysisV0Lam::CalculateKstar(TVector3 p1, TVector3 p2, double mass1, double mass2)
{
  Double_t e1 = sqrt(pow(mass1, 2) + p1.Mag2());
  Double_t e2 = sqrt(pow(mass2, 2) + p2.Mag2());

  TLorentzVector vec1(p1, e1);
  TLorentzVector vec2(p2, e2);

  TLorentzVector pDiff = vec1-vec2;
  TLorentzVector pSum = vec1+vec2;
  TLorentzVector kstarVec = 0.5*(pDiff - (pDiff * pSum/pSum.Mag2())*pSum);

  return -1.*kstarVec.Mag();
}

TVector3 AliAnalysisV0Lam::GetEmissionPoint(const AliAODMCParticle * const track, TVector3 primVertex)
{
  // Get a TVector3 for the emission point of an MC particle

  // Store global position into a TVector3
  TVector3 emisV(track->Xv(), track->Yv(), track->Zv());
  // Subtract away the primary vertex to get the emission point
  emisV -= primVertex;
  return emisV;
}

//________________________________________________________________________
void AliAnalysisV0Lam::DoPairStudies(const AliAnalysisV0LamEvent * const event, const Int_t centralityBin)
{
  //Loop over same- and mixed-event pairs and bin correlation function
  //numerators and denominators.

  for (Int_t topCutIndex = 0; topCutIndex < fNumberOfTopologicalCutValues; topCutIndex++)
  { // Start looping over all variable cut values
    for (Int_t i=0; i < event->fNumberCandidateV0; i++)
    { //Start looping over reconstructed V0s in this event
      AliReconstructedV0 &v01 = event->fReconstructedV0[i];

      Bool_t center1Lam  = v01.isLamCenter[topCutIndex];
      Bool_t center1ALam = v01.isALamCenter[topCutIndex];
      // Disregard V0 if it wasn't reconstructed as a center (anti)Lambda
      if (!(center1Lam || center1ALam)) continue;
      // Disregard V0 if it was removed via the judgment cuts
      if (v01.isDeemedUnworthy[topCutIndex]) continue;

      for (Int_t eventNumber=0; eventNumber<nEventsToMix+1; eventNumber++)
      { // Event buffer loop: eventNumber=0 is the current event, all other eventNumbers are past events
	Int_t startBin=0;
	// For same event pairs, start 2nd V0 loop at i+1 V0 to avoid
	// double counting
	if (eventNumber==0) startBin=i+1;
	for (Int_t j=startBin; j<(event+eventNumber)->fNumberCandidateV0; j++)
	{ // Second V0 loop (from past or current event)
	  AliReconstructedV0 &v02 = (event+eventNumber)->fReconstructedV0[j];
	  if (eventNumber==0)
	  { // Don't make pairs of V0s if they shared daughter tracks.
	    // This is redundant if the judgment cut is already employed
	    if (v01.daughter1ID == v02.daughter1ID) continue;
	    if (v01.daughter1ID == v02.daughter2ID) continue;
	    if (v01.daughter2ID == v02.daughter1ID) continue;
	    if (v01.daughter2ID == v02.daughter2ID) continue;
	  }
	  //Disregard second V0 if it was removed via judgment cuts
	  if (v02.isDeemedUnworthy[topCutIndex]) continue;
	  // A central V0 has a mass that falls within the accepted inv
	  // mass range.  Only make pairs with central V0s
	  Bool_t center2Lam  = v02.isLamCenter[topCutIndex];
	  Bool_t center2ALam = v02.isALamCenter[topCutIndex];
	  if (!(center2Lam || center2ALam)) continue;

	  // Determine the pair type
	  PairType pairType;
	  if (center1Lam && center2Lam) {
	    pairType = kLamLam;
	  } else if (center1ALam && center2ALam) {
	    pairType = kALamALam;
	  } else if (center1Lam && center2ALam) {
	    pairType = kLamALam;
	  } else if (center1ALam && center2Lam) {
	    pairType = kALamLam;
	  } else {
	    cout<<"Error: AliAnalysisV0Lam::DoPairStudies - Not a valid pair type"<<endl;
	    continue;
	  }

	  // Now we calculate a bunch of values that are used later during
	  // histogramming
	  // Double_t pairKt = (v01.v0Momentum + v02.v0Momentum).Pt()/2;
	  //Calculate k* for V0s and daughters using different mass assumptions
    // Double_t pairKstarLam =
    CalculateKstar(v01.v0Momentum, v02.v0Momentum, fPDGLambda,fPDGLambda);

	  // Check that the pair passes our pair cuts
	  vector<Bool_t> avgSepCutResults = CheckAvgSepCut(pairType, v01, v02);

	  if (topCutIndex == fNominalTopCutIndex) {
	    FillDecayLengthDiffHists(pairType, v01, v02, eventNumber);
	    FillAvgSepHists(pairType, v01, v02, eventNumber);
	  }

	  // If this passes the two-track cuts, fill correlation functions
	  for (UInt_t iTTC = 0; iTTC < avgSepCutResults.size(); iTTC++) {
	    if (avgSepCutResults[iTTC]) {
	      Int_t cutBin = 0;
	      if (kTopologicalStudy == fSysStudyType) {
		cutBin = topCutIndex;
	      } else if (kTwoTrackStudy == fSysStudyType) {
		cutBin = iTTC;
	      }
	      FillCorrelationHists(pairType, v01, v02, eventNumber, cutBin, centralityBin);

	      // Fill the momentum resolution matrices. We only
	      // want to fill this once, so make sure that there is
	      // no systematic cut study going on, or that we are
	      // using the default cut index (1)
	      // Do the same for the k* vs kT histograms
	      if (kNoStudy == fSysStudyType || cutBin == 1) {
		FillMomentumResolutionMatrix(pairType, v01, v02, eventNumber);
		FillKtVsKstarHists(pairType, v01, v02, eventNumber, centralityBin);
	      }
	    }
	  }
	}//end past event
      }//end event buffer
    }//end current event
  }//end variable cut loop
}//end DoPairStudies()


void AliAnalysisV0Lam::FillDecayLengthDiffHists(const PairType pairType, const AliReconstructedV0 &v01, const AliReconstructedV0 &v02, Bool_t isMixedEvent)
{
  // Fill histograms with the difference in proper decay lengths of the two particles

  if (v01.lorentzGammaLam <= 0) return; // Avoid divide by zero errors if something is wrong
  if (v02.lorentzGammaLam <= 0) return;

  Double_t properLength1 = v01.decayLength/v01.lorentzGammaLam;
  Double_t properLength2 = v02.decayLength/v02.lorentzGammaLam;

  Double_t difference = fabs(properLength1 - properLength2);

  if (kLamLam == pairType) {
    if (!isMixedEvent) {
      fHistSignalProperDecayLengthDiffLamLam->Fill(difference);
    } else {
      fHistBkgProperDecayLengthDiffLamLam->Fill(difference);
    }
  } else if (kALamALam == pairType) {
    if (!isMixedEvent) {
      fHistSignalProperDecayLengthDiffALamALam->Fill(difference);
    } else {
      fHistBkgProperDecayLengthDiffALamALam->Fill(difference);
    }
  } else if ((kLamALam == pairType) || (kALamLam == pairType)) {
    if (!isMixedEvent) {
      fHistSignalProperDecayLengthDiffLamALam->Fill(difference);
    } else {
      fHistBkgProperDecayLengthDiffLamALam->Fill(difference);
    }
  }
}

void AliAnalysisV0Lam::FillAvgSepHists(const PairType pairType, const AliReconstructedV0 &v01, const AliReconstructedV0 &v02, Bool_t isMixedEvent)
{
  // Calculate avg sep for all pairs of daughters
  Double_t avgSepPos = GetAverageSeparation(v01.daughterPosCorrectedGlobalPositions, v02.daughterPosCorrectedGlobalPositions);
  Double_t avgSepNeg = GetAverageSeparation(v01.daughterNegCorrectedGlobalPositions, v02.daughterNegCorrectedGlobalPositions);
  Double_t avgSepNegPos = GetAverageSeparation(v01.daughterNegCorrectedGlobalPositions, v02.daughterPosCorrectedGlobalPositions);
  Double_t avgSepPosNeg = GetAverageSeparation(v01.daughterPosCorrectedGlobalPositions, v02.daughterNegCorrectedGlobalPositions);

  Double_t dau1PosPt = v01.daughterPosMomentum.Pt();
  Double_t dau1NegPt = v01.daughterNegMomentum.Pt();
  Double_t dau2PosPt = v02.daughterPosMomentum.Pt();
  Double_t dau2NegPt = v02.daughterNegMomentum.Pt();

  if (kLamLam == pairType) {
    if(!isMixedEvent) {
      fSignalLamLamProtSep->Fill(avgSepPos, dau1PosPt, dau2PosPt);
      fSignalLamLamPiMinusSep->Fill(avgSepNeg, dau1NegPt, dau2NegPt);
      fSignalLamLamPlusMinusSep->Fill(avgSepPosNeg, dau1PosPt, dau2NegPt);
      fSignalLamLamPlusMinusSep->Fill(avgSepNegPos, dau1NegPt, dau2PosPt);
    } else {
      fBkgLamLamProtSep->Fill(avgSepPos, dau1PosPt, dau2PosPt);
      fBkgLamLamPiMinusSep->Fill(avgSepNeg, dau1NegPt, dau2NegPt);
      fBkgLamLamPlusMinusSep->Fill(avgSepPosNeg, dau1PosPt, dau2NegPt);
      fBkgLamLamPlusMinusSep->Fill(avgSepNegPos, dau1NegPt, dau2PosPt);
    }
  } else if (kALamALam == pairType) {
    if(!isMixedEvent) {
      fSignalALamALamAntiProtSep->Fill(avgSepNeg, dau1NegPt, dau2NegPt);
      fSignalALamALamPiPlusSep->Fill(avgSepPos, dau1PosPt, dau2PosPt);
      fSignalALamALamPlusMinusSep->Fill(avgSepPosNeg, dau1PosPt, dau2NegPt);
      fSignalALamALamPlusMinusSep->Fill(avgSepNegPos, dau1NegPt, dau2PosPt);
    } else {
      fBkgALamALamAntiProtSep->Fill(avgSepNeg, dau1NegPt, dau2NegPt);
      fBkgALamALamPiPlusSep->Fill(avgSepPos, dau1PosPt, dau2PosPt);
      fBkgALamALamPlusMinusSep->Fill(avgSepPosNeg, dau1PosPt, dau2NegPt);
      fBkgALamALamPlusMinusSep->Fill(avgSepNegPos, dau1NegPt, dau2PosPt);
    }
  } else if (kLamALam == pairType) {
    if(!isMixedEvent) {
      fSignalLamALamProtPiPlusSep->Fill(avgSepPos, dau1PosPt, dau2PosPt);
      fSignalLamALamAntiProtPiMinusSep->Fill(avgSepNeg, dau1NegPt, dau2NegPt);
      fSignalLamALamProtSep->Fill(avgSepPosNeg, dau1PosPt, dau2NegPt);
      fSignalLamALamPionSep->Fill(avgSepNegPos, dau1NegPt, dau2PosPt);
    } else {
      fBkgLamALamProtPiPlusSep->Fill(avgSepPos, dau1PosPt, dau2PosPt);
      fBkgLamALamAntiProtPiMinusSep->Fill(avgSepNeg, dau1NegPt, dau2NegPt);
      fBkgLamALamProtSep->Fill(avgSepPosNeg, dau1PosPt, dau2NegPt);
      fBkgLamALamPionSep->Fill(avgSepNegPos, dau1NegPt, dau2PosPt);
    }
  } else if (kALamLam == pairType) {
    if(!isMixedEvent) {
      fSignalLamALamProtPiPlusSep->Fill(avgSepPos, dau1PosPt, dau2PosPt);
      fSignalLamALamAntiProtPiMinusSep->Fill(avgSepNeg, dau1NegPt, dau2NegPt);
      fSignalLamALamProtSep->Fill(avgSepNegPos, dau1NegPt, dau2PosPt);
      fSignalLamALamPionSep->Fill(avgSepPosNeg, dau1PosPt, dau2NegPt);
    } else {
      fBkgLamALamProtPiPlusSep->Fill(avgSepPos, dau1PosPt, dau2PosPt);
      fBkgLamALamAntiProtPiMinusSep->Fill(avgSepNeg, dau1NegPt, dau2NegPt);
      fBkgLamALamProtSep->Fill(avgSepNegPos, dau1NegPt, dau2PosPt);
      fBkgLamALamPionSep->Fill(avgSepPosNeg, dau1PosPt, dau2NegPt);
    }
  }
}

void AliAnalysisV0Lam::FillCorrelationHists(const PairType pairType, const AliReconstructedV0 &v01, const AliReconstructedV0 &v02, const Bool_t isMixedEvent, const Int_t cutBin, const Int_t centralityBin)
{
  // Fill same event and mixed event k* histograms
  Double_t pairKstarLam = CalculateKstar(v01.v0Momentum, v02.v0Momentum, fPDGLambda,fPDGLambda);

  if (kLamLam == pairType) {
    if(!isMixedEvent) {
      fSignalLamLam->Fill(cutBin, centralityBin+1, pairKstarLam);
    } else {
      fBkgLamLam->Fill(cutBin, centralityBin+1, pairKstarLam);
    }
  } else if (kALamALam == pairType) {
    if(!isMixedEvent) {
      fSignalALamALam->Fill(cutBin, centralityBin+1, pairKstarLam);
    } else {
      fBkgALamALam->Fill(cutBin, centralityBin+1, pairKstarLam);
    }
  } else if ((kLamALam == pairType) || (kALamLam == pairType)) {
    if(!isMixedEvent) {
      fSignalLamALam->Fill(cutBin, centralityBin+1, pairKstarLam);
    } else {
      fBkgLamALam->Fill(cutBin, centralityBin+1, pairKstarLam);
    }
  }
}

void AliAnalysisV0Lam::FillKtVsKstarHists(const PairType pairType, const AliReconstructedV0 &v01, const AliReconstructedV0 &v02, const Bool_t isMixedEvent, const Int_t centralityBin)
{
  // Fill kT vs k* for same and mixed event pairs
  Double_t pairKstarLam = CalculateKstar(v01.v0Momentum, v02.v0Momentum, fPDGLambda,fPDGLambda);
  Double_t pairKtLam = (v01.v0Pt + v02.v0Pt) / 2;
  if (kLamLam == pairType) {
    if(!isMixedEvent) {
      fSignalKtVsKstarLamLam->Fill(pairKtLam, centralityBin+1, pairKstarLam);
    } else {
      fBkgKtVsKstarLamLam->Fill(pairKtLam, centralityBin+1, pairKstarLam);
    }
  } else if (kALamALam == pairType) {
    if(!isMixedEvent) {
      fSignalKtVsKstarALamALam->Fill(pairKtLam, centralityBin+1, pairKstarLam);
    } else {
      fBkgKtVsKstarALamALam->Fill(pairKtLam, centralityBin+1, pairKstarLam);
    }
  } else if ((kLamALam == pairType) || (kALamLam == pairType)) {
    if(!isMixedEvent) {
      fSignalKtVsKstarLamALam->Fill(pairKtLam, centralityBin+1, pairKstarLam);
    } else {
      fBkgKtVsKstarLamALam->Fill(pairKtLam, centralityBin+1, pairKstarLam);
    }
  }
}

vector<Bool_t> AliAnalysisV0Lam::CheckAvgSepCut(const PairType type, const AliReconstructedV0 &v01, const AliReconstructedV0 &v02)
{
  // Calculate avg sep for all pairs of daughters and check if they pass cuts

  vector<Bool_t> cutResults;
  if (fTestNoTTC) { // In case we are running with cuts turned off
    cutResults.push_back(kTRUE);
    return cutResults;
  }

  Double_t avgSepPos = GetAverageSeparation(v01.daughterPosCorrectedGlobalPositions, v02.daughterPosCorrectedGlobalPositions);
  Double_t avgSepNeg = GetAverageSeparation(v01.daughterNegCorrectedGlobalPositions, v02.daughterNegCorrectedGlobalPositions);
  Double_t avgSepNegPos = GetAverageSeparation(v01.daughterNegCorrectedGlobalPositions, v02.daughterPosCorrectedGlobalPositions);
  Double_t avgSepPosNeg = GetAverageSeparation(v01.daughterPosCorrectedGlobalPositions, v02.daughterNegCorrectedGlobalPositions);


  Int_t nVariableCuts = 1;
  if (kTwoTrackStudy == fSysStudyType) {
    nVariableCuts = 3;
  }

  // We may be varying one of the cut values. Loop over the variations.
  for (Int_t iVar = 0; iVar < nVariableCuts; iVar++) {
    // Initialize the nominal avg sep cut values.
    vector<Double_t> nominalCutValues(kNumberTTCTypes);
    nominalCutValues[kSameProtProt] = 12.; // Same sign prot-prot cut
    nominalCutValues[kSamePiPi]     = 10.; // Same sign pi-pi
    nominalCutValues[kSameProtPi]   = 10.; // Same sign prot-pi
    nominalCutValues[kDiffProtProt] = 10.; // Diff sign prot-prot
    nominalCutValues[kDiffPiPi]     = 25.; // Diff sign pi-pi
    nominalCutValues[kDiffProtPi]   = 15.; // Diff sign prot-pi


    // Vary one of the cut values
    if ((kTwoTrackStudy == fSysStudyType) && (iVar == 0)) {
      nominalCutValues[fVariableCutType] *= 0.9;
    } else if ((kTwoTrackStudy == fSysStudyType) && (iVar == 2)) {
      nominalCutValues[fVariableCutType] *= 1.1;
    }

    // Do different checks for each pair type
    Bool_t doesPass = kFALSE;
    if (kLamLam == type) {
      if (    (nominalCutValues[kSameProtProt] < avgSepPos)
	   && (nominalCutValues[kSamePiPi]     < avgSepNeg)
	   && (nominalCutValues[kDiffProtPi]   < avgSepPosNeg)
	   && (nominalCutValues[kDiffProtPi]   < avgSepNegPos)) {
	doesPass = kTRUE;
      }
    } else if (kALamALam == type) {
      if (    (nominalCutValues[kSamePiPi]     < avgSepPos)
	   && (nominalCutValues[kSameProtProt] < avgSepNeg)
	   && (nominalCutValues[kDiffProtPi]   < avgSepPosNeg)
	   && (nominalCutValues[kDiffProtPi]   < avgSepNegPos)) {
	doesPass = kTRUE;
      }
    } else if (kLamALam == type) {
      if (    (nominalCutValues[kSameProtPi]   < avgSepPos)
	   && (nominalCutValues[kSameProtPi]   < avgSepNeg)
	   && (nominalCutValues[kDiffProtProt] < avgSepPosNeg)
	   && (nominalCutValues[kDiffPiPi]     < avgSepNegPos)) {
	doesPass = kTRUE;
      }
    } else if (kALamLam == type) {
      if (    (nominalCutValues[kSameProtPi]   < avgSepPos)
	   && (nominalCutValues[kSameProtPi]   < avgSepNeg)
	   && (nominalCutValues[kDiffPiPi]     < avgSepPosNeg)
	   && (nominalCutValues[kDiffProtProt] < avgSepNegPos)) {
	doesPass = kTRUE;
      }
    } else {
      cerr<<"Error: AliAnalysisV0Lam::CheckAvgSepCut - Not a valid pair type!"<<endl;
    }
    cutResults.push_back(doesPass);
  } // end variable cut loop

  // return vector of check results
  return cutResults;
}


void AliAnalysisV0Lam::FillMomentumResolutionMatrix(const PairType type, const AliReconstructedV0 &v01, const AliReconstructedV0 &v02, Bool_t isMixedEvent)
{

  if(v01.v0MomentumTruth.Mag() < 0.0001) return;  //Not a real V0
  if(v02.v0MomentumTruth.Mag() < 0.0001) return;  //Not a real V0


  Double_t kstarRec = CalculateKstar(v01.v0Momentum, v02.v0Momentum, fPDGLambda, fPDGLambda);
  Double_t kstarTruth = CalculateKstar(v01.v0MomentumTruth, v02.v0MomentumTruth, fPDGLambda, fPDGLambda);

  // Fill the matrices with all pairs that pass reconstruction cuts, even if they aren't primary lambdas
  if (!isMixedEvent) {
    if (type == kLamLam) {
      fResMatrixLLSameAll->Fill(kstarTruth, kstarRec);
    } else if (type == kALamALam) {
      fResMatrixAASameAll->Fill(kstarTruth, kstarRec);
    } else if (type == kLamALam || type == kALamLam) {
      fResMatrixLASameAll->Fill(kstarTruth, kstarRec);
    } else {
      cerr << "AliAnalysisV0Lam - Not a valid pair type. Cannot fill resolution matrix." << endl;
    }
  } else {
    if (type == kLamLam) {
      fResMatrixLLMixedAll->Fill(kstarTruth, kstarRec);
    } else if (type == kALamALam) {
      fResMatrixAAMixedAll->Fill(kstarTruth, kstarRec);
    } else if (type == kLamALam || type == kALamLam) {
      fResMatrixLAMixedAll->Fill(kstarTruth, kstarRec);
    } else {
      cerr << "AliAnalysisV0Lam - Not a valid pair type. Cannot fill resolution matrix." << endl;
    }
  }

  // Fill the matrices with pairs of true primary lambdas
  if ((v01.mcOriginType == AliReconstructedV0::kPrimaryLambda || v01.mcOriginType == AliReconstructedV0::kPrimaryAntiLambda) &&
      (v02.mcOriginType == AliReconstructedV0::kPrimaryLambda || v02.mcOriginType == AliReconstructedV0::kPrimaryAntiLambda)) {
    if (!isMixedEvent) {
      if (type == kLamLam) {
	fResMatrixLLSamePure->Fill(kstarTruth, kstarRec);
      } else if (type == kALamALam) {
	fResMatrixAASamePure->Fill(kstarTruth, kstarRec);
      } else if (type == kLamALam || type == kALamLam) {
	fResMatrixLASamePure->Fill(kstarTruth, kstarRec);
      } else {
	cerr << "AliAnalysisV0Lam - Not a valid pair type. Cannot fill resolution matrix." << endl;
      }
    } else {
      if (type == kLamLam) {
	fResMatrixLLMixedPure->Fill(kstarTruth, kstarRec);
      } else if (type == kALamALam) {
	fResMatrixAAMixedPure->Fill(kstarTruth, kstarRec);
      } else if (type == kLamALam || type == kALamLam) {
	fResMatrixLAMixedPure->Fill(kstarTruth, kstarRec);
      } else {
	cerr << "AliAnalysisV0Lam - Not a valid pair type. Cannot fill resolution matrix." << endl;
      }
    }
  }


}


bool AliAnalysisV0Lam::RejectEventCentFlat(float MagField, float CentPercent)
{ // to flatten centrality distribution
  bool RejectEvent = kFALSE;
  int weightBinSign;
  TRandom3* fRandomNumber = new TRandom3();  //for 3D, random sign switching
  fRandomNumber->SetSeed(0);

  if(MagField > 0) weightBinSign = 0;
  else weightBinSign = 1;
  float kCentWeight[2][9] = {{.878,.876,.860,.859,.859,.88,.873,.879,.894},
                             {.828,.793,.776,.772,.775,.796,.788,.804,.839}};
  int weightBinCent = (int) CentPercent;
  if(fRandomNumber->Rndm() > kCentWeight[weightBinSign][weightBinCent]) RejectEvent = kTRUE;
  delete fRandomNumber;
  return RejectEvent;
}
