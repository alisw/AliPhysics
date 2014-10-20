#include <iostream>
#include <math.h>
#include "TChain.h"
#include "TFile.h"
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
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODv0.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisV0Lam.h"
#include "AliAODMCHeader.h"
#include "AliGenHijingEventHeader.h"
#define PI 3.1415927

// Analysis task for studying lambda-lambda femtoscopic correlations
// Author: Jai Salzwedel, jai.salzwedel@cern.ch


ClassImp(AliAnalysisV0Lam)

//________________________________________________________________________
AliAnalysisV0Lam::AliAnalysisV0Lam():
AliAnalysisTaskSE(),
  fAOD(0),
  fOutputList(0),
  fpidAOD(0)
{
}
//________________________________________________________________________
AliAnalysisV0Lam::AliAnalysisV0Lam(const char *name) 
: AliAnalysisTaskSE(name), 
  fAOD(0), 
  fOutputList(0),
  fpidAOD(0)
{
  // Define output slots here 
  // Output slot #1
  DefineOutput(1, TList::Class());
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
  // Set global variables
  fEventCount=           0;
  fPDGLambda =    1.115683;
  fPDGK0     =     .497614;
  fPDGProton =     .938272;
  fPDGPion   =    .1395702;
  fEtaDaughter =       0.8;    //max eta of daughter particles - default 0.8
  fMassWindowK0 =    0.018;    //accept .480-.514
  fMassWindowLam = 0.00568;    //accept 1.11003-1.121363
  fTOFLow =            0.8;    // Lower |P| limit
  fSigmaCutTPCProton = 3.0;    // max Nsigma allowed 
  fSigmaCutTPCPion   = 3.0;
  fSigmaCutTOFProton = 4.0;
  fSigmaCutTOFPion   = 4.0;
  fIsUsingVariableAvgSepCut = kFALSE; //Relevant for two-track cuts in DoPairStudies().
  fMaxV0Mult = 700; 
  int numberVariableAvgSepCuts = 12;
  

  // Setup V0 cut processor
  fCutProcessor = new AliAnalysisV0LamCutProcessing(fOutputList);
  fNumberOfVariableCutValues = fCutProcessor->GetNumberOfVariableCutValues();
  fDefaultVariableCutIndex = fCutProcessor->GetVariableCutIndex();
  if(fIsUsingVariableAvgSepCut){
    fNumberOfCfVariableCutValues = numberVariableAvgSepCuts;
  }
  else fNumberOfCfVariableCutValues = fNumberOfVariableCutValues;
  cout<<"Number of variable cf cut values: "<<fNumberOfCfVariableCutValues<<endl;
  fTotalLambda = 	 0; //tabulates number of v0s found for locally run jobs
  fTotalAntiLambda =     0;

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
  
  fOutputList = new TList();
  fOutputList->SetOwner(); 
  MyInit();// Initialize my settings
  
  fMultDistRough = new TH1F("fMultDistRough","Multiplicity Distribution",301,-.5,3001-.5);
  fMultDistRough->GetXaxis()->SetTitle("Event Multiplicity (pions)");
  fMultDistRough->GetYaxis()->SetTitle("# of events");
  fOutputList->Add(fMultDistRough);

  fCentrality = new TH1F("fCentrality", "Centrality Percentage of Event", 100, 0., 100.);
  fOutputList->Add(fCentrality);

  //TPC signal vs track momentum
  fTPCVsPPosLam = new TH2F("fTPCVsPPosLam","",50,0.,2.,250,0.,500.);
  fOutputList->Add(fTPCVsPPosLam);
    //TPC signal vs track momentum
  fTPCVsPNegLam = new TH2F("fTPCVsPNegLam","",50,0.,2.,250,0.,500.);
  fOutputList->Add(fTPCVsPNegLam);
    //TPC signal vs track momentum
  fTPCVsPPosALam = new TH2F("fTPCVsPPosALam","",50,0.,2.,250,0.,500.);
  fOutputList->Add(fTPCVsPPosALam);
    //TPC signal vs track momentum
  fTPCVsPNegALam = new TH2F("fTPCVsPNegALam","",50,0.,2.,250,0.,500.);
  fOutputList->Add(fTPCVsPNegALam);

  //V0 Shared daughter culling statistics
  fDataCompeted = new TH1F("fDataCompeted","",26, -0.5, 25.5);
  fOutputList->Add(fDataCompeted);

  fDataCulled = new TH1F("fDataCulled","",26, -0.5, 25.5);
  fOutputList->Add(fDataCulled);

  fRemainingV0s = new TH1F("fRemainingV0s","",26, -0.5, 25.5);
  fOutputList->Add(fRemainingV0s);

  fRemainingFrac = new TH1F("fRemainingFrac","",101, -.005, 1.005);
  fOutputList->Add(fRemainingFrac);

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
  fRemainingFromBeginningToRecon = new TH3F("fRemainingFromBeginningToRecon", "Fraction Remaining After Reconstruction", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5,  AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, 21, -0.025, 1.025);
  fRemainingFromV0FinderToRecon = new TH3F("fRemainingFromV0FinderToRecon", "Fraction From V0Finder That Remain After Reconstruction Stage", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5,  AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, 21, -0.025, 1.025);
  fMCTruthOfReconstructedParticles = new TH2F("fMCTruthOfReconstructedParticles", "MC Truth of Reconstructed Particles in Event", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5);

  // Particle multiplicities
  fMultDistLambda = new TH2F("fMultDistLambda", "Lambda multiplicity", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5,  21, -0.5, 21-0.5);
  fMultDistAntiLambda = new TH2F("fMultDistAntiLambda", "AntiLambda multiplicity", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5,  21, -0.5, 21-0.5);
  fMultCentLambda = new TH3F("fMultCentLambda", "Lambda multiplicity vs centrality", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5,  nCentBins, .5, nCentBins+.5, 21, -0.5, 21-0.5);
  fMultCentAntiLambda =  new TH3F("fMultCentAntiLambda", "AntiLambda multiplicity vs centrality", fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5,  nCentBins, .5, nCentBins+.5, 21, -0.5, 21-0.5);
  SetBinsOnOriginHists(fRemainingFromBeginningToRecon);
  SetBinsOnOriginHists(fRemainingFromV0FinderToRecon);
  SetBinsOnOriginHists(fMCTruthOfReconstructedParticles);
  fMultDistLambda->GetXaxis()->SetTitle("Var Bin");
  fMultDistLambda->GetYaxis()->SetTitle("Event Multiplicity (Lambdas)");
  fMultDistAntiLambda->GetXaxis()->SetTitle("Var Bin");
  fMultDistAntiLambda->GetXaxis()->SetTitle("Event Multiplicity (AntiLambdas)");
  fMultCentLambda->GetXaxis()->SetTitle("Var Bin");
  fMultCentLambda->GetYaxis()->SetTitle("Centrality");
  fMultCentLambda->GetZaxis()->SetTitle("Event Multiplicity (Lambdas)");
  fMultCentAntiLambda->GetXaxis()->SetTitle("Var Bin");
  fMultCentAntiLambda->GetYaxis()->SetTitle("Centrality");
  fMultCentAntiLambda->GetZaxis()->SetTitle("Event Multiplicity (AntiLambdas)");
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
  int kTBins = 200;
  int maxKtBin = 10.;
  int kStarBins = 800;
  int maxKStar = 2.;
  //Pair kT Tracking: Centrality bins, kT bins, k* bins
  fKtLamLamSig = new TH3F("fKtLamLamSig", "LamLam Pair Kt Same Event", nCentBins, .5, nCentBins +.5, kTBins, 0., maxKtBin, kStarBins, 0., maxKStar);
  fOutputList->Add(fKtLamLamSig);
  fKtALamALamSig = new TH3F("fKtALamALamSig", "ALamALam Pair Kt Same Event", nCentBins, .5, nCentBins +.5, kTBins, 0., maxKtBin, kStarBins, 0., maxKStar);
  fOutputList->Add(fKtALamALamSig);
  fKtLamALamSig = new TH3F("fKtLamALamSig", "LamALam Pair Kt Same Event", nCentBins, .5, nCentBins +.5, kTBins, 0., maxKtBin, kStarBins, 0., maxKStar);
  fOutputList->Add(fKtLamALamSig);
    fKtLamLamBkg = new TH3F("fKtLamLamBkg", "LamLam Pair Kt Mixed Event", nCentBins, .5, nCentBins +.5, kTBins, 0., maxKtBin, kStarBins, 0., maxKStar);
  fOutputList->Add(fKtLamLamBkg);
  fKtALamALamBkg = new TH3F("fKtALamALamBkg", "ALamALam Pair Kt Mixed Event", nCentBins, .5, nCentBins +.5, kTBins, 0., maxKtBin, kStarBins, 0., maxKStar);
  fOutputList->Add(fKtALamALamBkg);
  fKtLamALamBkg = new TH3F("fKtLamALamBkg", "LamALam Pair Kt Mixed Event", nCentBins, .5, nCentBins +.5, kTBins, 0., maxKtBin, kStarBins, 0., maxKStar);
  fOutputList->Add(fKtLamALamBkg);

  //Momentum resolution (pre-)correction analysis
  fHistPsmearingKreconVsKtruthLL = new TH2F ("fHistPsmearingKreconVsKtruthLL","Relative momentum resolution, recon vs truth LL", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fOutputList->Add(fHistPsmearingKreconVsKtruthLL);
  fHistPsmearingKreconMinusKtruthLL = new TH1F ("fHistPsmearingKreconMinusKtruthLL","Relative momentum resolution, recon - truth LL", kStarBins, -0.5, 0.5);
  fOutputList->Add(fHistPsmearingKreconMinusKtruthLL);

    fHistPsmearingKreconVsKtruthAA = new TH2F ("fHistPsmearingKreconVsKtruthAA","Relative momentum resolution, recon vs truth AA", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fOutputList->Add(fHistPsmearingKreconVsKtruthAA);
  fHistPsmearingKreconMinusKtruthAA = new TH1F ("fHistPsmearingKreconMinusKtruthAA","Relative momentum resolution, recon - truth AA", kStarBins, -0.5, 0.5);
  fOutputList->Add(fHistPsmearingKreconMinusKtruthAA);

    fHistPsmearingKreconVsKtruthLA = new TH2F ("fHistPsmearingKreconVsKtruthLA","Relative momentum resolution, recon vs truth LA", kStarBins, 0., maxKStar, kStarBins, 0., maxKStar);
  fOutputList->Add(fHistPsmearingKreconVsKtruthLA);
  fHistPsmearingKreconMinusKtruthLA = new TH1F ("fHistPsmearingKreconMinusKtruthLA","Relative momentum resolution, recon - truth LA", kStarBins, -0.5, 0.5);
  fOutputList->Add(fHistPsmearingKreconMinusKtruthLA);
  
  /////////Signal Distributions///////////////////
  //First bin is variable cut value, second bin is centrality, third bin is Kstar
  fSignalLamLam = new TH3F("fSignalLamLam","Same Event Pair Distribution", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamLam);

  fBkgLamLam = new TH3F("fBkgLamLam","Mixed Event Pair Distribution", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamLam);

  fSignalALamALam = new TH3F("fSignalALamALam","Same Event Pair Distribution", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalALamALam);

  fBkgALamALam = new TH3F("fBkgALamALam","Mixed Event Pair Distribution", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgALamALam);

  fSignalLamALam = new TH3F("fSignalLamALam","Same Event Pair Distribution", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALam);

  fBkgLamALam = new TH3F("fBkgLamALam","Mixed Event Pair Distribution", fNumberOfCfVariableCutValues, -0.5, fNumberOfCfVariableCutValues -0.5, nCentBins, .5, nCentBins+.5, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALam);

  // Daughter pair separation distributions
  // Binned according to (average separation, qinv)
  int avgSepBins = 200;
  double avgSepMaxValue = 20.;
  fSignalLamLamProtSep = new TH2F ("fSignalLamLamProtSep","Proton pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamLamProtSep);

  fSignalLamLamPiMinusSep = new TH2F ("fSignalLamLamPiMinusSep","PiMinus pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamLamPiMinusSep);

  fSignalALamALamAntiProtSep = new TH2F ("fSignalALamALamAntiProtSep","AntiProton pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalALamALamAntiProtSep);

  fSignalALamALamPiPlusSep = new TH2F ("fSignalALamALamPiPlusSep","PiPlus pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalALamALamPiPlusSep);

  fSignalLamALamAntiProtPiMinusSep = new TH2F ("fSignalLamALamAntiProtPiMinusSep","Neg particle pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALamAntiProtPiMinusSep);

  fSignalLamALamProtPiPlusSep = new TH2F ("fSignalLamALamProtPiPlusSep","Pos pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALamProtPiPlusSep);
  
  fBkgLamLamProtSep = new TH2F ("fBkgLamLamProtSep","Proton pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamLamProtSep);

  fBkgLamLamPiMinusSep = new TH2F ("fBkgLamLamPiMinusSep","PiMinus pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamLamPiMinusSep);

  fBkgALamALamAntiProtSep = new TH2F ("fBkgALamALamAntiProtSep","AntiProton pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgALamALamAntiProtSep);

  fBkgALamALamPiPlusSep = new TH2F ("fBkgALamALamPiPlusSep","PiPlus pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgALamALamPiPlusSep);

  fBkgLamALamAntiProtPiMinusSep = new TH2F ("fBkgLamALamAntiProtPiMinusSep","Neg particle pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALamAntiProtPiMinusSep);

  fBkgLamALamProtPiPlusSep = new TH2F ("fBkgLamALamProtPiPlusSep","Pos pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALamProtPiPlusSep);



  //opposite charged pair separation
  fSignalLamLamPlusMinusSep = new TH2F ("fSignalLamLamPlusMinusSep","Proton Pion pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamLamPlusMinusSep);

  fSignalALamALamPlusMinusSep = new TH2F ("fSignalALamALamPlusMinusSep","Proton Pion pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalALamALamPlusMinusSep);

  fSignalLamALamProtSep = new TH2F ("fSignalLamALamProtSep","Proton pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALamProtSep);

  fSignalLamALamPionSep = new TH2F ("fSignalLamALamPionSep","Pion pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALamPionSep);


  fBkgLamLamPlusMinusSep = new TH2F ("fBkgLamLamPlusMinusSep","Proton Pion pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamLamPlusMinusSep);

  fBkgALamALamPlusMinusSep = new TH2F ("fBkgALamALamPlusMinusSep","Proton Pion pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgALamALamPlusMinusSep);

  fBkgLamALamProtSep = new TH2F ("fBkgLamALamProtSep","Proton pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALamProtSep);

  fBkgLamALamPionSep = new TH2F ("fBkgLamALamPionSep","Pion pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALamPionSep);




  //primary-vertex corrected pair separation

  fSignalLamLamProtSepCorrected = new TH2F ("fSignalLamLamProtSepCorrected","Proton pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamLamProtSepCorrected);

  fSignalLamLamPiMinusSepCorrected = new TH2F ("fSignalLamLamPiMinusSepCorrected","PiMinus pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamLamPiMinusSepCorrected);

  fSignalALamALamAntiProtSepCorrected = new TH2F ("fSignalALamALamAntiProtSepCorrected","AntiProton pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalALamALamAntiProtSepCorrected);

  fSignalALamALamPiPlusSepCorrected = new TH2F ("fSignalALamALamPiPlusSepCorrected","PiPlus pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalALamALamPiPlusSepCorrected);

  fSignalLamALamAntiProtPiMinusSepCorrected = new TH2F ("fSignalLamALamAntiProtPiMinusSepCorrected","Neg particle pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALamAntiProtPiMinusSepCorrected);

  fSignalLamALamProtPiPlusSepCorrected = new TH2F ("fSignalLamALamProtPiPlusSepCorrected","Pos pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALamProtPiPlusSepCorrected);

  fBkgLamLamProtSepCorrected = new TH2F ("fBkgLamLamProtSepCorrected","Proton pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamLamProtSepCorrected);

  fBkgLamLamPiMinusSepCorrected = new TH2F ("fBkgLamLamPiMinusSepCorrected","PiMinus pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamLamPiMinusSepCorrected);

  fBkgALamALamAntiProtSepCorrected = new TH2F ("fBkgALamALamAntiProtSepCorrected","AntiProton pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgALamALamAntiProtSepCorrected);

  fBkgALamALamPiPlusSepCorrected = new TH2F ("fBkgALamALamPiPlusSepCorrected","PiPlus pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgALamALamPiPlusSepCorrected);

  fBkgLamALamAntiProtPiMinusSepCorrected = new TH2F ("fBkgLamALamAntiProtPiMinusSepCorrected","Neg particle pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALamAntiProtPiMinusSepCorrected);

  fBkgLamALamProtPiPlusSepCorrected = new TH2F ("fBkgLamALamProtPiPlusSepCorrected","Pos pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALamProtPiPlusSepCorrected);

  //opposite charged pair separation corrected for primary vertex
  fSignalLamLamPlusMinusSepCorrected = new TH2F ("fSignalLamLamPlusMinusSepCorrected","Proton Pion pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamLamPlusMinusSepCorrected);

  fSignalALamALamPlusMinusSepCorrected = new TH2F ("fSignalALamALamPlusMinusSepCorrected","Proton Pion pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalALamALamPlusMinusSepCorrected);

  fSignalLamALamProtSepCorrected = new TH2F ("fSignalLamALamProtSepCorrected","Proton pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALamProtSepCorrected);

  fSignalLamALamPionSepCorrected = new TH2F ("fSignalLamALamPionSepCorrected","Pion pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fSignalLamALamPionSepCorrected);


  fBkgLamLamPlusMinusSepCorrected = new TH2F ("fBkgLamLamPlusMinusSepCorrected","Proton Pion pair sep for Lam-Lam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamLamPlusMinusSepCorrected);

  fBkgALamALamPlusMinusSepCorrected = new TH2F ("fBkgALamALamPlusMinusSepCorrected","Proton Pion pair sep for ALam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgALamALamPlusMinusSepCorrected);

  fBkgLamALamProtSepCorrected = new TH2F ("fBkgLamALamProtSepCorrected","Proton pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALamProtSepCorrected);

  fBkgLamALamPionSepCorrected = new TH2F ("fBkgLamALamPionSepCorrected","Pion pair sep for Lam-ALam", avgSepBins, 0., avgSepMaxValue, kStarBins, 0., maxKStar);
  fOutputList->Add(fBkgLamALamPionSepCorrected);


  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisV0Lam::Exec(Option_t *) 
{
  // Main loop
  // Called for each event

  //Make sure we are using the correct event triggers
  if(!IsCorrectEventTrigger()) return;
  //  cout<<"===========  Event # "<<fEventCount+1<<"  ==========="<<endl;
  fEventCount++;
  fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!fAOD) {Printf("ERROR: fAOD not available"); return;}
  AliAODVertex *primaryVertexAOD;
  double vertex[3]={0};
  int zBin=0;
  double zStep=2*10/double(zVertexBins), zStart=-10.;
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
  else if(centralityPercentile <= 5.) centralityBin=19;
  else if(centralityPercentile <= 10.) centralityBin=18;
  else if(centralityPercentile <= 15.) centralityBin=17;
  else if(centralityPercentile <= 20.) centralityBin=16;
  else if(centralityPercentile <= 25.) centralityBin=15;
  else if(centralityPercentile <= 30.) centralityBin=14;
  else if(centralityPercentile <= 35.) centralityBin=13;
  else if(centralityPercentile <= 40.) centralityBin=12;
  else if(centralityPercentile <= 45.) centralityBin=11;
  else if(centralityPercentile <= 50.) centralityBin=10;
  else if(centralityPercentile <= 55.) centralityBin=9;
  else if(centralityPercentile <= 60.) centralityBin=8;
  else if(centralityPercentile <= 65.) centralityBin=7;
  else if(centralityPercentile <= 70.) centralityBin=6;
  else if(centralityPercentile <= 75.) centralityBin=5;
  else if(centralityPercentile <= 80.) centralityBin=4;
  else if(centralityPercentile <= 85.) centralityBin=3;
  else if(centralityPercentile <= 90.) centralityBin=2;
  else if(centralityPercentile <= 95.) centralityBin=1;
  else if(centralityPercentile <= 100.) centralityBin=0;
  else {Printf("Skipping Peripheral Event"); return;}
   
  //Vertexing
  primaryVertexAOD = fAOD->GetPrimaryVertex();
  vertex[0]=primaryVertexAOD->GetX(); 
  vertex[1]=primaryVertexAOD->GetY(); 
  vertex[2]=primaryVertexAOD->GetZ();
  if(vertex[0]<10e-5 && vertex[1]<10e-5 &&  vertex[2]<10e-5) return;
  if(fabs(vertex[2]) > 10) return; // Z-Vertex Cut
  for(int i=0; i<zVertexBins; i++)
  {
    if((vertex[2] > zStart+i*zStep) && (vertex[2] < zStart+(i+1)*zStep))
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
  fCutProcessor->SetCentralityBin(centralityBin+1);
//////////////////////////////////////////////////////////////////  
  //v0 tester
////////////////////////////////////////////////////////////////
  
  int v0Count = 0;
  vector<int> lambdaCount(fNumberOfVariableCutValues,0);
  vector<int> antiLambdaCount(fNumberOfVariableCutValues,0);
  TH1F *mcTruthOriginHist = nullptr;
  if(fIsMCEvent) mcTruthOriginHist = CreateLambdaOriginHist(mcArray,numberOfLastHijingLabel); //Find MC truths of all the MC particles (before detector effects)
  TH1F *v0OriginHist = new TH1F("v0OriginHist", "Lambda Origins", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1);
  TH2F *v0PassedCutsOriginHist = new TH2F("v0PassedCutsOriginHist", "Lambda Origins", AliReconstructedV0::kOriginTypeMax+1, 0, AliReconstructedV0::kOriginTypeMax+1, fNumberOfVariableCutValues, -0.5, fNumberOfVariableCutValues -0.5);
  for(int i = 0; i < fAOD->GetNumberOfV0s(); i++)
  {
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
    if(v0->GetOnFlyStatus())				 continue;
    //Now look at daughter tracks
    AliAODTrack* daughterTrackPos = (AliAODTrack*)v0->GetDaughter(0);
    AliAODTrack* daughterTrackNeg = (AliAODTrack*)v0->GetDaughter(1);
    if(!daughterTrackPos) continue; //Daughter tracks must exist
    if(!daughterTrackNeg) continue;
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
    double v0MomentumTruth[3] = {0.,0.,0.};
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
    if (daughterTrackNeg->Pt() > 0.3){
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
    fEvt->fReconstructedV0[v0Count].v0Momentum[0]  = v0->Px();
    fEvt->fReconstructedV0[v0Count].v0Momentum[1]  = v0->Py();
    fEvt->fReconstructedV0[v0Count].v0Momentum[2]  = v0->Pz();
    for(int pIndex = 0; pIndex <3; pIndex++)
    {
      fEvt->fReconstructedV0[v0Count].v0MomentumTruth[pIndex]  = v0MomentumTruth[pIndex];
    }
    fEvt->fReconstructedV0[v0Count].v0Pt     = v0->Pt();
    fEvt->fReconstructedV0[v0Count].v0Eta    = v0->Eta();
    fEvt->fReconstructedV0[v0Count].v0Phi    = v0->Phi();
    fEvt->fReconstructedV0[v0Count].massLam  = v0->MassLambda();
    fEvt->fReconstructedV0[v0Count].massALam = v0->MassAntiLambda();
    fEvt->fReconstructedV0[v0Count].massLamDifference  = fabs(v0->MassLambda() - fPDGLambda);
    fEvt->fReconstructedV0[v0Count].massALamDifference = fabs(v0->MassAntiLambda() - fPDGLambda);
    fEvt->fReconstructedV0[v0Count].massK0  = v0->MassK0Short();
    fEvt->fReconstructedV0[v0Count].lorentzGammaLam = v0->ELambda()/fPDGLambda;
    fEvt->fReconstructedV0[v0Count].v0DCA   = v0->DcaV0ToPrimVertex();
    fEvt->fReconstructedV0[v0Count].decayLength = v0->DecayLength(primaryVertexAOD);
    fEvt->fReconstructedV0[v0Count].decayVertexPosition[0] = v0->DecayVertexV0X();
    fEvt->fReconstructedV0[v0Count].decayVertexPosition[1] = v0->DecayVertexV0Y();
    fEvt->fReconstructedV0[v0Count].decayVertexPosition[2] = v0->DecayVertexV0Z();
    fEvt->fReconstructedV0[v0Count].cosPointing = v0->CosPointingAngle(primaryVertexAOD);
    fEvt->fReconstructedV0[v0Count].hasProtonDaughter     = hasProtonDaughter;
    fEvt->fReconstructedV0[v0Count].hasAntiProtonDaughter = hasAntiProtonDaughter;
    fEvt->fReconstructedV0[v0Count].hasPiPlusDaughter     = hasPiPlusDaughter;
    fEvt->fReconstructedV0[v0Count].hasPiMinusDaughter    = hasPiMinusDaughter;
    fEvt->fReconstructedV0[v0Count].mcOriginType = mcV0Origin;
    //Save Daughter information
    fEvt->fReconstructedV0[v0Count].daughter1ID = v0->GetNegID();
    fEvt->fReconstructedV0[v0Count].daughter2ID = v0->GetPosID();
    fEvt->fReconstructedV0[v0Count].daughterPosMomentum[0] = v0->MomPosX();
    fEvt->fReconstructedV0[v0Count].daughterPosMomentum[1] = v0->MomPosY();
    fEvt->fReconstructedV0[v0Count].daughterPosMomentum[2] = v0->MomPosZ();
    fEvt->fReconstructedV0[v0Count].daughterPosProtonE = v0->EPosProton();
    fEvt->fReconstructedV0[v0Count].daughterPosPionE = v0->EPosPion();
    fEvt->fReconstructedV0[v0Count].daughterNegMomentum[0] = v0->MomNegX();
    fEvt->fReconstructedV0[v0Count].daughterNegMomentum[1] = v0->MomNegY();
    fEvt->fReconstructedV0[v0Count].daughterNegMomentum[2] = v0->MomNegZ();
    fEvt->fReconstructedV0[v0Count].daughterNegProtonE = v0->ENegProton();
    fEvt->fReconstructedV0[v0Count].daughterNegPionE = v0->ENegPion();
    fEvt->fReconstructedV0[v0Count].daughtersDCA = v0->DcaV0Daughters();
    fEvt->fReconstructedV0[v0Count].daughterPosDCAPrimaryVertex = v0->DcaPosToPrimVertex();
    fEvt->fReconstructedV0[v0Count].daughterNegDCAPrimaryVertex = v0->DcaNegToPrimVertex();
    GetGlobalPositionAtGlobalRadiiThroughTPC(daughterTrackPos, bfield, fEvt->fReconstructedV0[v0Count].daughterPosGlobalPositions);// used for merging cuts later
    GetGlobalPositionAtGlobalRadiiThroughTPC(daughterTrackNeg, bfield, fEvt->fReconstructedV0[v0Count].daughterNegGlobalPositions);
    daughterTrackPos->GetXYZ(fEvt->fReconstructedV0[v0Count].daughterPosPositionDCA);
    daughterTrackNeg->GetXYZ(fEvt->fReconstructedV0[v0Count].daughterNegPositionDCA);
    daughterTrackPos->GetPxPyPz(fEvt->fReconstructedV0[v0Count].daughterPosMomentumDCA);
    daughterTrackNeg->GetPxPyPz(fEvt->fReconstructedV0[v0Count].daughterNegMomentumDCA);
    daughterTrackPos->GetCovarianceXYZPxPyPz(fEvt->fReconstructedV0[v0Count].daughterPosCovariance);
    daughterTrackNeg->GetCovarianceXYZPxPyPz(fEvt->fReconstructedV0[v0Count].daughterNegCovariance);
    for(int coord = 0; coord <3; coord++){
      fEvt->fPrimaryVertex[coord] = vertex[coord];
      for(int location = 0; location <9; location++) {
	//find the track locations relative to the primary vertex location.
	fEvt->fReconstructedV0[v0Count].daughterPosCorrectedGlobalPositions[location][coord] = fEvt->fReconstructedV0[v0Count].daughterPosGlobalPositions[location][coord] - vertex[coord];
	fEvt->fReconstructedV0[v0Count].daughterNegCorrectedGlobalPositions[location][coord] = fEvt->fReconstructedV0[v0Count].daughterNegGlobalPositions[location][coord] - vertex[coord];
      }
    }

    //Now analyze and histogram the V0
    fCutProcessor->CheckIfV0PassesCuts(& fEvt->fReconstructedV0[v0Count]);
    fCutProcessor->DoV0Histogramming(& fEvt->fReconstructedV0[v0Count]);
    FillReconstructedV0MCOrigin(& fEvt->fReconstructedV0[v0Count], v0PassedCutsOriginHist);
    AddV0ToMultiplicityCounts(& fEvt->fReconstructedV0[v0Count], lambdaCount, antiLambdaCount);
    FillTPCSignalHists(& fEvt->fReconstructedV0[v0Count], daughterTrackPos->P(), daughterTrackPos->GetTPCsignal(), daughterTrackNeg->P(), daughterTrackNeg->GetTPCsignal());
    if(fIsMCEvent) CheckForFakeV0s(& fEvt->fReconstructedV0[v0Count], fMCFakeParticleIdentity, fMCOtherV0Identity, mcV0Origin);
    
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
    mcTruthOriginHist = nullptr;
  }
  if(v0OriginHist){
    delete v0OriginHist;
    v0OriginHist = nullptr;
  }
  if(v0PassedCutsOriginHist){
    delete v0PassedCutsOriginHist;
    v0PassedCutsOriginHist = nullptr;
  }

  
  DoV0JudgmentCuts(fEvt, v0Count);
  HistogramEventMultiplicities(lambdaCount, antiLambdaCount, centralityBin);
  fTotalLambda += lambdaCount[fDefaultVariableCutIndex];
  fTotalAntiLambda += antiLambdaCount[fDefaultVariableCutIndex];
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
  cout<<"Total Lambdas found:\t"<<fTotalLambda<<"."<<endl
      <<"Total AntiLambdas found:\t"<<fTotalAntiLambda<<"."<<endl
      <<"Done"<<endl;
}






//________________________________________________________________________
void AliAnalysisV0Lam::GetGlobalPositionAtGlobalRadiiThroughTPC(const AliAODTrack *track, const Float_t bfield, Float_t globalPositionsAtRadii[9][3])
{
  // Gets the global position of the track at nine different radii in the TPC
  // track is the track you want to propagate
  // bfield is the magnetic field of your event
  //globalPositionsAtRadii is the array of global positions in the radii and xyz
  // Initialize the array to something indicating there was no propagation
  for(Int_t i=0;i<9;i++){
    for(Int_t j=0;j<3;j++){
      globalPositionsAtRadii[i][j]=-9999.;
    }
  }
  // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp; 
  etp.CopyFromVTrack(track);
  //printf("\nAfter CopyFromVTrack\n");
  //etp.Print();
  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.}; 
  // Counter for which radius we want
  Int_t iR=0;
  // The radii at which we get the global positions
  // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
  Float_t Rwanted[9]={85.,105.,125.,145.,165.,185.,205.,225.,245.};
  // The global radius we are at
  Float_t globalRadius=0;
  // Propagation is done in local x of the track
  for (Float_t x = etp.GetX();x<247.;x+=1.){ // GetX returns local coordinates
    // Starts at the tracks fX and goes outwards. x = 245 is the outer radial
    // limit of the TPC when the track is straight, i.e. has inifinite pt
    // and doesn't get bent. If the track's momentum is smaller than infinite,
    // it will develop a y-component, which adds to the global radius
    // Stop if the propagation was not succesful. This can happen for low pt
    // tracks that don't reach outer radii
    if(!etp.PropagateTo(x,bfield))break;
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates
    globalRadius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]); //Idea to speed up: compare squared radii
    // Roughly reached the radius we want
    if(globalRadius > Rwanted[iR]){
      // Bigger loop has bad precision, we're nearly one centimeter too far, go back in small steps.
      while (globalRadius>Rwanted[iR]){
        x-=.1;
        //      printf("propagating to x %5.2f\n",x);
        if(!etp.PropagateTo(x,bfield))break;
        etp.GetXYZ(xyz); // GetXYZ returns global coordinates
        globalRadius = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]); //Idea to speed up: compare squared radii
      }
      //printf("At Radius:%05.2f (local x %5.2f). Setting position to x %4.1f y %4.1f z %4.1f\n",globalRadius,x,xyz[0],xyz[1],xyz[2]);
      globalPositionsAtRadii[iR][0]=xyz[0];
      globalPositionsAtRadii[iR][1]=xyz[1];
      globalPositionsAtRadii[iR][2]=xyz[2];
      // Indicate we want the next radius   
      iR+=1;
    }
    if(iR>=8){
      // TPC edge reached
      return;
    }
  }
}



//________________________________________________________________________
double AliAnalysisV0Lam::GetAverageSeparation(Float_t globalPositions1st[9][3], Float_t globalPositions2nd[9][3])
{
  //Compute the separation of two daughter tracks, averaged over 9 different positions
  double avgSeparation = 0;
  for(int RadiusNumber = 0; RadiusNumber <9; RadiusNumber++){
    double sumsquare=0;
    for(int Component = 0; Component <3; Component++){
      sumsquare+= pow(globalPositions1st[RadiusNumber][Component]-globalPositions2nd[RadiusNumber][Component],2);
    }
    avgSeparation+=sqrt(sumsquare);
  }
  return avgSeparation/9.; 
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
  for (int cutIndex = 0; cutIndex < fNumberOfVariableCutValues; cutIndex++){
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
    AliAODMCParticle *mcDaughter1 = (AliAODMCParticle*)mcArray->At(mcParticle->GetDaughter(0));
    AliAODMCParticle *mcDaughter2 = (AliAODMCParticle*)mcArray->At(mcParticle->GetDaughter(1));
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
  for(int i = 0; i < fNumberOfVariableCutValues; i++){
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
void AliAnalysisV0Lam::GetMCParticleMomentumTruth(double *v0Momentum, AliAODv0 *v0, TClonesArray *mcArray)
{
  // Get the MC truth of the 3D momentum of a V0.  That info is copied into
  // v0Momentum
  int v0Id = GetV0MCParticleID(v0,mcArray);
  if (v0Id > 0){
    AliAODMCParticle* mcV0 = (AliAODMCParticle*)mcArray->At(v0Id);
    if(!mcV0->PxPyPz(v0Momentum)) cout<<"Err copying momentum truth"<<endl;
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
  AliAODTrack* daughterTrackPos = (AliAODTrack*)v0->GetDaughter(0);
  AliAODTrack* daughterTrackNeg = (AliAODTrack*)v0->GetDaughter(1);
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
  for(int i = 0; i < fNumberOfVariableCutValues; i++){
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

//________________________________________________________________________
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
  for(int i = 0; i < fNumberOfVariableCutValues; i++){
    if(v0->isLamCenter[i]) lambdaCount[i]++;
    if(v0->isALamCenter[i]) antiLambdaCount[i]++;
  }
}

//________________________________________________________________________
void AliAnalysisV0Lam::HistogramEventMultiplicities(vector<int> & lambdaCount, vector<int> & antiLambdaCount, int centralityBin)
{
  //Add the event yields to the output yield histograms
  for(int i = 0; i < fNumberOfVariableCutValues; i++){
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
  if(v0->isLamCenter[fDefaultVariableCutIndex])
  {
    fTPCVsPPosLam->Fill(posDaughterP,posDaughterTPCSignal);
    fTPCVsPNegLam->Fill(negDaughterP,negDaughterTPCSignal);
  }
  if(v0->isALamCenter[fDefaultVariableCutIndex])
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
  if(v0->isLamCenter[fDefaultVariableCutIndex] 
     || v0->isALamCenter[fDefaultVariableCutIndex])
  {
    if(AliReconstructedV0::kFake == mcV0Origin){
      //These V0s are fake
      if(v0->isLamCenter[fDefaultVariableCutIndex]) mcFakeParticleIdentity->Fill(0);
      if(v0->isALamCenter[fDefaultVariableCutIndex]) mcFakeParticleIdentity->Fill(1);
    }
    if(AliReconstructedV0::kUnassigned == mcV0Origin){
      // These V0s correspond to actual MC particles, but they aren't
      // Lambdas, Antilambdas, or K0s.
      if(v0->isLamCenter[fDefaultVariableCutIndex]) mcOtherV0Identity->Fill(0);
      if(v0->isALamCenter[fDefaultVariableCutIndex]) mcOtherV0Identity->Fill(1);
    }
  }
}

//________________________________________________________________________
double AliAnalysisV0Lam::CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2)
{
  // Calculate k* for any pair of particles, regardless of whether the
  // particles have the same mass.
  double kstar = 0.;
  double e1 = 0.;
  double e2 = 0.;
  for(int i = 0; i < 3; i++){
    kstar -= pow(momentum1[i]-momentum2[i],2);
    e1 += pow(momentum1[i],2);
    e2 += pow(momentum2[i],2);
  }
  e1 += pow(mass1,2);
  e1 = sqrt(e1);
  e2 += pow(mass2,2);
  e2 = sqrt(e2);
  
  kstar += pow(e1-e2,2);

  double totalMomentumSquared = 0;
  for(int i = 0; i < 3; i++){
    totalMomentumSquared -= pow(momentum1[i]+momentum2[i],2);
  }
  totalMomentumSquared += pow(e1+e2,2);
  kstar -= pow((pow(mass1,2)-pow(mass2,2)),2)/totalMomentumSquared;

  kstar *= -1.;
  kstar = sqrt(kstar); //At this point, we've actually calculated Qinv
  kstar *= 0.5; // kstar is 0.5*Qinv
  return kstar;
}

//________________________________________________________________________
void AliAnalysisV0Lam::BinMomentumSmearing(double *v01MomentumRecon, double *v01MomentumTruth, double *v02MomentumRecon, double *v02MomentumTruth, int pairType)
{
  // Used in MC studies to determine the amount of pair-wise relative-
  // momentum smearing.  For a given pair of V0s, fills output histograms
  // with reconstructed k* vs true k*.  Separate histograms are filled
  // for each pair type
  double v01PMag = sqrt( pow(v01MomentumTruth[0],2) + pow(v01MomentumTruth[1],2) + pow(v01MomentumTruth[2],2) );
  if(v01PMag < 0.0001) return;  //Not a real V0
  double v02PMag = sqrt( pow(v02MomentumTruth[0],2) + pow(v02MomentumTruth[1],2) + pow(v02MomentumTruth[2],2) );
  if(v02PMag < 0.0001) return;  //Not a real V0
  double reconKstar = CalculateKstar(v01MomentumRecon,v02MomentumRecon, fPDGLambda, fPDGLambda);
  double truthKstar = CalculateKstar(v01MomentumTruth,v02MomentumTruth, fPDGLambda, fPDGLambda);
  if(0 == pairType){ //Lambda-Lambda
    fHistPsmearingKreconVsKtruthLL->Fill(reconKstar,truthKstar);
    fHistPsmearingKreconMinusKtruthLL->Fill(reconKstar-truthKstar);
  }
  if(1 == pairType){ //Antilambda-Antilambda
    fHistPsmearingKreconVsKtruthAA->Fill(reconKstar,truthKstar);
    fHistPsmearingKreconMinusKtruthAA->Fill(reconKstar-truthKstar);
  }
  if(2 == pairType){ //Lambda-Antilambda
    fHistPsmearingKreconVsKtruthLA->Fill(reconKstar,truthKstar);
    fHistPsmearingKreconMinusKtruthLA->Fill(reconKstar-truthKstar);
  }
}
  
//________________________________________________________________________
void AliAnalysisV0Lam::DoPairStudies(const AliAnalysisV0LamEvent *event, const int centralityBin)
{
  //Loop over same- and mixed-event pairs and bin correlation function
  //numerators and denominators.

  //Default cut values used for the average separation cuts
  double avgSepIdenticalProtonCut = 3.; // used for same-charge prot avg sep
  double avgSepIdenticalPionCut = 4.; // used for same-charge pi-pi
  double avgSepNonIdenticalCut = 3.5; //used for same-charge pion-proton avg sep
  double &variableAvgSepValue = avgSepIdenticalPionCut; //the identical
  //pion cut will be the variable cut in this run.  This needs to be changed
  //by hand if one wants to study a different cut
  
  // These values are used for the average separation cuts of the
  // daughter tracks. When the average separation cut is being studied,
  // we can loop over many different cut values to see how those different
  // cuts affect the correlation function results
  double avgSepCutArray[12] = {0.,1.,2.,2.5,3.,3.5,4.,4.5,5.,6.,7.,10.};

  int numberOfAvgSepCutIndices; //this is the size of the avgSepCutArray
  if(fIsUsingVariableAvgSepCut) numberOfAvgSepCutIndices = sizeof(avgSepCutArray) / sizeof (avgSepCutArray[0]);
  else numberOfAvgSepCutIndices = 1;

  int defaultVariableAvgSepCutIndex = 0;
  if(fIsUsingVariableAvgSepCut) {
    for(int i = 0; i < numberOfAvgSepCutIndices; i++){
      // If variable avg sep cuts are being used, this find the default
      // value for the avg sep cut currently being studied.
      if (fabs(variableAvgSepValue - avgSepCutArray[i]) < 0.1) defaultVariableAvgSepCutIndex = i;
    }
  }
  for(int cutIndex = 0; cutIndex < fNumberOfVariableCutValues; cutIndex++)
  { // Start looping over all variable cut values
    for(int i=0; i < event->fNumberCandidateV0; i++) 
    { //Start looping over reconstructed V0s in this event
      bool center1Lam  = event->fReconstructedV0[i].isLamCenter[cutIndex];
      bool center1ALam = event->fReconstructedV0[i].isALamCenter[cutIndex];
      // Disregard V0 if it wasn't reconstructed as a center (anti)Lambda
      if(!(center1Lam || center1ALam)) continue;
      // Disregard V0 if it was removed via the judgment cuts
      if(event->fReconstructedV0[i].isDeemedUnworthy[cutIndex]) continue;
      for(int eventNumber=0; eventNumber<nEventsToMix+1; eventNumber++)
      { // Event buffer loop: eventNumber=0 is the current event, all other eventNumbers are past events
	int startBin=0;
	// For same event pairs, start 2nd V0 loop at i+1 V0 to avoid
	// double counting
	if(eventNumber==0) startBin=i+1;
	for(int j=startBin; j<(event+eventNumber)->fNumberCandidateV0; j++) 
	{ // Second V0 loop (from past or current event)
	  if(eventNumber==0)  
	  { // Don't make pairs of V0s if they shared daughter tracks.
	    // This is redundant if the judgment cut is already employed
	    if(event->fReconstructedV0[i].daughter1ID 
	       == (event+eventNumber)->fReconstructedV0[j].daughter1ID) continue;
	    if(event->fReconstructedV0[i].daughter1ID 
	       == (event+eventNumber)->fReconstructedV0[j].daughter2ID) continue;
	    if(event->fReconstructedV0[i].daughter2ID 
	       == (event+eventNumber)->fReconstructedV0[j].daughter1ID) continue;
	    if(event->fReconstructedV0[i].daughter2ID 
	       == (event+eventNumber)->fReconstructedV0[j].daughter2ID) continue;
	  }
	  //Disregard second V0 if it was removed via judgment cuts
	  if((event+eventNumber)->fReconstructedV0[j].isDeemedUnworthy[cutIndex]) continue;
	  // A central V0 has a mass that falls within the accepted inv
	  // mass range.  Only make pairs with central V0s
	  bool center2Lam  = (event+eventNumber)->fReconstructedV0[j].isLamCenter[cutIndex];
	  bool center2ALam = (event+eventNumber)->fReconstructedV0[j].isALamCenter[cutIndex];
	  if(!(center2Lam || center2ALam)) continue;
	  // Now we calculate a bunch of values that are used later during
	  // histogramming.
	  double pairKt = pow(event->fReconstructedV0[i].v0Momentum[0] + (event+eventNumber)->fReconstructedV0[j].v0Momentum[0],2.);
	  pairKt+= pow(event->fReconstructedV0[i].v0Momentum[1] + (event+eventNumber)->fReconstructedV0[j].v0Momentum[1],2.);
	  pairKt = sqrt(pairKt)/2.;
	  //Calculate k* for V0s and daughters using different mass assumptions
	  double pairKstarLam = CalculateKstar(event->fReconstructedV0[i].v0Momentum, (event+eventNumber)->fReconstructedV0[j].v0Momentum, fPDGLambda,fPDGLambda);
	  double pairKstarProtPlus = CalculateKstar(event->fReconstructedV0[i].daughterPosMomentum,(event+eventNumber)->fReconstructedV0[j].daughterPosMomentum, fPDGProton,fPDGProton);
	  double pairKstarProtMinus = CalculateKstar(event->fReconstructedV0[i].daughterNegMomentum,(event+eventNumber)->fReconstructedV0[j].daughterNegMomentum, fPDGProton,fPDGProton);
	  double pairKstarPiPlus = CalculateKstar(event->fReconstructedV0[i].daughterPosMomentum,(event+eventNumber)->fReconstructedV0[j].daughterPosMomentum, fPDGPion,fPDGPion);
	  double pairKstarPiMinus = CalculateKstar(event->fReconstructedV0[i].daughterNegMomentum,(event+eventNumber)->fReconstructedV0[j].daughterNegMomentum, fPDGPion,fPDGPion);
	  //used for lambda-antilambda daughter kstar
	  double pairKstarProtPlusPiPlus = 0; 
	  double pairKstarProtMinusPiMinus = 0;
	  double pairKstarProtPlusProtMinus = 0;
	  double pairKstarPiPlusPiMinus = 0;

	  // Need to be careful when calculating the k* of daughter tracks
	  // of non-identical V0s
	  if(center1Lam && center2ALam){
	    pairKstarProtPlusPiPlus = CalculateKstar(event->fReconstructedV0[i].daughterPosMomentum,(event+eventNumber)->fReconstructedV0[j].daughterPosMomentum, fPDGProton,fPDGPion);
	    pairKstarProtMinusPiMinus = CalculateKstar(event->fReconstructedV0[i].daughterNegMomentum,(event+eventNumber)->fReconstructedV0[j].daughterNegMomentum, fPDGPion,fPDGProton);
	    pairKstarProtPlusProtMinus = CalculateKstar(event->fReconstructedV0[i].daughterPosMomentum,(event+eventNumber)->fReconstructedV0[j].daughterNegMomentum, fPDGProton,fPDGProton); 
	    pairKstarPiPlusPiMinus = CalculateKstar(event->fReconstructedV0[i].daughterNegMomentum,(event+eventNumber)->fReconstructedV0[j].daughterPosMomentum, fPDGPion,fPDGPion); 
	  }
	  else if(center1ALam && center2Lam){
	    pairKstarProtPlusPiPlus = CalculateKstar(event->fReconstructedV0[i].daughterPosMomentum,(event+eventNumber)->fReconstructedV0[j].daughterPosMomentum, fPDGPion,fPDGProton);
	    pairKstarProtMinusPiMinus = CalculateKstar(event->fReconstructedV0[i].daughterNegMomentum,(event+eventNumber)->fReconstructedV0[j].daughterNegMomentum, fPDGProton,fPDGPion);
	    pairKstarProtPlusProtMinus = CalculateKstar(event->fReconstructedV0[i].daughterNegMomentum,(event+eventNumber)->fReconstructedV0[j].daughterPosMomentum, fPDGProton,fPDGProton);
	    pairKstarPiPlusPiMinus = CalculateKstar(event->fReconstructedV0[i].daughterPosMomentum,(event+eventNumber)->fReconstructedV0[j].daughterNegMomentum, fPDGPion,fPDGPion);
	  }
	  double pairKstarProtPlusPiMinus1 = CalculateKstar(event->fReconstructedV0[i].daughterPosMomentum,(event+eventNumber)->fReconstructedV0[j].daughterNegMomentum, fPDGProton,fPDGPion); // only relevant for LamLam
	  double pairKstarProtPlusPiMinus2 = CalculateKstar(event->fReconstructedV0[i].daughterNegMomentum,(event+eventNumber)->fReconstructedV0[j].daughterPosMomentum, fPDGPion,fPDGProton); // only relevant for LamLam

	  double pairKstarProtMinusPiPlus1 = CalculateKstar(event->fReconstructedV0[i].daughterNegMomentum,(event+eventNumber)->fReconstructedV0[j].daughterPosMomentum, fPDGProton,fPDGPion); // only relevant for ALamALam
	  double pairKstarProtMinusPiPlus2 = CalculateKstar(event->fReconstructedV0[i].daughterPosMomentum,(event+eventNumber)->fReconstructedV0[j].daughterNegMomentum, fPDGPion, fPDGProton); // only relevant for ALamALam
	  //Now find the average separation distance between daughter pairs.  Used to make a merging/splitting cut.
	  double avgSepPos = GetAverageSeparation(event->fReconstructedV0[i].daughterPosGlobalPositions, (event+eventNumber)->fReconstructedV0[j].daughterPosGlobalPositions);
	  double avgSepNeg = GetAverageSeparation(event->fReconstructedV0[i].daughterNegGlobalPositions, (event+eventNumber)->fReconstructedV0[j].daughterNegGlobalPositions);
	  double avgSepNegPos = GetAverageSeparation(event->fReconstructedV0[i].daughterNegGlobalPositions, (event+eventNumber)->fReconstructedV0[j].daughterPosGlobalPositions);
	  double avgSepPosNeg = GetAverageSeparation(event->fReconstructedV0[i].daughterPosGlobalPositions, (event+eventNumber)->fReconstructedV0[j].daughterNegGlobalPositions);
	  double correctedAvgSepPos = GetAverageSeparation(event->fReconstructedV0[i].daughterPosCorrectedGlobalPositions, (event+eventNumber)->fReconstructedV0[j].daughterPosCorrectedGlobalPositions);
	  double correctedAvgSepNeg = GetAverageSeparation(event->fReconstructedV0[i].daughterNegCorrectedGlobalPositions, (event+eventNumber)->fReconstructedV0[j].daughterNegCorrectedGlobalPositions);
	  double correctedAvgSepNegPos = GetAverageSeparation(event->fReconstructedV0[i].daughterNegCorrectedGlobalPositions, (event+eventNumber)->fReconstructedV0[j].daughterPosCorrectedGlobalPositions);
	  double correctedAvgSepPosNeg = GetAverageSeparation(event->fReconstructedV0[i].daughterPosCorrectedGlobalPositions, (event+eventNumber)->fReconstructedV0[j].daughterNegCorrectedGlobalPositions);

	  //Now we get to the actual pair histogramming
	  if(eventNumber==0) //Same event pair histogramming
	  {
	    //We do separate binning for each pair type
	    if(center1Lam && center2Lam){
	      // Some histograms we only fill when default variable cut
	      // values have been used
	      if(cutIndex == fDefaultVariableCutIndex){
		//same sign tracks
		fSignalLamLamProtSep->Fill(avgSepPos, pairKstarProtPlus);
		fSignalLamLamPiMinusSep->Fill(avgSepNeg, pairKstarPiMinus);
		fSignalLamLamProtSepCorrected->Fill(correctedAvgSepPos, pairKstarProtPlus);
		fSignalLamLamPiMinusSepCorrected->Fill(correctedAvgSepNeg, pairKstarPiMinus);
		//opposite sign tracks
		fSignalLamLamPlusMinusSep->Fill(avgSepPosNeg,pairKstarProtPlusPiMinus1);
		fSignalLamLamPlusMinusSep->Fill(avgSepNegPos,pairKstarProtPlusPiMinus2);
		fSignalLamLamPlusMinusSepCorrected->Fill(correctedAvgSepPosNeg,pairKstarProtPlusPiMinus1);
		fSignalLamLamPlusMinusSepCorrected->Fill(correctedAvgSepNegPos,pairKstarProtPlusPiMinus2);
	      }
	      for(int sepCutIndex = 0; sepCutIndex < numberOfAvgSepCutIndices; sepCutIndex++){ //looping over different avg sep cut values (if applicable)
		if(fIsUsingVariableAvgSepCut) variableAvgSepValue = avgSepCutArray[sepCutIndex];
		if((avgSepIdenticalProtonCut < correctedAvgSepPos)
		   && (avgSepIdenticalPionCut < correctedAvgSepNeg))
		{
		  fSignalLamLam->Fill(sepCutIndex, centralityBin+1, pairKstarLam);
		  if(defaultVariableAvgSepCutIndex == sepCutIndex){
		    //This implementation doesn't work properly
		    fKtLamLamSig->Fill(centralityBin+1,pairKt,pairKstarLam);
		  }
		}
	      }
	    }
	    if(center1ALam && center2ALam){
	      if(cutIndex == fDefaultVariableCutIndex){
		//same sign tracks
		fSignalALamALamAntiProtSep->Fill(avgSepNeg, pairKstarProtMinus);
		fSignalALamALamPiPlusSep->Fill(avgSepPos, pairKstarPiPlus);
		fSignalALamALamAntiProtSepCorrected->Fill(correctedAvgSepNeg, pairKstarProtMinus);
		fSignalALamALamPiPlusSepCorrected->Fill(correctedAvgSepPos, pairKstarPiPlus);

		//opposite sign tracks
		fSignalALamALamPlusMinusSep->Fill(avgSepPosNeg,pairKstarProtMinusPiPlus1);
		fSignalALamALamPlusMinusSep->Fill(avgSepNegPos,pairKstarProtMinusPiPlus2);
		fSignalALamALamPlusMinusSepCorrected->Fill(correctedAvgSepPosNeg,pairKstarProtMinusPiPlus1);
		fSignalALamALamPlusMinusSepCorrected->Fill(correctedAvgSepNegPos,pairKstarProtMinusPiPlus2);
	      }
	      for(int sepCutIndex = 0; sepCutIndex < numberOfAvgSepCutIndices; sepCutIndex++){
		if(fIsUsingVariableAvgSepCut) variableAvgSepValue = avgSepCutArray[sepCutIndex];
		if((avgSepIdenticalPionCut < correctedAvgSepPos)
		   && (avgSepIdenticalProtonCut < correctedAvgSepNeg))
		{
		  fSignalALamALam->Fill(sepCutIndex, centralityBin+1, pairKstarLam);
		  if(defaultVariableAvgSepCutIndex == sepCutIndex){
		    //This implementation doesn't work properly
		    fKtALamALamSig->Fill(centralityBin+1,pairKt,pairKstarLam);
		  }
		}
	      }
	    }
	    if((center1Lam && center2ALam) || (center1ALam && center2Lam)){
	      if(cutIndex == fDefaultVariableCutIndex){
		fSignalLamALamProtPiPlusSep->Fill(avgSepPos, pairKstarProtPlusPiPlus);
		fSignalLamALamAntiProtPiMinusSep->Fill(avgSepNeg, pairKstarProtMinusPiMinus);
		fSignalLamALamProtPiPlusSepCorrected->Fill(correctedAvgSepPos, pairKstarProtPlusPiPlus);
		fSignalLamALamAntiProtPiMinusSepCorrected->Fill(correctedAvgSepNeg, pairKstarProtMinusPiMinus);
		//opposite charge tracks
		if(center1Lam)
		{
		  fSignalLamALamProtSep->Fill(avgSepPosNeg, pairKstarProtPlusProtMinus);
		  fSignalLamALamPionSep->Fill(avgSepNegPos, pairKstarPiPlusPiMinus);
		  fSignalLamALamProtSepCorrected->Fill(correctedAvgSepPosNeg, pairKstarProtPlusProtMinus);
		  fSignalLamALamPionSepCorrected->Fill(correctedAvgSepNegPos, pairKstarPiPlusPiMinus);
		}
		else
		{
		  fSignalLamALamProtSep->Fill(avgSepNegPos, pairKstarProtPlusProtMinus);
		  fSignalLamALamPionSep->Fill(avgSepPosNeg, pairKstarPiPlusPiMinus);
		  fSignalLamALamProtSepCorrected->Fill(correctedAvgSepNegPos, pairKstarProtPlusProtMinus);
		  fSignalLamALamPionSepCorrected->Fill(correctedAvgSepPosNeg, pairKstarPiPlusPiMinus);
		}
	      }
	      for(int sepCutIndex = 0; sepCutIndex < numberOfAvgSepCutIndices; sepCutIndex++){
		if(fIsUsingVariableAvgSepCut) variableAvgSepValue = avgSepCutArray[sepCutIndex];
		if((avgSepNonIdenticalCut < correctedAvgSepPos)
		   && (avgSepNonIdenticalCut < correctedAvgSepNeg))
		{
		  fSignalLamALam->Fill(sepCutIndex, centralityBin+1, pairKstarLam);
		  if(defaultVariableAvgSepCutIndex == sepCutIndex){
		   //This implementation doesn't work properly
		    fKtLamALamSig->Fill(centralityBin+1,pairKt,pairKstarLam);
		  }
		}
	      }
	    }
	  } //end same event pair histogramming
	  else //Mixed event pair histogramming
	  {
	    if(center1Lam && center2Lam){
	      if(cutIndex == fDefaultVariableCutIndex){
		fBkgLamLamProtSep->Fill(avgSepPos, pairKstarProtPlus);
		fBkgLamLamPiMinusSep->Fill(avgSepNeg, pairKstarPiMinus);
		fBkgLamLamProtSepCorrected->Fill(correctedAvgSepPos, pairKstarProtPlus);
		fBkgLamLamPiMinusSepCorrected->Fill(correctedAvgSepNeg, pairKstarPiMinus);
		//opposite sign tracks
		fBkgLamLamPlusMinusSep->Fill(avgSepPosNeg,pairKstarProtPlusPiMinus1);
		fBkgLamLamPlusMinusSep->Fill(avgSepNegPos,pairKstarProtPlusPiMinus2);
		fBkgLamLamPlusMinusSepCorrected->Fill(correctedAvgSepPosNeg,pairKstarProtPlusPiMinus1);
		fBkgLamLamPlusMinusSepCorrected->Fill(correctedAvgSepNegPos,pairKstarProtPlusPiMinus2);
		if(fIsMCEvent)
		{ //collect momentum smearing information
		  int pairType = 0;
		  BinMomentumSmearing(event->fReconstructedV0[i].v0Momentum, event->fReconstructedV0[i].v0MomentumTruth, (event+eventNumber)->fReconstructedV0[j].v0Momentum, (event+eventNumber)->fReconstructedV0[j].v0MomentumTruth,pairType);
		}
	      }
	      for(int sepCutIndex = 0; sepCutIndex < numberOfAvgSepCutIndices; sepCutIndex++){
		if(fIsUsingVariableAvgSepCut) variableAvgSepValue = avgSepCutArray[sepCutIndex];
		if((avgSepIdenticalProtonCut < correctedAvgSepPos)
		   && (avgSepIdenticalPionCut < correctedAvgSepNeg))
		{
		  fBkgLamLam->Fill(sepCutIndex, centralityBin+1, pairKstarLam);
		  if(defaultVariableAvgSepCutIndex == sepCutIndex){
		    //This implementation doesn't work properly
		    fKtLamLamBkg->Fill(centralityBin+1,pairKt,pairKstarLam);
		  }
		}
	      }
	    }
	    if(center1ALam && center2ALam){
	      if(cutIndex == fDefaultVariableCutIndex){
		fBkgALamALamAntiProtSep->Fill(avgSepNeg, pairKstarProtMinus);
		fBkgALamALamPiPlusSep->Fill(avgSepPos, pairKstarPiPlus);
		fBkgALamALamAntiProtSepCorrected->Fill(correctedAvgSepNeg, pairKstarProtMinus);
		fBkgALamALamPiPlusSepCorrected->Fill(correctedAvgSepPos, pairKstarPiPlus);
		//opposite sign tracks
		fBkgALamALamPlusMinusSep->Fill(avgSepPosNeg,pairKstarProtMinusPiPlus1);
		fBkgALamALamPlusMinusSep->Fill(avgSepNegPos,pairKstarProtMinusPiPlus2);
		fBkgALamALamPlusMinusSepCorrected->Fill(correctedAvgSepPosNeg,pairKstarProtMinusPiPlus1);
		fBkgALamALamPlusMinusSepCorrected->Fill(correctedAvgSepNegPos,pairKstarProtMinusPiPlus2);
		if(fIsMCEvent)
		{ //collect momentum smearing information
		  int pairType = 1;
		  BinMomentumSmearing(event->fReconstructedV0[i].v0Momentum, event->fReconstructedV0[i].v0MomentumTruth, (event+eventNumber)->fReconstructedV0[j].v0Momentum, (event+eventNumber)->fReconstructedV0[j].v0MomentumTruth,pairType);
		}
	      }
	      for(int sepCutIndex = 0; sepCutIndex < numberOfAvgSepCutIndices; sepCutIndex++){
		if(fIsUsingVariableAvgSepCut) variableAvgSepValue = avgSepCutArray[sepCutIndex];
		if((avgSepIdenticalPionCut < correctedAvgSepPos)
		   && (avgSepIdenticalProtonCut < correctedAvgSepNeg))
		{
		  fBkgALamALam->Fill(sepCutIndex, centralityBin+1, pairKstarLam);
		  if(defaultVariableAvgSepCutIndex == sepCutIndex){
		    //This implementation doesn't work properly
		    fKtALamALamBkg->Fill(centralityBin+1,pairKt,pairKstarLam);
		  }
		}
	      }
	    }
	    if((center1Lam && center2ALam) || (center1ALam && center2Lam)){
	      if(cutIndex == fDefaultVariableCutIndex){
		fBkgLamALamProtPiPlusSep->Fill(avgSepPos, pairKstarProtPlusPiPlus);
		fBkgLamALamAntiProtPiMinusSep->Fill(avgSepNeg, pairKstarProtMinusPiMinus);
		fBkgLamALamProtPiPlusSepCorrected->Fill(correctedAvgSepPos, pairKstarProtPlusPiPlus);
		fBkgLamALamAntiProtPiMinusSepCorrected->Fill(correctedAvgSepNeg, pairKstarProtMinusPiMinus);
		//opposite charge tracks
		if(center1Lam)
		{
		  fBkgLamALamProtSep->Fill(avgSepPosNeg, pairKstarProtPlusProtMinus);
		  fBkgLamALamPionSep->Fill(avgSepNegPos, pairKstarPiPlusPiMinus);
		  fBkgLamALamProtSepCorrected->Fill(correctedAvgSepPosNeg, pairKstarProtPlusProtMinus);
		  fBkgLamALamPionSepCorrected->Fill(correctedAvgSepNegPos, pairKstarPiPlusPiMinus);
		}
		else
		{
		  fBkgLamALamProtSep->Fill(avgSepNegPos, pairKstarProtPlusProtMinus);
		  fBkgLamALamPionSep->Fill(avgSepPosNeg, pairKstarPiPlusPiMinus);
		  fBkgLamALamProtSepCorrected->Fill(correctedAvgSepNegPos, pairKstarProtPlusProtMinus);
		  fBkgLamALamPionSepCorrected->Fill(correctedAvgSepPosNeg, pairKstarPiPlusPiMinus);
		}
		if(fIsMCEvent)
		{ //collect momentum smearing information
		  int pairType = 2;
		  BinMomentumSmearing(event->fReconstructedV0[i].v0Momentum, event->fReconstructedV0[i].v0MomentumTruth, (event+eventNumber)->fReconstructedV0[j].v0Momentum, (event+eventNumber)->fReconstructedV0[j].v0MomentumTruth,pairType);
		}
	      }
	      for(int sepCutIndex = 0; sepCutIndex < numberOfAvgSepCutIndices; sepCutIndex++){
		if(fIsUsingVariableAvgSepCut) variableAvgSepValue = avgSepCutArray[sepCutIndex];
		if((avgSepNonIdenticalCut < correctedAvgSepPos)
		   && (avgSepNonIdenticalCut < correctedAvgSepNeg))
		{
		  fBkgLamALam->Fill(sepCutIndex, centralityBin+1, pairKstarLam);
		  if(defaultVariableAvgSepCutIndex == sepCutIndex){
		    //This implementation doesn't work properly
		    fKtLamALamBkg->Fill(centralityBin+1,pairKt,pairKstarLam);
		  }
		}
	      }//end loop over sepCutIndex
	    }//end Lam-ALam pair binning
	  }//end mixed event pair histogramming
	}//end past event
      }//end event buffer
    }//end current event
  }//end variable cut loop
}//end DoPairStudies()


