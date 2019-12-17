
//--- Task for investigation of the DoubleHyperHydrogen4 ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+   4LLH decays in the first step to a 4LHe and a pion   +
//+   and in the second step the 4LHe decays to a 3He a    +
//+   proton and a pion.                                   +
//+   We implemented 3 different methods using V0 finder   +
//+   and AliVertexer to reconstruct the 4LLH              +
//+   The first method uses 3 V0s [(4LHe, pi), (3He, pi)   + 
//+   and (p, pi)] and a combination of the tracks.        +
//+   The second method only uses tracks and the           +
//+   AliVertexer Class to build a secondary and tertiary  +
//+   Vertex.                                              +
//+   The third class uses a combination of both: The      +
//+   (4LHe,pi) will be found by the V0 finder and the     +
//+   tertiary vertex will be build by finding the 3       +
//+   tracks and reconstruct the vertex with AliVertexer.  +
//+   Also the BR to Hypertriton(2B), p, pi is checked in  +
//+   the combined Analysis.                               +
//+   m, p, pt and ct of the 4LH are also evaluated.       +
//+                                                        +
//+   Furthermore there is a void for the check of 4Li     +
//+   via He3,p correlation.                               +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom2.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "TPDGCode.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskDoubleHypNucTree.h"
#include "TLorentzVector.h"
#include <TClonesArray.h>
#include "TObject.h"

using namespace std;

ClassImp(AliAnalysisTaskDoubleHypNucTree)

// Default Constructor
AliAnalysisTaskDoubleHypNucTree::AliAnalysisTaskDoubleHypNucTree()
:AliAnalysisTaskSE("AliAnalysisTaskDoubleHypNucTree"),
fPIDCheckOnly(kFALSE),
ftrackAnalysis(kFALSE),
fV0Analysis(kFALSE),
fV0Combination(kFALSE),
fLi4Analysis(kFALSE),
fInputHandler(0),
fPID(0),
fESDevent(0),
fStack(),
fV0(),
fV01(),
fV02(),
fHistdEdx(0),
fHistdEdxV0(0),
fHistNumEvents(0),
fHistTrigger(0),
fHistV0(0),
fTree(0),
fHistogramList(NULL),
fPrimaryVertex(),
fMagneticField(),
fNV0Cand(),
fMCtrue(0),
fEventCuts(),
fPeriod(00),
fTriggerMask(),
fBetheSplines(),
fBetheParamsHe(),
fBetheParamsT(),
fmHHe4Pi(-99), 
fm4LH(-99),
fmLi4(-99),
fmDaughterSum1(-99), 
fmDaughterSum2(-99), 
fVertDiff(-99),
fp4LH(-99),
fpt4LH(-99),
fct4LH(-99),
fpHHe4Pi(-99), 
fctHHe4Pi(-99), 
fptHHe4Pi(-99), 
fpSum1(-99), 
fctSum1(-99), 
fptSum1(-99), 
fpSum2(-99), 
fctSum2(-99), 
fptSum2(-99), 
fpLi4(-99),
fptLi4(-99),
fdcaHe4Pi(-99), 
fdcaHe3Pi(-99), 
fdcaPPi(-99), 
fcosHe4Pi(-99), 
fcosHe3Pi(-99), 
fcosPPi(-99), 
fyHHe4pi(-99), 
fySum1(-99), 
fySum2(-99), 
fhe4DcaSec(-99), 
fhe3DcaTert(-99), 
fpDcaTert(-99), 
fpiDcaSec(-99), 
fpi1DcaTert(-99), 
fpi2DcaTert(-99),  
fpiDca(-99), 
fpi1Dca(-99), 
fpi2Dca(-99), 
fhe4Ncls(-99), 
fpiNcls(-99), 
fhe3Ncls(-99), 
fpi1Ncls(-99), 
fpNcls(-99), 
fpi2Ncls(-99),
fhe4NclsITS(-99), 
fpiNclsITS(-99), 
fhe3NclsITS(-99), 
fpi1NclsITS(-99), 
fpNclsITS(-99), 
fpi2NclsITS(-99), 
fhe4Dca(-99), 
fhe3Dca(-99), 
fpDca(-99), 
fhe4DedxSigma(-99), 
fpiDedxSigma(-99), 
fhe3DedxSigma(-99), 
fpi1DedxSigma(-99), 
fpDedxSigma(-99), 
fpi2DedxSigma(-99), 
ftDedxSigma(-99), 
fhe4P(-99), 
fpiP(-99), 
fhe3P(-99), 
fpi1P(-99), 
fpP(-99), 
fpi2P(-99), 
fhe4Dedx(-99), 
fpiDedx(-99), 
fhe3Dedx(-99), 
fpi1Dedx(-99), 
fpDedx(-99), 
fpi2Dedx(-99), 
farmalpha(-99), 
farmpt(-99),
farmalpha1(-99), 
farmpt1(-99),
farmalpha2(-99), 
farmpt2(-99), 
ftrig(-99), 
fz(-99), 
fmc(-99), 
fthetaP(-99), 
fthetaN(-99),
fDcaHe3P(-99),          //< DCA between He3 Track from V0 and P Track from another V0
fDcaHe3Pi2(-99),       //< DCA between He3 Track from V0 and Pi Track from another V0
fDcaPPi1(-99),          //< DCA between P Track from V0 and Pi Track from another V0
fPAHe3He4(-99),          //< PA between V0(He3 + Pi) and V0(He4 + pi)
fPAPHe4(-99),            //< PA between V0(P + Pi) and V0(He4 + pi)
fPA4LHe(-99),
fPA4LHe1(-99),
fCharge(-99),           //< anti or particle
fParticleSpecies(-99),  //< particle species
fonTheFly(-99),
fonTheFly1(-99),
fonTheFly2(-99),
fVertexPosition(),
fNumberV0s(-99),      //< number of v0s in event
fCentrality(-99),     //< centrality of event
frunnumber(-99),      //< number of run
fTrigger(),        //< array of Triggers
fTriggerClasses(), //< fired trigger classes
fEtaHe4(-99),
fEtaHe3(-99),
fEtaP(-99),
fEtaPi1(-99),
fEtaPi2(-99),
fEtaPi(-99),
fPhiHe4(-99),
fPhiHe3(-99),
fPhiP(-99),
fPhiPi1(-99),
fPhiPi2(-99),
fPhiPi(-99),
fGeoLengthHe4(-99),
fGeoLengthHe3(-99),
fGeoLengthP(-99),
fGeoLengthPi(-99),
fGeoLengthPi1(-99),
fGeoLengthPi2(-99),
fTOFSignalHe4(-99),
fTOFSignalHe3(-99),
fTOFSignalP(-99),
fTOFSignalPi(-99),
fTOFSignalPi1(-99),
fTOFSignalPi2(-99),
fBR(-99)
{

}

// Constructor
AliAnalysisTaskDoubleHypNucTree::AliAnalysisTaskDoubleHypNucTree(const char *name)
:AliAnalysisTaskSE(name),
fPIDCheckOnly(kFALSE),
ftrackAnalysis(kFALSE),
fV0Analysis(kFALSE),
fV0Combination(kFALSE),
fLi4Analysis(kFALSE),
fInputHandler(0),
fPID(0),
fESDevent(0),
fStack(),
fV0(),
fV01(),
fV02(),
fHistdEdx(0),
fHistdEdxV0(0),
fHistNumEvents(0),
fHistTrigger(0),
fHistV0(0),
fTree(0),
fHistogramList(NULL),
fPrimaryVertex(),
fMagneticField(),
fNV0Cand(),
fMCtrue(0),
fEventCuts(),
fPeriod(00),
fTriggerMask(),
fBetheSplines(),
fBetheParamsHe(),
fBetheParamsT(),
fmHHe4Pi(-99), 
fm4LH(-99),
fmLi4(-99),
fmDaughterSum1(-99), 
fmDaughterSum2(-99), 
fVertDiff(-99),
fp4LH(-99),
fpt4LH(-99),
fct4LH(-99),
fpHHe4Pi(-99), 
fctHHe4Pi(-99), 
fptHHe4Pi(-99), 
fpSum1(-99), 
fctSum1(-99), 
fptSum1(-99), 
fpSum2(-99), 
fctSum2(-99), 
fptSum2(-99), 
fpLi4(-99),
fptLi4(-99),
fdcaHe4Pi(-99), 
fdcaHe3Pi(-99), 
fdcaPPi(-99), 
fcosHe4Pi(-99), 
fcosHe3Pi(-99), 
fcosPPi(-99), 
fyHHe4pi(-99), 
fySum1(-99), 
fySum2(-99), 
fhe4DcaSec(-99), 
fhe3DcaTert(-99), 
fpDcaTert(-99), 
fpiDcaSec(-99), 
fpi1DcaTert(-99), 
fpi2DcaTert(-99),  
fpiDca(-99), 
fpi1Dca(-99), 
fpi2Dca(-99), 
fhe4Ncls(-99), 
fpiNcls(-99), 
fhe3Ncls(-99), 
fpi1Ncls(-99), 
fpNcls(-99), 
fpi2Ncls(-99),
fhe4NclsITS(-99), 
fpiNclsITS(-99), 
fhe3NclsITS(-99), 
fpi1NclsITS(-99), 
fpNclsITS(-99), 
fpi2NclsITS(-99), 
fhe4Dca(-99), 
fhe3Dca(-99), 
fpDca(-99), 
fhe4DedxSigma(-99), 
fpiDedxSigma(-99), 
fhe3DedxSigma(-99), 
fpi1DedxSigma(-99), 
fpDedxSigma(-99), 
fpi2DedxSigma(-99), 
ftDedxSigma(-99), 
fhe4P(-99), 
fpiP(-99), 
fhe3P(-99), 
fpi1P(-99), 
fpP(-99), 
fpi2P(-99), 
fhe4Dedx(-99), 
fpiDedx(-99), 
fhe3Dedx(-99), 
fpi1Dedx(-99), 
fpDedx(-99), 
fpi2Dedx(-99), 
farmalpha(-99), 
farmpt(-99),
farmalpha1(-99), 
farmpt1(-99),
farmalpha2(-99), 
farmpt2(-99), 
ftrig(-99), 
fz(-99), 
fmc(-99), 
fthetaP(-99), 
fthetaN(-99),
fDcaHe3P(-99),          //< DCA between He3 Track from V0 and P Track from another V0
fDcaHe3Pi2(-99),       //< DCA between He3 Track from V0 and Pi Track from another V0
fDcaPPi1(-99),          //< DCA between P Track from V0 and Pi Track from another V0
fPAHe3He4(-99),          //< PA between V0(He3 + Pi) and V0(He4 + pi)
fPAPHe4(-99),            //< PA between V0(P + Pi) and V0(He4 + pi)
fPA4LHe(-99),
fPA4LHe1(-99),
fCharge(-99),           //< anti or particle
fParticleSpecies(-99),  //< particle species
fonTheFly(-99),
fonTheFly1(-99),
fonTheFly2(-99),
fVertexPosition(),
fNumberV0s(-99),      //< number of v0s in event
fCentrality(-99),     //< centrality of event
frunnumber(-99),      //< number of run
fTrigger(),        //< array of Triggers
fTriggerClasses(), //< fired trigger classes
fEtaHe4(-99),
fEtaHe3(-99),
fEtaP(-99),
fEtaPi1(-99),
fEtaPi2(-99),
fEtaPi(-99),
fPhiHe4(-99),
fPhiHe3(-99),
fPhiP(-99),
fPhiPi1(-99),
fPhiPi2(-99),
fPhiPi(-99),
fGeoLengthHe4(-99),
fGeoLengthHe3(-99),
fGeoLengthP(-99),
fGeoLengthPi(-99),
fGeoLengthPi1(-99),
fGeoLengthPi2(-99),
fTOFSignalHe4(-99),
fTOFSignalHe3(-99),
fTOFSignalP(-99),
fTOFSignalPi(-99),
fTOFSignalPi1(-99),
fTOFSignalPi2(-99),
fBR(-99)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

// Destructor
AliAnalysisTaskDoubleHypNucTree::~AliAnalysisTaskDoubleHypNucTree() {

}

void AliAnalysisTaskDoubleHypNucTree::UserCreateOutputObjects() {
  fInputHandler = dynamic_cast<AliESDInputHandler*>
  (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!fInputHandler) {
    AliError("Could not get ESD InputHandler.\n");
    return;
  }
  fPID = fInputHandler->GetESDpid();
  if (!fPID) {
    AliError("Could not get PID response.\n");
    return;
  }
  fHistdEdx = new TH2F("fHistdEdX","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);
  fHistdEdxV0 = new TH2F("fHistdEdXV0","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);

  fHistNumEvents = new TH1F("fHistNumEvents","Number of Events",2,0,2);
  fHistNumEvents->GetXaxis()->SetBinLabel(1,"before PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(2,"after PhysSel");

  fHistTrigger = new TH1F("fHistTrigger","Trigger",7,0,7);
  fHistTrigger->GetXaxis()->SetBinLabel(1,"other");
  fHistTrigger->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(5,"HNU");
  fHistTrigger->GetXaxis()->SetBinLabel(6,"HQU");
  fHistTrigger->GetXaxis()->SetBinLabel(7,"HJT");
  fHistV0 = new TH1F("fHistV0","Trigger V0s",7,0,7);
  fHistV0->GetXaxis()->SetBinLabel(1,"other");
  fHistV0->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistV0->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistV0->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistV0->GetXaxis()->SetBinLabel(5,"HNU");
  fHistV0->GetXaxis()->SetBinLabel(6,"HQU");
  fHistV0->GetXaxis()->SetBinLabel(7,"HJT");

  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  fHistogramList->SetName(GetName());
  fHistogramList->Add(fHistdEdx);
  fHistogramList->Add(fHistdEdxV0);
  fHistogramList->Add(fHistNumEvents);
  fHistogramList->Add(fHistTrigger);
  fHistogramList->Add(fHistV0);

  fEventCuts.AddQAplotsToList(fHistogramList);
  //TREE for only V0 and combined track V0 analysis
  fTree = new TTree("tree","fTree");
  //Masses
  fTree->Branch("fmHHe4Pi", &fmHHe4Pi, "fmHHe4Pi/F");
  fTree->Branch("fmDaughterSum1", &fmDaughterSum1, "fmDaughterSum1/F");
  fTree->Branch("fmDaughterSum2", &fmDaughterSum2, "fmDaughterSum2/F");
  fTree->Branch("fmLi4", &fmLi4, "fmLi4/F");
  fTree->Branch("fm4LH", &fm4LH, "fm4LH/F");
  //Corr
  fTree->Branch("fDcaHe3P", &fDcaHe3P, "fDcaHe3P/F");
  fTree->Branch("fDcaHe3Pi2", &fDcaHe3Pi2, "fDcaHe3Pi2/F");
  fTree->Branch("fDcaPPi1", &fDcaPPi1, "fDcaPPi1/F");
  fTree->Branch("fPAHe3He4", &fPAHe3He4, "fPAHe3He4/F");
  fTree->Branch("fPAPHe4", &fPAPHe4, "fPAPHe4/F");
  fTree->Branch("fPA4LHe", &fPA4LHe, "fPA4LHe/F");
  fTree->Branch("fPA4LHe1", &fPA4LHe1, "fPA4LHe1/F");
  fTree->Branch("fVertDiff", &fVertDiff, "fVertDiff/F");
  //P and Pt
  fTree->Branch("fpHHe4Pi", &fpHHe4Pi, "fpHHe4Pi/F");
  fTree->Branch("fptHHe4Pi", &fptHHe4Pi, "fptHHe4Pi/F");  
  fTree->Branch("fpSum1", &fpSum1, "fpSum1/F");
  fTree->Branch("fptSum1", &fptSum1, "fptSum1/F");
  fTree->Branch("fpSum2", &fpSum2, "fpSum2/F");
  fTree->Branch("fptSum2", &fptSum2, "fptSum2/F");
  fTree->Branch("fpLi4", &fpLi4, "fpLi4/F");
  fTree->Branch("fptLi4", &fptLi4, "fptLi4/F");
  fTree->Branch("fp4LH", &fp4LH, "fp4LH/F");
  fTree->Branch("fpt4LH", &fpt4LH, "fpt4LH/F");
  //Ct
  fTree->Branch("fctHHe4Pi", &fctHHe4Pi, "fctHHe4Pi/F");
  fTree->Branch("fctSum1", &fctSum1, "fctSum1/F");
  fTree->Branch("fctSum2", &fctSum2, "fctSum2/F");
  fTree->Branch("fct4LH", &fct4LH, "fct4LH/F");
  //Particle P
  fTree->Branch("fhe4P", &fhe4P, "fhe4P/F");
  fTree->Branch("fpiP", &fpiP, "fpiP/F");
  fTree->Branch("fhe3P", &fhe3P, "fhe3P/F");
  fTree->Branch("fpi1P", &fpi1P, "fpi1P/F");
  fTree->Branch("fpP", &fpP, "fpP/F");
  fTree->Branch("fpi2P", &fpi2P, "fpi2P/F");
  //DCA's
  fTree->Branch("fdcaHe4Pi", &fdcaHe4Pi, "fdcaHe4Pi/F");
  fTree->Branch("fdcaHe3Pi", &fdcaHe3Pi, "fdcaHe3Pi/F");
  fTree->Branch("fdcaPPi", &fdcaPPi, "fdcaPPi/F");
  //PA's
  fTree->Branch("fcosHe4Pi", &fcosHe4Pi, "fcosHe4Pi/F");
  fTree->Branch("fcosHe3Pi", &fcosHe3Pi, "fcosHe3Pi/F");
  fTree->Branch("fcosPPi", &fcosPPi, "fcosPPi/F");
  //Rapidity
  fTree->Branch("fyHHe4pi", &fyHHe4pi, "fyHHe4pi/F");
  fTree->Branch("fySum1", &fySum1, "fySum1/F");
  fTree->Branch("fySum2", &fySum2, "fySum2/F");
  //DCA From Primary Vertex
  fTree->Branch("fpiDca", &fpiDca, "fpiDca/F");
  fTree->Branch("fpi1Dca", &fpi1Dca, "fpi1Dca/F");
  fTree->Branch("fpi2Dca", &fpi2Dca, "fpi2Dca/F");
  fTree->Branch("fhe4Dca", &fhe4Dca ,"fhe4Dca/F");
  fTree->Branch("fhe3Dca", &fhe3Dca ,"fhe3Dca/F");
  fTree->Branch("fpDca", &fpDca ,"fpDca/F");
  //DCA From secondary/tertiary Vertex
  fTree->Branch("fpiDcaSec", &fpiDcaSec, "fpiDcaSec/F");
  fTree->Branch("fpi1DcaTert", &fpi1DcaTert, "fpi1DcaTert/F");
  fTree->Branch("fpi2DcaTert", &fpi2DcaTert, "fpi2DcaTert/F");
  fTree->Branch("fhe4DcaSec", &fhe4DcaSec ,"fhe4DcaSec/F");
  fTree->Branch("fhe3DcaTert", &fhe3DcaTert ,"fhe3DcaTert/F");
  fTree->Branch("fpDcaTert", &fpDcaTert ,"fpDcaTert/F");
  //Number of Clusters
  fTree->Branch("fhe4Ncls", &fhe4Ncls, "fhe4Ncls/F");
  fTree->Branch("fpiNcls", &fpiNcls, "fpiNcls/F");
  fTree->Branch("fhe3Ncls", &fhe3Ncls, "fhe3Ncls/F");
  fTree->Branch("fpi1Ncls", &fpi1Ncls, "fpi1Ncls/F");
  fTree->Branch("fpNcls", &fpNcls, "fpNcls/F");
  fTree->Branch("fpi2Ncls", &fpi2Ncls, "fpi2Ncls/F");
  //Number of Clusters ITS
  fTree->Branch("fhe4NclsITS", &fhe4NclsITS, "fhe4NclsITS/F");
  fTree->Branch("fpiNclsITS", &fpiNclsITS, "fpiNclsITS/F");
  fTree->Branch("fhe3NclsITS", &fhe3NclsITS, "fhe3NclsITS/F");
  fTree->Branch("fpi1NclsITS", &fpi1NclsITS, "fpi1NclsITS/F");
  fTree->Branch("fpNclsITS", &fpNclsITS, "fpNclsITS/F");
  fTree->Branch("fpi2NclsITS", &fpi2NclsITS, "fpi2NclsITS/F");
  //Number of Sigmas PID
  fTree->Branch("fhe4DedxSigma", &fhe4DedxSigma, "fhe4DedxSigma/F");
  fTree->Branch("fhe3DedxSigma", &fhe3DedxSigma, "fhe3DedxSigma/F");
  fTree->Branch("fpDedxSigma", &fpDedxSigma, "fpDedxSigma/F");
  fTree->Branch("fpiDedxSigma", &fpiDedxSigma, "fpiDedxSigma/F");
  fTree->Branch("fpi1DedxSigma", &fpi1DedxSigma, "fpi1DedxSigma/F");
  fTree->Branch("fpi2DedxSigma", &fpi2DedxSigma, "fpi2DedxSigma/F");
  //TPC Signal
  fTree->Branch("fhe4Dedx", &fhe4Dedx, "fhe4Dedx/F");
  fTree->Branch("fpiDedx", &fpiDedx, "fpiDedx/F");
  fTree->Branch("fhe3Dedx", &fhe3Dedx, "fhe3Dedx/F");
  fTree->Branch("fpi1Dedx", &fpi1Dedx, "fpi1Dedx/F");
  fTree->Branch("fpDedx", &fpDedx, "fpDedx/F");
  fTree->Branch("fpi2Dedx", &fpi2Dedx, "fpi2Dedx/F");
  //Armenteros
  fTree->Branch("farmalpha", &farmalpha, "farmalpha/F");
  fTree->Branch("farmpt", &farmpt, "farmpt/F");
  fTree->Branch("farmalpha1", &farmalpha1, "farmalpha1/F");
  fTree->Branch("farmpt1", &farmpt1, "farmpt1/F");
  //Eta
  fTree->Branch("fEtaHe4", &fEtaHe4, "fEtaHe4/F");
  fTree->Branch("fEtaHe3", &fEtaHe3, "fEtaHe3/F");
  fTree->Branch("fEtaP", &fEtaP, "fEtaP/F");
  fTree->Branch("fEtaPi", &fEtaPi, "fEtaPi/F");
  fTree->Branch("fEtaPi1", &fEtaPi1, "fEtaPi1/F");
  fTree->Branch("fEtaPi2", &fEtaPi2, "fEtaPi2/F");
  //Phi
  fTree->Branch("fPhiHe4", &fPhiHe4, "fPhiHe4/F");
  fTree->Branch("fPhiHe3", &fPhiHe3, "fPhiHe3/F");
  fTree->Branch("fPhiP", &fPhiP, "fPhiP/F");
  fTree->Branch("fPhiPi", &fPhiPi, "fPhiPi/F");
  fTree->Branch("fPhiPi1", &fPhiPi1, "fPhiPi1/F");
  fTree->Branch("fPhiPi2", &fPhiPi2, "fPhiPi2/F");
  //GeoLength
  fTree->Branch("fGeoLengthHe4", &fGeoLengthHe4, "fGeoLengthHe4/F");
  fTree->Branch("fGeoLengthHe3", &fGeoLengthHe3, "fGeoLengthHe3/F");
  fTree->Branch("fGeoLengthP", &fGeoLengthP, "fGeoLengthP/F");
  fTree->Branch("fGeoLengthPi", &fGeoLengthPi, "fGeoLengthPi/F");
  fTree->Branch("fGeoLengthPi1", &fGeoLengthPi1, "fGeoLengthPi1/F");
  fTree->Branch("fGeoLengthPi2", &fGeoLengthPi2, "fGeoLengthPi2/F");
  //TOF Signal
  fTree->Branch("fTOFSignalHe4", &fTOFSignalHe4, "fTOFSignalHe4/F");
  fTree->Branch("fTOFSignalHe3", &fTOFSignalHe3, "fTOFSignalHe3/F");
  fTree->Branch("fTOFSignalP", &fTOFSignalP, "fTOFSignalP/F");
  fTree->Branch("fTOFSignalPi", &fTOFSignalPi, "fTOFSignalPi/F");
  fTree->Branch("fTOFSignalPi1", &fTOFSignalPi1, "fTOFSignalPi1/F");
  fTree->Branch("fTOFSignalPi2", &fTOFSignalPi2, "fTOFSignalPi2/F");
  //fTree->Branch("trig", &trig, "trig/F");
  fTree->Branch("fParticleSpecies", &fParticleSpecies, "fParticleSpecies/I");
  fTree->Branch("fonTheFly", &fonTheFly, "fonTheFly/I");
  fTree->Branch("fonTheFly1", &fonTheFly1, "fonTheFly1/I");
  fTree->Branch("fonTheFly2", &fonTheFly2, "fonTheFly2/I");
  fTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
  fTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
  fTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
  fTree->Branch("fz", &fz, "fz/I");
  fTree->Branch("fmc", &fmc, "fmc/I");
  fTree->Branch("fBR", &fBR, "fBR/I");
  
  PostData(1, fHistogramList);
  PostData(2, fTree);

  //********--------******** Info ********--------********
  //  EParticleType :                                   //
  //  kElectron = 0, kMuon = 1, kPion = 2, kKaon = 3,   //  
  //  kProton = 4, kDeuteron = 5, kTriton = 6, kHe3 = 7 //
  //  kAlpha = 8, kPhoton = 9, kPi0 = 10, kNeutron = 11 //
  //  kKaon0 = 12, kEleCon = 13, kUnknown = 14          //
  //                                                    //
  //  Masses: 4Li: 3.74958 GeV/c^2                      //
  //          4He: 3.727379 GeV/c^2                     //
  //          3He: 2.80923 GeV/c^2                      //
  //            t: 2.80925 GeV/c^2                      //
  //            d: 1.875613 GeV/c^2                     //
  //            p: 0.93827 GeV/c^2                      //
  //           pi: 0.13957 GeV/c^2                      //
  //          3LH: 2.99131 GeV/c^2                      //
  //          4LH: 3.931   GeV/c^2                      //
  //         4LHe: 3.929   GeV/c^2                      //
  //         5LHe: 4.841   GeV/c^2                      //
  //         4LLH: 4.106   GeV/c^2                      //
  //*******-------********--------********--------********

}

void AliAnalysisTaskDoubleHypNucTree::UserExec(Option_t *) {
  // MC
  fMCtrue = kTRUE;
  AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>
  (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcEventHandler) {
    fMCtrue = kFALSE;
  }
  AliMCEvent* mcEvent = 0x0;
  if (mcEventHandler) mcEvent = mcEventHandler->MCEvent();
  if (!mcEvent) {
    if (fMCtrue) return;
  }
  if (fMCtrue) {
    fStack = mcEvent->Stack();
    if (!fStack) return;
  }
  // Data
  fESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESDevent) {
    AliError("Could not get ESD Event.\n");
    return;
  }
  if (!fPID) {
    AliError("Could not get PID response.\n");
    return;
  }
  
  fHistNumEvents->Fill(0);
  Float_t centrality = -1;
  const AliESDVertex *vertex = fESDevent->GetPrimaryVertexSPD();
  fEventCuts.OverrideAutomaticTriggerSelection(fTriggerMask);
  if (fPeriod == 2010 || fPeriod == 2011) {
    if (vertex->GetNContributors() < 1) {
      vertex = fESDevent->GetPrimaryVertexSPD();
      if (vertex->GetNContributors() < 1) {
       PostData(1,fHistogramList);
       return;
     }
   }
   if (TMath::Abs(vertex->GetZ()) > 10) {
    PostData(1, fHistogramList);
    return;
  }
  centrality = fESDevent->GetCentrality()->GetCentralityPercentile("V0M");
  if(!fMCtrue){
    if (centrality < 0.0 || centrality > 100.0 ) {
      return;
    }
  }
}
if (fPeriod == 2015) {
  if(!fEventCuts.AcceptEvent(fESDevent)) {
    PostData(1,fHistogramList);
    return;
  }
    // 0 = V0M
  centrality = fEventCuts.GetCentrality(0);
  if(!fMCtrue){
    if (centrality < 0.0 || centrality > 100.0 ) {
      return;
    }
  }
}
if (fPeriod == 2016 || fPeriod == 2017 || fPeriod == 2018) {
  Int_t r = fESDevent->GetRunNumber();
  if(r == 297219 || r == 297194 || r == 297029 || r == 296890 || r == 296849 || r == 296750 || r == 296749) fEventCuts.UseTimeRangeCut(); 
  if(!fEventCuts.AcceptEvent(fESDevent)) {
    PostData(1,fHistogramList);
    return;
  }
  // 0 = V0M
  centrality = fEventCuts.GetCentrality(0);
}
AliESDVertex *esdVer1 = new AliESDVertex(*vertex);
AliVertexerTracks *vertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
AliVertexerTracks *vertexer1 = new AliVertexerTracks(fESDevent->GetMagneticField());
vertexer->SetVtxStart(esdVer1);

Double_t *dn = new Double_t[3];
Double_t *dd = new Double_t[3];
dn[0] = esdVer1->GetX();
dn[1] = esdVer1->GetY();
dn[2] = esdVer1->GetZ();
if(esdVer1) delete esdVer1;

Int_t runNumber = fESDevent->GetRunNumber();
frunnumber = runNumber;
TriggerSelection();
//Number of Events
fHistNumEvents->Fill(1);
//Centrality
fCentrality = centrality;
//MagneticField
fMagneticField  = fESDevent->GetMagneticField();
//Primary Vertex Position
fPrimaryVertex.SetXYZ(vertex->GetX(),vertex->GetY(),vertex->GetZ());
fVertexPosition = fPrimaryVertex; 
//V0
fNV0Cand = 0;
//For DCA
Double_t xthiss(0.0);
Double_t xpp(0.0);

AliESDtrackCuts trackCutsV0("AlitrackCutsV0", "AlitrackCutsV0");
AliESDtrackCuts trackCutsSec("AlitrackCutsSec", "AlitrackCutsSec");
AliESDtrackCuts trackCutsTert("AlitrackCutsTert", "AlitrackCutsTert");
AliESDtrackCuts trackCutsPi("AlitrackCutsPi", "AlitrackCutsPi");

if(fV0Analysis){
  trackCutsV0.SetEtaRange(-0.9,0.9);
  trackCutsV0.SetAcceptKinkDaughters(kFALSE);
  trackCutsV0.SetRequireTPCRefit(kTRUE);
  trackCutsV0.SetMaxChi2PerClusterTPC(5);
  trackCutsV0.SetMinNClustersTPC(60);
}
if(ftrackAnalysis){
  trackCutsSec.SetEtaRange(-0.9,0.9);
  trackCutsSec.SetAcceptKinkDaughters(kFALSE);
  trackCutsSec.SetRequireTPCRefit(kTRUE);
  trackCutsSec.SetMaxChi2PerClusterTPC(5);
  trackCutsSec.SetMinNClustersTPC(60);
  trackCutsSec.SetMaxRel1PtUncertainty(0.2);
  trackCutsSec.SetPtRange(0.0, 10.0);

  trackCutsTert.SetEtaRange(-0.9,0.9);
  trackCutsTert.SetAcceptKinkDaughters(kFALSE);
  trackCutsTert.SetRequireTPCRefit(kTRUE);
  trackCutsTert.SetMaxChi2PerClusterTPC(5);
  trackCutsTert.SetMinNClustersTPC(60);
  trackCutsTert.SetMinDCAToVertexXY(0.3);
  trackCutsTert.SetMinDCAToVertexZ(0.3);
  trackCutsTert.SetMaxRel1PtUncertainty(0.2);
  trackCutsTert.SetPtRange(0.0, 10.0);

  trackCutsPi.SetEtaRange(-0.9,0.9);
  trackCutsPi.SetAcceptKinkDaughters(kFALSE);
  trackCutsPi.SetRequireTPCRefit(kTRUE);
  trackCutsPi.SetMaxChi2PerClusterTPC(5);
  trackCutsPi.SetMinNClustersTPC(60);
  trackCutsPi.SetMaxRel1PtUncertainty(0.2);
  trackCutsPi.SetPtRange(0.0, 10.0);
}
if(fV0Combination){

  trackCutsV0.SetEtaRange(-0.9,0.9);
  trackCutsV0.SetAcceptKinkDaughters(kFALSE);
  trackCutsV0.SetRequireTPCRefit(kTRUE);
  trackCutsV0.SetMaxChi2PerClusterTPC(5);
  trackCutsV0.SetMinNClustersTPC(60);
  trackCutsV0.SetMaxRel1PtUncertainty(0.1);
  trackCutsV0.SetPtRange(0.0, 10.0);

  trackCutsTert.SetEtaRange(-0.9,0.9);
  trackCutsTert.SetAcceptKinkDaughters(kFALSE);
  trackCutsTert.SetRequireTPCRefit(kTRUE);
  trackCutsTert.SetMaxChi2PerClusterTPC(5);
  trackCutsTert.SetMinNClustersTPC(60);
  trackCutsTert.SetMaxRel1PtUncertainty(0.1);
  trackCutsTert.SetPtRange(0.0, 10.0);

  trackCutsPi.SetEtaRange(-0.9,0.9);
  trackCutsPi.SetAcceptKinkDaughters(kFALSE);
  trackCutsPi.SetRequireTPCRefit(kTRUE);
  trackCutsPi.SetMaxChi2PerClusterTPC(5);
  trackCutsPi.SetMinNClustersTPC(60);
  trackCutsPi.SetMaxRel1PtUncertainty(0.1);
  trackCutsPi.SetPtRange(0.0, 10.0);
}

Bool_t ITSClusterCut = kFALSE;
Int_t SecClusters = 0;
Int_t TertClusters = 0;

if(ITSClusterCut){
    SecClusters = 5; //has to be >=
    TertClusters = 2; // has to be <=
  }
  // Pidqa loop
  if(fPIDCheckOnly){
    dEdxCheck();
  }
  //Only use V0s to reconstruct 4LLH
  if(fV0Analysis){
    V0Analysis(trackCutsV0, xthiss, xpp, ITSClusterCut, SecClusters, TertClusters);
  }
  //only use Tracks to reconstruct 4LLH
  if(ftrackAnalysis){
    TrackAnalysis(trackCutsSec, trackCutsTert, trackCutsPi, vertexer, vertexer1, dn, dd, xthiss, xpp, ITSClusterCut, SecClusters, TertClusters);
  }
  //Use a Combination of V0 and Tracks to reconstruct 4LLH
  if(fV0Combination){
    CombinedAnalysis(trackCutsV0, trackCutsTert, trackCutsPi, vertexer, vertexer1, dn, dd, xthiss, xpp, ITSClusterCut, SecClusters, TertClusters);
  }
  //Check for Li4 in He3,p Correlations
  if(fLi4Analysis){
    Li4Analysis(trackCutsSec);
  }

  PostData(1, fHistogramList);
  PostData(2, fTree);
}
void AliAnalysisTaskDoubleHypNucTree::dEdxCheck(){
  AliESDtrackCuts* trackCutsPid = new AliESDtrackCuts("trackCutsPid", "trackCutsPid");
  trackCutsPid = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  trackCutsPid->SetEtaRange(-0.9,0.9);
  for (Int_t itrack = 0; itrack < fESDevent->GetNumberOfTracks(); itrack++) {
    AliESDtrack* track = fESDevent->GetTrack(itrack);
    if (!trackCutsPid->AcceptTrack(track)) continue;
    Double_t momentum = track->GetInnerParam()->GetP();
    fHistdEdx->Fill(momentum * track->GetSign(), track->GetTPCsignal());
  }
  delete trackCutsPid;
}
void AliAnalysisTaskDoubleHypNucTree::Li4Analysis(AliESDtrackCuts trackCutsSec){
  Int_t count = 0;
  //++++++++++++++++++++++++++++ He3 Track +++++++++++++++++++++++++++++++++++++//
  for (Int_t ATracks = 0; ATracks < fESDevent->GetNumberOfTracks(); ATracks++) {
    
    AliESDtrack* trackA = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(ATracks));

    if (!trackA->GetInnerParam()) continue;

    if (!trackCutsSec.AcceptTrack(trackA)) continue;

    Double_t ptotA = trackA->GetInnerParam()->GetP();
    Double_t signA = trackA->GetSign();    
    
    fHistdEdx->Fill(ptotA*signA, trackA->GetTPCsignal());
    
    if(fBetheSplines){
      if(TMath::Abs(fPID->NumberOfSigmasTPC(trackA, AliPID::kHe3)) > 3) continue;
    }
    else {
      if(TMath::Abs(Bethe(*trackA, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) > 3) continue;
    }       
    //+++++++++++++++++++++++++++++++ p Track +++++++++++++++++++++++++++++++++++++++++++//
    for (Int_t BTracks = ATracks + 1; BTracks < fESDevent->GetNumberOfTracks(); BTracks++) {

      AliESDtrack* trackB = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(BTracks));

      if (!trackB->GetInnerParam()) continue;

      if (!trackCutsSec.AcceptTrack(trackB)) continue;

      if(ATracks == BTracks) continue;     

      Double_t ptotB = trackB->GetInnerParam()->GetP();
      Double_t signB = trackB->GetSign();

      if (signB != signA) continue;

      fHistdEdx->Fill(ptotB*signB, trackB->GetTPCsignal());

      if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(trackB, AliPID::kProton)) > 3) continue;
      }
      else {
        if(TMath::Abs(Bethe(*trackB, AliPID::ParticleMass(AliPID::kProton), 2, fBetheParamsT)) > 3) continue;
      }
      //------------------------------------------ --------------------------------------------------------------//
      if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(trackA, AliPID::kHe3)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*trackA, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 3)){
        if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(trackB, AliPID::kProton)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*trackB, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3)) {
          //He3
          TLorentzVector fd(0.,0.,0.,0.);
          fd.SetXYZM(2*trackA->Px(), 2*trackA->Py(), 2*trackA->Pz(), AliPID::ParticleMass(AliPID::kHe3));
          //p
          TLorentzVector sd(0.,0.,0.,0.);
          sd.SetXYZM(trackB->Px(), trackB->Py(), trackB->Pz(), AliPID::ParticleMass(AliPID::kProton));
          //Mother
          TLorentzVector mother = fd + sd;
          //Daughter P
          fhe3P = fd.P();
          fpP = sd.P();
          //Mother m, p, pt
          fmLi4 = mother.M();
          fpLi4 = mother.P();
          fptLi4 = mother.Pt();
          //track TOF Signal
          fTOFSignalHe3 = trackA->GetTOFsignal();
          fTOFSignalP = trackB->GetTOFsignal();
          //dEdx Sigma
          if (fBetheSplines) {
            fhe3DedxSigma = fPID->NumberOfSigmasTPC(trackA, AliPID::kHe3);
            fpDedxSigma = fPID->NumberOfSigmasTPC(trackB, AliPID::kProton);
          } else {
            fhe3DedxSigma = Bethe(*trackA, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
            fpDedxSigma = Bethe(*trackB, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
          }
          //N Clusters TPC
          fhe3Ncls = trackA->GetTPCNcls();
          fpNcls = trackB->GetTPCNcls();
          //NClusters ITS
          fhe3NclsITS = trackA->GetNumberOfITSClusters();
          fpNclsITS = trackB->GetNumberOfITSClusters();
          //BR means kind of a bool for a later analysis of the reduced tree
          fBR = 3;
          //charge
          if(signA >0 && signB >0) fz = 3;
          if(signA<0 && signB<0) fz = -3;
          //Counter
          count++;
        }
      }
    }
    if (count>0) fTree->Fill();
  } 
}
void AliAnalysisTaskDoubleHypNucTree::V0Analysis(AliESDtrackCuts trackCutsV0, Double_t xthiss, Double_t xpp, Bool_t ITSClusterCut, Int_t SecClusters, Int_t TertClusters){
  //+++++++++++++++++++++++++++V0(4LHe, pi)++++++++++++++++++++++++++++++++++
  for (Int_t ivertex = 0; ivertex < fESDevent->GetNumberOfV0s(); ivertex++) {

    fHistV0->Fill(fTrigger);
    fV0 = fESDevent->GetV0(ivertex);
    
    Bool_t v0ChargeCorrect = kTRUE; 

    AliESDtrack* trackN = fESDevent->GetTrack(fV0->GetIndex(0));
    AliESDtrack* trackP = fESDevent->GetTrack(fV0->GetIndex(1));

    if (trackN->GetSign() > 0 ) {
      trackN = fESDevent->GetTrack(fV0->GetIndex(1));
      trackP = fESDevent->GetTrack(fV0->GetIndex(0));
      v0ChargeCorrect = kFALSE;
    }
    
    Double_t sign5 = trackP->GetSign();
    Double_t sign4 = trackN->GetSign();
    
    if (!trackCutsV0.AcceptTrack(trackN)) continue;
    if (!trackCutsV0.AcceptTrack(trackP)) continue;

    fHistdEdxV0->Fill(trackP->GetInnerParam()->GetP() * trackP->GetSign(), trackP->GetTPCsignal());
    fHistdEdxV0->Fill(trackN->GetInnerParam()->GetP() * trackN->GetSign(), trackN->GetTPCsignal());

    if(fPIDCheckOnly) continue;
    
    //special ITS Cluster cut, defined in UserExec
    if(ITSClusterCut){
      if(trackP->GetNumberOfITSClusters() <= SecClusters || trackN->GetNumberOfITSClusters() <= SecClusters) continue;  
    }

    Bool_t pionPositive     = kFALSE;
    Bool_t pionNegative     = kFALSE;
    Bool_t helium4Positive  = kFALSE;
    Bool_t helium4Negative  = kFALSE;
    
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kPion)) < 3) {
      pionPositive = kTRUE;
    }
    else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kPion)) < 3) {
      pionNegative = kTRUE;
    }
    else continue;
    //Use Framework Splines for p, He3, Alpha
    if (fBetheSplines) {
      if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kAlpha)) < 4) {
        helium4Positive = kTRUE;
      } else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kAlpha)) < 4) {
        helium4Negative = kTRUE;
      }
      else continue;
    //Use own Splines for p, He3, Alpha
    } else {

      if (TMath::Abs(Bethe(*trackP, AliPID::ParticleMass(AliPID::kAlpha),  2, fBetheParamsHe)) < 3) {
        helium4Positive = kTRUE;
      } else if (TMath::Abs(Bethe(*trackN, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe)) < 3) {
        helium4Negative = kTRUE;
      }
      else continue;
    }
    //+++++++++++++++++++++++++++V0(3He, pi)+++++++++++++++++++++++++++++++++++
    for (Int_t jvertex = 0; jvertex < fESDevent->GetNumberOfV0s(); jvertex++) {

      fHistV0->Fill(fTrigger);
      fV01 = fESDevent->GetV0(jvertex);

      Bool_t v0ChargeCorrect = kTRUE; 

      AliESDtrack* trackN1 = fESDevent->GetTrack(fV01->GetIndex(0));
      AliESDtrack* trackP1 = fESDevent->GetTrack(fV01->GetIndex(1));

      if (trackN1->GetSign() > 0 ) {
        trackN1 = fESDevent->GetTrack(fV01->GetIndex(1));
        trackP1 = fESDevent->GetTrack(fV01->GetIndex(0));
        v0ChargeCorrect = kFALSE;
      }
      Double_t sign1 = trackP1->GetSign();
      Double_t sign31 = trackN1->GetSign();
      if (!trackCutsV0.AcceptTrack(trackN1)) continue;
      if (!trackCutsV0.AcceptTrack(trackP1)) continue;

      fHistdEdxV0->Fill(trackP1->GetInnerParam()->GetP() * trackP1->GetSign(), trackP1->GetTPCsignal());
      fHistdEdxV0->Fill(trackN1->GetInnerParam()->GetP() * trackN1->GetSign(), trackN1->GetTPCsignal());

      if(fPIDCheckOnly) continue;
      
      //special ITS Cluster cut, defined in UserExec
      if(ITSClusterCut){
        if(trackP1->GetNumberOfITSClusters() >= TertClusters || trackN1->GetNumberOfITSClusters() >= TertClusters) continue;  
      }
      Bool_t pionPositive1     = kFALSE;
      Bool_t pionNegative1     = kFALSE;
      Bool_t helium3Positive  = kFALSE;
      Bool_t helium3Negative  = kFALSE;
      if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP1, AliPID::kPion)) < 3) {
        pionPositive1 = kTRUE;
      }
      else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN1, AliPID::kPion)) < 3) {
        pionNegative1 = kTRUE;
      }
      else continue;
      //Use Framework Splines for p, He3, Alpha
      if (fBetheSplines) {  
        if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP1, AliPID::kHe3)) < 4) {
          helium3Positive = kTRUE;
        } else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN1, AliPID::kHe3)) < 4) {
          helium3Negative = kTRUE;
        }
        else continue;
      //Use own Splines for p, He3, Alpha
      } else {
        if (TMath::Abs(Bethe(*trackP1, AliPID::ParticleMass(AliPID::kHe3),  2, fBetheParamsHe)) < 3) {
          helium3Positive = kTRUE;
        } else if (TMath::Abs(Bethe(*trackN1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 3) {
          helium3Negative = kTRUE;
        }
        else continue;
      }
      //+++++++++++++++++++++++++++V0(p, pi)+++++++++++++++++++++++++++++++++++++
      for (Int_t kvertex = 0; kvertex < fESDevent->GetNumberOfV0s(); kvertex++) {

        fHistV0->Fill(fTrigger);
        fV02 = fESDevent->GetV0(kvertex);
        Bool_t v0ChargeCorrect = kTRUE; 
        AliESDtrack* trackN2 = fESDevent->GetTrack(fV02->GetIndex(0));
        AliESDtrack* trackP2 = fESDevent->GetTrack(fV02->GetIndex(1));
        
        if (trackN2->GetSign() > 0 ) {
          trackN2 = fESDevent->GetTrack(fV02->GetIndex(1));
          trackP2 = fESDevent->GetTrack(fV02->GetIndex(0));
          v0ChargeCorrect = kFALSE;
        }
        Double_t sign2 = trackP2->GetSign();
        Double_t sign32 = trackN2->GetSign();
        if (!trackCutsV0.AcceptTrack(trackN2)) continue;
        if (!trackCutsV0.AcceptTrack(trackP2)) continue;

        fHistdEdxV0->Fill(trackP2->GetInnerParam()->GetP() * trackP2->GetSign(), trackP2->GetTPCsignal());
        fHistdEdxV0->Fill(trackN2->GetInnerParam()->GetP() * trackN2->GetSign(), trackN2->GetTPCsignal());

        if(fPIDCheckOnly) continue;
        
        //special ITS Cluster cut, defined in UserExec
        if(ITSClusterCut){
          if(trackP2->GetNumberOfITSClusters() >= TertClusters || trackN2->GetNumberOfITSClusters() >= TertClusters) continue;  
        }
        
        Bool_t pionPositive2     = kFALSE;
        Bool_t pionNegative2     = kFALSE;
        Bool_t protonPositive   = kFALSE;
        Bool_t protonNegative   = kFALSE;
        if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP2, AliPID::kPion)) < 3) {
          pionPositive2 = kTRUE;
        }
        else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN2, AliPID::kPion)) < 3) {
          pionNegative2 = kTRUE;
        }
        else continue;
        //Use Framework Splines for p, He3, Alpha
        if (fBetheSplines) {
          if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP2, AliPID::kProton)) < 4) {
            protonPositive = kTRUE;
          } else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN2, AliPID::kProton)) < 4) {
            protonNegative = kTRUE;
          }
          else continue;
        //Use own Splines for p, He3, Alpha
        } else {
          if (TMath::Abs(Bethe(*trackP2, AliPID::ParticleMass(AliPID::kProton),  1, fBetheParamsT)) < 3) {
            protonPositive = kTRUE;
          } else if (TMath::Abs(Bethe(*trackN2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3) {
            protonNegative = kTRUE;
          }
          else continue;
        }
        
        if (helium4Positive && pionNegative && helium3Positive && pionNegative1 && protonPositive && pionNegative2) {
          //BR for separation of different analysis: 1 = 4LLH --> (4LHe, pi), 2 = 4LLH -->(3LH, p, pi), 3 = Li4
          fBR = 1;
          fmc = -1; 
          fz = 1;
          fParticleSpecies = 1020010040;
          //V0 DCA
          fdcaHe4Pi = fV0->GetDcaV0Daughters();
          fdcaHe3Pi = fV01->GetDcaV0Daughters();
          fdcaPPi = fV02->GetDcaV0Daughters();
          //V0 Pointing-Angle
          fcosHe4Pi = fV0->GetV0CosineOfPointingAngle();
          fcosHe3Pi = fV01->GetV0CosineOfPointingAngle();
          fcosPPi = fV02->GetV0CosineOfPointingAngle();
          //V0-OnFlyStatus
          if(fV0->GetOnFlyStatus()) fonTheFly = 1;
          else if(!fV0->GetOnFlyStatus()) fonTheFly = 0;
          if(fV01->GetOnFlyStatus()) fonTheFly1 = 1;
          else if(!fV01->GetOnFlyStatus()) fonTheFly1 = 0;
          if(fV02->GetOnFlyStatus()) fonTheFly2 = 1;
          else if(!fV02->GetOnFlyStatus()) fonTheFly2 = 0;
          //Track TPC Signal
          fhe4Dedx = trackP->GetTPCsignal();
          fhe3Dedx = trackP1->GetTPCsignal();
          fpDedx = trackP2->GetTPCsignal();
          fpiDedx = trackN->GetTPCsignal();
          fpi1Dedx = trackN1->GetTPCsignal();
          fpi2Dedx = trackN2->GetTPCsignal();
          //Eta
          fEtaHe4 = trackP->Eta();
          fEtaHe3 = trackP1->Eta();
          fEtaP = trackP2->Eta();
          fEtaPi = trackN->Eta();
          fEtaPi1 = trackN1->Eta();
          fEtaPi2 = trackN2->Eta();
          //Phi
          fPhiHe4 = trackP->Phi();
          fPhiHe3 = trackP1->Phi();
          fPhiP = trackP2->Phi();
          fPhiPi = trackN->Phi();
          fPhiPi1 = trackN1->Phi();
          fPhiPi2 = trackN2->Phi();
          //GeoLength
          fGeoLengthHe4 = GeoLength(*trackP);
          fGeoLengthHe3 = GeoLength(*trackP1);
          fGeoLengthP = GeoLength(*trackP2);
          fGeoLengthPi = GeoLength(*trackN);
          fGeoLengthPi1 = GeoLength(*trackN1);
          fGeoLengthPi2 = GeoLength(*trackN2);
          //TOF Signal
          fTOFSignalHe4 = trackP->GetTOFsignal();
          fTOFSignalHe3 = trackP1->GetTOFsignal();
          fTOFSignalP = trackP2->GetTOFsignal();
          fTOFSignalPi = trackN->GetTOFsignal();
          fTOFSignalPi1 = trackN1->GetTOFsignal();
          fTOFSignalPi2 = trackN2->GetTOFsignal();
          //Track DeDx Sigma
          if (fBetheSplines) {
            fhe4DedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kAlpha);
            fhe3DedxSigma = fPID->NumberOfSigmasTPC(trackP1, AliPID::kHe3);
            fpDedxSigma = fPID->NumberOfSigmasTPC(trackP2, AliPID::kProton);
            fpiDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kPion);
            fpi1DedxSigma = fPID->NumberOfSigmasTPC(trackN1, AliPID::kPion);
            fpi2DedxSigma = fPID->NumberOfSigmasTPC(trackN2, AliPID::kPion);
          } else {
            fhe4DedxSigma = Bethe(*trackP, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
            fhe3DedxSigma = Bethe(*trackP1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
            fpDedxSigma = Bethe(*trackP2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
            fpiDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kPion);
            fpi1DedxSigma = fPID->NumberOfSigmasTPC(trackN1, AliPID::kPion);
            fpi2DedxSigma = fPID->NumberOfSigmasTPC(trackN2, AliPID::kPion);
          }
          //DCA From primary Vertex
          fhe4Dca = TMath::Abs(trackP->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
          fhe3Dca = TMath::Abs(trackP1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
          fpDca = TMath::Abs(trackP2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
          fpiDca = TMath::Abs(trackN->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
          fpi1Dca = TMath::Abs(trackN1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
          fpi2Dca = TMath::Abs(trackN2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
          //N Clusters TPC
          fhe4Ncls = trackP->GetTPCNcls();
          fhe3Ncls = trackP1->GetTPCNcls();
          fpNcls = trackP2->GetTPCNcls();
          fpiNcls = trackN->GetTPCNcls();
          fpi1Ncls = trackN1->GetTPCNcls();
          fpi2Ncls = trackN2->GetTPCNcls();
          //NClusters ITS
          fhe4NclsITS = trackP->GetNumberOfITSClusters();
          fhe3NclsITS = trackP1->GetNumberOfITSClusters();
          fpNclsITS = trackP2->GetNumberOfITSClusters();
          fpiNclsITS = trackN->GetNumberOfITSClusters();
          fpi1NclsITS = trackN1->GetNumberOfITSClusters();
          fpi2NclsITS = trackN2->GetNumberOfITSClusters();

          //Vertices
          TVector3 He4Vertex(fV0->Xv(), fV0->Yv(), fV0->Zv());
          AliESDVertex *He4Vert = new AliESDVertex(fV0->GetVertex());
          if(!trackP->PropagateToDCA(He4Vert, fMagneticField, 10)) continue;//10
          if(!trackN->PropagateToDCA(He4Vert, fMagneticField, 10)) continue;//10
          if(He4Vert) delete He4Vert;

          TVector3 He3Vertex(fV01->Xv(), fV01->Yv(), fV01->Zv());
          if(He4Vertex == He3Vertex) continue;
          AliESDVertex *He3Vert = new AliESDVertex(fV01->GetVertex());
          if(!trackP1->PropagateToDCA(He3Vert, fMagneticField, 10)) continue;//10
          if(!trackN1->PropagateToDCA(He3Vert, fMagneticField, 10)) continue;//10
          if(He3Vert) delete He3Vert;

          TVector3 PVertex(fV02->Xv(), fV02->Yv(), fV02->Zv());
          TVector3 Vertexdiff(TMath::Abs(He3Vertex.X() - PVertex.X()), TMath::Abs(He3Vertex.Y() - PVertex.Y()), 0);
          fVertDiff = Vertexdiff.Mag();
          //if(Vertexdiff.Mag() > 1) continue;
          AliESDVertex *PVert = new AliESDVertex(fV02->GetVertex());
          if(!trackP2->PropagateToDCA(PVert, fMagneticField, 10)) continue;//10
          if(!trackN2->PropagateToDCA(PVert, fMagneticField, 10)) continue;//10
          if(PVert) delete PVert;
          //V0 momenta
          TLorentzVector He4(0.,0.,0.,0.);
          He4.SetXYZM(2*trackP->Px(), 2*trackP->Py(), 2*trackP->Pz(), 3.929);
          TLorentzVector pi(0.,0.,0.,0.);
          pi.SetXYZM(trackN->Px(), trackN->Py(), trackN->Pz(), AliPID::ParticleMass(AliPID::kPion));
          TLorentzVector He4Mother = He4 + pi;
          fmHHe4Pi = He4Mother.M();
          fpHHe4Pi = He4Mother.P();
          fptHHe4Pi = He4Mother.Pt();
          fhe4P = He4.Pt();
          fpiP = pi.Pt();
          fyHHe4pi = He4Mother.Rapidity();

          TLorentzVector He3(0.,0.,0.,0.);
          He3.SetXYZM(2*trackP1->Px(), 2*trackP1->Py(), 2*trackP1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
          TLorentzVector pi1(0.,0.,0.,0.);
          pi1.SetXYZM(trackN1->Px(), trackN1->Py(), trackN1->Pz(), AliPID::ParticleMass(AliPID::kPion));
          TLorentzVector He3Mother = He3 + pi1; 
          fhe3P = He3.Pt();
          fpi1P = pi1.Pt();

          TLorentzVector Prot(0.,0.,0.,0.);
          Prot.SetXYZM(trackP2->Px(), trackP2->Py(), trackP2->Pz(), AliPID::ParticleMass(AliPID::kProton));
          TLorentzVector pi2(0.,0.,0.,0.);
          pi2.SetXYZM(trackN2->Px(), trackN2->Py(), trackN2->Pz(), AliPID::ParticleMass(AliPID::kPion));
          TLorentzVector PMother = Prot + pi2;
          fpP = Prot.Pt();
          fpi2P = pi2.Pt();

          TLorentzVector L4HeMother(0.,0.,0.,0.);
          L4HeMother = pi1 + He3 + Prot; 
          TLorentzVector L4HeMother1(0.,0.,0.,0.);
          L4HeMother1 = pi2 + He3 + Prot; 
          fmDaughterSum1 = L4HeMother.M();
          fmDaughterSum2 = L4HeMother1.M();
          fpSum1 = L4HeMother.P();
          fptSum1 = L4HeMother.Pt();
          fpSum2 = L4HeMother1.P();
          fptSum2 = L4HeMother1.Pt();
          fySum1 = L4HeMother.Rapidity();
          fySum2 = L4HeMother1.Rapidity();

          TVector3 secondaryVertex = He4Vertex;
          TVector3 secondaryVertex1 = He3Vertex;
          secondaryVertex = secondaryVertex - fPrimaryVertex;
          secondaryVertex1 = secondaryVertex1 - He4Vertex;
          fctHHe4Pi = secondaryVertex.Mag() * fmHHe4Pi / fpHHe4Pi;
          fctSum1 = secondaryVertex1.Mag() * fmDaughterSum1 / fpSum1;
          fctSum2 = secondaryVertex1.Mag() * fmDaughterSum2 / fpSum2;
          //Pointing Angle
          Double_t PA1 = He3Mother.Angle(-(He3Vertex - He4Vertex));
          fPAHe3He4 = TMath::Cos(PA1);

          Double_t PA2 = PMother.Angle(-(PVertex - He4Vertex));
          fPAPHe4 = TMath::Cos(PA2);

          Double_t PA3 = L4HeMother.Angle(-(He3Vertex - He4Vertex));
          fPA4LHe = TMath::Cos(PA3);

          Double_t PA4 = L4HeMother1.Angle(-(He3Vertex - He4Vertex));
          fPA4LHe1 = TMath::Cos(PA4);
          //DCA's
          Double_t DCA1 = trackP2->GetDCA(trackP1,fMagneticField,xthiss,xpp);
          fDcaHe3P = DCA1;

          Double_t DCA2 = trackN2->GetDCA(trackP1,fMagneticField,xthiss,xpp);
          fDcaHe3Pi2 = DCA2;

          Double_t DCA3 = trackP2->GetDCA(trackN1,fMagneticField,xthiss,xpp);
          fDcaPPi1 = DCA3;
          //DCA from sec/tert Vertices
          fhe4DcaSec = TMath::Abs(trackP->GetD(He4Vertex.X(), He4Vertex.Y(), fMagneticField));
          fhe3DcaTert = TMath::Abs(trackP1->GetD(He3Vertex.X(), He3Vertex.Y(), fMagneticField));
          fpDcaTert = TMath::Abs(trackP2->GetD(He3Vertex.X(), He3Vertex.Y(), fMagneticField));
          fpiDcaSec = TMath::Abs(trackN->GetD(He4Vertex.X(), He4Vertex.Y(), fMagneticField));
          fpi1DcaTert = TMath::Abs(trackN1->GetD(He3Vertex.X(), He3Vertex.Y(), fMagneticField));
          fpi2DcaTert = TMath::Abs(trackN2->GetD(He3Vertex.X(), He3Vertex.Y(), fMagneticField));
          //Armenteros Podolanski
          TVector3 vecN = pi.Vect();
          TVector3 vecP = He4.Vect();

          TVector3 vecM = He4Mother.Vect();

          fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
          fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

          farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
          farmpt = vecP.Mag()*sin(fthetaP);

          vecN = pi2.Vect();
          vecP = Prot.Vect();
          vecM = PMother.Vect();
          
          fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
          fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

          farmalpha1 = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
          farmpt1 = vecP.Mag()*sin(fthetaP);

          fNV0Cand = fNV0Cand + 1;

      }//end if helium4Positive && pionNegative
      fNumberV0s = (fNV0Cand);
      if (fNV0Cand) fTree->Fill();

      if (helium4Negative && pionPositive && helium3Negative && pionPositive1 && protonNegative && pionPositive2) {
        //BR for separation of different analysis: 1 = 4LLH --> (4LHe, pi), 2 = 4LLH -->(3LH, p, pi), 3 = Li4
        fBR = 1;
        fmc = -1; 
        fz = 1;
        fParticleSpecies = 1020010040;
        //V0 DCA
        fdcaHe4Pi = fV0->GetDcaV0Daughters();
        fdcaHe3Pi = fV01->GetDcaV0Daughters();
        fdcaPPi = fV02->GetDcaV0Daughters();
        //V0 Pointing-Angle
        fcosHe4Pi = fV0->GetV0CosineOfPointingAngle();
        fcosHe3Pi = fV01->GetV0CosineOfPointingAngle();
        fcosPPi = fV02->GetV0CosineOfPointingAngle();
        //V0-OnFlyStatus
        if(fV0->GetOnFlyStatus()) fonTheFly = 1;
        else if(!fV0->GetOnFlyStatus()) fonTheFly = 0;
        if(fV01->GetOnFlyStatus()) fonTheFly1 = 1;
        else if(!fV01->GetOnFlyStatus()) fonTheFly1 = 0;
        if(fV02->GetOnFlyStatus()) fonTheFly2 = 1;
        else if(!fV02->GetOnFlyStatus()) fonTheFly2 = 0;
        //Track TPC Signal
        fhe4Dedx = trackN->GetTPCsignal();
        fhe3Dedx = trackN1->GetTPCsignal();
        fpDedx = trackN2->GetTPCsignal();
        fpiDedx = trackP->GetTPCsignal();
        fpi1Dedx = trackP1->GetTPCsignal();
        fpi2Dedx = trackP2->GetTPCsignal();
        //Eta
        fEtaHe4 = trackN->Eta();
        fEtaHe3 = trackN1->Eta();
        fEtaP = trackN2->Eta();
        fEtaPi = trackP->Eta();
        fEtaPi1 = trackP1->Eta();
        fEtaPi2 = trackP2->Eta();
        //Phi
        fPhiHe4 = trackN->Phi();
        fPhiHe3 = trackN1->Phi();
        fPhiP = trackN2->Phi();
        fPhiPi = trackP->Phi();
        fPhiPi1 = trackP1->Phi();
        fPhiPi2 = trackP2->Phi();
        //GeoLength
        fGeoLengthHe4 = GeoLength(*trackN);
        fGeoLengthHe3 = GeoLength(*trackN1);
        fGeoLengthP = GeoLength(*trackN2);
        fGeoLengthPi = GeoLength(*trackP);
        fGeoLengthPi1 = GeoLength(*trackP1);
        fGeoLengthPi2 = GeoLength(*trackP2);
        //TOF Signal
        fTOFSignalHe4 = trackN->GetTOFsignal();
        fTOFSignalHe3 = trackN1->GetTOFsignal();
        fTOFSignalP = trackN2->GetTOFsignal();
        fTOFSignalPi = trackP->GetTOFsignal();
        fTOFSignalPi1 = trackP1->GetTOFsignal();
        fTOFSignalPi2 = trackP2->GetTOFsignal();
        //Track DeDx Sigma
        if (fBetheSplines) {
          fhe4DedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kAlpha);
          fhe3DedxSigma = fPID->NumberOfSigmasTPC(trackN1, AliPID::kHe3);
          fpDedxSigma = fPID->NumberOfSigmasTPC(trackN2, AliPID::kProton);
          fpiDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kPion);
          fpi1DedxSigma = fPID->NumberOfSigmasTPC(trackP1, AliPID::kPion);
          fpi2DedxSigma = fPID->NumberOfSigmasTPC(trackP2, AliPID::kPion);
        } else {
          fhe4DedxSigma = Bethe(*trackN, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
          fhe3DedxSigma = Bethe(*trackN1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
          fpDedxSigma = Bethe(*trackN2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
          fpiDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kPion);
          fpi1DedxSigma = fPID->NumberOfSigmasTPC(trackP1, AliPID::kPion);
          fpi2DedxSigma = fPID->NumberOfSigmasTPC(trackP2, AliPID::kPion);
        }
        //DCA From primary Vertex
        fhe4Dca = TMath::Abs(trackN->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
        fhe3Dca = TMath::Abs(trackN1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
        fpDca = TMath::Abs(trackN2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
        fpiDca = TMath::Abs(trackP->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
        fpi1Dca = TMath::Abs(trackP1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
        fpi2Dca = TMath::Abs(trackP2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
        //N Clusters TPC
        fhe4Ncls = trackN->GetTPCNcls();
        fhe3Ncls = trackN1->GetTPCNcls();
        fpNcls = trackN2->GetTPCNcls();
        fpiNcls = trackP->GetTPCNcls();
        fpi1Ncls = trackP1->GetTPCNcls();
        fpi2Ncls = trackP2->GetTPCNcls();
        //N ITS Cluster
        fhe4NclsITS = trackP->GetNumberOfITSClusters();
        fhe3NclsITS = trackP1->GetNumberOfITSClusters();
        fpNclsITS = trackP2->GetNumberOfITSClusters();
        fpiNclsITS = trackN->GetNumberOfITSClusters();
        fpi1NclsITS = trackN1->GetNumberOfITSClusters();
        fpi2NclsITS = trackN2->GetNumberOfITSClusters();
        //Vertices
        TVector3 He4Vertex(fV0->Xv(), fV0->Yv(), fV0->Zv());
        AliESDVertex *He4Vert = new AliESDVertex(fV0->GetVertex());
        if(!trackP->PropagateToDCA(He4Vert, fMagneticField, 10)) continue;//10
        if(!trackN->PropagateToDCA(He4Vert, fMagneticField, 10)) continue;//10
        if(He4Vert) delete He4Vert;

        TVector3 He3Vertex(fV01->Xv(), fV01->Yv(), fV01->Zv());
        if(He4Vertex == He3Vertex) continue;
        AliESDVertex *He3Vert = new AliESDVertex(fV01->GetVertex());
        if(!trackP1->PropagateToDCA(He3Vert, fMagneticField, 10)) continue;//10
        if(!trackN1->PropagateToDCA(He3Vert, fMagneticField, 10)) continue;//10
        if(He3Vert) delete He3Vert;

        TVector3 PVertex(fV02->Xv(), fV02->Yv(), fV02->Zv());
        TVector3 Vertexdiff(TMath::Abs(He3Vertex.X() - PVertex.X()), TMath::Abs(He3Vertex.Y() - PVertex.Y()), 0);
        fVertDiff = Vertexdiff.Mag();
        //if(Vertexdiff.Mag() > 1) continue;
        AliESDVertex *PVert = new AliESDVertex(fV02->GetVertex());
        if(!trackP2->PropagateToDCA(PVert, fMagneticField, 10)) continue;//10
        if(!trackN2->PropagateToDCA(PVert, fMagneticField, 10)) continue;//10
        if(PVert) delete PVert;
        //V0 momenta
        TLorentzVector He4(0.,0.,0.,0.);
        He4.SetXYZM(2*trackN->Px(), 2*trackN->Py(), 2*trackN->Pz(), 3.929);
        TLorentzVector pi(0.,0.,0.,0.);
        pi.SetXYZM(trackP->Px(), trackP->Py(), trackP->Pz(), AliPID::ParticleMass(AliPID::kPion));
        TLorentzVector He4Mother = He4 + pi;
        fmHHe4Pi = He4Mother.M();
        fpHHe4Pi = He4Mother.P();
        fptHHe4Pi = He4Mother.Pt();
        fhe4P = He4.Pt();
        fpiP = pi.Pt();
        fyHHe4pi = He4Mother.Rapidity();

        TLorentzVector He3(0.,0.,0.,0.);
        He3.SetXYZM(2*trackN1->Px(), 2*trackN1->Py(), 2*trackN1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
        TLorentzVector pi1(0.,0.,0.,0.);
        pi1.SetXYZM(trackP1->Px(), trackP1->Py(), trackP1->Pz(), AliPID::ParticleMass(AliPID::kPion));
        TLorentzVector He3Mother = He3 + pi1; 
        fhe3P = He3.Pt();
        fpi1P = pi1.Pt();

        TLorentzVector Prot(0.,0.,0.,0.);
        Prot.SetXYZM(trackN2->Px(), trackN2->Py(), trackN2->Pz(), AliPID::ParticleMass(AliPID::kProton));
        TLorentzVector pi2(0.,0.,0.,0.);
        pi2.SetXYZM(trackP2->Px(), trackP2->Py(), trackP2->Pz(), AliPID::ParticleMass(AliPID::kPion));
        TLorentzVector PMother = Prot + pi2;
        fpP = Prot.Pt();
        fpi2P = pi2.Pt();

        TLorentzVector L4HeMother(0.,0.,0.,0.);
        L4HeMother = pi1 + He3 + Prot; 
        TLorentzVector L4HeMother1(0.,0.,0.,0.);
        L4HeMother1 = pi2 + He3 + Prot; 
        fmDaughterSum1 = L4HeMother.M();
        fmDaughterSum2 = L4HeMother1.M();
        fpSum1 = L4HeMother.P();
        fptSum1 = L4HeMother.Pt();
        fpSum2 = L4HeMother1.P();
        fptSum2 = L4HeMother1.Pt();
        fySum1 = L4HeMother.Rapidity();
        fySum2 = L4HeMother1.Rapidity();
        //ct
        TVector3 secondaryVertex = He4Vertex;
        TVector3 secondaryVertex1 = He3Vertex;
        secondaryVertex = secondaryVertex - fPrimaryVertex;
        secondaryVertex1 = secondaryVertex1 - He4Vertex;
        fctHHe4Pi = secondaryVertex.Mag() * fmHHe4Pi / fpHHe4Pi;
        fctSum1 = secondaryVertex1.Mag() * fmDaughterSum1 / fpSum1;
        fctSum2 = secondaryVertex1.Mag() * fmDaughterSum2 / fpSum2;
        //PA
        Double_t PA1 = He3Mother.Angle(-(He3Vertex - He4Vertex));
        fPAHe3He4 = TMath::Cos(PA1);

        Double_t PA2 = PMother.Angle(-(PVertex - He4Vertex));
        fPAPHe4 = TMath::Cos(PA2);

        Double_t PA3 = L4HeMother.Angle(-(He3Vertex - He4Vertex));
        fPA4LHe = TMath::Cos(PA3);

        Double_t PA4 = L4HeMother1.Angle(-(He3Vertex - He4Vertex));
        fPA4LHe1 = TMath::Cos(PA4);
        //DCA
        Double_t DCA1 = trackN2->GetDCA(trackN1,fMagneticField,xthiss,xpp);
        fDcaHe3P = DCA1;

        Double_t DCA2 = trackP2->GetDCA(trackN1,fMagneticField,xthiss,xpp);
        fDcaHe3Pi2 = DCA2;

        Double_t DCA3 = trackN2->GetDCA(trackP1,fMagneticField,xthiss,xpp);
        fDcaPPi1 = DCA3;
        //DCA sec / tert Vertex
        fhe4DcaSec = TMath::Abs(trackN->GetD(He4Vertex.X(), He4Vertex.Y(), fMagneticField));
        fhe3DcaTert = TMath::Abs(trackN1->GetD(He3Vertex.X(), He3Vertex.Y(), fMagneticField));
        fpDcaTert = TMath::Abs(trackN2->GetD(He3Vertex.X(), He3Vertex.Y(), fMagneticField));
        fpiDcaSec = TMath::Abs(trackP->GetD(He4Vertex.X(), He4Vertex.Y(), fMagneticField));
        fpi1DcaTert = TMath::Abs(trackP1->GetD(He3Vertex.X(), He3Vertex.Y(), fMagneticField));
        fpi2DcaTert = TMath::Abs(trackP2->GetD(He3Vertex.X(), He3Vertex.Y(), fMagneticField));
        //Armenteros Podolanski
        TVector3 vecP = pi.Vect();
        TVector3 vecN = He4.Vect();

        TVector3 vecM = He4Mother.Vect();

        fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
        fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

        farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
        farmpt = vecP.Mag()*sin(fthetaP);

        vecP = pi2.Vect();
        vecN = Prot.Vect();
        vecM = PMother.Vect();
        
        fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
        fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

        farmalpha1 = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
        farmpt1 = vecP.Mag()*sin(fthetaP);

        fNV0Cand = fNV0Cand + 1;
        }//end helium4Negativ && pionPositive 
        fNumberV0s = (fNV0Cand);
        if (fNV0Cand) fTree->Fill();    
      }
    }
  }
}

void AliAnalysisTaskDoubleHypNucTree::TrackAnalysis(AliESDtrackCuts trackCutsSec, AliESDtrackCuts trackCutsTert, AliESDtrackCuts trackCutsPi, AliVertexerTracks *vertexer, AliVertexerTracks *vertexer1, Double_t *dn, Double_t *dd, Double_t xthiss, Double_t xpp, Bool_t ITSClusterCut, Int_t SecClusters, Int_t TertClusters){
  //Counter
  Int_t found = 0;
  //LorentzVectors
  TLorentzVector *lorentzVectorHe3 = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *lorentzVectorPionMinus = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *lorentzVectorProton = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *lorentzVectorPionMinus2 = new TLorentzVector(0.,0.,0.,0.);
  TLorentzVector *lorentzVectorDoubleHyperHydrogen4= new TLorentzVector(0.,0.,0.,0.);
  TVector3 *h = new TVector3(0., 0., 0.);
  //++++++++++++++++++++++++++++ 4LHe Track ++++++++++++++++++++++++++++++++++//
  for (Int_t mTracks = 0; mTracks < fESDevent->GetNumberOfTracks(); mTracks++) {

    AliESDtrack* track5 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(mTracks));

    if (!track5->GetInnerParam()) continue;

    if (!trackCutsSec.AcceptTrack(track5)) continue;
    //ITS Cluster Cut defined in User exec
    if(ITSClusterCut && track5->GetNumberOfITSClusters() <= SecClusters) continue;

    Double_t ptot5 = track5->GetInnerParam()->GetP();
    Double_t sign5 = track5->GetSign();    

    fHistdEdx->Fill(ptot5*sign5, track5->GetTPCsignal());

    if(fBetheSplines){
      if(TMath::Abs(fPID->NumberOfSigmasTPC(track5, AliPID::kAlpha)) > 3) continue;
    }
    else {
      if(TMath::Abs(Bethe(*track5, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe)) > 3) continue;
    }       
    //++++++++++++++++++++++++++++ 3He Track ++++++++++++++++++++++++++++++++++//
    for (Int_t iTracks = mTracks + 1; iTracks < fESDevent->GetNumberOfTracks(); iTracks++) {

      AliESDtrack* track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(iTracks));

      if (!track1->GetInnerParam()) continue;

      if (!trackCutsTert.AcceptTrack(track1)) continue;
      //ITS Cluster Cut defined in User exec
      if(ITSClusterCut && track1->GetNumberOfITSClusters() >= TertClusters) continue;

      if(iTracks == mTracks) continue;     

      Double_t ptot1 = track1->GetInnerParam()->GetP();
      Double_t sign1 = track1->GetSign();

      if (sign5 != sign1) continue;

      fHistdEdx->Fill(ptot1*sign1, track1->GetTPCsignal());

      if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) > 3) continue;
      }
      else {
        if(TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) > 3) continue;
      }

      //++++++++++++++++++++++++++++++++ p Track +++++++++++++++++++++++++++++++++++++++//
      for (Int_t jTracks = iTracks + 1; jTracks < fESDevent->GetNumberOfTracks(); jTracks++) {

        AliESDtrack* track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(jTracks));              

        if (!track2->GetInnerParam()) continue;

        if (!trackCutsTert.AcceptTrack(track2)) continue;
        //ITS Cluster Cut defined in User exec
        if(ITSClusterCut && track2->GetNumberOfITSClusters() >= TertClusters) continue;

        Double_t ptot2 = track2->GetInnerParam()->GetP();
        Double_t sign2 = track2->GetSign();

        if(mTracks == jTracks || iTracks == jTracks) continue;

        if(sign1 != sign2) continue; 

        fHistdEdx->Fill(ptot2*sign2, track2->GetTPCsignal());

        if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) > 3) continue;
        }
        else {
          if(TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) > 3) continue;
        }

        //++++++++++++++++++++++++++++++++ pi tert Track +++++++++++++++++++++++++++++++++++++++//
        for (Int_t kTracks = jTracks + 1; kTracks < fESDevent->GetNumberOfTracks(); kTracks++) {

          AliESDtrack* track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(kTracks));        

          if (!track3->GetInnerParam()) continue;

          if (!trackCutsPi.AcceptTrack(track3)) continue;
          //ITS Cluster Cut defined in User exec
          if(ITSClusterCut && track3->GetNumberOfITSClusters() >= TertClusters) continue;

          Double_t ptot3 = track3->GetInnerParam()->GetP();
          Double_t sign3 = track3->GetSign();

          if (sign2 == sign3) continue;

          if (jTracks == kTracks || iTracks == kTracks || mTracks == kTracks) continue; //reject using the same track twice

          fHistdEdx->Fill(ptot3*sign3, track3->GetTPCsignal());

          if(TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) > 3) continue;
 
          //++++++++++++++++++++++++++++++++ pi sec Track +++++++++++++++++++++++++++++++++++++++//
          for (Int_t lTracks = kTracks + 1; lTracks < fESDevent->GetNumberOfTracks(); lTracks++) {

            AliESDtrack* track4 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(lTracks));//Track for Double Hyper Hydrogen 4 Decay (Pion Track)

            if (!track4->GetInnerParam()) continue;

            if (!trackCutsPi.AcceptTrack(track4)) continue;
            //ITS Cluster Cut defined in User exec
            if(ITSClusterCut && track4->GetNumberOfITSClusters() <= SecClusters) continue;

            Double_t ptot4 = track4->GetInnerParam()->GetP();
            Double_t sign4 = track4->GetSign();

            if (sign4 != sign3) continue;// Data + + - -

            if (jTracks == lTracks || iTracks == lTracks || kTracks == lTracks || mTracks == lTracks) continue; //reject using the same track twice          

            fHistdEdx->Fill(ptot4*sign4, track4->GetTPCsignal());

              if(TMath::Abs(fPID->NumberOfSigmasTPC(track4, AliPID::kPion)) > 3) continue;
                  
            //--------------------------------------------------------------------------------------------------//
            if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 3)){
              if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3)) {
                if (TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) < 3) {
                  if (TMath::Abs(fPID->NumberOfSigmasTPC(track4, AliPID::kPion)) < 3) {
                    if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track5, AliPID::kAlpha)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track5, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe)) < 3)){      

                      //======= Vertex Reconstruction =======
                      TObjArray *trkArray = new TObjArray(3);
                      TObjArray *trkArray1 = new TObjArray(2);

                      trkArray->AddAt(track1,0);
                      trkArray->AddAt(track2,1);
                      trkArray->AddAt(track3,2);

                      trkArray1->AddAt(track5,0);
                      trkArray1->AddAt(track4,1);
                      AliESDVertex *secVertex1 = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray1);
                      if(trkArray1) delete trkArray1;

                      AliESDVertex *esdVer2 = new AliESDVertex(*secVertex1);     
                      vertexer1->SetVtxStart(esdVer2);
                      if(esdVer2) delete esdVer2;
                      AliESDVertex *tertVertex1 = (AliESDVertex*)vertexer1->VertexForSelectedESDTracks(trkArray);
                      if(trkArray) delete trkArray;
                      //HyperHelium4 Decay-->tertVtx
                      if(!track1->PropagateToDCA(tertVertex1, fMagneticField, 10)) continue;
                      if(!track2->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;
                      if(!track3->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;
                      //DoubleHyperHydrogen4 Decay-->seVtx
                      if(!track4->PropagateToDCA(secVertex1, fMagneticField, 10))continue;
                      if(!track5->PropagateToDCA(secVertex1, fMagneticField, 10))continue;
                      //DCAs
                      Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
                      fDcaHe3P = *dca;

                      Double_t *dca1 = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
                      fDcaHe3Pi2 = *dca1;

                      Double_t *dca3 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
                      fDcaPPi1 = *dca3;

                      Double_t *dca5 = new Double_t(track4->GetDCA(track5,fMagneticField,xthiss,xpp));
                      fdcaHe4Pi = *dca5;

                      if(dca) delete dca;
                      if(dca1) delete dca1;
                      if(dca3) delete dca3;
                      if(dca5) delete dca5;
                      //PA Array
                      dd[0]=dn[0]-secVertex1->GetX();
                      dd[1]=dn[1]-secVertex1->GetY();
                      dd[2]=dn[2]-secVertex1->GetZ();
                      //PA Vectors
                      TVector3 g(0.,0.,0.);
                      h->SetXYZ(-dd[0],-dd[1],-dd[2]);
                      g.SetXYZ(-(secVertex1->GetX()-tertVertex1->GetX()), -(secVertex1->GetY()-tertVertex1->GetY()), -(secVertex1->GetZ()-tertVertex1->GetZ()));

                      TVector3 secVtx(secVertex1->GetX(),secVertex1->GetY(),secVertex1->GetZ());
                      TVector3 tertVtx(tertVertex1->GetX(),tertVertex1->GetY(),tertVertex1->GetZ());

                      if(secVertex1) delete secVertex1;
                      if(tertVertex1) delete tertVertex1;

                      //============================ 4LLH positive ===========================
                      if(sign1>0 && sign2>0 && sign3<0 && sign4 <0 && sign5>0){
                        //BR for separation of different analysis: 1 = 4LLH --> (4LHe, pi), 2 = 4LLH -->(3LH, p, pi), 3 = Li4
                        fBR = 1;
                        fmc = -1; 
                        fz = 1;
                        fParticleSpecies = 1020010040;

                        lorentzVectorHe3->SetXYZM(2*track1->Px(),2*track1->Py(),2*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
                        lorentzVectorPionMinus->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
                        lorentzVectorProton->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
                        lorentzVectorPionMinus2->SetXYZM(track4->Px(),track4->Py(),track4->Pz(),AliPID::ParticleMass(AliPID::kPion));
                        TLorentzVector Helium4(2*track5->Px(), 2*track5->Py(), 2*track5->Pz(), 3.929);

                        *lorentzVectorDoubleHyperHydrogen4   = *lorentzVectorHe3 + *lorentzVectorPionMinus + *lorentzVectorProton + *lorentzVectorPionMinus2;    

                        TLorentzVector HHe4PiMother = Helium4 + *lorentzVectorPionMinus2; 
                        //Daughter momenta
                        fhe4P = Helium4.Pt();
                        fpiP = lorentzVectorPionMinus2->Pt();          
                        fhe3P =lorentzVectorHe3->Pt();
                        fpi1P = lorentzVectorPionMinus->Pt();
                        fpP = lorentzVectorProton->Pt();
                        //N Cls TPC
                        fhe4Ncls = track5->GetTPCNcls();
                        fhe3Ncls = track1->GetTPCNcls();
                        fpNcls = track2->GetTPCNcls();
                        fpiNcls = track4->GetTPCNcls();
                        fpi1Ncls = track3->GetTPCNcls();
                        //NCls ITS
                        fhe4NclsITS = track5->GetNumberOfITSClusters();
                        fhe3NclsITS = track1->GetNumberOfITSClusters();
                        fpNclsITS = track2->GetNumberOfITSClusters();
                        fpiNclsITS = track4->GetNumberOfITSClusters();
                        fpi1NclsITS = track3->GetNumberOfITSClusters();
                        //dEdx Tracks
                        fhe4Dedx = track5->GetTPCsignal();
                        fhe3Dedx = track1->GetTPCsignal();
                        fpDedx = track2->GetTPCsignal();
                        fpiDedx = track4->GetTPCsignal();
                        fpi1Dedx = track3->GetTPCsignal();
                        //Eta
                        fEtaHe4 = track5->Eta();
                        fEtaHe3 = track1->Eta();
                        fEtaP = track2->Eta();
                        fEtaPi = track3->Eta();
                        fEtaPi1 = track4->Eta();
                        //Phi
                        fPhiHe4 = track5->Phi();
                        fPhiHe3 = track1->Phi();
                        fPhiP = track2->Phi();
                        fPhiPi = track3->Phi();
                        fPhiPi1 = track4->Phi();
                        //GeoLength
                        fGeoLengthHe4 = GeoLength(*track5);
                        fGeoLengthHe3 = GeoLength(*track1);
                        fGeoLengthP = GeoLength(*track2);
                        fGeoLengthPi = GeoLength(*track3);
                        fGeoLengthPi1 = GeoLength(*track4);
                        //TOF Signal
                        fTOFSignalHe4 = track5->GetTOFsignal();
                        fTOFSignalHe3 = track1->GetTOFsignal();
                        fTOFSignalP = track2->GetTOFsignal();
                        fTOFSignalPi = track3->GetTOFsignal();
                        fTOFSignalPi1 = track4->GetTOFsignal();
                        //dEdx Sigma
                        if (fBetheSplines) {
                          fhe4DedxSigma = fPID->NumberOfSigmasTPC(track5, AliPID::kAlpha);
                          fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
                          fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
                          fpiDedxSigma = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
                          fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                        } else {
                          fhe4DedxSigma = Bethe(*track5, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
                          fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
                          fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
                          fpiDedxSigma = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
                          fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                        }
                        //DCA prim Vertex
                        fhe4Dca = TMath::Abs(track5->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                        fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                        fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                        fpiDca = TMath::Abs(track4->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                        fpi1Dca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField)); 
                        //Mother M (daughter sum)
                        fmDaughterSum1 = lorentzVectorDoubleHyperHydrogen4->M();
                        fpSum1 = lorentzVectorDoubleHyperHydrogen4->P();
                        fptSum1 = lorentzVectorDoubleHyperHydrogen4->Pt();
                        fySum1 = lorentzVectorDoubleHyperHydrogen4->Rapidity();
                        //Mother M (4LHe, pi - tracks)
                        fmHHe4Pi = HHe4PiMother.M();
                        fpHHe4Pi = HHe4PiMother.P();
                        fptHHe4Pi = HHe4PiMother.Pt();
                        fyHHe4pi = HHe4PiMother.Rapidity();
                        //DCA sec / tert Vertex
                        fhe4DcaSec = TMath::Abs(track5->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
                        fhe3DcaTert = TMath::Abs(track1->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                        fpDcaTert = TMath::Abs(track2->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                        fpiDcaSec = TMath::Abs(track4->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
                        fpi1DcaTert = TMath::Abs(track3->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                        //4LHe Vect
                        TLorentzVector HHe4(0.,0.,0.,0.);
                        HHe4=*lorentzVectorHe3 + *lorentzVectorPionMinus + *lorentzVectorProton;
                        //Mother M 4LHe
                        fmDaughterSum2 = HHe4.M();
                        fpSum2 = HHe4.P();
                        fptSum2 = HHe4.Pt();
                        fySum2 = HHe4.Rapidity();
                        //PA
                        Double_t pointingAngle = lorentzVectorDoubleHyperHydrogen4->Angle(*h);
                        fPA4LHe = TMath::Cos(pointingAngle);
                        fPA4LHe1 = TMath::Cos(HHe4.Angle(g));
                        //ct
                        fctHHe4Pi = h->Mag() * fmHHe4Pi / fpHHe4Pi;
                        fctSum1 = h->Mag() * fmDaughterSum1 / fpSum1;
                        fctSum2 = g.Mag() * fmDaughterSum2 / fpSum2;
                        //Armenteros Podolanski
                        TVector3 vecP = lorentzVectorPionMinus2->Vect();
                        TVector3 vecN = Helium4.Vect();

                        TVector3 vecM = HHe4PiMother.Vect();

                        fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                        fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                        farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                        farmpt = vecP.Mag()*sin(fthetaP);

                        vecN = lorentzVectorPionMinus->Vect();
                        vecP = lorentzVectorProton->Vect();
                        vecM = vecN + vecP;

                        fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                        fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                        farmalpha1 = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                        farmpt1 = vecP.Mag()*sin(fthetaP);

                        found++;

                      }//end DoubleHyperHydrogen4
                      //======================= 4LLH negative ================================
                      else if(sign1<0 && sign2<0 && sign3>0 && sign4>0 && sign5<0){
                        //BR for separation of different analysis: 1 = 4LLH --> (4LHe, pi), 2 = 4LLH -->(3LH, p, pi), 3 = Li4
                        fBR = 1;
                        fmc = -1; 
                        fz = -1;
                        fParticleSpecies = -1020010040;

                        lorentzVectorHe3->SetXYZM(2*track1->Px(),2*track1->Py(),2*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
                        lorentzVectorPionMinus->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
                        lorentzVectorProton->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
                        lorentzVectorPionMinus2->SetXYZM(track4->Px(),track4->Py(),track4->Pz(),AliPID::ParticleMass(AliPID::kPion));
                        TLorentzVector Helium4(2*track5->Px(), 2*track5->Py(), 2*track5->Pz(), 3.929);
                        *lorentzVectorDoubleHyperHydrogen4   = *lorentzVectorHe3 + *lorentzVectorPionMinus + *lorentzVectorProton + *lorentzVectorPionMinus2;    
                        TLorentzVector HHe4PiMother = Helium4 + *lorentzVectorPionMinus2; 
                        //Daughter P
                        fhe4P = Helium4.Pt();
                        fpiP = lorentzVectorPionMinus2->Pt();          
                        fhe3P =lorentzVectorHe3->Pt();
                        fpi1P = lorentzVectorPionMinus->Pt();
                        fpP = lorentzVectorProton->Pt();
                        //NCls TPC
                        fhe4Ncls = track5->GetTPCNcls();
                        fhe3Ncls = track1->GetTPCNcls();
                        fpNcls = track2->GetTPCNcls();
                        fpiNcls = track4->GetTPCNcls();
                        fpi1Ncls = track3->GetTPCNcls();
                        //NCls ITS
                        fhe4NclsITS = track5->GetNumberOfITSClusters();
                        fhe3NclsITS = track1->GetNumberOfITSClusters();
                        fpNclsITS = track2->GetNumberOfITSClusters();
                        fpiNclsITS = track4->GetNumberOfITSClusters();
                        fpi1NclsITS = track3->GetNumberOfITSClusters();
                        //track dEdx
                        fhe4Dedx = track5->GetTPCsignal();
                        fhe3Dedx = track1->GetTPCsignal();
                        fpDedx = track2->GetTPCsignal();
                        fpiDedx = track4->GetTPCsignal();
                        fpi1Dedx = track3->GetTPCsignal();
                        //Eta
                        fEtaHe4 = track5->Eta();
                        fEtaHe3 = track1->Eta();
                        fEtaP = track2->Eta();
                        fEtaPi = track3->Eta();
                        fEtaPi1 = track4->Eta();
                        //Phi
                        fPhiHe4 = track5->Phi();
                        fPhiHe3 = track1->Phi();
                        fPhiP = track2->Phi();
                        fPhiPi = track3->Phi();
                        fPhiPi1 = track4->Phi();
                        //GeoLength
                        fGeoLengthHe4 = GeoLength(*track5);
                        fGeoLengthHe3 = GeoLength(*track1);
                        fGeoLengthP = GeoLength(*track2);
                        fGeoLengthPi = GeoLength(*track3);
                        fGeoLengthPi1 = GeoLength(*track4);
                        //TOF Signal
                        fTOFSignalHe4 = track5->GetTOFsignal();
                        fTOFSignalHe3 = track1->GetTOFsignal();
                        fTOFSignalP = track2->GetTOFsignal();
                        fTOFSignalPi = track3->GetTOFsignal();
                        fTOFSignalPi1 = track4->GetTOFsignal();
                        //dEdx sigma
                        if (fBetheSplines) {
                          fhe4DedxSigma = fPID->NumberOfSigmasTPC(track5, AliPID::kAlpha);
                          fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
                          fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
                          fpiDedxSigma = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
                          fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                        } else {
                          fhe4DedxSigma = Bethe(*track5, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
                          fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
                          fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
                          fpiDedxSigma = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
                          fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                        }
                        //DCA prim Vertex
                        fhe4Dca = TMath::Abs(track5->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                        fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                        fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                        fpiDca = TMath::Abs(track4->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                        fpi1Dca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField)); 
                        //Mother M (daughter sum)
                        fmDaughterSum1 = lorentzVectorDoubleHyperHydrogen4->M();
                        fpSum1 = lorentzVectorDoubleHyperHydrogen4->P();
                        fptSum1 = lorentzVectorDoubleHyperHydrogen4->Pt();
                        fySum1 = lorentzVectorDoubleHyperHydrogen4->Rapidity();
                        //Mother M (4LHe, pi- track)
                        fmHHe4Pi = HHe4PiMother.M();
                        fpHHe4Pi = HHe4PiMother.P();
                        fptHHe4Pi = HHe4PiMother.Pt();
                        fyHHe4pi = HHe4PiMother.Rapidity();
                        //DCA sec / tert Vtx
                        fhe4DcaSec = TMath::Abs(track5->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
                        fhe3DcaTert = TMath::Abs(track1->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                        fpDcaTert = TMath::Abs(track2->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                        fpiDcaSec = TMath::Abs(track4->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
                        fpi1DcaTert = TMath::Abs(track3->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                        //4LHe Vec
                        TLorentzVector HHe4(0.,0.,0.,0.);
                        HHe4=*lorentzVectorHe3 + *lorentzVectorPionMinus + *lorentzVectorProton;
                        //4LHe m, p, pt, y
                        fmDaughterSum2 = HHe4.M();
                        fpSum2 = HHe4.P();
                        fptSum2 = HHe4.Pt();
                        fySum2 = HHe4.Rapidity();
                        //PA
                        Double_t pointingAngle = lorentzVectorDoubleHyperHydrogen4->Angle(*h);
                        fPA4LHe = TMath::Cos(pointingAngle);
                        fPA4LHe1 = TMath::Cos(HHe4.Angle(g));
                        //ct
                        fctHHe4Pi = h->Mag() * fmHHe4Pi / fpHHe4Pi;
                        fctSum1 = h->Mag() * fmDaughterSum1 / fpSum1;
                        fctSum2 = g.Mag() * fmDaughterSum2 / fpSum2;
                        //Armenteros Podolanski
                        TVector3 vecP = lorentzVectorPionMinus2->Vect();
                        TVector3 vecN = Helium4.Vect();

                        TVector3 vecM = HHe4PiMother.Vect();

                        fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                        fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                        farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                        farmpt = vecP.Mag()*sin(fthetaP);

                        vecP = lorentzVectorPionMinus->Vect();
                        vecN = lorentzVectorProton->Vect();
                        vecM = vecN + vecP;

                        fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                        fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                        farmalpha1 = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                        farmpt1 = vecP.Mag()*sin(fthetaP);

                        found++;
                      }//end AntiDoubleHyperHydrogen4 
                    }//PID Check
                  }//PID Check
                }//PID Check
              }//PID Check
            }//End PID Check 
          }//end loop track4
          if(found>0) fTree->Fill();
          found = 0; 
        }//end loop track3
      }//end loop track2
    }//end loop track1
  }//end loop track5
}

void AliAnalysisTaskDoubleHypNucTree::CombinedAnalysis(AliESDtrackCuts trackCutsV0, AliESDtrackCuts trackCutsTert, AliESDtrackCuts trackCutsPi, AliVertexerTracks *vertexer, AliVertexerTracks *vertexer1, Double_t *dn, Double_t *dd, Double_t xthiss, Double_t xpp, Bool_t ITSClusterCut, Int_t SecClusters, Int_t TertClusters){
  TLorentzVector *He3 = new TLorentzVector(0., 0., 0., 0.);
  TLorentzVector *Prot = new TLorentzVector(0., 0. ,0. ,0.);
  TLorentzVector *pi1 = new TLorentzVector(0., 0. ,0. ,0.);
  TLorentzVector *pi = new TLorentzVector(0., 0. ,0. ,0.);
  TLorentzVector *He4 = new TLorentzVector(0., 0. ,0. ,0.);
  TLorentzVector *He41 = new TLorentzVector(0., 0. ,0. ,0.);
  TLorentzVector *HypTrit = new TLorentzVector(0., 0. ,0. ,0.);
    //+++++++++++++++++++++++++++V0(4LHe, pi)++++++++++++++++++++++++++++++++++
  for (Int_t ivertex = 0; ivertex < fESDevent->GetNumberOfV0s(); ivertex++) {

    fHistV0->Fill(fTrigger);
    fV0 = fESDevent->GetV0(ivertex);

    Bool_t v0ChargeCorrect = kTRUE;
    AliESDtrack* trackN = fESDevent->GetTrack(fV0->GetIndex(0));
    AliESDtrack* trackP = fESDevent->GetTrack(fV0->GetIndex(1));

    if (trackN->GetSign() > 0 ) {
      trackN = fESDevent->GetTrack(fV0->GetIndex(1));
      trackP = fESDevent->GetTrack(fV0->GetIndex(0));
      v0ChargeCorrect = kFALSE;
    }

    Double_t sign5 = trackP->GetSign();
    Double_t sign4 = trackN->GetSign();

    if (!trackCutsV0.AcceptTrack(trackN)) continue;
    if (!trackCutsV0.AcceptTrack(trackP)) continue;

    fHistdEdxV0->Fill(trackP->GetInnerParam()->GetP() * trackP->GetSign(), trackP->GetTPCsignal());
    fHistdEdxV0->Fill(trackN->GetInnerParam()->GetP() * trackN->GetSign(), trackN->GetTPCsignal());

    if(fPIDCheckOnly) continue;

    //ITS Cluster Cut defined in User exec
    if(ITSClusterCut){
      if(trackP->GetNumberOfITSClusters() >= SecClusters || trackN->GetNumberOfITSClusters() >= SecClusters) continue;
    }
    //two different BR!
    Bool_t pionPositive     = kFALSE;
    Bool_t pionNegative     = kFALSE;
    Bool_t helium4Positive  = kFALSE;
    Bool_t helium4Negative  = kFALSE;
    Bool_t helium3Positive  = kFALSE;
    Bool_t helium3Negative  = kFALSE;

    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kPion)) < 3) {
      pionPositive = kTRUE;
    }
    else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kPion)) < 3) {
      pionNegative = kTRUE;
    }
    else continue;
    //Use Framework Splines for p, He3, Alpha
    if (fBetheSplines) {
      if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kAlpha)) < 4) {
        helium4Positive = kTRUE;
      } else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kAlpha)) < 4) {
        helium4Negative = kTRUE;
      }
      else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kHe3)) < 4) {
        helium3Positive = kTRUE;
      } else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kHe3)) < 4) {
        helium3Negative = kTRUE;
      }
      else continue;
    //Use own Splines for p, He3, Alpha
    } else {

      if (TMath::Abs(Bethe(*trackP, AliPID::ParticleMass(AliPID::kAlpha),  2, fBetheParamsHe)) < 3) {
        helium4Positive = kTRUE;
      } else if (TMath::Abs(Bethe(*trackN, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe)) < 3) {
        helium4Negative = kTRUE;
      }
      else if (TMath::Abs(Bethe(*trackP, AliPID::ParticleMass(AliPID::kHe3),  2, fBetheParamsHe)) < 3) {
        helium3Positive = kTRUE;
      } else if (TMath::Abs(Bethe(*trackN, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 3) {
        helium3Negative = kTRUE;
      }
      else continue;
    }
    //++++ 4LLH --> 4LHe + pi (pos) ++++//
    if (helium4Positive && pionNegative) {
      //++++++++++++++++++++++++++++++ He3 Track +++++++++++++++++++++++++++++++++++
      for (Int_t iTracks = 0; iTracks < fESDevent->GetNumberOfTracks(); iTracks++) {

        AliESDtrack* track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(iTracks));

        if (!track1->GetInnerParam()) continue;

        if (!trackCutsTert.AcceptTrack(track1)) continue;
        //ITS Cluster Cut defined in User Exec
        if(ITSClusterCut && track1->GetNumberOfITSClusters() >= TertClusters) continue;  

        Double_t ptot1 = track1->GetInnerParam()->GetP();
        Double_t sign1 = track1->GetSign();

        if(sign1 < 0) continue;

        fHistdEdx->Fill(ptot1*sign1, track1->GetTPCsignal());

        if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) > 3) continue;
        }
        else {
          if(TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) > 3) continue;
        }

        //+++++++++++++++++++++++++++++++++++ p Track ++++++++++++++++++++++++++++++++++++++++++
        for (Int_t jTracks = iTracks + 1; jTracks < fESDevent->GetNumberOfTracks(); jTracks++) {

          AliESDtrack* track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(jTracks));              

          if (!track2->GetInnerParam()) continue;

          if (!trackCutsTert.AcceptTrack(track2)) continue;
          //ITS Cluster Cut defined in User Exec
          if(ITSClusterCut && track2->GetNumberOfITSClusters() >= TertClusters) continue;

          Double_t ptot2 = track2->GetInnerParam()->GetP();
          Double_t sign2 = track2->GetSign();

          if(iTracks == jTracks) continue;

          if(sign2 < 0) continue; 

          fHistdEdx->Fill(ptot2*sign2, track2->GetTPCsignal());

        if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) > 3) continue;
        }
        else {
          if(TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) > 3) continue;
        }

          //+++++++++++++++++++++++++++++++++++ pi Track +++++++++++++++++++++++++++++++++++++++++
          for (Int_t kTracks = jTracks + 1; kTracks < fESDevent->GetNumberOfTracks(); kTracks++) {

            AliESDtrack* track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(kTracks));        

            if (!track3->GetInnerParam()) continue;

            if (!trackCutsPi.AcceptTrack(track3)) continue;
            //ITS Cluster Cut defined in User Exec
            if(ITSClusterCut && track3->GetNumberOfITSClusters() >= TertClusters) continue;

            Double_t ptot3 = track3->GetInnerParam()->GetP();
            Double_t sign3 = track3->GetSign();

            if (sign3 > 0) continue;

            if (jTracks == kTracks || iTracks == kTracks) continue; //reject using the same track twice

            fHistdEdx->Fill(ptot3*sign3, track3->GetTPCsignal());

            if(TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) > 3) continue;
            
            //------------------------------------------------------------------------------------------------//
            if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 3)){
              if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3)) {
                if (TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) < 3) {

                  //======= Vertex Reconstruction ======
                  TObjArray *trkArray = new TObjArray(3);

                  trkArray->AddAt(track1,0);
                  trkArray->AddAt(track2,1);
                  trkArray->AddAt(track3,2);

                  TVector3 He4Vertex(fV0->Xv(), fV0->Yv(), fV0->Zv());
                  AliESDVertex *He4Vert = new AliESDVertex(fV0->GetVertex());
                  vertexer1->SetVtxStart(He4Vert);
                  //V0 = 4LHe, pi --> sec Vertex
                  if(!trackP->PropagateToDCA(He4Vert, fMagneticField, 10)) continue;
                  if(!trackN->PropagateToDCA(He4Vert, fMagneticField, 10)) continue;
                  if(He4Vert) delete He4Vert;

                  AliESDVertex *tertVertex1 = (AliESDVertex*)vertexer1->VertexForSelectedESDTracks(trkArray);
                  if(trkArray) delete trkArray;
                  //HyperHelium4 Decay --> tert Vtx
                  if(!track1->PropagateToDCA(tertVertex1, fMagneticField, 10)) continue;
                  if(!track2->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;
                  if(!track3->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;
                  //DCA
                  Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
                  fDcaHe3P = *dca;

                  Double_t *dca1 = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
                  fDcaHe3Pi2 = *dca1;

                  Double_t *dca3 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
                  fDcaPPi1 = *dca3;

                  if(dca) delete dca;
                  if(dca1) delete dca1;
                  if(dca3) delete dca3;

                  TVector3 tertVtx(tertVertex1->GetX(),tertVertex1->GetY(),tertVertex1->GetZ());

                  if(tertVertex1) delete tertVertex1;

                  //===================================================================================================//
                  //BR for separation of different analysis: 1 = 4LLH --> (4LHe, pi), 2 = 4LLH -->(3LH, p, pi), 3 = Li4
                  fBR = 1;
                  fmc = -1; 
                  fz = 1;
                  fParticleSpecies = 1020010040;
                  //V0 DCA
                  fdcaHe4Pi = fV0->GetDcaV0Daughters();
                  //V0 Pointing-Angle
                  fcosHe4Pi = fV0->GetV0CosineOfPointingAngle();
                  //V0-OnFlyStatus
                  if(fV0->GetOnFlyStatus()) fonTheFly = 1;
                  else if(!fV0->GetOnFlyStatus()) fonTheFly = 0;
                  //Track TPC Signal
                  fhe4Dedx = trackP->GetTPCsignal();
                  fhe3Dedx = track1->GetTPCsignal();
                  fpDedx = track2->GetTPCsignal();
                  fpiDedx = trackN->GetTPCsignal();
                  fpi1Dedx = track3->GetTPCsignal();
                  //Eta
                  fEtaHe4 = trackP->Eta();
                  fEtaHe3 = track1->Eta();
                  fEtaP = track2->Eta();
                  fEtaPi = trackN->Eta();
                  fEtaPi1 = track3->Eta();
                  //Phi
                  fPhiHe4 = trackP->Phi();
                  fPhiHe3 = track1->Phi();
                  fPhiP = track2->Phi();
                  fPhiPi = trackN->Phi();
                  fPhiPi1 = track3->Phi();
                  //GeoLength
                  fGeoLengthHe4 = GeoLength(*trackP);
                  fGeoLengthHe3 = GeoLength(*track1);
                  fGeoLengthP = GeoLength(*track2);
                  fGeoLengthPi = GeoLength(*trackN);
                  fGeoLengthPi1 = GeoLength(*track3);
                  //TOF Signal
                  Float_t mass = 0;
                  Float_t time = -1;
                  Float_t beta = 0;
                  Float_t gamma = 0;
                  Float_t length = 0;
                  Float_t time0 = 0;
                  //TofSignal He4
                  length = length = trackP->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(trackP->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = trackP->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (trackP->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalHe4 = mass;
                  }
                  //TOFSignal He3
                  length = length = track1->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track1->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track1->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalHe3 = mass;
                  }
                  //TOFSignal P
                  length = length = track2->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track2->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track2->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalP = mass;
                  }
                  //TOFSignal Pi
                  length = length = trackN->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(trackN->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = trackN->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (trackN->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalPi = mass;
                  }
                  //TOFSignal Pi1
                  length = length = track3->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track3->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track3->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalPi1 = mass;
                  }
                  //Track DeDx Sigma
                  if (fBetheSplines) {
                    fhe4DedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kAlpha);
                    fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
                    fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
                    fpiDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kPion);
                    fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                  } else {
                    fhe4DedxSigma = Bethe(*trackP, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
                    fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
                    fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
                    fpiDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kPion);
                    fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                  }
                  //DCA From primary Vertex
                  fhe4Dca = TMath::Abs(trackP->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpiDca = TMath::Abs(trackN->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpi1Dca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  //N Clusters TPC
                  fhe4Ncls = trackP->GetTPCNcls();
                  fhe3Ncls = track1->GetTPCNcls();
                  fpNcls = track2->GetTPCNcls();
                  fpiNcls = trackN->GetTPCNcls();
                  fpi1Ncls = track3->GetTPCNcls();
                  //NCls ITS
                  fhe4NclsITS = trackP->GetNumberOfITSClusters();
                  fhe3NclsITS = track1->GetNumberOfITSClusters();
                  fpNclsITS = track2->GetNumberOfITSClusters();
                  fpiNclsITS = trackN->GetNumberOfITSClusters();
                  fpi1NclsITS = track3->GetNumberOfITSClusters();
                  //V0 Vectors with 4LHe Mass
                  He4->SetXYZM(2*trackP->Px(), 2*trackP->Py(), 2*trackP->Pz(), 3.929);
                  pi->SetXYZM(trackN->Px(), trackN->Py(), trackN->Pz(), AliPID::ParticleMass(AliPID::kPion));
                  TLorentzVector He4Mother = *He4 + *pi;
                  //V0 M, p, pt, y
                  fmHHe4Pi = He4Mother.M();
                  fpHHe4Pi = He4Mother.P();
                  fptHHe4Pi = He4Mother.Pt();
                  fhe4P = He4->Pt();
                  fpiP = pi->Pt();
                  fyHHe4pi = He4Mother.Rapidity();
                  //V0 M, p, pt with Alpha Mass
                  He41->SetXYZM(2*trackP->Px(), 2*trackP->Py(), 2*trackP->Pz(), AliPID::ParticleMass(AliPID::kAlpha));
                  TLorentzVector L4HMother = *He41 + *pi;
                  fm4LH = L4HMother.M();
                  fp4LH = L4HMother.P();
                  fpt4LH = L4HMother.Pt();
                  //4LHe Daughter P and Vect
                  He3->SetXYZM(2*track1->Px(), 2*track1->Py(), 2*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
                  Prot->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
                  pi1->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
                  fhe3P = He3->Pt();
                  fpi1P = pi1->Pt();
                  fpP = Prot->Pt();
                  //Mother M (Daughter Sum)
                  TLorentzVector L4HeMother(0.,0.,0.,0.);
                  L4HeMother = *pi1 + *He3 + *Prot + *pi; 
                  TLorentzVector L4HeMother1(0.,0.,0.,0.);
                  //4LHe Mother
                  L4HeMother1 = *pi1 + *He3 + *Prot; 
                  fmDaughterSum1 = L4HeMother.M();
                  fmDaughterSum2 = L4HeMother1.M();
                  fpSum1 = L4HeMother.P();
                  fptSum1 = L4HeMother.Pt();
                  fpSum2 = L4HeMother1.P();
                  fptSum2 = L4HeMother1.Pt();
                  fySum1 = L4HeMother.Rapidity();
                  fySum2 = L4HeMother1.Rapidity();
                  //ct
                  TVector3 secondaryVertex = He4Vertex;
                  secondaryVertex = secondaryVertex - fPrimaryVertex;
                  TVector3 tertdiff = tertVtx - He4Vertex;
                  fVertDiff = TMath::Sqrt(secondaryVertex.X()*secondaryVertex.X() + secondaryVertex.Y()*secondaryVertex.Y());
                  fctHHe4Pi = secondaryVertex.Mag() * fmHHe4Pi / fpHHe4Pi;
                  fctSum1 = secondaryVertex.Mag() * fmDaughterSum1 / fpSum1;
                  fctSum2 = tertdiff.Mag() * fmDaughterSum2 / fpSum2;
                  fct4LH = secondaryVertex.Mag() * fm4LH / fp4LH;
                  //PA
                  Double_t PA3 = L4HeMother1.Angle(-(tertVtx - He4Vertex));
                  fPA4LHe = TMath::Cos(PA3);

                  Double_t PA2 = L4HeMother.Angle(-(He4Vertex - fPrimaryVertex));
                  fPA4LHe1 = TMath::Cos(PA2);
                  //DCA sec Vtx
                  fhe4DcaSec = TMath::Abs(trackP->GetD(He4Vertex.X(), He4Vertex.Y(), fMagneticField));
                  fhe3DcaTert = TMath::Abs(track1->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                  fpDcaTert = TMath::Abs(track2->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                  fpiDcaSec = TMath::Abs(trackN->GetD(He4Vertex.X(), He4Vertex.Y(), fMagneticField));
                  fpi1DcaTert = TMath::Abs(track3->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                  //Armenteros Podolanski
                  TVector3 vecN = pi1->Vect();
                  TVector3 vecP = He4->Vect();

                  TVector3 vecM = He4Mother.Vect();

                  fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                  fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                  farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                  farmpt = vecP.Mag()*sin(fthetaP);

                  vecN = pi1->Vect();
                  vecP = Prot->Vect();
                  vecM = vecN + vecP;

                  fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                  fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                  farmalpha1 = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                  farmpt1 = vecP.Mag()*sin(fthetaP);

                  fNV0Cand = fNV0Cand + 1;
                }//end PID
              }//end PID
            }//end PID
          }//end track3
          fNumberV0s = (fNV0Cand);
          if (fNV0Cand) fTree->Fill();
        }//end track2
      }//end track1
    }//end if helium4Positive && pionNegative
    
    //++++++++++++|| 4LLH --> 4LHe + pi (neg) ||++++++++++++//
    if (helium4Negative && pionPositive) {
      //++++++++++++++++++++++++++++++ He3 Track +++++++++++++++++++++++++++++++++++
      for (Int_t iTracks = 0; iTracks < fESDevent->GetNumberOfTracks(); iTracks++) {

        AliESDtrack* track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(iTracks));

        if (!track1->GetInnerParam()) continue;

        if (!trackCutsTert.AcceptTrack(track1)) continue;
        //ITS Cluster Cut defined in User Exec
        if(ITSClusterCut && track1->GetNumberOfITSClusters() > TertClusters) continue;  

        Double_t ptot1 = track1->GetInnerParam()->GetP();
        Double_t sign1 = track1->GetSign();

        if(sign1 > 0) continue;

        fHistdEdx->Fill(ptot1*sign1, track1->GetTPCsignal());

        if(fBetheSplines){
          if(TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) > 3) continue;
        }
        else {
          if(TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) > 3) continue;
        }

        //+++++++++++++++++++++++++++++++++++ p Track ++++++++++++++++++++++++++++++++++++++++++
        for (Int_t jTracks = iTracks + 1; jTracks < fESDevent->GetNumberOfTracks(); jTracks++) {

          AliESDtrack* track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(jTracks));              

          if (!track2->GetInnerParam()) continue;

          if (!trackCutsTert.AcceptTrack(track2)) continue;
          //ITS Cluster Cut defined in User Exec
          if(ITSClusterCut && track2->GetNumberOfITSClusters() > TertClusters) continue;

          Double_t ptot2 = track2->GetInnerParam()->GetP();
          Double_t sign2 = track2->GetSign();

          if(iTracks == jTracks) continue;

          if(sign2 > 0) continue; 

          fHistdEdx->Fill(ptot2*sign2, track2->GetTPCsignal());

          if(fBetheSplines){
            if(TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) > 3) continue;
          }
          else {
            if(TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) > 3) continue;
          }

          //+++++++++++++++++++++++++++++++++++ pi Track +++++++++++++++++++++++++++++++++++++++++
          for (Int_t kTracks = jTracks + 1; kTracks < fESDevent->GetNumberOfTracks(); kTracks++) {

            AliESDtrack* track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(kTracks));        

            if (!track3->GetInnerParam()) continue;

            if (!trackCutsPi.AcceptTrack(track3)) continue;
            //ITS Cluster Cut defined in User Exec
            if(ITSClusterCut && track3->GetNumberOfITSClusters() > TertClusters) continue;

            Double_t ptot3 = track3->GetInnerParam()->GetP();
            Double_t sign3 = track3->GetSign();

            if (sign3 < 0) continue;

            if (jTracks == kTracks || iTracks == kTracks) continue; //reject using the same track twice

            fHistdEdx->Fill(ptot3*sign3, track3->GetTPCsignal());

            if(TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) > 3) continue;
           
            if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 3)){
              if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3)) {
                if (TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) < 3) {

                  //====== Vertex Reconstruction ========
                  TObjArray *trkArray = new TObjArray(3);

                  trkArray->AddAt(track1,0);
                  trkArray->AddAt(track2,1);
                  trkArray->AddAt(track3,2);
                  //V0 Vtx --> sec Vtx
                  TVector3 He4Vertex(fV0->Xv(), fV0->Yv(), fV0->Zv());
                  AliESDVertex *He4Vert = new AliESDVertex(fV0->GetVertex());
                  vertexer1->SetVtxStart(He4Vert);
                  if(!trackP->PropagateToDCA(He4Vert, fMagneticField, 10)) continue;//10
                  if(!trackN->PropagateToDCA(He4Vert, fMagneticField, 10)) continue;//10
                  if(He4Vert) delete He4Vert;
                  
                  AliESDVertex *tertVertex1 = (AliESDVertex*)vertexer1->VertexForSelectedESDTracks(trkArray);
                  if(trkArray) delete trkArray;
                  //HyperHelium4 Decay
                  if(!track1->PropagateToDCA(tertVertex1, fMagneticField, 10)) continue;//10
                  if(!track2->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;//10
                  if(!track3->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;//10
                  //DCA
                  Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
                  fDcaHe3P = *dca;

                  Double_t *dca1 = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
                  fDcaHe3Pi2 = *dca1;

                  Double_t *dca3 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
                  fDcaPPi1 = *dca3;

                  if(dca) delete dca;
                  if(dca1) delete dca1;
                  if(dca3) delete dca3;

                  TVector3 tertVtx(tertVertex1->GetX(),tertVertex1->GetY(),tertVertex1->GetZ());
                  if(tertVertex1) delete tertVertex1;
                  //==================================================================================================//
                  //BR for separation of different analysis: 1 = 4LLH --> (4LHe, pi), 2 = 4LLH -->(3LH, p, pi), 3 = Li4
                  fBR = 1;
                  fmc = -1; 
                  fz = -1;
                  fParticleSpecies = -1020010040;
                  //V0 DCA
                  fdcaHe4Pi = fV0->GetDcaV0Daughters();
                  //V0 Pointing-Angle
                  fcosHe4Pi = fV0->GetV0CosineOfPointingAngle();
                  //V0-OnFlyStatus
                  if(fV0->GetOnFlyStatus()) fonTheFly = 1;
                  else if(!fV0->GetOnFlyStatus()) fonTheFly = 0;
                  //Track TPC Signal
                  fhe4Dedx = trackN->GetTPCsignal();
                  fhe3Dedx = track1->GetTPCsignal();
                  fpDedx = track2->GetTPCsignal();
                  fpiDedx = trackP->GetTPCsignal();
                  fpi1Dedx = track3->GetTPCsignal();
                  //Eta
                  fEtaHe4 = trackN->Eta();
                  fEtaHe3 = track1->Eta();
                  fEtaP = track2->Eta();
                  fEtaPi = trackP->Eta();
                  fEtaPi1 = track3->Eta();
                  //Phi
                  fPhiHe4 = trackN->Phi();
                  fPhiHe3 = track1->Phi();
                  fPhiP = track2->Phi();
                  fPhiPi = trackP->Phi();
                  fPhiPi1 = track3->Phi();
                  //GeoLength
                  fGeoLengthHe4 = GeoLength(*trackN);
                  fGeoLengthHe3 = GeoLength(*track1);
                  fGeoLengthP = GeoLength(*track2);
                  fGeoLengthPi = GeoLength(*trackP);
                  fGeoLengthPi1 = GeoLength(*track3);
                  //TOF Signal
                  Float_t mass = 0;
                  Float_t time = -1;
                  Float_t beta = 0;
                  Float_t gamma = 0;
                  Float_t length = 0;
                  Float_t time0 = 0;
                  //TofSignal He4
                  length = length = trackN->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(trackN->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = trackN->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (trackN->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalHe4 = mass;
                  }
                  //TOFSignal He3
                  length = length = track1->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track1->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track1->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalHe3 = mass;
                  }
                  //TOFSignal P
                  length = length = track2->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track2->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track2->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalP = mass;
                  }
                  //TOFSignal Pi
                  length = length = trackP->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(trackP->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = trackP->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (trackP->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalPi = mass;
                  }
                  //TOFSignal Pi1
                  length = length = track3->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track3->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track3->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalPi1 = mass;
                  }
                  //Track DeDx Sigma
                  if (fBetheSplines) {
                    fhe4DedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kAlpha);
                    fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
                    fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
                    fpiDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kPion);
                    fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                  } else {
                    fhe4DedxSigma = Bethe(*trackN, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
                    fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
                    fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
                    fpiDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kPion);
                    fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                  }
                  //DCA From primary Vertex
                  fhe4Dca = TMath::Abs(trackN->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpiDca = TMath::Abs(trackP->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpi1Dca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  //N Clusters TPC
                  fhe4Ncls = trackN->GetTPCNcls();
                  fhe3Ncls = track1->GetTPCNcls();
                  fpNcls = track2->GetTPCNcls();
                  fpiNcls = trackP->GetTPCNcls();
                  fpi1Ncls = track3->GetTPCNcls();
                  //NCls ITS
                  fhe4NclsITS = trackN->GetNumberOfITSClusters();
                  fhe3NclsITS = track1->GetNumberOfITSClusters();
                  fpNclsITS = track2->GetNumberOfITSClusters();
                  fpiNclsITS = trackP->GetNumberOfITSClusters();
                  fpi1NclsITS = track3->GetNumberOfITSClusters();
                  //V0 Vect with 4LHe Mass
                  He4->SetXYZM(2*trackN->Px(), 2*trackN->Py(), 2*trackN->Pz(), 3.929);
                  pi->SetXYZM(trackP->Px(), trackP->Py(), trackP->Pz(), AliPID::ParticleMass(AliPID::kPion));
                  TLorentzVector He4Mother = *He4 + *pi;
                  //V0 M, p, pt, y
                  fmHHe4Pi = He4Mother.M();
                  fpHHe4Pi = He4Mother.P();
                  fptHHe4Pi = He4Mother.Pt();
                  fhe4P = He4->Pt();
                  fpiP = pi->Pt();
                  fyHHe4pi = He4Mother.Rapidity();
                  //V0 with Alpha Mass
                  He41->SetXYZM(2*trackN->Px(), 2*trackN->Py(), 2*trackN->Pz(), AliPID::ParticleMass(AliPID::kAlpha));
                  TLorentzVector L4HMother = *He41 + *pi;
                  fm4LH = L4HMother.M();
                  fp4LH = L4HMother.P();
                  fpt4LH = L4HMother.Pt();
                  //4LHe Daughters
                  He3->SetXYZM(2*track1->Px(), 2*track1->Py(), 2*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
                  pi1->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
                  Prot->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
                  fhe3P = He3->Pt();
                  fpi1P = pi1->Pt();
                  fpP = Prot->Pt();
                  //Mother Sum Daughters
                  TLorentzVector L4HeMother(0.,0.,0.,0.);
                  L4HeMother = *pi1 + *He3 + *Prot + *pi; 
                  //4LHe Sum Daughters
                  TLorentzVector L4HeMother1(0.,0.,0.,0.);
                  L4HeMother1 = *pi1 + *He3 + *Prot; 
                  fmDaughterSum1 = L4HeMother.M();
                  fmDaughterSum2 = L4HeMother1.M();
                  fpSum1 = L4HeMother.P();
                  fptSum1 = L4HeMother.Pt();
                  fpSum2 = L4HeMother1.P();
                  fptSum2 = L4HeMother1.Pt();
                  fySum1 = L4HeMother.Rapidity();
                  fySum2 = L4HeMother1.Rapidity();
                  //ct
                  TVector3 secondaryVertex = He4Vertex;
                  secondaryVertex = secondaryVertex - fPrimaryVertex;
                  TVector3 tertdiff = tertVtx - He4Vertex;
                  fVertDiff = TMath::Sqrt(secondaryVertex.X()*secondaryVertex.X() + secondaryVertex.Y()*secondaryVertex.Y());
                  fctHHe4Pi = secondaryVertex.Mag() * fmHHe4Pi / fpHHe4Pi;
                  fctSum1 = secondaryVertex.Mag() * fmDaughterSum1 / fpSum1;
                  fctSum2 = tertdiff.Mag() * fmDaughterSum2 / fpSum2;
                  fct4LH = secondaryVertex.Mag() * fm4LH / fp4LH;
                  //PA
                  Double_t PA3 = L4HeMother1.Angle(-(tertVtx - He4Vertex));
                  fPA4LHe = TMath::Cos(PA3);

                  Double_t PA2 = L4HeMother.Angle(-(He4Vertex - fPrimaryVertex));
                  fPA4LHe1 = TMath::Cos(PA2);
                  //DCA sec / tert
                  fhe4DcaSec = TMath::Abs(trackN->GetD(He4Vertex.X(), He4Vertex.Y(), fMagneticField));
                  fhe3DcaTert = TMath::Abs(track1->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                  fpDcaTert = TMath::Abs(track2->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                  fpiDcaSec = TMath::Abs(trackP->GetD(He4Vertex.X(), He4Vertex.Y(), fMagneticField));
                  fpi1DcaTert = TMath::Abs(track3->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
                  //Armenteros Podolanski
                  TVector3 vecP = pi->Vect();
                  TVector3 vecN = He4->Vect();

                  TVector3 vecM = He4Mother.Vect();

                  fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                  fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                  farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                  farmpt = vecP.Mag()*sin(fthetaP);

                  vecP = pi1->Vect();
                  vecN = Prot->Vect();
                  vecM = vecN + vecP;

                  fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                  fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                  farmalpha1 = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                  farmpt1 = vecP.Mag()*sin(fthetaP);

                  fNV0Cand = fNV0Cand + 1;  
                }//end PID
              }//end PID
            }//end PID
          }//end track3
          fNumberV0s = (fNV0Cand);
          if (fNV0Cand) fTree->Fill();
        }//end track2
      }//end track1
    }//end helium4Negativ && pionPositive 
    
    //++++++++++++++++++++++++++ 4LLH-->(3LH, p, pi) pos ++++++++++++++++++++++++++++++//      
    if (helium3Positive && pionNegative) {
      //++++++++++++++++++++++++++++++ He3 Track +++++++++++++++++++++++++++++++++++
      for (Int_t iTracks = 0; iTracks < fESDevent->GetNumberOfTracks(); iTracks++) {

        AliESDtrack* track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(iTracks));

        if (!track1->GetInnerParam()) continue;

        if (!trackCutsTert.AcceptTrack(track1)) continue;
        //ITS Cluster Cut defined in User Exec
        if(ITSClusterCut && track1->GetNumberOfITSClusters() >= TertClusters) continue;  

        Double_t ptot1 = track1->GetInnerParam()->GetP();
        Double_t sign1 = track1->GetSign();

        if(sign1 < 0) continue;

        fHistdEdx->Fill(ptot1*sign1, track1->GetTPCsignal());

        if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) > 3) continue;
        }
        else {
          if(TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) > 3) continue;
        }

        //!! special conditions for He3 Track --> track must be able to be constrained to V0 Vertex !!
        AliESDVertex *HyperTritonV = new AliESDVertex(fV0->GetVertex());
        Double_t b[3];
        track1->GetBxByBz(b);
        if(!track1->ConstrainToVertex(HyperTritonV, b)) continue;
        if(HyperTritonV) delete HyperTritonV;


        //+++++++++++++++++++++++++++++++++++ p Track ++++++++++++++++++++++++++++++++++++++++++
        for (Int_t jTracks = iTracks + 1; jTracks < fESDevent->GetNumberOfTracks(); jTracks++) {

          AliESDtrack* track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(jTracks));              

          if (!track2->GetInnerParam()) continue;

          if (!trackCutsTert.AcceptTrack(track2)) continue;
          //ITS Cluster Cut defined in User Exec
          if(ITSClusterCut && track2->GetNumberOfITSClusters() >= TertClusters) continue;

          Double_t ptot2 = track2->GetInnerParam()->GetP();
          Double_t sign2 = track2->GetSign();

          if(iTracks == jTracks) continue;

          if(sign2 < 0) continue; 

          fHistdEdx->Fill(ptot2*sign2, track2->GetTPCsignal());

        if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) > 3) continue;
        }
        else {
          if(TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) > 3) continue;
        }
        //+++++++++++++++++++++++++++++++++++ pi Track +++++++++++++++++++++++++++++++++++++++++
          for (Int_t kTracks = jTracks + 1; kTracks < fESDevent->GetNumberOfTracks(); kTracks++) {

            AliESDtrack* track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(kTracks));        

            if (!track3->GetInnerParam()) continue;

            if (!trackCutsPi.AcceptTrack(track3)) continue;
            //ITS Cluster Cut defined in User Exec
            if(ITSClusterCut && track3->GetNumberOfITSClusters() > TertClusters) continue;

            Double_t ptot3 = track3->GetInnerParam()->GetP();
            Double_t sign3 = track3->GetSign();

            if (sign3 > 0) continue;

            if (jTracks == kTracks || iTracks == kTracks) continue; //reject using the same track twice

            fHistdEdx->Fill(ptot3*sign3, track3->GetTPCsignal());

            if(TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) > 3) continue;
            
            //-------------------------------------------------------------------------------------------------//
            if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 3)){
              if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3)) {
                if (TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) < 3) {
                  //======= Vertex Reconstruction ======
                  TObjArray *trkArray = new TObjArray(3);

                  trkArray->AddAt(track1,0);
                  trkArray->AddAt(track2,1);
                  trkArray->AddAt(track3,2);

                  //HyperTriton V0 Vertex
                  TVector3 HypTritVertex(fV0->Xv(), fV0->Yv(), fV0->Zv());
                  AliESDVertex *HypTritVert = new AliESDVertex(fV0->GetVertex());
                  if(!trackP->PropagateToDCA(HypTritVert, fMagneticField, 10)) continue;
                  if(!trackN->PropagateToDCA(HypTritVert, fMagneticField, 10)) continue;
                  if(HypTritVert) delete HypTritVert;

                  AliESDVertex *SecVertex1 = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
                  if(trkArray) delete trkArray;
                  //HyperHelium4 Decay
                  if(!track1->PropagateToDCA(SecVertex1, fMagneticField, 10)) continue;
                  if(!track2->PropagateToDCA(SecVertex1, fMagneticField, 10))continue;
                  if(!track3->PropagateToDCA(SecVertex1, fMagneticField, 10))continue;
                  //DCA
                  Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
                  fDcaHe3P = *dca;

                  Double_t *dca1 = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
                  fDcaHe3Pi2 = *dca1;

                  Double_t *dca3 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
                  fDcaPPi1 = *dca3;

                  if(dca) delete dca;
                  if(dca1) delete dca1;
                  if(dca3) delete dca3;

                  TVector3 SecVtx(SecVertex1->GetX(),SecVertex1->GetY(),SecVertex1->GetZ());
                  if(SecVertex1) delete SecVertex1;

                  //========================================================================//
                  //BR for separation of different analysis: 1 = 4LLH --> (4LHe, pi), 2 = 4LLH -->(3LH, p, pi), 3 = Li4
                  fBR = 2;
                  fmc = -1; 
                  fz = 1;
                  fParticleSpecies = 1020010040;
                  //V0 DCA
                  fdcaHe4Pi = fV0->GetDcaV0Daughters();
                  //V0 Pointing-Angle
                  fcosHe4Pi = fV0->GetV0CosineOfPointingAngle();
                  //V0-OnFlyStatus
                  if(fV0->GetOnFlyStatus()) fonTheFly = 1;
                  else if(!fV0->GetOnFlyStatus()) fonTheFly = 0;
                  //Track TPC Signal
                  fhe4Dedx = track1->GetTPCsignal();
                  fhe3Dedx = trackP->GetTPCsignal();
                  fpDedx = track2->GetTPCsignal();
                  fpiDedx = trackN->GetTPCsignal();
                  fpi1Dedx = track3->GetTPCsignal();
                  //Eta
                  fEtaHe4 = track1->Eta();
                  fEtaHe3 = trackP->Eta();
                  fEtaP = track2->Eta();
                  fEtaPi = trackN->Eta();
                  fEtaPi1 = track3->Eta();
                  //Phi
                  fPhiHe4 = track1->Phi();
                  fPhiHe3 = trackP->Phi();
                  fPhiP = track2->Phi();
                  fPhiPi = trackN->Phi();
                  fPhiPi1 = track3->Phi();
                  //GeoLength
                  fGeoLengthHe4 = GeoLength(*track1);
                  fGeoLengthHe3 = GeoLength(*trackP);
                  fGeoLengthP = GeoLength(*track2);
                  fGeoLengthPi = GeoLength(*trackN);
                  fGeoLengthPi1 = GeoLength(*track3);
                  //TOF Signal
                  Float_t mass = 0;
                  Float_t time = -1;
                  Float_t beta = 0;
                  Float_t gamma = 0;
                  Float_t length = 0;
                  Float_t time0 = 0;
                  //TofSignal He3
                  length = length = trackP->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(trackP->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = trackP->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (trackP->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalHe3 = mass;
                  }
                  //TOFSignal Hypertriton
                  length = length = track1->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track1->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track1->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalHe4 = mass;
                  }
                  //TOFSignal P
                  length = length = track2->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track2->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track2->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalP = mass;
                  }
                  //TOFSignal Pi
                  length = length = trackN->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(trackN->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = trackN->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (trackN->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalPi = mass;
                  }
                  //TOFSignal Pi1
                  length = length = track3->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track3->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track3->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalPi1 = mass;
                  }
                  //Track DeDx Sigma
                  if (fBetheSplines) {
                    fhe4DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
                    fhe3DedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kHe3);
                    fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
                    fpiDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kPion);
                    fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                  } else {
                    fhe4DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
                    fhe3DedxSigma = Bethe(*trackP, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
                    fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
                    fpiDedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kPion);
                    fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                  }
                  //DCA From primary Vertex
                  fhe4Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fhe3Dca = TMath::Abs(trackP->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpiDca = TMath::Abs(trackN->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpi1Dca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  //N Clusters TPC
                  fhe4Ncls = track1->GetTPCNcls();
                  fhe3Ncls = trackP->GetTPCNcls();
                  fpNcls = track2->GetTPCNcls();
                  fpiNcls = trackN->GetTPCNcls();
                  fpi1Ncls = track3->GetTPCNcls();
                  //ITS NCls
                  fhe4NclsITS = track1->GetNumberOfITSClusters();
                  fhe3NclsITS = trackP->GetNumberOfITSClusters();
                  fpNclsITS = track2->GetNumberOfITSClusters();
                  fpiNclsITS = trackN->GetNumberOfITSClusters();
                  fpi1NclsITS = track3->GetNumberOfITSClusters();
                  //He3 Track as Hypertriton!! + Proton + pi = 4LLH
		              HypTrit->SetXYZM(2*track1->Px(), 2*track1->Py(), 2*track1->Pz(), 2.99131);
		              Prot->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
		              pi1->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
                  fhe4P = HypTrit->Pt();
                  fpi1P = pi1->Pt();
                  fpP = Prot->Pt();
                  //V0 Daughters
                  He3->SetXYZM(2*trackP->Px(), 2*trackP->Py(), 2*trackP->Pz(), AliPID::ParticleMass(AliPID::kHe3));
                  pi->SetXYZM(trackN->Px(), trackN->Py(), trackN->Pz(), AliPID::ParticleMass(AliPID::kPion));
                  fhe3P = He3->Pt();
                  fpiP = pi->Pt();
                  //Mother (tracks!)
                  TLorentzVector HypTritMother = *HypTrit + *pi1 + *Prot;
                  fmHHe4Pi = HypTritMother.M();
                  fpHHe4Pi = HypTritMother.P();
                  fptHHe4Pi = HypTritMother.Pt();
                  fyHHe4pi = HypTritMother.Rapidity();
                  //Mother (V0 tracks + p, pi tracks)
                  TLorentzVector L4HeMother(0.,0.,0.,0.);
                  L4HeMother = *pi1 + *He3 + *Prot + *pi; 
                  fmDaughterSum1 = L4HeMother.M();
                  fpSum1 = L4HeMother.P();
                  fptSum1 = L4HeMother.Pt();
                  fySum1 = L4HeMother.Rapidity();
                  //Hypertriton Mother
                  TLorentzVector L4HeMother1(0.,0.,0.,0.);
                  L4HeMother1 = *pi + *He3; 
                  fmDaughterSum2 = L4HeMother1.M();
                  fpSum2 = L4HeMother1.P();
                  fptSum2 = L4HeMother1.Pt();
                  fySum2 = L4HeMother1.Rapidity();
                  //ct
                  TVector3 secondaryVertex = SecVtx;
                  secondaryVertex = secondaryVertex - fPrimaryVertex;
                  TVector3 tertdiff = SecVtx - HypTritVertex;
                  fVertDiff = TMath::Sqrt(secondaryVertex.X()*secondaryVertex.X() + secondaryVertex.Y()*secondaryVertex.Y());
                  fctHHe4Pi = secondaryVertex.Mag() * fmHHe4Pi / fpHHe4Pi;
                  fctSum1 = secondaryVertex.Mag() * fmDaughterSum1 / fpSum1;
                  //HyperTriton ct
                  fctSum2 = tertdiff.Mag() * fmDaughterSum2 / fpSum2;
                  //PA
                  Double_t PA3 = L4HeMother1.Angle(-(HypTritVertex - SecVtx));
                  fPA4LHe1 = TMath::Cos(PA3);

                  Double_t PA2 = L4HeMother.Angle(-(SecVtx - fPrimaryVertex));
                  fPA4LHe = TMath::Cos(PA2);
                  //DCA sec vtx
                  fhe4DcaSec = TMath::Abs(track1->GetD(SecVtx.X(), SecVtx.Y(), fMagneticField));
                  fhe3DcaTert = TMath::Abs(trackP->GetD(HypTritVertex.X(), HypTritVertex.Y(), fMagneticField));
                  fpDcaTert = TMath::Abs(track2->GetD(SecVtx.X(), SecVtx.Y(), fMagneticField));
                  fpiDcaSec = TMath::Abs(trackN->GetD(HypTritVertex.X(), HypTritVertex.Y(), fMagneticField));
                  fpi1DcaTert = TMath::Abs(track3->GetD(SecVtx.X(), SecVtx.Y(), fMagneticField));
                  //Armenteros Podolanski
                  TVector3 vecN = pi->Vect();
                  TVector3 vecP = He3->Vect();

                  TVector3 vecM = vecP + vecN;

                  fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                  fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                  farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                  farmpt = vecP.Mag()*sin(fthetaP);

                  vecN = pi1->Vect();
                  vecP = Prot->Vect();
                  vecM = vecN + vecP;

                  fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                  fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                  farmalpha1 = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                  farmpt1 = vecP.Mag()*sin(fthetaP);

                  fNV0Cand = fNV0Cand + 1;
                  
                }//end PID
              }//end PID
            }//end PID
          }//end track3
          fNumberV0s = (fNV0Cand);
          if (fNV0Cand) fTree->Fill();
        }//end track2
      }//end track1
    }//end if helium3Positive && pionNegative
    if (helium3Negative && pionPositive) {
    //++++++++++++++++++++++++++++++ He3 Track +++++++++++++++++++++++++++++++++++
      for (Int_t iTracks = 0; iTracks < fESDevent->GetNumberOfTracks(); iTracks++) {

        AliESDtrack* track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(iTracks));

        if (!track1->GetInnerParam()) continue;

        if (!trackCutsTert.AcceptTrack(track1)) continue;

        if(ITSClusterCut && track1->GetNumberOfITSClusters() >= TertClusters) continue;  

        Double_t ptot1 = track1->GetInnerParam()->GetP();
        Double_t sign1 = track1->GetSign();

        if(sign1 > 0) continue;

        fHistdEdx->Fill(ptot1*sign1, track1->GetTPCsignal());

        if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) > 3) continue;
        }
        else {
          if(TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) > 3) continue;
        }

        AliESDVertex *HyperTritonV = new AliESDVertex(fV0->GetVertex());
        Double_t b[3];
        track1->GetBxByBz(b);
        if(!track1->ConstrainToVertex(HyperTritonV, b)) continue;
        if(HyperTritonV) delete HyperTritonV;


        //+++++++++++++++++++++++++++++++++++ p Track ++++++++++++++++++++++++++++++++++++++++++
        for (Int_t jTracks = iTracks + 1; jTracks < fESDevent->GetNumberOfTracks(); jTracks++) {

          AliESDtrack* track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(jTracks));              

          if (!track2->GetInnerParam()) continue;

          if (!trackCutsTert.AcceptTrack(track2)) continue;

          if(ITSClusterCut && track2->GetNumberOfITSClusters() >= TertClusters) continue;

          Double_t ptot2 = track2->GetInnerParam()->GetP();
          Double_t sign2 = track2->GetSign();

          if(iTracks == jTracks) continue;

          if(sign2 > 0) continue; 

          fHistdEdx->Fill(ptot2*sign2, track2->GetTPCsignal());

        if(fBetheSplines){
        if(TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) > 3) continue;
        }
        else {
          if(TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) > 3) continue;
        }
        //+++++++++++++++++++++++++++++++++++ pi Track +++++++++++++++++++++++++++++++++++++++++
          for (Int_t kTracks = jTracks + 1; kTracks < fESDevent->GetNumberOfTracks(); kTracks++) {

            AliESDtrack* track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(kTracks));        

            if (!track3->GetInnerParam()) continue;

            if (!trackCutsPi.AcceptTrack(track3)) continue;

            if(ITSClusterCut && track3->GetNumberOfITSClusters() > TertClusters) continue;

            Double_t ptot3 = track3->GetInnerParam()->GetP();
            Double_t sign3 = track3->GetSign();

            if (sign3 < 0) continue;

            if (jTracks == kTracks || iTracks == kTracks) continue; //reject using the same track twice

            fHistdEdx->Fill(ptot3*sign3, track3->GetTPCsignal());

            if(TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) > 3) continue;
            
            if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 3)){
              if ((fBetheSplines && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) < 3) || (!fBetheSplines && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) < 3)) {
                if (TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) < 3) {
                  //======= Vertex Reconstruction ======
                  TObjArray *trkArray = new TObjArray(3);

                  trkArray->AddAt(track1,0);
                  trkArray->AddAt(track2,1);
                  trkArray->AddAt(track3,2);

                  TVector3 HypTritVertex(fV0->Xv(), fV0->Yv(), fV0->Zv());
                  AliESDVertex *HypTritVert = new AliESDVertex(fV0->GetVertex());
                  
                  if(!trackP->PropagateToDCA(HypTritVert, fMagneticField, 10)) continue;
                  if(!trackN->PropagateToDCA(HypTritVert, fMagneticField, 10)) continue;
                  if(HypTritVert) delete HypTritVert;

                  AliESDVertex *SecVertex1 = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
                  if(trkArray) delete trkArray;
                  //HyperHelium4 Decay
                  if(!track1->PropagateToDCA(SecVertex1, fMagneticField, 10)) continue;
                  if(!track2->PropagateToDCA(SecVertex1, fMagneticField, 10))continue;
                  if(!track3->PropagateToDCA(SecVertex1, fMagneticField, 10))continue;

                  Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
                  fDcaHe3P = *dca;

                  Double_t *dca1 = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
                  fDcaHe3Pi2 = *dca1;

                  Double_t *dca3 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
                  fDcaPPi1 = *dca3;

                  if(dca) delete dca;
                  if(dca1) delete dca1;
                  if(dca3) delete dca3;

                  TVector3 SecVtx(SecVertex1->GetX(),SecVertex1->GetY(),SecVertex1->GetZ());

                  if(SecVertex1) delete SecVertex1;

                  //==================================================
                  fBR = 2;
                  fmc = -1; 
                  fz = -1;
                  fParticleSpecies = -1020010040;
                  //V0 DCA
                  fdcaHe4Pi = fV0->GetDcaV0Daughters();
                  //V0 Pointing-Angle
                  fcosHe4Pi = fV0->GetV0CosineOfPointingAngle();
                  //V0-OnFlyStatus
                  if(fV0->GetOnFlyStatus()) fonTheFly = 1;
                  else if(!fV0->GetOnFlyStatus()) fonTheFly = 0;
                  //Track TPC Signal
                  fhe4Dedx = track1->GetTPCsignal();
                  fhe3Dedx = trackN->GetTPCsignal();
                  fpDedx = track2->GetTPCsignal();
                  fpiDedx = trackP->GetTPCsignal();
                  fpi1Dedx = track3->GetTPCsignal();
                  //Eta
                  fEtaHe4 = track1->Eta();
                  fEtaHe3 = trackN->Eta();
                  fEtaP = track2->Eta();
                  fEtaPi = trackP->Eta();
                  fEtaPi1 = track3->Eta();
                  //Phi
                  fPhiHe4 = track1->Phi();
                  fPhiHe3 = trackN->Phi();
                  fPhiP = track2->Phi();
                  fPhiPi = trackP->Phi();
                  fPhiPi1 = track3->Phi();
                  //GeoLength
                  fGeoLengthHe4 = GeoLength(*track1);
                  fGeoLengthHe3 = GeoLength(*trackN);
                  fGeoLengthP = GeoLength(*track2);
                  fGeoLengthPi = GeoLength(*trackP);
                  fGeoLengthPi1 = GeoLength(*track3);
                  //TOF Signal
                  Float_t mass = 0;
                  Float_t time = -1;
                  Float_t beta = 0;
                  Float_t gamma = 0;
                  Float_t length = 0;
                  Float_t time0 = 0;
                  //TofSignal He3
                  length = length = trackN->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(trackN->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = trackN->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (trackN->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalHe3 = mass;
                  }
                  //TOFSignal Hypertriton
                  length = length = track1->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track1->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track1->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalHe4 = mass;
                  }
                  //TOFSignal P
                  length = length = track2->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track2->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track2->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalP = mass;
                  }
                  //TOFSignal Pi
                  length = length = trackP->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(trackP->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = trackP->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (trackP->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalPi = mass;
                  }
                  //TOFSignal Pi1
                  length = length = track3->GetIntegratedLength();
                  time0 = fPID->GetTOFResponse().GetStartTime(track3->P());//fESDpid->GetTOFResponse().GetTimeZero();
                  time = track3->GetTOFsignal() - time0;
                  if (time > 0) {
                    beta = length / (2.99792457999999984e-02 * time);
                    gamma = 1/TMath::Sqrt(1 - beta*beta);
                    mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
                    fTOFSignalPi1 = mass;
                  }
                  //Track DeDx Sigma
                  if (fBetheSplines) {
                    fhe4DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
                    fhe3DedxSigma = fPID->NumberOfSigmasTPC(trackN, AliPID::kHe3);
                    fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
                    fpiDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kPion);
                    fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                  } else {
                    fhe4DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
                    fhe3DedxSigma = Bethe(*trackN, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
                    fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
                    fpiDedxSigma = fPID->NumberOfSigmasTPC(trackP, AliPID::kPion);
                    fpi1DedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
                  }
                  //DCA From primary Vertex
                  fhe4Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fhe3Dca = TMath::Abs(trackN->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpiDca = TMath::Abs(trackP->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  fpi1Dca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
                  //N Clusters TPC
                  fhe4Ncls = track1->GetTPCNcls();
                  fhe3Ncls = trackN->GetTPCNcls();
                  fpNcls = track2->GetTPCNcls();
                  fpiNcls = trackP->GetTPCNcls();
                  fpi1Ncls = track3->GetTPCNcls();
                  //ITS N Cls
                  fhe4NclsITS = track1->GetNumberOfITSClusters();
                  fhe3NclsITS = trackN->GetNumberOfITSClusters();
                  fpNclsITS = track2->GetNumberOfITSClusters();
                  fpiNclsITS = trackP->GetNumberOfITSClusters();
                  fpi1NclsITS = track3->GetNumberOfITSClusters();
                  //Vectors
		              HypTrit->SetXYZM(2*track1->Px(), 2*track1->Py(), 2*track1->Pz(), 2.99131);
		              Prot->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
		              pi1->SetXYZM(track3->Px(), track3->Py(), track3->Pz(), AliPID::ParticleMass(AliPID::kPion));
                  fhe4P = HypTrit->Pt();
                  fpi1P = pi1->Pt();
                  fpP = Prot->Pt();
                  //V0 Vec
                  pi->SetXYZM(trackP->Px(), trackP->Py(), trackP->Pz(), AliPID::ParticleMass(AliPID::kPion));
                  He3->SetXYZM(2*trackN->Px(), 2*trackN->Py(), 2*trackN->Pz(), AliPID::ParticleMass(AliPID::kHe3));
                  fhe3P = He3->Pt();
                  fpiP = pi->Pt();
                  //Mother Sum tracks
                  TLorentzVector HypTritMother = *HypTrit + *pi1 + *Prot;
                  fmHHe4Pi = HypTritMother.M();
                  fpHHe4Pi = HypTritMother.P();
                  fptHHe4Pi = HypTritMother.Pt();
                  fyHHe4pi = HypTritMother.Rapidity();
                  //Mother sum daughters
                  TLorentzVector L4HeMother(0.,0.,0.,0.);
                  L4HeMother = *pi1 + *He3 + *Prot + *pi; 
                  //Hypertriton
                  TLorentzVector L4HeMother1(0.,0.,0.,0.);
                  L4HeMother1 = *pi + *He3; 
                  fmDaughterSum1 = L4HeMother.M();
                  fmDaughterSum2 = L4HeMother1.M();
                  fpSum1 = L4HeMother.P();
                  fptSum1 = L4HeMother.Pt();
                  fpSum2 = L4HeMother1.P();
                  fptSum2 = L4HeMother1.Pt();
                  fySum1 = L4HeMother.Rapidity();
                  fySum2 = L4HeMother1.Rapidity();
                  //ct
                  TVector3 secondaryVertex = SecVtx;
                  secondaryVertex = secondaryVertex - fPrimaryVertex;
                  TVector3 tertdiff = SecVtx - HypTritVertex;
                  fVertDiff = TMath::Sqrt(secondaryVertex.X()*secondaryVertex.X() + secondaryVertex.Y()*secondaryVertex.Y());
                  fctHHe4Pi = secondaryVertex.Mag() * fmHHe4Pi / fpHHe4Pi;
                  fctSum1 = secondaryVertex.Mag() * fmDaughterSum1 / fpSum1;
                  fctSum2 = tertdiff.Mag() * fmDaughterSum2 / fpSum2;
                  //PA
                  Double_t PA3 = L4HeMother1.Angle(-(HypTritVertex - SecVtx));
                  fPA4LHe1 = TMath::Cos(PA3);

                  Double_t PA2 = L4HeMother.Angle(-(SecVtx - fPrimaryVertex));
                  fPA4LHe = TMath::Cos(PA2);
                  //DCA sec / tert Vtx
                  fhe4DcaSec = TMath::Abs(track1->GetD(SecVtx.X(), SecVtx.Y(), fMagneticField));
                  fhe3DcaTert = TMath::Abs(trackN->GetD(HypTritVertex.X(), HypTritVertex.Y(), fMagneticField));
                  fpDcaTert = TMath::Abs(track2->GetD(SecVtx.X(), SecVtx.Y(), fMagneticField));
                  fpiDcaSec = TMath::Abs(trackP->GetD(HypTritVertex.X(), HypTritVertex.Y(), fMagneticField));
                  fpi1DcaTert = TMath::Abs(track3->GetD(SecVtx.X(), SecVtx.Y(), fMagneticField));
                  //Armenteros Podolanski
                  TVector3 vecP = pi->Vect();
                  TVector3 vecN = He3->Vect();

                  TVector3 vecM = vecP + vecN;

                  fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                  fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                  farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                  farmpt = vecP.Mag()*sin(fthetaP);

                  vecP = pi1->Vect();
                  vecN = Prot->Vect();
                  vecM = vecN + vecP;

                  fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
                  fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

                  farmalpha1 = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
                  farmpt1 = vecP.Mag()*sin(fthetaP);

                  fNV0Cand = fNV0Cand + 1;
                  
                }//end PID
              }//end PID
            }//end PID
          }//end track3
          fNumberV0s = (fNV0Cand);
          if (fNV0Cand) fTree->Fill();
        }//end track2
      }//end track1
    }//end if helium3Negative && pionPositive
  }//end V0 Loop
  
}

void AliAnalysisTaskDoubleHypNucTree::Terminate(const Option_t*) {
  if (!GetOutputData(0)) return;
}

/// Set trigger information in reduced event
/// \return returns kTRUE is successful.
Bool_t AliAnalysisTaskDoubleHypNucTree::TriggerSelection() {
  fTrigger = 0;
  TString classes = fESDevent->GetFiredTriggerClasses();
  fTriggerClasses = classes;
  if (classes.Contains("CINT7")) fTrigger = 1;
  if (classes.Contains("CVHMV0M")) fTrigger = 2;
  if (classes.Contains("CVHMSH2") || classes.Contains("CSHM")) fTrigger = 3;
  if (classes.Contains("HNU")) fTrigger = 4;
  if (classes.Contains("HQU")) fTrigger = 5;
  if (classes.Contains("HJT")) fTrigger = 6;
  fHistTrigger->Fill(fTrigger);
  Bool_t isTriggered = kTRUE;
  if (fTrigger == 0) isTriggered = kFALSE;
  return isTriggered;
}
/// Calculates number of sigma deviation from expected dE/dx in TPC
/// \param track particle track
/// \param mass mass hypothesis of particle
/// \param charge particle charge hypothesis
/// \param params Parameters of Aleph parametrization of Bethe Energy-loss
Double_t AliAnalysisTaskDoubleHypNucTree::Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params){
  Double_t expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass,params[0],params[1],params[2],params[3],params[4]);
  Double_t sigma = expected*params[5];
  if (TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;
}

Double_t AliAnalysisTaskDoubleHypNucTree::GeoLength(const AliESDtrack& track) {
  Double_t deadZoneWidth = 3.0;
  Double_t lengthInActiveZone = track.GetLengthInActiveZone(1, deadZoneWidth, 220, track.GetESDEvent()->GetMagneticField(),0,0);
  return lengthInActiveZone;
}
