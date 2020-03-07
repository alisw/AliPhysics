
//--- Task for investigation of the DoubleHyperHydrogen4 ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

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
fInputHandler(0),
fPID(0),
fESDevent(0),
fStack(),
fV0(),
fHistdEdx(0),
fHistdEdxV0(0),
fHistNumEvents(0),
fHistTrigger(0),
fHistV0(0),
aTree(0),
bTree(0),
cTree(0),
dTree(0),
eTree(0),
fTree(0),
gTree(0),
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
//++ DCAs ++
fDCAHe4Pi(-99),
fDCAHe4P(-99),
fDCAPPi(-99),
fDCAHe3P(-99),
fDCAHe3Pi(-99),
fDCAHe3Pi1(-99),
fDCAHe3d(-99),
fDCAdPi(-99),
fDCATPi(-99),
fDCATP(-99),
//++ 5LHe ++
fz5LHe(-99),
fm5LHe(-99),
fp5LHe(-99),
fpt5LHe(-99),
fct5LHe(-99),
fy5LHe(-99),
fPA5LHe(-99),
//++ 4LHe ++
fz4LHe(-99),
fm4LHe(-99),
fp4LHe(-99),
fpt4LHe(-99),
fct4LHe(-99),
fy4LHe(-99),
fPA4LHe(-99),
//++ 4LH ++
fz4LH(-99),
fm4LH(-99),
fp4LH(-99),
fpt4LH(-99),
fct4LH(-99),
fy4LH(-99),
fPA4LH(-99),
//++ 4Li ++
fzLi4(-99),
fmLi4(-99),
fpLi4(-99),
fptLi4(-99),
fyLi4(-99),
//++ 4LLH ++
fmDaughterSum(-99), 
fpSum(-99), 
fctSum(-99), 
fptSum(-99), 
fySum(-99), 
fz4LLH(-99),
fPA4LLH(-99),
fAngle4LLH(-99),
//++ 4he ++
fhe4P(-99), 
fhe4DcaSec(-99), 
fhe4Ncls(-99), 
fhe4NclsITS(-99), 
fhe4Dca(-99), 
fhe4DedxSigma(-99), 
fhe4Dedx(-99), 
fEtaHe4(-99),
fPhiHe4(-99),
fGeoLengthHe4(-99),
fTOFSignalHe4(-99),
//++ 3he ++
fhe3P(-99), 
fhe3DcaSec(-99), 
fhe3Ncls(-99), 
fhe3NclsITS(-99), 
fhe3Dca(-99), 
fhe3DedxSigma(-99), 
fhe3Dedx(-99), 
fEtaHe3(-99),
fPhiHe3(-99),
fGeoLengthHe3(-99),
fTOFSignalHe3(-99),
//++ 3H ++
ftP(-99),
ftDcaSec(-99), 
ftDca(-99), 
ftNcls(-99), 
ftNclsITS(-99),
ftDedxSigma(-99),  
ftDedx(-99), 
fEtaT(-99),
fPhiT(-99),
fGeoLengthT(-99),
fTOFSignalT(-99),
//++ 2H ++
fdP(-99),
fdDcaSec(-99), 
fdDca(-99), 
fdNcls(-99), 
fdNclsITS(-99),
fdDedxSigma(-99),  
fdDedx(-99), 
fEtaD(-99),
fPhiD(-99),
fGeoLengthD(-99),
fTOFSignalD(-99),
//++ p ++
fpP(-99),
fpDcaSec(-99), 
fpDca(-99), 
fpNcls(-99), 
fpNclsITS(-99),
fpDedxSigma(-99),  
fpDedx(-99), 
fEtaP(-99),
fPhiP(-99),
fGeoLengthP(-99),
fTOFSignalP(-99),
//++ pi ++
fpiP(-99),
fpiDcaSec(-99), 
fpiDca(-99), 
fpiNcls(-99), 
fpiNclsITS(-99),
fpiDedxSigma(-99),  
fpiDedx(-99), 
fEtaPi(-99),
fPhiPi(-99),
fGeoLengthPi(-99),
fTOFSignalPi(-99),
//++ pi ++
fpi1P(-99),
fpi1DcaSec(-99), 
fpi1Dca(-99), 
fpi1Ncls(-99), 
fpi1NclsITS(-99),
fpi1DedxSigma(-99),  
fpi1Dedx(-99), 
fEtaPi1(-99),
fPhiPi1(-99),
fGeoLengthPi1(-99),
fTOFSignalPi1(-99),
//
farmalpha(-99), 
farmpt(-99),
ftrig(-99), 
fmc(-99), 
fthetaP(-99), 
fthetaN(-99),
fonTheFly(-99),
fVertexPosition(),
fNumberV0s(-99),      //< number of v0s in event
fCentrality(-99),     //< centrality of event
frunnumber(-99),      //< number of run
fTrigger(),        //< array of Triggers
fTriggerClasses() //< fired trigger classes
{

}

// Constructor
AliAnalysisTaskDoubleHypNucTree::AliAnalysisTaskDoubleHypNucTree(const char *name)
:AliAnalysisTaskSE(name),
fPIDCheckOnly(kFALSE),
ftrackAnalysis(kFALSE),
fInputHandler(0),
fPID(0),
fESDevent(0),
fStack(),
fV0(),
fHistdEdx(0),
fHistdEdxV0(0),
fHistNumEvents(0),
fHistTrigger(0),
fHistV0(0),
aTree(0),
bTree(0),
cTree(0),
dTree(0),
eTree(0),
fTree(0),
gTree(0),
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
//++ DCAs ++
fDCAHe4Pi(-99),
fDCAHe4P(-99),
fDCAPPi(-99),
fDCAHe3P(-99),
fDCAHe3Pi(-99),
fDCAHe3Pi1(-99),
fDCAHe3d(-99),
fDCAdPi(-99),
fDCATPi(-99),
fDCATP(-99),
//++ 5LHe ++
fz5LHe(-99),
fm5LHe(-99),
fp5LHe(-99),
fpt5LHe(-99),
fct5LHe(-99),
fy5LHe(-99),
fPA5LHe(-99),
//++ 4LHe ++
fz4LHe(-99),
fm4LHe(-99),
fp4LHe(-99),
fpt4LHe(-99),
fct4LHe(-99),
fy4LHe(-99),
fPA4LHe(-99),
//++ 4LH ++
fz4LH(-99),
fm4LH(-99),
fp4LH(-99),
fpt4LH(-99),
fct4LH(-99),
fy4LH(-99),
fPA4LH(-99),
//++ 4Li ++
fzLi4(-99),
fmLi4(-99),
fpLi4(-99),
fptLi4(-99),
fyLi4(-99),
//++ 4LLH ++
fmDaughterSum(-99), 
fpSum(-99), 
fctSum(-99), 
fptSum(-99), 
fySum(-99), 
fz4LLH(-99),
fPA4LLH(-99),
fAngle4LLH(-99),
//++ 4he ++
fhe4P(-99), 
fhe4DcaSec(-99), 
fhe4Ncls(-99), 
fhe4NclsITS(-99), 
fhe4Dca(-99), 
fhe4DedxSigma(-99), 
fhe4Dedx(-99), 
fEtaHe4(-99),
fPhiHe4(-99),
fGeoLengthHe4(-99),
fTOFSignalHe4(-99),
//++ 3he ++
fhe3P(-99), 
fhe3DcaSec(-99), 
fhe3Ncls(-99), 
fhe3NclsITS(-99), 
fhe3Dca(-99), 
fhe3DedxSigma(-99), 
fhe3Dedx(-99), 
fEtaHe3(-99),
fPhiHe3(-99),
fGeoLengthHe3(-99),
fTOFSignalHe3(-99),
//++ 3H ++
ftP(-99),
ftDcaSec(-99), 
ftDca(-99), 
ftNcls(-99), 
ftNclsITS(-99),
ftDedxSigma(-99),  
ftDedx(-99), 
fEtaT(-99),
fPhiT(-99),
fGeoLengthT(-99),
fTOFSignalT(-99),
//++ 2H ++
fdP(-99),
fdDcaSec(-99), 
fdDca(-99), 
fdNcls(-99), 
fdNclsITS(-99),
fdDedxSigma(-99),  
fdDedx(-99), 
fEtaD(-99),
fPhiD(-99),
fGeoLengthD(-99),
fTOFSignalD(-99),
//++ p ++
fpP(-99),
fpDcaSec(-99), 
fpDca(-99), 
fpNcls(-99), 
fpNclsITS(-99),
fpDedxSigma(-99),  
fpDedx(-99), 
fEtaP(-99),
fPhiP(-99),
fGeoLengthP(-99),
fTOFSignalP(-99),
//++ pi ++
fpiP(-99),
fpiDcaSec(-99), 
fpiDca(-99), 
fpiNcls(-99), 
fpiNclsITS(-99),
fpiDedxSigma(-99),  
fpiDedx(-99), 
fEtaPi(-99),
fPhiPi(-99),
fGeoLengthPi(-99),
fTOFSignalPi(-99),
//++ pi ++
fpi1P(-99),
fpi1DcaSec(-99), 
fpi1Dca(-99), 
fpi1Ncls(-99), 
fpi1NclsITS(-99),
fpi1DedxSigma(-99),  
fpi1Dedx(-99), 
fEtaPi1(-99),
fPhiPi1(-99),
fGeoLengthPi1(-99),
fTOFSignalPi1(-99),
//
farmalpha(-99), 
farmpt(-99),
ftrig(-99), 
fmc(-99), 
fthetaP(-99), 
fthetaN(-99),
fonTheFly(-99),
fVertexPosition(),
fNumberV0s(-99),      //< number of v0s in event
fCentrality(-99),     //< centrality of event
frunnumber(-99),      //< number of run
fTrigger(),        //< array of Triggers
fTriggerClasses() //< fired trigger classes
{
	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
	DefineOutput(2, TTree::Class());
	DefineOutput(3, TTree::Class());
	DefineOutput(4, TTree::Class());
	DefineOutput(5, TTree::Class());
	DefineOutput(6, TTree::Class());
	DefineOutput(7, TTree::Class());
	DefineOutput(8, TTree::Class());
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
	//Histograms
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
  	//++ TREE for 4LH ++
	aTree = new TTree("tree4LH","aTree");
	aTree->Branch("fm4LH", &fm4LH, "fm4LH/F");
	aTree->Branch("fz4LH", &fz4LH, "fz4LH/F");
	aTree->Branch("fp4LH", &fp4LH, "f4LH/F");
	aTree->Branch("fpt4LH", &fpt4LH, "fpt4LH/F");
	aTree->Branch("fct4LH", &fct4LH, "fct4LH/F");
	aTree->Branch("fy4LH", &fy4LH, "fy4LH/F");
	aTree->Branch("fPA4LH", &fPA4LH, "fPA4LH/F");
	aTree->Branch("fDCAHe4Pi", &fDCAHe4Pi, "fDCAHe4Pi/F");

	aTree->Branch("fhe4P", &fhe4P, "fhe4P/F");
	aTree->Branch("fhe4Dca", &fhe4Dca ,"fhe4Dca/F");
	aTree->Branch("fhe4DcaSec", &fhe4DcaSec ,"fhe4DcaSec/F");
	aTree->Branch("fhe4Ncls", &fhe4Ncls, "fhe4Ncls/F");
	aTree->Branch("fhe4NclsITS", &fhe4NclsITS, "fhe4NclsITS/F");
	aTree->Branch("fhe4DedxSigma", &fhe4DedxSigma, "fhe4DedxSigma/F");
	aTree->Branch("fhe4Dedx", &fhe4Dedx, "fhe4Dedx/F");
	aTree->Branch("fEtaHe4", &fEtaHe4, "fEtaHe4/F");
	aTree->Branch("fPhiHe4", &fPhiHe4, "fPhiHe4/F");
	aTree->Branch("fGeoLengthHe4", &fGeoLengthHe4, "fGeoLengthHe4/F");
	aTree->Branch("fTOFSignalHe4", &fTOFSignalHe4, "fTOFSignalHe4/F");

	aTree->Branch("fpiP", &fpiP, "fpiP/F");
	aTree->Branch("fpiDca", &fpiDca, "fpiDca/F");
	aTree->Branch("fpiDcaSec", &fpiDcaSec, "fpiDcaSec/F");
	aTree->Branch("fpiNcls", &fpiNcls, "fpiNcls/F");
	aTree->Branch("fpiNclsITS", &fpiNclsITS, "fpiNclsITS/F");
	aTree->Branch("fpiDedxSigma", &fpiDedxSigma, "fpiDedxSigma/F");
	aTree->Branch("fpiDedx", &fpiDedx, "fpiDedx/F");
	aTree->Branch("fEtaPi", &fEtaPi, "fEtaPi/F");
	aTree->Branch("fPhiPi", &fPhiPi, "fPhiPi/F");
	aTree->Branch("fGeoLengthPi", &fGeoLengthPi, "fGeoLengthPi/F");
	aTree->Branch("fTOFSignalPi", &fTOFSignalPi, "fTOFSignalPi/F");

	aTree->Branch("farmalpha", &farmalpha, "farmalpha/F");
	aTree->Branch("farmpt", &farmpt, "farmpt/F");
	aTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
	aTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
	aTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
	aTree->Branch("fmc", &fmc, "fmc/I");
	//++ Tree for 4Li ++
	bTree = new TTree("tree4Li","bTree");
	bTree->Branch("fmLi4", &fmLi4, "fmLi4/F");
	bTree->Branch("fzLi4", &fzLi4, "fzLi4/F");
	bTree->Branch("fpLi4", &fpLi4, "fpLi4/F");
	bTree->Branch("fptLi4", &fptLi4, "fptLi4/F");
	bTree->Branch("fyLi4", &fyLi4, "fyLi4/F");
	bTree->Branch("fDCAHe3P", &fDCAHe3P, "fDCAHe3P/F");

	bTree->Branch("fhe3P", &fhe3P, "fhe3P/F");
	bTree->Branch("fhe3Dca", &fhe3Dca ,"fhe3Dca/F");
	bTree->Branch("fhe3Ncls", &fhe3Ncls, "fhe3Ncls/F");
	bTree->Branch("fhe3NclsITS", &fhe3NclsITS, "fhe3NclsITS/F");
	bTree->Branch("fhe3DedxSigma", &fhe3DedxSigma, "fhe3DedxSigma/F");
	bTree->Branch("fhe3Dedx", &fhe3Dedx, "fhe3Dedx/F");
	bTree->Branch("fEtaHe3", &fEtaHe3, "fEtaHe3/F");
	bTree->Branch("fPhiHe3", &fPhiHe3, "fPhiHe3/F");
	bTree->Branch("fGeoLengthHe3", &fGeoLengthHe3, "fGeoLengthHe3/F");
	bTree->Branch("fTOFSignalHe3", &fTOFSignalHe3, "fTOFSignalHe3/F");

	bTree->Branch("fpP", &fpP, "fpP/F");
	bTree->Branch("fpDca", &fpDca ,"fpDca/F");
	bTree->Branch("fpNcls", &fpNcls, "fpNcls/F");
	bTree->Branch("fpNclsITS", &fpNclsITS, "fpNclsITS/F");
	bTree->Branch("fpDedxSigma", &fpDedxSigma, "fpDedxSigma/F");
	bTree->Branch("fpDedx", &fpDedx, "fpDedx/F");
	bTree->Branch("fEtaP", &fEtaP, "fEtaP/F");
	bTree->Branch("fPhiP", &fPhiP, "fPhiP/F");
	bTree->Branch("fGeoLengthP", &fGeoLengthP, "fGeoLengthP/F");
	bTree->Branch("fTOFSignalP", &fTOFSignalP, "fTOFSignalP/F");

	bTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
	bTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
	bTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
	bTree->Branch("fmc", &fmc, "fmc/I");
	//++ Tree for 5LHe (4He,p,pi) ++
	cTree = new TTree("tree5LHe","cTree");
	cTree->Branch("fm5LHe", &fm5LHe, "fm5LHe/F");
	cTree->Branch("fz5LHe", &fz5LHe, "fz5LHe/F");
	cTree->Branch("fp5LHe", &fp5LHe, "fp5LHe/F");
	cTree->Branch("fpt5LHe", &fpt5LHe, "fpt5LHe/F");
	cTree->Branch("fy5LHe", &fy5LHe, "fy5LHe/F");
	cTree->Branch("fct5LHe", &fct5LHe, "fct5LHe/F");

	cTree->Branch("fPA5LHe", &fPA5LHe, "fPA5LHe/F");
	cTree->Branch("fDCAHe4P", &fDCAHe4P, "fDCAHe4P/F");
	cTree->Branch("fDCAHe4Pi", &fDCAHe4Pi, "fDCAHe4Pi/F");
	cTree->Branch("fDCAPPi", &fDCAPPi, "fDCAPPi/F");

	cTree->Branch("fhe4P", &fhe4P, "fhe4P/F");
	cTree->Branch("fhe4Dca", &fhe4Dca ,"fhe4Dca/F");
	cTree->Branch("fhe4DcaSec", &fhe4DcaSec ,"fhe4DcaSec/F");
	cTree->Branch("fhe4Ncls", &fhe4Ncls, "fhe4Ncls/F");
	cTree->Branch("fhe4NclsITS", &fhe4NclsITS, "fhe4NclsITS/F");
	cTree->Branch("fhe4DedxSigma", &fhe4DedxSigma, "fhe4DedxSigma/F");
	cTree->Branch("fhe4Dedx", &fhe4Dedx, "fhe4Dedx/F");
	cTree->Branch("fEtaHe4", &fEtaHe4, "fEtaHe4/F");
	cTree->Branch("fPhiHe4", &fPhiHe4, "fPhiHe4/F");
	cTree->Branch("fGeoLengthHe4", &fGeoLengthHe4, "fGeoLengthHe4/F");
	cTree->Branch("fTOFSignalHe4", &fTOFSignalHe4, "fTOFSignalHe4/F");

	cTree->Branch("fpP", &fpP, "fpP/F");
	cTree->Branch("fpDca", &fpDca ,"fpDca/F");
	cTree->Branch("fpDcaSec", &fpDcaSec ,"fpDcaSec/F");
	cTree->Branch("fpNcls", &fpNcls, "fpNcls/F");
	cTree->Branch("fpNclsITS", &fpNclsITS, "fpNclsITS/F");
	cTree->Branch("fpDedxSigma", &fpDedxSigma, "fpDedxSigma/F");
	cTree->Branch("fpDedx", &fpDedx, "fpDedx/F");
	cTree->Branch("fEtaP", &fEtaP, "fEtaP/F");
	cTree->Branch("fPhiP", &fPhiP, "fPhiP/F");
	cTree->Branch("fGeoLengthP", &fGeoLengthP, "fGeoLengthP/F");
	cTree->Branch("fTOFSignalP", &fTOFSignalP, "fTOFSignalP/F");

	cTree->Branch("fpiP", &fpiP, "fpiP/F");
	cTree->Branch("fpiDca", &fpiDca, "fpiDca/F");
	cTree->Branch("fpiDcaSec", &fpiDcaSec, "fpiDcaSec/F");
	cTree->Branch("fpiNcls", &fpiNcls, "fpiNcls/F");
	cTree->Branch("fpiNclsITS", &fpiNclsITS, "fpiNclsITS/F");
	cTree->Branch("fpiDedxSigma", &fpiDedxSigma, "fpiDedxSigma/F");
	cTree->Branch("fpiDedx", &fpiDedx, "fpiDedx/F");
	cTree->Branch("fEtaPi", &fEtaPi, "fEtaPi/F");
	cTree->Branch("fPhiPi", &fPhiPi, "fPhiPi/F");
	cTree->Branch("fGeoLengthPi", &fGeoLengthPi, "fGeoLengthPi/F");
	cTree->Branch("fTOFSignalPi", &fTOFSignalPi, "fTOFSignalPi/F");

	cTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
	cTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
	cTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
	cTree->Branch("fmc", &fmc, "fmc/I");
	//++ Tree for 5LHe (3He,d,pi) ++
	dTree = new TTree("tree5LHe2","dTree");
	dTree->Branch("fm5LHe", &fm5LHe, "fm5LHe/F");
	dTree->Branch("fz5LHe", &fz5LHe, "fz5LHe/F");
	dTree->Branch("fp5LHe", &fp5LHe, "fp5LHe/F");
	dTree->Branch("fpt5LHe", &fpt5LHe, "fpt5LHe/F");
	dTree->Branch("fy5LHe", &fy5LHe, "fy5LHe/F");
	dTree->Branch("fct5LHe", &fct5LHe, "fct5LHe/F");

	dTree->Branch("fPA5LHe", &fPA5LHe, "fPA5LHe/F");
	dTree->Branch("fDCAHe3d", &fDCAHe3d, "fDCAHe3d/F");
	dTree->Branch("fDCAHe3Pi", &fDCAHe3Pi, "fDCAHe3Pi/F");
	dTree->Branch("fDCAdPi", &fDCAdPi, "fDCAdPi/F");

	dTree->Branch("fhe3P", &fhe3P, "fhe3P/F");
	dTree->Branch("fhe3Dca", &fhe3Dca ,"fhe3Dca/F");
	dTree->Branch("fhe3DcaSec", &fhe3DcaSec ,"fhe3DcaSec/F");
	dTree->Branch("fhe3Ncls", &fhe3Ncls, "fhe3Ncls/F");
	dTree->Branch("fhe3NclsITS", &fhe3NclsITS, "fhe3NclsITS/F");
	dTree->Branch("fhe3DedxSigma", &fhe3DedxSigma, "fhe3DedxSigma/F");
	dTree->Branch("fhe3Dedx", &fhe3Dedx, "fhe3Dedx/F");
	dTree->Branch("fEtaHe3", &fEtaHe3, "fEtaHe3/F");
	dTree->Branch("fPhiHe3", &fPhiHe3, "fPhiHe3/F");
	dTree->Branch("fGeoLengthHe3", &fGeoLengthHe3, "fGeoLengthHe3/F");
	dTree->Branch("fTOFSignalHe3", &fTOFSignalHe3, "fTOFSignalHe3/F");

	dTree->Branch("fdP", &fdP, "fdP/F");
	dTree->Branch("fdDca", &fdDca ,"fdDca/F");
	dTree->Branch("fdDcaSec", &fdDcaSec ,"fdDcaSec/F");
	dTree->Branch("fdNcls", &fdNcls, "fdNcls/F");
	dTree->Branch("fdNclsITS", &fdNclsITS, "fdNclsITS/F");
	dTree->Branch("fdDedxSigma", &fdDedxSigma, "fdDedxSigma/F");
	dTree->Branch("fdDedx", &fdDedx, "fdDedx/F");
	dTree->Branch("fEtaD", &fEtaD, "fEtaD/F");
	dTree->Branch("fPhiD", &fPhiD, "fPhiD/F");
	dTree->Branch("fGeoLengthD", &fGeoLengthD, "fGeoLengthD/F");
	dTree->Branch("fTOFSignalD", &fTOFSignalD, "fTOFSignalD/F");

	dTree->Branch("fpiP", &fpiP, "fpiP/F");
	dTree->Branch("fpiDca", &fpiDca, "fpiDca/F");
	dTree->Branch("fpiDcaSec", &fpiDcaSec, "fpiDcaSec/F");
	dTree->Branch("fpiNcls", &fpiNcls, "fpiNcls/F");
	dTree->Branch("fpiNclsITS", &fpiNclsITS, "fpiNclsITS/F");
	dTree->Branch("fpiDedxSigma", &fpiDedxSigma, "fpiDedxSigma/F");
	dTree->Branch("fpiDedx", &fpiDedx, "fpiDedx/F");
	dTree->Branch("fEtaPi", &fEtaPi, "fEtaPi/F");
	dTree->Branch("fPhiPi", &fPhiPi, "fPhiPi/F");
	dTree->Branch("fGeoLengthPi", &fGeoLengthPi, "fGeoLengthPi/F");
	dTree->Branch("fTOFSignalPi", &fTOFSignalPi, "fTOFSignalPi/F");

	dTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
	dTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
	dTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
	dTree->Branch("fmc", &fmc, "fmc/I");
	//++ Tree for 4LHe (3He,p,pi) ++
	eTree = new TTree("tree4LHe","eTree");
	eTree->Branch("fm4LHe", &fm4LHe, "fm4LHe/F");
	eTree->Branch("fz4LHe", &fz4LHe, "fz4LHe/F");
	eTree->Branch("fp4LHe", &fp4LHe, "fp4LHe/F");
	eTree->Branch("fpt4LHe", &fpt4LHe, "fpt4LHe/F");
	eTree->Branch("fy4LHe", &fy4LHe, "fy4LHe/F");
	eTree->Branch("fct4LHe", &fct4LHe, "fct4LHe/F");

	eTree->Branch("fPA4LHe", &fPA4LHe, "fPA4LHe/F");
	eTree->Branch("fDCAHe3P", &fDCAHe3P, "fDCAHe3P/F");
	eTree->Branch("fDCAHe3Pi", &fDCAHe3Pi, "fDCAHe3Pi/F");
	eTree->Branch("fDCAPPi", &fDCAPPi, "fDCAPPi/F");

	eTree->Branch("fhe3P", &fhe3P, "fhe3P/F");
	eTree->Branch("fhe3Dca", &fhe3Dca ,"fhe3Dca/F");
	eTree->Branch("fhe3DcaSec", &fhe3DcaSec ,"fhe3DcaSec/F");
	eTree->Branch("fhe3Ncls", &fhe3Ncls, "fhe3Ncls/F");
	eTree->Branch("fhe3NclsITS", &fhe3NclsITS, "fhe3NclsITS/F");
	eTree->Branch("fhe3DedxSigma", &fhe3DedxSigma, "fhe3DedxSigma/F");
	eTree->Branch("fhe3Dedx", &fhe3Dedx, "fhe3Dedx/F");
	eTree->Branch("fEtaHe3", &fEtaHe3, "fEtaHe3/F");
	eTree->Branch("fPhiHe3", &fPhiHe3, "fPhiHe3/F");
	eTree->Branch("fGeoLengthHe3", &fGeoLengthHe3, "fGeoLengthHe3/F");
	eTree->Branch("fTOFSignalHe3", &fTOFSignalHe3, "fTOFSignalHe3/F");	

	eTree->Branch("fpP", &fpP, "fpP/F");
	eTree->Branch("fpDca", &fpDca ,"fpDca/F");
	eTree->Branch("fpDcaSec", &fpDcaSec ,"fpDcaSec/F");
	eTree->Branch("fpNcls", &fpNcls, "fpNcls/F");
	eTree->Branch("fpNclsITS", &fpNclsITS, "fpNclsITS/F");
	eTree->Branch("fpDedxSigma", &fpDedxSigma, "fpDedxSigma/F");
	eTree->Branch("fpDedx", &fpDedx, "fpDedx/F");
	eTree->Branch("fEtaP", &fEtaP, "fEtaP/F");
	eTree->Branch("fPhiP", &fPhiP, "fPhiP/F");
	eTree->Branch("fGeoLengthP", &fGeoLengthP, "fGeoLengthP/F");
	eTree->Branch("fTOFSignalP", &fTOFSignalP, "fTOFSignalP/F");

	eTree->Branch("fpiP", &fpiP, "fpiP/F");
	eTree->Branch("fpiDca", &fpiDca, "fpiDca/F");
	eTree->Branch("fpiDcaSec", &fpiDcaSec, "fpiDcaSec/F");
	eTree->Branch("fpiNcls", &fpiNcls, "fpiNcls/F");
	eTree->Branch("fpiNclsITS", &fpiNclsITS, "fpiNclsITS/F");
	eTree->Branch("fpiDedxSigma", &fpiDedxSigma, "fpiDedxSigma/F");
	eTree->Branch("fpiDedx", &fpiDedx, "fpiDedx/F");
	eTree->Branch("fEtaPi", &fEtaPi, "fEtaPi/F");
	eTree->Branch("fPhiPi", &fPhiPi, "fPhiPi/F");
	eTree->Branch("fGeoLengthPi", &fGeoLengthPi, "fGeoLengthPi/F");
	eTree->Branch("fTOFSignalPi", &fTOFSignalPi, "fTOFSignalPi/F");

	eTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
	eTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
	eTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
	eTree->Branch("fmc", &fmc, "fmc/I");
	//++ Tree for 4LH (t,p,pi) ++
	fTree = new TTree("tree4LH3B","fTree");
	fTree->Branch("fm4LH", &fm4LH, "fm4LH/F");
	fTree->Branch("fz4LH", &fz4LH, "fz4LH/F");
	fTree->Branch("fp4LH", &fp4LH, "f4LH/F");
	fTree->Branch("fpt4LH", &fpt4LH, "fpt4LH/F");
	fTree->Branch("fct4LH", &fct4LH, "fct4LH/F");
	fTree->Branch("fy4LH", &fy4LH, "fy4LH/F");

	fTree->Branch("fPA4LH", &fPA4LH, "fPA4LH/F");
	fTree->Branch("fDCATP", &fDCATP, "fDCATP/F");
	fTree->Branch("fDCATPi", &fDCAHe3Pi, "fDCATPi/F");
	fTree->Branch("fDCAPPi", &fDCAPPi, "fDCAPPi/F");

	fTree->Branch("ftP", &ftP, "ftP/F");
	fTree->Branch("ftDca", &ftDca ,"ftDca/F");
	fTree->Branch("ftDcaSec", &ftDcaSec ,"ftDcaSec/F");
	fTree->Branch("ftNcls", &ftNcls, "ftNcls/F");
	fTree->Branch("ftNclsITS", &ftNclsITS, "ftNclsITS/F");
	fTree->Branch("ftDedxSigma", &ftDedxSigma, "ftDedxSigma/F");
	fTree->Branch("ftDedx", &ftDedx, "ftDedx/F");
	fTree->Branch("fEtaT", &fEtaT, "fEtaT/F");
	fTree->Branch("fPhiT", &fPhiT, "fPhiT/F");
	fTree->Branch("fGeoLengthT", &fGeoLengthT, "fGeoLengthT/F");
	fTree->Branch("fTOFSignalT", &fTOFSignalT, "fTOFSignalT/F");

	fTree->Branch("fpP", &fpP, "fpP/F");
	fTree->Branch("fpDca", &fpDca ,"fpDca/F");
	fTree->Branch("fpDcaSec", &fpDcaSec ,"fpDcaSec/F");
	fTree->Branch("fpNcls", &fpNcls, "fpNcls/F");
	fTree->Branch("fpNclsITS", &fpNclsITS, "fpNclsITS/F");
	fTree->Branch("fpDedxSigma", &fpDedxSigma, "fpDedxSigma/F");
	fTree->Branch("fpDedx", &fpDedx, "fpDedx/F");
	fTree->Branch("fEtaP", &fEtaP, "fEtaP/F");
	fTree->Branch("fPhiP", &fPhiP, "fPhiP/F");
	fTree->Branch("fGeoLengthP", &fGeoLengthP, "fGeoLengthP/F");
	fTree->Branch("fTOFSignalP", &fTOFSignalP, "fTOFSignalP/F");

	fTree->Branch("fpiP", &fpiP, "fpiP/F");
	fTree->Branch("fpiDca", &fpiDca, "fpiDca/F");
	fTree->Branch("fpiDcaSec", &fpiDcaSec, "fpiDcaSec/F");
	fTree->Branch("fpiNcls", &fpiNcls, "fpiNcls/F");
	fTree->Branch("fpiNclsITS", &fpiNclsITS, "fpiNclsITS/F");
	fTree->Branch("fpiDedxSigma", &fpiDedxSigma, "fpiDedxSigma/F");
	fTree->Branch("fpiDedx", &fpiDedx, "fpiDedx/F");
	fTree->Branch("fEtaPi", &fEtaPi, "fEtaPi/F");
	fTree->Branch("fPhiPi", &fPhiPi, "fPhiPi/F");
	fTree->Branch("fGeoLengthPi", &fGeoLengthPi, "fGeoLengthPi/F");
	fTree->Branch("fTOFSignalPi", &fTOFSignalPi, "fTOFSignalPi/F");

	fTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
	fTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
	fTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
	fTree->Branch("fmc", &fmc, "fmc/I");
	//++ Tree for 4LLH (3He,p,pi,pi) ++
	gTree = new TTree("tree4LLH","gTree");
	gTree->Branch("fmDaughterSum", &fmDaughterSum, "fmDaughterSum/F");
	gTree->Branch("fpSum", &fpSum, "fpSum/F");
	gTree->Branch("fptSum", &fptSum, "fptSum/F");
	gTree->Branch("fctSum", &fctSum, "fctSum/F");
	gTree->Branch("fySum", &fySum, "fySum/F");
	gTree->Branch("fz4LLH", &fz4LLH, "fz4LLH/F");

	gTree->Branch("fm4LHe", &fm4LHe, "fm4LHe/F");
	gTree->Branch("fp4LHe", &fp4LHe, "fp4LHe/F");
	gTree->Branch("fpt4LHe", &fpt4LHe, "fpt4LHe/F");
	gTree->Branch("fy4LHe", &fy4LHe, "fy4LHe/F");
	gTree->Branch("fct4LHe", &fct4LHe, "fct4LHe/F");

	gTree->Branch("fDCAHe3P", &fDCAHe3P, "fDCAHe3P/F");
	gTree->Branch("fDCAHe3Pi", &fDCAHe3Pi, "fDCAHe3Pi/F");
	gTree->Branch("fDCAPPi", &fDCAPPi, "fDCAPPi/F");
	gTree->Branch("fDCAHe3Pi1", &fDCAHe3Pi1, "fDCAHe3Pi1/F");
	gTree->Branch("fPA4LLH", &fPA4LLH, "fPA4LLH/F");
	gTree->Branch("fAngle4LLH", &fAngle4LLH, "fAngle4LLH/F");
	
	gTree->Branch("fhe3P", &fhe3P, "fhe3P/F");
	gTree->Branch("fhe3Dca", &fhe3Dca ,"fhe3Dca/F");
	gTree->Branch("fhe3DcaSec", &fhe3DcaSec ,"fhe3DcaSec/F");
	gTree->Branch("fhe3Ncls", &fhe3Ncls, "fhe3Ncls/F");
	gTree->Branch("fhe3NclsITS", &fhe3NclsITS, "fhe3NclsITS/F");
	gTree->Branch("fhe3DedxSigma", &fhe3DedxSigma, "fhe3DedxSigma/F");
	gTree->Branch("fhe3Dedx", &fhe3Dedx, "fhe3Dedx/F");
	gTree->Branch("fEtaHe3", &fEtaHe3, "fEtaHe3/F");
	gTree->Branch("fPhiHe3", &fPhiHe3, "fPhiHe3/F");
	gTree->Branch("fGeoLengthHe3", &fGeoLengthHe3, "fGeoLengthHe3/F");
	gTree->Branch("fTOFSignalHe3", &fTOFSignalHe3, "fTOFSignalHe3/F");

	gTree->Branch("fpP", &fpP, "fpP/F");
	gTree->Branch("fpDca", &fpDca ,"fpDca/F");
	gTree->Branch("fpDcaSec", &fpDcaSec ,"fpDcaSec/F");
	gTree->Branch("fpNcls", &fpNcls, "fpNcls/F");
	gTree->Branch("fpNclsITS", &fpNclsITS, "fpNclsITS/F");
	gTree->Branch("fpDedxSigma", &fpDedxSigma, "fpDedxSigma/F");
	gTree->Branch("fpDedx", &fpDedx, "fpDedx/F");
	gTree->Branch("fEtaP", &fEtaP, "fEtaP/F");
	gTree->Branch("fPhiP", &fPhiP, "fPhiP/F");
	gTree->Branch("fGeoLengthP", &fGeoLengthP, "fGeoLengthP/F");
	gTree->Branch("fTOFSignalP", &fTOFSignalP, "fTOFSignalP/F");

	gTree->Branch("fpiP", &fpiP, "fpiP/F");
	gTree->Branch("fpiDca", &fpiDca, "fpiDca/F");
	gTree->Branch("fpiDcaSec", &fpiDcaSec, "fpiDcaSec/F");
	gTree->Branch("fpiNcls", &fpiNcls, "fpiNcls/F");
	gTree->Branch("fpiNclsITS", &fpiNclsITS, "fpiNclsITS/F");
	gTree->Branch("fpiDedxSigma", &fpiDedxSigma, "fpiDedxSigma/F");
	gTree->Branch("fpiDedx", &fpiDedx, "fpiDedx/F");
	gTree->Branch("fEtaPi", &fEtaPi, "fEtaPi/F");
	gTree->Branch("fPhiPi", &fPhiPi, "fPhiPi/F");
	gTree->Branch("fGeoLengthPi", &fGeoLengthPi, "fGeoLengthPi/F");
	gTree->Branch("fTOFSignalPi", &fTOFSignalPi, "fTOFSignalPi/F");

	gTree->Branch("fpi1P", &fpi1P, "fpi1P/F");
	gTree->Branch("fpi1Dca", &fpi1Dca, "fpi1Dca/F");
	gTree->Branch("fpi1DcaSec", &fpi1DcaSec, "fpi1DcaSec/F");
	gTree->Branch("fpi1Ncls", &fpi1Ncls, "fpi1Ncls/F");
	gTree->Branch("fpi1NclsITS", &fpi1NclsITS, "fpi1NclsITS/F");
	gTree->Branch("fpi1DedxSigma", &fpi1DedxSigma, "fpi1DedxSigma/F");
	gTree->Branch("fpi1Dedx", &fpi1Dedx, "fpi1Dedx/F");
	gTree->Branch("fEtaPi1", &fEtaPi1, "fEtaPi1/F");
	gTree->Branch("fPhiPi1", &fPhiPi1, "fPhiPi1/F");
	gTree->Branch("fGeoLengthPi1", &fGeoLengthPi1, "fGeoLengthPi1/F");
	gTree->Branch("fTOFSignalPi1", &fTOFSignalPi1, "fTOFSignalPi1/F");

	gTree->Branch("farmalpha", &farmalpha, "farmalpha/F");
	gTree->Branch("farmpt", &farmpt, "farmpt/F");
	gTree->Branch("frunnumber", &frunnumber,"frunnumber/I");
	gTree->Branch("fNumberV0s", &fNumberV0s, "fNumberV0s/I");
	gTree->Branch("fCentrality", &fCentrality, "fCentrality/I");
	gTree->Branch("fmc", &fmc, "fmc/I");

	PostData(1, fHistogramList);
	PostData(2, aTree);
	PostData(3, bTree);
	PostData(4, cTree);
	PostData(5, dTree);
	PostData(6, eTree);
	PostData(7, fTree);
	PostData(8, gTree);

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
  //++ MC ++
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
  //++ Data ++
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
    //++ 0 = V0M ++
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
  //++ 0 = V0M ++
		centrality = fEventCuts.GetCentrality(0);
	}
	//++ sec vertex reconstruction ++
	AliESDVertex *esdVer1 = new AliESDVertex(*vertex);
	AliVertexerTracks *vertexer = new AliVertexerTracks(fESDevent->GetMagneticField());
	vertexer->SetVtxStart(esdVer1);
	//++ arrays for pointing angle calculation ++
	Double_t *dn = new Double_t[3];
	Double_t *dd = new Double_t[3];
	dn[0] = esdVer1->GetX();
	dn[1] = esdVer1->GetY();
	dn[2] = esdVer1->GetZ();
	if(esdVer1) delete esdVer1;
	//++ runnumber ++
	Int_t runNumber = fESDevent->GetRunNumber();
	frunnumber = runNumber;
	TriggerSelection();
	//++ Number of Events ++
	fHistNumEvents->Fill(1);
	//++ Centrality ++
	fCentrality = centrality;
	//++ MagneticField ++
	fMagneticField  = fESDevent->GetMagneticField();
	//++ Primary Vertex Position ++
	fPrimaryVertex.SetXYZ(vertex->GetX(),vertex->GetY(),vertex->GetZ());
	fVertexPosition = fPrimaryVertex; 
	//++ V0 ++
	fNV0Cand = 0;
	//++ For DCA ++
	Double_t xthiss(0.0);
	Double_t xpp(0.0);

	fmc = fMCtrue;

	AliESDtrackCuts trackCutsNuc("AlitrackCutsNuc", "AlitrackCutsNuc");
	AliESDtrackCuts trackCutsP("AlitrackCutsP", "AlitrackCutsP");
	AliESDtrackCuts trackCutsPi("AlitrackCutsPi", "AlitrackCutsPi");
	
	if(ftrackAnalysis){
		//++ Track 1 ++
		trackCutsNuc.SetEtaRange(-0.9,0.9);
		trackCutsNuc.SetAcceptKinkDaughters(kFALSE);
		trackCutsNuc.SetRequireTPCRefit(kTRUE);
		trackCutsNuc.SetMaxChi2PerClusterTPC(5);
		trackCutsNuc.SetMinNClustersTPC(60);
		trackCutsNuc.SetMaxRel1PtUncertainty(0.1);
		trackCutsNuc.SetPtRange(0.0, 10.0);
		//++ Track 2 ++
		trackCutsP.SetEtaRange(-0.9,0.9);
		trackCutsP.SetAcceptKinkDaughters(kFALSE);
		trackCutsP.SetRequireTPCRefit(kTRUE);
		trackCutsP.SetMaxChi2PerClusterTPC(5);
		trackCutsP.SetMinNClustersTPC(60);
		trackCutsP.SetMaxRel1PtUncertainty(0.1);
		trackCutsP.SetPtRange(0.0, 5.0);
		//trackCutsP.SetMinDCAToVertexXY(0.05);
  		//trackCutsP.SetMinDCAToVertexZ(0.05);
		//++ Track 3 & Track 4 ++
		trackCutsPi.SetEtaRange(-0.9,0.9);
		trackCutsPi.SetAcceptKinkDaughters(kFALSE);
		trackCutsPi.SetRequireTPCRefit(kTRUE);
		trackCutsPi.SetMaxChi2PerClusterTPC(5);
		trackCutsPi.SetMinNClustersTPC(60);
		trackCutsPi.SetMaxRel1PtUncertainty(0.2);
		trackCutsPi.SetPtRange(0.0, 1.0);
		//trackCutsP.SetMinDCAToVertexXY(0.1);
  		//trackCutsP.SetMinDCAToVertexZ(0.1);
	}	
	//++ Pidqa loop ++
	if(fPIDCheckOnly){
		dEdxCheck();
	}
	//++ use Tracks to reconstruct 4LLH, 4LHe, 4LH, 5LHe ++
	if(ftrackAnalysis){
		TrackAnalysis(trackCutsNuc, trackCutsP, trackCutsPi, vertexer, dn, dd, xthiss, xpp);
	}

	PostData(1, fHistogramList);
	PostData(2, aTree);
	PostData(3, bTree);
	PostData(4, cTree);
	PostData(5, dTree);
	PostData(6, eTree);
	PostData(7, fTree);
	PostData(8, gTree);
	}
//------------------------------------------------
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
//------------------------------------------------
void AliAnalysisTaskDoubleHypNucTree::TrackAnalysis(AliESDtrackCuts trackCutsNuc, AliESDtrackCuts trackCutsP, AliESDtrackCuts trackCutsPi, AliVertexerTracks *vertexer, Double_t *dn, Double_t *dd, Double_t xthiss, Double_t xpp){
	//++ Lorentz Vectors ++
	TLorentzVector *vecFirstDaughter = new TLorentzVector(0.,0.,0.,0.);
	TLorentzVector *vecSecDaughter = new TLorentzVector(0.,0.,0.,0.);
	TLorentzVector *vecThirdDaughter = new TLorentzVector(0.,0.,0.,0.);
	TLorentzVector *vecFourthDaughter = new TLorentzVector(0.,0.,0.,0.);
	TLorentzVector *vecMother = new TLorentzVector(0.,0.,0.,0.);
	TVector3 *h = new TVector3(0., 0., 0.);
	//++ PID Bools ++
	Bool_t He4Pos1 = kFALSE;
	Bool_t He3Pos1 = kFALSE;
	Bool_t tPos1 = kFALSE;

	Bool_t tPos2 = kFALSE;
	Bool_t dPos2 = kFALSE;
	Bool_t pPos2 = kFALSE;
	Bool_t piPos2 = kFALSE;

	Bool_t piPos3 = kFALSE;
	Bool_t piPos4 = kFALSE;

	Bool_t He4Neg1 = kFALSE;
	Bool_t He3Neg1 = kFALSE;
	Bool_t tNeg1 = kFALSE;

	Bool_t tNeg2 = kFALSE;
	Bool_t dNeg2 = kFALSE;
	Bool_t pNeg2 = kFALSE;
	Bool_t piNeg2 = kFALSE;

	Bool_t piNeg3 = kFALSE;
	Bool_t piNeg4 = kFALSE;
	//+++++
	Bool_t Check4LH3B = kFALSE; //4LH 3-Body on/off
	Bool_t ShowProgress = kTRUE; //special progress for checks
   	//++ Variables for TOF ++
	Float_t mass = 0;
	Float_t time = -1;
	Float_t beta = 0;
	Float_t gamma = 0;
	Float_t length = 0;
	Float_t time0 = 0;
    //++++++++++++++++++++++++++++ first track ++++++++++++++++++++++++++++++++++//
	for (Int_t iTracks = 0; iTracks < fESDevent->GetNumberOfTracks(); iTracks++) {

		if(ShowProgress){
			if(iTracks == 0) cout<<"### Starting Analysis with "<<fESDevent->GetNumberOfTracks()<<" Tracks ###"<<endl;
		}

		AliESDtrack* track1 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(iTracks));

		if (!track1->GetInnerParam()) continue;

		if (!trackCutsNuc.AcceptTrack(track1)) continue;    

		Double_t ptot1 = track1->GetInnerParam()->GetP();
		Double_t sign1 = track1->GetSign();

		fHistdEdx->Fill(ptot1*sign1, track1->GetTPCsignal());
		//++ Reset PID Bools ++
		He4Pos1 = kFALSE;
		He3Pos1 = kFALSE;
		tPos1 = kFALSE;
		He4Neg1 = kFALSE;
		He3Neg1 = kFALSE;
		tNeg1 = kFALSE;
		//++ Check PID ++
		if(fBetheSplines){
			if(sign1 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kAlpha)) <= 3) He4Pos1 = kTRUE;
			else if(sign1 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kAlpha)) <= 3) He4Neg1 = kTRUE;
			else if(sign1 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) <= 3) He3Pos1 = kTRUE;
			else if(sign1 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kHe3)) <= 3) He3Neg1 = kTRUE;
			else if(Check4LH3B && sign1 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kTriton)) <= 3) tPos1 = kTRUE;
			else if(Check4LH3B && sign1 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kTriton)) <= 3) tNeg1 = kTRUE;
			else continue;
			//if(!He4Pos1 && !He4Neg1 && !He3Pos1 && !He3Neg1 && !tPos1 && !tNeg1) continue;
		}
		else {		
			if(sign1 >0 && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe)) <= 3) He4Pos1 = kTRUE;
			else if(sign1 <0 && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe)) <= 3) He4Neg1 = kTRUE;
			else if(sign1 >0 && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) <= 3) He3Pos1 = kTRUE;
			else if(sign1 <0 && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) <= 3) He3Neg1 = kTRUE;
			else if(Check4LH3B && sign1 >0 && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT)) <= 3) tPos1 = kTRUE;
			else if(Check4LH3B && sign1 <0 && TMath::Abs(Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT)) <= 3) tNeg1 = kTRUE;
			else continue;
			//if(!He4Pos1 && !He4Neg1 && !He3Pos1 && !He3Neg1 && !tPos1 && !tNeg1) continue;
		}
		if(ShowProgress){
			if(He4Pos1) cout<<"<--Alpha-->"<<endl;
			if(He4Neg1) cout<<"<--Anti-Alpha-->"<<endl;
			if(He3Pos1) cout<<"<-- 3He -->"<<endl;
			if(He3Neg1) cout<<"<-- Anti-3He -->"<<endl;
		}
      //++++++++++++++++++++++++++++++++ second track +++++++++++++++++++++++++++++++++++++++//
		for (Int_t jTracks = iTracks + 1; jTracks < fESDevent->GetNumberOfTracks(); jTracks++) {

			AliESDtrack* track2 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(jTracks));              

			if (!track2->GetInnerParam()) continue;

			if (!trackCutsP.AcceptTrack(track2)) continue;

			if(iTracks == jTracks) continue;

			Double_t ptot2 = track2->GetInnerParam()->GetP();
			Double_t sign2 = track2->GetSign();

			fHistdEdx->Fill(ptot2*sign2, track2->GetTPCsignal());
			//++ Reset PID Bools ++
			tPos2 = kFALSE;
			dPos2 = kFALSE;
			pPos2 = kFALSE;
			piPos2 = kFALSE;
			tNeg2 = kFALSE;
			dNeg2 = kFALSE;
			pNeg2 = kFALSE;
			piNeg2 = kFALSE;
			//++ Check PID ++
			if(fBetheSplines){		
				if(sign2 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kDeuteron)) <= 3) dPos2 = kTRUE;
				else if(sign2 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kDeuteron)) <= 3) dNeg2 = kTRUE;			
				else if(sign2 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) <= 3) pPos2 = kTRUE;
				else if(sign2 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kProton)) <= 3) pNeg2 = kTRUE;
				else if(sign2 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kPion)) <= 3) piPos2 = kTRUE;
				else if(sign2 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kPion)) <= 3) piNeg2 = kTRUE;
				else continue;
				//else if(!dPos2 && !dNeg2 && !pPos2 && !pNeg2 && !piPos2 && !piNeg2) continue;
			}
			else {
				if(sign2 >0 && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT)) <= 3) dPos2 = kTRUE;
				else if(sign2 <0 && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT)) <= 3) dNeg2 = kTRUE;
				else if(sign2 >0 && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) <= 3) pPos2 = kTRUE;
				else if(sign2 <0 && TMath::Abs(Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT)) <= 3) pNeg2 = kTRUE;
				else if(sign2 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kPion)) <= 3) piPos2 = kTRUE;
				else if(sign2 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kPion)) <= 3) piNeg2 = kTRUE;
				else continue;
				//else if(!dPos2 && !dNeg2 && !pPos2 && !pNeg2 && !piPos2 && !piNeg2) continue;
			}
			//++ sign control for different cases ++
			if((He3Pos1 || He3Neg1) && sign1 != sign2) continue;
			else if((tPos1 || tNeg1) && sign1 != sign2) continue;
			else if((piPos2 || piNeg2) && sign1 == sign2) continue;
			else if((!piPos2 && !piNeg2) && sign1 != sign2) continue;
			//++ Check for the different expected combinations ++
			if(!(He4Pos1 && piNeg2) && !(He4Neg1 && piPos2) && !(He3Pos1 && pPos2) && !(He3Neg1 && pNeg2) && !(He3Pos1 && dPos2) && !(He3Neg1 && dNeg2) && !(He4Pos1 && pPos2) && !(He4Neg1 && pNeg2) && !(tPos1 && pPos2) && !(tNeg1 && pNeg2)) continue;
			//++ pre DCA cut ++
			Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
			if(*dca > 1.0) continue;
			if(dca) delete dca;
        	//************ 4LH Pos **************
			if(He4Pos1 && piNeg2){
				//cout<<"4LH"<<endl;
				TObjArray *trkArray = new TObjArray(2);                  
				trkArray->AddAt(track1,0);
				trkArray->AddAt(track2,1);
				AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
				if(trkArray) delete trkArray;

				if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
				if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;

				Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
				fDCAHe4Pi = *dca;
				if(dca) delete dca;
	        	//PA Array
				dd[0]=dn[0]-secVertex->GetX();
				dd[1]=dn[1]-secVertex->GetY();
				dd[2]=dn[2]-secVertex->GetZ();
				h->SetXYZ(-dd[0],-dd[1],-dd[2]);
				TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
				if(secVertex) delete secVertex;
				//++ set momentum ++
				vecFirstDaughter->SetXYZM(2.*track1->Px(),2.*track1->Py(),2.*track1->Pz(),AliPID::ParticleMass(AliPID::kAlpha));
				vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kPion));
				*vecMother = *vecFirstDaughter + *vecSecDaughter;
				//++ 4he info ++
				fhe4P = vecFirstDaughter->P();
				fhe4Ncls = track1->GetTPCNcls();
				fhe4NclsITS = track1->GetNumberOfITSClusters();
				fhe4Dedx = track1->GetTPCsignal();
				fEtaHe4 = track1->Eta();
				fPhiHe4 = track1->Phi();
				fGeoLengthHe4 = GeoLength(*track1);
				fhe4Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
				fhe4DcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
	        	//++ TofSignal He4 ++
				length = track1->GetIntegratedLength();
		        time0 = fPID->GetTOFResponse().GetStartTime(track1->P());//fESDpid->GetTOFResponse().GetTimeZero();
		        time = track1->GetTOFsignal() - time0;
		        if (time > 0) {
		        	beta = length / (2.99792457999999984e-02 * time);
		        	gamma = 1/TMath::Sqrt(1 - beta*beta);
		            mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
		            fTOFSignalHe4 = mass;
		        }
		        //++ pi info ++
		        fpiP = vecSecDaughter->P();
		        fpiNcls = track2->GetTPCNcls();
		        fpiNclsITS = track2->GetNumberOfITSClusters();
		        fpiDedx = track2->GetTPCsignal();
		        fEtaPi = track2->Eta();
		        fPhiPi = track2->Phi();
		        fGeoLengthPi = GeoLength(*track2);
		        fpiDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
		        fpiDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
	        	//++ TofSignal Pi ++
		        length = track2->GetIntegratedLength();
		        time0 = fPID->GetTOFResponse().GetStartTime(track2->P());//fESDpid->GetTOFResponse().GetTimeZero();
		        time = track2->GetTOFsignal() - time0;
		        if (time > 0) {
		        	beta = length / (2.99792457999999984e-02 * time);
		        	gamma = 1/TMath::Sqrt(1 - beta*beta);
		            mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
		            fTOFSignalPi = mass;
		        }
		        //++ dEdx ++
		        if (fBetheSplines) {
		        	fhe4DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kAlpha);
		        	fpiDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
		        } else {
		        	fhe4DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
		        	fpiDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
		        }
		        //++ mother info ++
		        fz4LH = 1;
		        fm4LH = vecMother->M();
		        fp4LH = vecMother->P();
		        fpt4LH = vecMother->Pt();
		        fy4LH = vecMother->Rapidity();
		        fPA4LH = TMath::Cos(vecMother->Angle(*h));
		        fct4LH = h->Mag() * fm4LH / fp4LH;
		        //++ armenteros podolanski ++
		        TVector3 vecP = vecFirstDaughter->Vect();
		        TVector3 vecN = vecSecDaughter->Vect();
		        TVector3 vecM = vecMother->Vect();

		        fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
		        fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

		        farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
		        farmpt = vecP.Mag()*sin(fthetaP);
		        //++ fill tree a ++
		        aTree->Fill();
		    }
        	//************ 4LH Neg **************
		    if(He4Neg1 && piPos2){
	    		//cout<<"4LH -"<<endl;
		    	TObjArray *trkArray = new TObjArray(2);                  
		    	trkArray->AddAt(track1,0);
		    	trkArray->AddAt(track2,1);
		    	AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
		    	if(trkArray) delete trkArray;

		    	if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
		    	if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;

		    	Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
		    	fDCAHe4Pi = *dca;
		    	if(dca) delete dca;
	        	//++ PA Array ++
		    	dd[0]=dn[0]-secVertex->GetX();
		    	dd[1]=dn[1]-secVertex->GetY();
		    	dd[2]=dn[2]-secVertex->GetZ();
		    	h->SetXYZ(-dd[0],-dd[1],-dd[2]);
		    	TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
		    	if(secVertex) delete secVertex;
		    	//++ set momentum ++
		    	vecFirstDaughter->SetXYZM(2.*track1->Px(),2.*track1->Py(),2.*track1->Pz(),AliPID::ParticleMass(AliPID::kAlpha));
		    	vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kPion));
		    	*vecMother = *vecFirstDaughter + *vecSecDaughter;
		    	//++ he4 info ++
		    	fhe4P = vecFirstDaughter->P();
		    	fhe4Ncls = track1->GetTPCNcls();
		    	fhe4NclsITS = track1->GetNumberOfITSClusters();
		    	fhe4Dedx = track1->GetTPCsignal();
		    	fEtaHe4 = track1->Eta();
		    	fPhiHe4 = track1->Phi();
		    	fGeoLengthHe4 = GeoLength(*track1);
		    	fhe4Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
		    	fhe4DcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
	        	//++ TofSignal He4 ++
		    	length = track1->GetIntegratedLength();
		        time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
		        time = track1->GetTOFsignal() - time0;
		        if (time > 0) {
		        	beta = length / (2.99792457999999984e-02 * time);
		        	gamma = 1/TMath::Sqrt(1 - beta*beta);
		            mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
		            fTOFSignalHe4 = mass;
		        }
		        //++ pi info ++
		        fpiP = vecSecDaughter->P();
		        fpiNcls = track2->GetTPCNcls();
		        fpiNclsITS = track2->GetNumberOfITSClusters();
		        fpiDedx = track2->GetTPCsignal();
		        fEtaPi = track2->Eta();
		        fPhiPi = track2->Phi();
		        fGeoLengthPi = GeoLength(*track2);
		        fpiDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
		        fpiDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
		        //++ TofSignal Pi ++
		        length = track2->GetIntegratedLength();
		        time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
		        time = track2->GetTOFsignal() - time0;
		        if (time > 0) {
		        	beta = length / (2.99792457999999984e-02 * time);
		        	gamma = 1/TMath::Sqrt(1 - beta*beta);
		            mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
		            fTOFSignalPi = mass;
		        }
		        //++ dEdx ++
		        if (fBetheSplines) {
		        	fhe4DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kAlpha);
		        	fpiDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
		        } else {
		        	fhe4DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
		        	fpiDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kPion);
		        }
		        //++ mother info ++
		        fz4LH = -1;
		        fm4LH = vecMother->M();
		        fp4LH = vecMother->P();
		        fpt4LH = vecMother->Pt();
		        fy4LH = vecMother->Rapidity();
		        fPA4LH = TMath::Cos(vecMother->Angle(*h));
		        fct4LH = h->Mag() * fm4LH / fp4LH;
		        //++ armenteros podolanski ++
		        TVector3 vecP = vecSecDaughter->Vect();
		        TVector3 vecN = vecFirstDaughter->Vect();
		        TVector3 vecM = vecMother->Vect();

		        fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
		        fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

		        farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
		        farmpt = vecP.Mag()*sin(fthetaP);
		        //++ fill tree a ++
		        aTree->Fill();
		    }
        	//************ 4Li Pos **************
		    if(He3Pos1 && pPos2){
	    		//cout<<"4Li"<<endl;
	    		Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
		    	fDCAHe3P = *dca;
		    	if(dca) delete dca;
	    		//++ set momentum ++
		    	vecFirstDaughter->SetXYZM(2*track1->Px(), 2*track1->Py(), 2*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
		    	vecSecDaughter->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
		    	*vecMother = *vecFirstDaughter + *vecSecDaughter;
	      		//++ Daughter P ++
		    	fhe3P = vecFirstDaughter->P();
		    	fpP = vecSecDaughter->P();
	      		//++ Mother m, p, pt ++
		    	fzLi4 = 3;
		    	fmLi4 = vecMother->M();
		    	fpLi4 = vecMother->P();
		    	fptLi4 = vecMother->Pt();
		    	fyLi4 = vecMother->Rapidity();
	      		//++ TofSignal He3 ++
			    length = track1->GetIntegratedLength();
		        time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
		        time = track1->GetTOFsignal() - time0;
		        if (time > 0) {
		        	beta = length / (2.99792457999999984e-02 * time);
		        	gamma = 1/TMath::Sqrt(1 - beta*beta);
		            mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); 
		            fTOFSignalHe3 = mass;
		        }
		        //++ TofSignal p ++
		        length = track2->GetIntegratedLength();
		        time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
		        time = track2->GetTOFsignal() - time0;
		        if (time > 0) {
		        	beta = length / (2.99792457999999984e-02 * time);
		        	gamma = 1/TMath::Sqrt(1 - beta*beta);
		            mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
		            fTOFSignalP = mass;
		        }
		      	//++ dEdx Sigma ++
		        if (fBetheSplines) {
		        	fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
		        	fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
		        } else {
		        	fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
		        	fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
		        }
		      	//++ N Clusters TPC ++
		        fhe3Ncls = track1->GetTPCNcls();
		        fpNcls = track2->GetTPCNcls();
		      	//++ NClusters ITS ++
		        fhe3NclsITS = track1->GetNumberOfITSClusters();
		        fpNclsITS = track2->GetNumberOfITSClusters();
		        //++ fill tree b ++
		        bTree->Fill();
		    }
        	//************ 4Li Neg **************
		    if(He3Neg1 && pNeg2){
	    		//cout<<"4Li -"<<endl;
	    		Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
		    	fDCAHe3P = *dca;
		    	if(dca) delete dca;
	    		//++ set momentum ++
		    	vecFirstDaughter->SetXYZM(2*track1->Px(), 2*track1->Py(), 2*track1->Pz(), AliPID::ParticleMass(AliPID::kHe3));
		    	vecSecDaughter->SetXYZM(track2->Px(), track2->Py(), track2->Pz(), AliPID::ParticleMass(AliPID::kProton));
		    	*vecMother = *vecFirstDaughter + *vecSecDaughter;
          		//++ Daughter P ++
		    	fhe3P = vecFirstDaughter->P();
		    	fpP = vecSecDaughter->P();
          		//++ Mother m, p, pt ++
		    	fzLi4 = -3;
		    	fmLi4 = vecMother->M();
		    	fpLi4 = vecMother->P();
		    	fptLi4 = vecMother->Pt();
		    	fyLi4 = vecMother->Rapidity();
         		 //++ TofSignal He3 ++
		    	length = track1->GetIntegratedLength();
		        time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
		        time = track1->GetTOFsignal() - time0;
		        if (time > 0) {
		        	beta = length / (2.99792457999999984e-02 * time);
		        	gamma = 1/TMath::Sqrt(1 - beta*beta);
		            mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
		            fTOFSignalHe3 = mass;
		        }
		        //++ TofSignal p ++
		        length = track2->GetIntegratedLength();
		        time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
		        time = track2->GetTOFsignal() - time0;
		        if (time > 0) {
		        	beta = length / (2.99792457999999984e-02 * time);
		        	gamma = 1/TMath::Sqrt(1 - beta*beta);
		            mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
		            fTOFSignalP = mass;
		        }
	          	//++ dEdx Sigma ++
		        if (fBetheSplines) {
		        	fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
		        	fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
		        } else {
		        	fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
		        	fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
		        }
	          	//++ N Clusters TPC ++
		        fhe3Ncls = track1->GetTPCNcls();
		        fpNcls = track2->GetTPCNcls();
	          	//++ NClusters ITS ++
		        fhe3NclsITS = track1->GetNumberOfITSClusters();
		        fpNclsITS = track2->GetNumberOfITSClusters();
		        //++ fill tree b ++
		        bTree->Fill();
		    }
        	//++++++++++++++++++++++++++++++++ third track +++++++++++++++++++++++++++++++++++++++//
		    for (Int_t kTracks = jTracks + 1; kTracks < fESDevent->GetNumberOfTracks(); kTracks++) {

				AliESDtrack* track3 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(kTracks));        

				if (!track3->GetInnerParam()) continue;

				if (!trackCutsPi.AcceptTrack(track3)) continue;

				if (jTracks == kTracks || iTracks == kTracks) continue;

				Double_t ptot3 = track3->GetInnerParam()->GetP();
				Double_t sign3 = track3->GetSign();

				if(sign2 == sign3) continue;

				fHistdEdx->Fill(ptot3*sign3, track3->GetTPCsignal());
				//++ reset PID Bools ++
				piPos3 = kFALSE;
				piNeg3 = kFALSE;
				//++ check PID ++
				if(sign3 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) <= 3) piPos3 = kTRUE;
				else if(sign3 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track3, AliPID::kPion)) <= 3) piNeg3 = kTRUE;
				else continue;
				//++ check expected combinations ++
				if(!(He3Pos1 && pPos2 && piNeg3) && !(He3Neg1 && pNeg2 && piPos3) && !(He3Pos1 && dPos2 && piNeg3) && !(He3Neg1 && dNeg2 && piPos3) && !(He4Pos1 && pPos2 && piNeg3) && !(He4Neg1 && pNeg2 && piPos3) && !(tPos1 && pPos2 && piNeg3) && !(tNeg1 && pNeg2 && piPos3)) continue;
				//++ pre DCA cut ++
				Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
				if(*dca > 1.0) continue;
				Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
				if(*dca2 > 1.0) continue;
				if(dca) delete dca;
				if(dca2) delete dca2;
				//************ 4LH Pos (t + p + pi) **************
				if(tPos1 && pPos2 && piNeg3 && Check4LH3B){
					//cout<<"4LH3B"<<endl;
					TObjArray *trkArray = new TObjArray(3);                  
					trkArray->AddAt(track1,0);
					trkArray->AddAt(track2,1);
					trkArray->AddAt(track3,2);
					AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
					if(trkArray) delete trkArray;

					if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track3->PropagateToDCA(secVertex, fMagneticField, 10))continue;

					Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCATPi = *dca;
					Double_t *dca1 = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCATP = *dca1;
					Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
					fDCAPPi = *dca2;
					if(dca) delete dca;
					if(dca1) delete dca1;
					if(dca2) delete dca2;
					//++ PA Array ++
					dd[0]=dn[0]-secVertex->GetX();
					dd[1]=dn[1]-secVertex->GetY();
					dd[2]=dn[2]-secVertex->GetZ();
					h->SetXYZ(-dd[0],-dd[1],-dd[2]);
					TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
					if(secVertex) delete secVertex;
					//++ set momentum ++
					vecFirstDaughter->SetXYZM(track1->Px(),track1->Py(),track1->Pz(),AliPID::ParticleMass(AliPID::kTriton));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ t info ++
					ftP = vecFirstDaughter->P();
					ftNcls = track1->GetTPCNcls();
					ftNclsITS = track1->GetNumberOfITSClusters();
					ftDedx = track1->GetTPCsignal();
					fEtaT = track1->Eta();
					fPhiT = track1->Phi();
					fGeoLengthT = GeoLength(*track1);
					ftDca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					ftDcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal t ++
					length = track1->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
					time = track1->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalT = mass;
					}
					//++ p info ++
					fpP = vecSecDaughter->P();
					fpNcls = track2->GetTPCNcls();
					fpNclsITS = track2->GetNumberOfITSClusters();
					fpDedx = track2->GetTPCsignal();
					fEtaP = track2->Eta();
					fPhiP = track2->Phi();
					fGeoLengthP = GeoLength(*track2);
					fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal p ++
					length = track2->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
					time = track2->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalP = mass;
					}
					//++ pi info ++
					fpiP = vecThirdDaughter->P();
					fpiNcls = track3->GetTPCNcls();
					fpiNclsITS = track3->GetNumberOfITSClusters();
					fpiDedx = track3->GetTPCsignal();
					fEtaPi = track3->Eta();
					fPhiPi = track3->Phi();
					fGeoLengthPi = GeoLength(*track3);
					fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpiDcaSec = TMath::Abs(track3->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal Pi ++
					length = track3->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
					time = track3->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalPi = mass;
					}
					//++ dEdx ++
					if (fBetheSplines) {
						ftDedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kTriton);
						fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					} else {
						ftDedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 2, fBetheParamsT);
						fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					}
					//++ mother info ++
					fz4LH = 1;
					fm4LH = vecMother->M();
					fp4LH = vecMother->P();
					fpt4LH = vecMother->Pt();
					fy4LH = vecMother->Rapidity();
					fPA4LH = TMath::Cos(vecMother->Angle(*h));
					fct4LH = h->Mag() * fm4LH / fp4LH;
					//++ fill tree f ++
					fTree->Fill();
				}
				//************ 4LH Neg (t + p + pi) **************
				if(tNeg1 && pNeg2 && piPos3 && Check4LH3B){
					//cout<<"4LH3B -"<<endl;
					TObjArray *trkArray = new TObjArray(3);                  
					trkArray->AddAt(track1,0);
					trkArray->AddAt(track2,1);
					trkArray->AddAt(track3,2);
					AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
					if(trkArray) delete trkArray;

					if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track3->PropagateToDCA(secVertex, fMagneticField, 10))continue;

					Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCATPi = *dca;
					Double_t *dca1 = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCATP = *dca1;
					Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
					fDCAPPi = *dca2;
					if(dca) delete dca;
					if(dca1) delete dca1;
					if(dca2) delete dca2;
					//++ PA Array ++
					dd[0]=dn[0]-secVertex->GetX();
					dd[1]=dn[1]-secVertex->GetY();
					dd[2]=dn[2]-secVertex->GetZ();
					h->SetXYZ(-dd[0],-dd[1],-dd[2]);
					TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
					if(secVertex) delete secVertex;
					//++ set momentum ++
					vecFirstDaughter->SetXYZM(track1->Px(),track1->Py(),track1->Pz(),AliPID::ParticleMass(AliPID::kTriton));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ t info ++
					ftP = vecFirstDaughter->P();
					ftNcls = track1->GetTPCNcls();
					ftNclsITS = track1->GetNumberOfITSClusters();
					ftDedx = track1->GetTPCsignal();
					fEtaT = track1->Eta();
					fPhiT = track1->Phi();
					fGeoLengthT = GeoLength(*track1);
					ftDca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					ftDcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal t ++
					length = track1->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
					time = track1->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalT = mass;
					}
					//++ p info ++
					fpP = vecSecDaughter->P();
					fpNcls = track2->GetTPCNcls();
					fpNclsITS = track2->GetNumberOfITSClusters();
					fpDedx = track2->GetTPCsignal();
					fEtaP = track2->Eta();
					fPhiP = track2->Phi();
					fGeoLengthP = GeoLength(*track2);
					fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal p ++
					length = track2->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
					time = track2->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalP = mass;
					}
					//++ pi info ++
					fpiP = vecThirdDaughter->P();
					fpiNcls = track3->GetTPCNcls();
					fpiNclsITS = track3->GetNumberOfITSClusters();
					fpiDedx = track3->GetTPCsignal();
					fEtaPi = track3->Eta();
					fPhiPi = track3->Phi();
					fGeoLengthPi = GeoLength(*track3);
					fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpiDcaSec = TMath::Abs(track3->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal Pi ++
					length = track3->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
					time = track3->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalPi = mass;
					}
					//++ dEdx ++
					if (fBetheSplines) {
						ftDedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kTriton);
						fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					} else {
						ftDedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kTriton), 2, fBetheParamsT);
						fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					}
					//++ mother info ++
					fz4LH = -1;
					fm4LH = vecMother->M();
					fp4LH = vecMother->P();
					fpt4LH = vecMother->Pt();
					fy4LH = vecMother->Rapidity();
					fPA4LH = TMath::Cos(vecMother->Angle(*h));
					fct4LH = h->Mag() * fm4LH / fp4LH;
					//++ fill tree f ++
					fTree->Fill();
				}
				//************ 4LHe Pos (3He + p + pi) **************
				if(He3Pos1 && pPos2 && piNeg3){
					//cout<<"4LHe"<<endl;
					TObjArray *trkArray = new TObjArray(3);                  
					trkArray->AddAt(track1,0);
					trkArray->AddAt(track2,1);
					trkArray->AddAt(track3,2);
					AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
					if(trkArray) delete trkArray;

					if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track3->PropagateToDCA(secVertex, fMagneticField, 10))continue;

					Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe3Pi = *dca;
					Double_t *dca1 = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe3P = *dca1;
					Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
					fDCAPPi = *dca2;
					if(dca) delete dca;
					if(dca1) delete dca1;
					if(dca2) delete dca2;
					//++ PA Array ++
					dd[0]=dn[0]-secVertex->GetX();
					dd[1]=dn[1]-secVertex->GetY();
					dd[2]=dn[2]-secVertex->GetZ();
					h->SetXYZ(-dd[0],-dd[1],-dd[2]);
					TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
					if(secVertex) delete secVertex;
					//++ set momentum ++
					vecFirstDaughter->SetXYZM(2.*track1->Px(),2.*track1->Py(),2.*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ he3 info ++
					fhe3P = vecFirstDaughter->P();
					fhe3Ncls = track1->GetTPCNcls();
					fhe3NclsITS = track1->GetNumberOfITSClusters();
					fhe3Dedx = track1->GetTPCsignal();
					fEtaHe3 = track1->Eta();
					fPhiHe3 = track1->Phi();
					fGeoLengthHe3 = GeoLength(*track1);
					fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fhe3DcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal He3 ++
					length = track1->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
					time = track1->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalHe3 = mass;
					}
					//++ p info ++
					fpP = vecSecDaughter->P();
					fpNcls = track2->GetTPCNcls();
					fpNclsITS = track2->GetNumberOfITSClusters();
					fpDedx = track2->GetTPCsignal();
					fEtaP = track2->Eta();
					fPhiP = track2->Phi();
					fGeoLengthP = GeoLength(*track2);
					fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal p ++
					length = track2->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
					time = track2->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalP = mass;
					}
					//++ pi info ++
					fpiP = vecThirdDaughter->P();
					fpiNcls = track3->GetTPCNcls();
					fpiNclsITS = track3->GetNumberOfITSClusters();
					fpiDedx = track3->GetTPCsignal();
					fEtaPi = track3->Eta();
					fPhiPi = track3->Phi();
					fGeoLengthPi = GeoLength(*track3);
					fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpiDcaSec = TMath::Abs(track3->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal Pi ++
					length = track3->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
					time = track3->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalPi = mass;
					}
					//++ dEdx ++
					if (fBetheSplines) {
						fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
						fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					} else {
						fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
						fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					}
					//++ mother info ++
					fz4LHe = 2;
					fm4LHe = vecMother->M();
					fp4LHe = vecMother->P();
					fpt4LHe = vecMother->Pt();
					fy4LHe = vecMother->Rapidity();
					fPA4LHe = TMath::Cos(vecMother->Angle(*h));
					fct4LHe = h->Mag() * fm4LHe / fp4LHe;
					//++ fill tree e ++
					eTree->Fill();
				}
				//************ 4LHe Neg (3He + p + pi) **************
				if(He3Neg1 && pNeg2 && piPos3){
					//cout<<"4LHe -"<<endl;
					TObjArray *trkArray = new TObjArray(3);                  
					trkArray->AddAt(track1,0);
					trkArray->AddAt(track2,1);
					trkArray->AddAt(track3,2);
					AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
					if(trkArray) delete trkArray;

					if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track3->PropagateToDCA(secVertex, fMagneticField, 10))continue;

					Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe3Pi = *dca;
					Double_t *dca1 = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe3P = *dca1;
					Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
					fDCAPPi = *dca2;
					if(dca) delete dca;
					if(dca1) delete dca1;
					if(dca2) delete dca2;
					//++ PA Array ++
					dd[0]=dn[0]-secVertex->GetX();
					dd[1]=dn[1]-secVertex->GetY();
					dd[2]=dn[2]-secVertex->GetZ();
					h->SetXYZ(-dd[0],-dd[1],-dd[2]);
					TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
					if(secVertex) delete secVertex;
					//++ set momentum ++
					vecFirstDaughter->SetXYZM(2.*track1->Px(),2.*track1->Py(),2.*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ he3 info ++
					fhe3P = vecFirstDaughter->P();
					fhe3Ncls = track1->GetTPCNcls();
					fhe3NclsITS = track1->GetNumberOfITSClusters();
					fhe3Dedx = track1->GetTPCsignal();
					fEtaHe3 = track1->Eta();
					fPhiHe3 = track1->Phi();
					fGeoLengthHe3 = GeoLength(*track1);
					fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fhe3DcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal He3 ++
					length = track1->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
					time = track1->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalHe3 = mass;
					}
					//++ p info ++
					fpP = vecSecDaughter->P();
					fpNcls = track2->GetTPCNcls();
					fpNclsITS = track2->GetNumberOfITSClusters();
					fpDedx = track2->GetTPCsignal();
					fEtaP = track2->Eta();
					fPhiP = track2->Phi();
					fGeoLengthP = GeoLength(*track2);
					fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal p ++
					length = track2->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
					time = track2->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalP = mass;
					}
					//++ pi info ++
					fpiP = vecThirdDaughter->P();
					fpiNcls = track3->GetTPCNcls();
					fpiNclsITS = track3->GetNumberOfITSClusters();
					fpiDedx = track3->GetTPCsignal();
					fEtaPi = track3->Eta();
					fPhiPi = track3->Phi();
					fGeoLengthPi = GeoLength(*track3);
					fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpiDcaSec = TMath::Abs(track3->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal Pi ++
					length = track3->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
					time = track3->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalPi = mass;
					}
					//++ dEdx ++
					if (fBetheSplines) {
						fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
						fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					} else {
						fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
						fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					}	
					//++ mother info ++
					fz4LHe = -2;
					fm4LHe = vecMother->M();
					fp4LHe = vecMother->P();
					fpt4LHe = vecMother->Pt();
					fy4LHe = vecMother->Rapidity();
					fPA4LHe = TMath::Cos(vecMother->Angle(*h));
					fct4LHe = h->Mag() * fm4LHe / fp4LHe;
					//++ fill tree e ++
					eTree->Fill();
				}
				//************ 5LHe Pos (4He + p + pi) **************
				if(He4Pos1 && pPos2 && piNeg3){
					//cout<<"5LHe1"<<endl;
					TObjArray *trkArray = new TObjArray(3);                  
					trkArray->AddAt(track1,0);
					trkArray->AddAt(track2,1);
					trkArray->AddAt(track3,2);
					AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
					if(trkArray) delete trkArray;

					if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track3->PropagateToDCA(secVertex, fMagneticField, 10))continue;

					Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe4Pi = *dca;
					Double_t *dca1 = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe4P = *dca1;
					Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
					fDCAPPi = *dca2;
					if(dca) delete dca;
					if(dca1) delete dca1;
					if(dca2) delete dca2;
					//++ PA Array ++
					dd[0]=dn[0]-secVertex->GetX();
					dd[1]=dn[1]-secVertex->GetY();
					dd[2]=dn[2]-secVertex->GetZ();
					h->SetXYZ(-dd[0],-dd[1],-dd[2]);
					TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
					if(secVertex) delete secVertex;
					//++ set momentum ++
					vecFirstDaughter->SetXYZM(2.*track1->Px(),2.*track1->Py(),2.*track1->Pz(),AliPID::ParticleMass(AliPID::kAlpha));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ he4 info ++
					fhe4P = vecFirstDaughter->P();
					fhe4Ncls = track1->GetTPCNcls();
					fhe4NclsITS = track1->GetNumberOfITSClusters();
					fhe4Dedx = track1->GetTPCsignal();
					fEtaHe4 = track1->Eta();
					fPhiHe4 = track1->Phi();
					fGeoLengthHe4 = GeoLength(*track1);
					fhe4Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fhe4DcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal He4 ++
					length = track1->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
					time = track1->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalHe4 = mass;
					}
					//++ p info ++
					fpP = vecSecDaughter->P();
					fpNcls = track2->GetTPCNcls();
					fpNclsITS = track2->GetNumberOfITSClusters();
					fpDedx = track2->GetTPCsignal();
					fEtaP = track2->Eta();
					fPhiP = track2->Phi();
					fGeoLengthP = GeoLength(*track2);
					fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal p ++
					length = track2->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
					time = track2->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalP = mass;
					}
					//++ pi info ++
					fpiP = vecThirdDaughter->P();
					fpiNcls = track3->GetTPCNcls();
					fpiNclsITS = track3->GetNumberOfITSClusters();
					fpiDedx = track3->GetTPCsignal();
					fEtaPi = track3->Eta();
					fPhiPi = track3->Phi();
					fGeoLengthPi = GeoLength(*track3);
					fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpiDcaSec = TMath::Abs(track3->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal Pi ++
					length = track3->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
					time = track3->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalPi = mass;
					}
					//++ dEdx ++
					if (fBetheSplines) {
						fhe4DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kAlpha);
						fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					} else {
						fhe4DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
						fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					}
					//++ mother info ++
					fz5LHe = 2;
					fm5LHe = vecMother->M();
					fp5LHe = vecMother->P();
					fpt5LHe = vecMother->Pt();
					fy5LHe = vecMother->Rapidity();
					fPA5LHe = TMath::Cos(vecMother->Angle(*h));
					fct5LHe = h->Mag() * fm5LHe / fp5LHe;
					//++ fill tree c ++
					cTree->Fill();
				}
				//************ 5LHe Neg (4He + p + pi) **************
				if(He4Neg1 && pNeg2 && piPos3){
					//cout<<"5LHe1 -"<<endl;
					TObjArray *trkArray = new TObjArray(3);                  
					trkArray->AddAt(track1,0);
					trkArray->AddAt(track2,1);
					trkArray->AddAt(track3,2);
					AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
					if(trkArray) delete trkArray;

					if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track3->PropagateToDCA(secVertex, fMagneticField, 10))continue;

					Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe4Pi = *dca;
					Double_t *dca1 = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe4P = *dca1;
					Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
					fDCAPPi = *dca2;
					if(dca) delete dca;
					if(dca1) delete dca1;
					if(dca2) delete dca2;
					//++ PA Array ++
					dd[0]=dn[0]-secVertex->GetX();
					dd[1]=dn[1]-secVertex->GetY();
					dd[2]=dn[2]-secVertex->GetZ();
					h->SetXYZ(-dd[0],-dd[1],-dd[2]);
					TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
					if(secVertex) delete secVertex;
					//++ set momentum ++
					vecFirstDaughter->SetXYZM(2.*track1->Px(),2.*track1->Py(),2.*track1->Pz(),AliPID::ParticleMass(AliPID::kAlpha));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ he4 info ++
					fhe4P = vecFirstDaughter->P();
					fhe4Ncls = track1->GetTPCNcls();
					fhe4NclsITS = track1->GetNumberOfITSClusters();
					fhe4Dedx = track1->GetTPCsignal();
					fEtaHe4 = track1->Eta();
					fPhiHe4 = track1->Phi();
					fGeoLengthHe4 = GeoLength(*track1);
					fhe4Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fhe4DcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal He4 ++
					length = track1->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
					time = track1->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalHe4 = mass;
					}
					//++ p info ++
					fpP = vecSecDaughter->P();
					fpNcls = track2->GetTPCNcls();
					fpNclsITS = track2->GetNumberOfITSClusters();
					fpDedx = track2->GetTPCsignal();
					fEtaP = track2->Eta();
					fPhiP = track2->Phi();
					fGeoLengthP = GeoLength(*track2);
					fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal p ++
					length = track2->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
					time = track2->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalP = mass;
					}
					//++ pi info ++
					fpiP = vecThirdDaughter->P();
					fpiNcls = track3->GetTPCNcls();
					fpiNclsITS = track3->GetNumberOfITSClusters();
					fpiDedx = track3->GetTPCsignal();
					fEtaPi = track3->Eta();
					fPhiPi = track3->Phi();
					fGeoLengthPi = GeoLength(*track3);
					fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpiDcaSec = TMath::Abs(track3->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal Pi ++ 
					length = track3->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
					time = track3->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalPi = mass;
					}
					//++ dEdx ++
					if (fBetheSplines) {
						fhe4DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kAlpha);
						fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					} else {
						fhe4DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kAlpha), 2, fBetheParamsHe);
						fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					}
					//++ mother info ++
					fz5LHe = -2;
					fm5LHe = vecMother->M();
					fp5LHe = vecMother->P();
					fpt5LHe = vecMother->Pt();
					fy5LHe = vecMother->Rapidity();
					fPA5LHe = TMath::Cos(vecMother->Angle(*h));
					fct5LHe = h->Mag() * fm5LHe / fp5LHe;
					//++ fill tree c ++
					cTree->Fill();
				}
				//************ 5LHe Pos (3He + d + pi) **************
				if(He3Pos1 && dPos2 && piNeg3){
					//cout<<"5LHe2"<<endl;
					TObjArray *trkArray = new TObjArray(3);                  
					trkArray->AddAt(track1,0);
					trkArray->AddAt(track2,1);
					trkArray->AddAt(track3,2);
					AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
					if(trkArray) delete trkArray;

					if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track3->PropagateToDCA(secVertex, fMagneticField, 10))continue;

					Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe3Pi = *dca;
					Double_t *dca1 = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe3d = *dca1;
					Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
					fDCAdPi = *dca2;
					if(dca) delete dca;
					if(dca1) delete dca1;
					if(dca2) delete dca2;
					//++ PA Array ++
					dd[0]=dn[0]-secVertex->GetX();
					dd[1]=dn[1]-secVertex->GetY();
					dd[2]=dn[2]-secVertex->GetZ();
					h->SetXYZ(-dd[0],-dd[1],-dd[2]);
					TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
					if(secVertex) delete secVertex;
					//++ set momentum ++
					vecFirstDaughter->SetXYZM(2.*track1->Px(),2.*track1->Py(),2.*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kDeuteron));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ he3 info ++
					fhe3P = vecFirstDaughter->P();
					fhe3Ncls = track1->GetTPCNcls();
					fhe3NclsITS = track1->GetNumberOfITSClusters();
					fhe3Dedx = track1->GetTPCsignal();
					fEtaHe3 = track1->Eta();
					fPhiHe3 = track1->Phi();
					fGeoLengthHe3 = GeoLength(*track1);
					fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fhe3DcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal He3 ++
					length = track1->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
					time = track1->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalHe3 = mass;
					}
					//++ d info ++
					fdP = vecSecDaughter->P();
					fdNcls = track2->GetTPCNcls();
					fdNclsITS = track2->GetNumberOfITSClusters();
					fdDedx = track2->GetTPCsignal();
					fEtaD = track2->Eta();
					fPhiD = track2->Phi();
					fGeoLengthD = GeoLength(*track2);
					fdDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fdDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal d ++
					length = track2->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
					time = track2->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalD = mass;
					}
					//++ pi info ++
					fpiP = vecThirdDaughter->P();
					fpiNcls = track3->GetTPCNcls();
					fpiNclsITS = track3->GetNumberOfITSClusters();
					fpiDedx = track3->GetTPCsignal();
					fEtaPi = track3->Eta();
					fPhiPi = track3->Phi();
					fGeoLengthPi = GeoLength(*track3);
					fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpiDcaSec = TMath::Abs(track3->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal Pi ++
					length = track3->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
					time = track3->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalPi = mass;
					}
					//++ dEdx ++
					if (fBetheSplines) {
						fhe4DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
						fdDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kDeuteron);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					} else {
						fhe4DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
						fdDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					}
					//++ mother info ++
					fz5LHe = 2;
					fm5LHe = vecMother->M();
					fp5LHe = vecMother->P();
					fpt5LHe = vecMother->Pt();
					fy5LHe = vecMother->Rapidity();
					fPA5LHe = TMath::Cos(vecMother->Angle(*h));
					fct5LHe = h->Mag() * fm5LHe / fp5LHe;
					//++ fill tree d ++
					dTree->Fill();
				}
				//************ 5LHe Neg (3He + d + pi) **************
				if(He3Neg1 && dNeg2 && piPos3){
					//cout<<"5LHe2 - "<<endl;
					TObjArray *trkArray = new TObjArray(3);                  
					trkArray->AddAt(track1,0);
					trkArray->AddAt(track2,1);
					trkArray->AddAt(track3,2);
					AliESDVertex *secVertex = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
					if(trkArray) delete trkArray;

					if(!track1->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track2->PropagateToDCA(secVertex, fMagneticField, 10))continue;
					if(!track3->PropagateToDCA(secVertex, fMagneticField, 10))continue;

					Double_t *dca = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe3Pi = *dca;
					Double_t *dca1 = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
					fDCAHe3d = *dca1;
					Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
					fDCAdPi = *dca2;
					if(dca) delete dca;
					if(dca1) delete dca1;
					if(dca2) delete dca2;
					//++ PA Array ++
					dd[0]=dn[0]-secVertex->GetX();
					dd[1]=dn[1]-secVertex->GetY();
					dd[2]=dn[2]-secVertex->GetZ();
					h->SetXYZ(-dd[0],-dd[1],-dd[2]);
					TVector3 secVtx(secVertex->GetX(),secVertex->GetY(),secVertex->GetZ());
					if(secVertex) delete secVertex;
					//++ set momentum ++
					vecFirstDaughter->SetXYZM(2.*track1->Px(),2.*track1->Py(),2.*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kDeuteron));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ he3 info ++
					fhe3P = vecFirstDaughter->P();
					fhe3Ncls = track1->GetTPCNcls();
					fhe3NclsITS = track1->GetNumberOfITSClusters();
					fhe3Dedx = track1->GetTPCsignal();
					fEtaHe3 = track1->Eta();
					fPhiHe3 = track1->Phi();
					fGeoLengthHe3 = GeoLength(*track1);
					fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fhe3DcaSec = TMath::Abs(track1->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal He3 ++
					length = track1->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
					time = track1->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalHe3 = mass;
					}
					//++ d info ++
					fdP = vecSecDaughter->P();
					fdNcls = track2->GetTPCNcls();
					fdNclsITS = track2->GetNumberOfITSClusters();
					fdDedx = track2->GetTPCsignal();
					fEtaD = track2->Eta();
					fPhiD = track2->Phi();
					fGeoLengthD = GeoLength(*track2);
					fdDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fdDcaSec = TMath::Abs(track2->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal d ++
					length = track2->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
					time = track2->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalD = mass;
					}
					//++ pi info ++
					fpiP = vecThirdDaughter->P();
					fpiNcls = track3->GetTPCNcls();
					fpiNclsITS = track3->GetNumberOfITSClusters();
					fpiDedx = track3->GetTPCsignal();
					fEtaPi = track3->Eta();
					fPhiPi = track3->Phi();
					fGeoLengthPi = GeoLength(*track3);
					fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
					fpiDcaSec = TMath::Abs(track3->GetD(secVtx.X(), secVtx.Y(), fMagneticField));
					//++ TofSignal Pi ++
					length = track3->GetIntegratedLength();
					time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
					time = track3->GetTOFsignal() - time0;
					if (time > 0) {
						beta = length / (2.99792457999999984e-02 * time);
						gamma = 1/TMath::Sqrt(1 - beta*beta);
					mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
					fTOFSignalPi = mass;
					}
					//++ dEdx ++
					if (fBetheSplines) {
						fhe4DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
						fdDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kDeuteron);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					} else {
						fhe4DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
						fdDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kDeuteron), 1, fBetheParamsT);
						fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
					}
					//++ mother info ++
					fz5LHe = -2;
					fm5LHe = vecMother->M();
					fp5LHe = vecMother->P();
					fpt5LHe = vecMother->Pt();
					fy5LHe = vecMother->Rapidity();
					fPA5LHe = TMath::Cos(vecMother->Angle(*h));
					fct5LHe = h->Mag() * fm5LHe / fp5LHe;
					//++ fill tree d ++
					dTree->Fill();
				}
				//++++++++++++++++++++++++++++++++ fourth Track +++++++++++++++++++++++++++++++++++++++//
				for (Int_t lTracks = kTracks + 1; lTracks < fESDevent->GetNumberOfTracks(); lTracks++) {

					AliESDtrack* track4 = dynamic_cast<AliESDtrack*>(fESDevent->GetTrack(lTracks));

					if (!track4->GetInnerParam()) continue;

					if (!trackCutsPi.AcceptTrack(track4)) continue;

					if (jTracks == lTracks || iTracks == lTracks || kTracks == lTracks) continue; //reject using the same track twice 

					Double_t ptot4 = track4->GetInnerParam()->GetP();
					Double_t sign4 = track4->GetSign();

					if(sign3 != sign4) continue;
					      
					fHistdEdx->Fill(ptot4*sign4, track4->GetTPCsignal());
					//++ reset PID Bool ++
					piPos4 = kFALSE;
					piNeg4 = kFALSE;
					//++ PID check ++
					if(sign4 >0 && TMath::Abs(fPID->NumberOfSigmasTPC(track4, AliPID::kPion)) <= 3) piPos4 = kTRUE;
					else if(sign4 <0 && TMath::Abs(fPID->NumberOfSigmasTPC(track4, AliPID::kPion)) <= 3) piNeg4 = kTRUE;
					else continue;
					//++ check different combinations ++
					if(!(He3Pos1 && pPos2 && piNeg3 && piNeg4) && !(He3Neg1 && pNeg2 && piPos3 && piPos4)) continue;
					//++ set momentum for pre cuts ++
					vecFirstDaughter->SetXYZM(2*track1->Px(),2*track1->Py(),2*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
					vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
					vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
					vecFourthDaughter->SetXYZM(track4->Px(),track4->Py(),track4->Pz(),AliPID::ParticleMass(AliPID::kPion));
					TLorentzVector HHe4(0.,0.,0.,0.);
					HHe4=*vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
					//++ Calculate DCA of track 4 and 4LHe Mother ++
					Double_t dcaM4 = TMath::Abs(track4->GetD(HHe4.X(), HHe4.Y(), fMagneticField));
					fDCAHe3Pi1 = dcaM4;
					if(dcaM4 > 1.0) continue;
					//++ cut on opening angle between mother of first three tracks and track 4 ++
					Double_t angle = HHe4.Angle(vecFourthDaughter->Vect());
					if(angle > 0.75) continue;
					//++ cut on distance from primary vertex only for 4LLH ++
					if(TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField)) <0.02) continue;
					if(TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField)) <0.1) continue;
					if(TMath::Abs(track4->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField)) <0.1) continue;				
					//********** 4LLH Pos **********
					if(He3Pos1 && pPos2 && piNeg3 && piNeg4){
						//cout<<"4LLH"<<endl;
						TObjArray *trkArray = new TObjArray(3);

						trkArray->AddAt(track1,0);
						trkArray->AddAt(track2,1);
						trkArray->AddAt(track3,2);

						AliESDVertex *tertVertex1 = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
						if(trkArray) delete trkArray;
	
						if(!track1->PropagateToDCA(tertVertex1, fMagneticField, 10)) continue;
						if(!track2->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;
						if(!track3->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;

						Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
						fDCAHe3P = *dca;
						Double_t *dca1 = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
						fDCAHe3Pi = *dca1;
						Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
						fDCAPPi = *dca2;
						if(dca) delete dca;
						if(dca1) delete dca1;
						if(dca2) delete dca2;
						//++ PA Array ++
						dd[0]=dn[0]-tertVertex1->GetX();
						dd[1]=dn[1]-tertVertex1->GetY();
						dd[2]=dn[2]-tertVertex1->GetZ();
						h->SetXYZ(-dd[0],-dd[1],-dd[2]);
						TVector3 tertVtx(tertVertex1->GetX(),tertVertex1->GetY(),tertVertex1->GetZ());
						if(tertVertex1) delete tertVertex1;
						//++ set momentum ++
						vecFirstDaughter->SetXYZM(2*track1->Px(),2*track1->Py(),2*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
						vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
						vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
						vecFourthDaughter->SetXYZM(track4->Px(),track4->Py(),track4->Pz(),AliPID::ParticleMass(AliPID::kPion));
						*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter + *vecThirdDaughter;    
						//++ Daughter momenta ++
						fpiP = vecThirdDaughter->Pt();          
						fhe3P =vecFirstDaughter->Pt();
						fpi1P = vecFourthDaughter->Pt();
						fpP = vecSecDaughter->Pt();
						//++ Ncls TPC ++
						fhe3Ncls = track1->GetTPCNcls();
						fpNcls = track2->GetTPCNcls();
						fpiNcls = track3->GetTPCNcls();
						fpi1Ncls = track4->GetTPCNcls();
						//++ Ncls ITS ++
						fhe3NclsITS = track1->GetNumberOfITSClusters();
						fpNclsITS = track2->GetNumberOfITSClusters();
						fpiNclsITS = track3->GetNumberOfITSClusters();
						fpi1NclsITS = track4->GetNumberOfITSClusters();
						//++ dEdx Tracks ++
						fhe3Dedx = track1->GetTPCsignal();
						fpDedx = track2->GetTPCsignal();
						fpiDedx = track3->GetTPCsignal();
						fpi1Dedx = track4->GetTPCsignal();
						//++ Eta ++
						fEtaHe3 = track1->Eta();
						fEtaP = track2->Eta();
						fEtaPi = track3->Eta();
						fEtaPi1 = track4->Eta();
						//++ Phi ++
						fPhiHe3 = track1->Phi();
						fPhiP = track2->Phi();
						fPhiPi = track3->Phi();
						fPhiPi1 = track4->Phi();
						//++ GeoLength ++
						fGeoLengthHe3 = GeoLength(*track1);
						fGeoLengthP = GeoLength(*track2);
						fGeoLengthPi = GeoLength(*track3);
						fGeoLengthPi1 = GeoLength(*track4);
						//++ TOFSignal He3 ++
						length = length = track1->GetIntegratedLength();
						time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
						time = track1->GetTOFsignal() - time0;
						if (time > 0) {
							beta = length / (2.99792457999999984e-02 * time);
							gamma = 1/TMath::Sqrt(1 - beta*beta);
						mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
						fTOFSignalHe3 = mass;
						}
						//++ TOFSignal P ++
						length = length = track2->GetIntegratedLength();
						time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
						time = track2->GetTOFsignal() - time0;
						if (time > 0) {
							beta = length / (2.99792457999999984e-02 * time);
							gamma = 1/TMath::Sqrt(1 - beta*beta);
						mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); 
						fTOFSignalP = mass;
						}
						//++ TOFSignal Pi ++
						length = length = track3->GetIntegratedLength();
						time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
						time = track3->GetTOFsignal() - time0;
						if (time > 0) {
							beta = length / (2.99792457999999984e-02 * time);
							gamma = 1/TMath::Sqrt(1 - beta*beta);
						mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1); 
						fTOFSignalPi = mass;
						}
						//++ TOFSignal Pi1 ++
						length = length = track4->GetIntegratedLength();
						time0 = fPID->GetTOFResponse().GetStartTime(track4->P());
						time = track4->GetTOFsignal() - time0;
						if (time > 0) {
							beta = length / (2.99792457999999984e-02 * time);
							gamma = 1/TMath::Sqrt(1 - beta*beta);
						mass = (track4->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
						fTOFSignalPi1 = mass;
						}
						//++ dEdx Sigma ++
						if (fBetheSplines) {
							fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
							fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
							fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
							fpi1DedxSigma = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
						} else {
							fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
							fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
							fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
							fpi1DedxSigma = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
						}
						//++ DCA prim Vertex ++
						fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
						fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
						fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
						fpi1Dca = TMath::Abs(track4->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField)); 
						//++ DCA sec ++
						fhe3DcaSec = TMath::Abs(track1->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
						fpDcaSec = TMath::Abs(track2->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
						fpiDcaSec = TMath::Abs(track3->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
						fpi1DcaSec = TMath::Abs(track4->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
						//++ Mother M (daughter sum) ++
						fz4LLH = 1;
						fmDaughterSum = vecMother->M();
						fpSum = vecMother->P();
						fptSum = vecMother->Pt();
						fySum = vecMother->Rapidity();          
						//++ 4LHe Vect ++
						HHe4=*vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
						//++ Mother M 4LHe ++
						fm4LHe = HHe4.M();
						fp4LHe = HHe4.P();
						fpt4LHe = HHe4.Pt();
						fy4LHe = HHe4.Rapidity();
						//++ PA ++
						Double_t pointingAngle = HHe4.Angle(*h);
						fPA4LLH = TMath::Cos(pointingAngle);
						//++ Angle between 4LHe Vect and pi Vect ++
						fAngle4LLH = HHe4.Angle(vecFourthDaughter->Vect());
						//++ ct ++
						fctSum = h->Mag() * fmDaughterSum / fpSum;
						fct4LHe = h->Mag() * fm4LHe / fp4LHe;
						//++ Armenteros Podolanski ++
						TVector3 vecP = HHe4.Vect();
						TVector3 vecN = vecFourthDaughter->Vect();
						TVector3 vecM = vecMother->Vect();

						fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
						fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

						farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
						farmpt = vecP.Mag()*sin(fthetaP);
						//++ fill tree g ++
						gTree->Fill();
					}//end DoubleHyperHydrogen4

					//********** 4LLH Neg**********
					if(He3Neg1 && pNeg2 && piPos3 && piPos4){
						//cout<<"4LLH -"<<endl;
						TObjArray *trkArray = new TObjArray(3);
						trkArray->AddAt(track1,0);
						trkArray->AddAt(track2,1);
						trkArray->AddAt(track3,2);

						AliESDVertex *tertVertex1 = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(trkArray);
						if(trkArray) delete trkArray;
						if(!track1->PropagateToDCA(tertVertex1, fMagneticField, 10)) continue;
						if(!track2->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;
						if(!track3->PropagateToDCA(tertVertex1, fMagneticField, 10))continue;

						Double_t *dca = new Double_t(track2->GetDCA(track1,fMagneticField,xthiss,xpp));
						fDCAHe3P = *dca;
						Double_t *dca1 = new Double_t(track3->GetDCA(track1,fMagneticField,xthiss,xpp));
						fDCAHe3Pi = *dca1;
						Double_t *dca2 = new Double_t(track3->GetDCA(track2,fMagneticField,xthiss,xpp));
						fDCAPPi = *dca2;
						if(dca) delete dca;
						if(dca1) delete dca1;
						if(dca2) delete dca2;
						//++ PA Array ++
						dd[0]=dn[0]-tertVertex1->GetX();
						dd[1]=dn[1]-tertVertex1->GetY();
						dd[2]=dn[2]-tertVertex1->GetZ();
						h->SetXYZ(-dd[0],-dd[1],-dd[2]);
						TVector3 tertVtx(tertVertex1->GetX(),tertVertex1->GetY(),tertVertex1->GetZ());
						if(tertVertex1) delete tertVertex1;
						//++ set momentum ++
						vecFirstDaughter->SetXYZM(2*track1->Px(),2*track1->Py(),2*track1->Pz(),AliPID::ParticleMass(AliPID::kHe3));
						vecThirdDaughter->SetXYZM(track3->Px(),track3->Py(),track3->Pz(),AliPID::ParticleMass(AliPID::kPion));
						vecSecDaughter->SetXYZM(track2->Px(),track2->Py(),track2->Pz(),AliPID::ParticleMass(AliPID::kProton));
						vecFourthDaughter->SetXYZM(track4->Px(),track4->Py(),track4->Pz(),AliPID::ParticleMass(AliPID::kPion));
						*vecMother = *vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter + *vecThirdDaughter;    
						//++ Daughter momenta ++
						fpiP = vecThirdDaughter->Pt();          
						fhe3P =vecFirstDaughter->Pt();
						fpi1P = vecFourthDaughter->Pt();
						fpP = vecSecDaughter->Pt();
						//++ Ncls TPC ++
						fhe3Ncls = track1->GetTPCNcls();
						fpNcls = track2->GetTPCNcls();
						fpiNcls = track3->GetTPCNcls();
						fpi1Ncls = track4->GetTPCNcls();
						//++ Ncls ITS ++
						fhe3NclsITS = track1->GetNumberOfITSClusters();
						fpNclsITS = track2->GetNumberOfITSClusters();
						fpiNclsITS = track3->GetNumberOfITSClusters();
						fpi1NclsITS = track4->GetNumberOfITSClusters();
						//++ dEdx Tracks ++
						fhe3Dedx = track1->GetTPCsignal();
						fpDedx = track2->GetTPCsignal();
						fpiDedx = track3->GetTPCsignal();
						fpi1Dedx = track4->GetTPCsignal();
						//++ Eta ++
						fEtaHe3 = track1->Eta();
						fEtaP = track2->Eta();
						fEtaPi = track3->Eta();
						fEtaPi1 = track4->Eta();
					 	//++ Phi ++
						fPhiHe3 = track1->Phi();
						fPhiP = track2->Phi();
						fPhiPi = track3->Phi();
						fPhiPi1 = track4->Phi();
						//++ GeoLength ++
						fGeoLengthHe3 = GeoLength(*track1);
						fGeoLengthP = GeoLength(*track2);
						fGeoLengthPi = GeoLength(*track3);
						fGeoLengthPi1 = GeoLength(*track4);
						//++ TOFSignal He3 ++
						length = length = track1->GetIntegratedLength();
						time0 = fPID->GetTOFResponse().GetStartTime(track1->P());
						time = track1->GetTOFsignal() - time0;
						if (time > 0) {
							beta = length / (2.99792457999999984e-02 * time);
							gamma = 1/TMath::Sqrt(1 - beta*beta);
						mass = (track1->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
						fTOFSignalHe3 = mass;
						}
						//++ TOFSignal P ++
						length = length = track2->GetIntegratedLength();
						time0 = fPID->GetTOFResponse().GetStartTime(track2->P());
						time = track2->GetTOFsignal() - time0;
						if (time > 0) {
							beta = length / (2.99792457999999984e-02 * time);
							gamma = 1/TMath::Sqrt(1 - beta*beta);
						mass = (track2->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
						fTOFSignalP = mass;
						}
						//++ TOFSignal Pi ++
						length = length = track3->GetIntegratedLength();
						time0 = fPID->GetTOFResponse().GetStartTime(track3->P());
						time = track3->GetTOFsignal() - time0;
						if (time > 0) {
							beta = length / (2.99792457999999984e-02 * time);
							gamma = 1/TMath::Sqrt(1 - beta*beta);
						mass = (track3->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
						fTOFSignalPi = mass;
						}
						//++ TOFSignal Pi1 ++
						length = length = track4->GetIntegratedLength();
						time0 = fPID->GetTOFResponse().GetStartTime(track4->P());
						time = track4->GetTOFsignal() - time0;
						if (time > 0) {
							beta = length / (2.99792457999999984e-02 * time);
							gamma = 1/TMath::Sqrt(1 - beta*beta);
						mass = (track4->GetInnerParam()->GetP())/TMath::Sqrt(gamma*gamma - 1);
						fTOFSignalPi1 = mass;
						}
						//++ dEdx Sigma ++
						if (fBetheSplines) {
							fhe3DedxSigma = fPID->NumberOfSigmasTPC(track1, AliPID::kHe3);
							fpDedxSigma = fPID->NumberOfSigmasTPC(track2, AliPID::kProton);
							fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
							fpi1DedxSigma = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
						} else {
							fhe3DedxSigma = Bethe(*track1, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
							fpDedxSigma = Bethe(*track2, AliPID::ParticleMass(AliPID::kProton), 1, fBetheParamsT);
							fpiDedxSigma = fPID->NumberOfSigmasTPC(track3, AliPID::kPion);
							fpi1DedxSigma = fPID->NumberOfSigmasTPC(track4, AliPID::kPion);
						}
						//++ DCA prim Vertex ++
						fhe3Dca = TMath::Abs(track1->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
						fpDca = TMath::Abs(track2->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
						fpiDca = TMath::Abs(track3->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
						fpi1Dca = TMath::Abs(track4->GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField)); 
						//++ DCA sec ++
						fhe3DcaSec = TMath::Abs(track1->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
						fpDcaSec = TMath::Abs(track2->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
						fpiDcaSec = TMath::Abs(track3->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
						fpi1DcaSec = TMath::Abs(track4->GetD(tertVtx.X(), tertVtx.Y(), fMagneticField));
						//++ Mother M (daughter sum) ++
						fz4LLH = -1;
						fmDaughterSum = vecMother->M();
						fpSum = vecMother->P();
						fptSum = vecMother->Pt();
						fySum = vecMother->Rapidity();          
						//++ 4LHe Vect ++
						HHe4=*vecFirstDaughter + *vecSecDaughter + *vecThirdDaughter;
						//++ Mother M 4LHe ++
						fm4LHe = HHe4.M();
						fp4LHe = HHe4.P();
						fpt4LHe = HHe4.Pt();
						fy4LHe = HHe4.Rapidity();
						//++ PA ++
						Double_t pointingAngle = HHe4.Angle(*h);
						fPA4LLH = TMath::Cos(pointingAngle);
						fAngle4LLH = HHe4.Angle(vecFourthDaughter->Vect());
						//++ ct ++
						fctSum = h->Mag() * fmDaughterSum / fpSum;
						fct4LHe = h->Mag() * fm4LHe / fp4LHe;
						//++ Armenteros Podolanski ++
						TVector3 vecN = HHe4.Vect();
						TVector3 vecP = vecFourthDaughter->Vect();
						TVector3 vecM = vecMother->Vect();

						fthetaP = TMath::ACos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
						fthetaN = TMath::ACos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));

						farmalpha = ((vecP.Mag())*TMath::Cos(fthetaP)-(vecN.Mag())*TMath::Cos(fthetaN))/((vecP.Mag())*TMath::Cos(fthetaP)+(vecN.Mag())*TMath::Cos(fthetaN));
						farmpt = vecP.Mag()*sin(fthetaP);
						//++ fill tree g ++
						gTree->Fill();
					}//end Anti-DoubleHyperHydrogen4
				}//end loop track4
			}//end loop track3
		}//end loop track2
	}//end loop track1
}
//-----------------------------------------------------------------
void AliAnalysisTaskDoubleHypNucTree::Terminate(const Option_t*) {
	if (!GetOutputData(0)) return;
}
/// Set trigger information in reduced event
/// \return returns kTRUE is successful.
//-----------------------------------------------------------------
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
//-----------------------------------------------------------------
Double_t AliAnalysisTaskDoubleHypNucTree::Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params){
	Double_t expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass,params[0],params[1],params[2],params[3],params[4]);
	Double_t sigma = expected*params[5];
	if (TMath::IsNaN(expected)) return -999;
	return (track.GetTPCsignal() - expected) / sigma;
}
//-----------------------------------------------------------------
Double_t AliAnalysisTaskDoubleHypNucTree::GeoLength(const AliESDtrack& track) {
	Double_t deadZoneWidth = 3.0;
	Double_t lengthInActiveZone = track.GetLengthInActiveZone(1, deadZoneWidth, 220, track.GetESDEvent()->GetMagneticField(),0,0);
	return lengthInActiveZone;
}
