/****************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 *                                                                        *
 * Author: Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 * Anatask to compute flatenicity (arXiv:2204.13733)                      *
 **************************************************************************/

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TProfile.h"
#include "TVector3.h"
#include <AliAnalysisFilter.h>
#include <AliESDVertex.h>
#include <AliHeader.h>
#include <AliMultiplicity.h>
#include <Riostream.h>
#include <TBits.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>
using std::cout;
using std::endl;

#include "AliAnalysisTaskFlatenicityPiKp.h"

static const Int_t nCent = 9;
static const Int_t nEta = 4;

static Double_t centClass[nCent + 1] = {0.0,  1.0,  5.0,  10.0, 20.0,
	30.0, 40.0, 50.0, 70.0, 100.0};

static const Char_t* V0MClass[nCent] = {"0_1","1_5","5_10","10_20",
	"20_30","30_40","40_50","50_70","70_100"};

static const Char_t* etaClass[nEta] = {"02","24","46","68"};
static const Char_t* ParticleType[3] = {"Primaries","MaterialInt","WeakDecays"};

static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskFlatenicityPiKp) // classimp: necessary for root

AliAnalysisTaskFlatenicityPiKp::AliAnalysisTaskFlatenicityPiKp()
	: AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0Mindex(-1), fV0MBin("0_1"), fDetFlat("V0"), fIsMCclosure(kFALSE),
	fRemoveTrivialScaling(kFALSE), fnGen(-1), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), fdEdxCalibrated(kTRUE),
	fEtaCalibrationPos(0x0), fEtaCalibrationNeg(0x0), fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fPeriod("16l"),
	ftrackmult08(0), fv0mpercentile(0), fFlat(-1), fFlatMC(-1),
	fMultSelection(0x0), hPtPrimIn(0), hPtPrimOut(0), hPtSecOut(0), hPtOut(0),
	hFlatV0vsFlatTPC(0), hFlatenicityBefore(0), hFlatenicity(0),
	hFlatenicityMC(0), hFlatResponse(0), hFlatVsPt(0), hFlatVsPtMC(0),
	hActivityV0DataSect(0), hActivityV0McSect(0), hFlatVsNchMC(0),
	hFlatVsV0M(0), hCounter(0),
	hMCPtPionPos(0),hMCPtKaonPos(0),hMCPtProtonPos(0),
	hMCPtPionNeg(0),hMCPtKaonNeg(0),hMCPtProtonNeg(0),
	hTPCRecTracksPionPos(0), hTPCRecTracksKaonPos(0), hTPCRecTracksProtonPos(0),
	hTPCRecTracksPionNeg(0), hTPCRecTracksKaonNeg(0), hTPCRecTracksProtonNeg(0),
	hTOFRecTracksPionPos(0), hTOFRecTracksKaonPos(0), hTOFRecTracksProtonPos(0),
	hTOFRecTracksPionNeg(0), hTOFRecTracksKaonNeg(0), hTOFRecTracksProtonNeg(0),
	hrTPCRecTracksPion(0), hrTPCRecTracksKaon(0), hrTPCRecTracksProton(0),
	hPionTPCDCAxyNegData(0), hPionTPCDCAxyPosData(0), hProtonTPCDCAxyNegData(0), hProtonTPCDCAxyPosData(0),
	hPionTOFDCAxyNegData(0), hPionTOFDCAxyPosData(0), hProtonTOFDCAxyNegData(0), hProtonTOFDCAxyPosData(0),
	hMIPVsEta(0), pMIPVsEta(0), hPlateauVsEta(0), pPlateauVsEta(0)


{
	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hFlatVsPtV0M[i_c] = 0;

		for (Int_t i_eta = 0; i_eta < nEta; ++i_eta){

			hNsigmaPiPos[i_c][i_eta] = 0;
			hNsigmaKPos[i_c][i_eta] = 0;
			hNsigmaPPos[i_c][i_eta] = 0;
			hPtTPCEtaPos[i_c][i_eta] = 0;

			hNsigmaPiNeg[i_c][i_eta] = 0;
			hNsigmaKNeg[i_c][i_eta] = 0;
			hNsigmaPNeg[i_c][i_eta] = 0;
			hPtTPCEtaNeg[i_c][i_eta] = 0;

			hBetaPos[i_c][i_eta] = 0;
			hMomentumTOFEtaPos[i_c][i_eta] = 0;
			hPtTOFEtaPos[i_c][i_eta] = 0;

			hBetaNeg[i_c][i_eta] = 0;
			hMomentumTOFEtaNeg[i_c][i_eta] = 0;
			hPtTOFEtaNeg[i_c][i_eta] = 0;

			hdEdx[i_c][i_eta] = 0;
			hPtrTPC[i_c][i_eta] = 0;
			if (i_c==0){
				hPtVsP[i_eta] = 0;	
				hKaonContamination[i_eta] = 0;
				hKaondEdx[i_eta] = 0;
				hProtonContamination[i_eta] = 0;
				hProtondEdx[i_eta] = 0;
			}
		}
	}

	for(Int_t i = 0; i < 3; ++i){
		hPionTOFDCAxyNeg[i] = 0;
		hProtonTOFDCAxyNeg[i] = 0;
		hPionTOFDCAxyPos[i] = 0;
		hProtonTOFDCAxyPos[i] = 0;
		hPionTPCDCAxyNeg[i] = 0;
		hProtonTPCDCAxyNeg[i] = 0;	
		hPionTPCDCAxyPos[i] = 0;	
		hProtonTPCDCAxyPos[i] = 0;

	}

}
//_____________________________________________________________________________
AliAnalysisTaskFlatenicityPiKp::AliAnalysisTaskFlatenicityPiKp(const char *name)
	: AliAnalysisTaskSE(name), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0Mindex(-1), fV0MBin("0_1"), fDetFlat("V0"), fIsMCclosure(kFALSE),
	fRemoveTrivialScaling(kFALSE), fnGen(-1), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), fdEdxCalibrated(kTRUE),
	fEtaCalibrationPos(0x0), fEtaCalibrationNeg(0x0), fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fPeriod("16l"),
	ftrackmult08(0), fv0mpercentile(0), fFlat(-1), fFlatMC(-1),
	fMultSelection(0x0), hPtPrimIn(0), hPtPrimOut(0), hPtSecOut(0), hPtOut(0),
	hFlatV0vsFlatTPC(0), hFlatenicityBefore(0), hFlatenicity(0),
	hFlatenicityMC(0), hFlatResponse(0), hFlatVsPt(0), hFlatVsPtMC(0),
	hActivityV0DataSect(0), hActivityV0McSect(0), hFlatVsNchMC(0),
	hFlatVsV0M(0), hCounter(0), 
	hMCPtPionPos(0),hMCPtKaonPos(0),hMCPtProtonPos(0),
	hMCPtPionNeg(0),hMCPtKaonNeg(0),hMCPtProtonNeg(0),
	hTPCRecTracksPionPos(0), hTPCRecTracksKaonPos(0), hTPCRecTracksProtonPos(0),
	hTPCRecTracksPionNeg(0), hTPCRecTracksKaonNeg(0), hTPCRecTracksProtonNeg(0),
	hTOFRecTracksPionPos(0), hTOFRecTracksKaonPos(0), hTOFRecTracksProtonPos(0),
	hTOFRecTracksPionNeg(0), hTOFRecTracksKaonNeg(0), hTOFRecTracksProtonNeg(0),
	hrTPCRecTracksPion(0), hrTPCRecTracksKaon(0), hrTPCRecTracksProton(0),
	hPionTPCDCAxyNegData(0), hPionTPCDCAxyPosData(0), hProtonTPCDCAxyNegData(0), hProtonTPCDCAxyPosData(0),
	hPionTOFDCAxyNegData(0), hPionTOFDCAxyPosData(0), hProtonTOFDCAxyNegData(0), hProtonTOFDCAxyPosData(0),
	hMIPVsEta(0), pMIPVsEta(0), hPlateauVsEta(0), pPlateauVsEta(0)

{
	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hFlatVsPtV0M[i_c] = 0;

		for (Int_t i_eta = 0; i_eta < nEta; ++i_eta){

			hNsigmaPiPos[i_c][i_eta] = 0;
			hNsigmaKPos[i_c][i_eta] = 0;
			hNsigmaPPos[i_c][i_eta] = 0;
			hPtTPCEtaPos[i_c][i_eta] = 0;

			hNsigmaPiNeg[i_c][i_eta] = 0;
			hNsigmaKNeg[i_c][i_eta] = 0;
			hNsigmaPNeg[i_c][i_eta] = 0;
			hPtTPCEtaNeg[i_c][i_eta] = 0;

			hBetaPos[i_c][i_eta] = 0;
			hMomentumTOFEtaPos[i_c][i_eta] = 0;
			hPtTOFEtaPos[i_c][i_eta] = 0;

			hBetaNeg[i_c][i_eta] = 0;
			hMomentumTOFEtaNeg[i_c][i_eta] = 0;
			hPtTOFEtaNeg[i_c][i_eta] = 0;

			hdEdx[i_c][i_eta] = 0;
			hPtrTPC[i_c][i_eta] = 0;
			if (i_c==0){
				hPtVsP[i_eta] = 0;	
				hKaonContamination[i_eta] = 0;
				hKaondEdx[i_eta] = 0;
				hProtonContamination[i_eta] = 0;
				hProtondEdx[i_eta] = 0;
			}
		}
	}

	for(Int_t i = 0; i < 3; ++i){
		hPionTOFDCAxyNeg[i] = 0;
		hProtonTOFDCAxyNeg[i] = 0;
		hPionTOFDCAxyPos[i] = 0;
		hProtonTOFDCAxyPos[i] = 0;
		hPionTPCDCAxyNeg[i] = 0;
		hProtonTPCDCAxyNeg[i] = 0;	
		hPionTPCDCAxyPos[i] = 0;	
		hProtonTPCDCAxyPos[i] = 0;

	}

	DefineInput(0, TChain::Class()); // define the input of the analysis: in this
	// case you take a 'chain' of events
	DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
	// case it's a list of histograms
}
//_____________________________________________________________________________
AliAnalysisTaskFlatenicityPiKp::~AliAnalysisTaskFlatenicityPiKp() {
	// destructor
	if (fOutputList) {
		delete fOutputList; // at the end of your task, it is deleted from memory by
		// calling this function
		fOutputList = 0x0;
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::UserCreateOutputObjects() {

	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if(man){
		AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
		if(inputHandler)fPIDResponse = inputHandler->GetPIDResponse();
	}

	// create track filters
	fTrackFilter = new AliAnalysisFilter("trackFilter");
	AliESDtrackCuts *fCuts = new AliESDtrackCuts();
	fCuts->SetAcceptKinkDaughters(kFALSE);
	fCuts->SetRequireTPCRefit(kTRUE);
	fCuts->SetRequireITSRefit(kTRUE);
	fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
	fCuts->SetDCAToVertex2D(kFALSE);
	fCuts->SetRequireSigmaToVertex(kFALSE);
	fCuts->SetEtaRange(-0.8, 0.8);
	fCuts->SetMinNCrossedRowsTPC(70);
	fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fCuts->SetMaxChi2PerClusterTPC(4);
	fCuts->SetMaxDCAToVertexZ(2);
	fCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
	fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fCuts->SetMaxChi2PerClusterTPC(4);
	fCuts->SetMaxDCAToVertexZ(2);
	fCuts->SetMaxChi2PerClusterITS(36);
	fCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
	fCuts->SetMaxChi2PerClusterITS(36);
	fTrackFilter->AddCuts(fCuts);

	fTrackFilterPID = new AliAnalysisFilter("trackFilterPID");
	AliESDtrackCuts* fCutsPID = new AliESDtrackCuts;
	fCutsPID->SetMinNCrossedRowsTPC(70);
	fCutsPID->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	fCutsPID->SetMaxChi2PerClusterTPC(4);
	fCutsPID->SetAcceptKinkDaughters(kFALSE);
	fCutsPID->SetRequireTPCRefit(kTRUE);
	fCutsPID->SetRequireITSRefit(kTRUE);
	fCutsPID->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
	fCutsPID->SetMaxDCAToVertexZ(2);
	fCutsPID->SetDCAToVertex2D(kFALSE);
	fCutsPID->SetRequireSigmaToVertex(kFALSE);
	fCutsPID->SetMaxChi2PerClusterITS(36);
	fTrackFilterPID->AddCuts(fCutsPID);

	fEtaCalibrationPos = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	const Int_t nPtbins = 56;
	Double_t Ptbins[nPtbins+1] = {
		0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
		0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
		1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
		3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0,
		9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0,
		20.0,22.0,24.0,26.0,30.0};

	// create output objects
	float min_flat = -0.01;
	float max_flat = 1.01;
	int nbins_flat = 204;
	if (fRemoveTrivialScaling) {
		min_flat = -0.1;
		max_flat = 9.9;
		nbins_flat = 2000;
	}

	const Int_t nnSigmabins = 100;
	Double_t nSigmabins[nnSigmabins+1] = {0.0};

	for (Int_t i = 0; i <= nnSigmabins; ++i){
		nSigmabins[i] = -8.0 + i*0.16;
	}

	const Int_t nBetabins   = 100;
	Double_t Betabins[nBetabins+1] = { 0.0 };
	for( Int_t i = 0; i <= nBetabins; ++i){
		Betabins[i] = 0.2+((double)i)/100.0;
	}

	const Int_t dEdxHigh = 140;
	const Int_t dEdxLow = 40;

	const Int_t ndEdxbins = dEdxHigh-dEdxLow;
	Double_t dEdxbins[ndEdxbins+1];

	for(Int_t i = dEdxLow; i <= dEdxHigh; ++i){
		dEdxbins[i-dEdxLow] = i;
	}

	const Int_t nFlatenicitybins = 1000;
	Double_t Flatenicitybins[nFlatenicitybins+1] = {0.0};

	for (Int_t i = 0; i <= nFlatenicitybins; ++i){
		Flatenicitybins[i] = -0.1 + i*0.01;
	}

	OpenFile(1);
	fOutputList =
		new TList(); // this is a list which will contain all of your histograms
	fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

	hFlatV0vsFlatTPC =
		new TH2D("hFlatV0vsFlatTPC", "counter", nbins_flat, min_flat, max_flat,
				nbins_flat, min_flat, max_flat);
	//fOutputList->Add(hFlatV0vsFlatTPC);

	hFlatenicityBefore =
		new TH1F("hFlatenicityBefore", "Counter; Flatenicity; Entries; ", nbins_flat, min_flat, max_flat);
	fOutputList->Add(hFlatenicityBefore);

	hFlatenicity =
		new TH1D("hFlatenicity", "counter", nbins_flat, min_flat, max_flat);
	//fOutputList->Add(hFlatenicity);

	hFlatVsPt =
		new TH2D("hFlatVsPt", "Measured; Flatenicity; #it{p}_{T} (GeV/#it{c})",
				nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
	//fOutputList->Add(hFlatVsPt);

	Int_t SaveThisV0M = -1;
	if (fV0MBin == "0_1") SaveThisV0M = 0;
	else if (fV0MBin == "1_5") SaveThisV0M = 1;
	else if (fV0MBin == "5_10") SaveThisV0M = 2;
	else if (fV0MBin == "10_20") SaveThisV0M = 3;
	else if (fV0MBin == "20_30") SaveThisV0M = 4;
	else if (fV0MBin == "30_40") SaveThisV0M = 5;
	else if (fV0MBin == "40_50") SaveThisV0M = 6;
	else if (fV0MBin == "50_70") SaveThisV0M = 7;
	else SaveThisV0M = 8;

	hPionTPCDCAxyNegData = new TH2F("hPionTPCDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hPionTPCDCAxyPosData = new TH2F("hPionTPCDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hProtonTPCDCAxyNegData = new TH2F("hProtonTPCDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hProtonTPCDCAxyPosData = new TH2F("hProtonTPCDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	

	hPionTOFDCAxyNegData = new TH2F("hPionTOFDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hPionTOFDCAxyPosData = new TH2F("hPionTOFDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hProtonTOFDCAxyNegData = new TH2F("hProtonTOFDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
	hProtonTOFDCAxyPosData = new TH2F("hProtonTOFDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	

	hMIPVsEta = new TH2F("hMIPVsEta","; #eta; dE/dx_{MIP, primary tracks}", 50, -0.8, 0.8, 60-40, 40.0, 60.0);
	pMIPVsEta = new TProfile("pMIPVsEta","; #eta; #LT dE/dx #GT_{MIP, primary tracks}", 50, -0.8, 0.8, 40.0, 60.0);
	hPlateauVsEta = new TH2F("hPlateauVsEta","; #eta; dE/dx_{Plateau, primary tracks}",50, -0.8, 0.8, 50, 60.0, 110.0);
	pPlateauVsEta = new TProfile("pPlateauVsEta","; #eta; #LT dE/dx #GT_{Plateau, primary tracks}", 50, -0.8, 0.8, 60.0, 110.0);

	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		hFlatVsPtV0M[i_c] = new TH2D(Form("hFlatVsPtV0M_c%d", i_c),Form("Measured %1.0f-%1.0f%%V0M; Flatenicity; #it{p}_{T} (GeV/#it{c})",centClass[i_c], centClass[i_c + 1]),nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
		//fOutputList->Add(hFlatVsPtV0M[i_c]);

		for (Int_t i_eta = 0; i_eta < nEta; ++i_eta) {

			hNsigmaPiPos[i_c][i_eta] = new TH3F(Form("hNsigmaPiPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatenicitybins, Flatenicitybins);
			hNsigmaKPos[i_c][i_eta] = new TH3F(Form("hNsigmaKPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatenicitybins, Flatenicitybins);
			hNsigmaPPos[i_c][i_eta] = new TH3F(Form("hNsigmaPPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatenicitybins, Flatenicitybins);
			hPtTPCEtaPos[i_c][i_eta] = new TH2F(Form("hPtTPCEtaPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}; Flatenicity;)", nPtbins, Ptbins, nFlatenicitybins, Flatenicitybins);

			hNsigmaPiNeg[i_c][i_eta] = new TH3F(Form("hNsigmaPiNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatenicitybins, Flatenicitybins);
			hNsigmaKNeg[i_c][i_eta] = new TH3F(Form("hNsigmaKNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatenicitybins, Flatenicitybins);
			hNsigmaPNeg[i_c][i_eta] = new TH3F(Form("hNsigmaPNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins, Ptbins, nnSigmabins, nSigmabins, nFlatenicitybins, Flatenicitybins);
			hPtTPCEtaNeg[i_c][i_eta] = new TH2F(Form("hPtTPCEtaNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); Flatenicity;)", nPtbins, Ptbins, nFlatenicitybins, Flatenicitybins);

			hBetaPos[i_c][i_eta] = new TH3F(Form("hBetaPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p} (GeV/#it{c}); #beta; Flatenicity", nPtbins, Ptbins, nBetabins,Betabins, nFlatenicitybins, Flatenicitybins);
			hMomentumTOFEtaPos[i_c][i_eta] = new TH2F(Form("hMomentumTOFEtaPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatenicitybins, Flatenicitybins);
			hPtTOFEtaPos[i_c][i_eta] = new TH2F(Form("hPtTOFEtaPos_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatenicitybins, Flatenicitybins);

			hBetaNeg[i_c][i_eta] = new TH3F(Form("hBetaNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}; #beta; Flatenicity", nPtbins, Ptbins, nBetabins,Betabins, nFlatenicitybins, Flatenicitybins);
			hMomentumTOFEtaNeg[i_c][i_eta] = new TH2F(Form("hMomentumTOFEtaNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatenicitybins, Flatenicitybins);
			hPtTOFEtaNeg[i_c][i_eta] = new TH2F(Form("hPtTOFEtaNeg_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatenicitybins, Flatenicitybins);

			hdEdx[i_c][i_eta] = new TH3F(Form("hdEdx_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]), ";#it{p} (GeV/#it{c}); dE/dx; Flatenicity", nPtbins, Ptbins, ndEdxbins, dEdxbins, nFlatenicitybins, Flatenicitybins);
			hPtrTPC[i_c][i_eta] = new TH2F(Form("hPtrTPC_c%s_eta_%s",V0MClass[i_c],etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins, Ptbins, nFlatenicitybins, Flatenicitybins);

			if (i_c==0){

				hPtVsP[i_eta] = new TH2F(Form("hPtVsP_eta_%s",etaClass[i_eta]), ";#it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c})", nPtbins, Ptbins, nPtbins, Ptbins);
			}

			if (!fUseMC){ 
				if (SaveThisV0M == i_c){

					fOutputList->Add(hNsigmaPiPos[i_c][i_eta]);
					fOutputList->Add(hNsigmaKPos[i_c][i_eta]);
					fOutputList->Add(hNsigmaPPos[i_c][i_eta]);
					fOutputList->Add(hPtTPCEtaPos[i_c][i_eta]);
					fOutputList->Add(hNsigmaPiNeg[i_c][i_eta]);
					fOutputList->Add(hNsigmaKNeg[i_c][i_eta]);
					fOutputList->Add(hNsigmaPNeg[i_c][i_eta]);
					fOutputList->Add(hPtTPCEtaNeg[i_c][i_eta]);

					fOutputList->Add(hBetaPos[i_c][i_eta]);
					fOutputList->Add(hMomentumTOFEtaPos[i_c][i_eta]);
					fOutputList->Add(hPtTOFEtaPos[i_c][i_eta]);
					fOutputList->Add(hBetaNeg[i_c][i_eta]);
					fOutputList->Add(hMomentumTOFEtaNeg[i_c][i_eta]);
					fOutputList->Add(hPtTOFEtaNeg[i_c][i_eta]);
					fOutputList->Add(hdEdx[i_c][i_eta]);
					fOutputList->Add(hPtrTPC[i_c][i_eta]);

					fOutputList->Add(hPtVsP[i_eta]);

					if(i_eta == 0){
						fOutputList->Add(hPionTPCDCAxyNegData);
						fOutputList->Add(hPionTPCDCAxyPosData);
						fOutputList->Add(hProtonTPCDCAxyNegData);
						fOutputList->Add(hProtonTPCDCAxyPosData);
						fOutputList->Add(hPionTOFDCAxyNegData);
						fOutputList->Add(hPionTOFDCAxyPosData);
						fOutputList->Add(hProtonTOFDCAxyNegData);
						fOutputList->Add(hProtonTOFDCAxyPosData);
						fOutputList->Add(hMIPVsEta);
						fOutputList->Add(pMIPVsEta);
						fOutputList->Add(hPlateauVsEta);
						fOutputList->Add(pPlateauVsEta);
					}
				}
			}
		}
	}

	if (fUseMC) {
		hPtPrimIn =
			new TH1D("hPtPrimIn", "Prim In; #it{p}_{T} (GeV/#it{c}; counts)",
					nPtbins, Ptbins);
		fOutputList->Add(hPtPrimIn);

		hPtPrimOut =
			new TH1D("hPtPrimOut", "Prim Out; #it{p}_{T} (GeV/#it{c}; counts)",
					nPtbins, Ptbins);
		fOutputList->Add(hPtPrimOut);

		hPtSecOut =
			new TH1D("hPtSecOut", "Sec Out; #it{p}_{T} (GeV/#it{c}; counts)",
					nPtbins, Ptbins);
		fOutputList->Add(hPtSecOut);

		hPtOut = new TH1D("hPtOut", "all Out; #it{p}_{T} (GeV/#it{c}; counts)",
				nPtbins, Ptbins);
		fOutputList->Add(hPtOut);

		hFlatenicityMC =
			new TH1D("hFlatenicityMC", "counter", nbins_flat, min_flat, max_flat);
		fOutputList->Add(hFlatenicityMC);
		hFlatResponse =
			new TH2D("hFlatResponse", "; true flat; measured flat", nbins_flat,
					min_flat, max_flat, nbins_flat, min_flat, max_flat);
		fOutputList->Add(hFlatResponse);
		hFlatVsPtMC =
			new TH2D("hFlatVsPtMC", "MC true; Flatenicity; #it{p}_{T} (GeV/#it{c})",
					nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
		fOutputList->Add(hFlatVsPtMC);

		hFlatVsNchMC = new TH2D("hFlatVsNchMC", "; true flat; true Nch", nbins_flat,
				min_flat, max_flat, 100, -0.5, 99.5);
		fOutputList->Add(hFlatVsNchMC);

		hMCPtPionPos = new TH1F("hMCPtPionPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtPionPos);
		hMCPtKaonPos = new TH1F("hMCPtKaonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtKaonPos);
		hMCPtProtonPos = new TH1F("hMCPtProtonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtProtonPos);

		hMCPtPionNeg = new TH1F("hMCPtPionNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtPionNeg);
		hMCPtKaonNeg = new TH1F("hMCPtKaonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtKaonNeg);
		hMCPtProtonNeg = new TH1F("hMCPtProtonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hMCPtProtonNeg);

		hTPCRecTracksPionPos = new TH1F("hTPCRecTracksPionPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksPionPos);
		hTPCRecTracksKaonPos = new TH1F("hTPCRecTracksKaonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksKaonPos);
		hTPCRecTracksProtonPos = new TH1F("hTPCRecTracksProtonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksProtonPos);

		hTPCRecTracksPionNeg = new TH1F("hTPCRecTracksPionNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksPionNeg);
		hTPCRecTracksKaonNeg = new TH1F("hTPCRecTracksKaonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksKaonNeg);
		hTPCRecTracksProtonNeg = new TH1F("hTPCRecTracksProtonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTPCRecTracksProtonNeg);

		hTOFRecTracksPionPos = new TH1F("hTOFRecTracksPionPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksPionPos);
		hTOFRecTracksKaonPos = new TH1F("hTOFRecTracksKaonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksKaonPos);
		hTOFRecTracksProtonPos = new TH1F("hTOFRecTracksProtonPos",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksProtonPos);

		hTOFRecTracksPionNeg = new TH1F("hTOFRecTracksPionNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksPionNeg);
		hTOFRecTracksKaonNeg = new TH1F("hTOFRecTracksKaonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksKaonNeg);
		hTOFRecTracksProtonNeg = new TH1F("hTOFRecTracksProtonNeg",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hTOFRecTracksProtonNeg);

		hrTPCRecTracksPion = new TH1F("hrTPCRecTracksPion",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hrTPCRecTracksPion);
		hrTPCRecTracksKaon = new TH1F("hrTPCRecTracksKaon",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hrTPCRecTracksKaon);
		hrTPCRecTracksProton = new TH1F("hrTPCRecTracksProton",";#it{p}_{T} (GeV/#it{c}); Counts;", nPtbins, Ptbins);
		fOutputList->Add(hrTPCRecTracksProton);

		for(Int_t i = 0; i < 3; ++i){
			hPionTOFDCAxyNeg[i] = new TH2F(Form("hPion%sTOFDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hPionTOFDCAxyNeg[i]);
			hProtonTOFDCAxyNeg[i] = new TH2F(Form("hProton%sTOFDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hProtonTOFDCAxyNeg[i]);
			hPionTOFDCAxyPos[i] = new TH2F(Form("hPion%sTOFDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hPionTOFDCAxyPos[i]);
			hProtonTOFDCAxyPos[i] = new TH2F(Form("hProton%sTOFDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hProtonTOFDCAxyPos[i]);

			hPionTPCDCAxyNeg[i] = new TH2F(Form("hPion%sTPCDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hPionTPCDCAxyNeg[i]);
			hProtonTPCDCAxyNeg[i] = new TH2F(Form("hProton%sTPCDCAxyNeg",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hProtonTPCDCAxyNeg[i]);
			hPionTPCDCAxyPos[i] = new TH2F(Form("hPion%sTPCDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hPionTPCDCAxyPos[i]);
			hProtonTPCDCAxyPos[i] = new TH2F(Form("hProton%sTPCDCAxyPos",ParticleType[i]),"; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 1000, -3.5, 3.5 );	
			fOutputList->Add(hProtonTPCDCAxyPos[i]);

		}

		for (Int_t i_eta = 0; i_eta < nEta; ++i_eta) {

			hKaonContamination[i_eta] = new TH2F(Form("hKaonContamination_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); dEdx", nPtbins, Ptbins, nnSigmabins, nSigmabins);
			fOutputList->Add(hKaonContamination[i_eta]);

			hKaondEdx[i_eta] = new TH2F(Form("hKaondEdx_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); dEdx", nPtbins, Ptbins, nnSigmabins, nSigmabins);
			fOutputList->Add(hKaondEdx[i_eta]);

			hProtonContamination[i_eta] = new TH2F(Form("hProtonContamination_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); dEdx", nPtbins, Ptbins, nnSigmabins, nSigmabins);
			fOutputList->Add(hProtonContamination[i_eta]);

			hProtondEdx[i_eta] = new TH2F(Form("hProtondEdx_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); dEdx", nPtbins, Ptbins, nnSigmabins, nSigmabins);
			fOutputList->Add(hProtondEdx[i_eta]);
		}

	}

	hActivityV0DataSect =
		new TProfile("hActivityV0DataSect", "rec; V0 sector; #LTmultiplicity#GT",
				64, -0.5, 63.5);
	fOutputList->Add(hActivityV0DataSect);
	if (fUseMC) {
		hActivityV0McSect =
			new TProfile("hActivityV0McSect", "true; V0 sector; #LTmultiplicity#GT",
					64, -0.5, 63.5);
		fOutputList->Add(hActivityV0McSect);
	}

	hFlatVsV0M = new TH2F("hFlatVsV0M", "; V0M Percentile; Flatenicity;", nCent, centClass, nbins_flat,
			min_flat, max_flat);
	fOutputList->Add(hFlatVsV0M);

	hCounter = new TH1F("hCounter", "counter", 10, -0.5, 9.5);
	fOutputList->Add(hCounter);

	fEventCuts.AddQAplotsToList(fOutputList);
	PostData(1, fOutputList); // postdata will notify the analysis manager of
	// changes / updates to the
}
//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::UserExec(Option_t *) {

	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}

	fESD = dynamic_cast<AliESDEvent *>(event);

	if (!fESD) {
		Printf("%s:%d ESDEvent not found in Input Manager", (char *)__FILE__,
				__LINE__);
		this->Dump();
		return;
	}

	if (fUseMC) {

		//      E S D
		fMC = dynamic_cast<AliMCEvent *>(MCEvent());
		if (!fMC) {
			Printf("%s:%d MCEvent not found in Input Manager", (char *)__FILE__,
					__LINE__);
			this->Dump();
			return;
		}
		fMCStack = fMC->Stack();
	}

	AliHeader *headerMC;
	Bool_t isGoodVtxPosMC = kFALSE;

	if (fUseMC) {
		headerMC = fMC->Header();
		AliGenEventHeader *genHeader = headerMC->GenEventHeader();
		TArrayF vtxMC(3); // primary vertex  MC
		vtxMC[0] = 9999;
		vtxMC[1] = 9999;
		vtxMC[2] = 9999; // initialize with dummy
		if (genHeader) {
			genHeader->PrimaryVertex(vtxMC);
		}
		if (TMath::Abs(vtxMC[2]) <= 10)
			isGoodVtxPosMC = kTRUE;
	}

	// Trigger selection
	UInt_t fSelectMask = fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask & AliVEvent::kINT7;
	if (!isINT7selected)
		return;

	// Good events
	if (!fEventCuts.AcceptEvent(event)) {
		PostData(1, fOutputList);
		return;
	}

	// Good vertex
	Bool_t hasRecVertex = kFALSE;
	hasRecVertex = HasRecVertex();
	if (!hasRecVertex)
		return;

	// Multiplicity Estimation
	fv0mpercentile = -999;

	fMultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
	if (!fMultSelection)
		cout << "------- No AliMultSelection Object Found --------"
			<< fMultSelection << endl;
	fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
	hCounter->Fill(1);

	for (Int_t i_c = 0; i_c < nCent; ++i_c) {
		if (fv0mpercentile >= centClass[i_c] &&
				fv0mpercentile < centClass[i_c + 1]) {
			fV0Mindex = i_c;
		} else {
			continue;
		}
	}

	hCounter->Fill(3);

	Double_t flatenicity_v0 = GetFlatenicityV0();
	Double_t flatenicity_tpc = GetFlatenicityTPC();
	fFlat = flatenicity_v0; // default V0
	if (fDetFlat == "VO_TPC") {
		fFlat = (flatenicity_v0 + flatenicity_tpc) / 2.0;
	}
	if (fDetFlat == "TPC") {
		fFlat = flatenicity_tpc;
	}

	hFlatV0vsFlatTPC->Fill(flatenicity_tpc, flatenicity_v0);
	fFlatMC = -1;
	if (fUseMC) {
		fFlatMC = GetFlatenicityMC();
		if (fFlatMC >= 0) {
			hFlatenicityMC->Fill(fFlatMC);
			hFlatResponse->Fill(fFlatMC, fFlat);
			MakeMCanalysis();
			MakeMCanalysisPID();
			ElectronsContamination();
		}
	}
	if (fFlat >= 0) {

		hFlatenicityBefore->Fill(fFlat);
		if (flatenicity_v0 < 0.9 && flatenicity_tpc < 0.9) {
			hFlatenicity->Fill(fFlat);
		}
		if (fV0Mindex >= 0) {
			hCounter->Fill(5);
			hFlatVsV0M->Fill(fv0mpercentile, fFlat);
			MakeDataanalysis();
			MakePIDanalysis();
		}
	}

	if (fIsMCclosure) {
		Double_t randomUE = -1;
		gRandom->SetSeed(0);
		randomUE = gRandom->Uniform(0.0, 1.0);
		if (randomUE < 0.5) { // corrections (50% stat.)
			if (isGoodVtxPosMC) {
			}
		} else { // for testing the method
		}
	} else {
		if (fUseMC) {
			if (isGoodVtxPosMC) {
			}
		} else {
		}
	}

	PostData(1, fOutputList); // stream the result of this event to the output
	// manager which will write it to a file
}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::Terminate(Option_t *) {}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::MakeDataanalysis() {

	// rec
	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) {

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;
		hFlatVsPt->Fill(fFlat, esdtrack->Pt());
		hFlatVsPtV0M[fV0Mindex]->Fill(fFlat, esdtrack->Pt());
	}
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::MakePIDanalysis() {


	Int_t nESDTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = static_cast<AliESDtrack *>(fESD->GetTrack(iT));

		if (!esdTrack)
			continue;

		if (!fTrackFilterPID->IsSelected(esdTrack))
			continue;

		if (TMath::Abs(esdTrack->Eta()) > fEtaCut )
			continue;

		if (esdTrack->Pt() < fPtMin)
			continue;

		Double_t momentum = esdTrack->P();
		Double_t pt = esdTrack->Pt();
		Double_t eta = esdTrack->Eta();
		Double_t phi = esdTrack->Phi();
		Float_t dEdx = esdTrack->GetTPCsignal();
		UShort_t Ncl = esdTrack->GetTPCsignalN();
		Short_t charge = esdTrack->Charge();
		Float_t DCAxy = 0.0;
		Float_t dcaz = 0.0;
		esdTrack->GetImpactParameters(DCAxy,dcaz);

		Float_t nSigmaPi = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion);
		Float_t nSigmaK = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon);
		Float_t nSigmaP = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
		Float_t nSigmaPiTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion);
		Float_t nSigmaPTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton);

		Int_t nh = -1;
		if (TMath::Abs(eta)<0.2)
			nh = 0;
		else if (TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
			nh = 1;
		else if (TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
			nh = 2;
		else if (TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
			nh = 3;

		if (nh<0)
			continue;

		//! These are used to measure the Feed-Down correction

		if( charge < 0 ){
			if( (nSigmaPi >= -2.0) && (nSigmaPi <= 2.0) ) { hPionTPCDCAxyNegData->Fill(pt,DCAxy); }
			if( (nSigmaP >= -2.0) && (nSigmaP <= 2.0) ) { hProtonTPCDCAxyNegData->Fill(pt,DCAxy); }

			if( TOFPID(esdTrack) ){
				if( TMath::Sqrt(nSigmaPi*nSigmaPi + nSigmaPiTOF*nSigmaPiTOF) < 2.0) { hPionTOFDCAxyNegData->Fill(pt,DCAxy); }
				if( TMath::Sqrt(nSigmaP*nSigmaP + nSigmaPTOF*nSigmaPTOF) < 2.0) { hProtonTOFDCAxyNegData->Fill(pt,DCAxy); }
			}

		}else{
			if( (nSigmaPi >= -2.0) && (nSigmaPi <= 2.0) ) { hPionTPCDCAxyPosData->Fill(pt,DCAxy); }
			if( (nSigmaP >= -2.0) && (nSigmaP <= 2.0)) { hProtonTPCDCAxyPosData->Fill(pt,DCAxy); }

			if( TOFPID(esdTrack) ){
				if(TMath::Sqrt(nSigmaPi*nSigmaPi + nSigmaPiTOF*nSigmaPiTOF) < 2.0) { hPionTOFDCAxyPosData->Fill(pt,DCAxy); }
				if(TMath::Sqrt(nSigmaP*nSigmaP + nSigmaPTOF*nSigmaPTOF) < 2.0) { hProtonTOFDCAxyPosData->Fill(pt,DCAxy); }
			}
		}


		Double_t maxDCAxy = 10.0;
		maxDCAxy = fcutDCAxy->Eval(pt);

		if (TMath::Abs(DCAxy) > maxDCAxy )
			continue;

		if( charge < 0.0){
			hNsigmaPiNeg[fV0Mindex][nh]->Fill(pt,nSigmaPi,fFlat);
			hNsigmaKNeg[fV0Mindex][nh]->Fill(pt,nSigmaK,fFlat);
			hNsigmaPNeg[fV0Mindex][nh]->Fill(pt,nSigmaP,fFlat);
			hPtTPCEtaNeg[fV0Mindex][nh]->Fill(pt,fFlat);
		}else{
			hNsigmaPiPos[fV0Mindex][nh]->Fill(pt,nSigmaPi,fFlat);
			hNsigmaKPos[fV0Mindex][nh]->Fill(pt,nSigmaK,fFlat);
			hNsigmaPPos[fV0Mindex][nh]->Fill(pt,nSigmaP,fFlat);
			hPtTPCEtaPos[fV0Mindex][nh]->Fill(pt,fFlat);
		} 

		if( TOFPID(esdTrack) ){

			Double_t trkLength = esdTrack->GetIntegratedLength();
			Double_t beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);

			if(esdTrack->Charge() < 0.0){
				hBetaNeg[fV0Mindex][nh]->Fill(momentum,beta,fFlat);
				hMomentumTOFEtaNeg[fV0Mindex][nh]->Fill(momentum,fFlat);
				hPtTOFEtaNeg[fV0Mindex][nh]->Fill(pt,fFlat);
			}else{ 
				hBetaPos[fV0Mindex][nh]->Fill(momentum,beta,fFlat);
				hMomentumTOFEtaPos[fV0Mindex][nh]->Fill(momentum,fFlat);
				hPtTOFEtaPos[fV0Mindex][nh]->Fill(pt,fFlat);
			}
		}

		hPtVsP[nh]->Fill(momentum,pt);

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), 1.0, fcutLow, fcutHigh))
			continue;

		if (Ncl < fNcl)
			continue;

		if (fdEdxCalibrated) dEdx *= 50.0/EtaCalibration(eta);

		hdEdx[fV0Mindex][nh]->Fill(momentum,dEdx,fFlat);
		hPtrTPC[fV0Mindex][nh]->Fill(pt,fFlat);

		Bool_t IsTOFout=kFALSE;
		if ((esdTrack->GetStatus()&AliESDtrack::kTOFout)==0)
			IsTOFout=kTRUE;

		Float_t lengthtrack = esdTrack->GetIntegratedLength();
		Float_t timeTOF = esdTrack->GetTOFsignal();
		Double_t inttime[5]={0,0,0,0,0};
		esdTrack->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
		Float_t beta = -99;
		if ( !IsTOFout ){
			if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
				beta = inttime[0] / timeTOF;
		}

		if ( momentum <= 0.6 && momentum >= 0.4 ){ //only p:0.4-0.6 GeV, pion MIP
			if( dEdx < 60.0 && dEdx > 40.0 ){
				hMIPVsEta->Fill(eta,dEdx);
				pMIPVsEta->Fill(eta,dEdx);
			}
			if( dEdx > 70.0 && dEdx < 90.0 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta->Fill(eta,dEdx);
					pPlateauVsEta->Fill(eta,dEdx);
				}
			}
		}


	}//end of track loop


}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::MakeMCanalysis() {

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)
			continue;
		if (particle->Pt() < fPtMin)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;
		hFlatVsPtMC->Fill(fFlatMC, particle->Pt());
		hPtPrimIn->Fill(particle->Pt());
	}
	// rec
	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) {

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		if (TMath::Abs(esdtrack->Eta()) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;
		hPtOut->Fill(esdtrack->Pt());
		Int_t mcLabel = -1;
		mcLabel = TMath::Abs(esdtrack->GetLabel());
		if (fMC->IsPhysicalPrimary(mcLabel)) {
			hPtPrimOut->Fill(esdtrack->Pt());
		} else {
			hPtSecOut->Fill(esdtrack->Pt());
		}
	}
}

//______________________________________________________________________________

void AliAnalysisTaskFlatenicityPiKp::MakeMCanalysisPID() {

	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)
			continue;
		if (particle->Pt() < fPtMin)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;

		hFlatVsPtMC->Fill(fFlatMC, particle->Pt());
		hPtPrimIn->Fill(particle->Pt());

		Int_t pdgCode = particle->PdgCode();
		Int_t pidCode = GetPidCode(pdgCode);

		if(particle->Charge() < 0.0){
			if (pidCode==1) hMCPtPionNeg->Fill(particle->Pt());	
			if (pidCode==2) hMCPtKaonNeg->Fill(particle->Pt());	
			if (pidCode==3) hMCPtProtonNeg->Fill(particle->Pt());	
		}else{
			if (pidCode==1) hMCPtPionPos->Fill(particle->Pt());	
			if (pidCode==2) hMCPtKaonPos->Fill(particle->Pt());	
			if (pidCode==3) hMCPtProtonPos->Fill(particle->Pt());	

		}
	}

	// rec
	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; iT++){

		AliESDtrack *esdTrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdTrack)
			continue;

		if (!fTrackFilterPID->IsSelected(esdTrack))
			continue;

		if (TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		if (esdTrack->Pt() < fPtMin)
			continue;

		Int_t mcLabel = -1;
		mcLabel = TMath::Abs(esdTrack->GetLabel());

		AliMCParticle *mcTrack = 0;
		mcTrack = (AliMCParticle*)fMC->GetTrack(mcLabel);

		if (!mcTrack) 
			continue;

		if ( TMath::Abs(esdTrack->Charge())==0 )
			continue;

		Int_t pdgCode = mcTrack->PdgCode();
		Int_t pidCode = GetPidCode(pdgCode);
		Double_t pt = esdTrack->Pt();

		Double_t maxDCAxy = 10.0;
		maxDCAxy = fcutDCAxy->Eval(pt);
		Float_t DCAxy = 0.0;
		Float_t dcaz = 0.0;
		esdTrack->GetImpactParameters(DCAxy,dcaz);

		if( fMC->IsPhysicalPrimary(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hPionTOFDCAxyNeg[0]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyNeg[0]->Fill(pt,DCAxy);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[0]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyPos[0]->Fill(pt,DCAxy);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[0]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyNeg[0]->Fill(pt,DCAxy);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[0]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyPos[0]->Fill(pt,DCAxy);	
			}
		}

		if( fMC->IsSecondaryFromMaterial(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hPionTOFDCAxyNeg[1]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyNeg[1]->Fill(pt,DCAxy);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[1]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyPos[1]->Fill(pt,DCAxy);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[1]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyNeg[1]->Fill(pt,DCAxy);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[1]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyPos[1]->Fill(pt,DCAxy);	
			}
		}

		if( fMC->IsSecondaryFromWeakDecay(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hPionTOFDCAxyNeg[2]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyNeg[2]->Fill(pt,DCAxy);	
				}else{ 
					if (pidCode==1) hPionTOFDCAxyPos[2]->Fill(pt,DCAxy);	
					if (pidCode==3) hProtonTOFDCAxyPos[2]->Fill(pt,DCAxy);	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hPionTPCDCAxyNeg[2]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyNeg[2]->Fill(pt,DCAxy);	
			}else{
				if (pidCode==1) hPionTPCDCAxyPos[2]->Fill(pt,DCAxy);	
				if (pidCode==3) hProtonTPCDCAxyPos[2]->Fill(pt,DCAxy);	
			}
		}


		//! DCAxy cut to select primaries
		if (TMath::Abs(DCAxy) > maxDCAxy )
			continue;

		if( fMC->IsPhysicalPrimary(mcLabel) ){

			if( TOFPID(esdTrack) ){

				if(esdTrack->Charge() < 0.0){
					if (pidCode==1) hTOFRecTracksPionNeg->Fill(esdTrack->Pt());	
					if (pidCode==2)	hTOFRecTracksKaonNeg->Fill(esdTrack->Pt());	
					if (pidCode==3) hTOFRecTracksProtonNeg->Fill(esdTrack->Pt());	
				}else{ 
					if (pidCode==1) hTOFRecTracksPionPos->Fill(esdTrack->Pt());	
					if (pidCode==2)	hTOFRecTracksKaonPos->Fill(esdTrack->Pt());	
					if (pidCode==3) hTOFRecTracksProtonPos->Fill(esdTrack->Pt());	
				}
			}

			if(esdTrack->Charge() < 0.0){
				if (pidCode==1) hTPCRecTracksPionNeg->Fill(esdTrack->Pt());	
				if (pidCode==2) hTPCRecTracksKaonNeg->Fill(esdTrack->Pt());	
				if (pidCode==3) hTPCRecTracksProtonNeg->Fill(esdTrack->Pt());	
			}else{
				if (pidCode==1) hTPCRecTracksPionPos->Fill(esdTrack->Pt());	
				if (pidCode==2) hTPCRecTracksKaonPos->Fill(esdTrack->Pt());	
				if (pidCode==3) hTPCRecTracksProtonPos->Fill(esdTrack->Pt());	
			}

			if (esdTrack->GetTPCsignalN() < fNcl)
				continue;

			if (!PhiCut(esdTrack->Pt(), esdTrack->Phi(), esdTrack->Charge(), 1, fcutLow, fcutHigh))
				continue;

			if (pidCode==1) hrTPCRecTracksPion->Fill(esdTrack->Pt());	
			if (pidCode==2) hrTPCRecTracksKaon->Fill(esdTrack->Pt());	
			if (pidCode==3) hrTPCRecTracksProton->Fill(esdTrack->Pt());	
		}
	}

}

//______________________________________________________________________________

Double_t AliAnalysisTaskFlatenicityPiKp::GetFlatenicityTPC() {

	const int nRings2 = 4;
	const int nSectors2 = 8;
	const int nCells2 = nRings2 * nSectors2;
	float maxEta2[nRings2] = {-0.4, 0.0, +0.4, +0.8};
	float minEta2[nRings2] = {-0.8, -0.4, +0.0, +0.4};
	float maxPhi2[nSectors2] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
	float minPhi2[nSectors2] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
	float RhoLattice2[nCells2];
	for (int iCh = 0; iCh < nCells2; iCh++) {
		RhoLattice2[iCh] = 0.0;
	}
	int mult_glob = 0;
	Int_t nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) {

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		float eta_a = esdtrack->Eta();
		float phi_a = esdtrack->Phi();

		if (TMath::Abs(eta_a) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;
		int i_ch = 0;
		for (int ir = 0; ir < nRings2; ir++) {
			for (int is = 0; is < nSectors2; is++) {
				if (eta_a >= minEta2[ir] && eta_a < maxEta2[ir] &&
						phi_a >= minPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2) &&
						phi_a < maxPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2)) {
					RhoLattice2[i_ch]++;
					mult_glob++;
				}
				i_ch++;
			}
		}
	}

	double mRho_glob = 0;
	for (int iCell = 0; iCell < nCells2; ++iCell) {
		mRho_glob += 1.0 * RhoLattice2[iCell];
	}
	// average activity per cell
	mRho_glob /= (1.0 * nCells2);
	// get sigma
	double sRho_glob_tmp = 0;
	for (int iCell = 0; iCell < nCells2; ++iCell) {
		sRho_glob_tmp += TMath::Power(1.0 * RhoLattice2[iCell] - mRho_glob, 2);
	}
	sRho_glob_tmp /= (1.0 * nCells2 * nCells2);
	double sRho_glob = TMath::Sqrt(sRho_glob_tmp);
	float flatenicity_glob = 9999;
	if (mRho_glob > 0) {
		if (fRemoveTrivialScaling) {
			flatenicity_glob = TMath::Sqrt(mult_glob) * sRho_glob / mRho_glob;
		} else {
			flatenicity_glob = sRho_glob / mRho_glob;
		}
	}

	return flatenicity_glob;
}

//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityPiKp::GetFlatenicityV0() {

	AliVVZERO *lVV0 = 0x0;
	AliVEvent *lVevent = 0x0;
	lVevent = dynamic_cast<AliVEvent *>(InputEvent());
	if (!lVevent) {
		AliWarning("ERROR: ESD / AOD event not available \n");
		return -1;
	}
	// Get VZERO Information for multiplicity later
	lVV0 = lVevent->GetVZEROData();
	if (!lVV0) {
		AliError("AliVVZERO not available");
		return 9999;
	}
	// Flatenicity calculation
	const Int_t nRings = 4;
	const Int_t nSectors = 8;
	Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
	Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
	Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
	Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
	// Grid
	const Int_t nCells = nRings * 2 * nSectors;
	Float_t RhoLattice[nCells];
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		RhoLattice[iCh] = 0.0;
	}

	Int_t nringA = 0;
	Int_t nringC = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		Float_t detaV0 = -1;
		Float_t mult = lVV0->GetMultiplicity(iCh);
		if (iCh < 32) { // V0C
			if (iCh < 8) {
				nringC = 0;
			} else if (iCh >= 8 && iCh < 16) {
				nringC = 1;
			} else if (iCh >= 16 && iCh < 24) {
				nringC = 2;
			} else {
				nringC = 3;
			}
			detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
		} else { // V0A
			if (iCh < 40) {
				nringA = 0;
			} else if (iCh >= 40 && iCh < 48) {
				nringA = 1;
			} else if (iCh >= 48 && iCh < 56) {
				nringA = 2;
			} else {
				nringA = 3;
			}
			detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
		}
		RhoLattice[iCh] =
			mult / detaV0; // needed to consider the different eta coverage
	}
	// Filling histos with mult info
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]);
	}
	Float_t mRho = 0;
	Float_t flatenicity = -1;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		mRho += RhoLattice[iCh];
	}
	Float_t multiplicityV0M = mRho;
	// average activity per cell
	mRho /= (1.0 * nCells);
	// get sigma
	Double_t sRho_tmp = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
	}
	sRho_tmp /= (1.0 * nCells * nCells);
	Float_t sRho = TMath::Sqrt(sRho_tmp);
	if (mRho > 0) {
		if (fRemoveTrivialScaling) {
			flatenicity = TMath::Sqrt(multiplicityV0M) * sRho / mRho;
		} else {
			flatenicity = sRho / mRho;
		}
	} else {
		flatenicity = 9999;
	}
	return flatenicity;
}
//______________________________________________________________________________
Double_t AliAnalysisTaskFlatenicityPiKp::GetFlatenicityMC() {

	// Flatenicity calculation
	const Int_t nRings = 8;
	Float_t maxEta[nRings] = {-3.2, -2.7, -2.2, -1.7, 5.1, 4.5, 3.9, 3.4};
	Float_t minEta[nRings] = {-3.7, -3.2, -2.7, -2.2, 4.5, 3.9, 3.4, 2.8};

	const Int_t nSectors = 8;
	Float_t PhiBins[nSectors + 1];
	Float_t deltaPhi = (2.0 * TMath::Pi()) / (1.0 * nSectors);
	for (int i_phi = 0; i_phi < nSectors + 1; ++i_phi) {
		PhiBins[i_phi] = 0;
		if (i_phi < nSectors) {
			PhiBins[i_phi] = i_phi * deltaPhi;
		} else {
			PhiBins[i_phi] = 2.0 * TMath::Pi();
		}
	}

	// Grid
	const Int_t nCells = nRings * nSectors;
	Float_t RhoLattice[nCells];
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		RhoLattice[iCh] = 0.0;
	}

	Int_t nMult = 0;
	for (Int_t i = 0; i < fMC->GetNumberOfTracks(); ++i) {

		AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
		if (!particle)
			continue;
		if (!fMC->IsPhysicalPrimary(i))
			continue;
		if (particle->Pt() <= 0.0)
			continue;
		if (TMath::Abs(particle->Charge()) < 0.1)
			continue;
		Double_t phi = particle->Phi();
		Double_t eta = particle->Eta();

		Int_t i_segment = 0;
		for (int i_eta = 0; i_eta < nRings; ++i_eta) {

			for (int i_phi = 0; i_phi < nSectors; ++i_phi) {

				if (eta >= minEta[i_eta] && eta < maxEta[i_eta] &&
						phi >= PhiBins[i_phi] && phi < PhiBins[i_phi + 1]) {
					nMult++;
					RhoLattice[i_segment] += 1.0;
				}
				i_segment++;
			}
		}
	}

	Int_t i_segment = 0;
	for (int i_eta = 0; i_eta < nRings; ++i_eta) {
		for (int i_phi = 0; i_phi < nSectors; ++i_phi) {
			Float_t deltaEta = TMath::Abs(maxEta[i_eta] - minEta[i_eta]);
			RhoLattice[i_segment] /= deltaEta;
			// Filling histos with mult info
			hActivityV0McSect->Fill(i_segment, RhoLattice[i_segment]);
			i_segment++;
		}
	}

	Float_t mRho = 0;
	Float_t flatenicity = -1;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		mRho += RhoLattice[iCh];
	}
	// average activity per cell
	mRho /= (1.0 * nCells);
	// get sigma
	Float_t sRho_tmp = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
	}
	sRho_tmp /= (1.0 * nCells * nCells);
	Float_t sRho = TMath::Sqrt(sRho_tmp);
	if (mRho > 0) {
		if (fRemoveTrivialScaling) {
			flatenicity = TMath::Sqrt(1.0 * nMult) * sRho / mRho;
		} else {
			flatenicity = sRho / mRho;
		}
	} else {
		sRho = 9999;
	}
	hFlatVsNchMC->Fill(flatenicity, nMult);
	return flatenicity;
}

Bool_t AliAnalysisTaskFlatenicityPiKp::HasRecVertex() {

	float fMaxDeltaSpdTrackAbsolute = 0.5f;
	float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
	float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
	float fMaxResolutionSPDvertex = 0.25f;
	float fMaxDispersionSPDvertex = 1.e14f;

	Bool_t fRequireTrackVertex = true;
	unsigned long fFlag;
	fFlag = BIT(AliEventCuts::kNoCuts);

	const AliVVertex *vtTrc = fESD->GetPrimaryVertex();
	bool isTrackV = true;
	if (vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ())
		isTrackV = false;
	const AliVVertex *vtSPD = fESD->GetPrimaryVertexSPD();

	if (vtSPD->GetNContributors() > 0)
		fFlag |= BIT(AliEventCuts::kVertexSPD);

	if (vtTrc->GetNContributors() > 1 && isTrackV)
		fFlag |= BIT(AliEventCuts::kVertexTracks);

	if (((fFlag & BIT(AliEventCuts::kVertexTracks)) || !fRequireTrackVertex) &&
			(fFlag & BIT(AliEventCuts::kVertexSPD)))
		fFlag |= BIT(AliEventCuts::kVertex);

	const AliVVertex *&vtx =
		bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
	AliVVertex *fPrimaryVertex = const_cast<AliVVertex *>(vtx);
	if (!fPrimaryVertex)
		return kFALSE;

	/// Vertex quality cuts
	double covTrc[6], covSPD[6];
	vtTrc->GetCovarianceMatrix(covTrc);
	vtSPD->GetCovarianceMatrix(covSPD);
	double dz = bool(fFlag & AliEventCuts::kVertexSPD) &&
		bool(fFlag & AliEventCuts::kVertexTracks)
		? vtTrc->GetZ() - vtSPD->GetZ()
		: 0.; /// If one of the two vertices is not available this cut
	/// is always passed.
	double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
	double errTrc =
		bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
	double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
	/// vertex dispersion for run1, only for ESD, AOD code to be added here
	const AliESDVertex *vtSPDESD = dynamic_cast<const AliESDVertex *>(vtSPD);
	double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
	if ((TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute &&
				nsigTot <= fMaxDeltaSpdTrackNsigmaSPD &&
				nsigTrc <=
				fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
			(!vtSPD->IsFromVertexerZ() ||
			 TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
			(!vtSPD->IsFromVertexerZ() ||
			 vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for
			/// run1, only for ESD
	   ) // quality cut on vertexer SPD z
		fFlag |= BIT(AliEventCuts::kVertexQuality);

	Bool_t hasVtx = (TESTBIT(fFlag, AliEventCuts::kVertex)) &&
		(TESTBIT(fFlag, AliEventCuts::kVertexQuality));

	return hasVtx;
}

Bool_t AliAnalysisTaskFlatenicityPiKp::TOFPID(AliESDtrack * track) {
	UInt_t status;
	status=track->GetStatus();

	if (!(status & AliESDtrack::kTOFout) || !(status & AliESDtrack::kTIME))
		return kFALSE;

	if (track->GetIntegratedLength() < 350.)
		return kFALSE;

	if (TMath::Abs(track->GetTOFsignalDx()) > 10.0 || TMath::Abs(track->GetTOFsignalDz()) > 10.0)
		return kFALSE;

	return kTRUE;
}

Double_t AliAnalysisTaskFlatenicityPiKp::EtaCalibration(const Double_t &eta) {

	double aPos = 0.0;
	double bPos = 0.0;
	double cPos = 0.0;
	double dPos = 0.0;
	double ePos = 0.0;
	double fPos = 0.0;
	double gPos = 0.0;
	double hPos = 0.0;

	double aNeg = 0.0;
	double bNeg = 0.0;
	double cNeg = 0.0;
	double dNeg = 0.0;
	double eNeg = 0.0;
	double fNeg = 0.0;
	double gNeg = 0.0;
	double hNeg = 0.0;

	if (fPeriod == "16l"){
		aPos = 49.9216; bPos = 3.07252; cPos = -42.8044; dPos = 259.666; ePos = -910.432; fPos = 1776.09; gPos = -1740.65; hPos = 662.232;
		aNeg = 49.9732; bNeg = 4.03575; cNeg = 65.6189;  dNeg = 374.429; eNeg = 951.459;  fNeg = 1153.75; gNeg = 618.493;  hNeg = 100.499;
	} else if(fPeriod == "16k"){
		aPos = 49.9421; bPos = 2.3446; cPos = -41.2765; dPos = 279.695; ePos = -1027.73; fPos = 2022.84; gPos = -1967.79; hPos = 738.823;
		aNeg = 50.0477; bNeg = 8.27344; cNeg = 125.29;  dNeg = 736.8;   eNeg = 2057.75;  fNeg = 2935.38; gNeg = 2064.03;  hNeg = 565.983;
	} else if(fPeriod == "16deghijop"){
		aPos = 49.9743; bPos = 2.3388; cPos = -44.1496; dPos = 296.029; ePos = -1056.56; fPos = 2031.44; gPos = -1946.51; hPos = 723.89;
		aNeg = 50.0329; bNeg = 6.99747; cNeg = 107.168;  dNeg = 649.001; eNeg = 1875.17;  fNeg = 2785.78; gNeg = 2063.77;  hNeg = 606.868;
	} else if(fPeriod == "17data"){
		aPos = 49.6097; bPos = 0.922856; cPos = -6.57484; dPos = 65.3117; ePos = -372.142; fPos = 950.451; gPos = -1085.27; hPos = 458.144;
		aNeg = 49.6555; bNeg = 6.98696; cNeg = 102.734;  dNeg = 566.424; eNeg = 1513.64;  fNeg = 2092.01; gNeg = 1429.32;  hNeg = 375.642;
	} else{
		aPos = 49.6975; bPos = 2.32535; cPos = -42.6516; dPos = 283.058; ePos = -1009.58; fPos = 1945.89; gPos = -1871.23; hPos = 698.552;
		aNeg = 49.8071; bNeg = 9.78466; cNeg = 120.018;  dNeg = 603.325; eNeg = 1470.92;  fNeg = 1819.63; gNeg = 1073.82;  hNeg = 230.142;
	}

	for (Int_t i = 0; i < 8; ++i){
		fEtaCalibrationPos->SetParameter(i,0);
		fEtaCalibrationNeg->SetParameter(i,0);
	}

	Double_t Calibration = 0.0;

	if (eta<=0.0){
		fEtaCalibrationNeg->SetParameter(0,aNeg);
		fEtaCalibrationNeg->SetParameter(1,bNeg);
		fEtaCalibrationNeg->SetParameter(2,cNeg);
		fEtaCalibrationNeg->SetParameter(3,dNeg);
		fEtaCalibrationNeg->SetParameter(4,eNeg);
		fEtaCalibrationNeg->SetParameter(5,fNeg);
		fEtaCalibrationNeg->SetParameter(6,gNeg);
		fEtaCalibrationNeg->SetParameter(7,hNeg);

		Calibration = fEtaCalibrationNeg->Eval(eta);
	} else {
		fEtaCalibrationPos->SetParameter(0,aPos);
		fEtaCalibrationPos->SetParameter(1,bPos);
		fEtaCalibrationPos->SetParameter(2,cPos);
		fEtaCalibrationPos->SetParameter(3,dPos);
		fEtaCalibrationPos->SetParameter(4,ePos);
		fEtaCalibrationPos->SetParameter(5,fPos);
		fEtaCalibrationPos->SetParameter(6,gPos);
		fEtaCalibrationPos->SetParameter(7,hPos);

		Calibration = fEtaCalibrationPos->Eval(eta);
	}

	return Calibration;
}

void AliAnalysisTaskFlatenicityPiKp::ElectronsContamination() {

	Int_t iTracks(fESD->GetNumberOfTracks());          

	for (Int_t iT = 0; iT < iTracks; ++iT) {                

		AliESDtrack *esdTrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdTrack)
			continue;

		if (!fTrackFilterPID->IsSelected(esdTrack))
			continue;

		if (TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		if (esdTrack->Pt() < fPtMin)
			continue;

		Double_t maxDCAxy = 10.0;
		maxDCAxy = fcutDCAxy->Eval(esdTrack->Pt());
		Float_t DCAxy = 0.0;
		Float_t dcaz = 0.0;
		esdTrack->GetImpactParameters(DCAxy,dcaz);

		//! DCAxy cut to select primaries
		if (TMath::Abs(DCAxy) > maxDCAxy )
			continue;

		Int_t mcLabel = -1;
		mcLabel = TMath::Abs(esdTrack->GetLabel());

		if (!fMC->IsPhysicalPrimary(mcLabel) )
			continue;

		AliMCParticle *mcTrack = 0;
		mcTrack = (AliMCParticle*)fMC->GetTrack(mcLabel);

		if (!mcTrack) 
			continue;

		if (TMath::Abs(esdTrack->Charge())==0)
			continue;

		Int_t pdgCode = mcTrack->PdgCode();
		Int_t pidCode = GetPidCode(pdgCode);

		Int_t nh = -1;
		Double_t eta = esdTrack->Eta();
		if (TMath::Abs(eta)<0.2)
			nh = 0;
		else if (TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
			nh = 1;
		else if (TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
			nh = 2;
		else
			nh = 3;

		if ( nh < 0 )
			continue;

		Float_t nsigk = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon);
		Float_t nsigp = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
		Double_t pt = esdTrack->Pt();

		if (pidCode == 7){
			hKaonContamination[nh]->Fill(pt,nsigk);
			hProtonContamination[nh]->Fill(pt,nsigp);
		}

		if (pidCode == 2)
			hKaondEdx[nh]->Fill(pt,nsigk);

		if (pidCode == 3)
			hProtondEdx[nh]->Fill(pt,nsigp);

	} 

}


Bool_t AliAnalysisTaskFlatenicityPiKp::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t mag, TF1* phiCutLow, TF1* phiCutHigh) {

	if(pt < 2.0)
		return kTRUE;

	if(fESD->GetMagneticField() < 0)    // for negatve polarity field
		phi = TMath::TwoPi() - phi;
	if(q < 0) // for negatve charge
		phi = TMath::TwoPi()-phi;

	phi += TMath::Pi()/18.0; // to center gap in the middle
	phi = fmod(phi, TMath::Pi()/9.0);

	if(phi<phiCutHigh->Eval(pt)
			&& phi>phiCutLow->Eval(pt))
		return kFALSE; // reject track

	//	hPhi[fCentClass]->Fill(pt, phi);

	return kTRUE;
}

Int_t AliAnalysisTaskFlatenicityPiKp::GetPidCode(Int_t pdgCode) {
	// return our internal code for pions, kaons, and protons

	Int_t pidCode = 6;

	switch (TMath::Abs(pdgCode)) {
		case 211:
			pidCode = 1; // pion
			break;
		case 321:
			pidCode = 2; // kaon
			break;
		case 2212:
			pidCode = 3; // proton
			break;
		case 310:
			pidCode = 4; // K0s
			break;
		case 3122:
			pidCode = 5; // lambda
			break;
		case 333:
			pidCode = 6; // phi
			break;
		case 11:
			pidCode = 7; // electron
			break;
		default:
			pidCode = 8;  // something else?
	};

	return pidCode;
}

