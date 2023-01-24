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
#include <AliKFVertex.h>
#include <AliKFParticle.h>
using std::cout;
using std::endl;

#include "AliAnalysisTaskFlatenicityPiKp.h"

static const Int_t nCent = 9;
static const Int_t nEta = 4;

static Double_t centClass[nCent + 1] = {0.0,  1.0,  5.0,  10.0, 20.0,
	30.0, 40.0, 50.0, 70.0, 100.0};

static const Char_t* etaClass[nEta] = {"02","24","46","68"};
static const Char_t* ParticleType[3] = {"Primaries","MaterialInt","WeakDecays"};

static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskFlatenicityPiKp) // classimp: necessary for root

AliAnalysisTaskFlatenicityPiKp::AliAnalysisTaskFlatenicityPiKp()
	: AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
	fUseMC(kFALSE), fV0MMultiplicity(-1.0), fDetFlat("V0"), fV0MBin("0_1"), fIsMCclosure(kFALSE),
	fDeltaV0(kTRUE), fRemoveTrivialScaling(kFALSE), fnGen(-1), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), fdEdxCalibrated(kTRUE),
	fSaveDCAxyHistograms(kFALSE), fEtaCalibrationPos(0x0), fEtaCalibrationNeg(0x0),
	fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fPeriod("16l"),
	ftrackmult08(0), fv0mpercentile(0), fFlat(-1), fFlatTPC(-1.), fFlatMC(-1),
	fMultSelection(0x0), hFlat(0),
	hFlatenicityMC(0), hFlatenicityMCRec(0), hFlatResponse(0), hFlatVsPtMC(0),
	hActivityV0McSect(0), hFlatVsNchMC(0),
	hMCPtPionPos(0),hMCPtKaonPos(0),hMCPtProtonPos(0),
	hMCPtPionNeg(0),hMCPtKaonNeg(0),hMCPtProtonNeg(0),
	hTPCRecTracksPionPos(0), hTPCRecTracksKaonPos(0), hTPCRecTracksProtonPos(0),
	hTPCRecTracksPionNeg(0), hTPCRecTracksKaonNeg(0), hTPCRecTracksProtonNeg(0),
	hTOFRecTracksPionPos(0), hTOFRecTracksKaonPos(0), hTOFRecTracksProtonPos(0),
	hTOFRecTracksPionNeg(0), hTOFRecTracksKaonNeg(0), hTOFRecTracksProtonNeg(0),
	hrTPCRecTracksPion(0), hrTPCRecTracksKaon(0), hrTPCRecTracksProton(0),
	hPionTPCDCAxyNegData(0), hPionTPCDCAxyPosData(0), hProtonTPCDCAxyNegData(0), hProtonTPCDCAxyPosData(0),
	hPionTOFDCAxyNegData(0), hPionTOFDCAxyPosData(0), hProtonTOFDCAxyNegData(0), hProtonTOFDCAxyPosData(0),
	pMIPVsEta(0), pPlateauVsEta(0), pMIPVsEtaV0s(0)


{

	for (int i_eta = 0; i_eta < nEta; ++i_eta)
	{
		hNsigmaPiPos[i_eta] = 0;
		hNsigmaKPos[i_eta] = 0;
		hNsigmaPPos[i_eta] = 0;
		hNsigmaPiNeg[i_eta] = 0;
		hNsigmaKNeg[i_eta] = 0;
		hNsigmaPNeg[i_eta] = 0;
		hPtTPCEtaPos[i_eta] = 0;
		hPtTPCEtaNeg[i_eta] = 0;

		hNsigmaTOFKPos[i_eta] = 0;
		hNsigmaTOFPPos[i_eta] = 0;
		hNsigmaTOFKNeg[i_eta] = 0;
		hNsigmaTOFPNeg[i_eta] = 0;
		hBetaPos[i_eta] = 0;
		hMomentumTOFEtaPos[i_eta] = 0;
		hPtTOFEtaPos[i_eta] = 0;
		hBetaNeg[i_eta] = 0;
		hMomentumTOFEtaNeg[i_eta] = 0;
		hPtTOFEtaNeg[i_eta] = 0;

		hdEdx[i_eta] = 0;
		hPtrTPC[i_eta] = 0;
		hPtVsP[i_eta] = 0;	

		nsigma_kaon_h[i_eta] = 0;
		random_cont_in_kaon_h[i_eta] = 0;
		nsigma_proton_h[i_eta] = 0;
		random_cont_in_proton_h[i_eta] = 0;
		nsigma_pion_h[i_eta] = 0;
		random_cont_in_pion_h[i_eta] = 0;

		histPiV0[i_eta] = 0;
		histPV0[i_eta] = 0;
		histEV0[i_eta] = 0;
		histPiTof[i_eta] = 0;

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
	fUseMC(kFALSE), fV0MMultiplicity(-1.0), fDetFlat("V0"), fV0MBin("0_1"), fIsMCclosure(kFALSE),
	fDeltaV0(kTRUE), fRemoveTrivialScaling(kFALSE), fnGen(-1), fPIDResponse(0x0),
	fTrackFilter(0x0), fTrackFilterPID(0x0), fOutputList(0), fEtaCut(0.8), fPtMin(0.5), fNcl(70), fdEdxCalibrated(kTRUE), 
	fSaveDCAxyHistograms(kFALSE), fEtaCalibrationPos(0x0), fEtaCalibrationNeg(0x0),
	fcutLow(0x0), fcutHigh(0x0), fcutDCAxy(0x0), fPeriod("16l"),
	ftrackmult08(0), fv0mpercentile(0), fFlat(-1), fFlatTPC(-1.), fFlatMC(-1),
	fMultSelection(0x0), hFlat(0),
	hFlatenicityMC(0), hFlatenicityMCRec(0), hFlatResponse(0), hFlatVsPtMC(0),
	hActivityV0McSect(0), hFlatVsNchMC(0),
	hMCPtPionPos(0),hMCPtKaonPos(0),hMCPtProtonPos(0),
	hMCPtPionNeg(0),hMCPtKaonNeg(0),hMCPtProtonNeg(0),
	hTPCRecTracksPionPos(0), hTPCRecTracksKaonPos(0), hTPCRecTracksProtonPos(0),
	hTPCRecTracksPionNeg(0), hTPCRecTracksKaonNeg(0), hTPCRecTracksProtonNeg(0),
	hTOFRecTracksPionPos(0), hTOFRecTracksKaonPos(0), hTOFRecTracksProtonPos(0),
	hTOFRecTracksPionNeg(0), hTOFRecTracksKaonNeg(0), hTOFRecTracksProtonNeg(0),
	hrTPCRecTracksPion(0), hrTPCRecTracksKaon(0), hrTPCRecTracksProton(0),
	hPionTPCDCAxyNegData(0), hPionTPCDCAxyPosData(0), hProtonTPCDCAxyNegData(0), hProtonTPCDCAxyPosData(0),
	hPionTOFDCAxyNegData(0), hPionTOFDCAxyPosData(0), hProtonTOFDCAxyNegData(0), hProtonTOFDCAxyPosData(0),
	pMIPVsEta(0), pPlateauVsEta(0), pMIPVsEtaV0s(0)

{
	for (int i_eta = 0; i_eta < nEta; ++i_eta)
	{
		hNsigmaPiPos[i_eta] = 0;
		hNsigmaKPos[i_eta] = 0;
		hNsigmaPPos[i_eta] = 0;
		hNsigmaPiNeg[i_eta] = 0;
		hNsigmaKNeg[i_eta] = 0;
		hNsigmaPNeg[i_eta] = 0;
		hPtTPCEtaPos[i_eta] = 0;
		hPtTPCEtaNeg[i_eta] = 0;

		hNsigmaTOFKPos[i_eta] = 0;
		hNsigmaTOFPPos[i_eta] = 0;
		hNsigmaTOFKNeg[i_eta] = 0;
		hNsigmaTOFPNeg[i_eta] = 0;
		hBetaPos[i_eta] = 0;
		hMomentumTOFEtaPos[i_eta] = 0;
		hPtTOFEtaPos[i_eta] = 0;
		hBetaNeg[i_eta] = 0;
		hMomentumTOFEtaNeg[i_eta] = 0;
		hPtTOFEtaNeg[i_eta] = 0;

		hdEdx[i_eta] = 0;
		hPtrTPC[i_eta] = 0;
		hPtVsP[i_eta] = 0;	

		nsigma_kaon_h[i_eta] = 0;
		random_cont_in_kaon_h[i_eta] = 0;
		nsigma_proton_h[i_eta] = 0;
		random_cont_in_proton_h[i_eta] = 0;
		nsigma_pion_h[i_eta] = 0;
		random_cont_in_pion_h[i_eta] = 0;

		histPiV0[i_eta] = 0;
		histPV0[i_eta] = 0;
		histEV0[i_eta] = 0;
		histPiTof[i_eta] = 0;
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

	const int nPtbins = 56;
	double Ptbins[nPtbins+1] = {
		0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
		0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
		1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
		3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0,
		9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0,
		20.0,22.0,24.0,26.0,30.0 };

	const int nPtbins_low_pT = 17;
	double Ptbins_low_pT[nPtbins_low_pT+1] = {
		0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 
		0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2 };

	const int nPtbins_intermediate_pT = 30;
	double Ptbins_intermediate_pT[nPtbins_intermediate_pT+1] = {
		0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 
		1.1,  1.2, 1.3,  1.4, 1.5,  1.6, 1.7,  1.8, 1.9, 
		2.0,  2.2, 2.4,  2.6, 2.8,  3.0, 3.2,  3.4, 3.6, 
		3.8, 4.0, 4.5, 5.0 };

	const int nPtbins_nSigmaTOF_pT = 21;
	double Ptbins_nSigmaTOF_pT[nPtbins_nSigmaTOF_pT+1] = {
		0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 
		1.1,  1.2, 1.3,  1.4, 1.5,  1.6, 1.7,  1.8, 1.9, 
		2.0,  2.2, 2.4,  2.6 };

	const int nPtbins_high_pT = 31;
	double Ptbins_high_pT[nPtbins_high_pT+1] = {
		2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 
		4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 
		11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
		26.0, 30.0 };

	const int nPtBinsV0s = 22;
	double ptBinsV0s[nPtBinsV0s+1] = { 
		0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
		1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
		10.0};

	const int nDeltaPiBins = 80;
	const double deltaPiLow  = 20;
	const double deltaPiHigh = 100;

	// create output objects
	float min_flat = -0.01;
	float max_flat = 1.01;
	int nbins_flat = 1020;
	if (fRemoveTrivialScaling) {
		min_flat = -0.1;
		max_flat = 9.9;
		nbins_flat = 2000;
	}

	const int nFlatbins = 1020;
	double Flatbins[nFlatbins+1] = {0.0};
	for (int i = 0; i <= nFlatbins; ++i) {
		Flatbins[i] = -0.01 + (double)i * 0.001;
	}

	const int nMultbins = 100;
	double Multbins[nMultbins+1] = {0.0};
	for (int i = 0; i <= nMultbins; ++i) {
		Multbins[i] = -0.5 + (double)i;
	}

	const int nnSigmabins = 80;
	double nSigmabins[nnSigmabins+1] = {0.0};

	for (int i = 0; i <= nnSigmabins; ++i){
		nSigmabins[i] = -4.0 + i*0.1;
	}

	const int nBetabins   = 151;
	double Betabins[nBetabins+1] = { 0.0 };
	for(int i = 0; i <= nBetabins; ++i){
		Betabins[i] = 0.6 + 0.003333 * ((double)i);
	}

	const int dEdxHigh = 100;
	const int dEdxLow = 40;

	const int ndEdxbins = dEdxHigh-dEdxLow;
	double dEdxbins[ndEdxbins+1];

	for(int i = dEdxLow; i <= dEdxHigh; ++i){
		dEdxbins[i-dEdxLow] = i;
	}

	const int nV0Multbins = 51;
	double V0Multbins[nV0Multbins+1] = { 0.0 };

	for (Int_t i = 0; i <= nV0Multbins; ++i){
		V0Multbins[i] = -0.5 + (double)i;
	}

	const int nV0Sectorsbins = 64;
	double V0Sectorsbins[nV0Sectorsbins+1] = { 0.0 };

	for (int i = 0; i <= nV0Sectorsbins; ++i){
		V0Sectorsbins[i] = -0.5 + (double)i;
	}

	OpenFile(1);
	fOutputList = new TList(); // this is a list which will contain all of your histograms
	fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

	hPionTPCDCAxyNegData = new TH2F("hPionTPCDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 600, -3.0, 3.0);	
	hPionTPCDCAxyPosData = new TH2F("hPionTPCDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 600, -3.0, 3.0);	
	hProtonTPCDCAxyNegData = new TH2F("hProtonTPCDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 600, -3.0, 3.0);	
	hProtonTPCDCAxyPosData = new TH2F("hProtonTPCDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 600, -3.0, 3.0);	
	hPionTOFDCAxyNegData = new TH2F("hPionTOFDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 600, -3.0, 3.0);	
	hPionTOFDCAxyPosData = new TH2F("hPionTOFDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 600, -3.0, 3.0);	
	hProtonTOFDCAxyNegData = new TH2F("hProtonTOFDCAxyNeg","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 600, -3.0, 3.0);	
	hProtonTOFDCAxyPosData = new TH2F("hProtonTOFDCAxyPos","; #it{p}_{T} (GeV/#it{c}); DCA_{xy}", nPtbins, Ptbins, 600, -3.0, 3.0);	

	pMIPVsEta = new TProfile("pMIPVsEta","; #eta; #LT dE/dx #GT_{MIP, primary tracks}", 50, -0.8, 0.8, 40.0, 60.0);
	pPlateauVsEta = new TProfile("pPlateauVsEta","; #eta; #LT dE/dx #GT_{Plateau, primary tracks}", 50, -0.8, 0.8, 60.0, 110.0);
	pMIPVsEtaV0s = new TProfile("pMIPVsEtaV0s","; #eta; #LT dE/dx #GT_{MIP, secondary tracks}", 50, -0.8, 0.8, 40.0, 60.0);
	hFlat = new TH1F(Form("hFlat_V0_%s",fV0MBin.c_str()), "; Flatenicity; Counts", nFlatbins, Flatbins);

	if (!fUseMC && !fSaveDCAxyHistograms) { 
		fOutputList->Add(hFlat);
		fOutputList->Add(pMIPVsEta);
		fOutputList->Add(pPlateauVsEta);
		fOutputList->Add(pMIPVsEtaV0s);

		if (fSaveDCAxyHistograms) {
			fOutputList->Add(hPionTPCDCAxyNegData);
			fOutputList->Add(hPionTPCDCAxyPosData);
			fOutputList->Add(hProtonTPCDCAxyNegData);
			fOutputList->Add(hProtonTPCDCAxyPosData);
			fOutputList->Add(hPionTOFDCAxyNegData);
			fOutputList->Add(hPionTOFDCAxyPosData);
			fOutputList->Add(hProtonTOFDCAxyNegData);
			fOutputList->Add(hProtonTOFDCAxyPosData);
		}
	}

	for (int i_eta = 0; i_eta < nEta; ++i_eta) 
	{
		hNsigmaPiPos[i_eta] = new TH3F(Form("hNsigmaPiPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins_low_pT, Ptbins_low_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hNsigmaKPos[i_eta] = new TH3F(Form("hNsigmaKPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins_low_pT, Ptbins_low_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hNsigmaPPos[i_eta] = new TH3F(Form("hNsigmaPPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins_low_pT, Ptbins_low_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hNsigmaPiNeg[i_eta] = new TH3F(Form("hNsigmaPiNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins_low_pT, Ptbins_low_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hNsigmaKNeg[i_eta] = new TH3F(Form("hNsigmaKNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins_low_pT, Ptbins_low_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hNsigmaPNeg[i_eta] = new TH3F(Form("hNsigmaPNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity", nPtbins_low_pT, Ptbins_low_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hPtTPCEtaPos[i_eta] = new TH2F(Form("hPtTPCEtaPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}; Flatenicity;)", nPtbins_low_pT, Ptbins_low_pT, nFlatbins, Flatbins);
		hPtTPCEtaNeg[i_eta] = new TH2F(Form("hPtTPCEtaNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); Flatenicity;)", nPtbins_low_pT, Ptbins_low_pT, nFlatbins, Flatbins);
		hNsigmaTOFKPos[i_eta] = new TH3F(Form("hNsigmaTOFKPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity",nPtbins_nSigmaTOF_pT,Ptbins_nSigmaTOF_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hNsigmaTOFPPos[i_eta] = new TH3F(Form("hNsigmaTOFPPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity",nPtbins_nSigmaTOF_pT,Ptbins_nSigmaTOF_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hNsigmaTOFKNeg[i_eta] = new TH3F(Form("hNsigmaTOFKNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity",nPtbins_nSigmaTOF_pT,Ptbins_nSigmaTOF_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hNsigmaTOFPNeg[i_eta] = new TH3F(Form("hNsigmaTOFPNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); n#sigma; Flatenicity",nPtbins_nSigmaTOF_pT,Ptbins_nSigmaTOF_pT, nnSigmabins, nSigmabins, nFlatbins, Flatbins);
		hBetaPos[i_eta] = new TH3F(Form("hBetaPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p} (GeV/#it{c}); #beta; Flatenicity",nPtbins_intermediate_pT,Ptbins_intermediate_pT,nBetabins,Betabins,nFlatbins,Flatbins);
		hBetaNeg[i_eta] = new TH3F(Form("hBetaNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p} (GeV/#it{c}); #beta; Flatenicity", nPtbins_intermediate_pT, Ptbins_intermediate_pT, nBetabins,Betabins, nFlatbins, Flatbins);
		hMomentumTOFEtaPos[i_eta] = new TH2F(Form("hMomentumTOFEtaPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p} (GeV/#it{c}); Flatenicity", nPtbins_intermediate_pT, Ptbins_intermediate_pT, nFlatbins, Flatbins);
		hMomentumTOFEtaNeg[i_eta] = new TH2F(Form("hMomentumTOFEtaNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p} (GeV/#it{c}); Flatenicity", nPtbins_intermediate_pT, Ptbins_intermediate_pT, nFlatbins, Flatbins);
		hPtTOFEtaPos[i_eta] = new TH2F(Form("hPtTOFEtaPos_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins_intermediate_pT, Ptbins_intermediate_pT, nFlatbins, Flatbins);
		hPtTOFEtaNeg[i_eta] = new TH2F(Form("hPtTOFEtaNeg_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins_intermediate_pT, Ptbins_intermediate_pT, nFlatbins, Flatbins);

		hdEdx[i_eta] = new TH3F(Form("hdEdx_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]), ";#it{p} (GeV/#it{c}); dE/dx; Flatenicity", nPtbins_high_pT, Ptbins_high_pT, ndEdxbins, dEdxbins, nFlatbins, Flatbins);
		hPtrTPC[i_eta] = new TH2F(Form("hPtrTPC_c%s_eta_%s",fV0MBin.c_str(),etaClass[i_eta]),";#it{p}_{T} (GeV/#it{c}); Flatenicity", nPtbins_high_pT, Ptbins_high_pT, nFlatbins, Flatbins);
		hPtVsP[i_eta] = new TH2F(Form("hPtVsP_eta_%s",etaClass[i_eta]), ";#it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c})", nPtbins, Ptbins, nPtbins, Ptbins);

		histPiV0[i_eta] = new TH2F(Form("hPiV0_%s",etaClass[i_eta]), "Pions id by V0; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histPV0[i_eta] = new TH2F(Form("hPV0_%s",etaClass[i_eta]), "Protons id by V0; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histPiTof[i_eta] = new TH2F(Form("hPiTOF_%s",etaClass[i_eta]), "Primary Pions from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histEV0[i_eta] = new TH2F(Form("hEV0_%s",etaClass[i_eta]), "Electrons id by V0; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);

		if (!fUseMC && !fSaveDCAxyHistograms) { 

			fOutputList->Add(hNsigmaPiPos[i_eta]);
			fOutputList->Add(hNsigmaKPos[i_eta]);
			fOutputList->Add(hNsigmaPPos[i_eta]);
			fOutputList->Add(hNsigmaPiNeg[i_eta]);
			fOutputList->Add(hNsigmaKNeg[i_eta]);
			fOutputList->Add(hNsigmaPNeg[i_eta]);
			fOutputList->Add(hPtTPCEtaPos[i_eta]);
			fOutputList->Add(hPtTPCEtaNeg[i_eta]);
			fOutputList->Add(hNsigmaTOFKPos[i_eta]);
			fOutputList->Add(hNsigmaTOFPPos[i_eta]);
			fOutputList->Add(hNsigmaTOFKNeg[i_eta]);
			fOutputList->Add(hNsigmaTOFPNeg[i_eta]);
			fOutputList->Add(hBetaPos[i_eta]);
			fOutputList->Add(hMomentumTOFEtaPos[i_eta]);
			fOutputList->Add(hPtTOFEtaPos[i_eta]);
			fOutputList->Add(hBetaNeg[i_eta]);
			fOutputList->Add(hMomentumTOFEtaNeg[i_eta]);
			fOutputList->Add(hPtTOFEtaNeg[i_eta]);
			fOutputList->Add(hdEdx[i_eta]);
			fOutputList->Add(hPtrTPC[i_eta]);
			fOutputList->Add(hPtVsP[i_eta]);
			fOutputList->Add(histPiV0[i_eta]);
			fOutputList->Add(histPV0[i_eta]);
			fOutputList->Add(histPiTof[i_eta]);
			fOutputList->Add(histEV0[i_eta]);
		}
	}

	if (fUseMC) {
		hFlatenicityMC = new TH2F("hFlatenicityMC", ";True Flatenicity; V0M Percentile;", nbins_flat, min_flat, max_flat, nCent, centClass );
		fOutputList->Add(hFlatenicityMC);

		hFlatenicityMCRec = new TH2F("hFlatenicityMCRec",";rec Flatenicity;V0M Percentile;", nbins_flat, min_flat, max_flat, nCent, centClass );
		fOutputList->Add(hFlatenicityMCRec);

		hFlatResponse = new TH2F("hFlatResponse", "; true flat; measured flat", nbins_flat, min_flat, max_flat, nbins_flat, min_flat, max_flat);
		fOutputList->Add(hFlatResponse);

		hFlatVsPtMC = new TH2F("hFlatVsPtMC", "MC true; Flatenicity; #it{p}_{T} (GeV/#it{c})", nbins_flat, min_flat, max_flat, nPtbins, Ptbins);
		fOutputList->Add(hFlatVsPtMC);

		hFlatVsNchMC = new TH2F("hFlatVsNchMC", "; true flat; true Nch", nbins_flat, min_flat, max_flat, 100, -0.5, 99.5);
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

		for (int i_eta = 0; i_eta < nEta; ++i_eta) 
		{	
			nsigma_kaon_h[i_eta] = new TH2F(Form("nsigma_kaon_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(nsigma_kaon_h[i_eta]);
			random_cont_in_kaon_h[i_eta] = new TH2F(Form("random_cont_in_kaon_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(random_cont_in_kaon_h[i_eta]);

			nsigma_proton_h[i_eta] = new TH2F(Form("nsigma_proton_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(nsigma_proton_h[i_eta]);
			random_cont_in_proton_h[i_eta] = new TH2F(Form("random_cont_in_proton_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(random_cont_in_proton_h[i_eta]);

			nsigma_pion_h[i_eta] = new TH2F(Form("nsigma_pion_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(nsigma_pion_h[i_eta]);
			random_cont_in_pion_h[i_eta] = new TH2F(Form("random_cont_in_pion_h_%s",etaClass[i_eta]),"; #it{p}_{T} (GeV/#it{c}); n#sigma",nPtbins,Ptbins,nnSigmabins,nSigmabins);
			fOutputList->Add(random_cont_in_pion_h[i_eta]);
		}
	}

	if (fUseMC) {
		hActivityV0McSect = new TProfile("hActivityV0McSect", "true; V0 sector; #LTmultiplicity#GT", 64, -0.5, 63.5);
		fOutputList->Add(hActivityV0McSect);
	}

	/* fEventCuts.AddQAplotsToList(fOutputList); */
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
		/* PostData(1, fOutputList); */
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
	float v0multalice = fMultSelection->GetEstimator("V0M")->GetValue();

	// Analyze V0s for the MB sample
	AnalyzeV0s();

	if (fV0MBin=="0_1"){
		if (!(fv0mpercentile >= 0.0 && fv0mpercentile < 1.0)) { return; }
	}
	else if (fV0MBin=="1_5"){
		if (!(fv0mpercentile >= 1.0 && fv0mpercentile < 5.0)) { return; }
	}
	else if (fV0MBin=="5_10"){
		if (!(fv0mpercentile >= 5.0 && fv0mpercentile < 10.0)) { return; }
	}
	else if (fV0MBin=="10_20"){
		if (!(fv0mpercentile >= 10.0 && fv0mpercentile < 20.0)) { return; }
	}
	else if (fV0MBin=="20_30"){
		if (!(fv0mpercentile >= 20.0 && fv0mpercentile < 30.0)) { return; }
	}
	else if (fV0MBin=="30_40"){
		if (!(fv0mpercentile >= 30.0 && fv0mpercentile < 40.0)) { return; }
	}
	else if (fV0MBin=="40_50"){
		if (!(fv0mpercentile >= 40.0 && fv0mpercentile < 50.0)) { return; }
	}
	else if (fV0MBin=="50_70"){
		if (!(fv0mpercentile >= 50.0 && fv0mpercentile < 70.0)) { return; }
	}
	else{
		if (!(fv0mpercentile >= 70.0 && fv0mpercentile < 100.0)) { return; }
	}

	// This condition is added to reject the very few low-multiplicity events V0M Percentile 70-100%
	// that have a very large V0M Amplitude
	if ((fV0MBin=="70_100") && (fv0mpercentile >= 70.0 && fv0mpercentile < 100.0)){	
		if (v0multalice >= 400.0) { return; }
	}

	/* fMidRapidityMult = GetMidRapidityMultiplicity(); */
	/* fFlatTPC = GetFlatenicityTPC(); */ 
	fFlat = GetFlatenicityV0();

	fFlatMC = -1;
	if (fUseMC) {
		if (!isGoodVtxPosMC) { return; }
		fFlatMC = GetFlatenicityMC();
		//	if (fFlatMC >= 0) {
		hFlatenicityMC->Fill(fFlatMC,fv0mpercentile);
		hFlatenicityMCRec->Fill(fFlat,fv0mpercentile);
		hFlatResponse->Fill(fFlatMC, fFlat);
		MakeMCanalysis();
		MakeMCanalysisPID();
		nSigmaContamination();
		//	}
	}

	if (fFlat >= 0.0) {

		hFlat->Fill(fFlat);
		// piKp as a function of Flattenicity
		MakePIDanalysis();
		// Charged particle spectra as a function of Flattenicity
		/* MakeDataanalysis(); */
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
		/* hFlatVsPtV0M[fV0Mindex]->Fill(fFlat, esdtrack->Pt()); */
	}
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::MakePIDanalysis() 
{
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

		float nSigmaPi = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion);
		float nSigmaK = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon);
		float nSigmaP = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
		float nSigmaPiTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion);
		float nSigmaKTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon);
		float nSigmaPTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton);

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
			hNsigmaPiNeg[nh]->Fill(pt,nSigmaPi,fFlat);
			hNsigmaKNeg[nh]->Fill(pt,nSigmaK,fFlat);
			hNsigmaPNeg[nh]->Fill(pt,nSigmaP,fFlat);
			hPtTPCEtaNeg[nh]->Fill(pt,fFlat);
		}else{
			hNsigmaPiPos[nh]->Fill(pt,nSigmaPi,fFlat);
			hNsigmaKPos[nh]->Fill(pt,nSigmaK,fFlat);
			hNsigmaPPos[nh]->Fill(pt,nSigmaP,fFlat);
			hPtTPCEtaPos[nh]->Fill(pt,fFlat);
		} 

		if( TOFPID(esdTrack) ){

			double trkLength = esdTrack->GetIntegratedLength();
			double beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);

			if(esdTrack->Charge() < 0.0){
				hNsigmaTOFKNeg[nh]->Fill(pt,nSigmaKTOF,fFlat);
				hNsigmaTOFPNeg[nh]->Fill(pt,nSigmaPTOF,fFlat);
				hBetaNeg[nh]->Fill(momentum,beta,fFlat);
				hMomentumTOFEtaNeg[nh]->Fill(momentum,fFlat);
				hPtTOFEtaNeg[nh]->Fill(pt,fFlat);
			}else{ 
				hNsigmaTOFKPos[nh]->Fill(pt,nSigmaKTOF,fFlat);
				hNsigmaTOFPPos[nh]->Fill(pt,nSigmaPTOF,fFlat);
				hBetaPos[nh]->Fill(momentum,beta,fFlat);
				hMomentumTOFEtaPos[nh]->Fill(momentum,fFlat);
				hPtTOFEtaPos[nh]->Fill(pt,fFlat);
			}
			if (TMath::Sqrt(nSigmaPi*nSigmaPi + nSigmaPiTOF*nSigmaPiTOF) < 2.0) { histPiTof[nh]->Fill(momentum,dEdx); }
		}

		hPtVsP[nh]->Fill(momentum,pt);

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), 1.0, fcutLow, fcutHigh))
			continue;

		if (Ncl < fNcl)
			continue;

		if (fdEdxCalibrated) dEdx *= 50.0/EtaCalibration(eta);

		hdEdx[nh]->Fill(momentum,dEdx,fFlat);
		hPtrTPC[nh]->Fill(pt,fFlat);

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
				/* hMIPVsEta->Fill(eta,dEdx); */
				pMIPVsEta->Fill(eta,dEdx);
			}
			if( dEdx > 70.0 && dEdx < 90.0 ){
				if(TMath::Abs(beta-1)<0.1){
					/* hPlateauVsEta->Fill(eta,dEdx); */
					pPlateauVsEta->Fill(eta,dEdx);
				}
			}
		}


	}//end of track loop


}
//________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::AnalyzeV0s()
{

	Int_t nv0s = fESD->GetNumberOfV0s();

	const AliESDVertex *myBestPrimaryVertex = fESD->GetPrimaryVertex();

	if ( !myBestPrimaryVertex ) 
		return;
	if ( !(myBestPrimaryVertex->GetStatus()) )
		return;

	Double_t  lPrimaryVtxPosition[3];
	myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);
	Double_t  lPrimaryVtxCov[6];
	myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);
	Double_t  lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();

	AliAODVertex* myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);

	//
	// LOOP OVER V0s, K0s, L, AL
	//

	for (Int_t iV0=0; iV0<nv0s; iV0++)
	{

		AliESDv0 *esdV0 = fESD->GetV0(iV0);
		if (!esdV0) continue;

		//check onfly status
		//		if( !esdV0->GetOnFlyStatus() )
		//			continue;

		if (esdV0->GetOnFlyStatus()!=0)
			continue;

		// AliESDTrack (V0 Daughters)
		UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());

		AliESDtrack *pTrack = fESD->GetTrack(lKeyPos);
		AliESDtrack *nTrack = fESD->GetTrack(lKeyNeg);

		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retreive one of the daughter track");
			continue;
		}

		// Remove like-sign
		if (pTrack->GetSign() == nTrack->GetSign())
			continue;

		// Eta cut on decay products
		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
			continue;

		if (!fTrackFilterPID->IsSelected(pTrack)) { continue; }
		if (!fTrackFilterPID->IsSelected(nTrack)) { continue; }

		// Check if switch does anything!
		Bool_t isSwitched = kFALSE;
		if (pTrack->GetSign() < 0) { // switch
			isSwitched = kTRUE;
			AliESDtrack* helpTrack = nTrack;
			nTrack = pTrack;
			pTrack = helpTrack;
		}

		AliKFVertex primaryVtxKF( *myPrimaryVertex );
		AliKFParticle::SetField(fESD->GetMagneticField());

		// Also implement switch here!!!!!!
		AliKFParticle* negEKF  = 0; // e-
		AliKFParticle* posEKF  = 0; // e+
		AliKFParticle* negPiKF = 0; // pi -
		AliKFParticle* posPiKF = 0; // pi +
		AliKFParticle* posPKF  = 0; // p
		AliKFParticle* negAPKF = 0; // p-bar

		if(!isSwitched) {
			negEKF  = new AliKFParticle( *(esdV0->GetParamN()) , 11);
			posEKF  = new AliKFParticle( *(esdV0->GetParamP()) ,-11);
			negPiKF = new AliKFParticle( *(esdV0->GetParamN()) ,-211);
			posPiKF = new AliKFParticle( *(esdV0->GetParamP()) , 211);
			posPKF  = new AliKFParticle( *(esdV0->GetParamP()) , 2212);
			negAPKF = new AliKFParticle( *(esdV0->GetParamN()) ,-2212);
		} else { // switch + and -
			negEKF  = new AliKFParticle( *(esdV0->GetParamP()) , 11);
			posEKF  = new AliKFParticle( *(esdV0->GetParamN()) ,-11);
			negPiKF = new AliKFParticle( *(esdV0->GetParamP()) ,-211);
			posPiKF = new AliKFParticle( *(esdV0->GetParamN()) , 211);
			posPKF  = new AliKFParticle( *(esdV0->GetParamN()) , 2212);
			negAPKF = new AliKFParticle( *(esdV0->GetParamP()) ,-2212);
		}

		AliKFParticle v0GKF;  // Gamma e.g. from pi0
		v0GKF+=(*negEKF);
		v0GKF+=(*posEKF);
		v0GKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0K0sKF; // K0 short
		v0K0sKF+=(*negPiKF);
		v0K0sKF+=(*posPiKF);
		v0K0sKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0LambdaKF; // Lambda
		v0LambdaKF+=(*negPiKF);
		v0LambdaKF+=(*posPKF);
		v0LambdaKF.SetProductionVertex(primaryVtxKF);

		AliKFParticle v0AntiLambdaKF; // Lambda-bar
		v0AntiLambdaKF+=(*posPiKF);
		v0AntiLambdaKF+=(*negAPKF);
		v0AntiLambdaKF.SetProductionVertex(primaryVtxKF);

		Double_t dmassG     = TMath::Abs(v0GKF.GetMass());
		Double_t dmassK     = TMath::Abs(v0K0sKF.GetMass()-0.498);
		Double_t dmassL     = TMath::Abs(v0LambdaKF.GetMass()-1.116);
		Double_t dmassAL    = TMath::Abs(v0AntiLambdaKF.GetMass()-1.116);

		if( dmassG  > 0.1 &&   
				dmassK  > 0.1 &&   
				dmassL  > 0.1 &&   
				dmassAL > 0.1
		  )    
			continue;

		for( Int_t case_v0 = 0; case_v0 < 2; ++case_v0 ){

			switch(case_v0){
				case 0:{

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if(dmassG < 0.1)
						       continue;

					       if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01){
						       continue;
					       }

					       if(dmassL<0.01){
						       fillPos = kTRUE;
					       }

					       if(dmassAL<0.01) {
						       if(fillPos)
							       continue;
						       fillNeg = kTRUE;
					       }

					       if(dmassK<0.01) {
						       if(fillPos||fillNeg)
							       continue;
						       fillPos = kTRUE;
						       fillNeg = kTRUE;
					       }

					       for(Int_t j = 0; j < 2; j++) {

						       AliESDtrack* track = 0;

						       if(j==0) {

							       if(fillNeg)
								       track = nTrack;
							       else
								       continue;
						       } else {

							       if(fillPos)
								       track = pTrack;
							       else
								       continue;
						       }

						       if(track->GetTPCsignalN()<fNcl) { continue; }
						       Double_t phi = track->Phi();

						       /* if(!PhiCut(track->Pt(), phi, track->Charge(), Magf, fcutLow, fcutHigh)) */
						       if(!PhiCut(track->Pt(), phi, track->Charge(), 1.0, fcutLow, fcutHigh)) { continue; }

						       Double_t eta      = track->Eta();
						       Double_t momentum = track->P();
						       Double_t dedx     = track->GetTPCsignal();
						       Double_t dedxUnc  = track->GetTPCsignal();

						       if(fdEdxCalibrated)
							       dedx *= 50.0/EtaCalibration(eta);

						       if(fillPos&&fillNeg){
							       if(dedxUnc < 60.0 && dedxUnc > 40.0){
								       if(momentum<0.6&&momentum>0.4){
									       pMIPVsEtaV0s->Fill(eta,dedx);
								       }
							       }
						       }

						       Int_t nh = -1;

						       if(TMath::Abs(eta)<0.2)
							       nh = 0; 
						       else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
							       nh = 1; 
						       else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
							       nh = 2; 
						       else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
							       nh = 3; 

						       if(nh<0) { continue; }

						       if(fillPos&&fillNeg){
							       float nSigmaPiTPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
							       float nSigmaPiTOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
							       if ((momentum < 1.0) && (TMath::Abs(nSigmaPiTPC) < 3.0)) { histPiV0[nh]->Fill(momentum, dedx); }
							       if (!TOFPID(track)) { continue; }
							       if ((momentum >= 1.0 ) && (TMath::Abs(nSigmaPiTOF) < 3.0)) { histPiV0[nh]->Fill(momentum, dedx); }
						       }
						       else{
							       float nSigmaPTPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);
							       float nSigmaPTOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);
							       if ((momentum < 1.0) && (TMath::Abs(nSigmaPTPC) < 3.0)) { histPV0[nh]->Fill(momentum, dedx); }
							       if (!TOFPID(track)) { continue; }
							       if ((momentum >= 1.0) && (TMath::Abs(nSigmaPTOF) < 3.0)) { histPV0[nh]->Fill(momentum, dedx); }
						       }
					       }//end loop over two tracks

				       };
				       break;

				case 1:{//gammas

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if( dmassK>0.01 && dmassL>0.01 && dmassAL>0.01 ) {
						       if( dmassG<0.01 && dmassG>0.0001 ) {
							       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationEl(nTrack->Eta())) < 5.0)
								       fillPos = kTRUE;
						       } else {
							       continue;
						       }
					       }

					       if(fillPos == kTRUE && fillNeg == kTRUE)
						       continue;

					       AliESDtrack* track = 0;
					       if(fillNeg)
						       track = nTrack;
					       else if(fillPos)
						       track = pTrack;
					       else
						       continue;

					       Double_t dedx     = track->GetTPCsignal();
					       Double_t eta      = track->Eta();
					       Double_t phi      = track->Phi();
					       Double_t momentum = track->P();

					       if(fdEdxCalibrated)
						       dedx *= 50.0/EtaCalibration(eta);

					       if(track->GetTPCsignalN()<fNcl) { continue; }

					       if(!PhiCut(track->Pt(), phi, track->Charge(), 1.0, fcutLow, fcutHigh)) { continue; }

					       Int_t nh = -1;

					       if(TMath::Abs(eta)<0.2)
						       nh = 0; 
					       else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
						       nh = 1; 
					       else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
						       nh = 2; 
					       else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
						       nh = 3; 

					       if(nh<0) { continue; }

					       histEV0[nh]->Fill(momentum, dedx);

				       };
				       break;


			}//end switch

		}//end loop over case V0

		// clean up loop over v0

		delete negPiKF;
		delete posPiKF;
		delete posPKF;
		delete negAPKF;



	}

	delete myPrimaryVertex;

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
		/* hPtPrimIn->Fill(particle->Pt()); */
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
		/* hPtOut->Fill(esdtrack->Pt()); */
		Int_t mcLabel = -1;
		mcLabel = TMath::Abs(esdtrack->GetLabel());
		if (fMC->IsPhysicalPrimary(mcLabel)) {
			/* hPtPrimOut->Fill(esdtrack->Pt()); */
		} else {
			/* hPtSecOut->Fill(esdtrack->Pt()); */
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
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
			continue;

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

	const int nRings2 = 8;
	const int nSectors2 = 8;
	const int nCells2 = nRings2 * nSectors2;
	//	float maxEta2[nRings2] = {-0.4, 0.0, +0.4, +0.8};
	//	float minEta2[nRings2] = {-0.8, -0.4, +0.0, +0.4};
	float maxEta2[nRings2] = {-0.6, -0.4, -0.2, +0.0, +0.2, +0.4, +0.6, +0.8};
	float minEta2[nRings2] = {-0.8, -0.6, -0.4, -0.2, +0.0, +0.2, +0.4, +0.6};
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

	// only events with tracks in the central region
	if (mRho_glob <= 0.0) { return 9999; }

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
	float RhoLattice[nCells];
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		RhoLattice[iCh] = 0.0;
	}

	// before calibration
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		Float_t mult = lVV0->GetMultiplicity(iCh);
		RhoLattice[iCh] = mult;
		/* hActivityV0DataSectBefore->Fill(iCh, RhoLattice[iCh]); */
	}
	// after calibration
	/* if (fIsCalib) { */
	/* 	for (Int_t iCh = 0; iCh < nCells; iCh++) { */
	/* 		RhoLattice[iCh] *= fParVtx->Eval(0.0) / fParVtx->Eval(fVtxz); */
	/* 	} */
	/* } */

	// Filling histos with mult info
	/* float total_v0_tmp = 0.0; */
	/* for (Int_t iCh = 0; iCh < nCells; iCh++) { */
	/* 	hActivityV0DataSect->Fill(iCh, RhoLattice[iCh]); */
	/* total_v0_tmp += RhoLattice[iCh]; */
	/* } */
	/* float total_v0 = total_v0_tmp; */
	/* cout << "total_v0 = " << total_v0 << endl; */
	/* hV0vsVtxz->Fill(fVtxz, total_v0); */

	Int_t nringA = 0;
	Int_t nringC = 0;
	for (Int_t iCh = 0; iCh < nCells; iCh++) {
		Float_t detaV0 = -1;
		// Float_t mult = lVV0->GetMultiplicity(iCh);
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
		// consider the different eta coverage
		RhoLattice[iCh] /= detaV0;
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
Double_t AliAnalysisTaskFlatenicityPiKp::GetMidRapidityMultiplicity() {

	int mult_glob = 0;
	int nTracks = fESD->GetNumberOfTracks();
	for (Int_t iT = 0; iT < nTracks; ++iT) 
	{

		AliESDtrack *esdtrack = static_cast<AliESDtrack *>(
				fESD->GetTrack(iT)); // get a track (type AliesdTrack)
		if (!esdtrack)
			continue;
		if (!fTrackFilter->IsSelected(esdtrack))
			continue;
		float eta_a = esdtrack->Eta();

		if (TMath::Abs(eta_a) > fEtaCut)
			continue;
		if (esdtrack->Pt() < fPtMin)
			continue;

		mult_glob++;
	}

	return (double)mult_glob;
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
		if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMC))
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
			float mult = RhoLattice[i_segment];
			RhoLattice[i_segment] /= deltaEta;
			// Filling histos with mult info
			if (fDeltaV0) { mult /= deltaEta; }
			else { mult *= 1.0; }

			hActivityV0McSect->Fill(i_segment,mult);
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
//______________________________________________________________________________
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
//______________________________________________________________________________
Bool_t AliAnalysisTaskFlatenicityPiKp::TOFPID(AliESDtrack * track) 
{
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
//______________________________________________________________________________
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
	} else if(fPeriod == "16d" || fPeriod == "16e" || fPeriod == "16g" || fPeriod == "16h" || fPeriod == "16i" || fPeriod == "16j" || fPeriod == "16o" || fPeriod == "16p"){
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
//________________________________________________________________________
double AliAnalysisTaskFlatenicityPiKp::EtaCalibrationEl(const double &eta)
{

	double aPos = 0.0;
	double bPos = 0.0;
	double cPos = 0.0;
	double dPos = 0.0;
	double ePos = 0.0;

	double aNeg = 0.0;
	double bNeg = 0.0;
	double cNeg = 0.0;
	double dNeg = 0.0;
	double eNeg = 0.0;

	if(fPeriod=="16l"){
		aPos = 79.4195; bPos = 7.82459; cPos = -23.3466; dPos = 26.5577; ePos = -8.27151;
		aNeg = 79.8571; bNeg = -14.2921; cNeg = -66.6972; dNeg = -103.794; eNeg = -50.5771;
	}else if(fPeriod=="16k"){
		aPos = 80.254; bPos = 6.37076; cPos = -50.9878; dPos = 116.611; ePos = -79.0483;
		aNeg = 79.8728; bNeg = -3.08265; cNeg = -11.3778; dNeg = -20.6605; eNeg = -12.3861;
	}else if(fPeriod=="16deghijop"){
		aPos = 80.0719; bPos = 7.10053; cPos = -42.4788; dPos = 86.1074; ePos = -54.0891;
		aNeg = 79.6155; bNeg = -12.1254; cNeg = -66.2488; dNeg = -132.426; eNeg = -85.0155;
	}else if(fPeriod=="17data"){
		aPos = 82.4621; bPos = 5.20353; cPos = -32.2608; dPos = 63.4788; ePos = -39.3277;
		aNeg = 82.306; bNeg = -4.04076; cNeg = -22.133; dNeg = -40.5782; eNeg = -23.8157;
	}else{
		aPos = 79.7726; bPos = 6.83744; cPos = -40.0469; dPos = 78.987; ePos = -50.1373;
		aNeg = 79.4863; bNeg = -5.00403; cNeg = -21.6184;  dNeg = -39.1295; eNeg = -24.8757;
	}

	TF1* felededxfitPos = new TF1("felededxfitPos", "pol4", 0.0, 1.0);
	TF1* felededxfitNeg = new TF1("felededxfitNeg", "pol4", -1.0, 0.0);
	double dedx_electrons = 999;

	for(int i=0; i<5; ++i){
		felededxfitNeg->SetParameter(i,0.0);
		felededxfitPos->SetParameter(i,0.0);
	}

	if(eta<0.0){
		felededxfitNeg->SetParameter(0,aNeg);
		felededxfitNeg->SetParameter(1,bNeg);
		felededxfitNeg->SetParameter(2,cNeg);
		felededxfitNeg->SetParameter(3,dNeg);
		felededxfitNeg->SetParameter(4,eNeg);

		dedx_electrons = felededxfitNeg->Eval(eta);
	}
	else{
		felededxfitPos->SetParameter(0,aPos);
		felededxfitPos->SetParameter(1,bPos);
		felededxfitPos->SetParameter(2,cPos);
		felededxfitPos->SetParameter(3,dPos);
		felededxfitPos->SetParameter(4,ePos);

		dedx_electrons = felededxfitPos->Eval(eta);
	}

	delete felededxfitPos;
	delete felededxfitNeg;

	return dedx_electrons;

}
//______________________________________________________________________________
void AliAnalysisTaskFlatenicityPiKp::nSigmaContamination() {

	Int_t iTracks(fESD->GetNumberOfTracks());          

	for (Int_t iT = 0; iT < iTracks; ++iT) 
	{

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

		/* if (!fMC->IsPhysicalPrimary(mcLabel) ) */
		/* 	continue; */

		AliMCParticle *mc_particle = nullptr;
		mc_particle = (AliMCParticle*)fMC->GetTrack(mcLabel);

		if (!mc_particle) 
			continue;

		if (TMath::Abs(esdTrack->Charge())==0)
			continue;

		Int_t pdgCode = mc_particle->PdgCode();
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

		Float_t nsigma_pi = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion);
		Float_t nsigma_k = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon);
		Float_t nsigma_p = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
		Double_t pt = esdTrack->Pt();

		//	if (TMath::Abs(nsigma_k) <= 4.0){
		if (pidCode==2) { nsigma_kaon_h[nh]->Fill(pt,nsigma_k); }
		if (!(pidCode==2)) { random_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); }
		/* else if (pidCode==1) { pion_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); } */
		/* else if (pidCode==7) { electron_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); } */
		/* else { random_cont_in_kaon_h[nh]->Fill(pt,nsigma_k); } */
		//	}

		//	if (TMath::Abs(nsigma_p) <= 4.0){	
		if (pidCode==3) { nsigma_proton_h[nh]->Fill(pt,nsigma_p); }
		if (!(pidCode==3)) { random_cont_in_proton_h[nh]->Fill(pt,nsigma_p); }
		/* else if (pidCode==1) { pion_cont_in_proton_h[nh]->Fill(pt,nsigma_p); } */
		/* else if (pidCode==7) { electron_cont_in_proton_h[nh]->Fill(pt,nsigma_p); } */
		/* else { random_cont_in_proton_h[nh]->Fill(pt,nsigma_p); } */
		/* //	} */

		//	if (TMath::Abs(nsigma_pi) <= 3.0){	
		if (pidCode==1) { nsigma_pion_h[nh]->Fill(pt,nsigma_pi); }
		if (!(pidCode==1)) { random_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); }
		/* else if (pidCode==2) { kaon_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); } */
		/* else if (pidCode==7) { electron_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); } */
		/* else { random_cont_in_pion_h[nh]->Fill(pt,nsigma_pi); } */
		/* //	} */
	} 

}
//______________________________________________________________________________
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
//______________________________________________________________________________
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

