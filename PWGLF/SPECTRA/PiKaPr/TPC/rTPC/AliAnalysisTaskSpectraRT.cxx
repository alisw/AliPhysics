/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
 *************************************************************************/

#include "AliAnalysisTaskSpectraRT.h"

// ROOT includes
#include <TList.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TProfile.h>
#include <TParticle.h>
#include <TFile.h>

// AliRoot includes
#include "AliAnalysisTask.h"
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliPIDResponse.h>
#include "AliTPCPIDResponse.h"
#include "AliAnalysisUtils.h"

#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include "AliMCParticle.h"
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include <AliCentrality.h>
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h>
#include <AliVParticle.h>
#include <AliPID.h>
#include <AliAODPid.h>
#include <AliAODMCHeader.h>

#include <iostream>
#include <vector>
////class AliAnalysisTaskSpectraRT;    // your analysis class
using namespace std;


//
// Responsible:
// Antonio Ortiz (Lund)
// Peter Christiansen (Lund)
// Omar Vazquez (Lund)


static float Magf                = 1;
static const int nPid              = 4;
static const int nRegion           = 4;
static const int nHists            = 4;
static const int nRt               = 5;
const char* Region[4] = {"Toward","Away","Transverse","FullAzimuth"};
const char* Pid[nPid] = {"Charged","Pion","Kaon","Proton"};
const char* Charge[2] = {"Pos","Neg"};
static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps

ClassImp(AliAnalysisTaskSpectraRT)
	AliAnalysisTaskSpectraRT::AliAnalysisTaskSpectraRT():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fEventCuts(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fMCArray(0x0),
		fPIDResponse(0x0),
		fTrackFilterGolden(0x0),
		fTrackFilter(0x0),
		fHybridTrackCuts1(0x0),
		fHybridTrackCuts2(0x0),
		utils(0x0),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fIsMCclosure(kTRUE),
		fRandom(0x0),
		fNcl(70),
		fEtaCut(0.9),
		fdEdxCalibrated(kTRUE),
		fDeDxMIPMin(40),
		fDeDxMIPMax(60),
		fdEdxHigh(200),
		fdEdxLow(40),
		fPeriod("l"),
		fSetTPConlyTrkCuts(kFALSE),
		fSelectHybridTracks(kTRUE),
		fLeadPtCutMin(5.0),
		fLeadPtCutMax(40.0),
		fGenLeadPhi(0.0),
		fGenLeadPt(0.0),
		fGenLeadIn(0.0),
		fRecLeadPhi(0.0),
		fRecLeadPt(0.0),
		fRecLeadIn(0.0),
		fPtMin(0.15),
		fListOfObjects(0),
		fEvents(0x0),
		hNchTSData(0x0),
		hPhiTotal(0x0),
		hPhiStandard(0x0),
		hPhiHybrid1(0x0),
		fEtaCalibration(0x0),
		fEtaCalibrationEl(0x0),
		fcutDCAxy(0x0),
		fcutLow(0x0),
		fcutHigh(0x0),
		hMIPVsEta(0x0),
		pMIPVsEta(0x0),
		hPlateauVsEta(0x0),
		pPlateauVsEta(0x0),
		hMIPVsEtaV0s(0x0),
		pMIPVsEtaV0s(0x0)

{


	for(int r = 0; r < nRegion; r++){

		if(r < (nRegion-1)){
			hPhiData[r] = 0;
			hNchVsPtPosTPC[r] = 0;
			hNchVsPtNegTPC[r] = 0;
			hNchVsPtPosTOF[r] = 0;
			hNchVsPtNegTOF[r] = 0;
		}

		for(int j = 0; j < nHists; j++){
			hNchVsPtDataPosTOF[r][j] = 0;
			hNchVsPtDataNegTOF[r][j] = 0;
			hDeDxVsP[r][j] = 0;

		}	// ending

	}	// region        


	for(int r = 0; r < (nRegion-1); r++){

		for(int j = 0; j < nHists; j++){

			hNchVsPtDataPosPionTPC[r][j] = 0;
			hNchVsPtDataNegPionTPC[r][j] = 0;
			hNchVsPtDataPosKaonTPC[r][j] = 0;
			hNchVsPtDataNegKaonTPC[r][j] = 0;
			hNchVsPtDataPosProtonTPC[r][j] = 0;
			hNchVsPtDataNegProtonTPC[r][j] = 0;
		}	// ending

	}	// region        


	for(int j = 0; j < nHists; ++j){
		histPiTof[j] = 0;
		histEV0[j] = 0;
		histPV0[j] = 0;
		histPiV0[j] = 0;
		//	   hMIPVsPhi[j] = 0;
		//	   pMIPVsPhi[j] = 0;
		//	   hPlateauVsPhi[j] = 0;
		//	   pPlateauVsPhi[j] = 0;
		hPtVsP[j] = 0;

	}

}


AliAnalysisTaskSpectraRT::AliAnalysisTaskSpectraRT(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fEventCuts(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fMCArray(0x0),
	fPIDResponse(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilter(0x0),
	fHybridTrackCuts1(0x0),
	fHybridTrackCuts2(0x0),
	utils(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fIsMCclosure(kTRUE),
	fRandom(0x0),
	fNcl(70),
	fEtaCut(0.9),
	fdEdxCalibrated(kTRUE),
	fDeDxMIPMin(40),
	fDeDxMIPMax(60),
	fdEdxHigh(200),
	fdEdxLow(40),
	fPeriod("l"),
	fSetTPConlyTrkCuts(kFALSE),
	fSelectHybridTracks(kTRUE),
	fLeadPtCutMin(5.0),
	fLeadPtCutMax(40.0),
	fGenLeadPhi(0.0),
	fGenLeadPt(0.0),
	fGenLeadIn(0.0),
	fRecLeadPhi(0.0),
	fRecLeadPt(0.0),
	fRecLeadIn(0.0),
	fPtMin(0.15),
	fListOfObjects(0),
	fEvents(0x0),
	hNchTSData(0x0),
	hPhiTotal(0x0),
	hPhiStandard(0x0),
	hPhiHybrid1(0x0),
	fEtaCalibration(0x0),
	fEtaCalibrationEl(0x0),
	fcutDCAxy(0x0),
	fcutLow(0x0),
	fcutHigh(0x0),
	hMIPVsEta(0x0),
	pMIPVsEta(0x0),
	hPlateauVsEta(0x0),
	pPlateauVsEta(0x0),
	hMIPVsEtaV0s(0x0),
	pMIPVsEtaV0s(0x0)

{

	for(int r = 0; r < nRegion; r++){

		if(r < (nRegion-1)){
			hPhiData[r] = 0;
			hNchVsPtPosTPC[r] = 0;
			hNchVsPtNegTPC[r] = 0;
			hNchVsPtPosTOF[r] = 0;
			hNchVsPtNegTOF[r] = 0;
		}

		for(int j = 0; j < nHists; j++){

			hNchVsPtDataPosTOF[r][j] = 0;
			hNchVsPtDataNegTOF[r][j] = 0;
			hDeDxVsP[r][j] = 0 ;
		}

	}        

	for(int r = 0; r < (nRegion-1); r++){

		for(int j = 0; j < nHists; j++){

			hNchVsPtDataPosPionTPC[r][j] = 0;
			hNchVsPtDataNegPionTPC[r][j] = 0;
			hNchVsPtDataPosKaonTPC[r][j] = 0;
			hNchVsPtDataNegKaonTPC[r][j] = 0;
			hNchVsPtDataPosProtonTPC[r][j] = 0;
			hNchVsPtDataNegProtonTPC[r][j] = 0;
		}	// ending

	}	// region        


	for(int j = 0; j < nHists; ++j){
		histPiTof[j] = 0;
		histEV0[j] = 0;
		histPV0[j] = 0;
		histPiV0[j] = 0;
		//	   hMIPVsPhi[j] = 0;
		//	   pMIPVsPhi[j] = 0;
		//	   hPlateauVsPhi[j] = 0;
		//	   pPlateauVsPhi[j] = 0;
		hPtVsP[j] = 0;

	}


	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());//esto es nuevo
}

AliAnalysisTaskSpectraRT::~AliAnalysisTaskSpectraRT() {
	//
	// Destructor
	//

}
//______________________________________________________________________________
void AliAnalysisTaskSpectraRT::UserCreateOutputObjects()
{
	// This method is called once per worker node
	// Here we define the output: histograms and debug tree if requested
	// We also create the random generator here so it might get different seeds...
	fRandom = new TRandom(0); // 0 means random seed


	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if(man){
		AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
		if(inputHandler)fPIDResponse = inputHandler->GetPIDResponse();
	}

	//	fCuts *** leading particle ***
	if(!fTrackFilterGolden){
		fTrackFilterGolden = new AliAnalysisFilter("trackFilter2011");
		AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
		fTrackFilterGolden->AddCuts(esdTrackCutsGolden);
	}

	//	Track Cuts for Nch in the Transverse Side
	if(!fTrackFilter){
		fTrackFilter = new AliAnalysisFilter("trackFilterTPCOnly");
		SetTrackCuts(fTrackFilter);
	}

	if(!fHybridTrackCuts1){
		fHybridTrackCuts1 = new AliESDtrackCuts("fHybridTrackCuts1");	

		// TPC
		//		if(clusterCut == 0)  esdTrackCuts->SetMinNClustersTPC(50);
		//		else if (clusterCut == 1) {
		fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
		fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		//		}
		//		else {
		//			AliWarningClass(Form("Wrong value of the clusterCut parameter (%d), using cut on Nclusters",clusterCut));
		//			esdTrackCuts->SetMinNClustersTPC(50);
		//		}
		fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
		fHybridTrackCuts1->SetAcceptKinkDaughters(kFALSE);
		fHybridTrackCuts1->SetRequireTPCRefit(kTRUE);
		// ITS
		fHybridTrackCuts1->SetRequireITSRefit(kFALSE);
		fHybridTrackCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
				AliESDtrackCuts::kNone);
		//		if(selPrimaries) {
		// 7*(0.0015+0.0050/pt^1.1)
		fHybridTrackCuts1->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
		fHybridTrackCuts1->SetMaxChi2TPCConstrainedGlobal(36);
		//		}
		fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
		fHybridTrackCuts1->SetDCAToVertex2D(kFALSE);
		fHybridTrackCuts1->SetRequireSigmaToVertex(kFALSE);

		fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);

	} 

	if(!fHybridTrackCuts2){
		fHybridTrackCuts2 = new AliESDtrackCuts("fHybridTrackCuts2");	
		///	fHybridTrackCuts2->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
		fHybridTrackCuts2->SetRequireITSRefit(kFALSE);
	} 



	//OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	//
	// Histograms
	//


	fEvents = new TH1F( "fEvents", "; Evt. Sel.",12,0,12);
	fEvents->GetXaxis()->SetBinLabel(1, "Processed");
	fEvents->GetXaxis()->SetBinLabel(2, "PhysSel+Trigger");
	fEvents->GetXaxis()->SetBinLabel(3, "INEL>0");
	fEvents->GetXaxis()->SetBinLabel(4, "BG");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(5, "IsPileUpFromSPDinMultBins");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(6, "Incom DAQ");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(7, "Res&Proximity");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(8, "|Vtz|<10cm");//NotinVertexcut");
	fListOfObjects->Add(fEvents);

	const int nPtBinsV0s = 21;
	double ptBinsV0s[nPtBinsV0s+1] = { 
		0.1, 0.2, 0.3, 0.4,
		0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 
		1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 9.0};

	const int nPtBins = 25;
	double ptBins[nPtBins+1] = {
		0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
		0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 
		3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0};

	/*
	   const int nPtBins = 54;
	   double ptBins[nPtBins+1] = {
	   0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35,
	   0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
	   0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
	   1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40,
	   3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00,
	   9.00, 10.00, 12.00, 15.00, 20.00};
	   */

	const int nBinsNsigma = 50;
	double binsNsigma[nBinsNsigma+1] = {0};

	for(int i = 0; i <= nBinsNsigma; ++i){
		binsNsigma[i] = -10.0+i*0.4;
	}

	const int nDeltaPiBins   = 55;
	double DeltaPiBins[nDeltaPiBins+1] = { 0 };
	for(int i = 0; i <= nDeltaPiBins; ++i){
		DeltaPiBins[i] = 20.0+i*1.0;
	}

	const int ndEdxBins   = 100;
	double dEdxBins[ndEdxBins+1] = { 0 };
	for(int i = 0; i <= ndEdxBins; ++i){
		dEdxBins[i] = fdEdxLow+i*1.0;
	}

	const int nBetaBins   = 100;
	double BetaBins[nBetaBins+1] = { 0 };
	for(int i = 0; i <= nBetaBins; ++i){
		BetaBins[i] = 0.2+((double)i)/100.0;
	}


	/*
	   const Int_t ndcaBins = 100;
	   Double_t dcaBins[ndcaBins+1] = {
	   -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1,
	   -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1,
	   -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1,
	   -1.0, -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6,
	   -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15,
	   -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
	   0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
	   1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2,
	   2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4,
	   3.5, 3.6, 3.7, 3.8, 3.9, 4.0};
	   */ 

	/*
h: fMeanChT = 7.269
i: fMeanChT = 7.257
j: fMeanChT = 7.265
k: fMeanChT = 7.261
l: fMeanChT = 7.266
o: fMeanChT = 7.211
p: fMeanChT = 7.216


*/

	const int nBinsRT = 50;
	double binsRT[nBinsRT+1] = {0};

	for(int i = 0; i <= nBinsRT; ++i)
		binsRT[i] = (((double)i)-0.5)/1.0;

	const char* ending[nHists] = {"02", "24", "46", "68"};
	//const char* ending2[nHists] = {"02", "68"};

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+TMath::Pi()/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+TMath::Pi()/18.0+0.035", 0, 50);

	fEtaCalibration   = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	fEtaCalibrationEl = new TF1("fDeDxVsEtaEl", "pol4", 0.0, 1.0);

	hNchTSData = new TH1F("hMultTSData",";#it{N}_{acc}^{TS}; Entries",nBinsRT,binsRT);
	fListOfObjects->Add(hNchTSData);

	hPhiTotal = new TH2F("hPhiSum","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiTotal);

	hPhiStandard = new TH2F("hPhiStandard","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiStandard);

	hPhiHybrid1 = new TH2F("hPhiHybrid1","; #eta; #varphi",50,-0.8,0.8,100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0);
	fListOfObjects->Add(hPhiHybrid1);

	// Histos rTPC

	hMIPVsEta = new TH2F("hMIPVsEta","; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin,fDeDxMIPMin,fDeDxMIPMax);
	pMIPVsEta = new TProfile("pMIPVsEta","; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMin,fDeDxMIPMax);

	hPlateauVsEta = new TH2F("hPlateauVsEta","; #eta; dE/dx_{Plateau, primary tracks}",50,-0.8,0.8,50, 60, 110);
	pPlateauVsEta = new TProfile("pPlateauVsEta","; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",50,-0.8,0.8, 60, 110);

	hMIPVsEtaV0s = new TH2F("hMIPVsEtaV0s","; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin,fDeDxMIPMin,fDeDxMIPMax);
	pMIPVsEtaV0s = new TProfile("pMIPVsEtaV0s","; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMin,fDeDxMIPMax);

	fListOfObjects->Add(hMIPVsEta);
	fListOfObjects->Add(pMIPVsEta);
	fListOfObjects->Add(hPlateauVsEta);
	fListOfObjects->Add(pPlateauVsEta);
	fListOfObjects->Add(hMIPVsEtaV0s);
	fListOfObjects->Add(pMIPVsEtaV0s);

	for( int j = 0; j < nHists; j++ ){

		histEV0[j] = new TH2F(Form("histEV0_%s",ending[j]),"Electronsons from V0s; #it{p} (GeV/#it{c}); d#it{e}d#it{x}",nPtBinsV0s,ptBinsV0s,nDeltaPiBins,DeltaPiBins);
		histPV0[j] = new TH2F(Form("histPV0_%s",ending[j]),"Electronsons from V0s; #it{p} (GeV/#it{c}); d#it{e}d#it{x}",nPtBinsV0s,ptBinsV0s,nDeltaPiBins,DeltaPiBins);
		histPiV0[j] = new TH2F(Form("histPiV0_%s",ending[j]),"Electronsons from V0s; #it{p} (GeV/#it{c}); d#it{e}d#it{x}",nPtBinsV0s,ptBinsV0s,nDeltaPiBins,DeltaPiBins);
		histPiTof[j] = new TH2F(Form("hPiTOF_%s",ending[j]),"Pions from TOF;#it{p} (GeV/#it{c});d#it{e}d#it{x}",nPtBinsV0s,ptBinsV0s,nDeltaPiBins,DeltaPiBins);

		/*
		   hMIPVsPhi[j] = new TH2D(Form("hMIPVsPhi_%s",ending[j]),";#phi (rad); dE/dx MIP",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0,fDeDxMIPMax-fDeDxMIPMin,fDeDxMIPMin,fDeDxMIPMax);
		   hMIPVsPhi[j]->Sumw2();

		   pMIPVsPhi[j] = new TProfile(Form("pMIPVsPhi_%s",ending[j]),";#phi (rad); dE/dx MIP",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0,fDeDxMIPMin, fDeDxMIPMax);
		   pMIPVsPhi[j]->Sumw2();

		   hPlateauVsPhi[j] = new TH2D(Form("hPlateauVsPhi_%s",ending[j]),";#phi (rad); dE/dx Plateau",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0,20, 70, 90);
		   hPlateauVsPhi[j]->Sumw2();

		   pPlateauVsPhi[j] = new TProfile(Form("pPlateauVsPhi_%s",ending[j]),";#phi (rad); dE/dx Plateau",100,-TMath::Pi()/2.0,5.0*TMath::Pi()/2.0,fDeDxMIPMax, 95);
		   pPlateauVsPhi[j]->Sumw2();
		   */
		if(!fAnalysisMC){
			fListOfObjects->Add(histEV0[j]);
			fListOfObjects->Add(histPV0[j]);
			fListOfObjects->Add(histPiV0[j]);
			fListOfObjects->Add(histPiTof[j]);
			//			fListOfObjects->Add(hMIPVsPhi[j]);
			//			fListOfObjects->Add(pMIPVsPhi[j]);
			//			fListOfObjects->Add(hPlateauVsPhi[j]);
			//			fListOfObjects->Add(pPlateauVsPhi[j]);
		}
	}	// Only ending


	for( int j = 0; j < nHists; j++ ){

		hPtVsP[j] = new TH2F(Form("hPtVsP_%s",ending[j]),";#it{p} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c})",nPtBins,ptBins,nPtBins,ptBins);
		hPtVsP[j]->Sumw2();
		fListOfObjects->Add(hPtVsP[j]);

	}	// Only ending


	for(int r = 0; r < nRegion; ++r){

		if(r<(nRegion-1)){

			hPhiData[r] = new TH1F(Form("hPhiData_%s",Region[r]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);

			hNchVsPtPosTPC[r] = new TH2F(Form("hNchVsPtPosTPC_%s",Region[r]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

			hNchVsPtNegTPC[r] = new TH2F(Form("hNchVsPtNegTPC_%s",Region[r]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

			hNchVsPtPosTOF[r] = new TH2F(Form("hNchVsPtPosTOF_%s",Region[r]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

			hNchVsPtNegTOF[r] = new TH2F(Form("hNchVsPtNegTOF_%s",Region[r]),";#it{p}_{T} (GeV/#it{c}); #it{N}_{acc}^{TS}",nPtBins,ptBins,nBinsRT,binsRT);

		}

		for( int j = 0; j < nHists; j++ ){

			if(r<(nRegion-1)){
				hNchVsPtDataPosPionTPC[r][j] = new TH3F(Form("hNchVsPtData_Pos_Pion_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataNegPionTPC[r][j] = new TH3F(Form("hNchVsPtData_Neg_Pion_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataPosKaonTPC[r][j] = new TH3F(Form("hNchVsPtData_Pos_Kaon_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataNegKaonTPC[r][j] = new TH3F(Form("hNchVsPtData_Neg_Kaon_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataPosProtonTPC[r][j] = new TH3F(Form("hNchVsPtData_Pos_Proton_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);

				hNchVsPtDataNegProtonTPC[r][j] = new TH3F(Form("hNchVsPtData_Neg_Proton_TPC_%s_%s",Region[r],ending[j]),";#it{p}_{T}^{rec};n#sigma;#it{N}_{acc}",nPtBins,ptBins,nBinsNsigma,binsNsigma,nBinsRT,binsRT);
			}

			hNchVsPtDataPosTOF[r][j] = new TH3F(Form("hNchVsPtData_Pos_TOF_%s_%s",Region[r],ending[j]),";#it{p}^{rec};#beta;#it{N}_{acc}",nPtBins,ptBins,nBetaBins,BetaBins,nBinsRT,binsRT);

			hNchVsPtDataNegTOF[r][j] = new TH3F(Form("hNchVsPtData_Neg_TOF_%s_%s",Region[r],ending[j]),";#it{p}^{rec};#beta;#it{N}_{acc}",nPtBins,ptBins,nBetaBins,BetaBins,nBinsRT,binsRT);

			hDeDxVsP[r][j] = new TH3F(Form("hDeDxVsP_%s_%s",Region[r],ending[j]),";#it{p} (GeV/#it{c});#it{d}E/#it{d}x;#it{N}_{acc}^{TS}",nPtBins,ptBins,ndEdxBins,dEdxBins,nBinsRT,binsRT);

			if(!fAnalysisMC){

				if(r<(nRegion-1)){
					fListOfObjects->Add(hNchVsPtDataPosPionTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataNegPionTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataPosKaonTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataNegKaonTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataPosProtonTPC[r][j]);
					fListOfObjects->Add(hNchVsPtDataNegProtonTPC[r][j]);
				}

				fListOfObjects->Add(hNchVsPtDataPosTOF[r][j]);
				fListOfObjects->Add(hNchVsPtDataNegTOF[r][j]);
				fListOfObjects->Add(hDeDxVsP[r][j]);
			}

		}	// ending 

		if(!fAnalysisMC){
			if(r<(nRegion-1)){
				fListOfObjects->Add(hPhiData[r]);	
				fListOfObjects->Add(hNchVsPtPosTPC[r]);
				fListOfObjects->Add(hNchVsPtNegTPC[r]);
				fListOfObjects->Add(hNchVsPtPosTOF[r]);
				fListOfObjects->Add(hNchVsPtNegTOF[r]);
			}
		}

	}

	////	fEventCuts.AddQAplotsToList(fListOfObjects);
	PostData(1, fListOfObjects);

}
//______________________________________________________________________________
void AliAnalysisTaskSpectraRT::UserExec(Option_t *)
{
	// Main loop

	//
	// First we make sure that we have valid input(s)!
	//

	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}

	fESD = dynamic_cast<AliESDEvent*>(event);
	if(!fESD){
		Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
		this->Dump();
		return;
	}

	if (fAnalysisMC){

		fMC = dynamic_cast<AliMCEvent*>(MCEvent());
		if(!fMC){
			Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}

		fMCStack = fMC->Stack();

		if(!fMCStack){
			Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}

	}

	utils = new AliAnalysisUtils();
	if (!utils)
	{
		cout<<"------- No AnalysisUtils Object Found --------"<<utils<<endl;
		return;
	}

	fEvents->Fill(0.5);

	AliHeader* headerMC = fMC->Header();
	bool isGoodVtxPosMC = kFALSE;
	if (fAnalysisMC){
		AliGenEventHeader* genHeader = headerMC->GenEventHeader();
		TArrayF vtxMC(3); // primary vertex  MC
		vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
		if (genHeader) {
			genHeader->PrimaryVertex(vtxMC);
		}

		if(TMath::Abs(vtxMC[2])<=10){
			isGoodVtxPosMC = kTRUE;
		}

		// Before trigger selection
		GetLeadingObject(kTRUE);// leading particle at gen level
	}

	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	bool isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;
	fEvents->Fill(1.5);

	int INEL = -1;
	INEL = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0);
	if( INEL < 1 )
		return;
	fEvents->Fill(2.5);

	if( utils->IsSPDClusterVsTrackletBG(fESD) )
		return;
	fEvents->Fill(3.5);

	if( fESD->IsPileupFromSPDInMultBins() )
		return;
	fEvents->Fill(4.5);

	if( fESD->IsIncompleteDAQ())
		return;
	fEvents->Fill(5.5);

	if( !selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE) )
		return;
	fEvents->Fill(6.5);

	if( !IsGoodZvertexPos(fESD) )
		return;
	fEvents->Fill(7.5);

	GetLeadingObject(kFALSE);


	if(fIsMCclosure){
		double randomUE = gRandom->Uniform(0.0,1.0);
		if(randomUE<0.5){// corrections (50% stat.)
			if(isGoodVtxPosMC){
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
					GetDetectorResponse();
					GetMCCorrections();
				}
				//				if(fGenLeadPt>=fPtMin){
				//				}
			}
		}
		else{// for testing the method
			if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
				GetMultiplicityDistributions();
			}
		}
	}
	else{
		if(fAnalysisMC){
			if(isGoodVtxPosMC){
				if( (fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax)&&(fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax)){
					GetMultiplicityDistributions();
					GetDetectorResponse();
					GetMCCorrections();
				}

				//				if(fGenLeadPt>=fPtMin){
				//				}
			}
		}
		else{
			if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
				ProduceArrayTrksESD();
			}

			ProduceArrayV0ESD();

		}
	}

	PostData(1, fListOfObjects);
}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::GetLeadingObject(bool isMC) {

	Double_t flPt = 0;// leading pT
	Double_t flPhi = 0;
	Int_t flIndex = 0;

	if(isMC){
		for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

			AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
			if (!particle) continue;

			if (!fMC->IsPhysicalPrimary(i)) continue;  
			if (particle->Charge() == 0) continue;
			if ( TMath::Abs(particle->Eta()) > fEtaCut )continue;
			if( particle->Pt() < fPtMin)continue;

			if (flPt<particle->Pt()){
				flPt = particle->Pt();
				flPhi = particle->Phi();
				flIndex = i;
			}
		}

		fGenLeadPhi = flPhi;
		fGenLeadPt  = flPt;
		fGenLeadIn  = flIndex;
	}
	else{

		int iTracks(fESD->GetNumberOfTracks());          
		for(int i=0; i < iTracks; i++) {                

			AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 

			if(!track) continue;
			if(!fTrackFilterGolden->IsSelected(track)) continue;
			if(TMath::Abs(track->Eta()) > fEtaCut) continue;
			if(track->Pt() < fPtMin) continue;

			if (flPt<track->Pt()){
				flPt  = track->Pt();
				flPhi = track->Phi();
				flIndex = i;
			}

		}

		fRecLeadPhi = flPhi;
		fRecLeadPt  = flPt;
		fRecLeadIn  = flIndex;

	}

}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::GetMultiplicityDistributions(){

	int multTSgen=0;
	int multTSrec=0;
	const double pi = TMath::Pi();

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (i==fGenLeadIn) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)continue;
		if (particle->Pt() < fPtMin)continue;

		double DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;
		}
		else{
			multTSgen++;
		}
	}

	int iTracks = 0;           
	iTracks = fESD->GetNumberOfTracks();           

	for(int i = 0; i < iTracks; i++) {                 

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i));  
		if(!esdtrack) continue;
		if(esdtrack->Charge() == 0 ) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->Pt() < fPtMin) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->IsSelected(esdtrack)) { continue; } 
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kTRUE,kTRUE,kTRUE);
			if(!track) { continue; }
		}

		hPhiTotal->Fill(track->Eta(),track->Phi());
		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;
		}
		else{
			multTSrec++;
		}

	}
}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::GetDetectorResponse() {

	int multTSgen = 0;
	int multTSrec = 0;
	double pi = TMath::Pi();

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

		if(i==fGenLeadIn) continue;
		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)continue;
		if (particle->Pt() < fPtMin)continue;

		double DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;		
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;		
		}
		else{
			multTSgen++;
		}
	}

	int iTracks = 0;          
	iTracks = fESD->GetNumberOfTracks();          

	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;
		if(esdtrack->Charge() == 0 ) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->Pt() < fPtMin) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->IsSelected(esdtrack)) { continue; } 
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kTRUE,kTRUE,kTRUE);
			if(!track) { continue; }
		}

		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;		
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;		
		}
		else{
			multTSrec++;
		}
	}

}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::GetMCCorrections(){

	int multTSgen = 0;
	int multTSrec = 0;
	double pi = TMath::Pi();

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++){

		if (i==fGenLeadIn) continue;
		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)continue;
		if (particle->Pt() < fPtMin) continue;

		double DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;
		}
		else{
			multTSgen++;
		}
	}

	int iTracks(fESD->GetNumberOfTracks());          
	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->Pt() < fPtMin) continue;
		if(esdtrack->Charge()==0 ) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->IsSelected(esdtrack)) { continue; } 
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kTRUE,kTRUE,kTRUE);
			if(!track) { continue; }
		}

		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);
		if(TMath::Abs(DPhi)<pi/3.0){
			continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			continue;
		}
		else{
			multTSrec++;
		}
	}
}
//_____________________________________________________________________________
double AliAnalysisTaskSpectraRT::DeltaPhi(Double_t phi, Double_t Lphi,
		Double_t rangeMin, Double_t rangeMax)
{

	Double_t dphi = -999;
	Double_t pi = TMath::Pi();
//	if(Lphi > 2*pi || Lphi < 0)cout << "Lphi :: " << Lphi << endl;
//	if(phi  > 2*pi || phi < 0)cout << "phi = " << phi << endl;

	if(phi < 0)          phi += 2*pi;
	else if(phi > 2*pi)  phi -= 2*pi;
	if(Lphi < 0)         Lphi += 2*pi;
	else if(Lphi > 2*pi) Lphi -= 2*pi;
	dphi = Lphi - phi;
	if (dphi < rangeMin)      dphi += 2*pi;
	else if (dphi > rangeMax) dphi -= 2*pi;

	return dphi;
}
//_____________________________________________________________________________
short AliAnalysisTaskSpectraRT::GetPidCode(Int_t pdgCode) const
{
	// return our internal code for pions, kaons, and protons

	short pidCode = 6;

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
		case 11:
			pidCode = 4; // electron
			break;
		case 13:
			pidCode = 5; // muon
			break;
		default:
			pidCode = 6;  // something else?
	};

	return pidCode;
}
//_____________________________________________________________________________
void AliAnalysisTaskSpectraRT::ProduceArrayTrksESD(){

	const double pi = TMath::Pi();
	int multTSdata = 0;

	int iTracks(fESD->GetNumberOfTracks());          
	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!esdtrack) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->Pt() < fPtMin) continue;

		AliESDtrack* track = 0x0;
		if(!fSelectHybridTracks){
			if(!fTrackFilter->IsSelected(esdtrack)) { continue; } 
			else{ track = esdtrack; }
		}else{
			track = SetHybridTrackCuts(esdtrack,kTRUE,kTRUE,kTRUE);
			if(!track) { continue; }
		}

		hPhiTotal->Fill(track->Eta(),track->Phi());
		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			hPhiData[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hPhiData[1]->Fill(DPhi);
		}
		else{
			hPhiData[2]->Fill(DPhi);
			multTSdata++;
		}

	}

	hNchTSData->Fill(multTSdata);

	for(int iT = 0; iT < iTracks; iT++) {

		if(iT==fRecLeadIn) continue;
		AliESDtrack* esdTrack = (AliESDtrack*)fESD->GetTrack(iT);

		if(TMath::Abs(esdTrack->Eta()) > fEtaCut) continue;
		if(esdTrack->GetTPCsignalN() < fNcl) continue;
		if(esdTrack->Pt() < fPtMin) continue;
		if(!fTrackFilterGolden->IsSelected(esdTrack)) continue;

		double DPhi = DeltaPhi(esdTrack->Phi(), fRecLeadPhi);

		int nh = -1;
		double eta = esdTrack->Eta();
		if(TMath::Abs(eta)<0.2)
			nh = 0;
		else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
			nh = 1;
		else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
			nh = 2;
		else
			nh = 3;

		if( nh < 0 )
			continue;

		if(TMath::Abs(DPhi)<pi/3.0){
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTPC[0]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosPionTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataPosKaonTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataPosProtonTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTPC[0]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegPionTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataNegKaonTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataNegProtonTPC[0][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTPC[1]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosPionTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataPosKaonTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataPosProtonTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTPC[1]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegPionTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataNegKaonTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataNegProtonTPC[1][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
		}
		else{
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTPC[2]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosPionTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataPosKaonTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataPosProtonTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTPC[2]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegPionTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),multTSdata);	
				hNchVsPtDataNegKaonTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),multTSdata);	
				hNchVsPtDataNegProtonTPC[2][nh]->Fill(esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),multTSdata);	
			}
		}

		hPtVsP[nh]->Fill(esdTrack->P(),esdTrack->Pt());

		//
		//_______________________________ TOF PID
		//


		bool IsTOFout = kFALSE;
		IsTOFout = TOFPID(esdTrack);
		if(!IsTOFout) continue;

		double trkLength = esdTrack->GetIntegratedLength();
		double beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);

		if(TMath::Abs(DPhi)<pi/3.0){
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTOF[0]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosTOF[0][nh]->Fill(esdTrack->P(),beta,multTSdata);
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTOF[0]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegTOF[0][nh]->Fill(esdTrack->P(),beta,multTSdata);	
			}
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTOF[1]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosTOF[1][nh]->Fill(esdTrack->P(),beta,multTSdata);
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTOF[1]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegTOF[1][nh]->Fill(esdTrack->P(),beta,multTSdata);	
			}
		}
		else{
			if(esdTrack->Charge() > 0){
				hNchVsPtPosTOF[2]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataPosTOF[2][nh]->Fill(esdTrack->P(),beta,multTSdata);
			}
			if(esdTrack->Charge() < 0){
				hNchVsPtNegTOF[2]->Fill(esdTrack->Pt(),multTSdata);
				hNchVsPtDataNegTOF[2][nh]->Fill(esdTrack->P(),beta,multTSdata);	
			}
		}

		if(esdTrack->Charge() > 0)
			hNchVsPtDataPosTOF[3][nh]->Fill(esdTrack->P(),beta,multTSdata);

		if(esdTrack->Charge() < 0)
			hNchVsPtDataNegTOF[3][nh]->Fill(esdTrack->P(),beta,multTSdata);	


		//
		//_______________________________ rTPC PID
		//

		double momentum = esdTrack->P();
		float  dedx     = esdTrack->GetTPCsignal();

		if(!PhiCut(esdTrack->Pt(), esdTrack->Phi(), esdTrack->Charge(), Magf, fcutLow, fcutHigh))
			continue;

		if(fdEdxCalibrated){
			int index = -1;
			index = GetIndex();
			dedx *= 50/EtaCalibration(index,eta);
		}

		if( (momentum <= 0.6)&&(momentum >= 0.4) ){//only p:0.4-0.6 GeV, pion MIP
			if( (esdTrack->GetTPCsignal() < fDeDxMIPMax) && (esdTrack->GetTPCsignal() > fDeDxMIPMin) ){
				hMIPVsEta->Fill(eta,dedx);
				pMIPVsEta->Fill(eta,dedx);
			}
			if( (esdTrack->GetTPCsignal() > 70.0) && (esdTrack->GetTPCsignal() < 90.0) ){
				if( TMath::Abs(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kElectron))<2.0 ){
					hPlateauVsEta->Fill(eta,dedx);
					pPlateauVsEta->Fill(eta,dedx);
				}
			}
		}

		if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion))<2.0 )
			histPiTof[nh]->Fill(momentum,dedx);

		/*if( momentum <= 0.6 && momentum >= 0.4  ){
		  if( dedx < fDeDxMIPMax && dedx > fDeDxMIPMin ){
		  hMIPVsPhi[nh]->Fill(phi,dedx);
		  pMIPVsPhi[nh]->Fill(phi,dedx);
		  }
		  if( dedx > 70 && dedx < 90 ){
		  if( TMath::Abs(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kElectron))<2.0 ){
		  hPlateauVsPhi[nh]->Fill(phi,dedx);
		  pPlateauVsPhi[nh]->Fill(phi,dedx);
		  }
		  }
		  }*/

		if(TMath::Abs(DPhi)<pi/3.0){
			hDeDxVsP[0][nh]->Fill(momentum,dedx,multTSdata);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hDeDxVsP[1][nh]->Fill(momentum,dedx,multTSdata);
		}
		else{
			hDeDxVsP[2][nh]->Fill(momentum,dedx,multTSdata);
		}

		hDeDxVsP[3][nh]->Fill(momentum,dedx,multTSdata);


	}//end of track loop
}
//________________________________________________________________________
void AliAnalysisTaskSpectraRT::ProduceArrayV0ESD(){


	int index = -1;
	index = GetIndex();

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

	for (Int_t iV0=0; iV0<nv0s; iV0++) {

		AliESDv0 *esdV0 = fESD->GetV0(iV0);
		if ( !esdV0 ) continue;

		//check onfly status
		//              if( !esdV0->GetOnFlyStatus() )
		//                      continue;

		if( esdV0->GetOnFlyStatus()!=0 )
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
		//              if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
		//                      continue;

		UInt_t selectDebug_p = 0;
		if ( fTrackFilterGolden ) {
			selectDebug_p = fTrackFilterGolden->IsSelected(pTrack);
			if (!selectDebug_p) {
				continue;
			}
		}

		UInt_t selectDebug_n = 0;
		if ( fTrackFilterGolden ) {
			selectDebug_n = fTrackFilterGolden->IsSelected(nTrack);
			if (!selectDebug_n) {
				continue;
			}
		}


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

						       if(track->GetTPCsignalN()< fNcl )
							       continue;

						       if(!PhiCut(track->Pt(),track->Phi(), track->Charge(), Magf, fcutLow, fcutHigh))
							       continue;

						       Double_t eta      = track->Eta();
						       Double_t momentum = track->P();
						       Double_t dedx     = track->GetTPCsignal();

						       if(fdEdxCalibrated)
							       dedx *= 50/EtaCalibration(index,eta);						      

						       if(fillPos&&fillNeg){
							       if( (track->GetTPCsignal() < fDeDxMIPMax) && (track->GetTPCsignal() > fDeDxMIPMin) ){
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s->Fill(eta,dedx);
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

						       if(nh<0)
							       continue;

						       if(fillPos&&fillNeg){
							       histPiV0[nh]->Fill(track->P(), dedx);
						       }
						       else{
							       histPV0[nh]->Fill(track->P(), dedx);
						       }

					       }//end loop over two tracks

				       };
				       break;

				case 1:{//gammas

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if( dmassK>0.01 && dmassL>0.01 && dmassAL>0.01 ) {
						       if( dmassG<0.01 && dmassG>0.0001 ) {
							       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationEl(index,nTrack->Eta())) < 5)
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

					       if(fdEdxCalibrated)
						       dedx *= 50/EtaCalibration(index,track->Eta());						      

					       if(track->GetTPCsignalN() <= fNcl)
						       continue;

					       if(!PhiCut(track->Pt(), track->Phi(), track->Charge(), Magf, fcutLow, fcutHigh))
						       continue;

					       Int_t nh = -1;

					       if(TMath::Abs(eta)<0.2)
						       nh = 0;
					       else if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)<0.4)
						       nh = 1;
					       else if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)<0.6)
						       nh = 2;
					       else if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)<0.8)
						       nh = 3;

					       if(nh<0)
						       continue;

					       histEV0[nh]->Fill(track->P(),track->GetTPCsignal());

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
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::selectVertex2015pp(AliESDEvent *esd,
		Bool_t checkSPDres, //enable check on vtx resolution
		Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
		Bool_t checkProximity) //apply cut on relative position of spd and trk verteces
{

	if (!esd) return kFALSE;

	const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
	const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
	bool hasSPD = spdVertex->GetStatus();
	bool hasTrk = trkVertex->GetStatus();

	//Note that AliVertex::GetStatus checks that N_contributors is > 0
	//reject events if both are explicitly requested and none is available
	if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;

	//reject events if none between the SPD or track verteces are available
	//if no trk vertex, try to fall back to SPD vertex;
	if (!hasTrk) {
		if (!hasSPD) return kFALSE;
		//on demand check the spd vertex resolution and reject if not satisfied
		if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
	} else {
		if (hasSPD) {
			//if enabled check the spd vertex resolution and reject if not satisfied
			//if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
			if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
			if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
		}
	}

	//Cut on the vertex z position
	//const AliESDVertex * vertex = esd->GetPrimaryVertex();
	//if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::IsGoodSPDvertexRes(const AliESDVertex* spdVertex)
{

	if( !spdVertex ) return kFALSE;
	if( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::IsGoodZvertexPos(AliESDEvent *esd)
{

	if( !esd ) return kFALSE;
	//Cut on the vertex z position
	const AliESDVertex * vertex = esd->GetPrimaryVertex();
	if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::PhiCut(const double& pt, double phi, const double& q, const float& mag, TF1* phiCutLow, TF1* phiCutHigh)
{
	if(pt < 2.0)
		return kTRUE;

	if(mag < 0)    // for negatve polarity field
		phi = TMath::TwoPi() - phi;
	if(q < 0) // for negatve charge
		phi = TMath::TwoPi()-phi;

	phi += TMath::Pi()/18.0; // to center gap in the middle
	phi = fmod(phi, TMath::Pi()/9.0);

	if(phi<phiCutHigh->Eval(pt)
			&& phi>phiCutLow->Eval(pt))
		return kFALSE; // reject track

	//    hPhi[4]->Fill(pt, phi);

	return kTRUE;
}
//________________________________________________________________________
float AliAnalysisTaskSpectraRT::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){

	double maxDCAxy = 10;
	maxDCAxy = fMaxDCAxy->Eval(ptI);
	return maxDCAxy;

}
//________________________________________________________________________
void AliAnalysisTaskSpectraRT::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

	AliESDtrackCuts* esdTrackCuts = 0x0;
	if(fSetTPConlyTrkCuts){
		esdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
		esdTrackCuts->SetRequireTPCRefit(kTRUE);
		esdTrackCuts->SetRequireITSRefit(kTRUE);
		esdTrackCuts->SetEtaRange(-0.8,0.8);
	}
	else{
		esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
		esdTrackCuts->SetEtaRange(-0.8,0.8);
	}

	fTrackFilter->AddCuts(esdTrackCuts);
}
//________________________________________________________________________
AliESDtrack* AliAnalysisTaskSpectraRT::SetHybridTrackCuts(AliESDtrack *esdtrack, const bool fillPhiStand, const bool fillPhHyb1, const bool fillPhHyb2){

	// 
	// 	Get the Hybrid Tracks 
	// 	

	AliESDtrack *newTrack = 0x0;

	//	if(fTrackCutsType==0 || fTrackCutsType==3)
	//	{
	if(fTrackFilter->IsSelected(esdtrack))
	{
		newTrack = new AliESDtrack(*esdtrack);
		if(fillPhiStand) hPhiStandard->Fill(newTrack->Eta(),newTrack->Phi());
		////			newTrack->SetTRDQuality(0);
	}
	else if(fHybridTrackCuts1->AcceptTrack(esdtrack))
	{
		if(esdtrack->GetConstrainedParam())
		{
			newTrack = new AliESDtrack(*esdtrack);
			const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
			newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
			////				newTrack->SetTRDQuality(1);
			if(fillPhHyb1) hPhiHybrid1->Fill(newTrack->Eta(),newTrack->Phi());
		}
		else
			return 0x0;
	}
	/*else if(fHybridTrackCuts2->AcceptTrack(esdtrack))
	  {
	  if(esdtrack->GetConstrainedParam())
	  {
	  newTrack = new AliESDtrack(*esdtrack);
	  const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
	  newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
	/////				newTrack->SetTRDQuality(2);
	}
	else
	return 0x0;
	}*/
	else
	{
		return 0x0;
	}
	//	}

	return newTrack;

}
//________________________________________________________________________
double AliAnalysisTaskSpectraRT::EtaCalibration( const int &indx, const double &eta){
	//    h        i         j        l      k        o          p
	const Double_t aPos[nRt+2]      = {49.9044 ,50.0841  ,49.8419 ,49.9799 ,49.9659 ,50.0535 ,50.0649};
	const Double_t bPos[nRt+2]      = {4.05075 ,-0.743724,10.3952 ,2.99619 ,2.91366 ,-2.87404,-4.3589};
	const Double_t cPos[nRt+2]      = {-58.1027,-13.8508 ,-151.227,-45.718 ,-45.5994,31.159  ,66.268 };
	const Double_t dPos[nRt+2]      = {342.297 ,141.269  ,900.83  ,290.013 ,290.042 ,-151.257,-435.715};
	const Double_t ePos[nRt+2]      = {-1098.19,-567.054 ,-2833.25,-1018.42,-1014.49,282.703 ,1325.24};
	const Double_t fPos[nRt+2]      = {1944.32 ,1105.63  ,4884.59 ,1948.68 ,1931.84 ,-98.8756,-2030.25};
	const Double_t gPos[nRt+2]      = {-1749.14,-1022.19 ,-4318.99,-1864.06,-1839.36,-230.114,1542.09};
	const Double_t hPos[nRt+2]      = {617.929 ,355.158  ,1521.66 ,692.752 ,680.421 ,172.854 ,-468.577};

	const Double_t aNeg[nRt+2]      = {49.9261,50.0561,49.9583,50.078 ,50.046 ,49.9496,50.1258};
	const Double_t bNeg[nRt+2]      = {2.92422,4.68965,3.38038,6.67199,6.79992,2.45301,12.7977};
	const Double_t cNeg[nRt+2]      = {61.6661,65.891 ,53.1256,103.662,109.86 ,53.654 ,190.076};
	const Double_t dNeg[nRt+2]      = {421.545,394.542,314.489,611.034,668.241,363.689,1144.11};
	const Double_t eNeg[nRt+2]      = {1283.04,1100.56,825.296,1695.63,1916.44,1115.13,3411.98};
	const Double_t fNeg[nRt+2]      = {1944.85,1516.21,1021.01,2395.88,2815.04,1762.18,5402.99};
	const Double_t gNeg[nRt+2]      = {1442.98,989.24 ,548.44 ,1669.22,2057.21,1421.46,4379.16};
	const Double_t hNeg[nRt+2]      = {419.491,238.333,84.7945,455.362,595.391,469.45 ,1436.76};


	for(Int_t i=0; i<8; ++i)
		fEtaCalibration->SetParameter(i,0);

	if(eta<0){
		fEtaCalibration->SetParameter(0,aNeg[indx]);
		fEtaCalibration->SetParameter(1,bNeg[indx]);
		fEtaCalibration->SetParameter(2,cNeg[indx]);
		fEtaCalibration->SetParameter(3,dNeg[indx]);
		fEtaCalibration->SetParameter(4,eNeg[indx]);
		fEtaCalibration->SetParameter(5,fNeg[indx]);
		fEtaCalibration->SetParameter(6,gNeg[indx]);
		fEtaCalibration->SetParameter(7,hNeg[indx]);
	}
	else{
		fEtaCalibration->SetParameter(0,aPos[indx]);
		fEtaCalibration->SetParameter(1,bPos[indx]);
		fEtaCalibration->SetParameter(2,cPos[indx]);
		fEtaCalibration->SetParameter(3,dPos[indx]);
		fEtaCalibration->SetParameter(4,ePos[indx]);
		fEtaCalibration->SetParameter(5,fPos[indx]);
		fEtaCalibration->SetParameter(6,gPos[indx]);
		fEtaCalibration->SetParameter(7,hPos[indx]);
	}

	return fEtaCalibration->Eval(eta);

}
//________________________________________________________________________
double AliAnalysisTaskSpectraRT::EtaCalibrationEl(const int &indx, const double &eta){

	const Double_t aPosEl[nRt+2]    = {79.8647 ,79.6737 ,80.3915 ,80.1263 ,79.9957 ,79.6537 ,80.6434 };
	const Double_t bPosEl[nRt+2]    = {6.50512 ,16.0745 ,9.53925 ,5.28525 ,7.03079 ,15.0221 ,0.40293 };
	const Double_t cPosEl[nRt+2]    = {-35.9277,-80.5639,-69.3773,-32.7731,-42.9098,-83.6391,-21.8162};
	const Double_t dPosEl[nRt+2]    = {73.1535 ,148.866 ,143.956 ,68.4524 ,88.7057 ,168.5   ,61.9147 };
	const Double_t ePosEl[nRt+2]    = {-47.1041,-90.3376,-89.5518,-44.1566,-56.6554,-107.999,-44.6593};

	const Double_t aNegEl[nRt+2]    = {79.6366 ,80.0767 ,79.6157 ,79.8351 ,79.7387 ,79.3638 ,79.9111 };
	const Double_t bNegEl[nRt+2]    = {-11.3437,-2.51009,-16.2468,-8.46921,-8.60021,-17.1977,-1.66066};
	const Double_t cNegEl[nRt+2]    = {-65.1353,-23.6188,-92.0783,-44.5947,-44.1718,-82.7998,-6.96109};
	const Double_t dNegEl[nRt+2]    = {-134.447,-65.5053,-180.753,-86.2242,-84.4984,-143.394,-16.0465};
	const Double_t eNegEl[nRt+2]    = {-87.7848,-51.1463,-112.997,-53.6285,-51.945 ,-81.3439,-10.3587};


	for(Int_t i=0; i<5; ++i)
		fEtaCalibrationEl->SetParameter(i,0);

	if(eta<0){
		fEtaCalibrationEl->SetParameter(0,aNegEl[indx]);
		fEtaCalibrationEl->SetParameter(1,bNegEl[indx]);
		fEtaCalibrationEl->SetParameter(2,cNegEl[indx]);
		fEtaCalibrationEl->SetParameter(3,dNegEl[indx]);
		fEtaCalibrationEl->SetParameter(4,eNegEl[indx]);
	}
	else{
		fEtaCalibrationEl->SetParameter(0,aPosEl[indx]);
		fEtaCalibrationEl->SetParameter(1,bPosEl[indx]);
		fEtaCalibrationEl->SetParameter(2,cPosEl[indx]);
		fEtaCalibrationEl->SetParameter(3,dPosEl[indx]);
		fEtaCalibrationEl->SetParameter(4,ePosEl[indx]);
	}

	return fEtaCalibrationEl->Eval(eta);

}
//________________________________________________________________________
int AliAnalysisTaskSpectraRT::GetIndex()
{

	Int_t indx = -1;

	if(fPeriod=="h")
		indx = 0;
	else if(fPeriod=="i")
		indx = 1;
	else if(fPeriod=="j")
		indx = 2;
	else if(fPeriod=="l")
		indx = 3;
	else if(fPeriod=="k")
		indx = 4;
	else if(fPeriod=="o")
		indx = 5;
	else
		indx = 6;

	return indx;

}
//________________________________________________________________________
bool AliAnalysisTaskSpectraRT::TOFPID(AliESDtrack * track)
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

