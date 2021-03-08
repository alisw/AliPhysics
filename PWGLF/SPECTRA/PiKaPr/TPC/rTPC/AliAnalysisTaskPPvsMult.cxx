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

#include "AliAnalysisTaskPPvsMult.h"

// ROOT includes
#include <TList.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
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
#include <AliStack.h>

#include <TTreeStream.h>

#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenDPMjetEventHeader.h>

#include <AliCentrality.h>
#include <AliESDv0.h>
#include <AliKFVertex.h>
#include <AliAODVertex.h>

#include <AliAODTrack.h>
#include <AliVParticle.h>
#include <AliPID.h>
#include <AliAODPid.h>
#include <AliAODMCHeader.h>


// STL includes
#include <iostream>
using namespace std;


//
// Responsible:
// Antonio Ortiz (Lund)
// Peter Christiansen (Lund)
//


const Double_t AliAnalysisTaskPPvsMult::fgkClight = 2.99792458e-2;
Float_t Magf                    = 1;
const Int_t nHists              = 4;
const Int_t nCent               = 11;
const Double_t CentMin[nCent]   = {0.0,1.0,5.0 ,10.0,15.0,20.0,30.0,40.0,50.0,70.0,0.0};
const Double_t CentMax[nCent]   = {1.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,70.0,100.0,100.0};
static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps

ClassImp(AliAnalysisTaskPPvsMult)
	AliAnalysisTaskPPvsMult::AliAnalysisTaskPPvsMult():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fAOD(0x0),
		fEventCuts(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fMCArray(0x0),
		fPIDResponse(0x0),
		fTrackFilterGolden(0x0),
		fTrackFilterTPC(0x0),
		utils(0x0),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		///fAnalysisPbPb(kFALSE),
		fRandom(0x0),
		fVtxCut(0x0),
		fNcl(70),
		fEtaCut(0.9),
		cent(3),
		//fMinCent(0.0),
		//fMaxCent(100.0),
		fDeDxMIPMin(40),
		fDeDxMIPMax(60),
		fdEdxHigh(200),
		fdEdxLow(40),
		fPeriod("16l"),
		fMcProcessType(-999),
		fTriggeredEventMB(-999),
		fVtxStatus(-999),
		fZvtx(-999),
		fZvtxMC(-999),
		fRun(-999),
		fEventId(-999),
		fListOfObjects(0),
		fEvents(0x0),
		fVtxMC(0x0),
		fdEdxCalibrated(0x0),
		fTrackCuts(0),
		////		fMakePid(0x0),
		////		fLowPt(0x0),
		fLHC16l(0x0),
		fcent(0x0),
		fcentAfterPrimaries(0x0),
		fcentAfterV0s(0x0),
		hDCAxyVsPtPiNeg_TPC(0x0),
		hDCAxyVsPtPNeg_TPC(0x0),
		hDCAxyVsPtPiNeg_TOF(0x0),
		hDCAxyVsPtPNeg_TOF(0x0),
		hDCAxyVsPtPiPos_TPC(0x0),
		hDCAxyVsPtPPos_TPC(0x0),
		hDCAxyVsPtPiPos_TOF(0x0),
		hDCAxyVsPtPPos_TOF(0x0),
		fEtaCalibrationPos(0x0),
		fEtaCalibrationNeg(0x0),
		felededxfitPos(0x0),
		felededxfitNeg(0x0),
		fcutDCAxy(0x0),
		fcutLow(0x0),
		fcutHigh(0x0)

{

	for(Int_t j=0; j<nHists; ++j){

		hMIPVsV0M[j]=0;//TH2D, MIP vs V0M Mult. for different eta intervals
		pMIPVsV0M[j]=0;//TProfile, MIP vs V0M Mult. for different eta intervals
		hMIPVsNch[j]=0;//TH2D, MIP vs Nch for different eta intervals
		pMIPVsNch[j]=0;//TProfile, MIP vs Nch for different eta intervals

	}

	for(Int_t i = 0; i<nCent; ++i){


		// Histograms for PreCalibration

		hMIPVsEta[i]=0;
		pMIPVsEta[i]=0;
		hMIPVsEtaV0s[i]=0;
		pMIPVsEtaV0s[i]=0;
		hPlateauVsEta[i]=0;
		pPlateauVsEta[i]=0;
		hPhi[i]=0;

		// Histograms for PostCalibration

		hPt_TPC[i]=0;
		hPtpos_TPC[i]=0;
		hPtneg_TPC[i]=0;

		hPtpos_TOF[i]=0;
		hPtneg_TOF[i]=0;

		for(Int_t j=0; j<nHists; ++j){

			hPtpos_TPC_Eta[i][j]=0;//TH1D, Transverse momentum distribution  positive-charged hadrons
			hPtneg_TPC_Eta[i][j]=0;//TH1D, Transverse momentum distribution  negative-charged hadrons
			hPtpos_TOF_Eta[i][j]=0;
			hPtneg_TOF_Eta[i][j]=0;
			hPpos_TOF_Eta[i][j]=0;
			hPneg_TOF_Eta[i][j]=0;
			hPtVsP[i][j]=0;//TH2D, Transverse momentum Vs momentum

			hMIPVsPhi[i][j]=0;//TH2D, MIP vs phi for different eta intervals
			pMIPVsPhi[i][j]=0;//TProfile, MIP vs phi for different eta intervals
			hPlateauVsPhi[i][j]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
			pPlateauVsPhi[i][j]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c

			hDeDxVsP[i][j]=0;//TH2D, DeDx vs P

			hnSigmaPiPos[i][j]=0;//TH2D, nSigmas vs Pt Pions
			hnSigmaKPos[i][j]=0;//TH2D, nSigmas vs Pt Kaons
			hnSigmaPPos[i][j]=0;//TH2D, nSigmas vs Pt Protons

			hnSigmaPiNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiPions
			hnSigmaKNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiKaons
			hnSigmaPNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiProtons

			hBetavsPpos[i][j]=0;
			hBetavsPneg[i][j]=0;

			histPiV0[i][j]=0;//TH2D, dE/dx vs p, pi id by V0s
			histpPiV0[i][j]=0;//TH1D, pi id by V0s
			histPV0[i][j]=0;// TH2D, dE/dx vs p, p id by V0s
			histpPV0[i][j]=0;// TH1D, p id by V0s
			histPiTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
			histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
			histEV0[i][j]=0;

		}
	}

	for(Int_t pid=0;pid<3;++pid){
		hMcInNeg[pid]  = 0;
		hMcInPos[pid]  = 0;
		hMcOutNeg[pid] = 0;
		hMcOutPos[pid] = 0;
		hMcOutNegTOF[pid] = 0;
		hMcOutPosTOF[pid] = 0;
	}

	for(Int_t pid=0;pid<2;++pid){
		for(Int_t q=0;q<2;++q){
			hDCApTPrim[pid][q]  = 0;
			hDCApTWDec[pid][q]  = 0;
			hDCApTMate[pid][q]  = 0;

			hDCApTPrim_TOF[pid][q]  = 0;
			hDCApTWDec_TOF[pid][q]  = 0;
			hDCApTMate_TOF[pid][q]  = 0;
		}
	}

}


AliAnalysisTaskPPvsMult::AliAnalysisTaskPPvsMult(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fAOD(0x0),
	fEventCuts(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fMCArray(0x0),
	fPIDResponse(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilterTPC(0x0),
	utils(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	////fAnalysisPbPb(kFALSE),
	fRandom(0x0),
	fVtxCut(0x0),
	fNcl(70),
	fEtaCut(0.9),
	cent(3),
	//fMinCent(0.0),
	//fMaxCent(100.0),
	fDeDxMIPMin(40),
	fDeDxMIPMax(60),
	fdEdxHigh(200),
	fdEdxLow(40),
	fPeriod("16l"),
	fMcProcessType(-999),
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtxMC(-999),
	fRun(-999),
	fEventId(-999),
	fListOfObjects(0), 
	fEvents(0x0),
	fVtxMC(0x0),
	fdEdxCalibrated(0x0),
	fTrackCuts(0),
	////fMakePid(0x0),
	////fLowPt(0x0),
	fLHC16l(0x0),
	fcent(0x0),
	fcentAfterPrimaries(0x0),
	fcentAfterV0s(0x0),
	hDCAxyVsPtPiNeg_TPC(0x0),
	hDCAxyVsPtPNeg_TPC(0x0),
	hDCAxyVsPtPiNeg_TOF(0x0),
	hDCAxyVsPtPNeg_TOF(0x0),
	hDCAxyVsPtPiPos_TPC(0x0),
	hDCAxyVsPtPPos_TPC(0x0),
	hDCAxyVsPtPiPos_TOF(0x0),
	hDCAxyVsPtPPos_TOF(0x0),
	fEtaCalibrationPos(0x0),
	fEtaCalibrationNeg(0x0),
	felededxfitPos(0x0),
	felededxfitNeg(0x0),
	fcutDCAxy(0x0),
	fcutLow(0x0),
	fcutHigh(0x0)

{

	for(Int_t j=0; j<nHists; ++j){

		hMIPVsV0M[j]=0;//TH2D, MIP vs V0M Mult. for different eta intervals
		pMIPVsV0M[j]=0;//TProfile, MIP vs V0M Mult. for different eta intervals

		hMIPVsNch[j]=0;//TH2D, MIP vs Nch for different eta intervals
		pMIPVsNch[j]=0;//TProfile, MIP vs Nch for different eta intervals
	}


	for(Int_t i = 0; i<nCent; ++i){

		hMIPVsEta[i] = 0;
		hMIPVsEta[i]=0;
		pMIPVsEta[i]=0;
		hMIPVsEtaV0s[i]=0;
		pMIPVsEtaV0s[i]=0;
		hPlateauVsEta[i]=0;
		pPlateauVsEta[i]=0;
		hPhi[i]=0;

		hPt_TPC[i]=0;
		hPtpos_TPC[i]=0;
		hPtneg_TPC[i]=0;

		hPtpos_TOF[i]=0;
		hPtneg_TOF[i]=0;

		for(Int_t j=0; j<nHists; ++j){

			hPtpos_TPC_Eta[i][j]=0;//TH1D, Transverse momentum distribution  positive-charged hadrons
			hPtneg_TPC_Eta[i][j]=0;//TH1D, Transverse momentum distribution  negative-charged hadrons
			hPtpos_TOF_Eta[i][j]=0;
			hPtneg_TOF_Eta[i][j]=0;
			hPpos_TOF_Eta[i][j]=0;
			hPneg_TOF_Eta[i][j]=0;
			hPtVsP[i][j]=0;//TH2D, Transverse momentum Vs momentum

			hMIPVsPhi[i][j]=0;//TH2D, MIP vs phi for different eta intervals
			pMIPVsPhi[i][j]=0;//TProfile, MIP vs phi for different eta intervals
			hPlateauVsPhi[i][j]=0;//TH2D, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c
			pPlateauVsPhi[i][j]=0;//TProfile, dE/dx vs Phi, electrons 0.4<p<0.6 GeV/c

			hDeDxVsP[i][j]=0;//TH2D, DeDx vs P

			hnSigmaPiPos[i][j]=0;//TH2D, nSigmas vs Pt Pions
			hnSigmaKPos[i][j]=0;//TH2D, nSigmas vs Pt Kaons
			hnSigmaPPos[i][j]=0;//TH2D, nSigmas vs Pt Protons

			hnSigmaPiNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiPions
			hnSigmaKNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiKaons
			hnSigmaPNeg[i][j]=0;//TH2D, nSigmas vs Pt AntiProtons

			hBetavsPpos[i][j]=0;
			hBetavsPneg[i][j]=0;

			histPiV0[i][j]=0;//TH2D, dE/dx vs p, pi id by V0s
			histpPiV0[i][j]=0;//TH1D, pi id by V0s
			histPV0[i][j]=0;// TH2D, dE/dx vs p, p id by V0s
			histpPV0[i][j]=0;// TH1D, p id by V0s
			histPiTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
			histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
			histEV0[i][j]=0;

		}

	}

	for(Int_t pid=0; pid<3; ++pid){
		hMcInNeg[pid]=0;
		hMcInPos[pid]=0;
		hMcOutNeg[pid]=0;
		hMcOutPos[pid]=0;
		hMcOutNegTOF[pid] = 0;
		hMcOutPosTOF[pid] = 0;
	}

	for(Int_t pid=0;pid<2;++pid){
		for(Int_t q=0;q<2;++q){
			hDCApTPrim[pid][q]=0;
			hDCApTWDec[pid][q]=0;
			hDCApTMate[pid][q]=0;

			hDCApTPrim_TOF[pid][q]  = 0;
			hDCApTWDec_TOF[pid][q]  = 0;
			hDCApTMate_TOF[pid][q]  = 0;
		}
	}

	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());//esto es nuevo
}




AliAnalysisTaskPPvsMult::~AliAnalysisTaskPPvsMult() {

	//
	// Destructor
	//

	if(fListOfObjects) {
		delete fListOfObjects;
		fListOfObjects = 0x0;
	}

}



//______________________________________________________________________________



void AliAnalysisTaskPPvsMult::UserCreateOutputObjects()
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

	// Track selection for the spectra: ITSTPC2011
	// The DCAxy cut is applied inside the track loop
	if(!fTrackFilterGolden){

		fTrackFilterGolden = new AliESDtrackCuts("fTrackFilterGolden");

		// TPC
		//fTrackFilterGolden->SetMinNCrossedRowsTPC(70); // Varied for track selection
		fTrackFilterGolden->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		//fTrackFilterGolden->SetMaxChi2PerClusterTPC(4); // Varied for track selection
		fTrackFilterGolden->SetAcceptKinkDaughters(kFALSE);
		fTrackFilterGolden->SetRequireTPCRefit(kTRUE);
		// ITS
		fTrackFilterGolden->SetRequireITSRefit(kTRUE);
		fTrackFilterGolden->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
		//fTrackFilterGolden->SetMaxChi2TPCConstrainedGlobal(36);
		//fTrackFilterGolden->SetMaxDCAToVertexZ(2); // Varied for track selection
		fTrackFilterGolden->SetDCAToVertex2D(kFALSE);
		fTrackFilterGolden->SetRequireSigmaToVertex(kFALSE);
		fTrackFilterGolden->SetMaxChi2PerClusterITS(36);

		printf("fTrackCuts - %d\n",fTrackCuts);

		if(fTrackCuts==0){
			fTrackFilterGolden->SetMinNCrossedRowsTPC(70);
			fTrackFilterGolden->SetMaxChi2PerClusterTPC(4);
			fTrackFilterGolden->SetMaxDCAToVertexZ(2);
		}else if(fTrackCuts==1){
			fTrackFilterGolden->SetMinNCrossedRowsTPC(80);
			fTrackFilterGolden->SetMaxChi2PerClusterTPC(4);
			fTrackFilterGolden->SetMaxDCAToVertexZ(2);
		}else if(fTrackCuts==2){
			fTrackFilterGolden->SetMinNCrossedRowsTPC(70);
			fTrackFilterGolden->SetMaxChi2PerClusterTPC(5);
			fTrackFilterGolden->SetMaxDCAToVertexZ(2);
		}else if(fTrackCuts==3){
			fTrackFilterGolden->SetMinNCrossedRowsTPC(70);
			fTrackFilterGolden->SetMaxChi2PerClusterTPC(4);
			fTrackFilterGolden->SetMaxDCAToVertexZ(3);
		}else if(fTrackCuts==4){
			fTrackFilterGolden->SetMinNCrossedRowsTPC(60);
			fTrackFilterGolden->SetMaxChi2PerClusterTPC(4);
			fTrackFilterGolden->SetMaxDCAToVertexZ(2);
		}else{
			fTrackFilterGolden->SetMinNCrossedRowsTPC(70);
			fTrackFilterGolden->SetMaxChi2PerClusterTPC(4);
			fTrackFilterGolden->SetMaxDCAToVertexZ(1);
		}

	}

	OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	//
	// Histograms
	//

	fEvents = new TH2F( "fEvents", ";Evt. Sel.;Mult. Per.",12,0,12,13,0,13);
	fEvents->GetXaxis()->SetBinLabel(1, "Processed");
	fEvents->GetXaxis()->SetBinLabel(2, "PhysSel+Trigger");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(3, "INEL>0");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(4, "BG");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(5, "IsPileUpFromSPDinMultBins");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(6, "Incom DAQ");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(7, "Res&Proximity");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(8, "|Vtz|<10cm");//NotinVertexcut");
	fEvents->GetYaxis()->SetBinLabel(2,Form("%.2f-%.2f",CentMin[0],CentMax[0]));
	fEvents->GetYaxis()->SetBinLabel(3,Form("%.2f-%.2f",CentMin[1],CentMax[1]));
	fEvents->GetYaxis()->SetBinLabel(4,Form("%.2f-%.2f",CentMin[2],CentMax[2]));
	fEvents->GetYaxis()->SetBinLabel(5,Form("%.2f-%.2f",CentMin[3],CentMax[3]));
	fEvents->GetYaxis()->SetBinLabel(6,Form("%.2f-%.2f",CentMin[4],CentMax[4]));
	fEvents->GetYaxis()->SetBinLabel(7,Form("%.2f-%.2f",CentMin[5],CentMax[5]));
	fEvents->GetYaxis()->SetBinLabel(8,Form("%.2f-%.2f",CentMin[6],CentMax[6]));
	fEvents->GetYaxis()->SetBinLabel(9,Form("%.2f-%.2f",CentMin[7],CentMax[7]));
	fEvents->GetYaxis()->SetBinLabel(10,Form("%.2f-%.2f",CentMin[8],CentMax[8]));
	fEvents->GetYaxis()->SetBinLabel(11,Form("%.2f-%.2f",CentMin[9],CentMax[9]));
	fEvents->GetYaxis()->SetBinLabel(12,"0.0-100.0");
	fListOfObjects->Add(fEvents);

	fV0M = new TH1F("fV0M", ";V0M Percentile",100,0,100);
	fListOfObjects->Add(fV0M);

	fcent=new TH1F("fcent","fcent",13,0,13);
	fcentAfterPrimaries =new TH1F("fcentAfterPrimaries","fcentAfterPrimaries",13,0,13);
	fcentAfterV0s =new TH1F("fcentAfterV0s","fcentAfterV0s",13,0,13);
	fListOfObjects->Add(fcent);
	fListOfObjects->Add(fcentAfterPrimaries);
	fListOfObjects->Add(fcentAfterV0s);

	const Int_t nDeltaPiBins   = 80;
	const Double_t deltaPiLow  = 20;
	const Double_t deltaPiHigh = 100;

	///const Char_t *Pid[7]       = {"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};
	const Char_t *Q[2]         = {"Neg", "Pos"};

	const Int_t nPtBins = 63;
	Double_t ptBins[nPtBins+1] = {

		0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35,
		0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
		0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
		1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40,
		3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00,
		9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 18.00, 
		20.00,22.00,24.00,26.00,30.00

	};

	/*
	   const Int_t nPtBins = 53;
	   Double_t ptBins[nPtBins+1] = {
	   0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
	   0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,
	   1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
	   4.5,5,5.5,6,6.5,7,8,10,13,20,30};
	   */

	const int nV0MMultBins   = 100;
	double V0MMultBins[nV0MMultBins+1] = { 0 };
	for(int i = 0; i <= nV0MMultBins; ++i){
		V0MMultBins[i] = i;
	}

	const int nDCABins   = 800;
	double DCABins[nDCABins+1] = { 0 };
	for(int i = 0; i <= nDCABins; ++i){
		DCABins[i] = -4.0 + ((double)i)/100.0;
	}

	const Int_t nPtBinsV0s = 25;
	Double_t ptBinsV0s[nPtBinsV0s+1] = { 0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
		1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
		9.0 , 12.0, 15.0, 20.0 };

	const Char_t* ending[nHists] = {"02", "24", "46", "68"};
	const Char_t* LatexEta[nHists] = {"|#eta|<0.2", "0.2<|#eta|<0.4", "0.4<|#eta|<0.6", "0.6<|#eta|<0.8" };

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+TMath::Pi()/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+TMath::Pi()/18.0+0.035", 0, 50);


	fEtaCalibrationPos = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);

	felededxfitPos     = new TF1("felededxfitPos", "pol4", 0.0, 1.0);
	felededxfitNeg     = new TF1("felededxfitNeg", "pol4", -1.0, 0.0);

	hDCAxyVsPtPiNeg_TPC = new TH3F(Form("hDCAxyVsPtPionNeg_TPC"), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
	fListOfObjects->Add(hDCAxyVsPtPiNeg_TPC);
	hDCAxyVsPtPNeg_TPC = new TH3F(Form("hDCAxyVsPtProtonNeg_TPC"), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
	fListOfObjects->Add(hDCAxyVsPtPNeg_TPC);
	hDCAxyVsPtPiNeg_TOF = new TH3F(Form("hDCAxyVsPtPionNeg_TOF"), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
	fListOfObjects->Add(hDCAxyVsPtPiNeg_TOF);
	hDCAxyVsPtPNeg_TOF = new TH3F(Form("hDCAxyVsPtProtonNeg_TOF"), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
	fListOfObjects->Add(hDCAxyVsPtPNeg_TOF);

	hDCAxyVsPtPiPos_TPC = new TH3F(Form("hDCAxyVsPtPionPos_TPC"), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
	fListOfObjects->Add(hDCAxyVsPtPiPos_TPC);
	hDCAxyVsPtPPos_TPC = new TH3F(Form("hDCAxyVsPtProtonPos_TPC"), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
	fListOfObjects->Add(hDCAxyVsPtPPos_TPC);
	hDCAxyVsPtPiPos_TOF = new TH3F(Form("hDCAxyVsPtPionPos_TOF"), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
	fListOfObjects->Add(hDCAxyVsPtPiPos_TOF);
	hDCAxyVsPtPPos_TOF = new TH3F(Form("hDCAxyVsPtProtonPos_TOF"), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
	fListOfObjects->Add(hDCAxyVsPtPPos_TOF);

	Int_t nPhiBins = 36;

	for(Int_t i = 0; i<nCent; ++i){

		hMIPVsEta[i] = new TH2D(Form("hMIPVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsEta[i] = new TProfile(Form("pMIPVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8, fDeDxMIPMin, fDeDxMIPMax);
		hMIPVsEtaV0s[i] = new TH2D(Form("hMIPVsEtaV0s%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, secondary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsEtaV0s[i] = new TProfile(Form("pMIPVsEtaV0s%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, secondary tracks}",50,-0.8,0.8,fDeDxMIPMin, fDeDxMIPMax);
		hPlateauVsEta[i] = new TH2D(Form("hPlateauVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{Plateau, primary tracks}",50,-0.8,0.8,50, 60, 110);
		pPlateauVsEta[i] = new TProfile(Form("pPlateauVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",50,-0.8,0.8, 60, 110);
		hPhi[i] = new TH2D(Form("histPhi%.2f-%.2f",CentMin[i],CentMax[i]), ";pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4);

		hPt_TPC[i]     = new TH1D(Form("hPt_rTPC_%.2f_%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtneg_TPC[i]  = new TH1D(Form("hPt_TPC_neg_%.2f_%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtpos_TPC[i]  = new TH1D(Form("hPt_TPC_pos_%.2f_%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtpos_TOF[i] = new TH1D(Form("hPt_TOF_pos_%.2f_%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtneg_TOF[i] = new TH1D(Form("hPt_TOF_neg_%.2f_%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);

		for(Int_t j=0; j<nHists; j++) {

			hDeDxVsP[i][j] = new TH2D( Form("hDeDxVsP%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p} [GeV/c]; dE/dx", nPtBins, ptBins, fdEdxHigh-fdEdxLow, fdEdxLow, fdEdxHigh);

			hnSigmaPiPos[i][j] = new TH2D(Form("hnSigma_Pion_pos_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};n#sigma", nPtBins, ptBins, 40, -10, 10);

			hnSigmaKPos[i][j] = new TH2D(Form("hnSigma_Kaon_pos_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};n#sigma", nPtBins, ptBins, 40, -10, 10);

			hnSigmaPPos[i][j] = new TH2D(Form("hnSigma_Proton_pos_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};n#sigma", nPtBins, ptBins, 40, -10, 10);

			hnSigmaPiNeg[i][j] = new TH2D(Form("hnSigma_Pion_neg_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};n#sigma", nPtBins, ptBins, 40, -10, 10);

			hnSigmaKNeg[i][j] = new TH2D(Form("hnSigma_Kaon_neg_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};n#sigma", nPtBins, ptBins, 40, -10, 10);

			hnSigmaPNeg[i][j] = new TH2D(Form("hnSigma_Proton_neg_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};n#sigma", nPtBins, ptBins, 40, -10, 10);

			hBetavsPpos[i][j] = new TH2D(Form("hBetavsP_pos_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};BetavsP", nPtBins, ptBins, 300, 0.2, 1.3);

			hBetavsPneg[i][j] = new TH2D(Form("hBetavsP_neg_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};BetavsP", nPtBins, ptBins, 300, 0.2, 1.3);

			hMIPVsPhi[i][j] = new TH2D(Form("hMIPVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);

			pMIPVsPhi[i][j] = new TProfile(Form("pMIPVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMin, fDeDxMIPMax);

			hPlateauVsPhi[i][j]  = new TH2D(Form("hPlateauVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),20, 70, 90);

			pPlateauVsPhi[i][j] = new TProfile(Form("pPlateauVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMax, 95);

			hPtpos_TPC_Eta[i][j] = new TH1D(Form("hPt_TPC_pos_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);

			hPtneg_TPC_Eta[i][j] = new TH1D(Form("hPt_TPC_neg_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);

			hPtpos_TOF_Eta[i][j] = new TH1D(Form("hPt_TOF_pos_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);

			hPtneg_TOF_Eta[i][j] = new TH1D(Form("hPt_TOF_neg_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);

			hPpos_TOF_Eta[i][j] = new TH1D(Form("hP_TOF_pos_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);

			hPneg_TOF_Eta[i][j] = new TH1D(Form("hP_TOF_neg_%.2f_%.2f_%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);

			hPtVsP[i][j] = new TH2D(Form("hPtVsP%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p} [GeV/c]; #it{p}_{T}", nPtBins, ptBins, nPtBins, ptBins);

			histPiV0[i][j]  = new TH2D(Form("hPiV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Pions id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);

			histpPiV0[i][j]  = new TH1D(Form("hPiV01D%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Pions id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);

			histPV0[i][j]   = new TH2D(Form("hPV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Protons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);

			histpPV0[i][j]  = new TH1D(Form("hPV01D%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Protons id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);

			histPiTof[i][j] = new TH2D(Form("hPiTOF%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Primary Pions from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);

			histpPiTof[i][j]  = new TH1D(Form("hPTOF%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Primary Pions from TOF ; #it{p} (GeV/#it{c}); counts", 200, 0, 20);

			histEV0[i][j]   = new TH2D(Form("hEV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Electrons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);

		}// eta loop
	} // centrality loop


	for(Int_t i = 0; i<nHists; ++i ){

		hMIPVsV0M[i] = new TH2D(Form("hMIPVsV0M-%s",ending[i]), Form("%s; V0M mult.; dE/dx MIP",LatexEta[i]), 100, 0, 100, fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsV0M[i] = new TProfile(Form("pMIPVsV0M-%s",ending[i]), Form("%s; V0M mult.; dE/dx MIP",LatexEta[i]), 100, 0, 100, fDeDxMIPMin, fDeDxMIPMax);
		hMIPVsNch[i] = new TH2D(Form("hMIPVsNch-%s",ending[i]),"; TPC track mult. |#eta|<0.8; dE/dx MIP", 100, 1, 101, fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsNch[i] = new TProfile(Form("pMIPVsNch-%s",ending[i]),"; TPC track mult. |#eta|<0.8; dE/dx MIP", 100, 1, 101, fDeDxMIPMin, fDeDxMIPMax);
	}



	if(!fAnalysisMC){
		/*for(Int_t j=0; j<nHists; ++j){
		  fListOfObjects->Add(hMIPVsV0M[j]);
		  fListOfObjects->Add(pMIPVsV0M[j]);
		  fListOfObjects->Add(hMIPVsNch[j]);
		  fListOfObjects->Add(pMIPVsNch[j]);
		  }*/

		for(Int_t i=0; i<nCent; ++i ){

			fListOfObjects->Add(hMIPVsEta[i]);
			fListOfObjects->Add(pMIPVsEta[i]);
			//			fListOfObjects->Add(hMIPVsEtaV0s[i]);
			//			fListOfObjects->Add(pMIPVsEtaV0s[i]);
			fListOfObjects->Add(hPlateauVsEta[i]);
			fListOfObjects->Add(pPlateauVsEta[i]);
			fListOfObjects->Add(hPhi[i]);

			/*for(Int_t j=0; j<nHists; ++j){
			  fListOfObjects->Add(hMIPVsPhi[i][j]);
			  fListOfObjects->Add(pMIPVsPhi[i][j]);
			  fListOfObjects->Add(hPlateauVsPhi[i][j]);
			  fListOfObjects->Add(pPlateauVsPhi[i][j]);
			  }*/

			////			if(fMakePid){
			fListOfObjects->Add(hPt_TPC[i]);	
			////				if(fLowPt){
			fListOfObjects->Add(hPtneg_TPC[i]);
			fListOfObjects->Add(hPtpos_TPC[i]);
			fListOfObjects->Add(hPtneg_TOF[i]);
			fListOfObjects->Add(hPtpos_TOF[i]);

			///				}

			for(Int_t j=0; j<nHists; ++j){
				//					if(fLowPt){
				fListOfObjects->Add(hnSigmaPiPos[i][j]);
				fListOfObjects->Add(hnSigmaPiNeg[i][j]);	
				fListOfObjects->Add(hnSigmaKPos[i][j]);
				fListOfObjects->Add(hnSigmaKNeg[i][j]);
				fListOfObjects->Add(hnSigmaPPos[i][j]);
				fListOfObjects->Add(hnSigmaPNeg[i][j]);
				fListOfObjects->Add(hBetavsPpos[i][j]);
				fListOfObjects->Add(hBetavsPneg[i][j]);
				fListOfObjects->Add(hPtpos_TPC_Eta[i][j]);
				fListOfObjects->Add(hPtneg_TPC_Eta[i][j]);
				fListOfObjects->Add(hPtpos_TOF_Eta[i][j]);
				fListOfObjects->Add(hPtneg_TOF_Eta[i][j]);
				fListOfObjects->Add(hPpos_TOF_Eta[i][j]);
				fListOfObjects->Add(hPneg_TOF_Eta[i][j]);
				///					}

				fListOfObjects->Add(hPtVsP[i][j]);
				fListOfObjects->Add(histPiV0[i][j]);
				//					fListOfObjects->Add(histpPiV0[i][j]);
				//					fListOfObjects->Add(histpPV0[i][j]);
				//					fListOfObjects->Add(histpPiTof[i][j]);
				fListOfObjects->Add(histPiTof[i][j]);
				fListOfObjects->Add(histEV0[i][j]);
				fListOfObjects->Add(histPV0[i][j]);
				fListOfObjects->Add(hDeDxVsP[i][j]);
			}
			////			} //	if(MakePID) 
		} //	Cent
	} //	!fAnalysisMC


	else{
		for(Int_t pid=0; pid<3; pid++) {
			if(pid==0){
				hMcInNeg[pid] = new TH1D("hInNeg_Pion","",nPtBins,ptBins);
				hMcInPos[pid] = new TH1D("hInPos_Pion","",nPtBins,ptBins);
				hMcOutNeg[pid] = new TH1D("hPtRecNegTPC_Pion","",nPtBins,ptBins);
				hMcOutPos[pid] = new TH1D("hPtRecPosTPC_Pion","",nPtBins,ptBins);
				hMcOutNegTOF[pid] = new TH1D("hPtRecNegTOF_Pion","",nPtBins,ptBins);
				hMcOutPosTOF[pid] = new TH1D("hPtRecPosTOF_Pion","",nPtBins,ptBins);
			}
			else if(pid==1){
				hMcInNeg[pid] = new TH1D("hInNeg_Kaon","",nPtBins,ptBins);
				hMcInPos[pid] = new TH1D("hInPos_Kaon","",nPtBins,ptBins);
				hMcOutNeg[pid] = new TH1D("hPtRecNegTPC_Kaon","",nPtBins,ptBins);
				hMcOutPos[pid] = new TH1D("hPtRecPosTPC_Kaon","",nPtBins,ptBins);
				hMcOutNegTOF[pid] = new TH1D("hPtRecNegTOF_Kaon","",nPtBins,ptBins);
				hMcOutPosTOF[pid] = new TH1D("hPtRecPosTOF_Kaon","",nPtBins,ptBins);
			}
			else{
				hMcInNeg[pid] = new TH1D("hInNeg_Proton","",nPtBins,ptBins);
				hMcInPos[pid] = new TH1D("hInPos_Proton","",nPtBins,ptBins);
				hMcOutNeg[pid] = new TH1D("hPtRecNegTPC_Proton","",nPtBins,ptBins);
				hMcOutPos[pid] = new TH1D("hPtRecPosTPC_Proton","",nPtBins,ptBins);
				hMcOutNegTOF[pid] = new TH1D("hPtRecNegTOF_Proton","",nPtBins,ptBins);
				hMcOutPosTOF[pid] = new TH1D("hPtRecPosTOF_Proton","",nPtBins,ptBins);
			}

			fListOfObjects->Add(hMcInNeg[pid]);
			fListOfObjects->Add(hMcInPos[pid]);
			fListOfObjects->Add(hMcOutNeg[pid]);
			fListOfObjects->Add(hMcOutPos[pid]);
			fListOfObjects->Add(hMcOutNegTOF[pid]);
			fListOfObjects->Add(hMcOutPosTOF[pid]);

		}	// pid Eff

		for(Int_t pid=0; pid<2; ++pid){
			for(Int_t q=0; q<2; ++q){
				if(pid==0){
					hDCApTPrim[pid][q] = new TH3F(Form("hDCA_Pion%sPrim_TPC",Q[q]),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
					hDCApTWDec[pid][q] = new TH3F(Form("hDCA_Pion%sWDec_TPC",Q[q]),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
					hDCApTMate[pid][q] = new TH3F(Form("hDCA_Pion%sMate_TPC",Q[q]),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);

					hDCApTPrim_TOF[pid][q] = new TH3F(Form("hDCA_Pion%sPrim_TOF",Q[q]),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
					hDCApTWDec_TOF[pid][q] = new TH3F(Form("hDCA_Pion%sWDec_TOF",Q[q]),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
					hDCApTMate_TOF[pid][q] = new TH3F(Form("hDCA_Pion%sMate_TOF",Q[q]),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
				}else{
					hDCApTPrim[pid][q] = new TH3F(Form("hDCA_Proton%sPrim_TPC",Q[q]),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
					hDCApTWDec[pid][q] = new TH3F(Form("hDCA_Proton%sWDec_TPC",Q[q]),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
					hDCApTMate[pid][q] = new TH3F(Form("hDCA_Proton%sMate_TPC",Q[q]),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);

					hDCApTPrim_TOF[pid][q] = new TH3F(Form("hDCA_Proton%sPrim_TOF",Q[q]),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
					hDCApTWDec_TOF[pid][q] = new TH3F(Form("hDCA_Proton%sWDec_TOF",Q[q]),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);
					hDCApTMate_TOF[pid][q] = new TH3F(Form("hDCA_Proton%sMate_TOF",Q[q]),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)",nPtBins,ptBins,nDCABins,DCABins,nV0MMultBins,V0MMultBins);

				}

				fListOfObjects->Add(hDCApTPrim[pid][q]);
				fListOfObjects->Add(hDCApTWDec[pid][q]);
				fListOfObjects->Add(hDCApTMate[pid][q]);

				fListOfObjects->Add(hDCApTPrim_TOF[pid][q]);
				fListOfObjects->Add(hDCApTWDec_TOF[pid][q]);
				fListOfObjects->Add(hDCApTMate_TOF[pid][q]);
			}	// charge MC
		}	// pid MC

		fVtxMC = new TH1F("fVtxMC","mc vtx", 120, -30, 30);
		fListOfObjects->Add(fVtxMC);

	}


	//fEventCuts.AddQAplotsToList(fListOfObjects);
	// Post output data.
	PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskPPvsMult::UserExec(Option_t *)
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

	if (fAnalysisType == "ESD"){
		fESD = dynamic_cast<AliESDEvent*>(event);
		if(!fESD){
			Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}

	} else{
		fAOD = dynamic_cast<AliAODEvent*>(event);
		if(!fAOD){
			Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}
	}

	if (fAnalysisMC){
		if (fAnalysisType == "ESD"){
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
		} else { // AOD

			fMC = dynamic_cast<AliMCEvent*>(MCEvent());
			if(fMC)
				fMC->Dump();

			fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
			if(!fMCArray){
				Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}
		}
	}

	utils = new AliAnalysisUtils();
	if (!utils)
	{
		cout<<"------- No AnalysisUtils Object Found --------"<<utils<<endl;
		return;
	}

	printf("flag0\n");
	fEvents->Fill(0.5,11);

	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;
	fEvents->Fill(1.5,11);

	printf("flag1\n");
	Float_t V0MPer  = -1;
	Float_t MultBin = -1;

	AliMultSelection *MultSelection = (AliMultSelection*) fESD -> FindListObject("MultSelection");
	if(MultSelection-> IsEventSelected()){
		V0MPer = MultSelection->GetMultiplicityPercentile("V0M",kFALSE);
	}

	if((V0MPer>100)||(V0MPer<0))return;
	for(Int_t icent = 0; icent < (nCent-1); ++icent){
		if(V0MPer > CentMin[icent] && V0MPer <= CentMax[icent]){
			MultBin = icent+1;
		}
	}

	printf("flag2\n");
	Int_t INEL = -1;
	INEL = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0);
	if( INEL < 1 )
		return;

	printf("flag3\n");

	fEvents->Fill(2.5,MultBin);
	fEvents->Fill(2.5,11);

	if( utils->IsSPDClusterVsTrackletBG(fESD) ) 
		return;

	printf("flag4\n");

	fEvents->Fill(3.5,MultBin);
	fEvents->Fill(3.5,11);

	if( fESD->IsPileupFromSPDInMultBins() )
		return;

	printf("flag5\n");

	fEvents->Fill(4.5,MultBin);
	fEvents->Fill(4.5,11);

	if( fESD->IsIncompleteDAQ()) 
		return;

	fEvents->Fill(5.5,MultBin);
	fEvents->Fill(5.5,11);


	if( !selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE) ){return;}
	fEvents->Fill(6.5,MultBin);
	fEvents->Fill(6.5,11);

	if( !IsGoodZvertexPos(fESD) ){return;}
	fEvents->Fill(7.5,MultBin);
	fEvents->Fill(7.5,11);

	fV0M->Fill(V0MPer);

	//	cout << "INEL ============  " << INEL << endl;


	//  	const AliESDVertex * trkVertex = fESD->GetPrimaryVertexTracks();
	//  	const AliESDVertex * spdVertex = fESD->GetPrimaryVertexSPD();
	//  	Bool_t hasSPD = spdVertex->GetStatus();
	//  	Bool_t hasTrk = trkVertex->GetStatus();

	//   	Double_t cov[6]={0};
	//   	spdVertex->GetCovarianceMatrix(cov);
	//   	Double_t zRes = TMath::Sqrt(cov[5]);


	/*

	   if (!hasTrk){
	   if (!hasSPD){ return; } 
	   else{
	   if ( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ){ return; }
	   }
	   }
	   else {
	   if(hasTrk){
	   if ( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ){ return; }
	   }
	   }

	   fEvents->Fill(6.5,MultBin);
	   fEvents->Fill(6.5,11);

	   if(!(TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())<=0.5))return;
	   fEvents->Fill(7.5,MultBin);
	   fEvents->Fill(7.5,11);


	   if (!hasTrk){
	   if (!hasSPD){ return; } 
	   else{
	   if ( !(TMath::Abs(spdVertex->GetZ())<fVtxCut)  ){ return; }
	   }
	   }
	   else {
	   if(hasTrk){
	   if ( !(TMath::Abs(trkVertex->GetZ())<fVtxCut)  ){ return; }
	   }
	   }

	   fEvents->Fill(8.5,MultBin);
	   fEvents->Fill(8.5,11);
	   */



	if (fAnalysisType == "ESD")
		AnalyzeESD(fESD);

	if (fAnalysisMC) {
		if (fAnalysisType == "ESD"){

			AliHeader* headerMC = fMC->Header();
			if (headerMC) {

				AliGenEventHeader* genHeader = headerMC->GenEventHeader();
				TArrayF vtxMC(3); // primary vertex  MC
				vtxMC[0]=9999; vtxMC[1]=9999;  vtxMC[2]=9999; //initialize with dummy
				if (genHeader) {
					genHeader->PrimaryVertex(vtxMC);
				}
				fZvtxMC = vtxMC[2];

				// PYTHIA:
				AliGenPythiaEventHeader* pythiaGenHeader =
					dynamic_cast<AliGenPythiaEventHeader*>(headerMC->GenEventHeader());
				if (pythiaGenHeader) {  //works only for pythia
					fMcProcessType =  GetPythiaEventProcessType(pythiaGenHeader->ProcessType());
				}
				// PHOJET:
				AliGenDPMjetEventHeader* dpmJetGenHeader =
					dynamic_cast<AliGenDPMjetEventHeader*>(headerMC->GenEventHeader());
				if (dpmJetGenHeader) {
					fMcProcessType = GetDPMjetEventProcessType(dpmJetGenHeader->ProcessType());
				}
			}
		} else { // AOD
			AliAODMCHeader* mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader"));
			if(mcHeader) {
				fZvtxMC = mcHeader->GetVtxZ();
				if(strstr(mcHeader->GetGeneratorName(), "Pythia")) {
					fMcProcessType =  GetPythiaEventProcessType(mcHeader->GetEventType());
				} else {
					fMcProcessType =  GetDPMjetEventProcessType(mcHeader->GetEventType());
				}
			}
		}
	}


	PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskPPvsMult::AnalyzeESD(AliESDEvent* esdEvent)
{

	//	Float_t centrality = -10;
	//	if(fAnalysisPbPb){
	//		centrality = fEventCuts.GetCentrality(); /// Centrality calculated with the default estimator (V0M for LHC15o)

	//		if((centrality>fMaxCent)||(centrality<fMinCent))return;
	//		for(Int_t icent = 0; icent < nCent; ++icent){
	//			if(centrality > CentMin[icent] && centrality <= CentMax[icent]){
	//				cent = icent;
	//				fcent->Fill(icent+1);
	//				ProduceArrayTrksESD( esdEvent, cent );
	//				ProduceArrayV0ESD( esdEvent, cent );

	//				if(fAnalysisMC)
	//					ProcessMCTruthESD(cent);
	//			}
	//		}

	//		fcent->Fill(11);

	//	}

	//	else{

	Float_t V0MPercentile = -1;
	Float_t V0APercentile = -1;
	Float_t ADMPercentile = -1;

	AliMultSelection *MultSelection = (AliMultSelection*) esdEvent -> FindListObject("MultSelection");
	if(MultSelection-> IsEventSelected()){
		V0MPercentile = MultSelection->GetMultiplicityPercentile("V0M",false);
		V0APercentile = MultSelection->GetMultiplicityPercentile("V0A",false);
		ADMPercentile = MultSelection->GetMultiplicityPercentile("ADM",false);
	}

	if((V0MPercentile>100)||(V0MPercentile<0))return;
	for(Int_t icent = 0; icent < (nCent-1); ++icent){
		if(V0MPercentile > CentMin[icent] && V0MPercentile <= CentMax[icent]){
			cent = icent;
			fcent->Fill(icent+1);
			ProduceArrayTrksESD( esdEvent, cent );
			ProduceArrayV0ESD( esdEvent, cent );

			if(fAnalysisMC)
				ProcessMCTruthESD(cent);
		}
	}

	fcent->Fill(11);


	//	}



}


//_____________________________________________________________________________
Float_t AliAnalysisTaskPPvsMult::GetVertex(const AliVEvent* event) const
{
	Float_t zvtx = -999;

	const AliVVertex* primaryVertex = event->GetPrimaryVertex();

	if(primaryVertex->GetNContributors()>0)
		zvtx = primaryVertex->GetZ();

	return zvtx;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskPPvsMult::GetPidCode(Int_t pdgCode) const
{
	// return our internal code for pions, kaons, and protons

	Short_t pidCode = 6;

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
void AliAnalysisTaskPPvsMult::ProcessMCTruthESD(const Int_t Cent)
{
	// Fill the special MC histogram with the MC truth info

	///cout<<"Cent Inside ProcessMCTruth ::: "<<Cent<<endl;
	const Int_t nTracksMC = fMCStack->GetNtrack();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		TParticle* trackMC = fMCStack->Particle(iTracks);
		if( !trackMC )
			continue;

		if( !(fMCStack->IsPhysicalPrimary(iTracks)) )
			continue;

		TParticlePDG* pdgPart = trackMC ->GetPDG();
		Double_t chargeMC = pdgPart->Charge();

		if(chargeMC==0)
			continue;

		if ( TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		if( chargeMC < 0.0 ){
			if(pidCodeMC==1)hMcInNeg[0]->Fill(trackMC->Pt());
			if(pidCodeMC==2)hMcInNeg[1]->Fill(trackMC->Pt());
			if(pidCodeMC==3)hMcInNeg[2]->Fill(trackMC->Pt());
		}
		else{
			if(pidCodeMC==1)hMcInPos[0]->Fill(trackMC->Pt());
			if(pidCodeMC==2)hMcInPos[1]->Fill(trackMC->Pt());
			if(pidCodeMC==3)hMcInPos[2]->Fill(trackMC->Pt());
		}


	}//MC track loop



}

//_____________________________________________________________________________

Short_t AliAnalysisTaskPPvsMult::GetPythiaEventProcessType(Int_t pythiaType) {
	//
	// Get the process type of the event.  PYTHIA
	//
	// source PWG0   dNdpt

	Short_t globalType = -1; //init

	if(pythiaType==92||pythiaType==93){
		globalType = 2; //single diffractive
	}
	else if(pythiaType==94){
		globalType = 3; //double diffractive
	}
	//else if(pythiaType != 91){ // also exclude elastic to be sure... CKB??
	else {
		globalType = 1;  //non diffractive
	}
	return globalType;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskPPvsMult::GetDPMjetEventProcessType(Int_t dpmJetType) {
	//
	// get the process type of the event.  PHOJET
	//
	//source PWG0   dNdpt
	// can only read pythia headers, either directly or from cocktalil header
	Short_t globalType = -1;

	if (dpmJetType == 1 || dpmJetType == 4) { // explicitly inelastic plus central diffraction
		globalType = 1;
	}
	else if (dpmJetType==5 || dpmJetType==6) {
		globalType = 2;
	}
	else if (dpmJetType==7) {
		globalType = 3;
	}
	return globalType;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskPPvsMult::GetEventIdAsLong(AliVHeader* header) const
{
	// To have a unique id for each event in a run!
	// Modified from AliRawReader.h
	return ((ULong64_t)header->GetBunchCrossNumber()+
			(ULong64_t)header->GetOrbitNumber()*3564+
			(ULong64_t)header->GetPeriodNumber()*16777215*3564);
}


//____________________________________________________________________
TParticle* AliAnalysisTaskPPvsMult::FindPrimaryMother(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	Int_t motherLabel = FindPrimaryMotherLabel(stack, label);
	if (motherLabel < 0)
		return 0;

	return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskPPvsMult::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// returns its label
	//
	// Taken from AliPWG0Helper class
	//
	const Int_t nPrim  = stack->GetNprimary();

	while (label >= nPrim) {

		//printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

		TParticle* particle = stack->Particle(label);
		if (!particle) {

			AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
			return -1;
		}

		// find mother
		if (particle->GetMother(0) < 0) {

			AliDebugGeneral("FindPrimaryMotherLabel", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
			return -1;
		}

		label = particle->GetMother(0);
	}

	return label;
}

//____________________________________________________________________
TParticle* AliAnalysisTaskPPvsMult::FindPrimaryMotherV0(AliStack* stack, Int_t label)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	Int_t nSteps = 0;

	Int_t motherLabel = FindPrimaryMotherLabelV0(stack, label, nSteps);
	if (motherLabel < 0)
		return 0;

	return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskPPvsMult::FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// returns its label
	//
	// Taken from AliPWG0Helper class
	//
	nSteps = 0;
	const Int_t nPrim  = stack->GetNprimary();

	while (label >= nPrim) {

		//printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

		nSteps++; // 1 level down

		TParticle* particle = stack->Particle(label);
		if (!particle) {

			AliDebugGeneral("FindPrimaryMotherLabelV0", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
			return -1;
		}

		// find mother
		if (particle->GetMother(0) < 0) {

			AliDebugGeneral("FindPrimaryMotherLabelV0", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
			return -1;
		}

		label = particle->GetMother(0);
	}

	return label;
}
//__________________________________________________________________
void AliAnalysisTaskPPvsMult::ProduceArrayTrksESD( AliESDEvent *ESDevent, const Int_t Cent ){

	const Int_t nESDTracks = ESDevent->GetNumberOfTracks();

	fcentAfterPrimaries->Fill(Cent+1);
	fcentAfterPrimaries->Fill(11);

	Float_t V0MPer  = -1;

	AliMultSelection *MultSelection = (AliMultSelection*)ESDevent -> FindListObject("MultSelection");
	if(MultSelection-> IsEventSelected())
		V0MPer = MultSelection->GetMultiplicityPercentile("V0M",false);


	Int_t multTPC = 0;
	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);

		if( TMath::Abs(esdTrack->Eta()) > fEtaCut )
			continue;

		UInt_t selectDebug = 0;
		if (fTrackFilterTPC) {
			selectDebug = fTrackFilterTPC->IsSelected(esdTrack);
			if (!selectDebug) {
				continue;
			}
		}

		multTPC++;

	}

	for(Int_t iT = 0; iT < nESDTracks; iT++){

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
		if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		if ( fTrackFilterGolden ) {
			if (!fTrackFilterGolden->AcceptTrack(esdTrack)) {
				continue;
			}
		}

		if(esdTrack->GetTPCsignalN() < fNcl) continue;

		Double_t eta      = esdTrack->Eta();
		Double_t phi      = esdTrack->Phi();
		Double_t momentum = esdTrack->P();
		Double_t pt       = esdTrack->Pt();
		Float_t  dedx     = esdTrack->GetTPCsignal();

		double nSigmaPionTPC = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion);
		double nSigmaPionTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion);
		double nSigmaProtonTPC = fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton);
		double nSigmaProtonTOF = fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton);

		Float_t dcaxy = 0.;
		Float_t dcaz = 0.;
		esdTrack->GetImpactParameters(dcaxy,dcaz);

		Short_t pidCode     = 0;
		if(fAnalysisMC) {

			const Int_t label = TMath::Abs(esdTrack->GetLabel());
			TParticle* mcTrack = 0;
			mcTrack = fMCStack->Particle(label);

			if (!mcTrack) continue; ///{

			if( esdTrack->Charge()==0 )
				continue;

			Int_t pdgCode = mcTrack->GetPdgCode();
			pidCode = GetPidCode(pdgCode);

			if( fMCStack->IsPhysicalPrimary(label) ){

				if( esdTrack->Charge() < 0.0 ){

					if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) )  {

						if(PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), Magf, fcutLow, fcutHigh)){
							if(pidCode==1) hMcOutNeg[0]->Fill(esdTrack->Pt());
							if(pidCode==2) hMcOutNeg[1]->Fill(esdTrack->Pt());
							if(pidCode==3) hMcOutNeg[2]->Fill(esdTrack->Pt());
						}
					}

					if(pidCode==1) hDCApTPrim[0][0]->Fill(pt,dcaxy,V0MPer);
					if(pidCode==3) hDCApTPrim[1][0]->Fill(pt,dcaxy,V0MPer);

				}
				else{

					if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) ){

						if(PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), Magf, fcutLow, fcutHigh)){
							if(pidCode==1) hMcOutPos[0]->Fill(esdTrack->Pt());
							if(pidCode==2) hMcOutPos[1]->Fill(esdTrack->Pt());
							if(pidCode==3) hMcOutPos[2]->Fill(esdTrack->Pt());
						}

					}

					if(pidCode==1) hDCApTPrim[0][1]->Fill(pt,dcaxy,V0MPer);
					if(pidCode==3) hDCApTPrim[1][1]->Fill(pt,dcaxy,V0MPer);

				}
			}	// Primary particles MC

			if( fMCStack->IsSecondaryFromWeakDecay(label) ){

				if( esdTrack->Charge() < 0.0 ){
					if(pidCode==1) hDCApTWDec[0][0]->Fill(pt,dcaxy,V0MPer);
					if(pidCode==3) hDCApTWDec[1][0]->Fill(pt,dcaxy,V0MPer);
				}
				else{
					if(pidCode==1) hDCApTWDec[0][1]->Fill(pt,dcaxy,V0MPer);
					if(pidCode==3) hDCApTWDec[1][1]->Fill(pt,dcaxy,V0MPer);
				}
			}	// Weak Decay MC

			if( fMCStack->IsSecondaryFromMaterial(label) ){

				if( esdTrack->Charge() < 0.0 ){
					if(pidCode==1) hDCApTMate[0][0]->Fill(pt,dcaxy,V0MPer);
					if(pidCode==3) hDCApTMate[1][0]->Fill(pt,dcaxy,V0MPer);
				}
				else{
					if(pidCode==1) hDCApTMate[0][1]->Fill(pt,dcaxy,V0MPer);
					if(pidCode==3) hDCApTMate[1][1]->Fill(pt,dcaxy,V0MPer);
				}
			}	// Material Inte MC

			bool TOF_out = kFALSE;
			TOF_out = TOFPID(esdTrack);
			if(TOF_out){
				if( fMCStack->IsPhysicalPrimary(label) ){

					if( esdTrack->Charge() < 0.0 ){

						if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) )  {
							if(pidCode==1) hMcOutNegTOF[0]->Fill(esdTrack->Pt());
							if(pidCode==2) hMcOutNegTOF[1]->Fill(esdTrack->Pt());
							if(pidCode==3) hMcOutNegTOF[2]->Fill(esdTrack->Pt());
						}

						if(pidCode==1) hDCApTPrim_TOF[0][0]->Fill(pt,dcaxy,V0MPer);
						if(pidCode==3) hDCApTPrim_TOF[1][0]->Fill(pt,dcaxy,V0MPer);
					}
					else{

						if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) )  {
							if(pidCode==1) hMcOutPosTOF[0]->Fill(esdTrack->Pt());
							if(pidCode==2) hMcOutPosTOF[1]->Fill(esdTrack->Pt());
							if(pidCode==3) hMcOutPosTOF[2]->Fill(esdTrack->Pt());
						}

						if(pidCode==1) hDCApTPrim_TOF[0][1]->Fill(pt,dcaxy,V0MPer);
						if(pidCode==3) hDCApTPrim_TOF[1][1]->Fill(pt,dcaxy,V0MPer);
					}
				}	// Primary particles MC

				if( fMCStack->IsSecondaryFromWeakDecay(label) ){

					if( esdTrack->Charge() < 0.0 ){
						if(pidCode==1) hDCApTWDec_TOF[0][0]->Fill(pt,dcaxy,V0MPer);
						if(pidCode==3) hDCApTWDec_TOF[1][0]->Fill(pt,dcaxy,V0MPer);
					}
					else{
						if(pidCode==1) hDCApTWDec_TOF[0][1]->Fill(pt,dcaxy,V0MPer);
						if(pidCode==3) hDCApTWDec_TOF[1][1]->Fill(pt,dcaxy,V0MPer);
					}
				}	// Weak Decay MC

				if( fMCStack->IsSecondaryFromMaterial(label) ){

					if( esdTrack->Charge() < 0.0 ){
						if(pidCode==1) hDCApTMate_TOF[0][0]->Fill(pt,dcaxy,V0MPer);
						if(pidCode==3) hDCApTMate_TOF[1][0]->Fill(pt,dcaxy,V0MPer);
					}
					else{
						if(pidCode==1) hDCApTMate_TOF[0][1]->Fill(pt,dcaxy,V0MPer);
						if(pidCode==3) hDCApTMate_TOF[1][1]->Fill(pt,dcaxy,V0MPer);
					}
				}	// Material Inte MC

			} // DCA for TOF analysis

			///			}	//mcTrack
		}	//fAnalysis MC

		if(esdTrack->Charge() < 0.0){

			if(TMath::Abs(nSigmaPionTPC)<2.0)
				hDCAxyVsPtPiNeg_TPC->Fill(pt,dcaxy,V0MPer);

			if(TMath::Abs(nSigmaProtonTPC)<2.0)
				hDCAxyVsPtPNeg_TPC->Fill(pt,dcaxy,V0MPer);

		}else{

			if(TMath::Abs(nSigmaPionTPC)<2.0)
				hDCAxyVsPtPiPos_TPC->Fill(pt,dcaxy,V0MPer);

			if(TMath::Abs(nSigmaProtonTPC)<2.0)
				hDCAxyVsPtPPos_TPC->Fill(pt,dcaxy,V0MPer);
		}

		bool IsTOFin = kFALSE;
		IsTOFin = TOFPID(esdTrack);
		if(IsTOFin){
			if(esdTrack->Charge() < 0.0){

				if(TMath::Sqrt(TMath::Power(nSigmaPionTPC,2.0)+TMath::Power(nSigmaPionTOF,2.0))<2.0)
					hDCAxyVsPtPiNeg_TOF->Fill(pt,dcaxy,V0MPer);

				if(TMath::Sqrt(TMath::Power(nSigmaProtonTPC,2.0)+TMath::Power(nSigmaProtonTOF,2.0))<2.0)
					hDCAxyVsPtPNeg_TOF->Fill(pt,dcaxy,V0MPer);

			}else{
				if(TMath::Sqrt(TMath::Power(nSigmaPionTPC,2.0)+TMath::Power(nSigmaPionTOF,2.0))<2.0)
					hDCAxyVsPtPiPos_TOF->Fill(pt,dcaxy,V0MPer);

				if(TMath::Sqrt(TMath::Power(nSigmaProtonTPC,2.0)+TMath::Power(nSigmaProtonTOF,2.0))<2.0)
					hDCAxyVsPtPPos_TOF->Fill(pt,dcaxy,V0MPer);

			} // Negative charge tracks

		} // DCA for TOF analysis

		// DCAxy cut
		if( TMath::Abs(dcaxy) > GetMaxDCApTDep(fcutDCAxy,pt) )
			continue;

		int nh = -1;        
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

		if(esdTrack->Charge() < 0.0){
			hPtneg_TPC[Cent]->Fill(pt);
			hPtneg_TPC[10]->Fill(pt);
			hPtneg_TPC_Eta[Cent][nh]->Fill(pt);
			hPtneg_TPC_Eta[10][nh]->Fill(pt);
			hnSigmaPiNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaPiNeg[10][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaKNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaKNeg[10][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaPNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
			hnSigmaPNeg[10][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));

		}else{
			hPtpos_TPC[Cent]->Fill(pt);
			hPtpos_TPC[10]->Fill(pt);
			hPtpos_TPC_Eta[Cent][nh]->Fill(pt);
			hPtpos_TPC_Eta[10][nh]->Fill(pt);
			hnSigmaPiPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaPiPos[10][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaKPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaKPos[10][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaPPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
			hnSigmaPPos[10][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));

		}

		if(IsTOFin){

			double trkLength = esdTrack->GetIntegratedLength();
			double Beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);

			if(esdTrack->Charge() < 0.0){
				hBetavsPneg[Cent][nh]->Fill(esdTrack->P(),Beta);
				hBetavsPneg[10][nh]->Fill(esdTrack->P(),Beta);
				hPneg_TOF_Eta[Cent][nh]->Fill(esdTrack->P());
				hPneg_TOF_Eta[10][nh]->Fill(esdTrack->P());
				hPtneg_TOF_Eta[Cent][nh]->Fill(pt);
				hPtneg_TOF_Eta[10][nh]->Fill(pt);
				hPtneg_TOF[Cent]->Fill(pt);
				hPtneg_TOF[10]->Fill(pt);

			}else{
				hBetavsPpos[Cent][nh]->Fill(esdTrack->P(),Beta);
				hBetavsPpos[10][nh]->Fill(esdTrack->P(),Beta);
				hPpos_TOF_Eta[Cent][nh]->Fill(esdTrack->P());
				hPpos_TOF_Eta[10][nh]->Fill(esdTrack->P());
				hPtpos_TOF_Eta[Cent][nh]->Fill(pt);
				hPtpos_TOF_Eta[10][nh]->Fill(pt);
				hPtpos_TOF[Cent]->Fill(pt);
				hPtpos_TOF[10]->Fill(pt);

			} // Negative charge tracks

		} // TOF analysis

		//TOF
		Bool_t IsTOFout=kFALSE;
		if ((esdTrack->GetStatus()&AliESDtrack::kTOFout)==0)
			IsTOFout=kTRUE;
		Float_t lengthtrack=esdTrack->GetIntegratedLength();
		Float_t timeTOF=esdTrack->GetTOFsignal();
		Double_t inttime[5]={0,0,0,0,0};
		esdTrack->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis
		Float_t beta = -99;
		if ( !IsTOFout ){
			if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
				beta = inttime[0] / timeTOF;
		}

		hPtVsP[Cent][nh]->Fill(momentum,pt);
		hPtVsP[10][nh]->Fill(momentum,pt);

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), Magf, fcutLow, fcutHigh))
			continue;

		if(fdEdxCalibrated)
			dedx *= 50.0/EtaCalibration(eta);

		if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
			if( esdTrack->GetTPCsignal() < fDeDxMIPMax && esdTrack->GetTPCsignal() > fDeDxMIPMin ){
				hMIPVsEta[Cent]->Fill(eta,dedx);
				pMIPVsEta[Cent]->Fill(eta,dedx);
				hMIPVsEta[10]->Fill(eta,dedx);
				pMIPVsEta[10]->Fill(eta,dedx);
			}
			if( esdTrack->GetTPCsignal() > 70.0 && esdTrack->GetTPCsignal() < 90.0 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta[Cent]->Fill(eta,dedx);
					pPlateauVsEta[Cent]->Fill(eta,dedx);
					hPlateauVsEta[10]->Fill(eta,dedx);
					pPlateauVsEta[10]->Fill(eta,dedx);
				}
			}
		}

		if(beta>1){
			histPiTof[Cent][nh]->Fill(momentum, dedx);
			histPiTof[10][nh]->Fill(momentum, dedx);
		}

		if( momentum <= 0.6 && momentum >= 0.4  ){

			if( esdTrack->GetTPCsignal() < fDeDxMIPMax && esdTrack->GetTPCsignal() > fDeDxMIPMin ){
				hMIPVsPhi[Cent][nh]->Fill(phi,dedx);
				pMIPVsPhi[Cent][nh]->Fill(phi,dedx);
				hMIPVsPhi[10][nh]->Fill(phi,dedx);
				pMIPVsPhi[10][nh]->Fill(phi,dedx);

				hMIPVsNch[nh]->Fill(multTPC,dedx);
				pMIPVsNch[nh]->Fill(multTPC,dedx);

				hMIPVsV0M[nh]->Fill(V0MPer,dedx);
				pMIPVsV0M[nh]->Fill(V0MPer,dedx);

			}
			if( esdTrack->GetTPCsignal() > 70.0 && esdTrack->GetTPCsignal() < 90.0 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsPhi[Cent][nh]->Fill(phi,dedx);
					pPlateauVsPhi[Cent][nh]->Fill(phi,dedx);
					hPlateauVsPhi[10][nh]->Fill(phi,dedx);
					pPlateauVsPhi[10][nh]->Fill(phi,dedx);

				}
			}
		}

		hPt_TPC[Cent]->Fill(pt);
		hPt_TPC[10]->Fill(pt);
		hDeDxVsP[Cent][nh]->Fill(momentum,dedx);
		hDeDxVsP[10][nh]->Fill(momentum,dedx);


	}//end of track loop
}
//________________________________________________________________________
void AliAnalysisTaskPPvsMult::ProduceArrayV0ESD( AliESDEvent *ESDevent, const Int_t Cent ){

	Int_t nv0s = ESDevent->GetNumberOfV0s();

	fcentAfterV0s->Fill(Cent+1);
	fcentAfterV0s->Fill(11);

	const AliESDVertex *myBestPrimaryVertex = ESDevent->GetPrimaryVertex();

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

		AliESDv0 *esdV0 = ESDevent->GetV0(iV0);
		if ( !esdV0 ) continue;

		//check onfly status
		//		if( !esdV0->GetOnFlyStatus() )
		//			continue;

		if( esdV0->GetOnFlyStatus()!=0 )
			continue;

		// AliESDTrack (V0 Daughters)
		UInt_t lKeyPos = (UInt_t)TMath::Abs(esdV0->GetPindex());
		UInt_t lKeyNeg = (UInt_t)TMath::Abs(esdV0->GetNindex());

		AliESDtrack *pTrack = ESDevent->GetTrack(lKeyPos);
		AliESDtrack *nTrack = ESDevent->GetTrack(lKeyNeg);

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

		if ( fTrackFilterGolden ) {

			if (!fTrackFilterGolden->AcceptTrack(pTrack)) continue;

			if (!fTrackFilterGolden->AcceptTrack(nTrack)) continue;

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
		AliKFParticle::SetField(ESDevent->GetMagneticField());

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

						       if(track->GetTPCsignalN()<fNcl)continue;
						       Double_t phi     = track->Phi();

						       if(!PhiCut(track->Pt(), phi, track->Charge(), Magf, fcutLow, fcutHigh))
							       continue;

						       Double_t eta      = track->Eta();
						       Double_t momentum = track->P();
						       Double_t dedx     = track->GetTPCsignal();
						       Double_t dedxUnc  = track->GetTPCsignal();

						       if(fdEdxCalibrated)
							       dedx *= 50.0/EtaCalibration(eta);

						       if(fillPos&&fillNeg){
							       if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s[Cent]->Fill(eta,dedx);
									       pMIPVsEtaV0s[Cent]->Fill(eta,dedx);
									       hMIPVsEtaV0s[10]->Fill(eta,dedx);
									       pMIPVsEtaV0s[10]->Fill(eta,dedx);
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

							       histPiV0[Cent][nh]->Fill(momentum, dedx);
							       //histpPiV0[Cent][nh]->Fill(momentum);
							       histPiV0[10][nh]->Fill(momentum, dedx);
							       //histpPiV0[Cent][nh]->Fill(momentum);

						       }
						       else{
							       histPV0[Cent][nh]->Fill(momentum, dedx);
							       histPV0[10][nh]->Fill(momentum, dedx);
							       //histpPV0[Cent][nh]->Fill(momentum);

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

					       if(track->GetTPCsignalN()<fNcl)
						       continue;

					       if(!PhiCut(track->Pt(), phi, track->Charge(), Magf, fcutLow, fcutHigh))
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

					       //					       cout<<"Cent :: "<<Cent<<"Eta  :: "<<eta<<"dedx :: "<<dedx<<endl;
					       histEV0[Cent][nh]->Fill(momentum, dedx);
					       histEV0[10][nh]->Fill(momentum, dedx);

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
Bool_t AliAnalysisTaskPPvsMult::selectVertex2015pp(AliESDEvent *esd,
		Bool_t checkSPDres, //enable check on vtx resolution
		Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
		Bool_t checkProximity) //apply cut on relative position of spd and trk verteces
{

	if (!esd) return kFALSE;

	const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
	const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
	Bool_t hasSPD = spdVertex->GetStatus();
	Bool_t hasTrk = trkVertex->GetStatus();

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
Bool_t AliAnalysisTaskPPvsMult::IsGoodSPDvertexRes(const AliESDVertex* spdVertex)
{

	if( !spdVertex ) return kFALSE;
	if( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPPvsMult::IsGoodZvertexPos(AliESDEvent *esd)
{

	if( !esd ) return kFALSE;
	//Cut on the vertex z position
	const AliESDVertex * vertex = esd->GetPrimaryVertex();
	if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPPvsMult::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh)

{
	if(pt < 2.0)
		return kTRUE;

	//Double_t phi = track->Phi();
	if(fESD->GetMagneticField() < 0)    // for negatve polarity field
		phi = TMath::TwoPi() - phi;
	if(q < 0) // for negatve charge
		phi = TMath::TwoPi()-phi;

	phi += TMath::Pi()/18.0; // to center gap in the middle
	phi = fmod(phi, TMath::Pi()/9.0);

	if(phi<phiCutHigh->Eval(pt)
			&& phi>phiCutLow->Eval(pt))
		return kFALSE; // reject track

	hPhi[cent]->Fill(pt, phi);
	hPhi[10]->Fill(pt, phi);

	return kTRUE;
}
//________________________________________________________________________
Float_t AliAnalysisTaskPPvsMult::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){

	Double_t maxDCAxy = 10;
	maxDCAxy = fMaxDCAxy->Eval(ptI);
	return maxDCAxy;

}
//________________________________________________________________________
double AliAnalysisTaskPPvsMult::EtaCalibration(const double &eta){

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

	if(strcmp(fPeriod,"16l")==0){
		aPos = 49.9216; bPos = 3.07252; cPos = -42.8044; dPos = 259.666; ePos = -910.432; fPos = 1776.09; gPos = -1740.65; hPos = 662.232;
		aNeg = 49.9732; bNeg = 4.03575; cNeg = 65.6189;  dNeg = 374.429; eNeg = 951.459;  fNeg = 1153.75; gNeg = 618.493;  hNeg = 100.499;
	}else if(strcmp(fPeriod,"16k")==0){
		aPos = 49.9421; bPos = 2.3446; cPos = -41.2765; dPos = 279.695; ePos = -1027.73; fPos = 2022.84; gPos = -1967.79; hPos = 738.823;
		aNeg = 50.0477; bNeg = 8.27344; cNeg = 125.29;  dNeg = 736.8;   eNeg = 2057.75;  fNeg = 2935.38; gNeg = 2064.03;  hNeg = 565.983;
	}else if(strcmp(fPeriod,"16deghijop")==0){
		aPos = 49.9743; bPos = 2.3388; cPos = -44.1496; dPos = 296.029; ePos = -1056.56; fPos = 2031.44; gPos = -1946.51; hPos = 723.89;
		aNeg = 50.0329; bNeg = 6.99747; cNeg = 107.168;  dNeg = 649.001; eNeg = 1875.17;  fNeg = 2785.78; gNeg = 2063.77;  hNeg = 606.868;
	}else if(strcmp(fPeriod,"17data")==0){
		aPos = 49.6097; bPos = 0.922856; cPos = -6.57484; dPos = 65.3117; ePos = -372.142; fPos = 950.451; gPos = -1085.27; hPos = 458.144;
		aNeg = 49.6555; bNeg = 6.98696; cNeg = 102.734;  dNeg = 566.424; eNeg = 1513.64;  fNeg = 2092.01; gNeg = 1429.32;  hNeg = 375.642;
	}else{
		aPos = 49.6975; bPos = 2.32535; cPos = -42.6516; dPos = 283.058; ePos = -1009.58; fPos = 1945.89; gPos = -1871.23; hPos = 698.552;
		aNeg = 49.8071; bNeg = 9.78466; cNeg = 120.018;  dNeg = 603.325; eNeg = 1470.92;  fNeg = 1819.63; gNeg = 1073.82;  hNeg = 230.142;
	}

	for(int i=0; i<8; ++i){
		fEtaCalibrationNeg->SetParameter(i,0);
		fEtaCalibrationPos->SetParameter(i,0);
	}

	if(eta<0.0){
		fEtaCalibrationNeg->SetParameter(0,aNeg);
		fEtaCalibrationNeg->SetParameter(1,bNeg);
		fEtaCalibrationNeg->SetParameter(2,cNeg);
		fEtaCalibrationNeg->SetParameter(3,dNeg);
		fEtaCalibrationNeg->SetParameter(4,eNeg);
		fEtaCalibrationNeg->SetParameter(5,fNeg);
		fEtaCalibrationNeg->SetParameter(6,gNeg);
		fEtaCalibrationNeg->SetParameter(7,hNeg);

		return fEtaCalibrationNeg->Eval(eta);
	}
	else{
		fEtaCalibrationPos->SetParameter(0,aPos);
		fEtaCalibrationPos->SetParameter(1,bPos);
		fEtaCalibrationPos->SetParameter(2,cPos);
		fEtaCalibrationPos->SetParameter(3,dPos);
		fEtaCalibrationPos->SetParameter(4,ePos);
		fEtaCalibrationPos->SetParameter(5,fPos);
		fEtaCalibrationPos->SetParameter(6,gPos);
		fEtaCalibrationPos->SetParameter(7,hPos);

		return fEtaCalibrationPos->Eval(eta);
	}


}
//________________________________________________________________________
double AliAnalysisTaskPPvsMult::EtaCalibrationEl(const double &eta){

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

	if(strcmp(fPeriod,"16l")==0){
		aPos = 79.4195; bPos = 7.82459; cPos = -23.3466; dPos = 26.5577; ePos = -8.27151;
		aNeg = 79.8571; bNeg = -14.2921; cNeg = -66.6972; dNeg = -103.794; eNeg = -50.5771;
	}else if(strcmp(fPeriod,"16k")==0){
		aPos = 80.254; bPos = 6.37076; cPos = -50.9878; dPos = 116.611; ePos = -79.0483;
		aNeg = 79.8728; bNeg = -3.08265; cNeg = -11.3778; dNeg = -20.6605; eNeg = -12.3861;
	}else if(strcmp(fPeriod,"16deghijop")==0){
		aPos = 80.0719; bPos = 7.10053; cPos = -42.4788; dPos = 86.1074; ePos = -54.0891;
		aNeg = 79.6155; bNeg = -12.1254; cNeg = -66.2488; dNeg = -132.426; eNeg = -85.0155;
	}else if(strcmp(fPeriod,"17data")==0){
		aPos = 82.4621; bPos = 5.20353; cPos = -32.2608; dPos = 63.4788; ePos = -39.3277;
		aNeg = 82.306; bNeg = -4.04076; cNeg = -22.133; dNeg = -40.5782; eNeg = -23.8157;
	}else{
		aPos = 79.7726; bPos = 6.83744; cPos = -40.0469; dPos = 78.987; ePos = -50.1373;
		aNeg = 79.4863; bNeg = -5.00403; cNeg = -21.6184;  dNeg = -39.1295; eNeg = -24.8757;
	}

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

		return felededxfitNeg->Eval(eta);
	}
	else{
		felededxfitPos->SetParameter(0,aPos);
		felededxfitPos->SetParameter(1,bPos);
		felededxfitPos->SetParameter(2,cPos);
		felededxfitPos->SetParameter(3,dPos);
		felededxfitPos->SetParameter(4,ePos);

		return felededxfitPos->Eval(eta);
	}


}
//________________________________________________________________________
/*Double_t AliAnalysisTaskPPvsMult::EtaCalibrationNeg( const Int_t Cent, const Double_t eta){


  for(Int_t i=0; i<8; ++i)
  fEtaCalibrationNeg->SetParameter(i,0);

  fEtaCalibrationNeg->SetParameter(0,aNeg[Cent]);
  fEtaCalibrationNeg->SetParameter(1,bNeg[Cent]);
  fEtaCalibrationNeg->SetParameter(2,cNeg[Cent]);
  fEtaCalibrationNeg->SetParameter(3,dNeg[Cent]);
  fEtaCalibrationNeg->SetParameter(4,eNeg[Cent]);
  fEtaCalibrationNeg->SetParameter(5,fNeg[Cent]);
  fEtaCalibrationNeg->SetParameter(6,gNeg[Cent]);
  fEtaCalibrationNeg->SetParameter(7,hNeg[Cent]);

  return fEtaCalibrationNeg->Eval(eta);


  }
//________________________________________________________________________
Double_t AliAnalysisTaskPPvsMult::EtaCalibrationPos( const Int_t Cent, const Double_t eta){


for(Int_t i=0; i<8; ++i)
fEtaCalibration->SetParameter(i,0);

fEtaCalibration->SetParameter(0,aPos[Cent]);
fEtaCalibration->SetParameter(1,bPos[Cent]);
fEtaCalibration->SetParameter(2,cPos[Cent]);
fEtaCalibration->SetParameter(3,dPos[Cent]);
fEtaCalibration->SetParameter(4,ePos[Cent]);
fEtaCalibration->SetParameter(5,fPos[Cent]);
fEtaCalibration->SetParameter(6,gPos[Cent]);
fEtaCalibration->SetParameter(7,hPos[Cent]);

return fEtaCalibration->Eval(eta);

}
//________________________________________________________________________
Double_t AliAnalysisTaskPPvsMult::EtaCalibrationNegEl(const Int_t Cent, const Double_t eta){


for(Int_t i=0; i<5; ++i)
felededxfitNeg->SetParameter(i,0);


felededxfitNeg->SetParameter(0,aNegEl[Cent]);
felededxfitNeg->SetParameter(1,bNegEl[Cent]);
felededxfitNeg->SetParameter(2,cNegEl[Cent]);
felededxfitNeg->SetParameter(3,dNegEl[Cent]);
felededxfitNeg->SetParameter(4,eNegEl[Cent]);


return felededxfitNeg->Eval(eta);

}
//________________________________________________________________________
Double_t AliAnalysisTaskPPvsMult::EtaCalibrationPosEl(const Int_t Cent, const Double_t eta){


for(Int_t i=0; i<5; ++i)
felededxfitPos->SetParameter(i,0);

felededxfitPos->SetParameter(0,aPosEl[Cent]);
felededxfitPos->SetParameter(1,bPosEl[Cent]);
felededxfitPos->SetParameter(2,cPosEl[Cent]);
felededxfitPos->SetParameter(3,dPosEl[Cent]);
felededxfitPos->SetParameter(4,ePosEl[Cent]);

return felededxfitPos->Eval(eta);

}*/
//________________________________________________________________________
bool AliAnalysisTaskPPvsMult::TOFPID(AliESDtrack * track)
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
