/*************************************************************************
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
 * provided "as is" without express or implied warranty.                  * **************************************************************************/


#include "AliAnalysisTaskQAHighPtDeDxTest.h"

// ROOT includes
#include <TList.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
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
#include "AliMultSelection.h"
#include <AliPIDResponse.h>
#include "AliTPCPIDResponse.h"

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




const Double_t AliAnalysisTaskQAHighPtDeDxTest::fgkClight = 2.99792458e-2;
Float_t magf = -1;
TF1* cutLow  = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
TF1* cutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);
Double_t DeDxMIPMin  = 45;
Double_t DeDxMIPMax  = 55;
const Int_t nHists = 4;
Float_t centralityGlobal = -10;
Int_t etaLow[nHists+5]  = {-8, -8, -6, -4, -2, 0, 2, 4, 6};
Int_t etaHigh[nHists+5] = { 8, -6, -4, -2,  0, 2, 4, 6, 8};

Int_t EtaLow[nHists]  = {6, 4, 2, 0};
Int_t EtaHigh[nHists] = {8, 6, 4, 2};
const Int_t nCent = 10;
const Double_t CentMin[nCent] = {0.0,5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0};
const Double_t CentMax[nCent] = {5.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};

const Double_t aPos[nCent] = {49.8594  ,49.8542  ,49.8594  ,49.8545  ,49.8473   ,49.8589  ,49.8542   ,49.8569  ,49.8461   ,49.8512};
const Double_t bPos[nCent] = {-0.19737 ,-0.157422,-0.185804,-0.140537,-0.0502213,-0.155494,-0.0575582,-0.149527,-0.0901908,-0.0571605};
const Double_t cPos[nCent] = {0.936921 ,0.840229 ,0.856222 ,0.789579 ,0.508899  ,0.782057 ,0.577245  ,0.788935 ,0.688114  ,0.591065};
const Double_t dPos[nCent] = {-0.387464,-0.312826,-0.302875,-0.260128,-0.0490301 ,-0.23182,-0.125936 ,-0.262908,-0.211082 ,-0.152547};

const Double_t aNeg[nCent] = {49.8292  ,49.8493 ,49.8388  ,49.8413   ,49.8384  ,49.8436   ,49.8476  ,49.8336   ,49.8438 ,49.8446};
const Double_t bNeg[nCent] = {-0.122292,0.11876 ,0.0060856,-0.0200841,0.0228002,0.00971686,0.0146653,-0.0481027,0.077598,-0.0379178};
const Double_t cNeg[nCent] = {0.748659 ,1.39121 ,1.05511  ,0.972488  ,1.09636  ,1.04173   ,1.06027  ,0.924111  ,1.30225 ,0.852475};
const Double_t dNeg[nCent] = {0.415989 ,0.894481,0.617493 ,0.550473  ,0.645483 ,0.599914  ,0.615319 ,0.522093  ,0.840085,0.424588};

const Double_t aPosEl[nCent] = {78.7401 ,78.583   ,78.6193 ,78.5886,78.8275 ,78.6249 ,78.6817  ,78.5971,78.433  ,78.433};
const Double_t bPosEl[nCent] = {2.98839 ,2.04009  ,3.36522 ,2.58982,0.595512,3.50901 ,1.66237  ,2.80526,4.38372 ,4.38372 };
const Double_t cPosEl[nCent] = {-7.83643,-0.972868,-7.10595,-5.1087,0.45588 ,-9.23427,-0.980859,-6.0588,-9.11186,-9.11186 };
const Double_t dPosEl[nCent] = {4.1571  ,-2.85664 ,2.94706 ,2.01079,-2.98313,5.55299 ,-2.39829 ,2.8525 ,4.2765  ,4.2765};

const Double_t aNegEl[nCent]={78.5643  ,78.6004  ,78.2853 ,78.3675 ,78.3426 ,78.5809    ,78.2088 ,78.3638 ,78.1972,78.1972};
const Double_t bNegEl[nCent]={-1.62327 ,-0.662359,-3.78407,-3.31185,-3.7269 ,-1.19542  ,-5.41726,-3.44326,-4.20979,-4.20979};
const Double_t cNegEl[nCent]={-2.77851 ,1.53235  ,-7.07994,-6.34969,-7.29292,-0.00395374,-11.227,-6.04591,-6.53041,-6.53041};
const Double_t dNegEl[nCent]={-0.787795,3.57675  ,-3.2976 ,-2.87936,-3.32353,2.74549    ,-6.12438,-1.991  ,-1.802  ,-1.802};

const Bool_t CloseDCAxy = kTRUE;
Int_t nDeltaPiBins   = 80;
Double_t deltaPiLow  = 20;
Double_t deltaPiHigh = 100;
const Double_t dEdxHigh = 200;
const Double_t dEdxLow  = 40;
const Char_t *Pid[7]={"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};
ClassImp(AliAnalysisTaskQAHighPtDeDxTest)
	//_____________________________________________________________________________
	//AliAnalysisTaskQAHighPtDeDx::AliAnalysisTaskQAHighPtDeDx(const char *name):
	AliAnalysisTaskQAHighPtDeDxTest::AliAnalysisTaskQAHighPtDeDxTest():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fAOD(0x0),
		fEventCuts(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fMCArray(0x0),
		fPIDResponse(0x0),
		fTrackFilter2015PbPb(0x0),
		fTrackFilterGolden(0x0),
		fTrackFilterTPC(0x0),
		fCentEst("V0M"),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fAnalysisPbPb(kFALSE),
		ftrigBit(0x0),
		fRandom(0x0),
		fPileUpRej(kFALSE),
		fEtaCut(0.9),
		cent(0x0),
		fMinCent(0.0),
		fMaxCent(100.0),
		fStoreMcIn(kFALSE),//
		fMcProcessType(-999),
		fTriggeredEventMB(-999),
		fVtxStatus(-999),
		fZvtx(-999),
		fZvtxMC(-999),
		fRun(-999),
		fEventId(-999),
		fListOfObjects(0),
		fVtxMC(0x0),
		fdEdxCalibrated(0x0),
		fMakePid(0x0),
		fcent(0x0),
		fEtaCalibrationNeg(0x0),
		fEtaCalibration(0x0),
		felededxfitPos(0x0),
		felededxfitNeg(0x0),
		fcutDCAxy(0x0)


{


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

		hPtAll[i]=0;
		hPtAllPos[i]=0;
		hPtAllNeg[i]=0;

		hDCAxyVsPtPiNeg[i]=0;
		hDCAxyVsPtKNeg[i]=0;
		hDCAxyVsPtPNeg[i]=0;
		hDCAxyVsPtPiNegC[i]=0;
		hDCAxyVsPtKNegC[i]=0;
		hDCAxyVsPtPNegC[i]=0;
		hDCAxyVsPtPiPos[i]=0;
		hDCAxyVsPtKPos[i]=0;
		hDCAxyVsPtPPos[i]=0;
		hDCAxyVsPtPiPosC[i]=0;
		hDCAxyVsPtKPosC[i]=0;
		hDCAxyVsPtPPosC[i]=0;

		for(Int_t j=0; j<nHists; ++j){

			hPtPos[i][j]=0;//TH1D, Transverse momentum distribution  positive-charged hadrons
			hPtNeg[i][j]=0;//TH1D, Transverse momentum distribution  negative-charged hadrons
			hPtVsP[i][j]=0;//TH2D, Transverse momentum Vs momentum

			hMIPVsNch[i][j]=0;//TH2D, MIP vs Nch for different eta intervals
			pMIPVsNch[i][j]=0;//TProfile, MIP vs Nch for different eta intervals
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

			histPiV0[i][j]=0;//TH2D, dE/dx vs p, pi id by V0s
			histpPiV0[i][j]=0;//TH1D, pi id by V0s
			histPV0[i][j]=0;// TH2D, dE/dx vs p, p id by V0s
			histpPV0[i][j]=0;// TH1D, p id by V0s
			histPiTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
			histElTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
			histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
			histEV0[i][j]=0;

		}

	}


	//default constructor
	for(Int_t i=0;i<9;++i){

		for(Int_t pid=0;pid<7;++pid){
			hMcIn[pid][i]=0;
			hMcOut[pid][i]=0;
		}

	}

}


AliAnalysisTaskQAHighPtDeDxTest::AliAnalysisTaskQAHighPtDeDxTest(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fAOD(0x0),
	fEventCuts(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fMCArray(0x0),
	fPIDResponse(0x0),
	fTrackFilter2015PbPb(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilterTPC(0x0),
	fCentEst("V0M"),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fAnalysisPbPb(kFALSE),
	ftrigBit(0x0),
	fRandom(0x0),
	fPileUpRej(kFALSE),
	fEtaCut(0.9),
	cent(0x0),
	fMinCent(0.0),
	fMaxCent(100.0),
	fStoreMcIn(kFALSE),
	fMcProcessType(-999),
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtxMC(-999),
	fRun(-999),
	fEventId(-999),
	fListOfObjects(0), 
	fVtxMC(0x0),
	fdEdxCalibrated(0x0),
	fMakePid(0x0),
	fcent(0x0),
	fEtaCalibrationNeg(0x0),
	fEtaCalibration(0x0),
	felededxfitPos(0x0),
	felededxfitNeg(0x0),
	fcutDCAxy(0x0)

{


	for(Int_t i = 0; i<nCent; ++i){

		hMIPVsEta[i] = 0;
		hMIPVsEta[i]=0;
		pMIPVsEta[i]=0;
		hMIPVsEtaV0s[i]=0;
		pMIPVsEtaV0s[i]=0;
		hPlateauVsEta[i]=0;
		pPlateauVsEta[i]=0;
		hPhi[i]=0;

		hPtAll[i]=0;
		hPtAllPos[i]=0;
		hPtAllNeg[i]=0;

		hDCAxyVsPtPiNeg[i]=0;
		hDCAxyVsPtKNeg[i]=0;
		hDCAxyVsPtPNeg[i]=0;
		hDCAxyVsPtPiNegC[i]=0;
		hDCAxyVsPtKNegC[i]=0;
		hDCAxyVsPtPNegC[i]=0;
		hDCAxyVsPtPiPos[i]=0;
		hDCAxyVsPtKPos[i]=0;
		hDCAxyVsPtPPos[i]=0;
		hDCAxyVsPtPiPosC[i]=0;
		hDCAxyVsPtKPosC[i]=0;
		hDCAxyVsPtPPosC[i]=0;

		for(Int_t j=0; j<nHists; ++j){

			hPtPos[i][j]=0;//TH1D, Transverse momentum distribution  positive-charged hadrons
			hPtNeg[i][j]=0;//TH1D, Transverse momentum distribution  negative-charged hadrons
			hPtVsP[i][j]=0;//TH2D, Transverse momentum Vs momentum

			hMIPVsNch[i][j]=0;//TH2D, MIP vs Nch for different eta intervals
			pMIPVsNch[i][j]=0;//TProfile, MIP vs Nch for different eta intervals
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

			histPiV0[i][j]=0;//TH2D, dE/dx vs p, pi id by V0s
			histpPiV0[i][j]=0;//TH1D, pi id by V0s
			histPV0[i][j]=0;// TH2D, dE/dx vs p, p id by V0s
			histpPV0[i][j]=0;// TH1D, p id by V0s
			histPiTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
			histElTof[i][j]=0;//TH2D, dE/dx vs p for a "clean" sample of pions, beta>1
			histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
			histEV0[i][j]=0;

		}

	}


	// Default constructor (should not be used)
	for(Int_t i=0;i<9;++i){

		for(Int_t pid=0;pid<7;++pid){
			hMcIn[pid][i]=0;
			hMcOut[pid][i]=0;
		}


	}
	DefineOutput(1, TList::Class());//esto es nuevo
}




AliAnalysisTaskQAHighPtDeDxTest::~AliAnalysisTaskQAHighPtDeDxTest() {
	//
	// Destructor
	//

}



//______________________________________________________________________________



void AliAnalysisTaskQAHighPtDeDxTest::UserCreateOutputObjects()
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


	//OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner();

	//
	// Histograms
	//

	fcent=new TH1F("fcent","fcent",13,0,13);
	fListOfObjects->Add(fcent);

	const Int_t nPtBins = 63;
	Double_t ptBins[nPtBins+1] = {

		0.01,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
		0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,
		1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
		4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,30

	};


	const Int_t nPtBinsV0s = 25;
	Double_t ptBinsV0s[nPtBinsV0s+1] = { 0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
		1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
		9.0 , 12.0, 15.0, 20.0 };

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


	const Char_t* ending[nHists] = {"02", "24", "46", "68"};
	const Char_t* LatexEta[nHists] = {"|#eta|<0.2", "0.2<|#eta|<0.4", "0.4<|#eta|<0.6", "0.6<|#eta|<0.8" };

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);


	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol3", -1.0, 0.0);
	fEtaCalibration    = new TF1("fDeDxVsEtaPos", "pol3", 0.0, 1.0);

	felededxfitPos     = new TF1("felededxfitPos", "pol3", 0.0, 1.0);
	felededxfitNeg     = new TF1("felededxfitNeg", "pol3", -1.0, 0.0);


	Int_t nPhiBins = 36;

	for(Int_t i = 0; i<nCent; ++i){

		hMIPVsEta[i] = new TH2D(Form("hMIPVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, primary tracks}",16,-0.8,0.8,DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
		pMIPVsEta[i] = new TProfile(Form("pMIPVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, primary tracks}",16,-0.8,0.8, DeDxMIPMin, DeDxMIPMax);
		hMIPVsEtaV0s[i] = new TH2D(Form("hMIPVsEtaV0s%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, secondary tracks}",16,-0.8,0.8,DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
		pMIPVsEtaV0s[i] = new TProfile(Form("pMIPVsEtaV0s%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, secondary tracks}",16,-0.8,0.8,DeDxMIPMin, DeDxMIPMax);
		hPlateauVsEta[i] = new TH2D(Form("hPlateauVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{Plateau, primary tracks}",16,-0.8,0.8,20, 70, 90);
		pPlateauVsEta[i] = new TProfile(Form("pPlateauVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",16,-0.8,0.8, DeDxMIPMax, 95);
		hPhi[i] = new TH2D(Form("histPhi%.2f-%.2f",CentMin[i],CentMax[i]), ";pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4);

		hPtAll[i] = new TH1D(Form("hPt%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAll[i]->Sumw2();
		hPtAllNeg[i] = new TH1D(Form("hPtNeg%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAllNeg[i]->Sumw2();
		hPtAllPos[i] = new TH1D(Form("hPtPos%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAllPos[i]->Sumw2();

		hDCAxyVsPtPiNeg[i] = new TH2D(Form("hDCAxyVsPtPiNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPiNegC[i] = new TH2D(Form("hDCAxyVsPtPiNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPiNeg[i]->Sumw2();
		hDCAxyVsPtPiNegC[i]->Sumw2();
		hDCAxyVsPtKNeg[i] = new TH2D(Form("hDCAxyVsPtKNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtKNegC[i] = new TH2D(Form("hDCAxyVsPtKNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtKNeg[i]->Sumw2();
		hDCAxyVsPtKNegC[i]->Sumw2();
		hDCAxyVsPtPNeg[i] = new TH2D(Form("hDCAxyVsPtPNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPNegC[i] = new TH2D(Form("hDCAxyVsPtPNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPNeg[i]->Sumw2();
		hDCAxyVsPtPNegC[i]->Sumw2();

		hDCAxyVsPtPiPos[i] = new TH2D(Form("hDCAxyVsPtPiPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPiPosC[i] = new TH2D(Form("hDCAxyVsPtPiPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPiPos[i]->Sumw2();
		hDCAxyVsPtPiPosC[i]->Sumw2();
		hDCAxyVsPtKPos[i] = new TH2D(Form("hDCAxyVsPtKPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtKPosC[i] = new TH2D(Form("hDCAxyVsPtKPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtKPos[i]->Sumw2();
		hDCAxyVsPtKPosC[i]->Sumw2();
		hDCAxyVsPtPPos[i] = new TH2D(Form("hDCAxyVsPtPPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPPosC[i] = new TH2D(Form("hDCAxyVsPtPPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPPos[i]->Sumw2();
		hDCAxyVsPtPPosC[i]->Sumw2();

		for(Int_t j=0; j<nHists; j++) {

			hDeDxVsP[i][j] = new TH2D( Form("hDeDxVsP%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p} [GeV/c]; dE/dx", nPtBins, ptBins, dEdxHigh-dEdxLow, dEdxLow, dEdxHigh);
			hDeDxVsP[i][j]->Sumw2();

			hnSigmaPiPos[i][j] = new TH2D(Form("hnSigmaPiPos%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPiPos", nPtBins, ptBins, 40, -10, 10);
			hnSigmaPiPos[i][j]->Sumw2();

			hnSigmaKPos[i][j] = new TH2D(Form("hnSigmaKPos%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaKPos", nPtBins, ptBins, 40, -10, 10);
			hnSigmaKPos[i][j]->Sumw2();

			hnSigmaPPos[i][j] = new TH2D(Form("hnSigmaPPos%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPPos", nPtBins, ptBins, 40, -10, 10);
			hnSigmaPPos[i][j]->Sumw2();

			hnSigmaPiNeg[i][j] = new TH2D(Form("hnSigmaPiNeg%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPiNeg", nPtBins, ptBins, 40, -10, 10);
			hnSigmaPiNeg[i][j]->Sumw2();

			hnSigmaKNeg[i][j] = new TH2D(Form("hnSigmaKNeg%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaKNeg", nPtBins, ptBins, 40, -10, 10);
			hnSigmaKNeg[i][j]->Sumw2();

			hnSigmaPNeg[i][j] = new TH2D(Form("hnSigmaPNeg%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p}_{T};nSigmaPNeg", nPtBins, ptBins, 40, -10, 10);
			hnSigmaPNeg[i][j]->Sumw2();

			hMIPVsPhi[i][j] = new TH2D(Form("hMIPVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
			hMIPVsPhi[i][j]->Sumw2();

			pMIPVsPhi[i][j] = new TProfile(Form("pMIPVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),DeDxMIPMin, DeDxMIPMax);
			pMIPVsPhi[i][j]->Sumw2();

			hPlateauVsPhi[i][j]  = new TH2D(Form("hPlateauVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),20, 70, 90);
			hPlateauVsPhi[i][j]->Sumw2();

			pPlateauVsPhi[i][j] = new TProfile(Form("pPlateauVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),DeDxMIPMax, 95);
			pPlateauVsPhi[i][j]->Sumw2();

			hMIPVsNch[i][j] = new TH2D(Form("hMIPVsNch%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; TPC track mult. |#eta|<0.8; dE/dx MIP",LatexEta[j]), 800, 1, 4001, DeDxMIPMax-DeDxMIPMin, DeDxMIPMin, DeDxMIPMax);
			hMIPVsNch[i][j]->Sumw2();

			pMIPVsNch[i][j] = new TProfile(Form("pMIPVsNch%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; TPC track mult. |#eta|<0.8; dE/dx MIP",LatexEta[j]), 800, 1, 4001, DeDxMIPMin, DeDxMIPMax);
			pMIPVsNch[i][j]->Sumw2();

			hPtPos[i][j] = new TH1D(Form("hPtPos%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);
			hPtPos[i][j]->Sumw2();

			hPtNeg[i][j] = new TH1D(Form("hPtNeg%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins);
			hPtNeg[i][j]->Sumw2();

			hPtVsP[i][j] = new TH2D(Form("hPtVsP%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p} [GeV/c]; #it{p}_{T}", nPtBins, ptBins, nPtBins, ptBins);
			hPtVsP[i][j]->Sumw2();

			histPiV0[i][j]  = new TH2D(Form("hPiV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Pions id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histPiV0[i][j]->Sumw2();

			histpPiV0[i][j]  = new TH1D(Form("hPiV01D%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Pions id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
			histpPiV0[i][j]->Sumw2();

			histPV0[i][j]   = new TH2D(Form("hPV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Protons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histPV0[i][j]->Sumw2();

			histpPV0[i][j]  = new TH1D(Form("hPV01D%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Protons id by V0; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
			histpPV0[i][j]->Sumw2();

			histPiTof[i][j] = new TH2D(Form("hPiTOF%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Primary Pions from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histPiTof[i][j]->Sumw2();

			histElTof[i][j] = new TH2D(Form("hElTOF%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Primary Electrons from TOF; #it{p} (GeV/#it{c}); d#it{e}d#it{x}", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histElTof[i][j]->Sumw2();

			histpPiTof[i][j]  = new TH1D(Form("hPTOF%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Primary Pions from TOF ; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
			histpPiTof[i][j]->Sumw2();

			histEV0[i][j]   = new TH2D(Form("hEV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Electrons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
			histEV0[i][j]->Sumw2();

			//					fListOfObjects->Add(hMIPVsPhi[i][j]);
			//					fListOfObjects->Add(pMIPVsPhi[i][j]);
			//					fListOfObjects->Add(hPlateauVsPhi[i][j]);
			//					fListOfObjects->Add(pPlateauVsPhi[i][j]);

			if(fMakePid){
				fListOfObjects->Add(hDeDxVsP[i][j]);
				//				fListOfObjects->Add(hnSigmaPiPos[i][j]);
				//				fListOfObjects->Add(hnSigmaPiNeg[i][j]);
				//				fListOfObjects->Add(hnSigmaKPos[i][j]);
				//				fListOfObjects->Add(hnSigmaKNeg[i][j]);
				//				fListOfObjects->Add(hnSigmaPPos[i][j]);
				//				fListOfObjects->Add(hnSigmaPPos[i][j]);
				//				fListOfObjects->Add(hPtPos[i][j]);
				//				fListOfObjects->Add(hPtNeg[i][j]);
				fListOfObjects->Add(hPtVsP[i][j]);

				fListOfObjects->Add(histPiV0[i][j]);
				//				fListOfObjects->Add(histpPiV0[i][j]);
				//				fListOfObjects->Add(histpPV0[i][j]);
				fListOfObjects->Add(histPV0[i][j]);
				//				fListOfObjects->Add(histpPiTof[i][j]);
				fListOfObjects->Add(histPiTof[i][j]);
				fListOfObjects->Add(histElTof[i][j]);
				fListOfObjects->Add(histEV0[i][j]);
			}

			//				fListOfObjects->Add(hMIPVsNch[i]);
			//				fListOfObjects->Add(pMIPVsNch[i]);

		}// eta loop


		fListOfObjects->Add(hMIPVsEta[i]);
		fListOfObjects->Add(pMIPVsEta[i]);
		fListOfObjects->Add(hMIPVsEtaV0s[i]);
		fListOfObjects->Add(pMIPVsEtaV0s[i]);
		fListOfObjects->Add(hPlateauVsEta[i]);
		fListOfObjects->Add(pPlateauVsEta[i]);
		//				fListOfObjects->Add(hPhi[i]);

		if(fMakePid){

			fListOfObjects->Add(hPtAll[i]);
			fListOfObjects->Add(hPtAllNeg[i]);
			fListOfObjects->Add(hPtAllPos[i]);

			//			fListOfObjects->Add(hDCAxyVsPtPiNeg[i]);
			//			fListOfObjects->Add(hDCAxyVsPtPiPos[i]);
			//			fListOfObjects->Add(hDCAxyVsPtKNeg[i]);
			//			fListOfObjects->Add(hDCAxyVsPtKPos[i]);
			//			fListOfObjects->Add(hDCAxyVsPtPNeg[i]);
			//			fListOfObjects->Add(hDCAxyVsPtPPos[i]);

			//			fListOfObjects->Add(hDCAxyVsPtPiNegC[i]);
			//			fListOfObjects->Add(hDCAxyVsPtPiPosC[i]);
			//			fListOfObjects->Add(hDCAxyVsPtKNegC[i]);
			//			fListOfObjects->Add(hDCAxyVsPtKPosC[i]);
			//			fListOfObjects->Add(hDCAxyVsPtPNegC[i]);
			//			fListOfObjects->Add(hDCAxyVsPtPPosC[i]);


		}



	} // centrality loop


	if (fAnalysisMC) {
		for(Int_t i = 0; i < nHists; i++) {
			for(Int_t pid = 0; pid < 7; pid++) {

				hMcIn[pid][i] = new TH1D(Form("hIn%s%s", Pid[pid],ending[i]), Form("MC in (pid %s)", Pid[pid]),
						nPtBinsV0s, ptBinsV0s);
				fListOfObjects->Add(hMcIn[pid][i]);

				hMcOut[pid][i] = new TH1D(Form("hMcOut%s%s", Pid[pid],ending[i]), Form("MC out (pid %s)", Pid[pid]),
						nPtBinsV0s, ptBinsV0s);
				fListOfObjects->Add(hMcOut[pid][i]);


			}
		}

		fVtxMC = new TH1F("fVtxMC","mc vtx", 120, -30, 30);
		fListOfObjects->Add(fVtxMC);

	}

	fEventCuts.AddQAplotsToList(fListOfObjects);
	// Post output data.
	PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::UserExec(Option_t *)
{
	// Main loop

	//
	// First we make sure that we have valid input(s)!
	//

	AliVEvent *event = InputEvent();

	//	hEvents->Fill(0);
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


	//--------------- Event Selection --------------------


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
void AliAnalysisTaskQAHighPtDeDxTest::AnalyzeESD(AliESDEvent* esdEvent)
{


	if (!fEventCuts.AcceptEvent(esdEvent)){
		PostData(1, fListOfObjects);
		return;
	}
	else{
		ProduceArrayTrksESD( esdEvent );
		ProduceArrayV0ESD( esdEvent );
	}

	Float_t centrality = -10;

	if(fAnalysisPbPb){
		AliMultSelection *MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
		if(fCentEst == "V0M")
			centrality = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
		if(fCentEst == "V0A")
			centrality = MultSelection->GetMultiplicityPercentile("V0A",kTRUE);

		if((centrality>fMaxCent)||(centrality<fMinCent))return;
		for(Int_t icent = 0; icent < nCent; ++icent){
			if(centrality > CentMin[icent] && centrality <= CentMax[icent]){
				cent = icent;
				fcent->Fill(icent+1);
			}
		}

		fcent->Fill(11);

	}



}

//________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::AnalyzeAOD(AliAODEvent* aodEvent)
{
	fRun  = aodEvent->GetRunNumber();
	fEventId = 0;
	if(aodEvent->GetHeader())
		fEventId = GetEventIdAsLong(aodEvent->GetHeader());

	//UInt_t    time      = 0; // Missing AOD info? aodEvent->GetTimeStamp();
	magf      = aodEvent->GetMagneticField();

	//Int_t     trackmult = 0; // no pt cuts
	//Int_t     nadded    = 0;

	Bool_t isPileup = aodEvent->IsPileupFromSPD();
	if(fPileUpRej)
		if(isPileup)
			return;
	//fn1->Fill(4);



	if(fTriggeredEventMB) {// Only MC case can we have not triggered events

		// accepted event
		//fEvents->Fill(0);

		//if(fVtxStatus!=1) return; // accepted vertex
		//Int_t nAODTracks = aodEvent->GetNumberOfTracks();

		ProduceArrayTrksAOD( aodEvent );
		ProduceArrayV0AOD( aodEvent );

		//fEvents->Fill(1);




	} // end if triggered

}

//_____________________________________________________________________________
Float_t AliAnalysisTaskQAHighPtDeDxTest::GetVertex(const AliVEvent* event) const
{
	Float_t zvtx = -999;

	const AliVVertex* primaryVertex = event->GetPrimaryVertex();

	if(primaryVertex->GetNContributors()>0)
		zvtx = primaryVertex->GetZ();

	return zvtx;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskQAHighPtDeDxTest::GetPidCode(Int_t pdgCode) const
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
void AliAnalysisTaskQAHighPtDeDxTest::ProcessMCTruthESD()
{
	// Fill the special MC histogram with the MC truth info

	const Int_t nTracksMC = fMCStack->GetNtrack();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		//Cuts
		if(!(fMCStack->IsPhysicalPrimary(iTracks)))
			continue;

		TParticle* trackMC = fMCStack->Particle(iTracks);

		TParticlePDG* pdgPart = trackMC ->GetPDG();
		Double_t chargeMC = pdgPart->Charge();

		if(chargeMC==0)
			continue;

		if (TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		Double_t etaMC = trackMC->Eta();
		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);


		for(Int_t nh = 0; nh < 9; nh++) {

			if( etaMC > etaHigh[nh]/10.0 || etaMC < etaLow[nh]/10.0 )
				continue;

			hMcIn[0][nh]->Fill(trackMC->Pt());
			hMcIn[pidCodeMC][nh]->Fill(trackMC->Pt());


		}

	}//MC track loop



}

//_____________________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::ProcessMCTruthAOD()
{
	// Fill the special MC histogram with the MC truth info


	const Int_t nTracksMC = fMCArray->GetEntriesFast();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));

		if(!trackMC){
			AliError("Cannot get MC particle");
			continue;
		}

		//Cuts
		if(!(trackMC->IsPhysicalPrimary()))
			continue;


		Double_t chargeMC = trackMC->Charge();
		if(chargeMC==0)
			continue;


		if (TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		Double_t etaMC = trackMC->Eta();
		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		//cout<<"pidcode="<<pidCodeMC<<endl;
		for(Int_t nh = 0; nh < 9; nh++) {

			if( etaMC > etaHigh[nh]/10.0 || etaMC < etaLow[nh]/10.0 )
				continue;

			hMcIn[0][nh]->Fill(trackMC->Pt());
			hMcIn[pidCodeMC][nh]->Fill(trackMC->Pt());


		}

	}//MC track loop


}


//_____________________________________________________________________________
Short_t AliAnalysisTaskQAHighPtDeDxTest::GetPythiaEventProcessType(Int_t pythiaType) {
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
Short_t AliAnalysisTaskQAHighPtDeDxTest::GetDPMjetEventProcessType(Int_t dpmJetType) {
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
ULong64_t AliAnalysisTaskQAHighPtDeDxTest::GetEventIdAsLong(AliVHeader* header) const
{
	// To have a unique id for each event in a run!
	// Modified from AliRawReader.h
	return ((ULong64_t)header->GetBunchCrossNumber()+
			(ULong64_t)header->GetOrbitNumber()*3564+
			(ULong64_t)header->GetPeriodNumber()*16777215*3564);
}


//____________________________________________________________________
TParticle* AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMother(AliStack* stack, Int_t label)
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
Int_t AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
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
AliAODMCParticle* AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherAOD(AliAODMCParticle* startParticle)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	AliAODMCParticle* mcPart = startParticle;

	while (mcPart)
	{

		if(mcPart->IsPrimary())
			return mcPart;

		Int_t mother = mcPart->GetMother();

		mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
	}

	return 0;
}


//V0______________________________________
//____________________________________________________________________
TParticle* AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherV0(AliStack* stack, Int_t label)
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
Int_t AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps)
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

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskQAHighPtDeDxTest::FindPrimaryMotherAODV0(AliAODMCParticle* startParticle, Int_t& nSteps)
{
	//
	// Finds the first mother among the primary particles of the particle identified by <label>,
	// i.e. the primary that "caused" this particle
	//
	// Taken from AliPWG0Helper class
	//

	nSteps = 0;

	AliAODMCParticle* mcPart = startParticle;

	while (mcPart)
	{

		if(mcPart->IsPrimary())
			return mcPart;

		Int_t mother = mcPart->GetMother();

		mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
		nSteps++; // 1 level down
	}

	return 0;
}



//__________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::ProduceArrayTrksESD( AliESDEvent *ESDevent ){

	const Int_t nESDTracks = ESDevent->GetNumberOfTracks();

	Int_t multTPC = 0;
	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);


		if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
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

	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
		if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		// TracCuts 2015 Pb-Pb Open DCA
		UInt_t selectDebug = 0;
		if (fTrackFilter2015PbPb) {
			selectDebug = fTrackFilter2015PbPb->IsSelected(esdTrack);
			if (!selectDebug) {
				continue;
			}
		}

		Short_t ncl = esdTrack->GetTPCsignalN();
		if(ncl<70)continue;

		Double_t eta  = esdTrack->Eta();
		Double_t phi  = esdTrack->Phi();
		Double_t momentum = esdTrack->P();
		Double_t pt = esdTrack->Pt();
		Float_t  dedx    = esdTrack->GetTPCsignal();

		Float_t dcaxy = 0.;
		Float_t dcaz = 0.;
		esdTrack->GetImpactParameters(dcaxy,dcaz);


		if(fdEdxCalibrated){
			if(eta < 0)
				dedx *= 50/EtaCalibrationNeg(cent,eta);
			else
				dedx *= 50/EtaCalibrationPos(cent,eta);
		}

		if( TMath::Abs(dcaxy) > GetMaxDCApTDep(fcutDCAxy,pt) )
			continue;

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), magf, cutLow, cutHigh))
			continue;

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

		Short_t pidCode     = 0;

		if(fAnalysisMC) {

			const Int_t label = TMath::Abs(esdTrack->GetLabel());
			TParticle* mcTrack = fMCStack->Particle(label);
			if (mcTrack){

				Int_t pdgCode = mcTrack->GetPdgCode();
				pidCode = GetPidCode(pdgCode);

			}

		}

		if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP

			if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
				if(momentum<0.6&&momentum>0.4){
					hMIPVsEta[cent]->Fill(eta,dedx);
					pMIPVsEta[cent]->Fill(eta,dedx);
				}
			}
			if( dedx > 70 && dedx < 90 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta[cent]->Fill(eta,dedx);
					pPlateauVsEta[cent]->Fill(eta,dedx);
				}
			}
		}

		hPtAll[cent]->Fill(pt);

		if(esdTrack->Charge() < 0.){

			hPtAllNeg[cent]->Fill(pt);

			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion))<2.0)
				hDCAxyVsPtPiNeg[cent]->Fill(pt,dcaxy);
			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon))<2.0)
				hDCAxyVsPtKNeg[cent]->Fill(pt,dcaxy);
			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton))<2.0)
				hDCAxyVsPtPNeg[cent]->Fill(pt,dcaxy);

			if(CloseDCAxy){
				Double_t DCAxy = GetMaxDCApTDep(fcutDCAxy,pt);
				if( TMath::Abs(dcaxy) > DCAxy )
					continue;

				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion))<2.0)
					hDCAxyVsPtPiNegC[cent]->Fill(pt,dcaxy);
				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon))<2.0)
					hDCAxyVsPtKNegC[cent]->Fill(pt,dcaxy);
				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton))<2.0)
					hDCAxyVsPtPNegC[cent]->Fill(pt,dcaxy);
			}
		}
		else{

			hPtAllPos[cent]->Fill(pt);

			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion))<3.0)
				hDCAxyVsPtPiPos[cent]->Fill(pt,dcaxy);
			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon))<3.0)
				hDCAxyVsPtKPos[cent]->Fill(pt,dcaxy);
			if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton))<3.0)
				hDCAxyVsPtPPos[cent]->Fill(pt,dcaxy);

			if(CloseDCAxy){
				Double_t DCAxy = GetMaxDCApTDep(fcutDCAxy,pt);
				if( TMath::Abs(dcaxy) > DCAxy )
					continue;

				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion))<3.0)
					hDCAxyVsPtPiPosC[cent]->Fill(pt,dcaxy);
				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon))<3.0)
					hDCAxyVsPtKPosC[cent]->Fill(pt,dcaxy);
				if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton))<3.0)
					hDCAxyVsPtPPosC[cent]->Fill(pt,dcaxy);
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


		if(fAnalysisMC){
			hMcOut[0][nh]->Fill(esdTrack->Pt());
			hMcOut[pidCode][nh]->Fill(esdTrack->Pt());
		}

		if(beta>1){
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion),2)+TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),2))<3.0){
				histPiTof[cent][nh]->Fill(momentum, dedx);
			}
		}

		if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kElectron),2)+TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kElectron),2)) < 3){
			histElTof[cent][nh]->Fill(momentum, dedx);
		}    


		if( momentum <= 0.6 && momentum >= 0.4  ){

			if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
				hMIPVsPhi[cent][nh]->Fill(phi,dedx);
				pMIPVsPhi[cent][nh]->Fill(phi,dedx);

				hMIPVsNch[cent][nh]->Fill(multTPC,dedx);
				pMIPVsNch[cent][nh]->Fill(multTPC,dedx);

			}
			if( dedx > 70 && dedx < 90 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsPhi[cent][nh]->Fill(phi,dedx);
					pPlateauVsPhi[cent][nh]->Fill(phi,dedx);

				}
			}
		}

		hPtVsP[cent][nh]->Fill(momentum,pt);
		hDeDxVsP[cent][nh]->Fill(momentum,dedx);

		if(esdTrack->Charge() < 0.){
			hPtNeg[cent][nh]->Fill(pt);

			hnSigmaPiNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaKNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaPNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
		}
		else{
			hPtPos[cent][nh]->Fill(pt);

			hnSigmaPiPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaKPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaPPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
		}



	}//end of track loop




}
//__________________________________________________________________
void AliAnalysisTaskQAHighPtDeDxTest::ProduceArrayTrksAOD( AliAODEvent *AODevent ){

	//Bool_t OfficialCalibration = kTRUE;
	Int_t nAODTracks = AODevent->GetNumberOfTracks();
	Int_t multTPC = 0;

	//get multiplicity tpc only track cuts
	for(Int_t iT = 0; iT < nAODTracks; iT++) {

		//AliAODTrack* aodTrack = AODevent->GetTrack(iT);
		AliAODTrack* aodTrack = static_cast<AliAODTrack*>(AODevent->GetTrack(iT));
		if(!aodTrack)printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");


		if(TMath::Abs(aodTrack->Eta()) > fEtaCut)
			continue;


		if (fTrackFilterTPC) {
			// TPC only cuts is bit 1
			if(!aodTrack->TestFilterBit(1))
				continue;
		}

		multTPC++;

	}


	for(Int_t iT = 0; iT < nAODTracks; iT++) {

		AliAODTrack* aodTrack = static_cast<AliAODTrack*>(AODevent->GetTrack(iT));
		AliVTrack* const trackPID = dynamic_cast<AliVTrack*>(aodTrack);
		if(!trackPID)cout << "NO trackPID    " << endl;

		if (fTrackFilterGolden) {
			// "Global track RAA analysis QM2011 + Chi2ITS<36"; bit 1024
			if(!aodTrack->TestFilterBit(1024))
				continue;
		}


		if(TMath::Abs(aodTrack->Eta()) > fEtaCut)
			continue;


		Double_t eta  = aodTrack->Eta();
		Double_t phi  = aodTrack->Phi();
		Double_t momentum = aodTrack->P();
		Double_t pt = aodTrack->Pt();


		if(!PhiCut(aodTrack->Pt(), phi, aodTrack->Charge(), magf, cutLow, cutHigh))
			continue;


		AliAODPid* aodPid = aodTrack->GetDetPid();
		Short_t ncl     = -10;
		Float_t dedx    = -10;

		//TOF
		Float_t beta = -99;


		if(aodPid) {
			ncl     = aodPid->GetTPCsignalN();
			//			if(!OfficialCalibration)
			dedx    = aodPid->GetTPCsignal();
			//			else
			// 			dedx = fPIDResponse->GetTPCResponse().GetEtaCorrectedTrackdEdx( aodTrack, AliPID::kPion,      AliTPCPIDResponse::kdEdxDefault );

			//TOF
			Bool_t IsTOFout=kFALSE;
			Float_t lengthtrack=-999;//in aod we do not have that information, beta must be: beta=inttime/timeTOF
			Float_t timeTOF=-999;

			if ((aodTrack->GetStatus()&AliESDtrack::kTOFout)==0)
				IsTOFout=kTRUE;

			lengthtrack=-999;//in aod we do not have that information, beta must be: beta=inttime/timeTOF

			timeTOF=aodPid->GetTOFsignal();

			Double_t inttime[5]={0,0,0,0,0};
			aodPid->GetIntegratedTimes(inttime);// Returns the array with integrated times for each particle hypothesis


			if ( !IsTOFout ){
				if ( ( lengthtrack != 0 ) && ( timeTOF != 0) )
					beta = inttime[0] / timeTOF;
			}

		}


		if(ncl<70)
			continue;


		Short_t pidCode     = 0;

		if(fAnalysisMC) {

			const Int_t label = TMath::Abs(aodTrack->GetLabel());
			AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label));

			if (mcTrack){

				Int_t pdgCode = mcTrack->GetPdgCode();
				pidCode = GetPidCode(pdgCode);

			}

		}

		if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP
			if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
				if(momentum<0.6&&momentum>0.4){
					hMIPVsEta[cent]->Fill(eta,dedx);
					pMIPVsEta[cent]->Fill(eta,dedx);
				}
			}
			if( dedx > DeDxMIPMax+1 && dedx < 95 ){
				if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta[cent]->Fill(eta,dedx);
					pPlateauVsEta[cent]->Fill(eta,dedx);
				}
			}
		}


		for(Int_t nh = 0; nh < 9; nh++){

			if( eta > etaHigh[nh]/10.0 || eta < etaLow[nh]/10.0 )
				continue;

			if(fAnalysisMC){
				hMcOut[0][nh]->Fill(aodTrack->Pt());
				hMcOut[pidCode][nh]->Fill(aodTrack->Pt());
			}

			if(beta>1){
				histPiTof[cent][nh]->Fill(momentum,dedx);
				histpPiTof[cent][nh]->Fill(momentum);
			}

			if( momentum <= 0.6 && momentum >= 0.4  ){//only p:0.4-0.6 GeV, pion MIP
				//Fill  pion MIP, before calibration
				if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
					hMIPVsPhi[cent][nh]->Fill(phi,dedx);
					pMIPVsPhi[cent][nh]->Fill(phi,dedx);

					hMIPVsNch[cent][nh]->Fill(multTPC,dedx);
					pMIPVsNch[cent][nh]->Fill(multTPC,dedx);

				}

				//Fill electrons, before calibration
				if( dedx > DeDxMIPMax+1 && dedx < 95 ){
					if(TMath::Abs(beta-1)<0.1){
						hPlateauVsPhi[cent][nh]->Fill(phi,dedx);
						pPlateauVsPhi[cent][nh]->Fill(phi,dedx);
					}
				}
			}

			//				hPt[nh]->Fill(pt);
			hPtVsP[cent][nh]->Fill(momentum,pt);
			//				histAllCh[nh]->Fill(momentum,dedx);
			hDeDxVsP[cent][nh]->Fill(momentum,dedx);


			if(aodTrack->Charge() < 0.){
				hPtNeg[cent][nh]->Fill(pt);
				hnSigmaPiNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion));
				hnSigmaKNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon));
				hnSigmaPNeg[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton));
			}

			else{
				hPtPos[cent][nh]->Fill(pt);
				hnSigmaPiPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion));
				hnSigmaKPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon));
				hnSigmaPPos[cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton));
			}

		}//end loop over eta intervals





	}//end of track loop




}



//----------------------------------------------------------------------------------



void AliAnalysisTaskQAHighPtDeDxTest::ProduceArrayV0ESD( AliESDEvent *ESDevent ){

	Int_t nv0s = ESDevent->GetNumberOfV0s();

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
		if (!esdV0) continue;

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

		UInt_t selectDebug_p = 0;
		if (fTrackFilterGolden) {
			selectDebug_p = fTrackFilterGolden->IsSelected(pTrack);
			if (!selectDebug_p) {
				continue;
			}
		}

		UInt_t selectDebug_n = 0;
		if (fTrackFilterGolden) {
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

		Double_t dmassG     = v0GKF.GetMass();
		Double_t dmassK     = v0K0sKF.GetMass()-0.498;
		Double_t dmassL     = v0LambdaKF.GetMass()-1.116;
		Double_t dmassAL    = v0AntiLambdaKF.GetMass()-1.116;

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

						       if(track->GetTPCsignalN()<=70)continue;
						       Double_t phi     = track->Phi();

						       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, cutLow, cutHigh))
							       continue;

						       Double_t eta      = track->Eta();
						       Double_t momentum = track->Pt();
						       Double_t dedx     = track->GetTPCsignal();

						       if(fdEdxCalibrated){
							       if(eta < 0)
								       dedx *= 50/EtaCalibrationNeg(cent,eta);
							       else
								       dedx *= 50/EtaCalibrationPos(cent,eta);
						       }

						       if(fillPos&&fillNeg){
							       if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s[cent]->Fill(eta,dedx);
									       pMIPVsEtaV0s[cent]->Fill(eta,dedx);
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

							       histPiV0[cent][nh]->Fill(momentum, dedx);
							       histpPiV0[cent][nh]->Fill(momentum);

						       }
						       else{

							       histPV0[cent][nh]->Fill(momentum, dedx);
							       histpPV0[cent][nh]->Fill(momentum);

						       }
					       }//end loop over two tracks

				       };
				       break;

				case 1:{//gammas

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;


					       if((TMath::Abs(dmassK)>0.1) && (TMath::Abs(dmassL)>0.1) && (TMath::Abs(dmassAL)>0.1)) {
						       if((dmassG<0.01) && (dmassG>0.0001)) {

							       if(nTrack->Eta() > 0)
								       if( TMath::Abs(nTrack->GetTPCsignal()-felededxfitPos->Eval(nTrack->Eta())) < 5)
									       fillPos = kTRUE;

							       if(nTrack->Eta() < 0)
								       if( TMath::Abs(nTrack->GetTPCsignal()-felededxfitNeg->Eval(nTrack->Eta())) < 5)
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

					       if(fdEdxCalibrated){
						       if(eta < 0)
							       dedx *= 50/EtaCalibrationNeg(cent,eta);
						       else
							       dedx *= 50/EtaCalibrationPos(cent,eta);
					       }


					       if(track->GetTPCsignalN()<=70)continue;

					       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, cutLow, cutHigh))
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


					       histEV0[cent][nh]->Fill(momentum, dedx);

					       cout << "Cent   "<<cent<<"DE/DX :: "<<dedx<<"Eta ::  "<<eta<< endl;


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



//__________________________________________________________________________



void AliAnalysisTaskQAHighPtDeDxTest::ProduceArrayV0AOD( AliAODEvent *AODevent ){
	Int_t nv0s = AODevent->GetNumberOfV0s();
	/*
	   if(nv0s<1){
	   return;
	   }*/

	AliAODVertex *myBestPrimaryVertex = AODevent->GetPrimaryVertex();
	if (!myBestPrimaryVertex) return;



	// ################################
	// #### BEGINNING OF V0 CODE ######
	// ################################
	// This is the begining of the V0 loop
	for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
		AliAODv0 *aodV0 = AODevent->GetV0(iV0);
		if (!aodV0) continue;


		//check onfly status
		if(aodV0->GetOnFlyStatus()!=0)
			continue;

		// AliAODTrack (V0 Daughters)
		AliAODVertex* vertex = aodV0->GetSecondaryVtx();
		if (!vertex) {
			Printf("ERROR: Could not retrieve vertex");
			continue;
		}

		AliAODTrack *pTrack = (AliAODTrack*)vertex->GetDaughter(0);
		AliAODTrack *nTrack = (AliAODTrack*)vertex->GetDaughter(1);
		if (!pTrack || !nTrack) {
			Printf("ERROR: Could not retrieve one of the daughter track");
			continue;
		}

		// Remove like-sign
		if (pTrack->Charge() == nTrack->Charge()) {
			//cout<< "like sign, continue"<< endl;
			continue;
		}

		// Make sure charge ordering is ok
		if (pTrack->Charge() < 0) {
			AliAODTrack* helpTrack = pTrack;
			pTrack = nTrack;
			nTrack = helpTrack;
		}

		// Eta cut on decay products
		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
			continue;


		Double_t dmassG  = aodV0->InvMass2Prongs(0,1,11,11);
		Double_t dmassK  = aodV0->MassK0Short()-0.498;
		Double_t dmassL  = aodV0->MassLambda()-1.116;
		Double_t dmassAL = aodV0->MassAntiLambda()-1.116;

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

						       AliAODTrack* track = 0;

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

						       if(track->GetTPCsignalN()<=70)continue;

						       Double_t phi     = track->Phi();

						       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, cutLow, cutHigh))
							       continue;


						       //if(!PhiCut(pt, phi, charge, magf, cutLow, cutHigh, hPhi))
						       //	continue;

						       Double_t eta  = track->Eta();
						       Double_t momentum = track->Pt();
						       Double_t dedx = track->GetTPCsignal();

						       if(fillPos&&fillNeg){


							       if( dedx < DeDxMIPMax && dedx > DeDxMIPMin ){
								       if(momentum<0.6&&momentum>0.4){
									       hMIPVsEtaV0s[cent]->Fill(eta,dedx);
									       pMIPVsEtaV0s[cent]->Fill(eta,dedx);
								       }
							       }


						       }

						       for(Int_t nh = 0; nh < nHists; nh++) {



							       if( eta < etaLow[nh]/10.0 || eta > etaHigh[nh]/10.0 )
								       continue;

							       if(fillPos&&fillNeg){

								       histPiV0[cent][nh]->Fill(momentum, dedx);
								       histpPiV0[cent][nh]->Fill(momentum);

							       }
							       else{

								       histPV0[cent][nh]->Fill(momentum, dedx);
								       histpPV0[cent][nh]->Fill(momentum);

							       }

						       }


					       }//end loop over two tracks
				       };
				       break;

				case 1:{//gammas

					       Bool_t fillPos = kFALSE;
					       Bool_t fillNeg = kFALSE;

					       if(dmassK>0.01 && dmassL>0.01 && dmassAL>0.01) {
						       if(dmassG<0.01 && dmassG>0.0001) {

							       if( TMath::Abs(nTrack->GetTPCsignal() - 85.0) < 5)
								       fillPos = kTRUE;
							       if( TMath::Abs(pTrack->GetTPCsignal() - 85.0) < 5)
								       fillNeg = kTRUE;

						       } else {
							       continue;
						       }
					       }


					       if(fillPos == kTRUE && fillNeg == kTRUE)
						       continue;


					       AliAODTrack* track = 0;
					       if(fillNeg)
						       track = nTrack;
					       else if(fillPos)
						       track = pTrack;
					       else
						       continue;

					       Double_t dedx  = track->GetTPCsignal();
					       Double_t eta  = track->Eta();
					       Double_t phi  = track->Phi();
					       Double_t momentum = track->P();

					       if(track->GetTPCsignalN()<=70)continue;

					       if(!PhiCut(track->Pt(), phi, track->Charge(), magf, cutLow, cutHigh))
						       continue;

					       for(Int_t nh = 0; nh < nHists; nh++) {

						       if( eta < etaLow[nh]/10.0 || eta > etaHigh[nh]/10.0 )
							       continue;

						       histEV0[cent][nh]->Fill(momentum, dedx);

					       }

				       };
				       break;


			}//end switch
		}//end loop over V0s cases

	}//end loop over v0's




}

//-------------------------------------------------------------------------
Bool_t AliAnalysisTaskQAHighPtDeDxTest::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh)

{
	if(pt < 2.0)
		return kTRUE;

	//Double_t phi = track->Phi();
	if(mag < 0)    // for negatve polarity field
		phi = TMath::TwoPi() - phi;
	if(q < 0) // for negatve charge
		phi = TMath::TwoPi()-phi;

	phi += TMath::Pi()/18.0; // to center gap in the middle
	phi = fmod(phi, TMath::Pi()/9.0);

	if(phi<phiCutHigh->Eval(pt)
			&& phi>phiCutLow->Eval(pt))
		return kFALSE; // reject track

	hPhi[cent]->Fill(pt, phi);

	return kTRUE;
}

//-------------------------------------------------------------------------

Float_t AliAnalysisTaskQAHighPtDeDxTest::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){

	Double_t maxDCAxy = 10;
	maxDCAxy = fMaxDCAxy->Eval(ptI);
	return maxDCAxy;

}


//-------------------------------------------------------------------------

Double_t AliAnalysisTaskQAHighPtDeDxTest::EtaCalibrationNeg( Int_t Cent, Double_t eta){


	for(Int_t i=0; i<4; ++i){
		fEtaCalibrationNeg->SetParameter(i,0);
		felededxfitNeg->SetParameter(i,0);
	}

	switch(Cent){

		case 0:
			fEtaCalibrationNeg->SetParameter(0,aNeg[0]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[0]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[0]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[0]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[0]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[0]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[0]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[0]);

			felededxfitNeg->SetParameter(0,aNegEl[0]);
			felededxfitNeg->FixParameter(0,aNegEl[0]);
			felededxfitNeg->SetParameter(1,bNegEl[0]);
			felededxfitNeg->FixParameter(1,bNegEl[0]);
			felededxfitNeg->SetParameter(2,cNegEl[0]);
			felededxfitNeg->FixParameter(2,cNegEl[0]);
			felededxfitNeg->SetParameter(3,dNegEl[0]);
			felededxfitNeg->FixParameter(3,dNegEl[0]);
			break;

		case 1:
			fEtaCalibrationNeg->SetParameter(0,aNeg[1]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[1]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[1]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[1]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[1]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[1]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[1]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[1]);

			felededxfitNeg->SetParameter(0,aNegEl[1]);
			felededxfitNeg->FixParameter(0,aNegEl[1]);
			felededxfitNeg->SetParameter(1,bNegEl[1]);
			felededxfitNeg->FixParameter(1,bNegEl[1]);
			felededxfitNeg->SetParameter(2,cNegEl[1]);
			felededxfitNeg->FixParameter(2,cNegEl[1]);
			felededxfitNeg->SetParameter(3,dNegEl[1]);
			felededxfitNeg->FixParameter(3,dNegEl[1]);
			break;

		case 2:
			fEtaCalibrationNeg->SetParameter(0,aNeg[2]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[2]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[2]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[2]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[2]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[2]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[2]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[2]);

			felededxfitNeg->SetParameter(0,aNegEl[2]);
			felededxfitNeg->FixParameter(0,aNegEl[2]);
			felededxfitNeg->SetParameter(1,bNegEl[2]);
			felededxfitNeg->FixParameter(1,bNegEl[2]);
			felededxfitNeg->SetParameter(2,cNegEl[2]);
			felededxfitNeg->FixParameter(2,cNegEl[2]);
			felededxfitNeg->SetParameter(3,dNegEl[2]);
			felededxfitNeg->FixParameter(3,dNegEl[2]);
			break;

		case 3:
			fEtaCalibrationNeg->SetParameter(0,aNeg[3]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[3]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[3]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[3]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[3]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[3]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[3]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[3]);

			felededxfitNeg->SetParameter(0,aNegEl[3]);
			felededxfitNeg->FixParameter(0,aNegEl[3]);
			felededxfitNeg->SetParameter(1,bNegEl[3]);
			felededxfitNeg->FixParameter(1,bNegEl[3]);
			felededxfitNeg->SetParameter(2,cNegEl[3]);
			felededxfitNeg->FixParameter(2,cNegEl[3]);
			felededxfitNeg->SetParameter(3,dNegEl[3]);
			felededxfitNeg->FixParameter(3,dNegEl[3]);
			break;

		case 4:
			fEtaCalibrationNeg->SetParameter(0,aNeg[4]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[4]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[4]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[4]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[4]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[4]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[4]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[4]);

			felededxfitNeg->SetParameter(0,aNegEl[4]);
			felededxfitNeg->FixParameter(0,aNegEl[4]);
			felededxfitNeg->SetParameter(1,bNegEl[4]);
			felededxfitNeg->FixParameter(1,bNegEl[4]);
			felededxfitNeg->SetParameter(2,cNegEl[4]);
			felededxfitNeg->FixParameter(2,cNegEl[4]);
			felededxfitNeg->SetParameter(3,dNegEl[4]);
			felededxfitNeg->FixParameter(3,dNegEl[4]);
			break;

		case 5:
			fEtaCalibrationNeg->SetParameter(0,aNeg[5]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[5]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[5]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[5]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[5]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[5]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[5]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[5]);

			felededxfitNeg->SetParameter(0,aNegEl[5]);
			felededxfitNeg->FixParameter(0,aNegEl[5]);
			felededxfitNeg->SetParameter(1,bNegEl[5]);
			felededxfitNeg->FixParameter(1,bNegEl[5]);
			felededxfitNeg->SetParameter(2,cNegEl[5]);
			felededxfitNeg->FixParameter(2,cNegEl[5]);
			felededxfitNeg->SetParameter(3,dNegEl[5]);
			felededxfitNeg->FixParameter(3,dNegEl[5]);
			break;

		case 6:
			fEtaCalibrationNeg->SetParameter(0,aNeg[6]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[6]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[6]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[6]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[6]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[6]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[6]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[6]);

			felededxfitNeg->SetParameter(0,aNegEl[6]);
			felededxfitNeg->FixParameter(0,aNegEl[6]);
			felededxfitNeg->SetParameter(1,bNegEl[6]);
			felededxfitNeg->FixParameter(1,bNegEl[6]);
			felededxfitNeg->SetParameter(2,cNegEl[6]);
			felededxfitNeg->FixParameter(2,cNegEl[6]);
			felededxfitNeg->SetParameter(3,dNegEl[6]);
			felededxfitNeg->FixParameter(3,dNegEl[6]);
			break;

		case 7:
			fEtaCalibrationNeg->SetParameter(0,aNeg[7]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[7]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[7]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[7]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[7]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[7]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[7]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[7]);

			felededxfitNeg->SetParameter(0,aNegEl[7]);
			felededxfitNeg->FixParameter(0,aNegEl[7]);
			felededxfitNeg->SetParameter(1,bNegEl[7]);
			felededxfitNeg->FixParameter(1,bNegEl[7]);
			felededxfitNeg->SetParameter(2,cNegEl[7]);
			felededxfitNeg->FixParameter(2,cNegEl[7]);
			felededxfitNeg->SetParameter(3,dNegEl[7]);
			felededxfitNeg->FixParameter(3,dNegEl[7]);
			break;

		case 8:
			fEtaCalibrationNeg->SetParameter(0,aNeg[8]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[8]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[8]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[8]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[8]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[8]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[8]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[8]);

			felededxfitNeg->SetParameter(0,aNegEl[8]);
			felededxfitNeg->FixParameter(0,aNegEl[8]);
			felededxfitNeg->SetParameter(1,bNegEl[8]);
			felededxfitNeg->FixParameter(1,bNegEl[8]);
			felededxfitNeg->SetParameter(2,cNegEl[8]);
			felededxfitNeg->FixParameter(2,cNegEl[8]);
			felededxfitNeg->SetParameter(3,dNegEl[8]);
			felededxfitNeg->FixParameter(3,dNegEl[8]);
			break;

		case 9:
			fEtaCalibrationNeg->SetParameter(0,aNeg[9]);
			fEtaCalibrationNeg->FixParameter(0,aNeg[9]);
			fEtaCalibrationNeg->SetParameter(1,bNeg[9]);
			fEtaCalibrationNeg->FixParameter(1,bNeg[9]);
			fEtaCalibrationNeg->SetParameter(2,cNeg[9]);
			fEtaCalibrationNeg->FixParameter(2,cNeg[9]);
			fEtaCalibrationNeg->SetParameter(3,dNeg[9]);
			fEtaCalibrationNeg->FixParameter(3,dNeg[9]);

			felededxfitNeg->SetParameter(0,aNegEl[9]);
			felededxfitNeg->FixParameter(0,aNegEl[9]);
			felededxfitNeg->SetParameter(1,bNegEl[9]);
			felededxfitNeg->FixParameter(1,bNegEl[9]);
			felededxfitNeg->SetParameter(2,cNegEl[9]);
			felededxfitNeg->FixParameter(2,cNegEl[9]);
			felededxfitNeg->SetParameter(3,dNegEl[9]);
			felededxfitNeg->FixParameter(3,dNegEl[9]);
			break;

		default:
			cout<<"Wrong Centrality Parameter"<<endl;

	}

	return fEtaCalibrationNeg->Eval(eta);

}


//-------------------------------------------------------------------------

Double_t AliAnalysisTaskQAHighPtDeDxTest::EtaCalibrationPos( Int_t Cent, Double_t eta){


	for(Int_t i=0; i<4; ++i)
		fEtaCalibration->SetParameter(i,0);

	switch(Cent){

		case 0:
			fEtaCalibration->SetParameter(0,aPos[0]);
			fEtaCalibration->FixParameter(0,aPos[0]);
			fEtaCalibration->SetParameter(1,bPos[0]);
			fEtaCalibration->FixParameter(1,bPos[0]);
			fEtaCalibration->SetParameter(2,cPos[0]);
			fEtaCalibration->FixParameter(2,cPos[0]);
			fEtaCalibration->SetParameter(3,dPos[0]);
			fEtaCalibration->FixParameter(3,dPos[0]);

			felededxfitPos->SetParameter(0,aPosEl[0]);
			felededxfitPos->FixParameter(0,aPosEl[0]);
			felededxfitPos->SetParameter(1,bPosEl[0]);
			felededxfitPos->FixParameter(1,bPosEl[0]);
			felededxfitPos->SetParameter(2,cPosEl[0]);
			felededxfitPos->FixParameter(2,cPosEl[0]);
			felededxfitPos->SetParameter(3,dPosEl[0]);
			felededxfitPos->FixParameter(3,dPosEl[0]);
			break;

		case 1:
			fEtaCalibration->SetParameter(0,aPos[1]);
			fEtaCalibration->FixParameter(0,aPos[1]);
			fEtaCalibration->SetParameter(1,bPos[1]);
			fEtaCalibration->FixParameter(1,bPos[1]);
			fEtaCalibration->SetParameter(2,cPos[1]);
			fEtaCalibration->FixParameter(2,cPos[1]);
			fEtaCalibration->SetParameter(3,dPos[1]);
			fEtaCalibration->FixParameter(3,dPos[1]);

			felededxfitPos->SetParameter(0,aPosEl[1]);
			felededxfitPos->FixParameter(0,aPosEl[1]);
			felededxfitPos->SetParameter(1,bPosEl[1]);
			felededxfitPos->FixParameter(1,bPosEl[1]);
			felededxfitPos->SetParameter(2,cPosEl[1]);
			felededxfitPos->FixParameter(2,cPosEl[1]);
			felededxfitPos->SetParameter(3,dPosEl[1]);
			felededxfitPos->FixParameter(3,dPosEl[1]);
			break;

		case 2:
			fEtaCalibration->SetParameter(0,aPos[2]);
			fEtaCalibration->FixParameter(0,aPos[2]);
			fEtaCalibration->SetParameter(1,bPos[2]);
			fEtaCalibration->FixParameter(1,bPos[2]);
			fEtaCalibration->SetParameter(2,cPos[2]);
			fEtaCalibration->FixParameter(2,cPos[2]);
			fEtaCalibration->SetParameter(3,dPos[2]);
			fEtaCalibration->FixParameter(3,dPos[2]);

			felededxfitPos->SetParameter(0,aPosEl[2]);
			felededxfitPos->FixParameter(0,aPosEl[2]);
			felededxfitPos->SetParameter(1,bPosEl[2]);
			felededxfitPos->FixParameter(1,bPosEl[2]);
			felededxfitPos->SetParameter(2,cPosEl[2]);
			felededxfitPos->FixParameter(2,cPosEl[2]);
			felededxfitPos->SetParameter(3,dPosEl[2]);
			felededxfitPos->FixParameter(3,dPosEl[2]);
			break;

		case 3:
			fEtaCalibration->SetParameter(0,aPos[3]);
			fEtaCalibration->FixParameter(0,aPos[3]);
			fEtaCalibration->SetParameter(1,bPos[3]);
			fEtaCalibration->FixParameter(1,bPos[3]);
			fEtaCalibration->SetParameter(2,cPos[3]);
			fEtaCalibration->FixParameter(2,cPos[3]);
			fEtaCalibration->SetParameter(3,dPos[3]);
			fEtaCalibration->FixParameter(3,dPos[3]);

			felededxfitPos->SetParameter(0,aPosEl[3]);
			felededxfitPos->FixParameter(0,aPosEl[3]);
			felededxfitPos->SetParameter(1,bPosEl[3]);
			felededxfitPos->FixParameter(1,bPosEl[3]);
			felededxfitPos->SetParameter(2,cPosEl[3]);
			felededxfitPos->FixParameter(2,cPosEl[3]);
			felededxfitPos->SetParameter(3,dPosEl[3]);
			felededxfitPos->FixParameter(3,dPosEl[3]);
			break;

		case 4:
			fEtaCalibration->SetParameter(0,aPos[4]);
			fEtaCalibration->FixParameter(0,aPos[4]);
			fEtaCalibration->SetParameter(1,bPos[4]);
			fEtaCalibration->FixParameter(1,bPos[4]);
			fEtaCalibration->SetParameter(2,cPos[4]);
			fEtaCalibration->FixParameter(2,cPos[4]);
			fEtaCalibration->SetParameter(3,dPos[4]);
			fEtaCalibration->FixParameter(3,dPos[4]);

			felededxfitPos->SetParameter(0,aPosEl[4]);
			felededxfitPos->FixParameter(0,aPosEl[4]);
			felededxfitPos->SetParameter(1,bPosEl[4]);
			felededxfitPos->FixParameter(1,bPosEl[4]);
			felededxfitPos->SetParameter(2,cPosEl[4]);
			felededxfitPos->FixParameter(2,cPosEl[4]);
			felededxfitPos->SetParameter(3,dPosEl[4]);
			felededxfitPos->FixParameter(3,dPosEl[4]);
			break;

		case 5:
			fEtaCalibration->SetParameter(0,aPos[5]);
			fEtaCalibration->FixParameter(0,aPos[5]);
			fEtaCalibration->SetParameter(1,bPos[5]);
			fEtaCalibration->FixParameter(1,bPos[5]);
			fEtaCalibration->SetParameter(2,cPos[5]);
			fEtaCalibration->FixParameter(2,cPos[5]);
			fEtaCalibration->SetParameter(3,dPos[5]);
			fEtaCalibration->FixParameter(3,dPos[5]);

			felededxfitPos->SetParameter(0,aPosEl[5]);
			felededxfitPos->FixParameter(0,aPosEl[5]);
			felededxfitPos->SetParameter(1,bPosEl[5]);
			felededxfitPos->FixParameter(1,bPosEl[5]);
			felededxfitPos->SetParameter(2,cPosEl[5]);
			felededxfitPos->FixParameter(2,cPosEl[5]);
			felededxfitPos->SetParameter(3,dPosEl[5]);
			felededxfitPos->FixParameter(3,dPosEl[5]);
			break;

		case 6:
			fEtaCalibration->SetParameter(0,aPos[6]);
			fEtaCalibration->FixParameter(0,aPos[6]);
			fEtaCalibration->SetParameter(1,bPos[6]);
			fEtaCalibration->FixParameter(1,bPos[6]);
			fEtaCalibration->SetParameter(2,cPos[6]);
			fEtaCalibration->FixParameter(2,cPos[6]);
			fEtaCalibration->SetParameter(3,dPos[6]);
			fEtaCalibration->FixParameter(3,dPos[6]);

			felededxfitPos->SetParameter(0,aPosEl[6]);
			felededxfitPos->FixParameter(0,aPosEl[6]);
			felededxfitPos->SetParameter(1,bPosEl[6]);
			felededxfitPos->FixParameter(1,bPosEl[6]);
			felededxfitPos->SetParameter(2,cPosEl[6]);
			felededxfitPos->FixParameter(2,cPosEl[6]);
			felededxfitPos->SetParameter(3,dPosEl[6]);
			felededxfitPos->FixParameter(3,dPosEl[6]);
			break;

		case 7:
			fEtaCalibration->SetParameter(0,aPos[7]);
			fEtaCalibration->FixParameter(0,aPos[7]);
			fEtaCalibration->SetParameter(1,bPos[7]);
			fEtaCalibration->FixParameter(1,bPos[7]);
			fEtaCalibration->SetParameter(2,cPos[7]);
			fEtaCalibration->FixParameter(2,cPos[7]);
			fEtaCalibration->SetParameter(3,dPos[7]);
			fEtaCalibration->FixParameter(3,dPos[7]);

			felededxfitPos->SetParameter(0,aPosEl[7]);
			felededxfitPos->FixParameter(0,aPosEl[7]);
			felededxfitPos->SetParameter(1,bPosEl[7]);
			felededxfitPos->FixParameter(1,bPosEl[7]);
			felededxfitPos->SetParameter(2,cPosEl[7]);
			felededxfitPos->FixParameter(2,cPosEl[7]);
			felededxfitPos->SetParameter(3,dPosEl[7]);
			felededxfitPos->FixParameter(3,dPosEl[7]);
			break;

		case 8:
			fEtaCalibration->SetParameter(0,aPos[8]);
			fEtaCalibration->FixParameter(0,aPos[8]);
			fEtaCalibration->SetParameter(1,bPos[8]);
			fEtaCalibration->FixParameter(1,bPos[8]);
			fEtaCalibration->SetParameter(2,cPos[8]);
			fEtaCalibration->FixParameter(2,cPos[8]);
			fEtaCalibration->SetParameter(3,dPos[8]);
			fEtaCalibration->FixParameter(3,dPos[8]);

			felededxfitPos->SetParameter(0,aPosEl[8]);
			felededxfitPos->FixParameter(0,aPosEl[8]);
			felededxfitPos->SetParameter(1,bPosEl[8]);
			felededxfitPos->FixParameter(1,bPosEl[8]);
			felededxfitPos->SetParameter(2,cPosEl[8]);
			felededxfitPos->FixParameter(2,cPosEl[8]);
			felededxfitPos->SetParameter(3,dPosEl[8]);
			felededxfitPos->FixParameter(3,dPosEl[8]);
			break;

		case 9:
			fEtaCalibration->SetParameter(0,aPos[9]);
			fEtaCalibration->FixParameter(0,aPos[9]);
			fEtaCalibration->SetParameter(1,bPos[9]);
			fEtaCalibration->FixParameter(1,bPos[9]);
			fEtaCalibration->SetParameter(2,cPos[9]);
			fEtaCalibration->FixParameter(2,cPos[9]);
			fEtaCalibration->SetParameter(3,dPos[9]);
			fEtaCalibration->FixParameter(3,dPos[9]);

			felededxfitPos->SetParameter(0,aPosEl[9]);
			felededxfitPos->FixParameter(0,aPosEl[9]);
			felededxfitPos->SetParameter(1,bPosEl[9]);
			felededxfitPos->FixParameter(1,bPosEl[9]);
			felededxfitPos->SetParameter(2,cPosEl[9]);
			felededxfitPos->FixParameter(2,cPosEl[9]);
			felededxfitPos->SetParameter(3,dPosEl[9]);
			felededxfitPos->FixParameter(3,dPosEl[9]);
			break;

		default:
			cout<<"Wrong Centrality Parameter"<<endl;

	}

	return fEtaCalibration->Eval(eta);

}
