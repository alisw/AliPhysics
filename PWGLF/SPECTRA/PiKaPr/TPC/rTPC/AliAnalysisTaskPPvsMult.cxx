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

const Double_t aPos[nCent]      = {49.9799,49.9659,0,0,0,0,0,0,0,0,0};
const Double_t bPos[nCent]      = {2.99619,2.91366,0,0,0,0,0,0,0,0,0};
const Double_t cPos[nCent]      = {-45.718,-45.5994,0,0,0,0,0,0,0,0,0};
const Double_t dPos[nCent]      = {290.013,290.042,0,0,0,0,0,0,0,0,0};
const Double_t ePos[nCent]      = {-1018.42,-1014.49,0,0,0,0,0,0,0,0,0};
const Double_t fPos[nCent]      = {1948.68,1931.84,0,0,0,0,0,0,0,0,0};
const Double_t gPos[nCent]      = {-1864.06,-1839.36,0,0,0,0,0,0,0,0,0};
const Double_t hPos[nCent]      = {692.752,680.421,0,0,0,0,0,0,0,0,0};

const Double_t aNeg[nCent]      = {50.078,50.046,0,0,0,0,0,0,0,0,0};
const Double_t bNeg[nCent]      = {6.67199,6.79992,0,0,0,0,0,0,0,0,0};
const Double_t cNeg[nCent]      = {103.662,109.86,0,0,0,0,0,0,0,0,0};
const Double_t dNeg[nCent]      = {611.034,668.241,0,0,0,0,0,0,0,0,0};
const Double_t eNeg[nCent]      = {1695.63,1916.44,0,0,0,0,0,0,0,0,0};
const Double_t fNeg[nCent]      = {2395.88,2815.04,0,0,0,0,0,0,0,0,0};
const Double_t gNeg[nCent]      = {1669.22,2057.21,0,0,0,0,0,0,0,0,0};
const Double_t hNeg[nCent]      = {455.362,595.391,0,0,0,0,0,0,0,0,0};

const Double_t aPosEl[nCent]    = {80.1263,79.9957,0,0,0,0,0,0,0,0,0};
const Double_t bPosEl[nCent]    = {5.28525,7.03079,0,0,0,0,0,0,0,0,0};
const Double_t cPosEl[nCent]    = {-32.7731,-42.9098,0,0,0,0,0,0,0,0,0};
const Double_t dPosEl[nCent]    = {68.4524,88.7057,0,0,0,0,0,0,0,0,0};
const Double_t ePosEl[nCent]    = {-44.1566,-56.6554,0,0,0,0,0,0,0,0,0};

const Double_t aNegEl[nCent]    = {79.8351,79.7387,0,0,0,0,0,0,0,0,0};
const Double_t bNegEl[nCent]    = {-8.46921,-8.60021,0,0,0,0,0,0,0,0,0};
const Double_t cNegEl[nCent]    = {-44.5947,-44.1718,0,0,0,0,0,0,0,0,0};
const Double_t dNegEl[nCent]    = {-86.2242,-84.4984,0,0,0,0,0,0,0,0,0};
const Double_t eNegEl[nCent]    = {-53.6285,-51.945,0,0,0,0,0,0,0,0,0};

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
	fTrackFilter2015PbPb(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilterTPC(0x0),
	utils(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fAnalysisPbPb(kFALSE),
	fRandom(0x0),
	fVtxCut(0x0),
        fNcl(70),
	fEtaCut(0.9),
	cent(3),
	fMinCent(0.0),
	fMaxCent(100.0),
        fDeDxMIPMin(40),
        fDeDxMIPMax(60),
        fdEdxHigh(200),
        fdEdxLow(40),
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
	fMakePid(0x0),
	fLowPt(0x0),
	fLHC16l(0x0),
	fcent(0x0),
	fcentAfterPrimaries(0x0),
	fcentAfterV0s(0x0),
	fEtaCalibrationNeg(0x0),
	fEtaCalibration(0x0),
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
		histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
		histEV0[i][j]=0;

		}
	}

	//default constructor
	for(Int_t cent=0;cent<nCent;++cent){
		for(Int_t pid=0;pid<7;++pid){
			hMcIn[cent][pid]     = 0;
			hMcOut[cent][pid]    = 0;
			hMcInNeg[cent][pid]  = 0;
			hMcInPos[cent][pid]  = 0;
			hMcOutNeg[cent][pid] = 0;
			hMcOutPos[cent][pid] = 0;
		}
  			hPiondEdx[cent]=0;
                        hKaondEdx[cent]=0;
                        hProtondEdx[cent]=0;

	}

	for(Int_t cent=0;cent<10;++cent){
		for(Int_t pid=0;pid<7;++pid){
			for(Int_t q=0;q<3;++q){
			hDCApTPrim[cent][pid][q]  = 0;
			hDCApTWDec[cent][pid][q]  = 0;
			hDCApTMate[cent][pid][q]  = 0;

			hDCApTPrim2[cent][pid][q] = 0;
			hDCApTWDec2[cent][pid][q] = 0;
			hDCApTMate2[cent][pid][q] = 0;
			}
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
	fTrackFilter2015PbPb(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilterTPC(0x0),
	utils(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fAnalysisPbPb(kFALSE),
	fRandom(0x0),
	fVtxCut(0x0),
        fNcl(70),
	fEtaCut(0.9),
	cent(3),
	fMinCent(0.0),
	fMaxCent(100.0),
        fDeDxMIPMin(40),
        fDeDxMIPMax(60),
        fdEdxHigh(200),
        fdEdxLow(40),
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
	fMakePid(0x0),
	fLowPt(0x0),
	fLHC16l(0x0),
	fcent(0x0),
	fcentAfterPrimaries(0x0),
	fcentAfterV0s(0x0),
	fEtaCalibrationNeg(0x0),
	fEtaCalibration(0x0),
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
			histpPiTof[i][j]=0;//TH1D, for a "clean" sample of pions, beta>1
			histEV0[i][j]=0;

		}

	}

	// Default constructor (should not be used)
	for(Int_t cent=0; cent<nCent; ++cent){
		for(Int_t pid=0; pid<7; ++pid){
			hMcIn[cent][pid]=0;
			hMcOut[cent][pid]=0;
			hMcInNeg[cent][pid]=0;
			hMcInPos[cent][pid]=0;
			hMcOutNeg[cent][pid]=0;
			hMcOutPos[cent][pid]=0;
		}
		
			hPiondEdx[cent]=0;
                        hKaondEdx[cent]=0;
                        hProtondEdx[cent]=0;
	}

	for(Int_t cent=0;cent<10;++cent){
		for(Int_t pid=0;pid<7;++pid){
			for(Int_t q=0;q<3;++q){
				hDCApTPrim[cent][pid][q]=0;
				hDCApTWDec[cent][pid][q]=0;
				hDCApTMate[cent][pid][q]=0;

				hDCApTPrim2[cent][pid][q]=0;
				hDCApTWDec2[cent][pid][q]=0;
				hDCApTMate2[cent][pid][q]=0;
			}
		}
	}


	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());//esto es nuevo
}




AliAnalysisTaskPPvsMult::~AliAnalysisTaskPPvsMult() {
	//
	// Destructor
	//

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


	//OpenFile(1);
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
	fV0M->Sumw2();
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

        const Char_t *Pid[7]       = {"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};
        const Char_t *Q[3]         = {"", "Neg", "Pos"};

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

        fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+TMath::Pi()/18.0-0.025", 0, 50);
        fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+TMath::Pi()/18.0+0.035", 0, 50);


	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);
	fEtaCalibration    = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);

	felededxfitPos     = new TF1("felededxfitPos", "pol4", 0.0, 1.0);
	felededxfitNeg     = new TF1("felededxfitNeg", "pol4", -1.0, 0.0);


	Int_t nPhiBins = 36;

	for(Int_t i = 0; i<nCent; ++i){

		hMIPVsEta[i] = new TH2D(Form("hMIPVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, primary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsEta[i] = new TProfile(Form("pMIPVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, primary tracks}",50,-0.8,0.8, fDeDxMIPMin, fDeDxMIPMax);
		hMIPVsEtaV0s[i] = new TH2D(Form("hMIPVsEtaV0s%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{MIP, secondary tracks}",50,-0.8,0.8,fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsEtaV0s[i] = new TProfile(Form("pMIPVsEtaV0s%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{MIP, secondary tracks}",50,-0.8,0.8,fDeDxMIPMin, fDeDxMIPMax);
		hPlateauVsEta[i] = new TH2D(Form("hPlateauVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; dE/dx_{Plateau, primary tracks}",50,-0.8,0.8,50, 60, 110);
		pPlateauVsEta[i] = new TProfile(Form("pPlateauVsEta%.2f-%.2f",CentMin[i],CentMax[i]),"; #eta; #LT dE/dx #GT_{Plateau, primary tracks}",50,-0.8,0.8, 60, 110);
		hPhi[i] = new TH2D(Form("histPhi%.2f-%.2f",CentMin[i],CentMax[i]), ";pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4);

		hPtAll[i] = new TH1D(Form("hPt%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAll[i]->Sumw2();
		hPtAllNeg[i] = new TH1D(Form("hPtNeg%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAllNeg[i]->Sumw2();
		hPtAllPos[i] = new TH1D(Form("hPtPos%.2f-%.2f",CentMin[i],CentMax[i]),";#it{p}_{T};Counts",nPtBins,ptBins);
		hPtAllPos[i]->Sumw2();

		hDCAxyVsPtPiNeg[i] = new TH2D(Form("hDCAxyVsPtPiNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPiNegC[i] = new TH2D(Form("hDCAxyVsPtPiNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtPiNeg[i]->Sumw2();
		hDCAxyVsPtPiNegC[i]->Sumw2();
		hDCAxyVsPtKNeg[i] = new TH2D(Form("hDCAxyVsPtKNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtKNegC[i] = new TH2D(Form("hDCAxyVsPtKNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtKNeg[i]->Sumw2();
		hDCAxyVsPtKNegC[i]->Sumw2();
		hDCAxyVsPtPNeg[i] = new TH2D(Form("hDCAxyVsPtPNeg%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPNeg; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPNegC[i] = new TH2D(Form("hDCAxyVsPtPNegC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPNeg Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtPNeg[i]->Sumw2();
		hDCAxyVsPtPNegC[i]->Sumw2();
		hDCAxyVsPtPiPos[i] = new TH2D(Form("hDCAxyVsPtPiPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPiPosC[i] = new TH2D(Form("hDCAxyVsPtPiPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPiPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtPiPos[i]->Sumw2();
		hDCAxyVsPtPiPosC[i]->Sumw2();
		hDCAxyVsPtKPos[i] = new TH2D(Form("hDCAxyVsPtKPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtKPosC[i] = new TH2D(Form("hDCAxyVsPtKPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtKPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtKPos[i]->Sumw2();
		hDCAxyVsPtKPosC[i]->Sumw2();
		hDCAxyVsPtPPos[i] = new TH2D(Form("hDCAxyVsPtPPos%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPPos; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
		hDCAxyVsPtPPosC[i] = new TH2D(Form("hDCAxyVsPtPPosC%.2f-%.2f",CentMin[i],CentMax[i]), "hDCAxyVsPtPPos Close; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
		hDCAxyVsPtPPos[i]->Sumw2();
		hDCAxyVsPtPPosC[i]->Sumw2();

		for(Int_t j=0; j<nHists; j++) {

		hDeDxVsP[i][j] = new TH2D( Form("hDeDxVsP%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), ";#it{p} [GeV/c]; dE/dx", nPtBins, ptBins, fdEdxHigh-fdEdxLow, fdEdxLow, fdEdxHigh);
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

		hMIPVsPhi[i][j] = new TH2D(Form("hMIPVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		hMIPVsPhi[i][j]->Sumw2();

		pMIPVsPhi[i][j] = new TProfile(Form("pMIPVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx MIP",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsPhi[i][j]->Sumw2();

		hPlateauVsPhi[i][j]  = new TH2D(Form("hPlateauVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]),  nPhiBins, 0, 2*TMath::Pi(),20, 70, 90);
		hPlateauVsPhi[i][j]->Sumw2();

		pPlateauVsPhi[i][j] = new TProfile(Form("pPlateauVsPhi%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), Form("%s; #phi (rad); dE/dx Plateau",LatexEta[j]), nPhiBins, 0, 2*TMath::Pi(),fDeDxMIPMax, 95);
		pPlateauVsPhi[i][j]->Sumw2();

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

		histpPiTof[i][j]  = new TH1D(Form("hPTOF%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Primary Pions from TOF ; #it{p} (GeV/#it{c}); counts", 200, 0, 20);
		histpPiTof[i][j]->Sumw2();

		histEV0[i][j]   = new TH2D(Form("hEV0%.2f-%.2f-%s",CentMin[i],CentMax[i],ending[j]), "Electrons id by V0", nPtBinsV0s, ptBinsV0s, nDeltaPiBins, deltaPiLow, deltaPiHigh);
		histEV0[i][j]->Sumw2();

		}// eta loop
	} // centrality loop


	for(Int_t i = 0; i<nHists; ++i ){

		hMIPVsV0M[i] = new TH2D(Form("hMIPVsV0M-%s",ending[i]), Form("%s; V0M mult.; dE/dx MIP",LatexEta[i]), 100, 0, 100, fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		hMIPVsV0M[i]->Sumw2();
		pMIPVsV0M[i] = new TProfile(Form("pMIPVsV0M-%s",ending[i]), Form("%s; V0M mult.; dE/dx MIP",LatexEta[i]), 100, 0, 100, fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsV0M[i]->Sumw2();
		hMIPVsNch[i] = new TH2D(Form("hMIPVsNch-%s",ending[i]),"; TPC track mult. |#eta|<0.8; dE/dx MIP", 100, 1, 101, fDeDxMIPMax-fDeDxMIPMin, fDeDxMIPMin, fDeDxMIPMax);
		hMIPVsNch[i]->Sumw2();
		pMIPVsNch[i] = new TProfile(Form("pMIPVsNch-%s",ending[i]),"; TPC track mult. |#eta|<0.8; dE/dx MIP", 100, 1, 101, fDeDxMIPMin, fDeDxMIPMax);
		pMIPVsNch[i]->Sumw2();
}



		if(!fAnalysisMC){
		for(Int_t j=0; j<nHists; ++j){
		fListOfObjects->Add(hMIPVsV0M[j]);
		fListOfObjects->Add(pMIPVsV0M[j]);
		fListOfObjects->Add(hMIPVsNch[j]);
		fListOfObjects->Add(pMIPVsNch[j]);
		}

		for(Int_t i=0; i<nCent; ++i ){
		fListOfObjects->Add(hMIPVsEta[i]);
		fListOfObjects->Add(pMIPVsEta[i]);
		fListOfObjects->Add(hMIPVsEtaV0s[i]);
		fListOfObjects->Add(pMIPVsEtaV0s[i]);
		fListOfObjects->Add(hPlateauVsEta[i]);
		fListOfObjects->Add(pPlateauVsEta[i]);
		fListOfObjects->Add(hPhi[i]);

		for(Int_t j=0; j<nHists; ++j){
		fListOfObjects->Add(hMIPVsPhi[i][j]);
		fListOfObjects->Add(pMIPVsPhi[i][j]);
		fListOfObjects->Add(hPlateauVsPhi[i][j]);
		fListOfObjects->Add(pPlateauVsPhi[i][j]);
			}

		if(fMakePid){
		fListOfObjects->Add(hPtAll[i]);	
		if(fLowPt){
		fListOfObjects->Add(hPtAllNeg[i]);
		fListOfObjects->Add(hPtAllPos[i]);

		fListOfObjects->Add(hDCAxyVsPtPiNeg[i]);
		fListOfObjects->Add(hDCAxyVsPtPiPos[i]);
		fListOfObjects->Add(hDCAxyVsPtPNeg[i]);
		fListOfObjects->Add(hDCAxyVsPtPPos[i]);

		fListOfObjects->Add(hDCAxyVsPtPiNegC[i]);
		fListOfObjects->Add(hDCAxyVsPtPiPosC[i]);
		fListOfObjects->Add(hDCAxyVsPtPNegC[i]);
		fListOfObjects->Add(hDCAxyVsPtPPosC[i]);
		}

		for(Int_t j=0; j<nHists; ++j){
		if(fLowPt){
		fListOfObjects->Add(hnSigmaPiPos[i][j]);
		fListOfObjects->Add(hnSigmaPiNeg[i][j]);	
		fListOfObjects->Add(hnSigmaKPos[i][j]);
		fListOfObjects->Add(hnSigmaKNeg[i][j]);
		fListOfObjects->Add(hnSigmaPPos[i][j]);
		fListOfObjects->Add(hnSigmaPNeg[i][j]);
		fListOfObjects->Add(hPtPos[i][j]);
		fListOfObjects->Add(hPtNeg[i][j]);
		}

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
			} //	if(MakePID) 
		} //	Cent
	} //	!fAnalysisMC


	else{
	for(Int_t cent=0; cent<nCent; cent++) {
	for(Int_t pid=0; pid<7; pid++) {
	hMcIn[cent][pid]=new TH1D(Form("hIn_%.2f-%.2f-%s",CentMin[cent],CentMax[cent],Pid[pid]), Form("MC in (pid %s)", Pid[pid]),nPtBins,ptBins);
	hMcIn[cent][pid]->Sumw2();
	hMcInNeg[cent][pid]=new TH1D(Form("hInNeg_%.2f-%.2f-%s",CentMin[cent],CentMax[cent],Pid[pid]),Form("MC in (pid %s)",Pid[pid]),nPtBins,ptBins);
	hMcInNeg[cent][pid]->Sumw2();
	hMcInPos[cent][pid]=new TH1D(Form("hInPos_%.2f-%.2f-%s",CentMin[cent],CentMax[cent],Pid[pid]),Form("MC in (pid %s)",Pid[pid]),nPtBins,ptBins);
	hMcInPos[cent][pid]->Sumw2();
	hMcOut[cent][pid]=new TH1D(Form("hMcOut_%.2f-%.2f-%s",CentMin[cent],CentMax[cent],Pid[pid]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
	hMcOut[cent][pid]->Sumw2();
	hMcOutNeg[cent][pid]=new TH1D(Form("hMcOutNeg_%.2f-%.2f-%s",CentMin[cent],CentMax[cent],Pid[pid]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
	hMcOutNeg[cent][pid]->Sumw2();
	hMcOutPos[cent][pid]=new TH1D(Form("hMcOutPos_%.2f-%.2f-%s",CentMin[cent],CentMax[cent],Pid[pid]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
	hMcOutPos[cent][pid]->Sumw2();

//				fListOfObjects->Add(hMcIn[cent][pid]);
//				fListOfObjects->Add(hMcInNeg[cent][pid]);
//				fListOfObjects->Add(hMcInPos[cent][pid]);
//				fListOfObjects->Add(hMcOut[cent][pid]);
//				fListOfObjects->Add(hMcOutNeg[cent][pid]);
//				fListOfObjects->Add(hMcOutPos[cent][pid]);

			}	// pid Eff

	hPiondEdx[cent]=new TH2D(Form("h%s_dEdx_%d",Pid[1],cent),
	";#it{p}_{T} (GeV/#it{c}); #it{d}E/#it{d}x",
	nPtBins,ptBins,fdEdxHigh-fdEdxLow,fdEdxLow,fdEdxHigh);
	hPiondEdx[cent]->Sumw2();
	hKaondEdx[cent]=new TH2D(Form("h%s_dEdx_%d",Pid[2],cent),
	";#it{p}_{T} (GeV/#it{c}); #it{d}E/#it{d}x",
	nPtBins,ptBins,fdEdxHigh-fdEdxLow,fdEdxLow,fdEdxHigh);
	hKaondEdx[cent]->Sumw2();
	hProtondEdx[cent]=new TH2D(Form("h%s_dEdx_%d",Pid[3],cent),
	";#it{p}_{T} (GeV/#it{c}); #it{d}E/#it{d}x",
	nPtBins,ptBins,fdEdxHigh-fdEdxLow,fdEdxLow,fdEdxHigh);
	hProtondEdx[cent]->Sumw2();
	fListOfObjects->Add(hPiondEdx[cent]);
	fListOfObjects->Add(hKaondEdx[cent]);
	fListOfObjects->Add(hProtondEdx[cent]);

		}	// cent Eff

	for(Int_t i_cent=0; i_cent<10; ++i_cent){
	for(Int_t pid=0; pid<7; ++pid){
	for(Int_t q=0; q<3; ++q){
	hDCApTPrim[i_cent][pid][q] = 0;
	hDCApTPrim[i_cent][pid][q] = new TH2D(Form("hDCA_%s%sPrimcent%d",Pid[pid],Q[q],i_cent),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
	hDCApTPrim[i_cent][pid][q]->Sumw2();
	hDCApTWDec[i_cent][pid][q] = 0;
	hDCApTWDec[i_cent][pid][q] = new TH2D(Form("hDCA_%s%sWDeccent%d",Pid[pid],Q[q],i_cent),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
	hDCApTWDec[i_cent][pid][q]->Sumw2();
	hDCApTMate[i_cent][pid][q] = 0;
	hDCApTMate[i_cent][pid][q] = new TH2D(Form("hDCA_%s%sMatecent%d",Pid[pid],Q[q],i_cent),"from material; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, ndcaBins, dcaBins );
	hDCApTMate[i_cent][pid][q]->Sumw2();

		// narrower bin size
	hDCApTPrim2[i_cent][pid][q] = 0;
	hDCApTPrim2[i_cent][pid][q] = new TH2D(Form("hDCA2_%s%sPrimcent%d",Pid[pid],Q[q],i_cent),"primaries; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
	hDCApTPrim2[i_cent][pid][q]->Sumw2();
	hDCApTWDec2[i_cent][pid][q] = 0;
	hDCApTWDec2[i_cent][pid][q] = new TH2D(Form("hDCA2_%s%sWDeccent%d",Pid[pid],Q[q],i_cent),"from weak decays; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
	hDCApTWDec2[i_cent][pid][q]->Sumw2();
	hDCApTMate2[i_cent][pid][q] = 0;
	hDCApTMate2[i_cent][pid][q] = new TH2D(Form("hDCA2_%s%sMatecent%d",Pid[pid],Q[q],i_cent),"from material; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", nPtBins, ptBins, 800, -4.0, 4.0 );
	hDCApTMate2[i_cent][pid][q]->Sumw2();

//					fListOfObjects->Add(hDCApTPrim[i_cent][pid][q]);
//					fListOfObjects->Add(hDCApTPrim2[i_cent][pid][q]);
//					fListOfObjects->Add(hDCApTWDec[i_cent][pid][q]);
//					fListOfObjects->Add(hDCApTWDec2[i_cent][pid][q]);
//					fListOfObjects->Add(hDCApTMate[i_cent][pid][q]);
//					fListOfObjects->Add(hDCApTMate2[i_cent][pid][q]);
				}	// charge MC
			}	// pid MC
		}	//cent DCA MC

		fVtxMC = new TH1F("fVtxMC","mc vtx", 120, -30, 30);
		fListOfObjects->Add(fVtxMC);

	}


	fEventCuts.AddQAplotsToList(fListOfObjects);
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

	fEvents->Fill(0.5,11);

	UInt_t fSelectMask= fInputHandler->IsEventSelected();
   	Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;
	fEvents->Fill(1.5,11);

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

	Int_t INEL = -1;
	      INEL = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0);
	if( INEL < 1 )
		return;

	fEvents->Fill(2.5,MultBin);
	fEvents->Fill(2.5,11);

  	if( utils->IsSPDClusterVsTrackletBG(fESD) ) 
		return;

	fEvents->Fill(3.5,MultBin);
	fEvents->Fill(3.5,11);

	if( fESD->IsPileupFromSPDInMultBins() )
		return;

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

	Float_t centrality = -10;
	if(fAnalysisPbPb){
		centrality = fEventCuts.GetCentrality(); /// Centrality calculated with the default estimator (V0M for LHC15o)

		if((centrality>fMaxCent)||(centrality<fMinCent))return;
		for(Int_t icent = 0; icent < nCent; ++icent){
			if(centrality > CentMin[icent] && centrality <= CentMax[icent]){
				cent = icent;
				fcent->Fill(icent+1);
				ProduceArrayTrksESD( esdEvent, cent );
				ProduceArrayV0ESD( esdEvent, cent );

				if(fAnalysisMC)
					ProcessMCTruthESD(cent);
				cout<<"Cent ::: "<<cent<<endl;
			}
		}

		fcent->Fill(11);

	}

	else{

      	Float_t V0MPercentile = -1;
        Float_t V0APercentile = -1;
        Float_t ADMPercentile = -1;

        AliMultSelection *MultSelection = (AliMultSelection*) esdEvent -> FindListObject("MultSelection");
        if(MultSelection-> IsEventSelected()){
                V0MPercentile = MultSelection->GetMultiplicityPercentile("V0M",false);
                V0APercentile = MultSelection->GetMultiplicityPercentile("V0A",false);
                ADMPercentile = MultSelection->GetMultiplicityPercentile("ADM",false);
        }

	cout << "V0MPercentile ~~~~~~~~~~~~ " << V0MPercentile<< endl;
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


	}



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

	cout<<"Cent Inside ProcessMCTruth ::: "<<Cent<<endl;
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

		//		if ( TMath::Abs(trackMC->Y()) > 0.5 )
		//			continue;

		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		hMcIn[Cent][0]->Fill(trackMC->Pt());
		hMcIn[Cent][pidCodeMC]->Fill(trackMC->Pt());

		hMcIn[10][0]->Fill(trackMC->Pt());
		hMcIn[10][pidCodeMC]->Fill(trackMC->Pt());

		if( chargeMC < 0 ){
			hMcInNeg[Cent][0]->Fill(trackMC->Pt());
			hMcInNeg[Cent][pidCodeMC]->Fill(trackMC->Pt());

			hMcInNeg[10][0]->Fill(trackMC->Pt());
			hMcInNeg[10][pidCodeMC]->Fill(trackMC->Pt());
		}
		else{
			hMcInPos[Cent][0]->Fill(trackMC->Pt());
			hMcInPos[Cent][pidCodeMC]->Fill(trackMC->Pt());

			hMcInPos[10][0]->Fill(trackMC->Pt());
			hMcInPos[10][pidCodeMC]->Fill(trackMC->Pt());
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

	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
		if(TMath::Abs(esdTrack->Eta()) > fEtaCut)
			continue;

		UInt_t selectDebug = 0;
		if ( fTrackFilterGolden ) {
			selectDebug = fTrackFilterGolden->IsSelected(esdTrack);
			if (!selectDebug) {
				continue;
			}
		}

		Short_t ncl = esdTrack->GetTPCsignalN();
		if(ncl<fNcl)
		continue;
		
		if(ncl < 90)cout<<" ¡¡¡¡¡¡¡¡¡¡¡     ncl < 90"<<endl;
		
		Double_t eta      = esdTrack->Eta();
		Double_t phi      = esdTrack->Phi();
		Double_t momentum = esdTrack->P();
		Double_t pt       = esdTrack->Pt();
		Float_t  dedx     = esdTrack->GetTPCsignal();
		Float_t  dedxUnc  = esdTrack->GetTPCsignal();

		Float_t dcaxy = 0.;
		Float_t dcaz = 0.;
		esdTrack->GetImpactParameters(dcaxy,dcaz);

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), Magf, fcutLow, fcutHigh))
			continue;

		if(fdEdxCalibrated){
			if(eta < 0){
				if(fLHC16l == 1)
				dedx *= 50/EtaCalibrationNeg(0,eta);
				else
				dedx *= 50/EtaCalibrationNeg(1,eta);
			}
			else{
				if(fLHC16l == 1)
				dedx *= 50/EtaCalibrationPos(0,eta);
				else
				dedx *= 50/EtaCalibrationPos(1,eta);
			}
		}

		if(esdTrack->Charge() < 0.0){
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion),2))<2.0){
				hDCAxyVsPtPiNeg[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtPiNegC[Cent]->Fill(pt,dcaxy);

			}
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon),2))<2.0){
				hDCAxyVsPtKNeg[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtKNegC[Cent]->Fill(pt,dcaxy);

			}
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton),2))<2.0){
				hDCAxyVsPtPNeg[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtPNegC[Cent]->Fill(pt,dcaxy);
			}
		}
		else{
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion),2))<2.0){
				hDCAxyVsPtPiPos[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtPiPosC[Cent]->Fill(pt,dcaxy);
			}
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon),2))<2.0){
				hDCAxyVsPtKPos[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtKPosC[Cent]->Fill(pt,dcaxy);
			}
			if(TMath::Sqrt(TMath::Power(fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),2)+TMath::Power(fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton),2))<2.0){
				hDCAxyVsPtPPos[Cent]->Fill(pt,dcaxy);
				hDCAxyVsPtPPosC[Cent]->Fill(pt,dcaxy);
			}
		}


		Short_t pidCode     = 0;
		if(fAnalysisMC) {

			const Int_t label = TMath::Abs(esdTrack->GetLabel());
			TParticle* mcTrack = 0;
			mcTrack = fMCStack->Particle(label);

			if (mcTrack){

				if( esdTrack->Charge()==0 )
					continue;

				Int_t pdgCode = mcTrack->GetPdgCode();
				pidCode = GetPidCode(pdgCode);

				if( fMCStack->IsPhysicalPrimary(label) ){

					if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) ){
						hMcOut[Cent][0]->Fill(esdTrack->Pt());
						hMcOut[Cent][pidCode]->Fill(esdTrack->Pt());

						hMcOut[10][0]->Fill(esdTrack->Pt());
						hMcOut[10][pidCode]->Fill(esdTrack->Pt());
 
						if(pidCode == 1){
                                                hPiondEdx[Cent]->Fill(esdTrack->Pt(),esdTrack->GetTPCsignal());
                                                hPiondEdx[10]->Fill(esdTrack->Pt(),esdTrack->GetTPCsignal());}
                                                else if(pidCode == 2){
                                                hKaondEdx[Cent]->Fill(esdTrack->Pt(),esdTrack->GetTPCsignal());
                                                hKaondEdx[10]->Fill(esdTrack->Pt(),esdTrack->GetTPCsignal());}
                                                else if(pidCode == 3){
                                                hProtondEdx[Cent]->Fill(esdTrack->Pt(),esdTrack->GetTPCsignal());
                                                hProtondEdx[10]->Fill(esdTrack->Pt(),esdTrack->GetTPCsignal());}
                                                else{
                                                continue;}
					}

					hDCApTPrim[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTPrim[Cent][pidCode][0]->Fill(pt,dcaxy);

					hDCApTPrim2[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTPrim2[Cent][pidCode][0]->Fill(pt,dcaxy);

					if( esdTrack->Charge() < 0.0 ){

						if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) )  {
							hMcOutNeg[Cent][0]->Fill(esdTrack->Pt());
							hMcOutNeg[Cent][pidCode]->Fill(esdTrack->Pt());

							hMcOutNeg[10][0]->Fill(esdTrack->Pt());
							hMcOutNeg[10][pidCode]->Fill(esdTrack->Pt());
						}

						hDCApTPrim[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTPrim[Cent][pidCode][1]->Fill(pt,dcaxy);

						hDCApTPrim2[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTPrim2[Cent][pidCode][1]->Fill(pt,dcaxy);
					}
					else{

						if( TMath::Abs(dcaxy) < GetMaxDCApTDep(fcutDCAxy,pt) ){
							hMcOutPos[Cent][0]->Fill(esdTrack->Pt());
							hMcOutPos[Cent][pidCode]->Fill(esdTrack->Pt());

							hMcOutPos[10][0]  ->Fill(esdTrack->Pt());
							hMcOutPos[10][pidCode]  ->Fill(esdTrack->Pt());
						}

						hDCApTPrim[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTPrim[Cent][pidCode][2]->Fill(pt,dcaxy);

						hDCApTPrim2[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTPrim2[Cent][pidCode][2]->Fill(pt,dcaxy);
					}
				}	// Primary particles MC

				if( fMCStack->IsSecondaryFromWeakDecay(label) ){

					hDCApTWDec[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTWDec[Cent][pidCode][0]->Fill(pt,dcaxy);

					hDCApTWDec2[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTWDec2[Cent][pidCode][0]->Fill(pt,dcaxy);

					if( esdTrack->Charge() < 0.0 ){

						hDCApTWDec[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTWDec[Cent][pidCode][1]->Fill(pt,dcaxy);

						hDCApTWDec2[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTWDec2[Cent][pidCode][1]->Fill(pt,dcaxy);

					}
					else{

						hDCApTWDec[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTWDec[Cent][pidCode][2]->Fill(pt,dcaxy);

						hDCApTWDec2[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTWDec2[Cent][pidCode][2]->Fill(pt,dcaxy);

					}
				}	// Weak Decay MC

				if( fMCStack->IsSecondaryFromMaterial(label) ){

					hDCApTMate[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTMate[Cent][pidCode][0]->Fill(pt,dcaxy);

					hDCApTMate2[Cent][0][0]->Fill(pt,dcaxy);
					hDCApTMate2[Cent][pidCode][0]->Fill(pt,dcaxy);

					if( esdTrack->Charge() < 0.0 ){

						hDCApTMate[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTMate[Cent][pidCode][1]->Fill(pt,dcaxy);

						hDCApTMate2[Cent][0][1]->Fill(pt,dcaxy);
						hDCApTMate2[Cent][pidCode][1]->Fill(pt,dcaxy);

					}
					else{

						hDCApTMate[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTMate[Cent][pidCode][2]->Fill(pt,dcaxy);

						hDCApTMate2[Cent][0][2]->Fill(pt,dcaxy);
						hDCApTMate2[Cent][pidCode][2]->Fill(pt,dcaxy);

					}
				}	// Material Inte MC
			}	//mcTrack
		}	//fAnalysis MC

		//==========================  DCAxy cut

		cout<<"DCA before::  "<<dcaxy<<endl;
		if( TMath::Abs(dcaxy) > GetMaxDCApTDep(fcutDCAxy,pt) )
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

		if(!fdEdxCalibrated){
			//			cout<<"PreCalibration"<<endl;
			if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
				if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
					hMIPVsEta[Cent]->Fill(eta,dedxUnc);
					pMIPVsEta[Cent]->Fill(eta,dedxUnc);
					hMIPVsEta[10]->Fill(eta,dedxUnc);
					pMIPVsEta[10]->Fill(eta,dedxUnc);
				}
				if( dedxUnc > 70 && dedxUnc < 90 ){
					if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta[Cent]->Fill(eta,dedxUnc);
					pPlateauVsEta[Cent]->Fill(eta,dedxUnc);
					hPlateauVsEta[10]->Fill(eta,dedxUnc);
					pPlateauVsEta[10]->Fill(eta,dedxUnc);
					}
				}
			}
		}
		else{
			//			cout<<"PostCalibration"<<endl;
			if( momentum <= 0.6 && momentum >= 0.4 ){//only p:0.4-0.6 GeV, pion MIP
				if( dedxUnc < fDeDxMIPMax && dedxUnc > fDeDxMIPMin ){
					hMIPVsEta[Cent]->Fill(eta,dedx);
					pMIPVsEta[Cent]->Fill(eta,dedx);
					hMIPVsEta[10]->Fill(eta,dedx);
					pMIPVsEta[10]->Fill(eta,dedx);
				}
				if( dedxUnc > 70 && dedxUnc < 90 ){
					if(TMath::Abs(beta-1)<0.1){
					hPlateauVsEta[Cent]->Fill(eta,dedx);
					pPlateauVsEta[Cent]->Fill(eta,dedx);
					hPlateauVsEta[10]->Fill(eta,dedx);
					pPlateauVsEta[10]->Fill(eta,dedx);
					}
				}
			}
		}

		hPtAll[Cent]->Fill(pt);
		hPtAll[10]->Fill(pt);

		if(esdTrack->Charge() < 0.)
			hPtAllNeg[Cent]->Fill(pt);
		
		else
			hPtAllPos[Cent]->Fill(pt);
		

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

		if(beta>1){
			histPiTof[Cent][nh]->Fill(momentum, dedx);
			histPiTof[10][nh]->Fill(momentum, dedx);
		}

		if( momentum <= 0.6 && momentum >= 0.4  ){

			if( dedx < fDeDxMIPMax && dedx > fDeDxMIPMin ){
				hMIPVsPhi[Cent][nh]->Fill(phi,dedx);
				pMIPVsPhi[Cent][nh]->Fill(phi,dedx);
				hMIPVsPhi[10][nh]->Fill(phi,dedx);
				pMIPVsPhi[10][nh]->Fill(phi,dedx);

				hMIPVsNch[nh]->Fill(multTPC,dedx);
				pMIPVsNch[nh]->Fill(multTPC,dedx);

				hMIPVsV0M[nh]->Fill(V0MPer,dedx);   // Put here V0M Multiplicity and Declare the Histos. !!!
				pMIPVsV0M[nh]->Fill(V0MPer,dedx);

			}
			if( dedx > 70 && dedx < 90 ){
				if(TMath::Abs(beta-1)<0.1){
				hPlateauVsPhi[Cent][nh]->Fill(phi,dedx);
				pPlateauVsPhi[Cent][nh]->Fill(phi,dedx);
				hPlateauVsPhi[10][nh]->Fill(phi,dedx);
				pPlateauVsPhi[10][nh]->Fill(phi,dedx);

				}
			}
		}

		hPtVsP[Cent][nh]->Fill(momentum,pt);
		hPtVsP[10][nh]->Fill(momentum,pt);
		hDeDxVsP[Cent][nh]->Fill(momentum,dedx);
		hDeDxVsP[10][nh]->Fill(momentum,dedx);

		if(esdTrack->Charge() < 0.){
			hPtNeg[Cent][nh]->Fill(pt);

			hnSigmaPiNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaKNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaPNeg[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
		}
		else{
			hPtPos[Cent][nh]->Fill(pt);

			hnSigmaPiPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hnSigmaKPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hnSigmaPPos[Cent][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
		}

	}//end of track loop


}

//----------------------------------------------------------------------------------

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
//		if(TMath::Abs(pTrack->Eta()) > fEtaCut || TMath::Abs(nTrack->Eta()) > fEtaCut)
//			continue;

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

						       if(track->GetTPCsignalN()<70)continue;
						       Double_t phi     = track->Phi();

						       if(!PhiCut(track->Pt(), phi, track->Charge(), Magf, fcutLow, fcutHigh))
							       continue;

						       Double_t eta      = track->Eta();
						       Double_t momentum = track->P();
						       Double_t dedx     = track->GetTPCsignal();
						       Double_t dedxUnc  = track->GetTPCsignal();


						       if(fdEdxCalibrated){
							       if(eta < 0){
									if(fLHC16l == 1)
								        dedx   *= 50.0/EtaCalibrationNeg(0,eta);
									else
								        dedx   *= 50.0/EtaCalibrationNeg(1,eta);
							       }
							       else{
									if(fLHC16l == 1)
								        dedx   *= 50.0/EtaCalibrationPos(0,eta);
									else
								        dedx   *= 50.0/EtaCalibrationPos(1,eta);
							       }
						       }


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
							       if(nTrack->Eta() > 0){
								if(fLHC16l == 1){
								       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationPosEl(0,nTrack->Eta())) < 5)
									       fillPos = kTRUE;
									}
								else{
								       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationPosEl(1,nTrack->Eta())) < 5)
									       fillPos = kTRUE;
		cout << "-------------------LHC16k"  << endl;
									}

								}
							       if(nTrack->Eta() < 0){
								if(fLHC16l == 1){
								       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationNegEl(0,nTrack->Eta())) < 5)
									       fillPos = kTRUE;
								}
								else{
								       if( TMath::Abs(nTrack->GetTPCsignal()-EtaCalibrationNegEl(1,nTrack->Eta())) < 5)
									       fillPos = kTRUE;
		cout << "-------------------LHC16k"  << endl;
								}
							}

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
						       if(eta < 0){
								if(fLHC16l == 1)
							        dedx *= 50/EtaCalibrationNeg(0,eta);
								else
							        dedx *= 50/EtaCalibrationNeg(1,eta);
							}
						       else{
								if(fLHC16l == 1)
							        dedx *= 50/EtaCalibrationPos(0,eta);
								else
							        dedx *= 50/EtaCalibrationPos(1,eta);
							}
					       }


//					       if(track->GetTPCsignalN()<=70)continue;

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


					       cout<<"Cent :: "<<Cent<<"Eta  :: "<<eta<<"dedx :: "<<dedx<<endl;
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
Double_t AliAnalysisTaskPPvsMult::EtaCalibrationNeg( const Int_t Cent, const Double_t eta){


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

}
