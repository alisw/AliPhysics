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

#include "AliAnalysisTaskSpherocity.h"

// ROOT includes
#include <TList.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
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
//#include <AliSpherocityUtils.h>

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


const Double_t AliAnalysisTaskSpherocity::fgkClight = 2.99792458e-2;
Float_t MAGF                    = 1;
const Int_t nHists              = 4;
const Int_t nCent               = 11;
const Double_t CentMin[nCent]   = {0.0,1.0,5.0 ,10.0,15.0,20.0,30.0,40.0,50.0,70.0,0.0};
const Double_t CentMax[nCent]   = {1.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,70.0,100.0,100.0};
static const double C_Value = TMath::C()*(1.e2/1.e12); // cm/ps

const Char_t *So[3]         = {"Jetty", "Isotropic", "Reference"};

ClassImp(AliAnalysisTaskSpherocity)
	AliAnalysisTaskSpherocity::AliAnalysisTaskSpherocity():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fAOD(0x0),
		fEventCuts(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fMCArray(0x0),
		fPIDResponse(0x0),
		fTrackFilterGolden(0x0),
		fTrackFilter(0x0),
		utils(0x0),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fisV0Mestimator(kFALSE),
		fRandom(0x0),
		fNcl(70),
		fEtaCut(0.9),
		fDeDxMIPMin(40),
		fDeDxMIPMax(60),
		fdEdxHigh(200),
		fdEdxLow(40),
		fJettyCutOff(0.503),
		fIsotrCutOff(0.774),
		fJettyCutOff_0(0.1),
		fJettyCutOff_1(0.3),
		fJettyCutOff_2(0.4),
		fIsotrCutOff_0(0.95),
		fIsotrCutOff_1(0.9),
		fIsotrCutOff_2(0.85),
		fMinMult(10),
		fNrec(0),
		fSizeStep(0.1),
		fMcProcessType(-999),
		fTriggeredEventMB(-999),
		fVtxStatus(-999),
		fZvtx(-999),
		fZvtxMC(-999),
		fListOfObjects(0),
		fEvents(0x0),
		fdEdxCalibrated(kFALSE),
		fPeriod("16l"),
		hphiso(0x0),
		hetaso(0x0),
		hPtTruthVsPtRec(0x0),
		hPtTruthVsPtRecJetty(0x0),
		hPtTruthVsPtRecIsotr(0x0),
		hTruthPhiSo(0x0),
		hTruthEtaSo(0x0),
		hSOtvsTrks(0x0),
		hSOtvsTrkst(0x0),
		hSOtvsV0M(0x0),
		hSOrvsV0M(0x0),
		hSOrvsTrks(0x0),
		hRefMultVsRefMultPer(0x0),
		fEtaCalibrationNeg(0x0),
		fEtaCalibrationPos(0x0),
		fcutDCAxy(0x0),
		fcutLow(0x0),
		fcutHigh(0x0)

{



	for(Int_t i = 0; i< 2; ++i){

		for(Int_t so = 0; so < 4; ++so){

			hPtAll[i][so]=0;
			hPtpos_TPC[i][so]=0;
			hPtneg_TPC[i][so]=0;
			hPtpos_TOF[i][so]=0;
			hPtneg_TOF[i][so]=0;

			for(Int_t j = 0; j < nHists; ++j){

				hDeDxVsP[i][so][j]=0;
				hnSigmaPiPos[i][so][j]=0;
				hnSigmaKPos[i][so][j]=0;
				hnSigmaPPos[i][so][j]=0;
				hnSigmaPiNeg[i][so][j]=0;
				hnSigmaKNeg[i][so][j]=0;
				hnSigmaPNeg[i][so][j]=0;

				hBetavsPpos[i][so][j]=0;
				hBetavsPneg[i][so][j]=0;

				hPtpos_TPC_Eta[i][so][j]=0;
				hPtneg_TPC_Eta[i][so][j]=0;
				hPtpos_TOF_Eta[i][so][j]=0;
				hPtneg_TOF_Eta[i][so][j]=0;
				hPpos_TOF_Eta[i][so][j]=0;
				hPneg_TOF_Eta[i][so][j]=0;
			}
		}
	}

	hPtAllSoInt=0;
	hPtpos_TPCSoInt=0;
	hPtneg_TPCSoInt=0;
	hPtpos_TOFSoInt=0;
	hPtneg_TOFSoInt=0;

	for(Int_t j=0; j<nHists; ++j){

		hDeDxVsPSoInt[j]=0;
		hnSigmaPiPosSoInt[j]=0;
		hnSigmaKPosSoInt[j]=0;
		hnSigmaPPosSoInt[j]=0;
		hnSigmaPiNegSoInt[j]=0;
		hnSigmaKNegSoInt[j]=0;
		hnSigmaPNegSoInt[j]=0;

		hBetavsPposSoInt[j]=0;
		hBetavsPnegSoInt[j]=0;

		hPtpos_TPC_EtaSoInt[j]=0;
		hPtneg_TPC_EtaSoInt[j]=0;
		hPtpos_TOF_EtaSoInt[j]=0;
		hPtneg_TOF_EtaSoInt[j]=0;
		hPpos_TOF_EtaSoInt[j]=0;
		hPneg_TOF_EtaSoInt[j]=0;
		hPtVsP[j]=0;
	}

	//	hTruthEtaSo = 0;
	//	hTruthPhiSo = 0;
	//default constructor
	for(Int_t cent=0;cent<nCent;++cent){
		//		hSOtVsSOm[cent]   = 0;
		for(Int_t pid=0;pid<7;++pid){
			for(Int_t so=0;so<3;++so){
				hMcIn[cent][pid][so]     = 0;
				hMcOut[cent][pid][so]    = 0;
				hMcInNeg[cent][pid][so]  = 0;
				hMcInPos[cent][pid][so]  = 0;
				hMcOutNeg[cent][pid][so] = 0;
				hMcOutPos[cent][pid][so] = 0;
			}
		}
	}

}

AliAnalysisTaskSpherocity::AliAnalysisTaskSpherocity(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fAOD(0x0),
	fEventCuts(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fMCArray(0x0),
	fPIDResponse(0x0),
	fTrackFilterGolden(0x0),
	fTrackFilter(0x0),
	utils(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fisV0Mestimator(kFALSE),
	fRandom(0x0),
	fNcl(70),
	fEtaCut(0.9),
	fDeDxMIPMin(40),
	fDeDxMIPMax(60),
	fdEdxHigh(200),
	fdEdxLow(40),
	fJettyCutOff(0.503),
	fIsotrCutOff(0.774),
	fJettyCutOff_0(0.1),
	fJettyCutOff_1(0.3),
	fJettyCutOff_2(0.4),
	fIsotrCutOff_0(0.95),
	fIsotrCutOff_1(0.9),
	fIsotrCutOff_2(0.85),
	fMinMult(10),
	fNrec(0),
	fSizeStep(0.1),
	fMcProcessType(-999),
	fTriggeredEventMB(-999),
	fVtxStatus(-999),
	fZvtx(-999),
	fZvtxMC(-999),
	fListOfObjects(0), 
	fEvents(0x0),
	fdEdxCalibrated(kFALSE),
	fPeriod("16l"),
	hphiso(0x0),
	hetaso(0x0),
	hPtTruthVsPtRec(0x0),
	hPtTruthVsPtRecJetty(0x0),
	hPtTruthVsPtRecIsotr(0x0),
	hTruthPhiSo(0x0),
	hTruthEtaSo(0x0),
	hSOtvsTrks(0x0),
	hSOtvsTrkst(0x0),
	hSOtvsV0M(0x0),
	hSOrvsV0M(0x0),
	hSOrvsTrks(0x0),
	hRefMultVsRefMultPer(0x0),
	fEtaCalibrationNeg(0x0),
	fEtaCalibrationPos(0x0),
	fcutDCAxy(0x0),
	fcutLow(0x0),
	fcutHigh(0x0)
{

	for(Int_t i = 0; i< 2; ++i){

		for(Int_t so = 0; so < 4; ++so){

			hPtAll[i][so]=0;
			hPtpos_TPC[i][so]=0;
			hPtneg_TPC[i][so]=0;
			hPtpos_TOF[i][so]=0;
			hPtneg_TOF[i][so]=0;

			for(Int_t j = 0; j < nHists; ++j){

				hDeDxVsP[i][so][j]=0;
				hnSigmaPiPos[i][so][j]=0;
				hnSigmaKPos[i][so][j]=0;
				hnSigmaPPos[i][so][j]=0;
				hnSigmaPiNeg[i][so][j]=0;
				hnSigmaKNeg[i][so][j]=0;
				hnSigmaPNeg[i][so][j]=0;

				hBetavsPpos[i][so][j]=0;
				hBetavsPneg[i][so][j]=0;

				hPtpos_TPC_Eta[i][so][j]=0;
				hPtneg_TPC_Eta[i][so][j]=0;
				hPtpos_TOF_Eta[i][so][j]=0;
				hPtneg_TOF_Eta[i][so][j]=0;
				hPpos_TOF_Eta[i][so][j]=0;
				hPneg_TOF_Eta[i][so][j]=0;
			}
		}
	}

	hPtAllSoInt=0;
	hPtpos_TPCSoInt=0;
	hPtneg_TPCSoInt=0;
	hPtpos_TOFSoInt=0;
	hPtneg_TOFSoInt=0;

	for(Int_t j=0; j<nHists; ++j){

		hDeDxVsPSoInt[j]=0;
		hnSigmaPiPosSoInt[j]=0;
		hnSigmaKPosSoInt[j]=0;
		hnSigmaPPosSoInt[j]=0;
		hnSigmaPiNegSoInt[j]=0;
		hnSigmaKNegSoInt[j]=0;
		hnSigmaPNegSoInt[j]=0;

		hBetavsPposSoInt[j]=0;
		hBetavsPnegSoInt[j]=0;

		hPtpos_TPC_EtaSoInt[j]=0;
		hPtneg_TPC_EtaSoInt[j]=0;
		hPtpos_TOF_EtaSoInt[j]=0;
		hPtneg_TOF_EtaSoInt[j]=0;
		hPpos_TOF_EtaSoInt[j]=0;
		hPneg_TOF_EtaSoInt[j]=0;
		hPtVsP[j]=0;
	}

	// Default constructor (should not be used)
	for(Int_t cent=0; cent<nCent; ++cent){
		//		hSOtVsSOm[cent] = 0;
		for(Int_t pid=0; pid<7; ++pid){
			for(Int_t so=0; so<3; ++so){
				hMcIn[cent][pid][so]=0;
				hMcOut[cent][pid][so]=0;
				hMcInNeg[cent][pid][so]=0;
				hMcInPos[cent][pid][so]=0;
				hMcOutNeg[cent][pid][so]=0;
				hMcOutPos[cent][pid][so]=0;
			}
		}
	}

	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());//esto es nuevo
}




AliAnalysisTaskSpherocity::~AliAnalysisTaskSpherocity() {
	//
	// Destructor
	//

}
//______________________________________________________________________________
void AliAnalysisTaskSpherocity::UserCreateOutputObjects()
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

	// Definition of trackcuts
	if(!fTrackFilter){
		fTrackFilter = new AliAnalysisFilter("trackFilter2015");
		SetTrackCutsSpherocity(fTrackFilter);
	}

	if(!fTrackFilterGolden){

		fTrackFilterGolden = new AliESDtrackCuts("fTrackFilterGolden");
		fTrackFilterGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
		fTrackFilterGolden->SetEtaRange(-0.8,0.8);

	}

	/*if(!fSpheroUtils){
	  fSpheroUtils = new AliSpherocityUtils();
	  fSpheroUtils->SetTrackFilter(fTrackFilter);
	  fSpheroUtils->SetMinMult(10);
	  fSpheroUtils->Init();

	  printf("JettyCutOff  =  %f       IsotrCutOff  =  %f\n",fJettyCutOff,fIsotrCutOff);
	  }*/

	//OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner(kTRUE);

	//
	// Histograms
	//

	const int nBinsPer = 10;
	float BinsPer[nBinsPer+1] = {0.0,1.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,70.0,100.0};

	const int nBins = 12;
	float Bins[nBins+1] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0}; 

	fEvents = new TH2F( "fEvents", ";Evt. Sel.;Mult. Per.",nBins,Bins,nBinsPer,BinsPer);
	fEvents->GetXaxis()->SetBinLabel(1, "Processed");
	fEvents->GetXaxis()->SetBinLabel(2, "Trigger");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(3, "IsPileUpFromSPDinMultBins");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(4, "DAQ");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(5, "BG");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(6, "INEL>0");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(7, "VtxRes&Proximity");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(8, "|Vtz|<10cm");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(9, "non-Selected So");//NotinVertexcut");
	fEvents->GetXaxis()->SetBinLabel(10, "Selected So");//NotinVertexcut");
	fEvents->GetYaxis()->SetBinLabel(1,Form("%.2f-%.2f",BinsPer[0],BinsPer[1]));
	fEvents->GetYaxis()->SetBinLabel(2,Form("%.2f-%.2f",BinsPer[1],BinsPer[2]));
	fEvents->GetYaxis()->SetBinLabel(3,Form("%.2f-%.2f",BinsPer[2],BinsPer[3]));
	fEvents->GetYaxis()->SetBinLabel(4,Form("%.2f-%.2f",BinsPer[3],BinsPer[4]));
	fEvents->GetYaxis()->SetBinLabel(5,Form("%.2f-%.2f",BinsPer[4],BinsPer[5]));
	fEvents->GetYaxis()->SetBinLabel(6,Form("%.2f-%.2f",BinsPer[5],BinsPer[6]));
	fEvents->GetYaxis()->SetBinLabel(7,Form("%.2f-%.2f",BinsPer[6],BinsPer[7]));
	fEvents->GetYaxis()->SetBinLabel(8,Form("%.2f-%.2f",BinsPer[7],BinsPer[8]));
	fEvents->GetYaxis()->SetBinLabel(9,Form("%.2f-%.2f",BinsPer[8],BinsPer[9]));
	fEvents->GetYaxis()->SetBinLabel(10,Form("%.2f-%.2f",BinsPer[9],BinsPer[10]));
	//	fEvents->GetYaxis()->SetBinLabel(11,"0.0-100.0");
	fListOfObjects->Add(fEvents);

	//	hMult = new TH1F("hMult","Mult of events with SO measured",100,0,100);
	//	fListOfObjects->Add(hMult);

	//	fTrcksVsTrklets = new TH2F("fTrcksVsTrklets",";SPD_{Tracklets};Global Tracks",100,0,100,100,0,100);
	//	fListOfObjects->Add(fTrcksVsTrklets);

	//fcent=new TH1F("fcent","fcent",13,0,13);
	//fcentAfterPrimaries =new TH1F("fcentAfterPrimaries","fcentAfterPrimaries",13,0,13);
	//	fcentAfterV0s =new TH1F("fcentAfterV0s","fcentAfterV0s",13,0,13);
	//	fListOfObjects->Add(fcent);
	//	fListOfObjects->Add(fcentAfterPrimaries);
	//	fListOfObjects->Add(fcentAfterV0s);

	//	const Int_t nDeltaPiBins   = 80;
	//	const Double_t deltaPiLow  = 20;
	//	const Double_t deltaPiHigh = 100;

	const Char_t *Pid[7]       = {"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};

	const int nPtBins = 63;
	float ptBins[nPtBins+1] = {
		0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35,
		0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
		0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
		1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40,
		3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00,
		9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 18.00,
		20.00,22.00,24.00,26.00,30.00};

	//	const int ndEdxBins   = 200;
	//	double dEdxBins[ndEdxBins+1] = { 0 };
	//	for(int i = 0; i <= ndEdxBins; ++i){
	//		dEdxBins[i] = fdEdxLow+i*1.0;
	//	}

	const int nBinsNsigma = 50;
	float BinsNsigma[nBinsNsigma+1] = {0};

	for(int i = 0; i <= nBinsNsigma; ++i){
		BinsNsigma[i] = -10.0+i*0.4;
	}

	const int nBetaBins   = 100;
	float BetaBins[nBetaBins+1] = { 0 };
	for(int i = 0; i <= nBetaBins; ++i){
		BetaBins[i] = 0.2+((double)i)/100.0;
	}

	//	const int nPtBinsV0s = 25;
	//	double ptBinsV0s[nPtBinsV0s+1] = {
	//		0.0 , 0.1 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 ,
	//		1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 5.0 , 7.0 ,
	//		9.0 , 12.0, 15.0, 20.0 };

	printf("==============================\n");
	printf("Running the Fixed Code\n");
	printf("Jetty cut: %f\n",fJettyCutOff);
	printf("Isotropic cut: %f\n",fIsotrCutOff);
	printf("==============================\n");

	const int ndEdxBins = fdEdxHigh-fdEdxLow;
	float dEdxBins[ndEdxBins+1];

	for(int i = fdEdxLow; i <= fdEdxHigh; ++i){
		dEdxBins[i-fdEdxLow] = i;
		//        printf("edges :: %f\n",binsdEdx[i-fdEdxLow]);
	}

	const Char_t* ending[nHists] = {"02", "24", "46", "68"};
	////	const Char_t* LatexEta[nHists] = {"|#eta|<0.2", "0.2<|#eta|<0.4", "0.4<|#eta|<0.6", "0.6<|#eta|<0.8" };

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);


	fEtaCalibrationPos = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	fEtaCalibrationNeg = new TF1("fDeDxVsEtaNeg", "pol7", -1.0, 0.0);

	//	felededxfitPos     = new TF1("felededxfitPos", "pol4", 0.0, 1.0);
	//	felededxfitNeg     = new TF1("felededxfitNeg", "pol4", -1.0, 0.0);

	hphiso = new TH1D("hphiso","spherocity; #phi; counts",64,0.0,2*TMath::Pi());
	fListOfObjects->Add(hphiso);

	hetaso = new TH1D("hetaso","spherocity; #eta; counts",40,-1.0,1.0);
	fListOfObjects->Add(hetaso);

	///	Int_t nPhiBins = 36;

	hSOrvsV0M  = new TH2F("hSOrVsV0M","Measured SO vs V0M Per.;V0M Per.;#it{S}_{O} Reconstructed",100,0,100,1000,0,1);
	fListOfObjects->Add(hSOrvsV0M);

	hSOrvsTrks  = new TH2F("hSOrVsTrks","Measured SO vs Measured Ref. Mult. |#eta|<0.8;Reference mult. (|#eta|<0.8);#it{S}_{O} Reconstructed",100, 0, 100, 1000, 0, 1);
	fListOfObjects->Add(hSOrvsTrks);

	hRefMultVsRefMultPer = new TH2F("hRefMultVsRefMultPer","Ref Mult. vs Ref. Mult. Per. |#eta|<0.8;Ref. Mult Per.;Ref. Mult",100, 0, 100, 100, 0, 100);
	fListOfObjects->Add(hRefMultVsRefMultPer);

	for(Int_t topology = 0; topology < 2; ++topology){
		for(Int_t so = 0; so < 4; ++so){
			hPtAll[topology][so] = new TH2F(Form("hPt_rTPC_%s_%d",So[topology],so),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
			hPtpos_TPC[topology][so] = new TH2F(Form("hPt_pos_TPC_%s_%d",So[topology],so),";#it{p}_{T};MultPer",nPtBins,ptBins,nBinsPer,BinsPer);
			hPtneg_TPC[topology][so] = new TH2F(Form("hPt_neg_TPC_%s_%d",So[topology],so),";#it{p}_{T};MultPer",nPtBins,ptBins,nBinsPer,BinsPer);
			hPtpos_TOF[topology][so] = new TH2F(Form("hPt_pos_TOF_%s_%d",So[topology],so),";#it{p}_{T};MultPer",nPtBins,ptBins,nBinsPer,BinsPer);
			hPtneg_TOF[topology][so] = new TH2F(Form("hPt_neg_TOF_%s_%d",So[topology],so),";#it{p}_{T};MultPer",nPtBins,ptBins,nBinsPer,BinsPer);

			for(Int_t j = 0; j < nHists; j++) {            

				hDeDxVsP[topology][so][j] = new TH3F(Form("hDeDxVsP_%s_%d_%s",So[topology],so,ending[j]), ";#it{p} [GeV/c]; dE/dx; MultPer",nPtBins,ptBins,ndEdxBins,dEdxBins,nBinsPer,BinsPer);
				hnSigmaPiPos[topology][so][j] = new TH3F(Form("hnSigma_Pion_pos_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};nSigmaPiPos;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
				hnSigmaKPos[topology][so][j] = new TH3F(Form("hnSigma_Kaon_pos_%s_%d_%s",So[topology],so,ending[j]), ";#it{p}_{T};nSigmaKPos;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
				hnSigmaPPos[topology][so][j] = new TH3F(Form("hnSigma_Proton_pos_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};nSigmaPPos;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
				hnSigmaPiNeg[topology][so][j] = new TH3F(Form("hnSigma_Pion_neg_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};nSigmaPiNeg;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
				hnSigmaKNeg[topology][so][j] = new TH3F(Form("hnSigma_Kaon_neg_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};nSigmaKNeg;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
				hnSigmaPNeg[topology][so][j] = new TH3F(Form("hnSigma_Proton_neg_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};nSigmaPNeg;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
				hBetavsPpos[topology][so][j] = new TH3F(Form("hBetavsP_pos_%s_%d_%s",So[topology],so,ending[j]),";#it{p};#beta;MultPer",nPtBins,ptBins,nBetaBins,BetaBins,nBinsPer,BinsPer);
				hBetavsPneg[topology][so][j] = new TH3F(Form("hBetavsP_neg_%s_%d_%s",So[topology],so,ending[j]),";#it{p};#beta;MultPer",nPtBins,ptBins,nBetaBins,BetaBins,nBinsPer,BinsPer);
				hPtpos_TPC_Eta[topology][so][j] = new TH2F(Form("hPt_pos_TPC_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
				hPtneg_TPC_Eta[topology][so][j] = new TH2F(Form("hPt_neg_TPC_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
				hPtpos_TOF_Eta[topology][so][j] = new TH2F(Form("hPt_pos_TOF_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
				hPtneg_TOF_Eta[topology][so][j] = new TH2F(Form("hPt_neg_TOF_%s_%d_%s",So[topology],so,ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
				hPpos_TOF_Eta[topology][so][j] = new TH2F(Form("hP_pos_TOF_%s_%d_%s",So[topology],so,ending[j]),";#it{p};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
				hPneg_TOF_Eta[topology][so][j] = new TH2F(Form("hP_neg_TOF_%s_%d_%s",So[topology],so,ending[j]),";#it{p};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
			}
		}
	} // topology loop

	for(Int_t j = 0; j < nHists; j++) {            

		hDeDxVsPSoInt[j] = new TH3F(Form("hDeDxVsP_%s_%s",So[2],ending[j]), ";#it{p} [GeV/c]; dE/dx; MultPer",nPtBins,ptBins,ndEdxBins,dEdxBins,nBinsPer,BinsPer);
		hnSigmaPiPosSoInt[j] = new TH3F(Form("hnSigma_Pion_pos_%s_%s",So[2],ending[j]),";#it{p}_{T};nSigmaPiPos;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
		hnSigmaKPosSoInt[j] = new TH3F(Form("hnSigma_Kaon_pos_%s_%s",So[2],ending[j]), ";#it{p}_{T};nSigmaKPos;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
		hnSigmaPPosSoInt[j] = new TH3F(Form("hnSigma_Proton_pos_%s_%s",So[2],ending[j]),";#it{p}_{T};nSigmaPPos;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
		hnSigmaPiNegSoInt[j] = new TH3F(Form("hnSigma_Pion_neg_%s_%s",So[2],ending[j]),";#it{p}_{T};nSigmaPiNeg;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
		hnSigmaKNegSoInt[j] = new TH3F(Form("hnSigma_Kaon_neg_%s_%s",So[2],ending[j]),";#it{p}_{T};nSigmaKNeg;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
		hnSigmaPNegSoInt[j] = new TH3F(Form("hnSigma_Proton_neg_%s_%s",So[2],ending[j]),";#it{p}_{T};nSigmaPNeg;MultPer",nPtBins,ptBins,nBinsNsigma,BinsNsigma,nBinsPer,BinsPer);
		hBetavsPposSoInt[j] = new TH3F(Form("hBetavsP_pos_%s_%s",So[2],ending[j]),";#it{p};#beta;MultPer",nPtBins,ptBins,nBetaBins,BetaBins,nBinsPer,BinsPer);
		hBetavsPnegSoInt[j] = new TH3F(Form("hBetavsP_neg_%s_%s",So[2],ending[j]),";#it{p};#beta;MultPer",nPtBins,ptBins,nBetaBins,BetaBins,nBinsPer,BinsPer);
		hPtpos_TPC_EtaSoInt[j] = new TH2F(Form("hPt_pos_TPC_%s_%s",So[2],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
		hPtneg_TPC_EtaSoInt[j] = new TH2F(Form("hPt_neg_TPC_%s_%s",So[2],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
		hPtpos_TOF_EtaSoInt[j] = new TH2F(Form("hPt_pos_TOF_%s_%s",So[2],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
		hPtneg_TOF_EtaSoInt[j] = new TH2F(Form("hPt_neg_TOF_%s_%s",So[2],ending[j]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
		hPpos_TOF_EtaSoInt[j] = new TH2F(Form("hP_pos_TOF_%s_%s",So[2],ending[j]),";#it{p};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
		hPneg_TOF_EtaSoInt[j] = new TH2F(Form("hP_neg_TOF_%s_%s",So[2],ending[j]),";#it{p};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
		hPtVsP[j] = new TH2D(Form("hPtVsP_%s",ending[j]), ";#it{p} [GeV/c]; #it{p}_{T}", nPtBins, ptBins, nPtBins, ptBins);
	}// eta loop

	hPtAllSoInt = new TH2F(Form("hPt_rTPC_%s",So[2]),";#it{p}_{T};Counts",nPtBins,ptBins,nBinsPer,BinsPer);
	hPtpos_TPCSoInt = new TH2F(Form("hPt_pos_TPC_%s",So[2]),";#it{p}_{T};MultPer",nPtBins,ptBins,nBinsPer,BinsPer);
	hPtneg_TPCSoInt = new TH2F(Form("hPt_neg_TPC_%s",So[2]),";#it{p}_{T};MultPer",nPtBins,ptBins,nBinsPer,BinsPer);
	hPtpos_TOFSoInt = new TH2F(Form("hPt_pos_TOF_%s",So[2]),";#it{p}_{T};MultPer",nPtBins,ptBins,nBinsPer,BinsPer);
	hPtneg_TOFSoInt = new TH2F(Form("hPt_neg_TOF_%s",So[2]),";#it{p}_{T};MultPer",nPtBins,ptBins,nBinsPer,BinsPer);

	if(!fAnalysisMC){

		for(Int_t i = 0; i< 2; ++i ){

			for(Int_t so = 0; so < 4; ++so){

				fListOfObjects->Add(hPtAll[i][so]);
				fListOfObjects->Add(hPtneg_TPC[i][so]);
				fListOfObjects->Add(hPtpos_TPC[i][so]);
				fListOfObjects->Add(hPtneg_TOF[i][so]);
				fListOfObjects->Add(hPtpos_TOF[i][so]);

				for(Int_t j = 0; j < nHists; ++j){
					fListOfObjects->Add(hDeDxVsP[i][so][j]);
					fListOfObjects->Add(hnSigmaPiPos[i][so][j]);
					fListOfObjects->Add(hnSigmaPiNeg[i][so][j]);
					fListOfObjects->Add(hnSigmaKPos[i][so][j]);
					fListOfObjects->Add(hnSigmaKNeg[i][so][j]);
					fListOfObjects->Add(hnSigmaPPos[i][so][j]);
					fListOfObjects->Add(hnSigmaPNeg[i][so][j]);
					fListOfObjects->Add(hBetavsPneg[i][so][j]);
					fListOfObjects->Add(hBetavsPpos[i][so][j]);
					fListOfObjects->Add(hPtneg_TPC_Eta[i][so][j]);
					fListOfObjects->Add(hPtpos_TPC_Eta[i][so][j]);
					fListOfObjects->Add(hPtneg_TOF_Eta[i][so][j]);
					fListOfObjects->Add(hPtpos_TOF_Eta[i][so][j]);
					fListOfObjects->Add(hPneg_TOF_Eta[i][so][j]);
					fListOfObjects->Add(hPpos_TOF_Eta[i][so][j]);
				} // eta
			} // so
		} // topology

		for(Int_t j = 0; j < nHists; ++j){
			fListOfObjects->Add(hDeDxVsPSoInt[j]);
			fListOfObjects->Add(hnSigmaPiPosSoInt[j]);
			fListOfObjects->Add(hnSigmaKPosSoInt[j]);
			fListOfObjects->Add(hnSigmaPPosSoInt[j]);
			fListOfObjects->Add(hnSigmaPiNegSoInt[j]);
			fListOfObjects->Add(hnSigmaKNegSoInt[j]);
			fListOfObjects->Add(hnSigmaPNegSoInt[j]);
			fListOfObjects->Add(hBetavsPposSoInt[j]);
			fListOfObjects->Add(hBetavsPnegSoInt[j]);
			fListOfObjects->Add(hPtpos_TPC_EtaSoInt[j]);
			fListOfObjects->Add(hPtneg_TPC_EtaSoInt[j]);
			fListOfObjects->Add(hPtpos_TOF_EtaSoInt[j]);
			fListOfObjects->Add(hPtneg_TOF_EtaSoInt[j]);
			fListOfObjects->Add(hPpos_TOF_EtaSoInt[j]);
			fListOfObjects->Add(hPneg_TOF_EtaSoInt[j]);
			fListOfObjects->Add(hPtVsP[j]);
		}

		fListOfObjects->Add(hPtAllSoInt);
		fListOfObjects->Add(hPtpos_TPCSoInt);
		fListOfObjects->Add(hPtneg_TPCSoInt);
		fListOfObjects->Add(hPtpos_TOFSoInt);
		fListOfObjects->Add(hPtneg_TOFSoInt);

	} //	!fAnalysisMC


	///	else{

	for(Int_t cent=0; cent<nCent; cent++) {
		for(Int_t pid=0; pid<7; pid++){
			for(Int_t so=0; so<3; so++){
				hMcIn[cent][pid][so]=new TH1D(Form("hIn_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]), Form("MC in (pid %s)", Pid[pid]),nPtBins,ptBins);
				hMcInNeg[cent][pid][so]=new TH1D(Form("hInNeg_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC in (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcInPos[cent][pid][so]=new TH1D(Form("hInPos_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC in (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcOut[cent][pid][so]=new TH1D(Form("hMcOut_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcOutNeg[cent][pid][so]=new TH1D(Form("hMcOutNeg_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);
				hMcOutPos[cent][pid][so]=new TH1D(Form("hMcOutPos_%.2f-%.2f-%s-%s",CentMin[cent],CentMax[cent],Pid[pid],So[so]),Form("MC out (pid %s)",Pid[pid]),nPtBins,ptBins);

				/*if(cent<3){
				  fListOfObjects->Add(hMcIn[cent][pid][so]);
				  fListOfObjects->Add(hMcInNeg[cent][pid][so]);
				  fListOfObjects->Add(hMcInPos[cent][pid][so]);
				  fListOfObjects->Add(hMcOut[cent][pid][so]);
				  fListOfObjects->Add(hMcOutNeg[cent][pid][so]);
				  fListOfObjects->Add(hMcOutPos[cent][pid][so]);
				  }*/
			}
		}	// pid Eff

		///			fListOfObjects->Add(hPhi[cent]);

		//			hSOtVsSOm[cent] = new TH2D(Form("hSOtVsSOm-%.2f-%.2f",CentMin[cent],CentMax[cent]),";#it{S}_{O} generated;#it{S}_{O} reconstructed",10, 0, 1, 10, 0, 1);
		//			hSOtVsSOm[cent]->Sumw2();
		//			fListOfObjects->Add(hSOtVsSOm[cent]);

	}	// cent Eff

	hPtTruthVsPtRec = new TH2D("hPtTruthVsPtRec","Pt Truth Vs Pt Rec;#it{p}_{T}^{Gen};#it{p}_{T}^{Rec} ",200, 0, 200, 200, 0, 200);
	//	fListOfObjects->Add(hPtTruthVsPtRec);

	hPtTruthVsPtRecJetty = new TH2D("hPtTruthVsPtRecJetty","Pt Truth Vs Pt Rec Jetty Events;#it{p}_{T}^{Gen};#it{p}_{T}^{Rec} ",200, 0, 200, 200, 0, 200);
	//	fListOfObjects->Add(hPtTruthVsPtRecJetty);

	hPtTruthVsPtRecIsotr = new TH2D("hPtTruthVsPtRecIsotr","Pt Truth Vs Pt Rec Isotropic Events;#it{p}_{T}^{Gen};#it{p}_{T}^{Rec} ",200, 0, 200, 200, 0, 200);
	//	fListOfObjects->Add(hPtTruthVsPtRecIsotr);

	hTruthEtaSo = new TH1D("hTruthEtaSo","spherocity; #eta; counts",40,-1.0,1.0);
	//	fListOfObjects->Add(hTruthEtaSo);

	hTruthPhiSo = new TH1D("hTruthPhiSo","spherocity; #phi; counts",64,0.0,2*TMath::Pi());
	//	fListOfObjects->Add(hTruthPhiSo);

	hSOtvsTrks  = new TH2D("hSOtVsTrks","Truth SO vs Measured Ref. Mult.;Ref. mult. (|#eta|<0.8);#it{S}_{O} Truth ",100, 0, 100, 1000, 0, 1);
	//	fListOfObjects->Add(hSOtvsTrks);

	hSOtvsTrkst  = new TH2D("hSOtVsTrkst","Truth SO vs Truth. Mult.;Truth Mult. (|#eta|<0.8);#it{S}_{O} Truth ",100, 0, 100, 1000, 0, 1);
	//	fListOfObjects->Add(hSOtvsTrkst);

	hSOtvsV0M  = new TH2D("hSOtVsV0M","Truth SO vs V0M Per.;V0M Per.;#it{S}_{O} Truth",100, 0, 100, 1000, 0, 1);;
	//	fListOfObjects->Add(hSOtvsV0M);

	///}

	PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskSpherocity::UserExec(Option_t *)
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
	if ( !utils ){return;}

	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7;
	if(!isINT7selected)
		return;

	////fEvents->Fill(1.5,10);

	//--------------- Event Selection --------------------

	float MultPercentile = -1;
	float V0MPercentile = -1;
	float RefPercentile = -1;
	int fnRefGlobal = -1;
	//	int IndxTrksMult = -1;
	//	int IndxV0MMult = -1;

	fnRefGlobal = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8 );

	AliMultSelection *MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
	if(!MultSelection){ return; }
	if( MultSelection-> IsEventSelected() ){
		V0MPercentile = MultSelection->GetMultiplicityPercentile("V0M",false);
		RefPercentile = MultSelection->GetMultiplicityPercentile("RefMult08",false);

		//		IndxV0MMult  = GetCentralityClass(V0MPercentile);
		//		IndxTrksMult = GetCentralityClass(RefPercentile);

		if(fisV0Mestimator) MultPercentile = MultSelection->GetMultiplicityPercentile("V0M",false);
		else MultPercentile = MultSelection->GetMultiplicityPercentile("RefMult08",false);
	}
	else{ return; }

	/*	if(!fisV0Mestimator) {
		if(IndxTrksMult>=0)
		fCentClass = IndxTrksMult;
		else
		return;
		} else {
		if(IndxV0MMult>=0)
		fCentClass = IndxV0MMult;
		else
		return;
		}
		*/
	if( fESD->IsPileupFromSPDInMultBins() ){return;}
	fEvents->Fill(2.5,MultPercentile);
	///fEvents->Fill(2.5,10);

	if( fESD->IsIncompleteDAQ() ){return;}
	fEvents->Fill(3.5,MultPercentile);
	///fEvents->Fill(3.5,10);

	if( utils->IsSPDClusterVsTrackletBG(fESD) ){return;}
	fEvents->Fill(4.5,MultPercentile);
	///fEvents->Fill(4.5,10);

	if( !MultSelection->GetThisEventINELgtZERO() ){return;}
	fEvents->Fill(5.5,MultPercentile);
	///fEvents->Fill(5.5,10);

	if( !selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE) ){return;}
	fEvents->Fill(6.5,MultPercentile);
	///fEvents->Fill(6.5,10);

	if( !IsGoodZvertexPos(fESD) ){return;}
	fEvents->Fill(7.5,MultPercentile);
	////fEvents->Fill(7.5,10);

	float SOm = -1.0;
	SOm = GetSpherocity(hphiso, hetaso);

	double SOt = -1.0;
	//    if(fAnalysisMC)
	//        SOt = fSpheroUtils->GetEventShapeTrue(fMCStack,hTruthPhiSo,hTruthEtaSo);

	// Events with non-measured spherocity
	if( SOm < 0.0 ){
		fEvents->Fill(8.5,MultPercentile);
		return;
	}

	hSOrvsV0M->Fill(V0MPercentile,SOm);
	hSOrvsTrks->Fill(RefPercentile,SOm);
	hRefMultVsRefMultPer->Fill(RefPercentile,fnRefGlobal);

	fEvents->Fill(9.5,MultPercentile);

	//	fTrcksVsTrklets->Fill(fnRefGlobal,nRec);

	ProduceArrayTrksESD(SOm);

	if(fAnalysisMC){

		int TruthMult = -1;
		TruthMult = GetMultiplicityParticles(0.8);

		if(SOt>0){
			hSOtvsTrks->Fill(fnRefGlobal,SOt);
			hSOtvsV0M->Fill(V0MPercentile,SOt);
			hSOtvsTrkst->Fill(TruthMult,SOt);
			//			hSOtVsSOm[fCentClass]->Fill(SOt,SOm);
			//			hSOtVsSOm[10]->Fill(SOt,SOm);
		}
	}

	PostData(1, fListOfObjects);
}
//_____________________________________________________________________________
/*Int_t AliAnalysisTaskSpherocity::GetCentralityClass(Float_t percentile)
  {

  if((percentile<=0))return -1;

  Int_t Index = -1;
  if((percentile>70) && (percentile<=100))Index=9;
  else if((percentile>50) && (percentile<=70))Index=8;
  else if((percentile>40) && (percentile<=50))Index=7;
  else if((percentile>30) && (percentile<=40))Index=6;
  else if((percentile>20) && (percentile<=30))Index=5;
  else if((percentile>15) && (percentile<=20))Index=4;
  else if((percentile>10) && (percentile<=15))Index=3;
  else if((percentile>5) && (percentile<=10))Index=2;
  else if((percentile>1) && (percentile<=5))Index=1;
  else if((percentile>0) && (percentile<=1))Index=0;
  else Index=-1;

  return Index;
  }
  */
//_____________________________________________________________________________
void AliAnalysisTaskSpherocity::PtRecVsPtTruth( AliESDEvent *ESDevent, const Bool_t IsJetty )
{
	const Int_t nESDTracks = ESDevent->GetNumberOfTracks();
	for(Int_t iT = 0; iT < nESDTracks; iT++){

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);
		if(!esdTrack){continue;}

		UInt_t selectDebug = 0;
		if(fTrackFilter){
			selectDebug = fTrackFilter->IsSelected(esdTrack);
			if (!selectDebug) {continue;}
		}

		if(TMath::Abs(esdTrack->Eta())>fEtaCut)
			printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      ==========  %f\n",esdTrack->Eta());


		Float_t dcaxy = 0.;
		Float_t dcaz = 0.;
		esdTrack->GetImpactParameters(dcaxy,dcaz);
		if(TMath::Abs(dcaxy)>GetMaxDCApTDep(fcutDCAxy,esdTrack->Pt())){continue;}

		Short_t ncl = esdTrack->GetTPCNcls();
		if(ncl<50)
			printf("--------------------------------------------------         --------------- Ncl      ==========  %d\n",ncl);

		const Int_t label = TMath::Abs(esdTrack->GetLabel());
		TParticle* mcTrack = fMCStack->Particle(label);

		TParticlePDG* pdgPart = mcTrack->GetPDG();
		Double_t chargeMC = pdgPart->Charge();

		if (mcTrack){

			if((esdTrack->Charge()==0) || (TMath::Abs(chargeMC) < 0.1))
				continue;

			if((TMath::Abs(esdTrack->Eta())>fEtaCut) || (TMath::Abs(mcTrack->Eta())>fEtaCut))
				continue;

			if( fMCStack->IsPhysicalPrimary(label) ){
				hPtTruthVsPtRec->Fill(mcTrack->Pt(),esdTrack->Pt());

				if(IsJetty)
					hPtTruthVsPtRecJetty->Fill(mcTrack->Pt(),esdTrack->Pt());
				if(!IsJetty)
					hPtTruthVsPtRecIsotr->Fill(mcTrack->Pt(),esdTrack->Pt());
			}
		}
	}
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSpherocity::GetMultiplicityParticles(Double_t etaCut)
{
	// Fill the special MC histogram with the MC truth info

	Int_t trackmult = 0;
	const Int_t nTracksMC = fMCStack->GetNtrack();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		TParticle* trackMC = fMCStack->Particle(iTracks);
		if(!trackMC)
			continue;

		if( !(fMCStack->IsPhysicalPrimary(iTracks)) )
			continue;


		TParticlePDG* pdgPart = trackMC ->GetPDG();
		Double_t chargeMC = pdgPart->Charge();

		if( TMath::Abs(chargeMC) < 0.1 )
			continue;

		if ( TMath::Abs( trackMC->Eta() ) > etaCut )
			continue;

		trackmult++;

	}//MC track loop

	return trackmult;

}
//_____________________________________________________________________________
Short_t AliAnalysisTaskSpherocity::GetPidCode(Int_t pdgCode) const
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
void AliAnalysisTaskSpherocity::ProcessMCTruthESD(const Int_t Cent, const Int_t Spherocity)
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

		Int_t pdgCode = trackMC->GetPdgCode();
		Short_t pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		hMcIn[Cent][0][Spherocity]->Fill(trackMC->Pt());
		hMcIn[Cent][pidCodeMC][Spherocity]->Fill(trackMC->Pt());

		if( chargeMC < 0 ){
			hMcInNeg[Cent][0][Spherocity]->Fill(trackMC->Pt());
			hMcInNeg[Cent][pidCodeMC][Spherocity]->Fill(trackMC->Pt());
		}
		else{
			hMcInPos[Cent][0][Spherocity]->Fill(trackMC->Pt());
			hMcInPos[Cent][pidCodeMC][Spherocity]->Fill(trackMC->Pt());
		}

	}//MC track loop
}
//____________________________________________________________________
TParticle* AliAnalysisTaskSpherocity::FindPrimaryMother(AliStack* stack, Int_t label)
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
Int_t AliAnalysisTaskSpherocity::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
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
TParticle* AliAnalysisTaskSpherocity::FindPrimaryMotherV0(AliStack* stack, Int_t label)
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
Int_t AliAnalysisTaskSpherocity::FindPrimaryMotherLabelV0(AliStack* stack, Int_t label, Int_t& nSteps)
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
void AliAnalysisTaskSpherocity::ProduceArrayTrksESD( const float& Spherocity ){


	float MultPer  = -1.0;
	AliMultSelection *MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
	if(!MultSelection){ return; }
	if( MultSelection-> IsEventSelected() ){

		if(fisV0Mestimator) MultPer = MultSelection->GetMultiplicityPercentile("V0M",false);
		else MultPer = MultSelection->GetMultiplicityPercentile("RefMult08",false);

	}

	if((MultPer < 0.0) || (MultPer > 100.0)) return;

	int Topology = -1;

	if(( 0.0 < Spherocity )&&( Spherocity <= fJettyCutOff)){

		Topology = 0;

	}
	if(( Spherocity >= fIsotrCutOff) && ( Spherocity < 1.0 )){

		Topology = 1;
	}

	const Int_t nESDTracks = fESD->GetNumberOfTracks();
	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = fESD->GetTrack(iT);
		if( TMath::Abs(esdTrack->Eta()) > fEtaCut )
			continue;

		if(!fTrackFilterGolden->AcceptTrack(esdTrack))
			continue;

		double eta      = esdTrack->Eta();
		double phi      = esdTrack->Phi();
		double momentum = esdTrack->P();
		double pt       = esdTrack->Pt();
		float  dedx     = esdTrack->GetTPCsignal();
		float dcaxy = 0.0;
		float dcaz = 0.0;
		esdTrack->GetImpactParameters(dcaxy,dcaz);

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

		short ncl = esdTrack->GetTPCsignalN();
		if( ncl < fNcl )
			continue;

		if( TMath::Abs(dcaxy) > GetMaxDCApTDep(fcutDCAxy,pt) )
			continue;

		if( TOFPID(esdTrack) ){

			double trkLength = esdTrack->GetIntegratedLength();
			double beta = trkLength/((esdTrack->GetTOFsignal()-fPIDResponse->GetTOFResponse().GetStartTime(esdTrack->P()))*C_Value);

			if(esdTrack->Charge() < 0.0){

				if( Topology == 0 ){ // Topology = 0 -> Jetty topology

					if( Spherocity < fJettyCutOff_0 ){
						hBetavsPneg[Topology][0][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtneg_TOF_Eta[Topology][0][nh]->Fill(pt,MultPer);
						hPneg_TOF_Eta[Topology][0][nh]->Fill(esdTrack->P(),MultPer);
						hPtneg_TOF[Topology][0]->Fill(pt,MultPer);
					}
					if( Spherocity < fJettyCutOff_1 ){
						hBetavsPneg[Topology][1][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtneg_TOF_Eta[Topology][1][nh]->Fill(pt,MultPer);
						hPneg_TOF_Eta[Topology][1][nh]->Fill(esdTrack->P(),MultPer);
						hPtneg_TOF[Topology][1]->Fill(pt,MultPer);
					}
					if( Spherocity < fJettyCutOff_2 ){
						hBetavsPneg[Topology][2][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtneg_TOF_Eta[Topology][2][nh]->Fill(pt,MultPer);
						hPneg_TOF_Eta[Topology][2][nh]->Fill(esdTrack->P(),MultPer);
						hPtneg_TOF[Topology][2]->Fill(pt,MultPer);
					}

					hBetavsPneg[Topology][3][nh]->Fill(esdTrack->P(),beta,MultPer);
					hPtneg_TOF_Eta[Topology][3][nh]->Fill(pt,MultPer);
					hPneg_TOF_Eta[Topology][3][nh]->Fill(esdTrack->P(),MultPer);
					hPtneg_TOF[Topology][3]->Fill(pt,MultPer);
				}
				if( Topology == 1 ){ // Topology = 1 -> Isotropic topology 

					if( Spherocity > fIsotrCutOff_0 ){
						hBetavsPneg[Topology][0][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtneg_TOF_Eta[Topology][0][nh]->Fill(pt,MultPer);
						hPneg_TOF_Eta[Topology][0][nh]->Fill(esdTrack->P(),MultPer);
						hPtneg_TOF[Topology][0]->Fill(pt,MultPer);
					}
					if( Spherocity > fIsotrCutOff_1 ){
						hBetavsPneg[Topology][1][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtneg_TOF_Eta[Topology][1][nh]->Fill(pt,MultPer);
						hPneg_TOF_Eta[Topology][1][nh]->Fill(esdTrack->P(),MultPer);
						hPtneg_TOF[Topology][1]->Fill(pt,MultPer);
					}
					if( Spherocity > fIsotrCutOff_2 ){
						hBetavsPneg[Topology][2][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtneg_TOF_Eta[Topology][2][nh]->Fill(pt,MultPer);
						hPneg_TOF_Eta[Topology][2][nh]->Fill(esdTrack->P(),MultPer);
						hPtneg_TOF[Topology][2]->Fill(pt,MultPer);
					}

					hBetavsPneg[Topology][3][nh]->Fill(esdTrack->P(),beta,MultPer);
					hPtneg_TOF_Eta[Topology][3][nh]->Fill(pt,MultPer);
					hPneg_TOF_Eta[Topology][3][nh]->Fill(esdTrack->P(),MultPer);
					hPtneg_TOF[Topology][3]->Fill(pt,MultPer);
				}

				hBetavsPnegSoInt[nh]->Fill(esdTrack->P(),beta,MultPer);
				hPtneg_TOF_EtaSoInt[nh]->Fill(pt,MultPer);
				hPneg_TOF_EtaSoInt[nh]->Fill(esdTrack->P(),MultPer);
				hPtneg_TOFSoInt->Fill(pt,MultPer);

			}
			else{ // Positive charge

				if( Topology == 0 ){ // Topology = 0 -> Jetty topology 

					if( Spherocity < fJettyCutOff_0 ){
						hBetavsPpos[Topology][0][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtpos_TOF_Eta[Topology][0][nh]->Fill(pt,MultPer);
						hPpos_TOF_Eta[Topology][0][nh]->Fill(esdTrack->P(),MultPer);
						hPtpos_TOF[Topology][0]->Fill(pt,MultPer);
					}
					if( Spherocity < fJettyCutOff_1 ){
						hBetavsPpos[Topology][1][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtpos_TOF_Eta[Topology][1][nh]->Fill(pt,MultPer);
						hPpos_TOF_Eta[Topology][1][nh]->Fill(esdTrack->P(),MultPer);
						hPtpos_TOF[Topology][1]->Fill(pt,MultPer);
					}
					if( Spherocity < fJettyCutOff_2 ){
						hBetavsPpos[Topology][2][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtpos_TOF_Eta[Topology][2][nh]->Fill(pt,MultPer);
						hPpos_TOF_Eta[Topology][2][nh]->Fill(esdTrack->P(),MultPer);
						hPtpos_TOF[Topology][2]->Fill(pt,MultPer);
					}

					hBetavsPpos[Topology][3][nh]->Fill(esdTrack->P(),beta,MultPer);
					hPtpos_TOF_Eta[Topology][3][nh]->Fill(pt,MultPer);
					hPpos_TOF_Eta[Topology][3][nh]->Fill(esdTrack->P(),MultPer);
					hPtpos_TOF[Topology][3]->Fill(pt,MultPer);
				}
				if( Topology == 1 ){ // Topology = 1 -> Isotropic topology 

					if( Spherocity > fIsotrCutOff_0 ){
						hBetavsPpos[Topology][0][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtpos_TOF_Eta[Topology][0][nh]->Fill(pt,MultPer);
						hPpos_TOF_Eta[Topology][0][nh]->Fill(esdTrack->P(),MultPer);
						hPtpos_TOF[Topology][0]->Fill(pt,MultPer);
					}
					if( Spherocity > fIsotrCutOff_1 ){
						hBetavsPpos[Topology][1][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtpos_TOF_Eta[Topology][1][nh]->Fill(pt,MultPer);
						hPpos_TOF_Eta[Topology][1][nh]->Fill(esdTrack->P(),MultPer);
						hPtpos_TOF[Topology][1]->Fill(pt,MultPer);
					}
					if( Spherocity > fIsotrCutOff_2 ){
						hBetavsPpos[Topology][2][nh]->Fill(esdTrack->P(),beta,MultPer);
						hPtpos_TOF_Eta[Topology][2][nh]->Fill(pt,MultPer);
						hPpos_TOF_Eta[Topology][2][nh]->Fill(esdTrack->P(),MultPer);
						hPtpos_TOF[Topology][2]->Fill(pt,MultPer);
					}

					hBetavsPpos[Topology][3][nh]->Fill(esdTrack->P(),beta,MultPer);
					hPtpos_TOF_Eta[Topology][3][nh]->Fill(pt,MultPer);
					hPpos_TOF_Eta[Topology][3][nh]->Fill(esdTrack->P(),MultPer);
					hPtpos_TOF[Topology][3]->Fill(pt,MultPer);
				}

				hBetavsPposSoInt[nh]->Fill(esdTrack->P(),beta,MultPer);
				hPtpos_TOF_EtaSoInt[nh]->Fill(pt,MultPer);
				hPpos_TOF_EtaSoInt[nh]->Fill(esdTrack->P(),MultPer);
				hPtpos_TOFSoInt->Fill(pt,MultPer);
			}
		}

		hPtVsP[nh]->Fill(momentum,pt);

		if(!PhiCut(esdTrack->Pt(), phi, esdTrack->Charge(), MAGF, fcutLow, fcutHigh))
			continue;

		if(fdEdxCalibrated) dedx *= 50/EtaCalibration(eta);

		if( Topology == 0 ){ // Topology = 0 -> Jetty topology

			if( Spherocity < fJettyCutOff_0 ){
				hPtAll[Topology][0]->Fill(pt,MultPer);
				hDeDxVsP[Topology][0][nh]->Fill(momentum,dedx,MultPer);
			}
			if( Spherocity < fJettyCutOff_1 ){
				hPtAll[Topology][1]->Fill(pt,MultPer);
				hDeDxVsP[Topology][1][nh]->Fill(momentum,dedx,MultPer);
			}
			if( Spherocity < fJettyCutOff_2 ){
				hPtAll[Topology][2]->Fill(pt,MultPer);
				hDeDxVsP[Topology][2][nh]->Fill(momentum,dedx,MultPer);
			}

			hPtAll[Topology][3]->Fill(pt,MultPer);
			hDeDxVsP[Topology][3][nh]->Fill(momentum,dedx,MultPer);
		}
		if( Topology == 1 ){ // Topology = 1 -> Isotropic topology 

			if( Spherocity > fIsotrCutOff_0 ){
				hPtAll[Topology][0]->Fill(pt,MultPer);
				hDeDxVsP[Topology][0][nh]->Fill(momentum,dedx,MultPer);
			}
			if( Spherocity > fIsotrCutOff_1 ){
				hPtAll[Topology][1]->Fill(pt,MultPer);
				hDeDxVsP[Topology][1][nh]->Fill(momentum,dedx,MultPer);
			}
			if( Spherocity > fIsotrCutOff_2 ){
				hPtAll[Topology][2]->Fill(pt,MultPer);
				hDeDxVsP[Topology][2][nh]->Fill(momentum,dedx,MultPer);
			}

			hPtAll[Topology][3]->Fill(pt,MultPer);
			hDeDxVsP[Topology][3][nh]->Fill(momentum,dedx,MultPer);
		}

		hPtAllSoInt->Fill(pt,MultPer);
		hDeDxVsPSoInt[nh]->Fill(momentum,dedx,MultPer);

		if(esdTrack->Charge() < 0.0){

			if( Topology == 0 ){ // Topology = 0 -> Jetty topology

				if( Spherocity < fJettyCutOff_0 ){
					hnSigmaPiNeg[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKNeg[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPNeg[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtneg_TPC_Eta[Topology][0][nh]->Fill(pt,MultPer);
					hPtneg_TPC[Topology][0]->Fill(pt,MultPer);
				}
				if( Spherocity < fJettyCutOff_1 ){
					hnSigmaPiNeg[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKNeg[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPNeg[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtneg_TPC_Eta[Topology][1][nh]->Fill(pt,MultPer);
					hPtneg_TPC[Topology][1]->Fill(pt,MultPer);
				}
				if( Spherocity < fJettyCutOff_2 ){
					hnSigmaPiNeg[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKNeg[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPNeg[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtneg_TPC_Eta[Topology][2][nh]->Fill(pt,MultPer);
					hPtneg_TPC[Topology][2]->Fill(pt,MultPer);
				}

				hnSigmaPiNeg[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
				hnSigmaKNeg[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
				hnSigmaPNeg[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
				hPtneg_TPC_Eta[Topology][3][nh]->Fill(pt,MultPer);
				hPtneg_TPC[Topology][3]->Fill(pt,MultPer);
			}
			if( Topology == 1 ){ // Topology = 1 -> Isotropic topology 

				if( Spherocity > fIsotrCutOff_0 ){
					hnSigmaPiNeg[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKNeg[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPNeg[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtneg_TPC_Eta[Topology][0][nh]->Fill(pt,MultPer);
					hPtneg_TPC[Topology][0]->Fill(pt,MultPer);
				}
				if( Spherocity > fIsotrCutOff_1 ){
					hnSigmaPiNeg[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKNeg[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPNeg[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtneg_TPC_Eta[Topology][1][nh]->Fill(pt,MultPer);
					hPtneg_TPC[Topology][1]->Fill(pt,MultPer);
				}
				if( Spherocity > fIsotrCutOff_2 ){
					hnSigmaPiNeg[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKNeg[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPNeg[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtneg_TPC_Eta[Topology][2][nh]->Fill(pt,MultPer);
					hPtneg_TPC[Topology][2]->Fill(pt,MultPer);
				}

				hnSigmaPiNeg[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
				hnSigmaKNeg[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
				hnSigmaPNeg[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
				hPtneg_TPC_Eta[Topology][3][nh]->Fill(pt,MultPer);
				hPtneg_TPC[Topology][3]->Fill(pt,MultPer);
			}

			hnSigmaPiNegSoInt[nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
			hnSigmaKNegSoInt[nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
			hnSigmaPNegSoInt[nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
			hPtneg_TPC_EtaSoInt[nh]->Fill(pt,MultPer);
			hPtneg_TPCSoInt->Fill(pt,MultPer);

		}else{

			if( Topology == 0 ){ // Topology = 0 -> Jetty topology

				if( Spherocity < fJettyCutOff_0 ){
					hnSigmaPiPos[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKPos[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPPos[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtpos_TPC_Eta[Topology][0][nh]->Fill(pt,MultPer);
					hPtpos_TPC[Topology][0]->Fill(pt,MultPer);
				}
				if( Spherocity < fJettyCutOff_1 ){
					hnSigmaPiPos[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKPos[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPPos[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtpos_TPC_Eta[Topology][1][nh]->Fill(pt,MultPer);
					hPtpos_TPC[Topology][1]->Fill(pt,MultPer);
				}
				if( Spherocity < fJettyCutOff_2 ){
					hnSigmaPiPos[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKPos[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPPos[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtpos_TPC_Eta[Topology][2][nh]->Fill(pt,MultPer);
					hPtpos_TPC[Topology][2]->Fill(pt,MultPer);
				}

				hnSigmaPiPos[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
				hnSigmaKPos[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
				hnSigmaPPos[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
				hPtpos_TPC_Eta[Topology][3][nh]->Fill(pt,MultPer);
				hPtpos_TPC[Topology][3]->Fill(pt,MultPer);

			}
			if( Topology == 1 ){ // Topology = 1 -> Isotropic topology 

				if( Spherocity > fIsotrCutOff_0 ){
					hnSigmaPiPos[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKPos[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPPos[Topology][0][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtpos_TPC_Eta[Topology][0][nh]->Fill(pt,MultPer);
					hPtpos_TPC[Topology][0]->Fill(pt,MultPer);
				}
				if( Spherocity > fIsotrCutOff_1 ){
					hnSigmaPiPos[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKPos[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPPos[Topology][1][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtpos_TPC_Eta[Topology][1][nh]->Fill(pt,MultPer);
					hPtpos_TPC[Topology][1]->Fill(pt,MultPer);
				}
				if( Spherocity > fIsotrCutOff_2 ){
					hnSigmaPiPos[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
					hnSigmaKPos[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
					hnSigmaPPos[Topology][2][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
					hPtpos_TPC_Eta[Topology][2][nh]->Fill(pt,MultPer);
					hPtpos_TPC[Topology][2]->Fill(pt,MultPer);
				}

				hnSigmaPiPos[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
				hnSigmaKPos[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
				hnSigmaPPos[Topology][3][nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
				hPtpos_TPC_Eta[Topology][3][nh]->Fill(pt,MultPer);
				hPtpos_TPC[Topology][3]->Fill(pt,MultPer);
			}

			hnSigmaPiPosSoInt[nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion),MultPer);
			hnSigmaKPosSoInt[nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon),MultPer);
			hnSigmaPPosSoInt[nh]->Fill(pt,fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton),MultPer);
			hPtpos_TPC_EtaSoInt[nh]->Fill(pt,MultPer);
			hPtpos_TPCSoInt->Fill(pt,MultPer);
		} // negative charge

	}//end of track loop


}

//----------------------------------------------------------------------------------
/*
   void AliAnalysisTaskSpherocity::ProduceArrayV0ESD( AliESDEvent *ESDevent, const Int_t Cent, const Int_t Spherocity ){

   Int_t nv0s = ESDevent->GetNumberOfV0s();

   fcentAfterV0s->Fill(Cent);
   fcentAfterV0s->Fill(10);

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

				       if(track->GetTPCsignalN()<fNcl)continue;
				       Double_t phi     = track->Phi();

				       if(!PhiCut(track->Pt(), phi, track->Charge(), MAGF, fcutLow, fcutHigh))
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
							       hMIPVsEtaV0s[Cent][Spherocity]->Fill(eta,dedx);
							       pMIPVsEtaV0s[Cent][Spherocity]->Fill(eta,dedx);
							       //hMIPVsEtaV0s[10][Spherocity]->Fill(eta,dedx);
							       //pMIPVsEtaV0s[10][Spherocity]->Fill(eta,dedx);
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
					       histPiV0[Cent][nh][Spherocity]->Fill(momentum, dedx);
					       //histPiV0[10][nh]->Fill(momentum, dedx);
				       }
				       else{
					       histPV0[Cent][nh][Spherocity]->Fill(momentum, dedx);
					       //histPV0[10][nh]->Fill(momentum, dedx);
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

			       if(!PhiCut(track->Pt(), phi, track->Charge(), MAGF, fcutLow, fcutHigh))
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
			       histEV0[Cent][nh][Spherocity]->Fill(momentum, dedx);
			       //histEV0[10][nh]->Fill(momentum, dedx);
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
*/

//________________________________________________________________________
Bool_t AliAnalysisTaskSpherocity::selectVertex2015pp(AliESDEvent *esd,
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
Bool_t AliAnalysisTaskSpherocity::IsGoodSPDvertexRes(const AliESDVertex* spdVertex)
{

	if( !spdVertex ) return kFALSE;
	if( spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25) ) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSpherocity::IsGoodZvertexPos(AliESDEvent *esd)
{

	if( !esd ) return kFALSE;
	//Cut on the vertex z position
	const AliESDVertex * vertex = esd->GetPrimaryVertex();
	if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
	return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskSpherocity::SetTrackCutsSpherocity(AliAnalysisFilter* fTrackFilter){

	AliESDtrackCuts* esdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	esdTrackCuts->SetRequireTPCRefit(kTRUE);
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	esdTrackCuts->SetEtaRange(-0.8,0.8);
	fTrackFilter->AddCuts(esdTrackCuts);
}
//________________________________________________________________________
Bool_t AliAnalysisTaskSpherocity::PhiCut(Double_t pt, Double_t phi, Double_t q, Float_t   mag, TF1* phiCutLow, TF1* phiCutHigh)
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

	//	hPhi[fCentClass]->Fill(pt, phi);

	return kTRUE;
}
//________________________________________________________________________
Float_t AliAnalysisTaskSpherocity::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI){

	Double_t maxDCAxy = 10;
	maxDCAxy = fMaxDCAxy->Eval(ptI);
	return maxDCAxy;

}
//________________________________________________________________________
float AliAnalysisTaskSpherocity::GetSpherocity( TH1D * hphi, TH1D *heta )
{
	vector<Float_t> pt;
	vector<Float_t> eta;
	vector<Float_t> phi;

	if(fESD)
		fNrec = ReadESDEvent( pt, eta, phi, hphi, heta );

	if( fNrec < fMinMult )
		return -0.5;

	float spherocity = AnalyseGetSpherocity( pt, eta, phi );

	return spherocity;

}
//_____________________________________________________________________
int AliAnalysisTaskSpherocity::ReadESDEvent( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TH1D * hphi, TH1D *heta ){

	ptArray.clear();
	etaArray.clear();
	phiArray.clear();

	int nTracks = fESD->GetNumberOfTracks();
	int nRec    = 0;

	for(int iT = 0; iT < nTracks; iT++) {
		AliESDtrack* Track = 0;
		Track = fESD->GetTrack(iT);
		if(!Track)
			continue;

		float eta  = Track->Eta();
		float pt   = Track->Pt();
		float phi  = Track->Phi();

		if( !(TMath::Abs(eta) < fEtaCut) )
			continue;

		//cuts in pt
		if( pt > 1E8 || pt < 0.15 )
			continue;

		//quality cuts
		if(!fTrackFilter->IsSelected(Track))
			continue;

		if(hphi)
			hphi->Fill(phi);
		if(heta)
			heta->Fill(eta);

		ptArray.push_back(pt);
		etaArray.push_back(eta);
		phiArray.push_back(phi);

		nRec++;

	} //close first loop on nTracks

	return nRec;

}
//________________________________________________________________________
float AliAnalysisTaskSpherocity::AnalyseGetSpherocity( const vector<Float_t> &pt, const vector<Float_t> &eta, const vector<Float_t> &phi ){


	float spherocity = -10.0;
	float pFull = 0;
	float Spherocity = 2;

	//computing total pt
	float sumapt = 0;
	for(int i1 = 0; i1 < fNrec; ++i1){
		//    sumapt += pt[i1];       
		sumapt++;       
	}

	//Getting thrust
	for(int i = 0; i < 360/(fSizeStep); ++i){
		float numerador = 0;
		float phiparam  = 0;
		float nx = 0;
		float ny = 0;
		phiparam=( (TMath::Pi()) * i * fSizeStep ) / 180; // parametrization of the angle
		nx = TMath::Cos(phiparam);            // x component of an unitary vector n
		ny = TMath::Sin(phiparam);            // y component of an unitary vector n
		for(int i1 = 0; i1 < fNrec; ++i1){

			//            float pxA = pt[i1] * TMath::Cos( phi[i1] );
			//            float pyA = pt[i1] * TMath::Sin( phi[i1] );
			float pxA = 1.0 * TMath::Cos( phi[i1] );
			float pyA = 1.0 * TMath::Sin( phi[i1] );

			numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
		}
		pFull=TMath::Power( (numerador / sumapt),2 );
		if(pFull < Spherocity)//maximization of pFull
		{
			Spherocity = pFull;
		}
	}

	spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;


	return spherocity;

}
//________________________________________________________________________
bool AliAnalysisTaskSpherocity::TOFPID(AliESDtrack * track)
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
//________________________________________________________________________
double AliAnalysisTaskSpherocity::EtaCalibration(const double &eta){

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
	}else{
		aPos = 49.6975; bPos = 2.32535; cPos = -42.6516; dPos = 283.058; ePos = -1009.58; fPos = 1945.89; gPos = -1871.23; hPos = 698.552;
		aNeg = 49.8071; bNeg = 9.78466; cNeg = 120.018;  dNeg = 603.325; eNeg = 1470.92;  fNeg = 1819.63; gNeg = 1073.82;  hNeg = 230.142;
	}

	for(int i=0; i<8; ++i){
		fEtaCalibrationPos->SetParameter(i,0);
		fEtaCalibrationNeg->SetParameter(i,0);
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
/*Double_t AliAnalysisTaskSpherocity::EtaCalibrationNeg( const Int_t Cent, const Double_t eta){


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
Double_t AliAnalysisTaskSpherocity::EtaCalibrationPos( const Int_t Cent, const Double_t eta){


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
Double_t AliAnalysisTaskSpherocity::EtaCalibrationNegEl(const Int_t Cent, const Double_t eta){


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
Double_t AliAnalysisTaskSpherocity::EtaCalibrationPosEl(const Int_t Cent, const Double_t eta){


for(Int_t i=0; i<5; ++i)
felededxfitPos->SetParameter(i,0);

felededxfitPos->SetParameter(0,aPosEl[Cent]);
felededxfitPos->SetParameter(1,bPosEl[Cent]);
felededxfitPos->SetParameter(2,cPosEl[Cent]);
felededxfitPos->SetParameter(3,dPosEl[Cent]);
felededxfitPos->SetParameter(4,ePosEl[Cent]);

return felededxfitPos->Eval(eta);

}
*/
