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


static const int nPid              = 4;
static const int nRegion           = 4;
static const int nHists            = 4;
static const int nRt               = 5;
const char* Region[4] = {"Toward","Away","Transverse","FullAzimuth"};
const char* Pid[nPid] = {"Charged","Pion","Kaon","Proton"};
const char* Charge[2] = {"Pos","Neg"};

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
		utils(0x0),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fIsMCclosure(kTRUE),
		fRandom(0x0),
		fNcl(70),
		fEtaCut(0.9),
		fDeDxMIPMin(40),
		fDeDxMIPMax(60),
		fdEdxHigh(200),
		fdEdxLow(40),
		fPeriod("l"),
		fMeanChT(7.11),
		fMeanMultTSMCGen(7.412),
		fMeanMultTSMCRec(7.412),
		fSetTPConlyTrkCuts(kFALSE),
		fLeadPtCutMin(5.0),
		fLeadPtCutMax(40.0),
		fGenLeadPhi(0.0),
		fGenLeadPt(0.0),
		fGenLeadIn(0.0),
		fRecLeadPhi(0.0),
		fRecLeadPt(0.0),
		fRecLeadIn(0.0),
		fPtMin(0.3),
		fListOfObjects(0),
		fEvents(0x0),
		hNchRecVsPtRecOut(0x0),
		hNchRecVsPtGenOut(0x0),
		hPtPriGen(0x0),
		hPtRec(0x0),
		hPtPriRec(0x0),
		hPtSecRec(0x0),
		fPtLVsNchGen(0x0),
		fPtLVsNchRec(0x0),
		hMultTSGen(0x0),
		hMultTSRec(0x0),
		hNchTSGen(0x0),
		hNchTSRec(0x0),
		hNchTSGen_1(0x0),
		hNchTSContamination(0x0),
		hNchTSRecAll(0x0),
		hNchResponse(0x0),
		hNchRMvsPt(0x0),
		hNchTSGenTest(0x0),
		hNchTSRecTest(0x0),
		hNchTSData(0x0),
		hRTData(0x0),
		hPtLVsRT(0x0),
		fEtaCalibration(0x0),
		fEtaCalibrationEl(0x0),
		fcutDCAxy(0x0),
		fcutLow(0x0),
		fcutHigh(0x0)

{

	for(int pid = 0; pid < nPid; ++pid){
		hNchGenVsPtGenIn[pid] = 0;
		hNchRecVsPtGenIn[pid] = 0;
		hNchGenVsPtGenPosIn[pid] = 0;
		hNchGenVsPtGenNegIn[pid] = 0;
		hNchGenVsPtRecIn[pid] = 0;
		hNchGenVsPtRecPosIn[pid] = 0;
		hNchGenVsPtRecNegIn[pid] = 0;
		hNchGenVsPtRecInTOF[pid] = 0;
		hNchGenVsPtRecPosInTOF[pid] = 0;
		hNchGenVsPtRecNegInTOF[pid] = 0;
		hPtResponsePID[pid] = 0;
	}

	for(int r = 0; r < nRegion; r++){

		hPhiGen[r] = 0;
		hPhiRec[r] = 0;
		hPhiData[r] = 0;
		hNchGenGTZVsPtGen[r] = 0;
		hNchGenGTZVsPtRec[r] = 0;

		for(int pid = 0; pid < nPid; ++pid){

			hNchGenVsPtGenPID[r][pid] = 0;
			hNchGenVsPtRec[r][pid] = 0;

			hNchVsPtDataTPC[r][pid] = 0;
			hNchVsPtDataPosTPC[r][pid] = 0;
			hNchVsPtDataNegTPC[r][pid] = 0;

			hNchVsPtDataTOF[r][pid] = 0;
			hNchVsPtDataPosTOF[r][pid] = 0;
			hNchVsPtDataNegTOF[r][pid] = 0;
		}
	}        

	for(int j = 0; j < nHists; ++j){
		hPtVsP[j]=0;
		///hBetavsPMB[j]=0;
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
	utils(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fIsMCclosure(kTRUE),
	fRandom(0x0),
	fNcl(70),
	fEtaCut(0.9),
	fDeDxMIPMin(40),
	fDeDxMIPMax(60),
	fdEdxHigh(200),
	fdEdxLow(40),
	fPeriod("l"),
	fMeanChT(7.11),
	fMeanMultTSMCGen(7.412),
	fMeanMultTSMCRec(7.412),
	fSetTPConlyTrkCuts(kFALSE),
	fLeadPtCutMin(5.0),
	fLeadPtCutMax(40.0),
	fGenLeadPhi(0.0),
	fGenLeadPt(0.0),
	fGenLeadIn(0.0),
	fRecLeadPhi(0.0),
	fRecLeadPt(0.0),
	fRecLeadIn(0.0),
	fPtMin(0.3),
	fListOfObjects(0),
	fEvents(0x0),
	hNchRecVsPtRecOut(0x0),
	hNchRecVsPtGenOut(0x0),
	hPtPriGen(0x0),
	hPtRec(0x0),
	hPtPriRec(0x0),
	hPtSecRec(0x0),
	fPtLVsNchGen(0x0),
	fPtLVsNchRec(0x0),
	hMultTSGen(0x0),
	hMultTSRec(0x0),
	hNchTSGen(0x0),
	hNchTSRec(0x0),
	hNchTSGen_1(0x0),
	hNchTSContamination(0x0),
	hNchTSRecAll(0x0),
	hNchResponse(0x0),
	hNchRMvsPt(0x0),
	hNchTSGenTest(0x0),
	hNchTSRecTest(0x0),
	hNchTSData(0x0),
	hRTData(0x0),
	hPtLVsRT(0x0),
	fEtaCalibration(0x0),
	fEtaCalibrationEl(0x0),
	fcutDCAxy(0x0),
	fcutLow(0x0),
	fcutHigh(0x0)

{

	for(int pid = 0; pid < nPid; ++pid){
		hNchGenVsPtGenIn[pid] = 0;
		hNchRecVsPtGenIn[pid] = 0;
		hNchGenVsPtGenPosIn[pid] = 0;
		hNchGenVsPtGenNegIn[pid] = 0;
		hNchGenVsPtRecIn[pid] = 0;
		hNchGenVsPtRecPosIn[pid] = 0;
		hNchGenVsPtRecNegIn[pid] = 0;
		hNchGenVsPtRecInTOF[pid] = 0;
		hNchGenVsPtRecPosInTOF[pid] = 0;
		hNchGenVsPtRecNegInTOF[pid] = 0;
		hPtResponsePID[pid] = 0;
	}

	for(int r = 0; r < nRegion; r++){

		hPhiGen[r] = 0;
		hPhiRec[r] = 0;
		hPhiData[r] = 0;
		hNchGenGTZVsPtGen[r] = 0;
		hNchGenGTZVsPtRec[r] = 0;

		for(int pid = 0; pid < 4; ++pid){

			hNchGenVsPtGenPID[r][pid] = 0;
			hNchGenVsPtRec[r][pid] = 0;

			hNchVsPtDataTPC[r][pid] = 0;
			hNchVsPtDataPosTPC[r][pid] = 0;
			hNchVsPtDataNegTPC[r][pid] = 0;

			hNchVsPtDataTOF[r][pid] = 0;
			hNchVsPtDataPosTOF[r][pid] = 0;
			hNchVsPtDataNegTOF[r][pid] = 0;
		}

	}        
	for(int j = 0; j < nHists; ++j){
		hPtVsP[j]=0;
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

	// Definition of trackcuts
	if(!fTrackFilter){
		fTrackFilter = new AliAnalysisFilter("trackFilter");
		SetTrackCuts(fTrackFilter);
	}
	
	printf("The min cut in Pt is: %f\n",fPtMin);


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

	const int nPtBins = 54;
	double ptBins[nPtBins+1] = {
		0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35,
		0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
		0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
		1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40,
		3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00,
		9.00, 10.00, 12.00, 15.00, 20.00};

	const int nBinsNsigma = 40;
	double binsNsigma[nBinsNsigma+1] = {0};

	for(int i = 0; i <= nBinsNsigma; ++i){
		binsNsigma[i] = -10.0+i*0.5;

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

	printf("<Nch>_{gen} = %f  <Nch>_{rec} = %f  <Nch>_{dat} = %f\n",fMeanMultTSMCGen,fMeanMultTSMCRec,fMeanChT);

	const int nBinsRT = 100;
	double binsRT[nBinsRT+1] = {0};

	for(int i = 0; i <= nBinsRT; ++i)
		binsRT[i] = (((double)i)-0.5)/1.0;

	const char* ending[nHists] = {"02", "24", "46", "68"};

	fcutDCAxy = new TF1("fMaxDCAxy","[0]+[1]/(x^[2])",0,1e10);
	fcutDCAxy->SetParameter(0,0.0105);
	fcutDCAxy->SetParameter(1,0.0350);
	fcutDCAxy->SetParameter(2,1.1);

	fcutLow = new TF1("StandardPhiCutLow",  "0.1/x/x+TMath::Pi()/18.0-0.025", 0, 50);
	fcutHigh = new TF1("StandardPhiCutHigh", "0.12/x+TMath::Pi()/18.0+0.035", 0, 50);

	fEtaCalibration   = new TF1("fDeDxVsEtaPos", "pol7", 0.0, 1.0);
	fEtaCalibrationEl = new TF1("fDeDxVsEtaEl", "pol4", 0.0, 1.0);

	hNchTSData = new TH1D("hMultTSData",";#it{N}_{acc} mult; Entries",100,-0.5,99.5);
	hNchTSData->Sumw2();
	fListOfObjects->Add(hNchTSData);

	hRTData = new TH1D("hRTData",";#it{R}_{T}; Entries",nBinsRT,binsRT);
	hRTData->Sumw2();
	fListOfObjects->Add(hRTData);

	hPtLVsRT = new TH2D("hPtLVsRTData", "; #it{p}^{L}_{rec} (GeV/#it{c}); #it{R}_{T}",nPtBins,ptBins,nBinsRT,binsRT);
	hPtLVsRT->Sumw2();
	fListOfObjects->Add(hPtLVsRT);

	for(int r = 0; r < nRegion; ++r){

		hPhiData[r] = new TH1D(Form("hPhiData_%s",Region[r]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
		hPhiData[r]->Sumw2();

		for(int pid = 0; pid < nPid; ++pid){

			hNchVsPtDataTPC[r][pid] = new TH3D(Form("hNchVsPtData_TPC_%s_%s",Region[r],Pid[pid]),";#it{R}_{T};#it{p}_{T}^{rec};n#sigma",nBinsRT,binsRT,nPtBins,ptBins,nBinsNsigma,binsNsigma);
			hNchVsPtDataTPC[r][pid]->Sumw2();
			fListOfObjects->Add(hNchVsPtDataTPC[r][pid]);

			hNchVsPtDataPosTPC[r][pid] = new TH3D(Form("hNchVsPtData_Pos_TPC_%s_%s",Region[r],Pid[pid]),";#it{R}_{T};#it{p}_{T}^{rec};n#sigma",nBinsRT,binsRT,nPtBins,ptBins,nBinsNsigma,binsNsigma);
			hNchVsPtDataPosTPC[r][pid]->Sumw2();
			fListOfObjects->Add(hNchVsPtDataPosTPC[r][pid]);

			hNchVsPtDataNegTPC[r][pid] = new TH3D(Form("hNchVsPtData_Neg_TPC_%s_%s",Region[r],Pid[pid]),";#it{R}_{T};#it{p}_{T}^{rec};n#sigma",nBinsRT,binsRT,nPtBins,ptBins,nBinsNsigma,binsNsigma);
			hNchVsPtDataNegTPC[r][pid]->Sumw2();
			fListOfObjects->Add(hNchVsPtDataNegTPC[r][pid]);

			hNchVsPtDataTOF[r][pid] = new TH3D(Form("hNchVsPtData_TOF_%s_%s",Region[r],Pid[pid]),";#it{R}_{T};#it{p}_{T}^{rec};n#sigma",nBinsRT,binsRT,nPtBins,ptBins,nBinsNsigma,binsNsigma);
			hNchVsPtDataTOF[r][pid]->Sumw2();
			fListOfObjects->Add(hNchVsPtDataTOF[r][pid]);

			hNchVsPtDataPosTOF[r][pid] = new TH3D(Form("hNchVsPtData_Pos_TOF_%s_%s",Region[r],Pid[pid]),";#it{R}_{T};#it{p}_{T}^{rec};n#sigma",nBinsRT,binsRT,nPtBins,ptBins,nBinsNsigma,binsNsigma);
			hNchVsPtDataPosTOF[r][pid]->Sumw2();
			fListOfObjects->Add(hNchVsPtDataPosTOF[r][pid]);

			hNchVsPtDataNegTOF[r][pid] = new TH3D(Form("hNchVsPtData_Neg_TOF_%s_%s",Region[r],Pid[pid]),";#it{R}_{T};#it{p}_{T}^{rec};n#sigma",nBinsRT,binsRT,nPtBins,ptBins,nBinsNsigma,binsNsigma);
			hNchVsPtDataNegTOF[r][pid]->Sumw2();
			fListOfObjects->Add(hNchVsPtDataNegTOF[r][pid]);

		}

		fListOfObjects->Add(hPhiData[r]);

	}

	for( int j = 0; j < nHists; j++ ){

		hPtVsP[j] = new TH2D(Form("hPtVsP_%s",ending[j]),";#it{p} [GeV/c]; #it{p}_{T}",nPtBins,ptBins,nPtBins,ptBins);
		hPtVsP[j]->Sumw2();
		////		fListOfObjects->Add(hPtVsP[j]);

	}

	if(fAnalysisMC){

		hMultTSRec = new TH1D("hMultTSRec",";#it{N}_{acc} mult; Entries",100,-0.5,99.5);
		hMultTSRec->Sumw2();
		fListOfObjects->Add(hMultTSRec);

		hMultTSGen = new TH1D("hMultTSGen",";#it{N}_{acc} mult; Entries",100,-0.5,99.5);
		hMultTSGen->Sumw2();
		fListOfObjects->Add(hMultTSGen);

		fPtLVsNchGen = new TH2D("fPtLVsNchGen", "; #it{p}^{L}_{gen} (GeV/#it{c}); #it{R}_{T}",nPtBins,ptBins,nBinsRT,binsRT);
		fPtLVsNchGen->Sumw2();
		fListOfObjects->Add(fPtLVsNchGen);

		fPtLVsNchRec = new TH2D("fPtLVsNchRec", "; #it{p}^{L}_{rec} (GeV/#it{c}); #it{R}_{T}",nPtBins,ptBins,nBinsRT,binsRT);
		fPtLVsNchRec->Sumw2();
		fListOfObjects->Add(fPtLVsNchRec);

		hNchTSGen = new TH1D("hNchTSGen","; #it{R}_{T}; Entries",nBinsRT,binsRT);
		hNchTSGen->Sumw2();
		fListOfObjects->Add(hNchTSGen);

		hNchTSRec = new TH1D("hNchTSRec","; #it{R}_{T}; Entries",nBinsRT,binsRT);
		hNchTSRec->Sumw2();
		fListOfObjects->Add(hNchTSRec);

		hNchTSGen_1 = new TH1D("hNchTSGen_RecGreaterThanZero","; #it{R}_{T}; Entries",nBinsRT,binsRT);
		hNchTSGen_1->Sumw2();
		fListOfObjects->Add(hNchTSGen_1);

		hNchTSContamination = new TH1D("hNchTSContamination","; #it{R}_{T}; Entries",nBinsRT,binsRT);
		hNchTSContamination->Sumw2();
		fListOfObjects->Add(hNchTSContamination);

		hNchTSRecAll = new TH1D("hNchTSRecAll","; #it{R}_{T}; Entries",nBinsRT,binsRT);
		hNchTSRecAll->Sumw2();
		fListOfObjects->Add(hNchTSRecAll);

		hNchTSGenTest = new TH1D("hNchTSGenTest","; #it{R}_{T}; Entries",nBinsRT,binsRT);
		hNchTSGenTest->Sumw2();
		fListOfObjects->Add(hNchTSGenTest);

		hNchTSRecTest = new TH1D("hNchTSRecTest","; #it{R}_{T}; Entries",nBinsRT,binsRT);
		hNchTSRecTest->Sumw2();
		fListOfObjects->Add(hNchTSRecTest);

		for( int i = 0; i < nRegion; ++i ){

			hPhiGen[i]= new TH1D(Form("hPhiGen_%s",Region[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			hPhiGen[i]->Sumw2();
			fListOfObjects->Add(hPhiGen[i]);

			hPhiRec[i] = new TH1D(Form("hPhiRec_%s",Region[i]),"",64,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			hPhiRec[i]->Sumw2();
			fListOfObjects->Add(hPhiRec[i]);

			hNchGenGTZVsPtGen[i] = new TH2D(Form("hNchGenGTZVsPtGen_%s_Charged",Region[i]),";#it{R}_{T};#it{p}_{T}^{gen};",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenGTZVsPtGen[i]->Sumw2();
			fListOfObjects->Add(hNchGenGTZVsPtGen[i]);

			hNchGenGTZVsPtRec[i] = new TH2D(Form("hNchGenGTZVsPtRec_%s_Charged",Region[i]),";#it{R}_{T};#it{p}_{T}^{gen};",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenGTZVsPtRec[i]->Sumw2();
			fListOfObjects->Add(hNchGenGTZVsPtRec[i]);

			for(int pid = 0; pid < nPid; ++pid){
				hNchGenVsPtGenPID[i][pid] = new TH2D(Form("hNchGenVsPtGen_%s_%s",Region[i],Pid[pid]),";#it{R}_{T};#it{p}_{T}^{gen};",nBinsRT,binsRT,nPtBins,ptBins);
				hNchGenVsPtGenPID[i][pid]->Sumw2();
				fListOfObjects->Add(hNchGenVsPtGenPID[i][pid]);

				hNchGenVsPtRec[i][pid] = new TH2D(Form("hNchGenVsPtRec_%s_%s",Region[i],Pid[pid]),";#it{R}_{T};#it{p}_{T}^{gen};",nBinsRT,binsRT,nPtBins,ptBins);
				hNchGenVsPtRec[i][pid]->Sumw2();
				fListOfObjects->Add(hNchGenVsPtRec[i][pid]);

			}
		}

		for(int pid = 0; pid < nPid; pid++){

			hPtResponsePID[pid] = new TH2D(Form("hPtResponse_%s",Pid[pid]),"; #it{p}_{T}^{rec}; #it{p}_{T}^{gen}",nPtBins,ptBins,nPtBins,ptBins);
			hPtResponsePID[pid]->Sumw2();
			fListOfObjects->Add(hPtResponsePID[pid]);

			hNchGenVsPtGenIn[pid] = new TH2D(Form("hNchGenVsPtGenIn_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{gen}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtGenIn[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtGenIn[pid]);

			hNchGenVsPtGenPosIn[pid] = new TH2D(Form("hNchGenVsPtGenPosIn_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{gen}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtGenPosIn[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtGenPosIn[pid]);

			hNchGenVsPtGenNegIn[pid] = new TH2D(Form("hNchGenVsPtGenNegIn_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{gen}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtGenNegIn[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtGenNegIn[pid]);

			hNchGenVsPtRecIn[pid] = new TH2D(Form("hNchGenVsPtRecIn_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtRecIn[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtRecIn[pid]);

			hNchGenVsPtRecPosIn[pid] = new TH2D(Form("hNchGenVsPtRecPosIn_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtRecPosIn[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtRecPosIn[pid]);

			hNchGenVsPtRecNegIn[pid] = new TH2D(Form("hNchGenVsPtRecNegIn_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtRecNegIn[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtRecNegIn[pid]);

			hNchRecVsPtGenIn[pid] = new TH2D(Form("hNchRecVsPtGenIn_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchRecVsPtGenIn[pid]->Sumw2();
			fListOfObjects->Add(hNchRecVsPtGenIn[pid]);

			hNchGenVsPtRecInTOF[pid] = new TH2D(Form("hNchGenVsPtRecInTOF_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtRecInTOF[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtRecInTOF[pid]);

			hNchGenVsPtRecPosInTOF[pid] = new TH2D(Form("hNchGenVsPtRecPosInTOF_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtRecPosInTOF[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtRecPosInTOF[pid]);

			hNchGenVsPtRecNegInTOF[pid] = new TH2D(Form("hNchGenVsPtRecNegInTOF_%s",Pid[pid]),"; #it{R}_{T}^{gen}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
			hNchGenVsPtRecNegInTOF[pid]->Sumw2();
			fListOfObjects->Add(hNchGenVsPtRecNegInTOF[pid]);

		}

		hNchResponse = new TH2D("hNchResponse","; #it{R}_{T}^{rec}; #it{R}_{T}^{gen}",nBinsRT,binsRT,nBinsRT,binsRT);
		hNchResponse->Sumw2();
		fListOfObjects->Add(hNchResponse);

		hNchRMvsPt = new TH3D("hNchRMvsPt","; #it{R}_{T}^{rec}; #it{R}_{T}^{gen}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nBinsRT,binsRT,nPtBins,ptBins);
		hNchRMvsPt->Sumw2();
		fListOfObjects->Add(hNchRMvsPt);

		hNchRecVsPtRecOut = new TH2D("hNchRecVsPtRecOut","; #it{R}_{T}^{rec}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
		hNchRecVsPtRecOut->Sumw2();
		fListOfObjects->Add(hNchRecVsPtRecOut);

		hNchRecVsPtGenOut = new TH2D("hNchRecVsPtGenOut","; #it{R}_{T}^{rec}; #it{p}_{T}^{rec}",nBinsRT,binsRT,nPtBins,ptBins);
		hNchRecVsPtGenOut->Sumw2();
		fListOfObjects->Add(hNchRecVsPtGenOut);

		hPtPriGen = new TH1D("hPtPriGen","; #it{p}_{T}; Entries",nPtBins,ptBins);
		hPtPriGen->Sumw2();
		fListOfObjects->Add(hPtPriGen);

		hPtRec = new TH1D("hPtRec","; #it{p}_{T}; Entries",nPtBins,ptBins);
		hPtRec->Sumw2();
		fListOfObjects->Add(hPtRec);

		hPtPriRec = new TH1D("hPtPriRec","; #it{p}_{T}; Entries",nPtBins,ptBins);
		hPtPriRec->Sumw2();
		fListOfObjects->Add(hPtPriRec);

		hPtSecRec = new TH1D("hPtSecRec","; #it{p}_{T}; Entries",nPtBins,ptBins);
		hPtSecRec->Sumw2();
		fListOfObjects->Add(hPtSecRec);

	}

	fEventCuts.AddQAplotsToList(fListOfObjects);
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

	printf("---- fGenLeadPt   =  %f     fRecLeadPt  =  %f\n",fGenLeadPt,fRecLeadPt);

	if(fIsMCclosure){
		double randomUE = gRandom->Uniform(0.0,1.0);
		if(randomUE<0.5){// corrections (50% stat.)
			if(isGoodVtxPosMC){
				if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
					GetDetectorResponse();
				}
				if(fGenLeadPt>=fPtMin){
					GetMCCorrections();
				}
			}
		}
		else{// for testing the method
			if( ( fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax ) && ( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax )){
				GetMultiplicityDistributions();
				////ProduceArrayTrksESD();
			}
		}
	}
	else{
		if(fAnalysisMC){
			if(isGoodVtxPosMC){
				if( (fGenLeadPt>=fLeadPtCutMin && fGenLeadPt<fLeadPtCutMax)&&(fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax)){
					GetMultiplicityDistributions();
					GetDetectorResponse();
				}

				// UE analysis
				if(fGenLeadPt>=fPtMin){
					GetMCCorrections();
				}
			}
		}
		else{
			// KNO scaling
			if(( fRecLeadPt>=fLeadPtCutMin && fRecLeadPt<fLeadPtCutMax ))
				ProduceArrayTrksESD();

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
			if( track->Pt() < fPtMin) continue;

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

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));  
		if (!track) continue;
		if (i==fRecLeadIn) continue;
		if (!fTrackFilter->IsSelected(track)) continue;
		if (track->Charge() == 0 ) continue;
		if (TMath::Abs(track->Eta()) > fEtaCut) continue;
		if (track->Pt() < fPtMin) continue;

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

	double RTgen = (double)multTSgen/fMeanMultTSMCGen;

	hMultTSGen->Fill(multTSgen);
	hNchTSGenTest->Fill(RTgen);
	fPtLVsNchGen->Fill(fGenLeadPt,RTgen);

	//	Filling Nch vs pT True
	for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (i==fGenLeadIn) continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;
		if (particle->Charge() == 0) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut) continue;
		if (particle->Pt() < fPtMin) continue;

		int pdgCode = particle->PdgCode();
		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		double DPhi = DeltaPhi(particle->Phi(), fGenLeadPhi);

		hNchGenVsPtGenPID[3][0]->Fill(RTgen,particle->Pt());
		if(multTSrec > 0)hNchGenGTZVsPtGen[3]->Fill(RTgen,particle->Pt());

		if(pidCodeMC==1)
			hNchGenVsPtGenPID[3][1]->Fill(RTgen,particle->Pt());
		if(pidCodeMC==2)
			hNchGenVsPtGenPID[3][2]->Fill(RTgen,particle->Pt());
		if(pidCodeMC==3)
			hNchGenVsPtGenPID[3][3]->Fill(RTgen,particle->Pt());

		if(TMath::Abs(DPhi)<pi/3.0){

			hNchGenVsPtGenPID[0][0]->Fill(RTgen,particle->Pt());
			if(multTSrec > 0)hNchGenGTZVsPtGen[0]->Fill(RTgen,particle->Pt());

			if(pidCodeMC==1)
				hNchGenVsPtGenPID[0][1]->Fill(RTgen,particle->Pt());
			else if(pidCodeMC==2)
				hNchGenVsPtGenPID[0][2]->Fill(RTgen,particle->Pt());
			else if(pidCodeMC==3)
				hNchGenVsPtGenPID[0][3]->Fill(RTgen,particle->Pt());
			else
				continue;
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){

			hNchGenVsPtGenPID[1][0]->Fill(RTgen,particle->Pt());
			if(multTSrec > 0)hNchGenGTZVsPtGen[1]->Fill(RTgen,particle->Pt());

			if(pidCodeMC==1)
				hNchGenVsPtGenPID[1][1]->Fill(RTgen,particle->Pt());
			else if(pidCodeMC==2)
				hNchGenVsPtGenPID[1][2]->Fill(RTgen,particle->Pt());
			else if(pidCodeMC==3)
				hNchGenVsPtGenPID[1][3]->Fill(RTgen,particle->Pt());
			else
				continue;
		}
		else{

			hNchGenVsPtGenPID[2][0]->Fill(RTgen,particle->Pt());
			if(multTSrec > 0)hNchGenGTZVsPtGen[2]->Fill(RTgen,particle->Pt());

			if(pidCodeMC==1)
				hNchGenVsPtGenPID[2][1]->Fill(RTgen,particle->Pt());
			else if(pidCodeMC==2)
				hNchGenVsPtGenPID[2][2]->Fill(RTgen,particle->Pt());
			else if(pidCodeMC==3)
				hNchGenVsPtGenPID[2][3]->Fill(RTgen,particle->Pt());
			else
				continue;
		}

	}

	double RTrec = (double)multTSrec/fMeanMultTSMCRec;

	hMultTSRec->Fill(multTSrec);
	hNchTSRecTest->Fill(RTrec);
	fPtLVsNchRec->Fill(fRecLeadPt,RTrec);

	if(multTSrec > 0)
		hNchTSGen_1->Fill(RTgen);

	// Filling rec pT vs UE (for pT I use 2015 track cuts, UE uses TPC-only)
	for(int i=0; i < iTracks; i++){  

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
		if (!track) continue;
		if (i==fRecLeadIn) continue;
		if (!fTrackFilterGolden->IsSelected(track)) continue;
		if (track->Charge() == 0 ) continue;
		if (TMath::Abs(track->Eta()) > fEtaCut) continue;
		if (track->Pt() < fPtMin) continue;

		const int label = TMath::Abs(track->GetLabel());
		AliMCParticle *trackMC = (AliMCParticle*)fMC->GetTrack(label);

		int pdgCode = trackMC->PdgCode();
		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		bool IsTOFout = kFALSE;
		IsTOFout = TOFPID(track);

		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		hNchGenVsPtRec[3][0]->Fill(RTgen,track->Pt());
		if(multTSrec > 0)hNchGenGTZVsPtRec[3]->Fill(RTgen,track->Pt());
		hNchVsPtDataTPC[3][0]->Fill(RTrec,track->Pt(),0.25);
		if(IsTOFout)hNchVsPtDataTOF[3][0]->Fill(RTrec,track->Pt(),0.25);

		if(pidCodeMC==1){
			hNchGenVsPtRec[3][1]->Fill(RTgen,track->Pt());
			hNchVsPtDataTPC[3][1]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
			if(IsTOFout)hNchVsPtDataTOF[3][1]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion));
		}
		if(pidCodeMC==2){
			hNchGenVsPtRec[3][2]->Fill(RTgen,track->Pt());
			hNchVsPtDataTPC[3][2]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
			if(IsTOFout)hNchVsPtDataTOF[3][2]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
		}
		if(pidCodeMC==3){
			hNchGenVsPtRec[3][3]->Fill(RTgen,track->Pt());
			hNchVsPtDataTPC[3][3]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
			if(IsTOFout)hNchVsPtDataTOF[3][3]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton));
		}

		if(TMath::Abs(DPhi)<pi/3.0){

			hNchGenVsPtRec[0][0]->Fill(RTgen,track->Pt());
			if(multTSrec > 0)hNchGenGTZVsPtRec[0]->Fill(RTgen,track->Pt());
			hNchVsPtDataTPC[0][0]->Fill(RTrec,track->Pt(),0.25);
			if(IsTOFout)hNchVsPtDataTOF[0][0]->Fill(RTrec,track->Pt(),0.25);

			if(pidCodeMC==1){
				hNchGenVsPtRec[0][1]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[0][1]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
				if(IsTOFout)hNchVsPtDataTOF[0][1]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion));
			}
			else if(pidCodeMC==2){
				hNchGenVsPtRec[0][2]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[0][2]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
				if(IsTOFout)hNchVsPtDataTOF[0][2]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
			}
			else if(pidCodeMC==3){
				hNchGenVsPtRec[0][3]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[0][3]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
				if(IsTOFout)hNchVsPtDataTOF[0][3]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton));
			}
			else 
				continue;
		}

		else if(TMath::Abs(DPhi-pi)<pi/3.0){

			hNchGenVsPtRec[1][0]->Fill(RTgen,track->Pt());
			if(multTSrec > 0)hNchGenGTZVsPtRec[1]->Fill(RTgen,track->Pt());
			hNchVsPtDataTPC[1][0]->Fill(RTrec,track->Pt(),0.25);
			if(IsTOFout)hNchVsPtDataTOF[1][0]->Fill(RTrec,track->Pt(),0.25);

			if(pidCodeMC==1){
				hNchGenVsPtRec[1][1]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[1][1]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
				if(IsTOFout)hNchVsPtDataTOF[1][1]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion));
			}
			else if(pidCodeMC==2){
				hNchGenVsPtRec[1][2]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[1][2]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
				if(IsTOFout)hNchVsPtDataTOF[1][2]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
			}
			else if(pidCodeMC==3){
				hNchGenVsPtRec[1][3]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[1][3]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
				if(IsTOFout)hNchVsPtDataTOF[1][3]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton));
			}
			else 
				continue;
		}

		else{
			hNchGenVsPtRec[2][0]->Fill(RTgen,track->Pt());
			if(multTSrec > 0)hNchGenGTZVsPtRec[2]->Fill(RTgen,track->Pt());
			hNchVsPtDataTPC[2][0]->Fill(RTrec,track->Pt(),0.25);
			if(IsTOFout)hNchVsPtDataTOF[2][0]->Fill(RTrec,track->Pt(),0.25);

			if(pidCodeMC==1){
				hNchGenVsPtRec[2][1]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[2][1]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
				if(IsTOFout)hNchVsPtDataTOF[2][1]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion));
			}
			else if(pidCodeMC==2){
				hNchGenVsPtRec[2][2]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[2][2]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
				if(IsTOFout)hNchVsPtDataTOF[2][2]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon));
			}
			else if(pidCodeMC==3){
				hNchGenVsPtRec[2][3]->Fill(RTgen,track->Pt());
				hNchVsPtDataTPC[2][3]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
				if(IsTOFout)hNchVsPtDataTOF[2][3]->Fill(RTrec,track->Pt(),fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton));
			}
			else 
				continue;
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
			hPhiGen[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hPhiGen[1]->Fill(DPhi);
		}
		else{
			multTSgen++;
			hPhiGen[2]->Fill(DPhi);
		}

		hPhiGen[3]->Fill(DPhi);
	}

	int iTracks = 0;          
	iTracks = fESD->GetNumberOfTracks();          

	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!track) continue;
		if(!fTrackFilter->IsSelected(track)) continue;
		if (track->Charge() == 0 ) continue;
		if(TMath::Abs(track->Eta()) > fEtaCut) continue;
		if( track->Pt() < fPtMin) continue;

		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			hPhiRec[0]->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hPhiRec[1]->Fill(DPhi);
		}
		else{
			multTSrec++;
			hPhiRec[2]->Fill(DPhi);
		}

		hPhiRec[3]->Fill(DPhi);
	}

	double RTgen = (double)multTSgen/fMeanMultTSMCGen;
	double RTrec = (double)multTSrec/fMeanMultTSMCRec;

	hNchResponse->Fill(RTrec,RTgen);

	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!track) continue;
		if(!fTrackFilterGolden->IsSelected(track)) continue;
		if (track->Charge() == 0 ) continue;
		if(TMath::Abs(track->Eta()) > fEtaCut) continue;
		if( track->Pt() < fPtMin) continue;

		int label = TMath::Abs(track->GetLabel());
		AliMCParticle *trackMC = (AliMCParticle*)fMC->GetTrack(label);

		int pdgCode = trackMC->PdgCode();
		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		hPtResponsePID[0]->Fill(track->Pt(),trackMC->Pt());

		if(pidCodeMC==1)	
			hPtResponsePID[1]->Fill(track->Pt(),trackMC->Pt());

		else if(pidCodeMC==2)	
			hPtResponsePID[2]->Fill(track->Pt(),trackMC->Pt());

		else if(pidCodeMC==3)	
			hPtResponsePID[3]->Fill(track->Pt(),trackMC->Pt());

		else
			continue;

	}

		printf(" iTracks = %d\n",iTracks);

	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!track) continue;
		if(!fTrackFilter->IsSelected(track)) continue;
		if (track->Charge() == 0 ) continue;
		if(TMath::Abs(track->Eta()) > fEtaCut) continue;
		if( track->Pt() < fPtMin) continue;

		hNchRMvsPt->Fill(RTrec,RTgen,track->Pt());

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

	double RTgen = (double)multTSgen/fMeanMultTSMCGen;
	hNchTSGen->Fill(RTgen);

	int iTracks(fESD->GetNumberOfTracks());          
	for(int i = 0; i < iTracks; i++){              

		if(i==fRecLeadIn) continue;
		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 
		if(!track) continue;
		if(!fTrackFilter->IsSelected(track)) continue;
		if(TMath::Abs(track->Eta()) > fEtaCut) continue;
		if(track->Pt() < fPtMin) continue;
		if(track->Charge()==0 ) continue;

		const int label = TMath::Abs(track->GetLabel());
		if( !fMC->IsPhysicalPrimary(label) ) continue;

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

	double RTrec = (double)multTSrec/fMeanMultTSMCRec;
	if(multTSrec > 0){
		hNchTSRec->Fill(RTrec);
		//		hNchTSGen_1->Fill(RTgen);
	}

	hNchTSRecAll->Fill(RTrec);
	if(RTgen==0)hNchTSContamination->Fill(RTrec);

	for (int i = 0; i < fMC->GetNumberOfTracks(); i++){

		AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
		if (!particle) continue;
		if (particle->Charge() == 0) continue;
		if (TMath::Abs(particle->Eta()) > fEtaCut)continue;
		if (particle->Pt() < fPtMin)continue;
		if (!fMC->IsPhysicalPrimary(i)) continue;

		hPtPriGen->Fill(particle->Pt());

		int pdgCode = particle->PdgCode();
		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		hNchGenVsPtGenIn[0]->Fill(RTgen,particle->Pt());
		hNchRecVsPtGenIn[0]->Fill(RTrec,particle->Pt());
		hNchRecVsPtGenOut->Fill(RTrec,particle->Pt());

		if(particle->Charge() > 0)
			hNchGenVsPtGenPosIn[0]->Fill(RTgen,particle->Pt());

		if(particle->Charge() < 0)
			hNchGenVsPtGenNegIn[0]->Fill(RTgen,particle->Pt());

		if(pidCodeMC==1){

			hNchGenVsPtGenIn[1]->Fill(RTgen,particle->Pt());
			hNchRecVsPtGenIn[1]->Fill(RTrec,particle->Pt());

			if(particle->Charge() > 0)
				hNchGenVsPtGenPosIn[1]->Fill(RTgen,particle->Pt());

			if(particle->Charge() < 0)
				hNchGenVsPtGenNegIn[1]->Fill(RTgen,particle->Pt());

		}
		else if(pidCodeMC==2){

			hNchGenVsPtGenIn[2]->Fill(RTgen,particle->Pt());
			hNchRecVsPtGenIn[2]->Fill(RTrec,particle->Pt());

			if(particle->Charge() > 0)
				hNchGenVsPtGenPosIn[2]->Fill(RTgen,particle->Pt());

			if(particle->Charge() < 0)
				hNchGenVsPtGenNegIn[2]->Fill(RTgen,particle->Pt());
		}
		else if(pidCodeMC==3){

			hNchGenVsPtGenIn[3]->Fill(RTgen,particle->Pt());
			hNchRecVsPtGenIn[3]->Fill(RTrec,particle->Pt());

			if(particle->Charge() > 0)
				hNchGenVsPtGenPosIn[3]->Fill(RTgen,particle->Pt());

			if(particle->Charge() < 0)
				hNchGenVsPtGenNegIn[3]->Fill(RTgen,particle->Pt());
		}
		else
			continue;
	}

	for(int i = 0; i < iTracks; i++){                

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 

		if (!track) continue;        
		if (track->Charge() == 0 ) continue;
		if (!fTrackFilterGolden->IsSelected(track)) continue;        
		if (TMath::Abs(track->Eta()) > fEtaCut) continue;        
		if (track->Pt() < fPtMin) continue;

		const int label = TMath::Abs(track->GetLabel());
		AliMCParticle *trackMC = (AliMCParticle*)fMC->GetTrack(label);

		int pdgCode = trackMC->PdgCode();
		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		hPtRec->Fill(track->Pt());

		if( fMC->IsPhysicalPrimary(label) ){
			hPtPriRec->Fill(track->Pt());
			hNchRecVsPtRecOut->Fill(RTrec,track->Pt());

			hNchGenVsPtRecIn[0]->Fill(RTgen,track->Pt());

			if(track->Charge() > 0)
				hNchGenVsPtRecPosIn[0]->Fill(RTgen,track->Pt());

			if(track->Charge() < 0)
				hNchGenVsPtRecNegIn[0]->Fill(RTgen,track->Pt());


			if(TOFPID(track)){

				hNchGenVsPtRecInTOF[0]->Fill(RTgen,track->Pt());

				if(track->Charge() > 0)
					hNchGenVsPtRecPosInTOF[0]->Fill(RTgen,track->Pt());

				if(track->Charge() < 0)
					hNchGenVsPtRecNegInTOF[0]->Fill(RTgen,track->Pt());

			}

			if(pidCodeMC==1){
				hNchGenVsPtRecIn[1]->Fill(RTgen,track->Pt());

				if(track->Charge() > 0)
					hNchGenVsPtRecPosIn[1]->Fill(RTgen,track->Pt());

				if(track->Charge() < 0)
					hNchGenVsPtRecNegIn[1]->Fill(RTgen,track->Pt());

				if(TOFPID(track)){

					hNchGenVsPtRecInTOF[1]->Fill(RTgen,track->Pt());

					if(track->Charge() > 0)
						hNchGenVsPtRecPosInTOF[1]->Fill(RTgen,track->Pt());

					if(track->Charge() < 0)
						hNchGenVsPtRecNegInTOF[1]->Fill(RTgen,track->Pt());

				}	
			}

			else if(pidCodeMC==2){

				hNchGenVsPtRecIn[2]->Fill(RTgen,track->Pt());

				if(track->Charge() > 0)
					hNchGenVsPtRecPosIn[2]->Fill(RTgen,track->Pt());

				if(track->Charge() < 0)
					hNchGenVsPtRecNegIn[2]->Fill(RTgen,track->Pt());

				if(TOFPID(track)){

					hNchGenVsPtRecInTOF[2]->Fill(RTgen,track->Pt());

					if(track->Charge() > 0)
						hNchGenVsPtRecPosInTOF[2]->Fill(RTgen,track->Pt());

					if(track->Charge() < 0)
						hNchGenVsPtRecNegInTOF[2]->Fill(RTgen,track->Pt());
				}	
			}

			else if(pidCodeMC==3){

				hNchGenVsPtRecIn[3]->Fill(RTgen,track->Pt());

				if(track->Charge() > 0)
					hNchGenVsPtRecPosIn[3]->Fill(RTgen,track->Pt());

				if(track->Charge() < 0)
					hNchGenVsPtRecNegIn[3]->Fill(RTgen,track->Pt());

				if(TOFPID(track)){

					hNchGenVsPtRecInTOF[3]->Fill(RTgen,track->Pt());

					if(track->Charge() > 0)
						hNchGenVsPtRecPosInTOF[3]->Fill(RTgen,track->Pt());

					if(track->Charge() < 0)
						hNchGenVsPtRecNegInTOF[3]->Fill(RTgen,track->Pt());
				}
			}

			else
				continue;
		}
		if( fMC->IsSecondaryFromWeakDecay(label) || fMC->IsSecondaryFromMaterial(label)){
			hPtSecRec->Fill(track->Pt());
		}

	}

}
//_____________________________________________________________________________
double AliAnalysisTaskSpectraRT::DeltaPhi(Double_t phi, Double_t Lphi,
		Double_t rangeMin, Double_t rangeMax)
{

	Double_t dphi = -999;
	Double_t pi = TMath::Pi();
	if(Lphi > 2*pi || Lphi < 0)cout << "Lphi :: " << Lphi << endl;
	if(phi  > 2*pi || phi < 0)cout << "phi = " << phi << endl;

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

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i)); 

		if(!track) continue;

		if(!fTrackFilter->IsSelected(track)) continue;

		if(TMath::Abs(track->Eta()) > fEtaCut) continue;

		if( track->Pt() < fPtMin) continue;

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

		hPhiData[3]->Fill(DPhi);
	}

	hNchTSData->Fill(multTSdata);

	double RT = multTSdata/fMeanChT;
	hPtLVsRT->Fill(fRecLeadPt,RT);
	hRTData->Fill(RT);

	for(int iT = 0; iT < iTracks; iT++) {

		if(iT==fRecLeadIn) continue;
		AliESDtrack* esdTrack = (AliESDtrack*)fESD->GetTrack(iT);

		if(TMath::Abs(esdTrack->Eta()) > fEtaCut) continue;
		////		if(esdTrack->GetTPCsignalN() < fNcl) continue;
		if(esdTrack->Pt() < fPtMin) continue;
		if(!fTrackFilterGolden->IsSelected(esdTrack)) continue;

		double DPhi = DeltaPhi(esdTrack->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){

			hNchVsPtDataTPC[0][0]->Fill(RT,esdTrack->Pt(),0.25);
			hNchVsPtDataTPC[0][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hNchVsPtDataTPC[0][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hNchVsPtDataTPC[0][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));

			if(esdTrack->Charge() > 0){
				hNchVsPtDataPosTPC[0][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataPosTPC[0][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));	
				hNchVsPtDataPosTPC[0][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));	
				hNchVsPtDataPosTPC[0][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
			}

			if(esdTrack->Charge() < 0){
				hNchVsPtDataNegTPC[0][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataNegTPC[0][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));	
				hNchVsPtDataNegTPC[0][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));	
				hNchVsPtDataNegTPC[0][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));	
			}
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){

			hNchVsPtDataTPC[1][0]->Fill(RT,esdTrack->Pt(),0.25);
			hNchVsPtDataTPC[1][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hNchVsPtDataTPC[1][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hNchVsPtDataTPC[1][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));

			if(esdTrack->Charge() > 0){
				hNchVsPtDataPosTPC[1][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataPosTPC[1][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));	
				hNchVsPtDataPosTPC[1][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));	
				hNchVsPtDataPosTPC[1][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
			}

			if(esdTrack->Charge() < 0){
				hNchVsPtDataNegTPC[1][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataNegTPC[1][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));	
				hNchVsPtDataNegTPC[1][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));	
				hNchVsPtDataNegTPC[1][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));	
			}
		}
		else{
			hNchVsPtDataTPC[2][0]->Fill(RT,esdTrack->Pt(),0.25);
			hNchVsPtDataTPC[2][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
			hNchVsPtDataTPC[2][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
			hNchVsPtDataTPC[2][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));

			if(esdTrack->Charge() > 0){
				hNchVsPtDataPosTPC[2][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataPosTPC[2][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));	
				hNchVsPtDataPosTPC[2][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));	
				hNchVsPtDataPosTPC[2][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
			}

			if(esdTrack->Charge() < 0){
				hNchVsPtDataNegTPC[2][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataNegTPC[2][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));	
				hNchVsPtDataNegTPC[2][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));	
				hNchVsPtDataNegTPC[2][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));	
			}
		}

		hNchVsPtDataTPC[3][0]->Fill(RT,esdTrack->Pt(),0.25);
		hNchVsPtDataTPC[3][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));
		hNchVsPtDataTPC[3][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));
		hNchVsPtDataTPC[3][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));

		if(esdTrack->Charge() > 0){
			hNchVsPtDataPosTPC[3][0]->Fill(RT,esdTrack->Pt(),0.25);	
			hNchVsPtDataPosTPC[3][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));	
			hNchVsPtDataPosTPC[3][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));	
			hNchVsPtDataPosTPC[3][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));
		}

		if(esdTrack->Charge() < 0){
			hNchVsPtDataNegTPC[3][0]->Fill(RT,esdTrack->Pt(),0.25);	
			hNchVsPtDataNegTPC[3][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kPion));	
			hNchVsPtDataNegTPC[3][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kKaon));	
			hNchVsPtDataNegTPC[3][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTPC(esdTrack,AliPID::kProton));	
		}

		//_______________________________ TOF PID

		bool IsTOFout = kFALSE;
		IsTOFout = TOFPID(esdTrack);
		if(!IsTOFout) continue;

		if(TMath::Abs(DPhi)<pi/3.0){

			hNchVsPtDataTOF[0][0]->Fill(RT,esdTrack->Pt(),0.25);
			hNchVsPtDataTOF[0][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));
			hNchVsPtDataTOF[0][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));
			hNchVsPtDataTOF[0][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));

			if(esdTrack->Charge() > 0){
				hNchVsPtDataPosTOF[0][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataPosTOF[0][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));	
				hNchVsPtDataPosTOF[0][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));	
				hNchVsPtDataPosTOF[0][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));
			}

			if(esdTrack->Charge() < 0){
				hNchVsPtDataNegTOF[0][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataNegTOF[0][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));	
				hNchVsPtDataNegTOF[0][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));	
				hNchVsPtDataNegTOF[0][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));	
			}
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){

			hNchVsPtDataTOF[1][0]->Fill(RT,esdTrack->Pt(),0.25);
			hNchVsPtDataTOF[1][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));
			hNchVsPtDataTOF[1][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));
			hNchVsPtDataTOF[1][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));

			if(esdTrack->Charge() > 0){
				hNchVsPtDataPosTOF[1][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataPosTOF[1][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));	
				hNchVsPtDataPosTOF[1][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));	
				hNchVsPtDataPosTOF[1][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));
			}

			if(esdTrack->Charge() < 0){
				hNchVsPtDataNegTOF[1][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataNegTOF[1][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));	
				hNchVsPtDataNegTOF[1][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));	
				hNchVsPtDataNegTOF[1][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));	
			}
		}
		else{
			hNchVsPtDataTOF[2][0]->Fill(RT,esdTrack->Pt(),0.25);
			hNchVsPtDataTOF[2][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));
			hNchVsPtDataTOF[2][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));
			hNchVsPtDataTOF[2][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));

			if(esdTrack->Charge() > 0){
				hNchVsPtDataPosTOF[2][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataPosTOF[2][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));	
				hNchVsPtDataPosTOF[2][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));	
				hNchVsPtDataPosTOF[2][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));
			}

			if(esdTrack->Charge() < 0){
				hNchVsPtDataNegTOF[2][0]->Fill(RT,esdTrack->Pt(),0.25);	
				hNchVsPtDataNegTOF[2][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));	
				hNchVsPtDataNegTOF[2][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));	
				hNchVsPtDataNegTOF[2][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));	
			}
		}

		hNchVsPtDataTOF[3][0]->Fill(RT,esdTrack->Pt(),0.25);
		hNchVsPtDataTOF[3][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));
		hNchVsPtDataTOF[3][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));
		hNchVsPtDataTOF[3][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));

		if(esdTrack->Charge() > 0){
			hNchVsPtDataPosTOF[3][0]->Fill(RT,esdTrack->Pt(),0.25);	
			hNchVsPtDataPosTOF[3][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));	
			hNchVsPtDataPosTOF[3][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));	
			hNchVsPtDataPosTOF[3][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));
		}

		if(esdTrack->Charge() < 0){
			hNchVsPtDataNegTOF[3][0]->Fill(RT,esdTrack->Pt(),0.25);	
			hNchVsPtDataNegTOF[3][1]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kPion));	
			hNchVsPtDataNegTOF[3][2]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kKaon));	
			hNchVsPtDataNegTOF[3][3]->Fill(RT,esdTrack->Pt(),fPIDResponse->NumberOfSigmasTOF(esdTrack,AliPID::kProton));	
		}

	}//end of track loop
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
		printf("Setting TPCOnly Track Cuts\n");
	}
	else{
		esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,1);
		printf("Setting ITSTPC2011 Track Cuts\n");
	}

	fTrackFilter->AddCuts(esdTrackCuts);
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
