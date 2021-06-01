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
		fPtMinCut(0.15),
		fPtMaxCut(50.0),
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
		fTrackCuts(0),
		fPeriod("16l"),
		hphiso(0x0),
		hetaso(0x0),
		hTruthPhiSo(0x0),
		hTruthEtaSo(0x0),
		hSOrvsV0M(0x0),
		hSOrvsTrks(0x0),
		hMultPercvsNch(0x0),
		hSOtvsSOrV0M(0x0),
		hSOtvsSOrCL1(0x0),
		fEtaCalibrationNeg(0x0),
		fEtaCalibrationPos(0x0),
		fcutDCAxy(0x0),
		fcutLow(0x0),
		fcutHigh(0x0),
		fRecLeadPhi(0.0),
		fRecLeadPt(0.0),
		fRecLeadIn(0),
		fNT(0.0),
		hPhiStandard(0x0),
		hPhiHybrid1(0x0),
		hPhiTotal(0x0),
		fGeometricalCut(0x0),
		fHybridTrackCuts1(0x0),
		hV0MvsSOvsNT(0x0),
		hCL1vsSOvsNT(0x0),
		hV0MvsSOvsNoNT(0x0),
		hCL1vsSOvsNoNT(0x0),
		hV0MvsNoSOvsNT(0x0),
		hCL1vsNoSOvsNT(0x0),
		hV0MvsNoSOvsNoNT(0x0),
		hCL1vsNoSOvsNoNT(0x0),
		hPhiToward(0x0),
		hPhiAway(0x0),
		hPhiTransverse(0x0)

{

	hPionSimvsSOV0M = 0x0;
	hPionGenvsSOV0M = 0x0;
	hKaonSimvsSOV0M = 0x0;
	hKaonGenvsSOV0M = 0x0;
	hProtonSimvsSOV0M = 0x0;
	hProtonGenvsSOV0M = 0x0;

	hPionSimvsSOV0M2 = 0x0;
	hKaonSimvsSOV0M2 = 0x0;
	hProtonSimvsSOV0M2 = 0x0;
	hPionSimvsSOCL12 = 0x0;
	hKaonSimvsSOCL12 = 0x0;
	hProtonSimvsSOCL12 = 0x0;

	hPionSimvsSOCL1 = 0x0;
	hPionGenvsSOCL1 = 0x0;
	hKaonSimvsSOCL1 = 0x0;
	hKaonGenvsSOCL1 = 0x0;
	hProtonSimvsSOCL1 = 0x0;
	hProtonGenvsSOCL1 = 0x0;

	hPhiSimV0M = 0x0;
	hPhiGenV0M = 0x0;
	hPhiSimCL1 = 0x0;
	hPhiGenCL1 = 0x0;

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
	fPtMinCut(0.15),
	fPtMaxCut(50.0),
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
	fTrackCuts(0),
	fPeriod("16l"),
	hphiso(0x0),
	hetaso(0x0),
	hTruthPhiSo(0x0),
	hTruthEtaSo(0x0),
	hSOrvsV0M(0x0),
	hSOrvsTrks(0x0),
	hMultPercvsNch(0x0),
	hSOtvsSOrV0M(0x0),
	hSOtvsSOrCL1(0x0),
	fEtaCalibrationNeg(0x0),
	fEtaCalibrationPos(0x0),
	fcutDCAxy(0x0),
	fcutLow(0x0),
	fcutHigh(0x0),
	fRecLeadPhi(0.0),
	fRecLeadPt(0.0),
	fRecLeadIn(0),
	fNT(0.0),
	hPhiStandard(0x0),
	hPhiHybrid1(0x0),
	hPhiTotal(0x0),
	fGeometricalCut(0x0),
	fHybridTrackCuts1(0x0),
	hV0MvsSOvsNT(0x0),
	hCL1vsSOvsNT(0x0),
	hV0MvsSOvsNoNT(0x0),
	hCL1vsSOvsNoNT(0x0),
	hV0MvsNoSOvsNT(0x0),
	hCL1vsNoSOvsNT(0x0),
	hV0MvsNoSOvsNoNT(0x0),
	hCL1vsNoSOvsNoNT(0x0),
	hPhiToward(0x0),
	hPhiAway(0x0),
	hPhiTransverse(0x0)
{

	hPionSimvsSOV0M = 0x0;
	hPionGenvsSOV0M = 0x0;
	hKaonSimvsSOV0M = 0x0;
	hKaonGenvsSOV0M = 0x0;
	hProtonSimvsSOV0M = 0x0;
	hProtonGenvsSOV0M = 0x0;

	hPionSimvsSOV0M2 = 0x0;
	hKaonSimvsSOV0M2 = 0x0;
	hProtonSimvsSOV0M2 = 0x0;
	hPionSimvsSOCL12 = 0x0;
	hKaonSimvsSOCL12 = 0x0;
	hProtonSimvsSOCL12 = 0x0;

	hPionSimvsSOCL1 = 0x0;
	hPionGenvsSOCL1 = 0x0;
	hKaonSimvsSOCL1 = 0x0;
	hKaonGenvsSOCL1 = 0x0;
	hProtonSimvsSOCL1 = 0x0;
	hProtonGenvsSOCL1 = 0x0;

	hPhiSimV0M = 0x0;
	hPhiGenV0M = 0x0;
	hPhiSimCL1 = 0x0;
	hPhiGenCL1 = 0x0;

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

	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
}

AliAnalysisTaskSpherocity::~AliAnalysisTaskSpherocity() {

	//
	// Destructor
	//

	if(fListOfObjects) {
		delete fListOfObjects;
		fListOfObjects = 0x0;
	}
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
		//fTrackFilterGolden->SetMinNCrossedRowsTPC(70); // Varied for track selection
		fTrackFilterGolden->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		//fTrackFilterGolden->SetMaxChi2PerClusterTPC(4); // Varied for track selection
		fTrackFilterGolden->SetAcceptKinkDaughters(kFALSE);
		fTrackFilterGolden->SetRequireTPCRefit(kTRUE);
		fTrackFilterGolden->SetRequireITSRefit(kTRUE);
		fTrackFilterGolden->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
		fTrackFilterGolden->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
		//                fTrackFilterGolden->SetMaxChi2TPCConstrainedGlobal(36);
		//fTrackFilterGolden->SetMaxDCAToVertexZ(2); // Varied for track selection
		fTrackFilterGolden->SetDCAToVertex2D(kFALSE);
		fTrackFilterGolden->SetRequireSigmaToVertex(kFALSE);
		fTrackFilterGolden->SetMaxChi2PerClusterITS(36);
		fTrackFilterGolden->SetEtaRange(-0.8,0.8);

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

	if(!fHybridTrackCuts1){
		fHybridTrackCuts1 = new AliESDtrackCuts("fHybridTrackCuts1");
		fHybridTrackCuts1->SetMinNCrossedRowsTPC(70);
		fHybridTrackCuts1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
		fHybridTrackCuts1->SetMaxChi2PerClusterTPC(4);
		fHybridTrackCuts1->SetAcceptKinkDaughters(kFALSE);
		fHybridTrackCuts1->SetRequireTPCRefit(kTRUE);
		fHybridTrackCuts1->SetRequireITSRefit(kFALSE);
		fHybridTrackCuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
		fHybridTrackCuts1->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
		//fHybridTrackCuts1->SetMaxChi2TPCConstrainedGlobal(36);
		fHybridTrackCuts1->SetMaxDCAToVertexZ(2);
		fHybridTrackCuts1->SetDCAToVertex2D(kFALSE);
		fHybridTrackCuts1->SetRequireSigmaToVertex(kFALSE);
		fHybridTrackCuts1->SetMaxChi2PerClusterITS(36);

	}

	if(!fGeometricalCut){
		fGeometricalCut = new AliESDtrackCuts("fGeometricalCut");
		fGeometricalCut->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7);
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

	//fcent=new TH1F("fcent","fcent",13,0,13);
	//fcentAfterPrimaries =new TH1F("fcentAfterPrimaries","fcentAfterPrimaries",13,0,13);
	//	fcentAfterV0s =new TH1F("fcentAfterV0s","fcentAfterV0s",13,0,13);
	//	fListOfObjects->Add(fcent);
	//	fListOfObjects->Add(fcentAfterPrimaries);
	//	fListOfObjects->Add(fcentAfterV0s);

	//	const Int_t nDeltaPiBins   = 80;
	//	const Double_t deltaPiLow  = 20;
	//	const Double_t deltaPiHigh = 100;

	////	const Char_t *Pid[7]       = {"Ch","Pion","Kaon","Proton","Electron","Muon","Oher"};

	const int nPtBins = 63;
	float ptBins[nPtBins+1] = {
		0.01, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.25, 0.30, 0.35,
		0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
		0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
		1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40,
		3.60, 3.80, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 8.00,
		9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 18.00,
		20.00,22.00,24.00,26.00,30.00};

	const int nBinsSO = 1000;
	float BinsSO[nBinsSO+1] = {0};

	for(int i = 0; i <= nBinsSO; ++i){
		BinsSO[i] = 0.0 + i*0.001;
	}

	const int nBinsSO0 = 100;
	float BinsSO0[nBinsSO0+1] = {0};

	for(int i = 0; i <= nBinsSO0; ++i){
		BinsSO0[i] = 0.0 + i*0.01;
	}

	const int nBinsNsigma = 50;
	float BinsNsigma[nBinsNsigma+1] = {0.0};

	for(int i = 0; i <= nBinsNsigma; ++i){
		BinsNsigma[i] = -10.0+i*0.4;
	}

	const int nBetaBins   = 100;
	float BetaBins[nBetaBins+1] = { 0.0 };
	for(int i = 0; i <= nBetaBins; ++i){
		BetaBins[i] = 0.2+((double)i)/100.0;
	}

	const int nPhiBins = 100; 
	float PhiBins[nPhiBins+1] = {0.0};
	for(int i = 0; i <= nPhiBins; ++i){
		PhiBins[i] = -1.0*TMath::Pi() + ((float)i)*(0.01)*(8.0*TMath::Pi()/2);
	}

	const int nEtaBins0 = 50; 
	float EtaBins0[nEtaBins0+1] = {0.0};
	for(int i = 0; i <= nEtaBins0; ++i){
		EtaBins0[i] = -0.8 + 1.6*((float)i)*(0.02);
	}

	const int nEtaBins = 100; 
	float EtaBins[nEtaBins+1] = {0.0};
	for(int i = 0; i <= nEtaBins; ++i){
		EtaBins[i] = -2.0 + 2.0*((float)i)*(0.01);
	}

	const int nBinsNT = 50;
	float BinsNT[nBinsNT+1] = {0};

	for(int i = 0; i <= nBinsNT; ++i){
		BinsNT[i] = ((float)i)-0.5;
	}

	const int nNchBins = 100; 
	float NchBins[nNchBins+1] = {0.0};
	for(int i = 0; i <= nNchBins; ++i){
		NchBins[i] = -0.5 + (float)i;
	}

	//printf("==============================\n");
	//printf("Running the Fixed Code\n");
	//printf("Jetty cut: %f\n",fJettyCutOff);
	//printf("Isotropic cut: %f\n",fIsotrCutOff);
	//printf("==============================\n");

	const int ndEdxBins = fdEdxHigh-fdEdxLow;
	float dEdxBins[ndEdxBins+1];

	for(int i = fdEdxLow; i <= fdEdxHigh; ++i){
		dEdxBins[i-fdEdxLow] = i;
	}

	const char* ending[nHists] = {"02", "24", "46", "68"};

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

	hSOrvsV0M  = new TH2F("hSOrVsV0M","Measured SO vs V0M Per.;V0M Per.;#it{S}_{O} Reconstructed",100,0,100,1000,0,1);
	fListOfObjects->Add(hSOrvsV0M);

	hSOrvsTrks  = new TH2F("hSOrVsTrks","Measured SO vs Measured Ref. Mult. |#eta|<0.8;Reference mult. (|#eta|<0.8);#it{S}_{O} Reconstructed",100, 0, 100, 1000, 0, 1);
	fListOfObjects->Add(hSOrvsTrks);

	hMultPercvsNch = new TH3F("hMultPercvsNch",";V0M; CL1; #it{N}_{ch}",nBinsPer,BinsPer,nBinsPer,BinsPer,nNchBins,NchBins);
	fListOfObjects->Add(hMultPercvsNch);

	hSOtvsSOrV0M = new TH3F("hSOtvsSOrV0M",";#it{S}_{O}^{r};#it{S}_{O}^{t};V0M",nBinsSO,BinsSO,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hSOtvsSOrCL1 = new TH3F("hSOtvsSOrCL1",";#it{S}_{O}^{r};#it{S}_{O}^{t};CL1",nBinsSO,BinsSO,nBinsSO,BinsSO,nBinsPer,BinsPer);

	hPionSimvsSOV0M = new TH3F("hPionSimvsSOV0M",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hPionGenvsSOV0M = new TH3F("hPionGenvsSOV0M",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hKaonSimvsSOV0M = new TH3F("hKaonSimvsSOV0M",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hKaonGenvsSOV0M = new TH3F("hKaonGenvsSOV0M",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hProtonSimvsSOV0M = new TH3F("hProtonSimvsSOV0M",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hProtonGenvsSOV0M = new TH3F("hProtonGenvsSOV0M",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);

	hPionSimvsSOV0M2 = new TH3F("hPionSimvsSOV0M2",";#it{p}_{T}^{gen};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hKaonSimvsSOV0M2 = new TH3F("hKaonSimvsSOV0M2",";#it{p}_{T}^{gen};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hProtonSimvsSOV0M2 = new TH3F("hProtonSimvsSOV0M2",";#it{p}_{T}^{gen};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);

	hPionSimvsSOCL1 = new TH3F("hPionSimvsSOTrks",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hPionGenvsSOCL1 = new TH3F("hPionGenvsSOTrks",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hKaonSimvsSOCL1 = new TH3F("hKaonSimvsSOTrks",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hKaonGenvsSOCL1 = new TH3F("hKaonGenvsSOTrks",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hProtonSimvsSOCL1 =new TH3F("hProtonSimvsSOTrks",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hProtonGenvsSOCL1 =new TH3F("hProtonGenvsSOTrks",";#it{p}_{T};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);

	hPionSimvsSOCL12 = new TH3F("hPionSimvsSOTrks2",";#it{p}_{T}^{gen};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hKaonSimvsSOCL12 = new TH3F("hKaonSimvsSOTrks2",";#it{p}_{T}^{gen};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);
	hProtonSimvsSOCL12 = new TH3F("hProtonSimvsSOTrks2",";#it{p}_{T}^{gen};#it{S}_{O}^{r};V0M",nPtBins,ptBins,nBinsSO,BinsSO,nBinsPer,BinsPer);

	hTruthEtaSo = new TH1D("hTruthEtaSo","spherocity; #eta; counts",40,-1.0,1.0);
	hTruthPhiSo = new TH1D("hTruthPhiSo","spherocity; #phi; counts",64,0.0,2*TMath::Pi());

	hPhiSimV0M = new TH3F("hPhiSimV0M",";#phi;#eta;#it{S}_{O}",nPhiBins,PhiBins,nEtaBins,EtaBins,nBinsSO,BinsSO);
	hPhiGenV0M = new TH3F("hPhiGenV0M",";#phi;#eta;#it{S}_{O}",nPhiBins,PhiBins,nEtaBins,EtaBins,nBinsSO,BinsSO);

	hPhiSimCL1 = new TH3F("hPhiSimCL1",";#phi;#eta;#it{S}_{O}",nPhiBins,PhiBins,nEtaBins,EtaBins,nBinsSO,BinsSO);
	hPhiGenCL1 = new TH3F("hPhiGenCL1",";#phi;#eta;#it{S}_{O}",nPhiBins,PhiBins,nEtaBins,EtaBins,nBinsSO,BinsSO);

	hPhiStandard = new TH2F("hPhiStandard",";#eta ;#phi",nEtaBins0,EtaBins0,nPhiBins,PhiBins);
	fListOfObjects->Add(hPhiStandard);
	hPhiHybrid1 = new TH2F("hPhiHybrid",";#eta ;#phi",nEtaBins0,EtaBins0,nPhiBins,PhiBins);
	fListOfObjects->Add(hPhiHybrid1);
	hPhiTotal = new TH2F("hPhi_Standard_Hybrid",";#eta ;#phi",nEtaBins0,EtaBins0,nPhiBins,PhiBins);
	fListOfObjects->Add(hPhiTotal);

	hV0MvsSOvsNT = new TH3F("hV0MvsSOvsNT",";NT;#it{S}_{O};V0M",nBinsNT,BinsNT,nBinsSO0,BinsSO0,nBinsPer,BinsPer);
	fListOfObjects->Add(hV0MvsSOvsNT);
	hCL1vsSOvsNT = new TH3F("hTrksvsSOvsNT",";NT;#it{S}_{O};CL1",nBinsNT,BinsNT,nBinsSO0,BinsSO0,nBinsPer,BinsPer);
	fListOfObjects->Add(hCL1vsSOvsNT);
	hV0MvsSOvsNoNT = new TH3F("hV0MvsSOvsNoNT",";NT;#it{S}_{O};V0M",nBinsNT,BinsNT,nBinsSO0,BinsSO0,nBinsPer,BinsPer);
	fListOfObjects->Add(hV0MvsSOvsNoNT);
	hCL1vsSOvsNoNT = new TH3F("hTrksvsSOvsNoNT",";NT;#it{S}_{O};CL1",nBinsNT,BinsNT,nBinsSO0,BinsSO0,nBinsPer,BinsPer);
	fListOfObjects->Add(hCL1vsSOvsNoNT);

	hV0MvsNoSOvsNT = new TH3F("hV0MvsNoSOvsNT",";NT;#it{S}_{O};V0M",nBinsNT,BinsNT,nBinsSO0,BinsSO0,nBinsPer,BinsPer);
	fListOfObjects->Add(hV0MvsNoSOvsNT);
	hCL1vsNoSOvsNT = new TH3F("hTrksvsNoSOvsNT",";NT;#it{S}_{O};CL1",nBinsNT,BinsNT,nBinsSO0,BinsSO0,nBinsPer,BinsPer);
	fListOfObjects->Add(hCL1vsNoSOvsNT);
	hV0MvsNoSOvsNoNT = new TH3F("hV0MvsNoSOvsNoNT",";NT;#it{S}_{O};V0M",nBinsNT,BinsNT,nBinsSO0,BinsSO0,nBinsPer,BinsPer);
	fListOfObjects->Add(hV0MvsNoSOvsNoNT);
	hCL1vsNoSOvsNoNT = new TH3F("hTrksvsNoSOvsNoNT",";NT;#it{S}_{O};CL1",nBinsNT,BinsNT,nBinsSO0,BinsSO0,nBinsPer,BinsPer);
	fListOfObjects->Add(hCL1vsNoSOvsNoNT);

	hPhiToward = new TH1F("hPhiToward",";#phi;Entries",nPhiBins,PhiBins);
	fListOfObjects->Add(hPhiToward);
	hPhiAway = new TH1F("hPhiAway",";#phi;Entries",nPhiBins,PhiBins);
	fListOfObjects->Add(hPhiAway);
	hPhiTransverse = new TH1F("hPhiTransverse",";#phi;Entries",nPhiBins,PhiBins);
	fListOfObjects->Add(hPhiTransverse);

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


	else{
		fListOfObjects->Add(hSOtvsSOrV0M);
		fListOfObjects->Add(hSOtvsSOrCL1);
		fListOfObjects->Add(hPionSimvsSOV0M);
		fListOfObjects->Add(hPionSimvsSOV0M2);
		fListOfObjects->Add(hPionGenvsSOV0M);
		fListOfObjects->Add(hKaonSimvsSOV0M);
		fListOfObjects->Add(hKaonSimvsSOV0M2);
		fListOfObjects->Add(hKaonGenvsSOV0M);
		fListOfObjects->Add(hProtonSimvsSOV0M);
		fListOfObjects->Add(hProtonSimvsSOV0M2);
		fListOfObjects->Add(hProtonGenvsSOV0M);
		fListOfObjects->Add(hPionSimvsSOCL1);
		fListOfObjects->Add(hPionSimvsSOCL12);
		fListOfObjects->Add(hPionGenvsSOCL1);
		fListOfObjects->Add(hKaonSimvsSOCL1);
		fListOfObjects->Add(hKaonSimvsSOCL12);
		fListOfObjects->Add(hKaonGenvsSOCL1);
		fListOfObjects->Add(hProtonSimvsSOCL1);
		fListOfObjects->Add(hProtonSimvsSOCL12);
		fListOfObjects->Add(hProtonGenvsSOCL1);
		fListOfObjects->Add(hTruthEtaSo);
		fListOfObjects->Add(hTruthPhiSo);
		fListOfObjects->Add(hPhiSimV0M);
		fListOfObjects->Add(hPhiSimCL1);
		fListOfObjects->Add(hPhiGenV0M);
		fListOfObjects->Add(hPhiGenCL1);

	}

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
	//	int fnRefGlobal = -1;
	//	int IndxTrksMult = -1;
	//	int IndxV0MMult = -1;

	////fnRefGlobal = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8 );

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

	if( fESD->IsPileupFromSPDInMultBins() ){return;}
	fEvents->Fill(2.5,MultPercentile);

	if( fESD->IsIncompleteDAQ() ){return;}
	fEvents->Fill(3.5,MultPercentile);

	if( utils->IsSPDClusterVsTrackletBG(fESD) ){return;}
	fEvents->Fill(4.5,MultPercentile);

	if( !MultSelection->GetThisEventINELgtZERO() ){return;}
	fEvents->Fill(5.5,MultPercentile);

	if( !selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE) ){return;}
	fEvents->Fill(6.5,MultPercentile);

	if( !IsGoodZvertexPos(fESD) ){return;}
	fEvents->Fill(7.5,MultPercentile);

	float SOm = -1.0;
	SOm = GetSpherocity(hphiso, hetaso);

	float SOt = -1.0;
	if(fAnalysisMC)
		SOt = GetEventShapeTrue(fMCStack,hTruthPhiSo,hTruthEtaSo);

	GetLeadingObject();
	GetNT();

	if( (fRecLeadPt >= 5.0) && (SOm < 0.0) ){
		hV0MvsNoSOvsNT->Fill(fNT,SOm,V0MPercentile);
		hCL1vsNoSOvsNT->Fill(fNT,SOm,RefPercentile);
	}

	if( (fRecLeadPt < 5.0) && (SOm < 0.0) ){
		hV0MvsNoSOvsNoNT->Fill(fNT,SOm,V0MPercentile);
		hCL1vsNoSOvsNoNT->Fill(fNT,SOm,RefPercentile);
	}

	// Events with non-measured spherocity
	if( SOm < 0.0 ){
		fEvents->Fill(8.5,MultPercentile);
		return;
	}

	hSOrvsV0M->Fill(V0MPercentile,SOm);
	hSOrvsTrks->Fill(RefPercentile,SOm);

	fEvents->Fill(9.5,MultPercentile);

	ProduceArrayTrksESD(SOm);
	int Nch = -1;
	Nch = ReadESDEvent();
	hMultPercvsNch->Fill(V0MPercentile,RefPercentile,Nch);

	if( (fRecLeadPt >= 5.0) && (fRecLeadPt < 40.0) ){
		hV0MvsSOvsNT->Fill(fNT,SOm,V0MPercentile);
		hCL1vsSOvsNT->Fill(fNT,SOm,RefPercentile);
	}

	if( fRecLeadPt < 5.0){
		hV0MvsSOvsNoNT->Fill(fNT,SOm,V0MPercentile);
		hCL1vsSOvsNoNT->Fill(fNT,SOm,RefPercentile);
	}

	if(fAnalysisMC){

		hSOtvsSOrV0M->Fill(SOm,SOt,V0MPercentile);
		hSOtvsSOrCL1->Fill(SOm,SOt,RefPercentile);

		ProcessMCTruthV0M(V0MPercentile,SOt,SOm);
		ProcessMCTruthCL1(RefPercentile,SOt,SOm);

		//		int TruthMult = -1;
		//		TruthMult = GetMultiplicityParticles(0.8);

		AnalyseSimDataV0M(SOm,SOt);
		AnalyseSimDataCL1(SOm,SOt);

	}

	PostData(1, fListOfObjects);
}
//_____________________________________________________________________________
void AliAnalysisTaskSpherocity::AnalyseSimDataV0M(const float& SOm, const float& SOt)
{

	float MultPer  = -1.0;
	AliMultSelection *MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
	if(!MultSelection){ return; }
	if( MultSelection-> IsEventSelected() ){

		MultPer = MultSelection->GetMultiplicityPercentile("V0M",false);

	}

	if((MultPer < 0.0) || (MultPer > 100.0)) return;

	int nESDTracks = fESD->GetNumberOfTracks();
	for(int iT = 0; iT < nESDTracks; iT++){

		AliESDtrack* esdTrack = fESD->GetTrack(iT);
		if(!esdTrack)
			continue;

		if( TMath::Abs(esdTrack->Eta()) > fEtaCut )
			continue;

		if(esdTrack->Pt() < fPtMinCut)
			continue;

		if(esdTrack->GetTPCsignalN() < fNcl)
			continue;

		if(!PhiCut(esdTrack->Pt(), esdTrack->Phi(), esdTrack->Charge(), MAGF, fcutLow, fcutHigh))
			continue;

		if(!fTrackFilterGolden->AcceptTrack(esdTrack))
			continue;

		float dcaxy = 0.0;
		float dcaz = 0.0;
		esdTrack->GetImpactParameters(dcaxy,dcaz);

		if(TMath::Abs(dcaxy)>GetMaxDCApTDep(fcutDCAxy,esdTrack->Pt()))
			continue;

		const int label = TMath::Abs(esdTrack->GetLabel());
		TParticle* mcTrack = 0;
		mcTrack = fMCStack->Particle(label);

		if(!mcTrack) 
			continue;

		if((esdTrack->Charge()==0) || (TMath::Abs(esdTrack->Charge()) < 0.1))
			continue;

		int pdgCode = mcTrack->GetPdgCode();
		int pidCode = GetPidCode(pdgCode);

		if( fMCStack->IsPhysicalPrimary(label) ){

			if(pidCode==1){ 
				hPionSimvsSOV0M->Fill(esdTrack->Pt(),SOm,MultPer);	
			}
			if(pidCode==2){ 
				hKaonSimvsSOV0M->Fill(esdTrack->Pt(),SOm,MultPer);	
			}
			if(pidCode==3){
				hProtonSimvsSOV0M->Fill(esdTrack->Pt(),SOm,MultPer);	
			}

			// Second loop 
			for(int jT = 0; jT < nESDTracks; jT++){

				if(iT==jT)
					continue;

				AliESDtrack* jTrack = fESD->GetTrack(jT);
				if(!esdTrack)
					continue;

				if( TMath::Abs(jTrack->Eta()) > fEtaCut )
					continue;

				if(jTrack->Pt() < fPtMinCut)
					continue;

				if(!fTrackFilterGolden->AcceptTrack(jTrack))
					continue;

				float dcaxy = 0.0;
				float dcaz = 0.0;
				jTrack->GetImpactParameters(dcaxy,dcaz);

				if(TMath::Abs(dcaxy)>GetMaxDCApTDep(fcutDCAxy,jTrack->Pt()))
					continue;

				const int label = TMath::Abs(jTrack->GetLabel());
				TParticle* j_mcTrack = 0;
				j_mcTrack = fMCStack->Particle(label);

				if(!j_mcTrack) 
					continue;

				if((jTrack->Charge()==0) || (TMath::Abs(jTrack->Charge()) < 0.1))
					continue;

				double dPhisim = -999.0;
				double dPhigen = -999.0;
				double dEtasim = -999.0;
				double dEtagen = -999.0;

				dPhisim = DeltaPhi(esdTrack->Phi(),jTrack->Phi());
				dPhigen = DeltaPhi(mcTrack->Phi(),j_mcTrack->Phi()); 

				dEtasim = esdTrack->Eta()-jTrack->Eta();
				dEtagen = mcTrack->Eta()-j_mcTrack->Eta(); 

				if( fMCStack->IsPhysicalPrimary(label) ){

					if(MultPer > 10.0) continue;

					hPhiSimV0M->Fill(dPhisim,dEtasim,SOm);
					hPhiGenV0M->Fill(dPhigen,dEtagen,SOt);

				}
			}
		}
	}
}
//_____________________________________________________________________________
void AliAnalysisTaskSpherocity::AnalyseSimDataCL1(const float& SOm, const float& SOt)
{

	float MultPer  = -1.0;
	AliMultSelection *MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
	if(!MultSelection){ return; }
	if( MultSelection-> IsEventSelected() ){

		MultPer = MultSelection->GetMultiplicityPercentile("RefMult08",false);

	}

	if((MultPer < 0.0) || (MultPer > 100.0)) return;

	int nESDTracks = fESD->GetNumberOfTracks();
	for(int iT = 0; iT < nESDTracks; iT++){

		AliESDtrack* esdTrack = fESD->GetTrack(iT);
		if(!esdTrack)
			continue;

		if( TMath::Abs(esdTrack->Eta()) > fEtaCut )
			continue;

		if(esdTrack->Pt() < fPtMinCut)
			continue;

		if(esdTrack->GetTPCsignalN() < fNcl)
			continue;

		if(!PhiCut(esdTrack->Pt(), esdTrack->Phi(), esdTrack->Charge(), MAGF, fcutLow, fcutHigh))
			continue;

		if(!fTrackFilterGolden->AcceptTrack(esdTrack))
			continue;

		float dcaxy = 0.0;
		float dcaz = 0.0;
		esdTrack->GetImpactParameters(dcaxy,dcaz);

		if(TMath::Abs(dcaxy)>GetMaxDCApTDep(fcutDCAxy,esdTrack->Pt()))
			continue;

		const int label = TMath::Abs(esdTrack->GetLabel());
		TParticle* mcTrack = 0;
		mcTrack = fMCStack->Particle(label);

		if(!mcTrack) 
			continue;

		if((esdTrack->Charge()==0) || (TMath::Abs(esdTrack->Charge()) < 0.1))
			continue;

		int pdgCode = mcTrack->GetPdgCode();
		int pidCode = GetPidCode(pdgCode);

		if( fMCStack->IsPhysicalPrimary(label) ){

			if(pidCode==1){ 
				hPionSimvsSOCL1->Fill(esdTrack->Pt(),SOm,MultPer);	
			}
			if(pidCode==2){ 
				hKaonSimvsSOCL1->Fill(esdTrack->Pt(),SOm,MultPer);	
			}
			if(pidCode==3){
				hProtonSimvsSOCL1->Fill(esdTrack->Pt(),SOm,MultPer);	
			}

			// Second loop 
			for(int jT = 0; jT < nESDTracks; jT++){

				if(iT==jT)
					continue;

				AliESDtrack* jTrack = fESD->GetTrack(jT);
				if(!esdTrack)
					continue;

				if( TMath::Abs(jTrack->Eta()) > fEtaCut )
					continue;

				if(jTrack->Pt() < fPtMinCut)
					continue;

				if(!fTrackFilterGolden->AcceptTrack(jTrack))
					continue;

				float dcaxy = 0.0;
				float dcaz = 0.0;
				jTrack->GetImpactParameters(dcaxy,dcaz);

				if(TMath::Abs(dcaxy)>GetMaxDCApTDep(fcutDCAxy,jTrack->Pt()))
					continue;

				const int label = TMath::Abs(jTrack->GetLabel());
				TParticle* j_mcTrack = 0;
				j_mcTrack = fMCStack->Particle(label);

				if(!j_mcTrack) 
					continue;

				if((jTrack->Charge()==0) || (TMath::Abs(jTrack->Charge()) < 0.1))
					continue;

				double dPhisim = -999.0;
				double dPhigen = -999.0;
				double dEtasim = -999.0;
				double dEtagen = -999.0;

				dPhisim = DeltaPhi(esdTrack->Phi(),jTrack->Phi());
				dPhigen = DeltaPhi(mcTrack->Phi(),j_mcTrack->Phi()); 

				dEtasim = esdTrack->Eta()-jTrack->Eta();
				dEtagen = mcTrack->Eta()-j_mcTrack->Eta(); 

				if( fMCStack->IsPhysicalPrimary(label) ){

					if(MultPer > 10.0) continue;

					hPhiSimCL1->Fill(dPhisim,dEtasim,SOm);
					hPhiGenCL1->Fill(dPhigen,dEtagen,SOt);

				}
			}
		}
	}
}
//_____________________________________________________________________________
int AliAnalysisTaskSpherocity::GetMultiplicityParticles(Double_t etaCut)
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
void AliAnalysisTaskSpherocity::ProcessMCTruthV0M(const float& MultPer, const float& SOt, const float& SOm)
{
	// Fill the special MC histogram with the MC truth info

	const int nTracksMC = fMCStack->GetNtrack();

	for (int iTracks = 0; iTracks < nTracksMC; iTracks++) {

		TParticle* trackMC = fMCStack->Particle(iTracks);
		if( !trackMC )
			continue;

		if( !(fMCStack->IsPhysicalPrimary(iTracks)) )
			continue;

		TParticlePDG* pdgPart = trackMC ->GetPDG();
		double chargeMC = pdgPart->Charge();

		if((chargeMC==0) || TMath::Abs(chargeMC) < 0.1)
			continue;

		if( TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		if(trackMC->Pt() < fPtMinCut)
			continue;

		int pdgCode = trackMC->GetPdgCode();
		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		if(pidCodeMC==1){

			hPionSimvsSOV0M2->Fill(trackMC->Pt(),SOm,MultPer);	
			hPionGenvsSOV0M->Fill(trackMC->Pt(),SOt,MultPer);	

		}else if(pidCodeMC==2){

			hKaonSimvsSOV0M2->Fill(trackMC->Pt(),SOm,MultPer);	
			hKaonGenvsSOV0M->Fill(trackMC->Pt(),SOt,MultPer);	

		}else if(pidCodeMC==3){

			hProtonSimvsSOV0M2->Fill(trackMC->Pt(),SOm,MultPer);	
			hProtonGenvsSOV0M->Fill(trackMC->Pt(),SOt,MultPer);	

		}else{ continue; }

	}//MC track loop
}
//_____________________________________________________________________________
void AliAnalysisTaskSpherocity::ProcessMCTruthCL1(const float& MultPer, const float& SOt, const float& SOm)
{
	// Fill the special MC histogram with the MC truth info

	const int nTracksMC = fMCStack->GetNtrack();

	for (int iTracks = 0; iTracks < nTracksMC; iTracks++) {

		TParticle* trackMC = fMCStack->Particle(iTracks);
		if( !trackMC )
			continue;

		if( !(fMCStack->IsPhysicalPrimary(iTracks)) )
			continue;

		TParticlePDG* pdgPart = trackMC ->GetPDG();
		double chargeMC = pdgPart->Charge();

		if((chargeMC==0) || TMath::Abs(chargeMC) < 0.1)
			continue;

		if(TMath::Abs(trackMC->Eta()) > fEtaCut )
			continue;

		if(trackMC->Pt() < fPtMinCut)
			continue;

		int pdgCode = trackMC->GetPdgCode();
		short pidCodeMC = 0;
		pidCodeMC = GetPidCode(pdgCode);

		if(pidCodeMC==1){

			hPionSimvsSOCL12->Fill(trackMC->Pt(),SOm,MultPer);	
			hPionGenvsSOCL1->Fill(trackMC->Pt(),SOt,MultPer);	

		}else if(pidCodeMC==2){

			hKaonSimvsSOCL12->Fill(trackMC->Pt(),SOm,MultPer);	
			hKaonGenvsSOCL1->Fill(trackMC->Pt(),SOt,MultPer);	

		}else if(pidCodeMC==3){

			hProtonSimvsSOCL12->Fill(trackMC->Pt(),SOm,MultPer);	
			hProtonGenvsSOCL1->Fill(trackMC->Pt(),SOt,MultPer);	

		}else{ continue; }

	}//MC track loop
}
//____________________________________________________________________
/*
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
*/
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
float AliAnalysisTaskSpherocity::GetMaxDCApTDep( TF1 *fMaxDCAxy, Double_t ptI)
{

	Double_t maxDCAxy = 10;
	maxDCAxy = fMaxDCAxy->Eval(ptI);
	return maxDCAxy;

}
//________________________________________________________________________
float AliAnalysisTaskSpherocity::GetEventShapeTrue( AliStack *event, TH1D *hphi, TH1D *heta )
{

	double lreturnval = 1.0;

	fMCStack = event;

	lreturnval = GetSpherocityMC( hphi, heta );

	return lreturnval;
}
//________________________________________________________________________
float AliAnalysisTaskSpherocity::GetSpherocityMC( TH1D * hphi, TH1D *heta )
{
	vector<Float_t> pt;
	vector<Float_t> eta;
	vector<Float_t> phi;

	fNrec = ReadMC(pt, eta, phi, hphi, heta);
	if( fNrec < fMinMult )
		return -0.5;

	float spherocity = AnalyseGetSpherocity( pt, eta, phi );

	return spherocity;
}
//________________________________________________________________________
int AliAnalysisTaskSpherocity::ReadMC( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, vector<Float_t> &phiArray, TH1D * hphi, TH1D *heta )
{

	ptArray.clear();
	etaArray.clear();
	phiArray.clear();

	int nTracks = fMCStack->GetNtrack();
	int nTrue = 0;

	for(int iT = 0; iT < nTracks; iT++) {
		//Cuts
		if(!(fMCStack->IsPhysicalPrimary(iT)))
			continue;

		TParticle* Track = fMCStack->Particle(iT);
		if(!Track)
			continue;

		float eta = Track->Eta();
		float pt  = Track->Pt();
		float phi = Track->Phi();

		TParticlePDG* pdgPart = Track->GetPDG();
		double chargeMC = pdgPart->Charge();

		if( TMath::Abs(chargeMC) < 0.1 )
			continue;

		if( TMath::Abs(eta) > fEtaCut )
			continue;

		//cuts in pt
		if( pt > fPtMaxCut || pt < fPtMinCut )
			continue;

		if(hphi)
			hphi->Fill(phi);
		if(heta)
			heta->Fill(eta);

		ptArray.push_back(pt);
		etaArray.push_back(eta);
		phiArray.push_back(phi);

		nTrue++;

	} //close first loop on nTracks

	return nTrue;

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
		if( pt > fPtMaxCut || pt < fPtMinCut )
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
//_____________________________________________________________________
int AliAnalysisTaskSpherocity::ReadESDEvent(){

	//	ptArray.clear();
	//	etaArray.clear();
	//	phiArray.clear();

	int nTracks = fESD->GetNumberOfTracks();
	int nRec    = 0;

	for(int iT = 0; iT < nTracks; iT++) {
		AliESDtrack* Track = 0;
		Track = fESD->GetTrack(iT);
		if(!Track)
			continue;

		float eta  = Track->Eta();
		float pt   = Track->Pt();
		//		float phi  = Track->Phi();

		if( !(TMath::Abs(eta) < fEtaCut) )
			continue;

		//cuts in pt
		if( pt > fPtMaxCut || pt < fPtMinCut )
			continue;

		//quality cuts
		if(!fTrackFilter->IsSelected(Track))
			continue;

		//		if(hphi)
		//			hphi->Fill(phi);
		//		if(heta)
		//			heta->Fill(eta);

		//		ptArray.push_back(pt);
		//		etaArray.push_back(eta);
		//		phiArray.push_back(phi);

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
	}else if(strcmp(fPeriod,"17data")==0){
		aPos = 49.6097; bPos = 0.922856; cPos = -6.57484; dPos = 65.3117; ePos = -372.142; fPos = 950.451; gPos = -1085.27; hPos = 458.144;
		aNeg = 49.6555; bNeg = 6.98696; cNeg = 102.734;  dNeg = 566.424; eNeg = 1513.64;  fNeg = 2092.01; gNeg = 1429.32;  hNeg = 375.642;
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
void AliAnalysisTaskSpherocity::GetNT()
{
	const double pi = TMath::Pi();
	int multTS = 0;

	int iTracks(fESD->GetNumberOfTracks());
	for(int i = 0; i < iTracks; i++){

		if(i==fRecLeadIn) continue;
		AliESDtrack* esdtrack = static_cast<AliESDtrack*>(fESD->GetTrack(i));
		if(!esdtrack) continue;
		if(TMath::Abs(esdtrack->Eta()) > fEtaCut) continue;
		if(esdtrack->Pt() < fPtMinCut) continue;

		AliESDtrack* track = 0x0;
		track = SetHybridTrackCuts(esdtrack,kTRUE,kTRUE);
		if(!track) continue;

		hPhiTotal->Fill(track->Eta(),track->Phi());
		double DPhi = DeltaPhi(track->Phi(), fRecLeadPhi);

		if(TMath::Abs(DPhi)<pi/3.0){
			hPhiToward->Fill(DPhi);
		}
		else if(TMath::Abs(DPhi-pi)<pi/3.0){
			hPhiAway->Fill(DPhi);
		}
		else{
			hPhiTransverse->Fill(DPhi);
			multTS++;
		}

		delete track;

	}

	fNT = multTS;	

}
//________________________________________________________________________
double AliAnalysisTaskSpherocity::DeltaPhi(double phi, double Lphi, double rangeMin, double rangeMax)
{

	double dphi = -999;
	double pi = TMath::Pi();
	//      if(Lphi > 2*pi || Lphi < 0)cout << "Lphi :: " << Lphi << endl;
	//      if(phi  > 2*pi || phi < 0)cout << "phi = " << phi << endl;

	if(phi < 0)          phi += 2*pi;
	else if(phi > 2*pi)  phi -= 2*pi;
	if(Lphi < 0)         Lphi += 2*pi;
	else if(Lphi > 2*pi) Lphi -= 2*pi;
	dphi = Lphi - phi;
	if (dphi < rangeMin)      dphi += 2*pi;
	else if (dphi > rangeMax) dphi -= 2*pi;

	return dphi;
}
//________________________________________________________________________
void AliAnalysisTaskSpherocity::GetLeadingObject()
{

	double flPt = 0.0;
	double flPhi = 0.0;
	int flIndex = 0;

	int iTracks(fESD->GetNumberOfTracks());
	for(int i=0; i < iTracks; i++) {

		AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
		if(!track) continue;
		if(TMath::Abs(track->Eta()) > fEtaCut) continue;
		if(track->Pt() < fPtMinCut) continue;

		AliESDtrack* track_hybrid = 0x0;
		track_hybrid = SetHybridTrackCuts(track,kFALSE,kFALSE);
		if(!track_hybrid) continue;

		if(!fGeometricalCut->AcceptTrack(track_hybrid)) continue;

		if(flPt < track_hybrid->Pt()){
			flPt  = track_hybrid->Pt();
			flPhi = track_hybrid->Phi();
			flIndex = i;
		}

		delete track_hybrid;

	}

	fRecLeadPhi = flPhi;
	fRecLeadPt  = flPt;
	fRecLeadIn  = flIndex;

}
//________________________________________________________________________
AliESDtrack* AliAnalysisTaskSpherocity::SetHybridTrackCuts(AliESDtrack *esdtrack, const bool fillPhiStand, const bool fillPhHyb1)
{

	AliESDtrack *newTrack = 0x0;

	if(fTrackFilter->IsSelected(esdtrack)){
		newTrack = new AliESDtrack(*esdtrack);
		if(fillPhiStand) hPhiStandard->Fill(newTrack->Eta(),newTrack->Phi());
		////                    newTrack->SetTRDQuality(0);
	}
	else if(fHybridTrackCuts1->AcceptTrack(esdtrack)){
		if(esdtrack->GetConstrainedParam()){
			newTrack = new AliESDtrack(*esdtrack);
			const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
			newTrack->Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
			////                            newTrack->SetTRDQuality(1);
			if(fillPhHyb1) hPhiHybrid1->Fill(newTrack->Eta(),newTrack->Phi());
		}
		else{ return 0x0; }
	}

	else{
		return 0x0;
	}

	return newTrack;

}
