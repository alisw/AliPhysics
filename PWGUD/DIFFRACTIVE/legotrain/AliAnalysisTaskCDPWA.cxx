/*************************************************************************
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
 **************************************************************************/
//
// Select events according to gap conditions, analyze two track events in pp
// collisions
//
// Author:
//  Xianguo Lu <lu@physi.uni-heidelberg.de>
//  continued by
//  Felix Reidt <Felix.Reidt@cern.ch>
//  continued by
//  Taesoo Kim <taesoo.kim@cern.ch>

#include <TArrayI.h>
#include <TList.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TParticle.h>

#include <TMath.h>

#include "AliVEvent.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliRawEventHeaderBase.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliESDVZERO.h"
#include "AliESDAD.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliStack.h"
#include "TGeoManager.h"

#include "AliMultiplicitySelectionCPPWA.h"
#include "AliAnalysisTaskCDPWA.h"

ClassImp(AliAnalysisTaskCDPWA);
ClassImp(AliAnalysisTaskCDPWA::EventInfo);
ClassImp(AliAnalysisTaskCDPWA::CombInfo);
ClassImp(AliAnalysisTaskCDPWA::TrackInfo);

//______________________________________________________________________________
void AliAnalysisTaskCDPWA::EventInfo::Fill(const AliESDEvent *esdEvent) 
{
	const AliESDHeader *esdHeader = esdEvent->GetHeader();
	if (NULL == esdHeader) // this is already dealt with in UserExec
		return;

	fRunNumber = esdEvent->GetRunNumber();
	fTimeStamp = esdHeader->GetTimeStamp();
	//My own period
	if (fRunNumber < 117224) fPeriod = 1;//10b
	else if (fRunNumber < 121041) fPeriod = 2;//10c
	else if (fRunNumber < 126438) fPeriod = 3;//10d
	else if (fRunNumber < 130851) fPeriod = 4;//10e
	else if (fRunNumber < 135032) fPeriod = 5;//10f
	else if (fRunNumber < 224891) fPeriod = 6;//12b (8 TeV)
	else if (fRunNumber < 227750) fPeriod = 7;//15f
	else if (fRunNumber < 232465) fPeriod = 8;//15g
	else if (fRunNumber < 235196) fPeriod = 9;//15h
	else if (fRunNumber < 236892) fPeriod = 10;//15i
	else if (fRunNumber < 238682) fPeriod = 11;//15j
	else if (fRunNumber < 239207) fPeriod = 12;//15k
	else if (fRunNumber < 243374) fPeriod = 13;//15l
	else if (fRunNumber < 244340) fPeriod = 14;//15m
	else if (fRunNumber < 244824) fPeriod = 15;//15n
	else if (fRunNumber <= 246994) fPeriod = 16;//15o

	return;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::CombInfo::Fill(const AliESDEvent *esdEvent) 
{
	if (NULL == esdEvent) return;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::TrackInfo::Fill(AliESDtrack *tr, AliPIDResponse *pidResponse, AliPIDCombined *pidCombined) 
{
	if (NULL == tr || NULL == pidResponse) {
		AliError(Form("tr=%p pidResponse=%p", tr, pidResponse));
		return;
	}
	fSign = tr->GetSign();
	fPx   = tr->Px();
	fPy   = tr->Py();  
	fPz   = tr->Pz();  
	fEnergy = tr->E();
	fIntegratedLength = tr->GetIntegratedLength();
	fITSSignal = tr->GetITSsignal();
	fTPCSignal = tr->GetTPCsignal();
	fTOFSignal = tr->GetTOFsignal();
	fTRDSignal = tr->GetTRDsignal();

	fPIDStatus[0] = pidResponse->CheckPIDStatus(AliPIDResponse::kITS, tr);
	fPIDStatus[1] = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, tr);
	fPIDStatus[2] = pidResponse->CheckPIDStatus(AliPIDResponse::kTOF, tr);
	fPIDStatus[3] = pidResponse->CheckPIDStatus(AliPIDResponse::kTRD, tr);

	for (Int_t i=0; i<AliPID::kSPECIES; ++i) {
		const AliPID::EParticleType particleType(static_cast<const AliPID::EParticleType>(i));
		fITSSigma[i] = pidResponse->NumberOfSigmasITS(tr, particleType);
		fTPCSigma[i] = pidResponse->NumberOfSigmasTPC(tr, particleType);
		fTOFSigma[i] = pidResponse->NumberOfSigmasTOF(tr, particleType);
		fITSSignalDelta[i] = pidResponse->GetSignalDelta(AliPIDResponse::kITS,tr,particleType);
	}
	pidResponse->ComputePIDProbability(AliPIDResponse::kTPC,tr,AliPID::kSPECIES,fTPCProb);
	pidResponse->ComputePIDProbability(AliPIDResponse::kTOF,tr,AliPID::kSPECIES,fTOFProb);
	pidResponse->ComputePIDProbability(AliPIDResponse::kITS,tr,AliPID::kSPECIES,fITSProb);

	pidCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
	fDetMask[0] = pidCombined->ComputeProbabilities(tr,pidResponse,fTPCBayesProb);
	pidCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
	fDetMask[1] = pidCombined->ComputeProbabilities(tr,pidResponse,fTOFBayesProb);
	pidCombined->SetDetectorMask(AliPIDResponse::kDetITS);
	fDetMask[2] = pidCombined->ComputeProbabilities(tr,pidResponse,fITSBayesProb);
	pidCombined->SetDetectorMask(AliPIDResponse::kDetTRD);
	fDetMask[3] = pidCombined->ComputeProbabilities(tr,pidResponse,fTRDBayesProb);
	pidCombined->SetDetectorMask(AliPIDResponse::kDetTPC | 
			AliPIDResponse::kDetTOF |
		       	AliPIDResponse::kDetITS |
			AliPIDResponse::kDetTRD);
	fDetMask[4] = pidCombined->ComputeProbabilities(tr,pidResponse,ftotBayesProb);
}
//______________________________________________________________________________
AliAnalysisTaskCDPWA::AliAnalysisTaskCDPWA(const char* name):
	AliAnalysisTaskSE(name)
	, fIsRun2(0)
	, fIsMC(0)
	, fCombmode(0)
	, fIsPythia8(0)
	, fIsSys(0)
	, fIsPWAMC(0)
	, fnSys(0)
	, fIsPythia(0)
	, fIsPhojet(0)
	, fIsEPOS(0)
	, fTree(NULL)
	, fTree_Comb(NULL)
	, fList(NULL)
	, fMCTrack("TLorentzVector",4)
	, fTrackInfo("AliAnalysisTaskCDPWA::TrackInfo",4*11)
	, fEventInfo()
	, fCombInfo()
	, fTriggerAnalysis()
	, fAnalysisUtils()
	, fPIDResponse(0x0)
	, fPIDCombined(0x0)
	, fTrackCuts(0x0)
	, fESDEvent(0x0)
	, fMult(0x0)
	, fRunNumber(0)
	//Hist
	, fHistEvent(0x0)
	, fHistTrigger(0x0)
	, fHistEventProcesses(0x0)
	, fHistPrimVtxX(0x0)
	, fHistPrimVtxY(0x0)
	, fHistPrimVtxZ(0x0)
	, fMultRegionsMC(0x0)
	, fSPDFiredChipsClvsHd(0x0)
	, fSPDwc(0x0)
	, fRunVsMBOR_V0(0x0)
	, fRunVsMBAND_V0(0x0)
	, fRunVsMBOR_AD(0x0)
	, fRunVsMBAND_AD(0x0)
	, fRunVsMBOR_Global(0x0)
	, fRunVsMBAND_Global(0x0)
	, fRunVs2t(0x0)
	, fRunVs2t_ITSSA(0x0)
	, fRunVs4t(0x0)
	, fRunVs4t_ITSSA(0x0)
	, fMultNG_ST(0x0)
	, fMultNG_MS(0x0)
	, fMultDG_ST(0x0)
	, fMultDG_MS(0x0)
	, fMassNG_ST_2t(0x0)
	, fMassNG_MS_2t(0x0)
	, fMassDG_ST_2t(0x0)
	, fMassDG_MS_2t(0x0)
	, fMassNG_ST_4t(0x0)
	, fMassNG_MS_4t(0x0)
	, fMassDG_ST_4t(0x0)
	, fMassDG_MS_4t(0x0)
	, fTrackCutsInfo(0x0)
	, fTrackCutsInfo_ITSSA(0x0)
	, fhClusterVsTracklets_bf(0x0)
	, fhClusterVsTracklets_af(0x0)
	, fMult_ITSSA(0x0)
	, fMult_DG_ITSSA(0x0)
	, fMult_NG_ITSSA(0x0)
	, fMult_ITSSA_MS(0x0)
	, fMult_DG_ITSSA_MS(0x0)
	, fMult_NG_ITSSA_MS(0x0)
	, fTPCSignal(0x0)
	, fTOFSignal(0x0)
	, fITSSignal(0x0)
	, fTRDSignal(0x0)
	, fMult_Gen(0x0)
	, fMult_Gen_Process(0x0)
	, fMult_Rec_DG_Process(0x0)
	, fMult_Rec_NG_Process(0x0)
{
	for (Int_t i = 0; i < 10; i++) {
		if (i < 2) {
			fV0Time[i] = 0x0;
			fADTime[i] = 0x0;
		}
		if (i < 6) fRunFiducial[i] = 0x0;
		fRunVsDG[i] = 0x0;
	}
	//
	// standard constructor (the one which should be used)
	//
	// slot in TaskSE must start from 1
	DefineOutput(1, TTree::Class());
	DefineOutput(2, TTree::Class());
	DefineOutput(3, TList::Class());
}
//______________________________________________________________________________
AliAnalysisTaskCDPWA::AliAnalysisTaskCDPWA():
	AliAnalysisTaskSE()
	, fIsRun2(0)
	, fIsMC(0)
	, fCombmode(0)
	, fIsPythia8(0)
	, fIsSys(0)
	, fIsPWAMC(0)
	, fnSys(0)
	, fIsPythia(0)
	, fIsPhojet(0)
	, fIsEPOS(0)
	, fTree(NULL)
	, fTree_Comb(NULL)
	, fList(NULL)
	, fMCTrack("TLorentzVector",4)
	, fTrackInfo("AliAnalysisTaskCDPWA::TrackInfo",4*11)
	, fEventInfo()
	, fCombInfo()
	, fTriggerAnalysis()
	, fAnalysisUtils()
	, fPIDResponse(0x0)
	, fPIDCombined(0x0)
	, fTrackCuts(0x0)
	, fESDEvent(0x0)
	, fMult(0x0)
	, fRunNumber(0)
	//Hist
	, fHistEvent(0x0)
	, fHistTrigger(0x0)
	, fHistEventProcesses(0x0)
	, fHistPrimVtxX(0x0)
	, fHistPrimVtxY(0x0)
	, fHistPrimVtxZ(0x0)
	, fMultRegionsMC(0x0)
	, fSPDFiredChipsClvsHd(0x0)
	, fSPDwc(0x0)
	, fRunVsMBOR_V0(0x0)
	, fRunVsMBAND_V0(0x0)
	, fRunVsMBOR_AD(0x0)
	, fRunVsMBAND_AD(0x0)
	, fRunVsMBOR_Global(0x0)
	, fRunVsMBAND_Global(0x0)
	, fRunVs2t(0x0)
	, fRunVs2t_ITSSA(0x0)
	, fRunVs4t(0x0)
	, fRunVs4t_ITSSA(0x0)
	, fMultNG_ST(0x0)
	, fMultNG_MS(0x0)
	, fMultDG_ST(0x0)
	, fMultDG_MS(0x0)
	, fMassNG_ST_2t(0x0)
	, fMassNG_MS_2t(0x0)
	, fMassDG_ST_2t(0x0)
	, fMassDG_MS_2t(0x0)
	, fMassNG_ST_4t(0x0)
	, fMassNG_MS_4t(0x0)
	, fMassDG_ST_4t(0x0)
	, fMassDG_MS_4t(0x0)
	, fTrackCutsInfo(0x0)
	, fTrackCutsInfo_ITSSA(0x0)
	, fhClusterVsTracklets_bf(0x0)
	, fhClusterVsTracklets_af(0x0)
	, fMult_ITSSA(0x0)
	, fMult_DG_ITSSA(0x0)
	, fMult_NG_ITSSA(0x0)
	, fMult_ITSSA_MS(0x0)
	, fMult_DG_ITSSA_MS(0x0)
	, fMult_NG_ITSSA_MS(0x0)
	, fTPCSignal(0x0)
	, fTOFSignal(0x0)
	, fITSSignal(0x0)
	, fTRDSignal(0x0)
	, fMult_Gen(0x0)
	, fMult_Gen_Process(0x0)
	, fMult_Rec_DG_Process(0x0)
	, fMult_Rec_NG_Process(0x0)
{
	for (Int_t i = 0; i < 10; i++) {
		if (i < 2) {
			fV0Time[i] = 0x0;
			fADTime[i] = 0x0;
		}
		if (i < 6) fRunFiducial[i] = 0x0;
		fRunVsDG[i] = 0x0;
	}
}
//______________________________________________________________________________
AliAnalysisTaskCDPWA::~AliAnalysisTaskCDPWA()
{
	//Destructor(pointer should be deleted)
	if (fTree) {delete fTree; fTree=0x0;}
	if (fTree_Comb) {delete fTree_Comb; fTree_Comb=0x0;}
	if (fList) {delete fList; fList=0x0;}
	if (fTrackCuts) {delete fTrackCuts; fTrackCuts=0x0;}
	if (fPIDResponse) {delete fPIDResponse; fPIDResponse=0x0;}
	if (fPIDCombined) {delete fPIDCombined; fPIDCombined=0x0;}
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::SetBranchPWA(TTree *t) {
	t->Branch("fEventInfo",&fEventInfo);
	t->Branch("fTrackInfo",&fTrackInfo);
	if (fIsPWAMC) {
		t->Branch("fMCTrack",&fMCTrack);
	}
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::SetBranchComb(TTree *t) {
	t->Branch("fCombInfo",&fCombInfo);
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::UserCreateOutputObjects()
{
	// TTree + TList should be defined-----------------------------------------
	OpenFile(1);
	fTree = new TTree("tree1","PWAtree");
	SetBranchPWA(fTree);

	OpenFile(2);
	fTree_Comb = new TTree("tree_Comb","tree_Comb");
	if (fCombmode) SetBranchComb(fTree_Comb);

	OpenFile(3);
	fList = new TList();
	fList->SetOwner();
	Int_t minRn = 114500;//min run number in Run1
	Int_t maxRn = 179300;//max run number in Run1
	if (fIsRun2) {//Run2
		minRn = 224800;
		maxRn = 260000;
	}
	Int_t diff = maxRn - minRn;
	if (!fHistEvent) {
		fHistEvent = new TH1D("fHistEvent","Number of Events for each selection",kAll,0,kAll);
		fHistEvent->GetXaxis()->SetBinLabel(kInput+1,"InputEvent");
		fHistEvent->GetXaxis()->SetBinLabel(kMCCheck+1,"MCCheck");
		fHistEvent->GetXaxis()->SetBinLabel(kOnlineTrigger+1,"OnlineTrigger");
		fHistEvent->GetXaxis()->SetBinLabel(kOfflineCut+1,"OfflineMBOR");
		fHistEvent->GetXaxis()->SetBinLabel(kVtxCut+1,"VertexCut");
		fHistEvent->GetXaxis()->SetBinLabel(kPileUpCut+1,"PileUpCut");
		fHistEvent->GetXaxis()->SetBinLabel(kClusterCut+1,"ClusterCut");
		fHistEvent->GetXaxis()->SetBinLabel(kMBOR_V0+1,"MBOR_V0");
		fHistEvent->GetXaxis()->SetBinLabel(kMBAND_V0+1,"MBAND_V0");
		fHistEvent->GetXaxis()->SetBinLabel(kMBOR_AD+1,"MBOR_AD");
		fHistEvent->GetXaxis()->SetBinLabel(kMBAND_AD+1,"MBAND_AD");
		fHistEvent->GetXaxis()->SetBinLabel(kMBOR_Global+1,"MBOR_Global");
		fHistEvent->GetXaxis()->SetBinLabel(kMBAND_Global+1,"MBAND_Global");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0SPD+1,"DG_V0SPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGADSPD+1,"DG_ADSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADSPD+1,"DG_V0ADSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADFMDSPD+1,"DG_V0ADFMDSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADFMDZDCSPD+1,"DG_V0ADFMDZDCSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0FMDSPD+1,"DG_V0FMDSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0FMDZDCSPD+1,"DG_V0FMDZDCSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kNG_Data+1,"No-gap");
		fHistEvent->GetXaxis()->SetBinLabel(kSGA_Data+1,"Single-gap A");
		fHistEvent->GetXaxis()->SetBinLabel(kSGC_Data+1,"Single-gap C");
		fHistEvent->GetXaxis()->SetBinLabel(k2Tracks+1,"2Tracks");
		fHistEvent->GetXaxis()->SetBinLabel(k4Tracks+1,"4Tracks");
		fHistEvent->GetXaxis()->SetBinLabel(k2Tracks_ITSSA+1,"2Tracks_ITSSA");
		fHistEvent->GetXaxis()->SetBinLabel(k4Tracks_ITSSA+1,"4Tracks_ITSSA");
		fList->Add(fHistEvent);
	}
	if (!fHistTrigger) {
		fHistTrigger = new TH1D("fHistTrigger","",kTriAll,0,kTriAll);
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CINT5+1,"CINT5-B-NOPF-ALLNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CINT7+1,"CINT7-B-NOPF-ALLNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CINT10+1,"CINT10-B-NOPF-ALLNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CINT11+1,"CINT11-B-NOPF-CENTNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CADAND+1,"CADAND-B-NOPF-ALLNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_C0SMB+1,"C0SMB-B-NOPF-ALLNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CDG6+1,"CDG6-B-NOPF-CENTNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CDG6_SPD2+1,"CDG6-B-SPD2-CENTNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CDG7_SPD2+1,"CDG7-B-SPD2-CENTNOTRD");
		fList->Add(fHistTrigger);
	}
	if (!fHistEventProcesses) {
		fHistEventProcesses = new TH1D("fHistEventProcesses","",kBinMCAll,0,kBinMCAll);
		fHistEventProcesses->GetXaxis()->SetBinLabel(kBinND+1,"ND");
		fHistEventProcesses->GetXaxis()->SetBinLabel(kBinCD+1,"CD");
		fHistEventProcesses->GetXaxis()->SetBinLabel(kBinSD1+1,"SD1");
		fHistEventProcesses->GetXaxis()->SetBinLabel(kBinSD2+1,"SD2");
		fHistEventProcesses->GetXaxis()->SetBinLabel(kBinDD+1,"DD");
		fHistEventProcesses->GetXaxis()->SetBinLabel(kBinEL+1,"Elastic");
		fList->Add(fHistEventProcesses);
	}
	if (!fHistPrimVtxX) {
		fHistPrimVtxX = new TH1D("fHistPrimVtxX","",1000,-2,2);
		fList->Add(fHistPrimVtxX);
	}
	if (!fHistPrimVtxY) {
		fHistPrimVtxY = new TH1D("fHistPrimVtxY","",1000,-2,2);
		fList->Add(fHistPrimVtxY);
	}
	if (!fHistPrimVtxZ) {
		fHistPrimVtxZ = new TH1D("fHistPrimVtxZ","",1000,-20,20);
		fList->Add(fHistPrimVtxZ);
	}
	if (!fMultRegionsMC) {
		fMultRegionsMC = new TH1D("fMultRegions","",4,0,4);
		fList->Add(fMultRegionsMC);
	}
	if (!fSPDFiredChipsClvsHd) {
		fSPDFiredChipsClvsHd = new TH2D("fSPDFiredChipsClvsHd","",1200,0,1200,1200,0,1200);
		fList->Add(fSPDFiredChipsClvsHd);
	}
	if (!fSPDwc) {
		fSPDwc = new TH1D("fSPDwc","",1200,0,1200);
		fList->Add(fSPDwc);
	}
	if (!fRunVsMBOR_V0) {
		fRunVsMBOR_V0 = new TH1D("fRunVsMBOR_V0","",diff,minRn,maxRn);
		fList->Add(fRunVsMBOR_V0);
	}
	if (!fRunVsMBAND_V0) {
		fRunVsMBAND_V0 = new TH1D("fRunVsMBAND_V0","",diff,minRn,maxRn);
		fList->Add(fRunVsMBAND_V0);
	}
	if (!fRunVsMBOR_AD) {
		fRunVsMBOR_AD = new TH1D("fRunVsMBOR_AD","",diff,minRn,maxRn);
		fList->Add(fRunVsMBOR_AD);
	}
	if (!fRunVsMBAND_AD) {
		fRunVsMBAND_AD = new TH1D("fRunVsMBAND_AD","",diff,minRn,maxRn);
		fList->Add(fRunVsMBAND_AD);
	}
	if (!fRunVsMBOR_Global) {
		fRunVsMBOR_Global = new TH1D("fRunVsMBOR_Global","",diff,minRn,maxRn);
		fList->Add(fRunVsMBOR_Global);
	}
	if (!fRunVsMBAND_Global) {
		fRunVsMBAND_Global = new TH1D("fRunVsMBAND_Global","",diff,minRn,maxRn);
		fList->Add(fRunVsMBAND_Global);
	}
	Char_t tname[50];
	const Int_t trigger_no = 8;
	const char *trigger_name[trigger_no] = {
		"CINT5-B-NOPF-ALLNOTRD",
		"CINT7-B-NOPF-ALLNOTRD",
		"CINT10-B-NOPF-ALLNOTRD",
		"C0SMB-B-NOPF-ALLNOTRD",
		"CADAND-B-NOPF-ALLNOTRD",
		"CDG6-B-NOPF-CENTNOTRD",
		"CDG6-B-SPD2-CENTNOTRD",
		"CDG7-B-SPD2-CENTNOTRD"
	};
	for (Int_t i = 0; i < trigger_no; i++) {
		sprintf(tname,"%s",trigger_name[i]);
		fRunOnline[i] = new TH1D(tname,"",diff,minRn,maxRn);
		fList->Add(fRunOnline[i]);
	}
	for (Int_t i = 0; i < 10; i++) {
		if (!fRunVsDG[i]) {
			fRunVsDG[i] = new TH1D(Form("fRunVsDG_%d",i),"",diff,minRn,maxRn);
			fList->Add(fRunVsDG[i]);
		}
	}
	if (!fRunVs2t) {
		fRunVs2t = new TH1D("fRunVs2t","",diff,minRn,maxRn);
		fList->Add(fRunVs2t);
	}
	if (!fRunVs2t_ITSSA) {
		fRunVs2t_ITSSA = new TH1D("fRunVs2t_ITSSA","",diff,minRn,maxRn);
		fList->Add(fRunVs2t_ITSSA);
	}
	if (!fRunVs4t) {
		fRunVs4t = new TH1D("fRunVs4t","",diff,minRn,maxRn);
		fList->Add(fRunVs4t);
	}
	if (!fRunVs4t_ITSSA) {
		fRunVs4t_ITSSA = new TH1D("fRunVs4t_ITSSA","",diff,minRn,maxRn);
		fList->Add(fRunVs4t_ITSSA);
	}
	if (!fMultNG_ST) {
		fMultNG_ST = new TH1D("fMultNG_ST","",100,0,100);
		fList->Add(fMultNG_ST);
	}
	if (!fMultNG_MS) {
		fMultNG_MS = new TH1D("fMultNG_MS","",100,0,100);
		fList->Add(fMultNG_MS);
	}
	if (!fMultDG_ST) {
		fMultDG_ST = new TH1D("fMultDG_ST","",100,0,100);
		fList->Add(fMultDG_ST);
	}
	if (!fMultDG_MS) {
		fMultDG_MS = new TH1D("fMultDG_MS","",100,0,100);
		fList->Add(fMultDG_MS);
	}
	if (!fMassNG_ST_2t) {
		fMassNG_ST_2t = new TH1D("fMassNG_ST_2t","",200,0,2);
		fList->Add(fMassNG_ST_2t);
	}
	if (!fMassNG_MS_2t) {
		fMassNG_MS_2t = new TH1D("fMassNG_MS_2t","",200,0,2);
		fList->Add(fMassNG_MS_2t);
	}
	if (!fMassDG_ST_2t) {
		fMassDG_ST_2t = new TH1D("fMassDG_ST_2t","",200,0,2);
		fList->Add(fMassDG_ST_2t);
	}
	if (!fMassDG_MS_2t) {
		fMassDG_MS_2t = new TH1D("fMassDG_MS_2t","",200,0,2);
		fList->Add(fMassDG_MS_2t);
	}
	if (!fMassNG_ST_4t) {
		fMassNG_ST_4t = new TH1D("fMassNG_ST_4t","",500,0,5);
		fList->Add(fMassNG_ST_4t);
	}
	if (!fMassNG_MS_4t) {
		fMassNG_MS_4t = new TH1D("fMassNG_MS_4t","",500,0,5);
		fList->Add(fMassNG_MS_4t);
	}
	if (!fMassDG_ST_4t) {
		fMassDG_ST_4t = new TH1D("fMassDG_ST_4t","",500,0,5);
		fList->Add(fMassDG_ST_4t);
	}
	if (!fMassDG_MS_4t) {
		fMassDG_MS_4t = new TH1D("fMassDG_MS_4t","",500,0,5);
		fList->Add(fMassDG_MS_4t);
	}
	if (!fTrackCutsInfo) {
		fTrackCutsInfo = new TH1D("fTrackCutsInfo","",20,-10,10);
		fList->Add(fTrackCutsInfo);
	}
	if (!fTrackCutsInfo_ITSSA) {
		fTrackCutsInfo_ITSSA = new TH1D("fTrackCutsInfo_ITSSA","",20,-10,10);
		fList->Add(fTrackCutsInfo_ITSSA);
	}
	if (!fhClusterVsTracklets_bf) {
		fhClusterVsTracklets_bf = new TH2D("fhClusterVsTracklets_bf","",200,0,200,1000,0,1000);
		fList->Add(fhClusterVsTracklets_bf);
	}
	if (!fhClusterVsTracklets_af) {
		fhClusterVsTracklets_af = new TH2D("fhClusterVsTracklets_af","",200,0,200,1000,0,1000);
		fList->Add(fhClusterVsTracklets_af);
	}
	if (!fMult_ITSSA) {
		fMult_ITSSA = new TH1D("fMult_ITSSA","",200,0,200);
		fList->Add(fMult_ITSSA);
	}
	if (!fMult_DG_ITSSA) {
		fMult_DG_ITSSA = new TH1D("fMult_DG_ITSSA","",200,0,200);
		fList->Add(fMult_DG_ITSSA);
	}
	if (!fMult_NG_ITSSA) {
		fMult_NG_ITSSA = new TH1D("fMult_NG_ITSSA","",200,0,200);
		fList->Add(fMult_NG_ITSSA);
	}
	if (!fMult_ITSSA_MS) {
		fMult_ITSSA_MS = new TH1D("fMult_ITSSA_MS","",200,0,200);
		fList->Add(fMult_ITSSA_MS);
	}
	if (!fMult_DG_ITSSA_MS) {
		fMult_DG_ITSSA_MS = new TH1D("fMult_DG_ITSSA_MS","",200,0,200);
		fList->Add(fMult_DG_ITSSA_MS);
	}
	if (!fMult_NG_ITSSA_MS) {
		fMult_NG_ITSSA_MS = new TH1D("fMult_NG_ITSSA_MS","",200,0,200);
		fList->Add(fMult_NG_ITSSA_MS);
	}
	if (!fTPCSignal) {
		fTPCSignal = new TH2D("fTPCSignal","",500,0,5,400,0,200);
		fList->Add(fTPCSignal);
	}
	if (!fTOFSignal) {
		fTOFSignal = new TH2D("fTOFSignal","",500,0,5,400,0,1.1);
		fList->Add(fTOFSignal);
	}
	if (!fITSSignal) {
		fITSSignal = new TH2D("fITSSignal","",500,0,5,400,0,200);
		fList->Add(fITSSignal);
	}
	if (!fTRDSignal) {
		fTRDSignal = new TH2D("fTRDSignal","",500,0,5,400,0,50);
		fList->Add(fTRDSignal);
	}
	for (Int_t i = 0; i < 6; i++) {
		if (!fRunFiducial[i]) {
			fRunFiducial[i] = new TH1D(Form("fRunFiducial_%d",i),"",diff,minRn,maxRn);
			fList->Add(fRunFiducial[i]);
		}
	}
	if (fIsMC) {
		for (Int_t i = 0; i < 7; i++) {
			if (!fMC_Eta[i]) {
				fMC_Eta[i] = new TH1D(Form("fMC_Eta_%d",i),"",300,-15,15);
				fList->Add(fMC_Eta[i]);
			}
			if (!fMC_DiffMass[i]) {
				fMC_DiffMass[i] = new TH1D(Form("fMC_DiffMass_%d",i),"",1000,0,1000);
				fList->Add(fMC_DiffMass[i]);
			}
			if (!fMC_DiffMass_PDG[i]) {
				fMC_DiffMass_PDG[i] = new TH1D(Form("fMC_DiffMass_PDG_%d",i),"",1000,0,1000);
				fList->Add(fMC_DiffMass_PDG[i]);
			}
		}
		if (!fMult_Gen) {
			fMult_Gen = new TH1D("fMult_Gen","",500,0,500);
			fList->Add(fMult_Gen);
		}
		if (!fMult_Gen_Process) {
			fMult_Gen_Process = new TH2D("fMult_Gen_Process","",kBinMCAll,0,kBinMCAll,500,0,500);
			fList->Add(fMult_Gen_Process);
		}
		if (!fMult_Rec_DG_Process) {
			fMult_Rec_DG_Process = new TH2D("fMult_Rec_DG_Process","",kBinMCAll,0,kBinMCAll,500,0,500);
			fList->Add(fMult_Rec_DG_Process);
		}
		if (!fMult_Rec_NG_Process) {
			fMult_Rec_NG_Process = new TH2D("fMult_Rec_NG_Process","",kBinMCAll,0,kBinMCAll,500,0,500);
			fList->Add(fMult_Rec_NG_Process);
		}
	}
	for (Int_t i = 0; i < 2; i++) {
		if (!fV0Time[i]) {
			fV0Time[i] = new TH2D(Form("fV0Time_%d",i),"",4000,-200,200,4000,-200,200);
			fList->Add(fV0Time[i]);
		}
		if (!fADTime[i]) {
			fADTime[i] = new TH2D(Form("fADTime_%d",i),"",4000,-200,200,4000,-200,200);
			fList->Add(fADTime[i]);
		}
	}


	// Track Cuts
	if (!fIsRun2) {
		fTrackCuts = new AliESDtrackCuts();
		fTrackCuts->GetStandardITSTPCTrackCuts2010(1,0);//0 for d and e
	}
	else {
		fTrackCuts = new AliESDtrackCuts();
		{
			fTrackCuts -> SetMaxDCAToVertexXYPtDep("(0.0182+0.0350/pt^1.01)");
			fTrackCuts -> SetMinNCrossedRowsTPC(120);
			fTrackCuts -> SetMaxDCAToVertexZ(2);
			fTrackCuts -> SetEtaRange(-0.9,0.9);
			fTrackCuts -> SetMaxChi2PerClusterTPC(4);
			fTrackCuts -> SetRequireTPCRefit(kTRUE);
			fTrackCuts -> SetRequireITSRefit(kTRUE);
			fTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
			fTrackCuts -> SetAcceptKinkDaughters(kFALSE);
			fTrackCuts -> SetMaxChi2PerClusterITS(36);
			fTrackCuts -> SetMaxChi2TPCConstrainedGlobal(36);
			fTrackCuts -> SetPtRange(0.15);
			fTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
			fTrackCuts -> SetMaxFractionSharedTPCClusters(0.4);
		}
	}

	//PID Combined
	fPIDCombined = new AliPIDCombined;
	fPIDCombined->SetDefaultTPCPriors();//Need more update..
	fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);//This is default
	fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);//Do we need??

	PostOutputs();
	//-------------------------------------------------------------------------
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::UserExec(Option_t *)
{
	//Initialize
	fEventInfo.Clear();
	fCombInfo.Clear();
	fTrackInfo.Delete();
	if (fIsPWAMC) fMCTrack.Delete();

	//Check corrupted file, load ESDevent and handler
	if (!CheckInput()) {
		AliFatal("Corrupted file");
		return;
	}

	//Begin analysis
	//Fill Event information
	fEventInfo.Fill(fESDEvent);
	fHistEvent->Fill(kInput);
	fRunFiducial[0]->Fill(fRunNumber);

	//MC studies
	Bool_t passMC = kFALSE;
	if (fIsMC) {
		if (fIsPWAMC) passMC = DoMCPWA();
		if (passMC) fHistEvent->Fill(kMCCheck);
	}
	
	//Online trigger only in Run2
	if (fIsRun2 && !fIsMC && !CheckOnlineTrigger(fESDEvent)) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kOnlineTrigger);
	fRunFiducial[1]->Fill(fRunNumber);

	//Offline MB_OR cut (Ground condition)
	//Time measurement before cut
	DoTimeMeasurements(fESDEvent, 0);
	Bool_t lstPass = kFALSE;
	Bool_t normalCut = kFALSE;
	AliInputEventHandler *inputHandler = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	if (!fIsRun2) {//For Run1
		if (inputHandler->IsEventSelected() & AliVEvent::kMB) {//passed
			fCombInfo.fComb_IsPassMBOR = kTRUE;
			fHistEvent->Fill(kOfflineCut);
			fRunFiducial[2]->Fill(fRunNumber);
			lstPass = kTRUE;
			normalCut = kTRUE;
		}
		else fCombInfo.fComb_IsPassMBOR = kFALSE;
	}
	else {//For Run2
		if (inputHandler->IsEventSelected() & AliVEvent::kUserDefined) {
			fCombInfo.fComb_IsPassMBOR = kTRUE;
			fHistEvent->Fill(kOfflineCut);
			fRunFiducial[2]->Fill(fRunNumber);
			lstPass = kTRUE;
			normalCut = kTRUE;
		}
		else fCombInfo.fComb_IsPassMBOR = kFALSE;
	}
	if (lstPass) DoTimeMeasurements(fESDEvent, 1);

	//Event Selection 1: Vertex Cut
	//10cm: default but 15cm is allowed to do systematic study
	fCombInfo.fComb_IsPassVertex = DoVertexCut(fESDEvent);
	lstPass &= fCombInfo.fComb_IsPassVertex;
	normalCut &= fEventInfo.fSysVertex[0];
	if (normalCut) {
		fHistEvent->Fill(kVtxCut);
		fRunFiducial[3]->Fill(fRunNumber);
	}

	//Event Selection 2: Pileup Cut
	Double_t pileup_Ncont = 2.;
	Double_t pileup_dist = 0.8;
	fCombInfo.fComb_IsPileUp = fESDEvent->IsPileupFromSPD(pileup_Ncont,pileup_dist,3.,2.,5.);
	lstPass &= !fCombInfo.fComb_IsPileUp;
	normalCut &= !fCombInfo.fComb_IsPileUp;
	//Systematics for pile up cut
	fEventInfo.fSysPileUp[0] = fCombInfo.fComb_IsPileUp;
	fEventInfo.fSysPileUp[1] = fESDEvent->IsPileupFromSPD(pileup_Ncont+1,pileup_dist,3.,2.,5.);
	fEventInfo.fSysPileUp[2] = fESDEvent->IsPileupFromSPD(pileup_Ncont+2,pileup_dist,3.,2.,5.);
	if (normalCut) {
		fHistEvent->Fill(kPileUpCut);
		fRunFiducial[4]->Fill(fRunNumber);
	}

	//Event Selection 3: Cluster cut
	//Normal first
	//Histogram before cutting
	Int_t nClustersLayer0 = fESDEvent->GetNumberOfITSClusters(0);
	Int_t nClustersLayer1 = fESDEvent->GetNumberOfITSClusters(1);
	Int_t nTracklets      = fMult->GetNumberOfTracklets();
	fCombInfo.fComb_SPDCluster = nClustersLayer0 + nClustersLayer1;
	fCombInfo.fComb_SPDTracklets = nTracklets;
	fhClusterVsTracklets_bf->Fill(nTracklets,nClustersLayer0+nClustersLayer1);

	fCombInfo.fComb_IsPassClusterCut = !fAnalysisUtils.IsSPDClusterVsTrackletBG(fESDEvent);
	normalCut &= fCombInfo.fComb_IsPassClusterCut;
	//Systematic
	fEventInfo.fSysCluster[0] = !fCombInfo.fComb_IsPassClusterCut;
	fAnalysisUtils.SetASPDCvsTCut(55); fAnalysisUtils.SetBSPDCvsTCut(4);
	fEventInfo.fSysCluster[1] = !fAnalysisUtils.IsSPDClusterVsTrackletBG(fESDEvent);//Tighter
	fAnalysisUtils.SetASPDCvsTCut(75); fAnalysisUtils.SetBSPDCvsTCut(4);
	fEventInfo.fSysCluster[2] = !fAnalysisUtils.IsSPDClusterVsTrackletBG(fESDEvent);//Looser
	fAnalysisUtils.SetASPDCvsTCut(65); fAnalysisUtils.SetBSPDCvsTCut(3);
	fEventInfo.fSysCluster[3] = !fAnalysisUtils.IsSPDClusterVsTrackletBG(fESDEvent);//Tighter
	fAnalysisUtils.SetASPDCvsTCut(65); fAnalysisUtils.SetBSPDCvsTCut(5);
	fEventInfo.fSysCluster[4] = !fAnalysisUtils.IsSPDClusterVsTrackletBG(fESDEvent);//Looser
	//For lst Pass
	fAnalysisUtils.SetASPDCvsTCut(75); fAnalysisUtils.SetBSPDCvsTCut(5);
	lstPass &= !fAnalysisUtils.IsSPDClusterVsTrackletBG(fESDEvent);

	if (normalCut) {
		fHistEvent->Fill(kClusterCut);
		fRunFiducial[5]->Fill(fRunNumber);
		fhClusterVsTracklets_af->Fill(nTracklets,nClustersLayer0+nClustersLayer1);
	}

	//Combinatorics and Fill eventinfo.Gap
	DoCombinatorics(fESDEvent);

	//Reject not lst events to reduce memory leak
	if (fCombmode) fTree_Comb->Fill();
	if (!lstPass) {
		PostOutputs();
		return;
	}

	//CEP ANALYSIS START FROM HERE!
	//Detector QA procedures
	Double_t fTOFlength = 0.;
	Double_t timeTOF = 0.;
	Double_t fTOFTime = 0.;
	Double_t fTOFPIDmom = 0.;
	Double_t fBeta = 0.;
	for (Int_t iTrk = 0; iTrk < fESDEvent->GetNumberOfTracks(); iTrk++) {
		AliESDtrack *trk = (AliESDtrack*)fESDEvent->GetTrack(iTrk);
		if(!trk) continue;
		if(!fTrackCuts->AcceptTrack(trk)) continue;
		if (normalCut) {
			fTPCSignal->Fill(trk->P(),trk->GetTPCsignal());
			fITSSignal->Fill(trk->P(),trk->GetITSsignal());
			fTRDSignal->Fill(trk->P(),trk->GetTRDsignal());
		}
		//// *** TOF *** //// for beta extraction
		fTOFlength = trk->GetIntegratedLength();
		if(fTOFlength<=0) continue;
		fTOFPIDmom = trk->P();
		if((fTOFPIDmom*fTOFPIDmom)==0) continue;
		fTOFlength = fTOFlength*0.01;
		timeTOF = trk->GetTOFsignal();
		fTOFTime = timeTOF*1E-3*TMath::C()*1.E-9;
		fBeta = fTOFlength/fTOFTime;
		if (normalCut) fTOFSignal->Fill(trk->P(),fBeta);
	}

	//Trigger analysis
	Bool_t IsMBOR_V0 = (fCombInfo.fComb_DetHit[1] || fCombInfo.fComb_DetHit[2] || fCombInfo.fComb_DetHit[0]) ? kTRUE : kFALSE;
	Bool_t IsMBAND_V0 = (fCombInfo.fComb_DetHit[1] && fCombInfo.fComb_DetHit[2]) ? kTRUE : kFALSE;
	Bool_t IsMBOR_AD = (fCombInfo.fComb_DetHit[7] || fCombInfo.fComb_DetHit[8] || fCombInfo.fComb_DetHit[0]) ? kTRUE : kFALSE;
	Bool_t IsMBAND_AD = (fCombInfo.fComb_DetHit[7] && fCombInfo.fComb_DetHit[8]) ? kTRUE : kFALSE;
	Bool_t IsMBOR_Global = (fCombInfo.fComb_DetHit[0] || fCombInfo.fComb_DetHit[1] || fCombInfo.fComb_DetHit[2] || fCombInfo.fComb_DetHit[7] || fCombInfo.fComb_DetHit[8]) ? kTRUE : kFALSE;
	Bool_t IsMBAND_Global = ((fCombInfo.fComb_DetHit[1] || fCombInfo.fComb_DetHit[7]) && (fCombInfo.fComb_DetHit[2] || fCombInfo.fComb_DetHit[8])) ? kTRUE : kFALSE;

	if (normalCut) {
		if (IsMBOR_V0) {fHistEvent->Fill(kMBOR_V0); fRunVsMBOR_V0->Fill(fRunNumber);}
		if (IsMBAND_V0) {fHistEvent->Fill(kMBAND_V0); fRunVsMBAND_V0->Fill(fRunNumber);}
		if (IsMBOR_AD) {fHistEvent->Fill(kMBOR_AD); fRunVsMBOR_AD->Fill(fRunNumber);}
		if (IsMBAND_AD) {fHistEvent->Fill(kMBAND_AD); fRunVsMBAND_AD->Fill(fRunNumber);}
		if (IsMBOR_Global) {fHistEvent->Fill(kMBOR_Global); fRunVsMBOR_Global->Fill(fRunNumber);}
		if (IsMBAND_Global) {fHistEvent->Fill(kMBAND_Global); fRunVsMBAND_Global->Fill(fRunNumber);}
	}

	//Gap analysis
	//V0+SPD, AD+SPD are basic gap def.
	Bool_t fDG_Det[10]={kFALSE,};
	fDG_Det[0] = (fEventInfo.fV0Gap & fEventInfo.fSPDFired);//V0SPD
	fDG_Det[1] = (fEventInfo.fADGap & fEventInfo.fSPDFired);//ADSPD
	fDG_Det[2] = (fEventInfo.fADGap & fDG_Det[0]);//V0ADSPD
	fDG_Det[3] = (fEventInfo.fFMDGap & fDG_Det[2]);//V0ADFMDSPD
	fDG_Det[4] = (fEventInfo.fZDCGap & fDG_Det[3]);//V0ADFMDZDCSPD
	fDG_Det[5] = (fEventInfo.fFMDGap & fDG_Det[0]);//V0FMDSPD
	fDG_Det[6] = (fEventInfo.fZDCGap & fDG_Det[5]);//V0FMDZDCSPD
	//No-gap and single-gap
	fDG_Det[8] = (!fCombInfo.fComb_DetHit[1] & fCombInfo.fComb_DetHit[2]);//Single-gap at A-side
	fDG_Det[9] = (!fCombInfo.fComb_DetHit[2] & fCombInfo.fComb_DetHit[1]);//Single-gap at C-side
	fDG_Det[7] = (!fDG_Det[0] & !fDG_Det[8] & !fDG_Det[9]);//No-gap

	//Fill histo
	if (normalCut) {
		for (Int_t i = 0; i < 10; i++) {
			if (fDG_Det[i]) {
				fHistEvent->Fill(kDGV0SPD+i);
				fRunVsDG[i]->Fill(fRunNumber);
			}
		}
	}

	//Standard track cuts in NG and DG
	Double_t pionmass = 0.139570;
	TLorentzVector lv_track1;
	TLorentzVector lv_sum_2t;
	TLorentzVector lv_sum_4t;
	Int_t nmbTrk_ST = 0;
	Int_t trkCharge = 0;
	//N-tracks
	for (Int_t iTrack = 0; iTrack < fESDEvent->GetNumberOfTracks(); iTrack++) {
		AliESDtrack *track1 = fESDEvent->GetTrack(iTrack);
		if (!track1) continue;
		if (!fTrackCuts->AcceptTrack(track1)) continue;
		nmbTrk_ST++;
	}
	//2-tracks
	if (nmbTrk_ST == 2) {
		for (Int_t iTrack = 0; iTrack < fESDEvent->GetNumberOfTracks(); iTrack++) {
			AliESDtrack *track1 = fESDEvent->GetTrack(iTrack);
			if (!track1) continue;
			if (!fTrackCuts->AcceptTrack(track1)) continue;
			lv_track1.SetXYZM(track1->Px(),track1->Py(),track1->Pz(),pionmass);
			lv_sum_2t += lv_track1;
			trkCharge += (Int_t)track1->GetSign();
		}
	}
	//4-tracks
	else if (nmbTrk_ST == 4) {
		trkCharge = 0;
		for (Int_t iTrack = 0; iTrack < fESDEvent->GetNumberOfTracks(); iTrack++) {
			AliESDtrack *track1 = fESDEvent->GetTrack(iTrack);
			if (!track1) continue;
			if (!fTrackCuts->AcceptTrack(track1)) continue;
			lv_track1.SetXYZM(track1->Px(),track1->Py(),track1->Pz(),pionmass);
			lv_sum_4t += lv_track1;
			trkCharge += (Int_t)track1->GetSign();
		}
	}
	if (normalCut && fDG_Det[7]) {//NG
		fMultNG_ST->Fill(nmbTrk_ST);
		if (nmbTrk_ST == 2 && trkCharge == 0) fMassNG_ST_2t->Fill(lv_sum_2t.M());
		if (nmbTrk_ST == 4 && trkCharge == 0) fMassNG_ST_4t->Fill(lv_sum_4t.M());
	}
	if (normalCut && fDG_Det[0]) {//DG_V0SPD
		fMultDG_ST->Fill(nmbTrk_ST);
		if (nmbTrk_ST == 2 && trkCharge == 0) fMassDG_ST_2t->Fill(lv_sum_2t.M());
		if (nmbTrk_ST == 4 && trkCharge == 0) fMassDG_ST_4t->Fill(lv_sum_4t.M());
	}

	fEventInfo.fNtrk_ST = nmbTrk_ST;

	//Martin's track cuts in NG and DG
	//two/four track event investigation in DG and NG
	Bool_t is10bc = kFALSE;
	if (fEventInfo.fPeriod <= 2) is10bc=kTRUE;
	AliMultiplicitySelectionCPPWA selec;
	//(clusterCut,useITSSA,isRun2,nCluster)
	if (!fIsRun2) selec.InitDefaultTrackCuts((Int_t)is10bc,0,kFALSE,0);// 1 = b,c and 0 = d,e
	else selec.InitDefaultTrackCuts(0,0,kTRUE,0);//Run2
	TArrayI indices;
	Int_t Nsel = selec.GetNumberOfITSTPCtracks(fESDEvent,indices);
	
	//2-tracks & 4-tracks
	lv_sum_2t.SetPxPyPzE(0,0,0,0);
	lv_sum_4t.SetPxPyPzE(0,0,0,0);
	trkCharge = 0;
	if (Nsel == 2) {
		for (Int_t iTrack = 0; iTrack < 2; iTrack++) {
			AliESDtrack *track1 = fESDEvent->GetTrack(indices.At(iTrack));
			if(!track1) continue;
			lv_track1.SetXYZM(track1->Px(),track1->Py(),track1->Pz(),pionmass);
			lv_sum_2t += lv_track1;
			trkCharge += (Int_t)track1->GetSign();
		}
	}
	else if (Nsel == 4) {
		trkCharge = 0;
		for (Int_t iTrack = 0; iTrack < 4; iTrack++) {
			AliESDtrack *track1 = fESDEvent->GetTrack(indices.At(iTrack));
			if(!track1) continue;
			lv_track1.SetXYZM(track1->Px(),track1->Py(),track1->Pz(),pionmass);
			lv_sum_4t += lv_track1;
			trkCharge += (Int_t)track1->GetSign();
		}
	}

	if (normalCut && fDG_Det[7]) {//NG
		fMultNG_MS->Fill(Nsel);
		if (Nsel == 2 && trkCharge == 0) fMassNG_MS_2t->Fill(lv_sum_2t.M());
		if (Nsel == 4 && trkCharge == 0) fMassNG_MS_4t->Fill(lv_sum_4t.M());
	}
	if (normalCut && fDG_Det[0]) {//DG_V0SPD
		fMultDG_MS->Fill(Nsel);
		if (Nsel == 2) {
			if (trkCharge == 0) fMassDG_MS_2t->Fill(lv_sum_2t.M());
			fHistEvent->Fill(k2Tracks);
			fRunVs2t->Fill(fRunNumber);
		}
		else if (Nsel == 4) {
			if (trkCharge == 0) fMassDG_MS_4t->Fill(lv_sum_4t.M());
			fHistEvent->Fill(k4Tracks);
			fRunVs4t->Fill(fRunNumber);
		}
		//Underflow
		fTrackCutsInfo->Fill(Nsel);
	}

	fEventInfo.fNtrk_MS = Nsel;

	//Move to the DG Event to reduce memory leak
	if (!fDG_Det[0] || !fDG_Det[1]) {
		PostOutputs();
		return;
	}

	//Store normal version
	Bool_t storeT = kFALSE;

	if (!fIsMC) {
		for (Int_t i = 0; i < 11; i++) {
			Nsel = 0;
			selec.Clear();
			indices = 0x0;
			if (!fIsRun2) selec.InitDefaultTrackCuts((Int_t)is10bc,0,kFALSE,i);
			else selec.InitDefaultTrackCuts(0,0,kTRUE,0);
			Nsel = selec.GetNumberOfITSTPCtracks(fESDEvent,indices);
			fEventInfo.fCheckTwoTrack[i] = (Nsel == 2) ? kTRUE : kFALSE;
			fEventInfo.fCheckFourTrack[i] = (Nsel == 4) ? kTRUE : kFALSE;

			if (Nsel == 2 || Nsel == 4) {
				storeT = kTRUE;
				for (Int_t iSel = 0; iSel < Nsel; iSel++) {
					new (fTrackInfo[Nsel*i+iSel]) TrackInfo(dynamic_cast<AliESDtrack*>(fESDEvent->GetTrack(indices.At(iSel))), fPIDResponse, fPIDCombined);
				}
			}
		}
	}
	else if (fIsMC && fIsPWAMC) {
		//By default it has to run it over with 'fnSys'
		selec.Clear();
		indices = 0x0;
		if (!fIsRun2) selec.InitDefaultTrackCuts((Int_t)is10bc,0,kFALSE,fnSys);
		else selec.InitDefaultTrackCuts(0,0,kTRUE,0);
		Nsel = selec.GetNumberOfITSTPCtracks(fESDEvent,indices);
		fEventInfo.fCheckTwoTrack[0] = (Nsel == 2) ? kTRUE : kFALSE;
		fEventInfo.fCheckFourTrack[0] = (Nsel == 4) ? kTRUE : kFALSE;
		for (Int_t iSel = 0; iSel < Nsel; iSel++) {
			new (fTrackInfo[iSel]) TrackInfo(dynamic_cast<AliESDtrack*>(fESDEvent->GetTrack(indices.At(iSel))), fPIDResponse, fPIDCombined);
		}

	}

	if (!fIsMC) {
		if (storeT) fTree->Fill();
	}

	PostOutputs();
	return;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::PostOutputs()
{
	if (fIsMC) {
		fTree->Fill();
		if (fCombmode) fTree_Comb->Fill();
	}
	PostData(1,fTree);
	PostData(2,fTree_Comb);
	PostData(3,fList);
	return;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::CheckInput()
{
	//Event
	AliVEvent *event = InputEvent();
	if (NULL == event) return 0;

	fESDEvent = dynamic_cast<AliESDEvent*>(InputEvent());
	if (NULL == fESDEvent) return 0; 
	
	//Input Handler
	const AliAnalysisManager* man(AliAnalysisManager::GetAnalysisManager());
	if (NULL == man) return 0;

	AliESDInputHandler* inputHandler(dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler()));  
	if (NULL == inputHandler) return 0;

	fPIDResponse = inputHandler->GetPIDResponse();
	if (NULL == fPIDResponse) return 0;

	const AliESDHeader *esdHeader = fESDEvent->GetHeader();
	if (NULL == esdHeader) return 0;

	if (kFALSE == fIsMC && esdHeader->GetEventType() != AliRawEventHeaderBase::kPhysicsEvent)
	return 0;

	fMult = fESDEvent->GetMultiplicity();
	if (NULL == fMult) return 0;

	//Magnetic field
	if(TMath::Abs(fESDEvent->GetMagneticField()) < 1) return 0;

	if (fIsMC) {
		fMCEvent = MCEvent();
		if (NULL == fMCEvent) return 0;
	}

	//AliTriggerAnalysis
	fTriggerAnalysis.SetDoFMD(kTRUE);
	fTriggerAnalysis.SetFMDThreshold(0.3,0.5);

	//RunNumber
	fRunNumber = fESDEvent->GetRunNumber();

	return kTRUE;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::CheckOnlineTrigger(
		const AliESDEvent *esd
		)
{
	const Int_t trigger_no = 9;
	const char *trigger_name[trigger_no] = {
		"CINT5-B-NOPF-ALLNOTRD",
		"CINT7-B-NOPF-ALLNOTRD",
		"CINT10-B-NOPF-ALLNOTRD",
		"CINT11-B-NOPF-CENTNOTRD",
		"CADAND-B-NOPF-ALLNOTRD",
		"C0SMB-B-NOPF-ALLNOTRD",
		"CDG6-B-NOPF-CENTNOTRD",
		"CDG6-B-SPD2-CENTNOTRD",
		"CDG7-B-SPD2-CENTNOTRD"
	};
	Char_t tname[50];
	for (Int_t i = 0; i < trigger_no; i++) {
		sprintf(tname,"%s",trigger_name[i]);
		if (esd->IsTriggerClassFired(tname)) {
			fHistTrigger->Fill(kTri_CINT5+i);
			fRunOnline[i]->Fill(fRunNumber);
		}
	}

	if (esd->IsTriggerClassFired("CINT10-B-NOPF-ALLNOTRD") ||
			esd->IsTriggerClassFired("C0SMB-B-NOPF-ALLNOTRD") ||
			esd->IsTriggerClassFired("CINT11-B-NOPF-CENTNOTRD"))
		return kTRUE;
	else return kFALSE;

	return kTRUE;

}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::DoTimeMeasurements(
		const AliESDEvent *esd,
		const Int_t seq
		)
{
	AliESDVZERO *vzero = (AliESDVZERO*)esd->GetVZEROData();
	fV0Time[seq]->Fill(vzero->GetV0ATime()-vzero->GetV0CTime(),vzero->GetV0ATime()+vzero->GetV0CTime());
	if (fIsRun2) {
		AliESDAD *ad = (AliESDAD*)esd->GetADData();
		fADTime[seq]->Fill(ad->GetADATime()-ad->GetADCTime(),ad->GetADATime()+ad->GetADCTime());
	}

	//SPD Check
	if (seq == 1) {
		fSPDFiredChipsClvsHd->Fill(fTriggerAnalysis.SPDFiredChips(esd,0),fTriggerAnalysis.SPDFiredChips(esd,1));
		for (Int_t i = 0; i < 1200; i++) {
			if (fMult->TestFastOrFiredChips(i)) {
				fSPDwc->Fill((Double_t)i);
			}
		}
	}
	return;
		
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::DoVertexCut(
		const AliESDEvent *esd
		)
{
	//Get Vertex first
	const AliESDVertex *vtxTracks = esd->GetPrimaryVertexTracks();
	const AliESDVertex *vtxSPD = esd->GetPrimaryVertexSPD();
	const AliESDVertex *vtxTPC = esd->GetPrimaryVertexTPC();

	//Fill eventinfo first
	fEventInfo.fVertexTracks[0] = vtxTracks->GetX();
	fEventInfo.fVertexTracks[1] = vtxTracks->GetY();
	fEventInfo.fVertexTracks[2] = vtxTracks->GetZ();
	fEventInfo.fVertexSPD[0] = vtxSPD->GetX();
	fEventInfo.fVertexSPD[1] = vtxSPD->GetY();
	fEventInfo.fVertexSPD[2] = vtxSPD->GetZ();
	fEventInfo.fVertexTPC[0] = vtxTPC->GetX();
	fEventInfo.fVertexTPC[1] = vtxTPC->GetY();
	fEventInfo.fVertexTPC[2] = vtxTPC->GetZ();

	//Normal primary vertex reconstruction
	if (vtxTracks->GetNContributors() < 1) {
		vtxTracks = esd->GetPrimaryVertexSPD();
		if (vtxTracks->GetNContributors() < 1) {//No vertex pass event
			return kFALSE;
		}
	}
	fEventInfo.fVertexUsed[0] = vtxTracks->GetX();
	fEventInfo.fVertexUsed[1] = vtxTracks->GetY();
	fEventInfo.fVertexUsed[2] = vtxTracks->GetZ();

	Double_t cutValue = 15.;
	Bool_t isPass = (TMath::Abs(vtxTracks->GetZ()) < cutValue) ? kTRUE : kFALSE;

	if (isPass) {
		fHistPrimVtxX->Fill(vtxTracks->GetX());
		fHistPrimVtxY->Fill(vtxTracks->GetY());
		fHistPrimVtxZ->Fill(vtxTracks->GetZ());
	}

	//Fill Systematic Info
	fEventInfo.fSysVertex[0] = (TMath::Abs(vtxTracks->GetZ()) < 10.) ? kTRUE : kFALSE;
	fEventInfo.fSysVertex[1] = (TMath::Abs(vtxTracks->GetZ()) < 15.) ? kTRUE : kFALSE;
	fEventInfo.fSysVertex[2] = (TMath::Abs(vtxTracks->GetZ()) < 5.) ? kTRUE : kFALSE;

	return isPass;

}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::DoCombinatorics(
		const AliESDEvent *esd
		)
{
	fCombInfo.fComb_DetHit[0] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kSPDGFO);
	fCombInfo.fComb_DetHit[1] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kV0A);
	fCombInfo.fComb_DetHit[2] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kV0C);
	fCombInfo.fComb_DetHit[3] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kFMDA);
	fCombInfo.fComb_DetHit[4] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kFMDC);
	fCombInfo.fComb_DetHit[5] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kZDCA);
	fCombInfo.fComb_DetHit[6] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kZDCC);
	if (fIsRun2) {
		fCombInfo.fComb_DetHit[7] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kADA);
		fCombInfo.fComb_DetHit[8] = fTriggerAnalysis.IsOfflineTriggerFired(esd,AliTriggerAnalysis::kADC);
	}

	//Event info
	fEventInfo.fSPDFired = fCombInfo.fComb_DetHit[0];
	fEventInfo.fV0Gap = !fCombInfo.fComb_DetHit[1] & !fCombInfo.fComb_DetHit[2];
	fEventInfo.fFMDGap = !fCombInfo.fComb_DetHit[3] & !fCombInfo.fComb_DetHit[4];
	fEventInfo.fZDCGap = !fCombInfo.fComb_DetHit[5] & !fCombInfo.fComb_DetHit[6];
	if (fIsRun2) fEventInfo.fADGap = !fCombInfo.fComb_DetHit[7] & !fCombInfo.fComb_DetHit[8];

	return;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::DoMCPWA()
{
	AliStack* stack = fMCEvent->Stack();
	if (!stack) return kFALSE;
	TParticle *part1 = NULL;
	Int_t n_proton = 0;
	Int_t n_pion = 0;
	Int_t nPrimaries = stack->GetNprimary();

	// Track Gap && Multiplicity in Region Bins -------------------------------
	Int_t gapA = 0;
	Int_t gapAv0 = 0;
	Int_t gapAv0fmd = 0;
	Int_t gapAad = 0;
	Int_t gapC = 0;
	Int_t gapCv0 = 0;
	Int_t gapCv0fmd = 0;
	Int_t gapCad = 0;
	Int_t central = 0;
	for (Int_t iTracks = 0; iTracks <  nPrimaries; ++iTracks) {
		TParticle *part = (TParticle*)stack->Particle(iTracks);
		if (!part) continue;
		if (stack->IsPhysicalPrimary(iTracks) && (part->GetPDG()->Charge() != 0.) && part->GetStatusCode() == 1) {
			if (part->Eta() > -0.9 && part->Eta() < 0.9) central++;
			if (part->Eta() > 0.9 && part->Eta() < 6.3) gapA++;
			if (part->Eta() > 2.8 && part->Eta() < 5.1) gapAv0++;
			if (part->Eta() > 1.7 && part->Eta() < 5.1) gapAv0fmd++;
			if (part->Eta() > 4.77 && part->Eta() < 6.30) gapAad++;
			if ((part->Eta() < -0.9 && part->Eta() > -3.7) || (part->Eta() < -4.92 && part->Eta() > -6.96)) gapC++;
			if (part->Eta() < -1.9 && part->Eta() > -3.7) gapCv0++;
			if (part->Eta() < -1.9 && part->Eta() > -3.7) gapCv0fmd++;
			if (part->Eta() < -4.92 && part->Eta() > -6.96) gapCad++;
		}
	}

	if (fMultRegionsMC) {
		// multiplicity distribution separated in A-side, central barrel and C-side
		fMultRegionsMC->Fill(gapA, 0);
		fMultRegionsMC->Fill(central, 1);
		fMultRegionsMC->Fill(gapC, 2);
		fMultRegionsMC->Fill(gapA+gapC, 3);
	}

	//Store generated protons and pions for PWA
	if (fIsPWAMC) {
		for (Int_t j = 0; j < nPrimaries; j++) {
			part1 = (TParticle*)stack->Particle(j);
			if (!part1) continue;
			if (stack->IsPhysicalPrimary(j) && (part1->GetPDG()->Charge() != 0.) && part1->GetStatusCode() == 1 && TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 2212) {
				n_proton++;
				continue;
			}
			if (stack->IsPhysicalPrimary(j) && (part1->GetPDG()->Charge() != 0.) && part1->GetStatusCode() == 1 && (TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 211 || TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 321)) {
				n_pion++;
				continue;
			}
		}
		if (n_proton != 2 || n_pion !=2) return kFALSE;

		//Store
		Int_t nCount = 0;
		for (Int_t j = 0; j < nPrimaries; j++) {
			part1 = (TParticle*)stack->Particle(j);
			if (!part1) continue;
			if (stack->IsPhysicalPrimary(j) && (part1->GetPDG()->Charge() != 0.) && part1->GetStatusCode() == 1 && TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 2212) {
				new (fMCTrack[nCount]) TLorentzVector(part1->Px(),part1->Py(),part1->Pz(),part1->Energy());
				nCount++;
			}
			if (stack->IsPhysicalPrimary(j) && (part1->GetPDG()->Charge() != 0.) && part1->GetStatusCode() == 1 && (TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 211 || TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 321)) {
				new (fMCTrack[nCount]) TLorentzVector(part1->Px(),part1->Py(),part1->Pz(),part1->Energy());
				nCount++;
			}
		}

	}

	return kTRUE;

}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::Terminate(Option_t *)
{
}
//______________________________________________________________________________
