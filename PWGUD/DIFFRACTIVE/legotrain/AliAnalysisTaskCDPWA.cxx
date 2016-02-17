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

#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TTree.h>
#include <THnSparse.h>
#include <TArrayI.h>
#include <TRandom3.h>
#include <TString.h>

#include "AliTriggerAnalysis.h"
#include "AliESDInputHandler.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliSPDUtils.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSAlignMille2Module.h"
#include "AliPhysicsSelection.h"
#include "AliESDtrackCuts.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliGeomManager.h"
#include "AliESDAD.h"
#include "AliESDVZERO.h"
#include "TGeoManager.h"

#include "AliMultiplicitySelectionCPPWA.h"
#include "AliAnalysisTaskCDPWA.h"

//______________________________________________________________________________
AliAnalysisTaskCDPWA::AliAnalysisTaskCDPWA(const char* name):
	AliAnalysisTaskSE(name)
	, fMaxVtxDst(0.5) // value to be checked with the vertex study histograms
	, fPIDResponse(0x0)
	, fPIDCombined(0x0)
	, fESDEvent(0x0)
	, fPhysicsSelection(0x0)
	, fTrackCuts(0x0)
	, fTrackCuts_ITSSA(0x0)
	, fVtxDst(-1.)
	, fVtxZ(-20)
	, fMCprocessType(0)
	, fMCprocess(-1)
	, fRun(-999)
	, fPIDmode(0)
	, fTree(0x0)
	, fCheckTwoPion(0)
	, fCheckFourPion(0)
	, fCheckTwoPion_ITSSA(0)
	, fDGV0SPD(0)
	, fDGADSPD(0)
	, fSPDFired(0)
	, fV0Gap(0)
	, fADGap(0)
	, fFMDGap(0)
	, fZDCGap(0)
	, fIsMC(0)
	, fRunNumber(0)
	, fList(0x0)
	, fHistEvent(0x0)
	, fHistTrigger(0x0)
	, fHistEventProcesses(0x0)
	, fHistPrimVtxX(0x0)
	, fHistPrimVtxY(0x0)
	, fHistPrimVtxZ(0x0)
	, fHitV0A(0x0)
	, fHitV0C(0x0)
	, fMultRegionsMC(0x0)
	, fSPDFiredChipsCluster(0x0)
	, fSPDFiredChipsHardware(0x0)
	, fSPDwc(0x0)
	, fRunVsMBOR_V0(0x0)
	, fRunVsMBAND_V0(0x0)
	, fRunVsMBOR_AD(0x0)
	, fRunVsMBAND_AD(0x0)
	, fRunVsMBOR_Global(0x0)
	, fRunVsMBAND_Global(0x0)
	, fRunVsDG_V0(0x0)
	, fRunVsDG_AD(0x0)
	, fRunVsDG_V0SPD(0x0)
	, fRunVsDG_ADSPD(0x0)
	, fRunVsDG_V0ADSPD(0x0)
	, fRunVsDG_V0ADFMDSPD(0x0)
	, fRunVsDG_V0ADFMDZDCSPD(0x0)
	, fRunVs2t(0x0)
	, fRunVs2t_ITSSA(0x0)
	, fMultNG(0x0)
	, fMultNG_MS(0x0)
	, fMultDG(0x0)
	, fMultDG_MS(0x0)
	, fMassNG(0x0)
	, fMassNG_MS(0x0)
	, fMassDG(0x0)
	, fTrackCutsInfo(0x0)
	, fhClusterVsTracklets_bf(0x0)
	, fhClusterVsTracklets_af(0x0)
	, fMult_ITSSA(0x0)
	, fMult_DG_ITSSA(0x0)
	, fADATime_bf(0x0)
	, fADATime_af(0x0)
	, fADCTime_bf(0x0)
	, fADCTime_af(0x0)
	, fV0ATime_bf(0x0)
	, fV0ATime_af(0x0)
	, fV0CTime_bf(0x0)
	, fV0CTime_af(0x0)
	, fADATime_V0SPD(0x0)
	, fADCTime_V0SPD(0x0)
	, fADATime_V0ADSPD(0x0)
	, fADCTime_V0ADSPD(0x0)

{
	//
	// standard constructor (the one which should be used)
	//
	// slot in TaskSE must start from 1
	
	for (Int_t i = 0; i < 3; i++) {
		fVertex[i] = -999.;
	}
	for (Int_t i = 0; i < 2; i++) {//Two Pion
		fTwoPionMask[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fTwoPionTrack[i][j] = 0;
			fTwoPionTPCSigma[i][j] = 0;
			fTwoPionTOFSigma[i][j] = 0;
			fTwoPionITSSigma[i][j] = 0;
			fTwoPionDetMask[i][j] = 0;
			fTwoPionTrack_ITSSA[i][j] = 0;
			if ( j < 5) fTwoPionBayesProb[i][j] = 0;
		}
	}
	for (Int_t i = 0; i < 4; i++) {//Four Pion
		fFourPionMask[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fFourPionTrack[i][j] = 0;
			fFourPionTPCSigma[i][j] = 0;
			fFourPionTOFSigma[i][j] = 0;
			fFourPionITSSigma[i][j] = 0;
			fFourPionDetMask[i][j] = 0;
			fFourPionTrack_ITSSA[i][j] = 0;
			if ( j < 5) fFourPionBayesProb[i][j] = 0;
		}
	}


	// MC variable
	for (Int_t i=0; i< 5; i++) {
		for (Int_t j = 0; j < 2; j++) {
			fMCGenProtonTrack[j][i] = 0.;
			fMCGenPionTrack[j][i] = 0.;
		}
	}
	DefineOutput(1, TTree::Class());
	DefineOutput(2, TList::Class());
}
//______________________________________________________________________________
AliAnalysisTaskCDPWA::AliAnalysisTaskCDPWA():
	AliAnalysisTaskSE()
	, fMaxVtxDst(0.5)
	, fPIDResponse(0x0)
	, fPIDCombined(0x0)
	, fESDEvent(0x0)
	, fPhysicsSelection(0x0)
	, fTrackCuts(0x0)
	, fTrackCuts_ITSSA(0x0)
	, fVtxDst(-1.)
	, fVtxZ(-20)
	, fMCprocessType(0)
	, fMCprocess(-1)
	, fRun(-999)
	, fPIDmode(0)
	, fTree(0x0)
	, fCheckTwoPion(0)
	, fCheckFourPion(0)
	, fCheckTwoPion_ITSSA(0)
	, fDGV0SPD(0)
	, fDGADSPD(0)
	, fSPDFired(0)
	, fV0Gap(0)
	, fADGap(0)
	, fFMDGap(0)
	, fZDCGap(0)
	, fIsMC(0)
	, fRunNumber(0)
	, fList(0x0)
	, fHistEvent(0x0)
	, fHistTrigger(0x0)
	, fHistEventProcesses(0x0)
	, fHistPrimVtxX(0x0)
	, fHistPrimVtxY(0x0)
	, fHistPrimVtxZ(0x0)
	, fHitV0A(0x0)
	, fHitV0C(0x0)
	, fMultRegionsMC(0x0)
	, fSPDFiredChipsCluster(0x0)
	, fSPDFiredChipsHardware(0x0)
	, fSPDwc(0x0)
	, fRunVsMBOR_V0(0x0)
	, fRunVsMBAND_V0(0x0)
	, fRunVsMBOR_AD(0x0)
	, fRunVsMBAND_AD(0x0)
	, fRunVsMBOR_Global(0x0)
	, fRunVsMBAND_Global(0x0)
	, fRunVsDG_V0(0x0)
	, fRunVsDG_AD(0x0)
	, fRunVsDG_V0SPD(0x0)
	, fRunVsDG_ADSPD(0x0)
	, fRunVsDG_V0ADSPD(0x0)
	, fRunVsDG_V0ADFMDSPD(0x0)
	, fRunVsDG_V0ADFMDZDCSPD(0x0)
	, fRunVs2t(0x0)
	, fRunVs2t_ITSSA(0x0)
	, fMultNG(0x0)
	, fMultNG_MS(0x0)
	, fMultDG(0x0)
	, fMultDG_MS(0x0)
	, fMassNG(0x0)
	, fMassNG_MS(0x0)
	, fMassDG(0x0)
	, fTrackCutsInfo(0x0)
	, fhClusterVsTracklets_bf(0x0)
	, fhClusterVsTracklets_af(0x0)
	, fMult_ITSSA(0x0)
	, fMult_DG_ITSSA(0x0)
	, fADATime_bf(0x0)
	, fADATime_af(0x0)
	, fADCTime_bf(0x0)
	, fADCTime_af(0x0)
	, fV0ATime_bf(0x0)
	, fV0ATime_af(0x0)
	, fV0CTime_bf(0x0)
	, fV0CTime_af(0x0)
	, fADATime_V0SPD(0x0)
	, fADCTime_V0SPD(0x0)
	, fADATime_V0ADSPD(0x0)
	, fADCTime_V0ADSPD(0x0)
{
	for (Int_t i = 0; i < 3; i++) {
		fVertex[i] = -999.;
	}
	for (Int_t i = 0; i < 2; i++) {//Two Pion
		fTwoPionMask[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fTwoPionTrack[i][j] = 0;
			fTwoPionTPCSigma[i][j] = 0;
			fTwoPionTOFSigma[i][j] = 0;
			fTwoPionITSSigma[i][j] = 0;
			fTwoPionDetMask[i][j] = 0;
			fTwoPionTrack_ITSSA[i][j] = 0;
			if ( j < 5) fTwoPionBayesProb[i][j] = 0;
		}
	}
	for (Int_t i = 0; i < 4; i++) {//Four Pion
		fFourPionMask[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fFourPionTrack[i][j] = 0;
			fFourPionTPCSigma[i][j] = 0;
			fFourPionTOFSigma[i][j] = 0;
			fFourPionITSSigma[i][j] = 0;
			fFourPionDetMask[i][j] = 0;
			fFourPionTrack_ITSSA[i][j] = 0;
			if ( j < 5) fFourPionBayesProb[i][j] = 0;
		}
	}

	for (Int_t i=0; i< 5; i++) {
		for (Int_t j = 0; j < 2; j++) {
			fMCGenProtonTrack[j][i] = 0.;
			fMCGenPionTrack[j][i] = 0.;
		}
	}
}
//______________________________________________________________________________
AliAnalysisTaskCDPWA::~AliAnalysisTaskCDPWA()
{
	//Destructor(pointer should be deleted)
	if (fTree) {
		delete fTree;
		fTree = 0x0;
	}
	if (fList) {
		delete fList;
		fList = 0x0;
	}
	if (fPhysicsSelection) {
		delete fPhysicsSelection;
		fPhysicsSelection = 0x0;
	}
	if (fTrackCuts) {
		delete fTrackCuts;
	}
	if (fTrackCuts_ITSSA) {
		delete fTrackCuts_ITSSA;
	}
	if (fHistEvent) {
		delete fHistEvent;
		fHistEvent = 0x0;
	}
	if (fPIDResponse) delete fPIDResponse;
	if (fPIDCombined) delete fPIDCombined;
	if (fHistEventProcesses) delete fHistEventProcesses;
	if (fHistTrigger) delete fHistTrigger;
	if (fHistPrimVtxX) delete fHistPrimVtxX;
	if (fHistPrimVtxY) delete fHistPrimVtxY;
	if (fHistPrimVtxZ) delete fHistPrimVtxZ;
	if (fHitV0A) delete fHitV0A;
	if (fHitV0C) delete fHitV0C;
	if (fMultRegionsMC) delete fMultRegionsMC;
	if (fSPDFiredChipsCluster) delete fSPDFiredChipsCluster;
	if (fSPDFiredChipsHardware) delete fSPDFiredChipsHardware;
	if (fSPDwc) delete fSPDwc;
	if (fRunVsMBOR_V0) delete fRunVsMBOR_V0;
	if (fRunVsMBAND_V0) delete fRunVsMBAND_V0;
	if (fRunVsMBOR_AD) delete fRunVsMBOR_AD;
	if (fRunVsMBAND_AD) delete fRunVsMBAND_AD;
	if (fRunVsMBOR_Global) delete fRunVsMBOR_Global;
	if (fRunVsMBAND_Global) delete fRunVsMBAND_Global;
	if (fRunVsDG_V0) delete fRunVsDG_V0;
	if (fRunVsDG_AD) delete fRunVsDG_AD;
	if (fRunVsDG_V0SPD) delete fRunVsDG_V0SPD;
	if (fRunVsDG_ADSPD) delete fRunVsDG_ADSPD;
	if (fRunVsDG_V0ADSPD) delete fRunVsDG_V0ADSPD;
	if (fRunVsDG_V0ADFMDSPD) delete fRunVsDG_V0ADFMDSPD;
	if (fRunVsDG_V0ADFMDZDCSPD) delete fRunVsDG_V0ADFMDZDCSPD;
	if (fRunVs2t) delete fRunVs2t;
	if (fRunVs2t_ITSSA) delete fRunVs2t_ITSSA;
	if (fMultNG) delete fMultNG;
	if (fMultNG_MS) delete fMultNG_MS;
	if (fMultDG) delete fMultDG;
	if (fMultDG_MS) delete fMultDG_MS;
	if (fMassNG) delete fMassNG;
	if (fMassNG_MS) delete fMassNG_MS;
	if (fMassDG) delete fMassDG;
	if (fTrackCutsInfo) delete fTrackCutsInfo;
	if (fhClusterVsTracklets_bf) delete fhClusterVsTracklets_bf;
	if (fhClusterVsTracklets_af) delete fhClusterVsTracklets_af;
	if (fMult_ITSSA) delete fMult_ITSSA;
	if (fMult_DG_ITSSA) delete fMult_DG_ITSSA;
	if (fADATime_bf) delete fADATime_bf;
	if (fADATime_af) delete fADATime_af;
	if (fADCTime_bf) delete fADCTime_bf;
	if (fADCTime_af) delete fADCTime_af;
	if (fV0ATime_bf) delete fV0ATime_bf;
	if (fV0ATime_af) delete fV0ATime_af;
	if (fV0CTime_bf) delete fV0CTime_bf;
	if (fV0CTime_af) delete fV0CTime_af;
	if (fADATime_V0SPD) delete fADATime_V0SPD;
	if (fADCTime_V0SPD) delete fADCTime_V0SPD;
	if (fADATime_V0ADSPD) delete fADATime_V0ADSPD;
	if (fADCTime_V0ADSPD) delete fADCTime_V0ADSPD;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::UserCreateOutputObjects()
{
	// TTree + TList should be defined-----------------------------------------
	OpenFile(1);
	fTree = new TTree("tree1","PWAtree");
	//For Data + Rec
	fTree->Branch("CheckTwoPion",&fCheckTwoPion);
	fTree->Branch("CheckFourPion",&fCheckFourPion);
	fTree->Branch("CheckTwoPion_ITSSA",&fCheckTwoPion_ITSSA);
	fTree->Branch("DGV0SPD",&fDGV0SPD);
	fTree->Branch("DGADSPD",&fDGADSPD);
	fTree->Branch("SPDFired",&fSPDFired);
	fTree->Branch("V0Gap",&fV0Gap);
	fTree->Branch("ADGap",&fADGap);
	fTree->Branch("FMDGap",&fFMDGap);
	fTree->Branch("ZDCGap",&fZDCGap);
	fTree->Branch("RunNumber",&fRunNumber);
	for (Int_t i = 0; i < 2; i++) {
		fTree->Branch(Form("TwoPionMask_%d",i),&fTwoPionMask[i]);
		for (Int_t j = 0; j < 9; j++) {
			fTree->Branch(Form("TwoPionTrack_%d_%d",i,j),&fTwoPionTrack[i][j]);
			fTree->Branch(Form("TwoPionTPCSigma_%d_%d",i,j),&fTwoPionTPCSigma[i][j]);
			fTree->Branch(Form("TwoPionTOFSigma_%d_%d",i,j),&fTwoPionTOFSigma[i][j]);
			fTree->Branch(Form("TwoPionITSSigma_%d_%d",i,j),&fTwoPionITSSigma[i][j]);
			fTree->Branch(Form("TwoPionDetMask_%d_%d",i,j),&fTwoPionDetMask[i][j]);
			fTree->Branch(Form("TwoPionTrack_ITSSA_%d_%d",i,j),&fTwoPionTrack_ITSSA[i][j]);
			if ( j < 5) fTree->Branch(Form("TwoPionBayesProb_%d_%d",i,j),&fTwoPionBayesProb[i][j]);
		}
	}
	for (Int_t i = 0; i < 4; i++) {
		fTree->Branch(Form("FourPionMask_%d",i),&fFourPionMask[i]);
		for (Int_t j = 0; j < 9; j++) {
			fTree->Branch(Form("FourPionTrack_%d_%d",i,j),&fFourPionTrack[i][j]);
			fTree->Branch(Form("FourPionTPCSigma_%d_%d",i,j),&fFourPionTPCSigma[i][j]);
			fTree->Branch(Form("FourPionTOFSigma_%d_%d",i,j),&fFourPionTOFSigma[i][j]);
			fTree->Branch(Form("FourPionITSSigma_%d_%d",i,j),&fFourPionITSSigma[i][j]);
			fTree->Branch(Form("FourPionDetMask_%d_%d",i,j),&fFourPionDetMask[i][j]);
			fTree->Branch(Form("FourPionTrack_ITSSA_%d_%d",i,j),&fFourPionTrack_ITSSA[i][j]);
			if ( j < 5) fTree->Branch(Form("FourPionBayesProb_%d_%d",i,j),&fFourPionBayesProb[i][j]);
		}
	}
	//For MC
	for (Int_t i = 0; i < 5; i++) {
		for (Int_t j = 0; j < 2; j++) {
			fTree->Branch(Form("MCGenProtonTrack_%d_%d",j,i),&fMCGenProtonTrack[j][i]);
			fTree->Branch(Form("MCGenPionTrack_%d_%d",j,i),&fMCGenPionTrack[j][i]);
		}
	}
	//For vertex
	for (Int_t i = 0; i < 3; i++) {
		fTree->Branch(Form("Vertex_%d",i),&fVertex[i]);
	}
	//fTree->Branch("test",&fPIDResponse);
	//-------------------------------------------------------------------------

	// TList-------------------------------------------------------------------
	Int_t minRn = 224800;
	Int_t maxRn = 240000;
	Int_t diff = maxRn - minRn;

	OpenFile(2);
	fList = new TList();
	fList->SetOwner();
	if (!fHistEvent) {
		fHistEvent = new TH1D("fHistEvent","Number of Events for each selection",kAll,0,kAll);
		fHistEvent->GetXaxis()->SetBinLabel(kInput+1,"InputEvent");
		fHistEvent->GetXaxis()->SetBinLabel(kMCCheck+1,"MCCheck");
		fHistEvent->GetXaxis()->SetBinLabel(kOnlineTrigger+1,"OnlineTrigger");
		fHistEvent->GetXaxis()->SetBinLabel(kMBOR+1,"OfflineMBOR");
		fHistEvent->GetXaxis()->SetBinLabel(kVtxCut+1,"VertexCut");
		fHistEvent->GetXaxis()->SetBinLabel(kPileUpCut+1,"PileUpCut");
		fHistEvent->GetXaxis()->SetBinLabel(kClusterCut+1,"ClusterCut");
		fHistEvent->GetXaxis()->SetBinLabel(kMBOR_V0+1,"MBOR_V0");
		fHistEvent->GetXaxis()->SetBinLabel(kMBAND_V0+1,"MBAND_V0");
		fHistEvent->GetXaxis()->SetBinLabel(kMBOR_AD+1,"MBOR_AD");
		fHistEvent->GetXaxis()->SetBinLabel(kMBAND_AD+1,"MBAND_AD");
		fHistEvent->GetXaxis()->SetBinLabel(kMBOR_Global+1,"MBOR_Global");
		fHistEvent->GetXaxis()->SetBinLabel(kMBAND_Global+1,"MBAND_Global");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0+1,"DG_V0");
		fHistEvent->GetXaxis()->SetBinLabel(kDGAD+1,"DG_AD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0SPD+1,"DG_V0SPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGADSPD+1,"DG_ADSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADSPD+1,"DG_V0ADSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADFMDSPD+1,"DG_V0ADFMDSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADFMDZDCSPD+1,"DG_V0ADFMDZDCSPD");
		fHistEvent->GetXaxis()->SetBinLabel(k2Tracks+1,"2Tracks");
		fHistEvent->GetXaxis()->SetBinLabel(k4Tracks+1,"4Tracks");
		fList->Add(fHistEvent);
	}
	if (!fHistTrigger) {
		fHistTrigger = new TH1D("fHistTrigger","",kTriAll,0,kTriAll);
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CINT5+1,"CINT5-B-NOPF-ALLNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CINT7+1,"CINT7-B-NOPF-ALLNOTRD");
		fHistTrigger->GetXaxis()->SetBinLabel(kTri_CINT10+1,"CINT10-B-NOPF-ALLNOTRD");
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
		fHistEventProcesses->GetXaxis()->SetBinLabel(kBinSD2+1,"SD1");
		fHistEventProcesses->GetXaxis()->SetBinLabel(kBinDD+1,"DD");
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
	if (!fHitV0A) {
		fHitV0A = new TH1D("fHitV0A","Hit on V0A(0:No hit, 1:Hit)",2,0,2);
		fList->Add(fHitV0A);
	}
	if (!fHitV0C) {
		fHitV0C = new TH1D("fHitV0C","Hit on V0C(0:No hit, 1:Hit)",2,0,2);
		fList->Add(fHitV0C);
	}
	if (!fMultRegionsMC) {
		fMultRegionsMC = new TH1D("fMultRegions","",4,0,4);
		fList->Add(fMultRegionsMC);
	}
	if (!fSPDFiredChipsCluster) {
		fSPDFiredChipsCluster = new TH1D("fSPDFiredChipsCluster","",1200,0,1200);
		fList->Add(fSPDFiredChipsCluster);
	}
	if (!fSPDFiredChipsHardware) {
		fSPDFiredChipsHardware = new TH1D("fSPDFiredChipsHardware","",1200,0,1200);
		fList->Add(fSPDFiredChipsHardware);
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
	if (!fRunVsDG_V0) {
		fRunVsDG_V0 = new TH1D("fRunVsDG_V0","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_V0);
	}
	if (!fRunVsDG_AD) {
		fRunVsDG_AD = new TH1D("fRunVsDG_AD","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_AD);
	}
	if (!fRunVsDG_V0SPD) {
		fRunVsDG_V0SPD = new TH1D("fRunVsDG_V0SPD","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_V0SPD);
	}
	if (!fRunVsDG_ADSPD) {
		fRunVsDG_ADSPD = new TH1D("fRunVsDG_ADSPD","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_ADSPD);
	}
	if (!fRunVsDG_V0ADSPD) {
		fRunVsDG_V0ADSPD = new TH1D("fRunVsDG_V0ADSPD","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_V0ADSPD);
	}
	if (!fRunVsDG_V0ADFMDSPD) {
		fRunVsDG_V0ADFMDSPD = new TH1D("fRunVsDG_V0ADFMDSPD","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_V0ADFMDSPD);
	}
	if (!fRunVsDG_V0ADFMDZDCSPD) {
		fRunVsDG_V0ADFMDZDCSPD = new TH1D("fRunVsDG_V0ADFMDZDCSPD","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_V0ADFMDZDCSPD);
	}
	if (!fRunVs2t) {
		fRunVs2t = new TH1D("fRunVs2t","",diff,minRn,maxRn);
		fList->Add(fRunVs2t);
	}
	if (!fRunVs2t_ITSSA) {
		fRunVs2t_ITSSA = new TH1D("fRunVs2t_ITSSA","",diff,minRn,maxRn);
		fList->Add(fRunVs2t_ITSSA);
	}
	if (!fMultNG) {
		fMultNG = new TH1D("fMultNG","",100,0,100);
		fList->Add(fMultNG);
	}
	if (!fMultNG_MS) {
		fMultNG_MS = new TH1D("fMultNG_MS","",100,0,100);
		fList->Add(fMultNG_MS);
	}
	if (!fMultDG) {
		fMultDG = new TH1D("fMultDG","",100,0,100);
		fList->Add(fMultDG);
	}
	if (!fMultDG_MS) {
		fMultDG_MS = new TH1D("fMultDG_MS","",100,0,100);
		fList->Add(fMultDG_MS);
	}
	if (!fMassNG) {
		fMassNG = new TH1D("fMassNG","",50,0,2);
		fList->Add(fMassNG);
	}
	if (!fMassNG_MS) {
		fMassNG_MS = new TH1D("fMassNG_MS","",50,0,2);
		fList->Add(fMassNG_MS);
	}
	if (!fMassDG) {
		fMassDG = new TH1D("fMassDG","",50,0,2);
		fList->Add(fMassDG);
	}
	if (!fTrackCutsInfo) {
		fTrackCutsInfo = new TH1D("fTrackCutsInfo","",20,-10,10);
		fList->Add(fTrackCutsInfo);
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
	if (!fADATime_bf) {
		fADATime_bf = new TH1D("fADATime_bf","",5000,-100,400);
		fList->Add(fADATime_bf);
	}
	if (!fADATime_af) {
		fADATime_af = new TH1D("fADATime_af","",5000,-100,400);
		fList->Add(fADATime_af);
	}
	if (!fADCTime_bf) {
		fADCTime_bf = new TH1D("fADCTime_bf","",5000,-100,400);
		fList->Add(fADCTime_bf);
	}
	if (!fADCTime_af) {
		fADCTime_af = new TH1D("fADCTime_af","",5000,-100,400);
		fList->Add(fADCTime_af);
	}
	if (!fV0ATime_bf) {
		fV0ATime_bf = new TH1D("fV0ATime_bf","",8000,-400,400);
		fList->Add(fV0ATime_bf);
	}
	if (!fV0ATime_af) {
		fV0ATime_af = new TH1D("fV0ATime_af","",8000,-400,400);
		fList->Add(fV0ATime_af);
	}
	if (!fV0CTime_bf) {
		fV0CTime_bf = new TH1D("fV0CTime_bf","",8000,-400,400);
		fList->Add(fV0CTime_bf);
	}
	if (!fV0CTime_af) {
		fV0CTime_af = new TH1D("fV0CTime_af","",8000,-400,400);
		fList->Add(fV0CTime_af);
	}
	if (!fADATime_V0SPD) {
		fADATime_V0SPD = new TH1D("fADATime_V0SPD","",5000,-100,400);
		fList->Add(fADATime_V0SPD);
	}
	if (!fADCTime_V0SPD) {
		fADCTime_V0SPD = new TH1D("fADCTime_V0SPD","",5000,-100,400);
		fList->Add(fADCTime_V0SPD);
	}
	if (!fADATime_V0ADSPD) {
		fADATime_V0ADSPD = new TH1D("fADATime_V0ADSPD","",5000,-100,400);
		fList->Add(fADATime_V0ADSPD);
	}
	if (!fADCTime_V0ADSPD) {
		fADCTime_V0ADSPD = new TH1D("fADCTime_V0ADSPD","",5000,-100,400);
		fList->Add(fADCTime_V0ADSPD);
	}

	// Track Cuts
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
	fTrackCuts_ITSSA = new AliESDtrackCuts();
	{
		fTrackCuts_ITSSA->GetStandardITSSATrackCuts2010(1,0);//0 for noPID
	}

	//PID Combined
	fPIDCombined = new AliPIDCombined;
	fPIDCombined->SetDefaultTPCPriors();
	fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
	fPIDCombined->SetDetectorMask(AliPIDResponse::kDetITS | 
			AliPIDResponse::kDetTPC | 
			AliPIDResponse::kDetTRD | 
			AliPIDResponse::kDetTOF);

	PostOutputs();
	//-------------------------------------------------------------------------
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::UserExec(Option_t *)
{
	// Initialize all tree varialbles------------------------------------------
	for (Int_t i = 0; i < 2; i++) {//Two Pion
		fTwoPionMask[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fTwoPionTrack[i][j] = 0;
			fTwoPionTPCSigma[i][j] = 0;
			fTwoPionTOFSigma[i][j] = 0;
			fTwoPionITSSigma[i][j] = 0;
			fTwoPionDetMask[i][j] = 0;
			fTwoPionTrack_ITSSA[i][j] = 0;
			if ( j < 5) fTwoPionBayesProb[i][j] = 0;
		}
	}
	for (Int_t i = 0; i < 4; i++) {//Four Pion
		fFourPionMask[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fFourPionTrack[i][j] = 0;
			fFourPionTPCSigma[i][j] = 0;
			fFourPionTOFSigma[i][j] = 0;
			fFourPionITSSigma[i][j] = 0;
			fFourPionDetMask[i][j] = 0;
			fFourPionTrack_ITSSA[i][j] = 0;
			if ( j < 5) fFourPionBayesProb[i][j] = 0;
		}
	}
	for (Int_t i = 0; i < 5; i++) {
		for (Int_t j = 0; j < 2; j++) {
			fMCGenProtonTrack[j][i] = 0.;
			fMCGenPionTrack[j][i] = 0.;
		}
	}
	for (Int_t i = 0; i < 3; i++) {
		fVertex[i] = -999.;
	}
	fCheckTwoPion = kFALSE;
	fCheckFourPion = kFALSE;
	fCheckTwoPion_ITSSA = kFALSE;

	fDGV0SPD = kFALSE;
	fDGADSPD = kFALSE;
	fSPDFired = kFALSE;
	fV0Gap = kFALSE;
	fADGap = kFALSE;
	fFMDGap = kFALSE;
	fZDCGap = kFALSE;
	fRunNumber = -999;
	fIsMC = kFALSE;
	//-------------------------------------------------------------------------

	// TODO:: Recommendtation by Christoph Mayer-------------------------------
	/*
	if (fESDEvent->IsIncompleteDAQ()) {
		PostOutputs();
		return;
	}
	*/
	//-------------------------------------------------------------------------

	//Check Input event is available-------------------------------------------
	//Load Input handler and MC
	if (!CheckInput()) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kInput);
	//-------------------------------------------------------------------------

	//MC TRUTH-----------------------------------------------------------------
	Int_t nMCprimaries = 0;
	DetermineMCprocessType();
	//Fill Histo for MC processes
	for (Int_t i = 0; i < kBinMCAll; i++) {
		if (fMCprocessType == i) 
			fHistEventProcesses->Fill(i);
	}
	if (fMCEvent) {
		nMCprimaries = DoMCTruth();
		if (nMCprimaries == -1) { // Something is wrong with MC
			PostOutputs();
			return;
		}
	}
	fHistEvent->Fill(kMCCheck);
	//-------------------------------------------------------------------------

	//At least we have to use CINT10 + C0SMB-----------------------------------
	if (!CheckOnlineTrigger(fESDEvent)) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kOnlineTrigger);
	//-------------------------------------------------------------------------

	//Offline MBOR cut(available only for pass2)-------------------------------
	//Timing check before cutting out
	AliESDAD *ad = (AliESDAD*)fESDEvent->GetADData();
	if(!ad) {
		PostOutputs();
		return;
	}
	DoTimeV0AD(fESDEvent, fV0ATime_bf, fV0CTime_bf, fADATime_bf, fADCTime_bf);

	AliInputEventHandler *inputHandler = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	Bool_t isMinimumBias = kFALSE;
	if (inputHandler->IsEventSelected() & AliVEvent::kUserDefined) {
		isMinimumBias = kTRUE;
	}
	else {
		PostOutputs();
		return;
	}
	//Time measurement after cutting out
	DoTimeV0AD(fESDEvent, fV0ATime_af, fV0CTime_af, fADATime_af, fADCTime_af);

	// SPD FastOR chip check to see noisy chips
	DoSPDCheck(fESDEvent, fSPDFiredChipsCluster, fSPDFiredChipsHardware, fSPDwc);
	fHistEvent->Fill(kMBOR);
	//-------------------------------------------------------------------------

	// Event selection :: Vertex cut-------------------------------------------
	const Bool_t eventIsValid = CutEvent(fESDEvent, fHistPrimVtxX, fHistPrimVtxY, fHistPrimVtxZ);
	if (!eventIsValid) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kVtxCut); //After Vertex cut
	//-------------------------------------------------------------------------

	// Event selection :: Pile-up cut -----------------------------------------
	// Which is same as dN/dpT 13TeV analysis
	const Bool_t isPileup = fESDEvent->IsPileupFromSPD(5, 0.8, 3., 2., 5.);
	if (isPileup) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kPileUpCut); // After pile up
	//-------------------------------------------------------------------------

	// Cluster Cut ------------------------------------------------------------
	// Note that this is responsible for USER!!
	Bool_t IsClusterCut = kFALSE;
	const AliMultiplicity *tmp_mult = fESDEvent->GetMultiplicity();
	Int_t nClustersLayer0 = fESDEvent->GetNumberOfITSClusters(0);
	Int_t nClustersLayer1 = fESDEvent->GetNumberOfITSClusters(1);
	Int_t nTracklets      = tmp_mult->GetNumberOfTracklets();

	fhClusterVsTracklets_bf->Fill(nTracklets,nClustersLayer0+nClustersLayer1);

	Double_t cut_slope = 4.;//4 is default
	Double_t cut_b = 65.;//65 is default
	if (nClustersLayer0 + nClustersLayer1 <= (cut_b+cut_slope*nTracklets)) IsClusterCut = kTRUE;
	if (IsClusterCut == kFALSE) {
		PostOutputs();
		return;
	}
	fhClusterVsTracklets_af->Fill(nTracklets,nClustersLayer0+nClustersLayer1);
	fHistEvent->Fill(kClusterCut); //Cluster Cut
	//-------------------------------------------------------------------------

	// TRIGGER ANALYSIS -------------------------------------------------------
	// MBOR and MBAND for Run2
	fRunNumber = (Int_t)fESDEvent->GetRunNumber();
	AliTriggerAnalysis fTrigger;
	fTrigger.SetFMDThreshold(0.3,0.5);//FMD Threshold for 7TeV is used.

	Bool_t IsV0A = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0A);
	Bool_t IsV0C = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0C);
	Bool_t IsADA = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kADA);
	Bool_t IsADC = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kADC);
	Bool_t IsSPD = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kSPDGFO);
	Bool_t IsFMDA = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kFMDA);
	Bool_t IsFMDC = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kFMDC);
//	Bool_t IsZDCA = kFALSE;
//	Bool_t IsZDCC = kFALSE;
	Bool_t IsZDCA = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kZDCA);//Not TDC info..
	Bool_t IsZDCC = fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kZDCC);//Not TDC info..

	Bool_t IsMBOR_V0 = (IsV0A || IsV0C || IsSPD) ? kTRUE : kFALSE;
	Bool_t IsMBAND_V0 = (IsV0A && IsV0C) ? kTRUE : kFALSE;
	Bool_t IsMBOR_AD = (IsADA || IsADC || IsSPD) ? kTRUE : kFALSE;
	Bool_t IsMBAND_AD = (IsADA && IsADC) ? kTRUE : kFALSE;
	Bool_t IsMBOR_Global = (IsV0A || IsV0C || IsADA || IsADC || IsSPD) ? kTRUE : kFALSE;
	Bool_t IsMBAND_Global = ((IsV0A || IsADA) && (IsV0C || IsADC)) ? kTRUE : kFALSE;

	if (IsMBOR_V0) {fHistEvent->Fill(kMBOR_V0); fRunVsMBOR_V0->Fill(fESDEvent->GetRunNumber());}
	if (IsMBAND_V0) {fHistEvent->Fill(kMBAND_V0); fRunVsMBAND_V0->Fill(fESDEvent->GetRunNumber());}
	if (IsMBOR_AD) {fHistEvent->Fill(kMBOR_AD); fRunVsMBOR_AD->Fill(fESDEvent->GetRunNumber());}
	if (IsMBAND_AD) {fHistEvent->Fill(kMBAND_AD); fRunVsMBAND_AD->Fill(fESDEvent->GetRunNumber());}
	if (IsMBOR_Global) {fHistEvent->Fill(kMBOR_Global); fRunVsMBOR_Global->Fill(fESDEvent->GetRunNumber());}
	if (IsMBAND_Global) {fHistEvent->Fill(kMBAND_Global); fRunVsMBAND_Global->Fill(fESDEvent->GetRunNumber());}
	//-------------------------------------------------------------------------

	// Determine GAP ----------------------------------------------------------
	// Below two conditions will be used as criteria
	fDGV0SPD = (!IsV0A && !IsV0C && IsSPD) ? kTRUE : kFALSE;
	fDGADSPD = (!IsADA && !IsADC && IsSPD) ? kTRUE : kFALSE;
	fSPDFired = IsSPD;

	// Sub-sample for double gap conditions
	fV0Gap = (!IsV0A && !IsV0C) ? kTRUE : kFALSE;
	fADGap = (!IsADA && !IsADC) ? kTRUE : kFALSE;
	fFMDGap = (!IsFMDA && !IsFMDC) ? kTRUE : kFALSE;
	fZDCGap = (!IsZDCA && !IsZDCC) ? kTRUE : kFALSE;

	if (fV0Gap) {fHistEvent->Fill(kDGV0); fRunVsDG_V0->Fill(fESDEvent->GetRunNumber());}
	if (fADGap) {fHistEvent->Fill(kDGAD); fRunVsDG_AD->Fill(fESDEvent->GetRunNumber());}
	if (fDGV0SPD) {fHistEvent->Fill(kDGV0SPD); fRunVsDG_V0SPD->Fill(fESDEvent->GetRunNumber());}
	if (fDGADSPD) {fHistEvent->Fill(kDGADSPD); fRunVsDG_ADSPD->Fill(fESDEvent->GetRunNumber());}

	// Local variable for tmp..
	Bool_t fDGV0ADSPD = (fDGV0SPD && fADGap) ? kTRUE : kFALSE;

	// Sub-sample histogram
	if (fDGV0SPD && fADGap) {fHistEvent->Fill(kDGV0ADSPD); fRunVsDG_V0ADSPD->Fill(fESDEvent->GetRunNumber());}
	if (fDGV0SPD && fADGap && fFMDGap) {fHistEvent->Fill(kDGV0ADFMDSPD); fRunVsDG_V0ADFMDSPD->Fill(fESDEvent->GetRunNumber());}
	if (fDGV0SPD && fADGap && fFMDGap && fZDCGap) {fHistEvent->Fill(kDGV0ADFMDZDCSPD); fRunVsDG_V0ADFMDZDCSPD->Fill(fESDEvent->GetRunNumber());}

	// Define no-gap events to compare(only using V0)
	Bool_t IsNG = kFALSE;
	Bool_t IsGapA = kFALSE;
	Bool_t IsGapC = kFALSE;

	if (IsV0A && !IsV0C) IsGapC = kTRUE;
	if (IsV0C && !IsV0A) IsGapA = kTRUE;
	if (!IsGapC && !IsGapA && !fV0Gap) IsNG = kTRUE;
	//-------------------------------------------------------------------------

	// Track cuts by Martin's selection----------------------------------------
	AliMultiplicitySelectionCPPWA *selec = new AliMultiplicitySelectionCPPWA();
	selec->InitDefaultTrackCuts(0);// 1 = 10b and 10c, 0 = 10d and 10e
	TArrayI indices;
	Int_t Nsel = selec->GetNumberOfITSTPCtracks(fESDEvent,indices);
	Double_t pionmass = 0.139570;
	TLorentzVector lv_track1;
	TLorentzVector lv_track2;
	TLorentzVector lv_sum;
	//-------------------------------------------------------------------------

	// Mass distribution for NG and DG events with Martin's trackcuts and standard trackcuts
	Int_t nmbTrk_NG = 0;
	Int_t nmbTrk_DG = 0;
	if (IsNG) {//For No-gap
		//2tracks with Standard ITSTPC cut
		for (Int_t iTrack = 0; iTrack < fESDEvent->GetNumberOfTracks(); iTrack++) {
			AliESDtrack *track1 = fESDEvent->GetTrack(iTrack);
			if (!track1) continue;
			if (!fTrackCuts->AcceptTrack(track1)) continue;
			for (Int_t jTrack = iTrack+1; jTrack < fESDEvent->GetNumberOfTracks(); jTrack++) {
				AliESDtrack *track2 = fESDEvent->GetTrack(jTrack);
				if (!track2) continue;
				if (!fTrackCuts->AcceptTrack(track2)) continue;

				if (track1->GetSign()>0 && track2->GetSign()>0) continue;
				if (track1->GetSign()<0 && track2->GetSign()<0) continue;

				lv_track1.SetXYZM(track1->Px(),track1->Py(),track1->Pz(),pionmass);
				lv_track2.SetXYZM(track2->Px(),track2->Py(),track2->Pz(),pionmass);
				lv_sum = lv_track1 + lv_track2;
				fMassNG->Fill(lv_sum.M());
			}
			nmbTrk_NG++;
		}
		//2tracks with Martin's selection
		if (Nsel == 2){
			for (Int_t iTrack = 0; iTrack < 2; iTrack++) {
				AliESDtrack *track1 = fESDEvent->GetTrack(indices.At(iTrack));
				for (Int_t jTrack = iTrack+1; jTrack < 2; jTrack++) {
					AliESDtrack *track2 = fESDEvent->GetTrack(indices.At(1));
					if(!track1 || !track2) continue;
					if (track1->GetSign()>0 && track2->GetSign()>0) continue;
					if (track1->GetSign()<0 && track2->GetSign()<0) continue;

					lv_track1.SetXYZM(track1->Px(),track1->Py(),track1->Pz(),pionmass);
					lv_track2.SetXYZM(track2->Px(),track2->Py(),track2->Pz(),pionmass);
					lv_sum = lv_track1 + lv_track2;
					fMassNG_MS->Fill(lv_sum.M());
				}
			}
		}
		fMultNG->Fill(nmbTrk_NG);
		fMultNG_MS->Fill(Nsel);
	}
	else if (fV0Gap) {//For DG only for V0 with standard trackcuts
		if (Nsel == 2) {
			fHistEvent->Fill(k2Tracks);
			fRunVs2t->Fill(fESDEvent->GetRunNumber());
		}
		else if (Nsel == 4) fHistEvent->Fill(k4Tracks);
		for (Int_t iTrack = 0; iTrack < fESDEvent->GetNumberOfTracks(); iTrack++) {
			AliESDtrack *track1 = fESDEvent->GetTrack(iTrack);
			if (!track1) continue;
			if (!fTrackCuts->AcceptTrack(track1)) continue;
			for (Int_t jTrack = iTrack+1; jTrack < fESDEvent->GetNumberOfTracks(); jTrack++) {
				AliESDtrack *track2 = fESDEvent->GetTrack(jTrack);
				if (!track2) continue;
				if (!fTrackCuts->AcceptTrack(track2)) continue;

				if (track1->GetSign()>0 && track2->GetSign()>0) continue;
				if (track1->GetSign()<0 && track2->GetSign()<0) continue;

				lv_track1.SetXYZM(track1->Px(),track1->Py(),track1->Pz(),pionmass);
				lv_track2.SetXYZM(track2->Px(),track2->Py(),track2->Pz(),pionmass);
				lv_sum = lv_track1 + lv_track2;
				fMassDG->Fill(lv_sum.M());
			}
			nmbTrk_DG++;
		}
		fMultDG_MS->Fill(Nsel);
		fMultDG->Fill(nmbTrk_DG);
	}

	Bool_t fCheck2tracks = kFALSE;
	Bool_t fCheck4tracks = kFALSE;
	Bool_t fCheck2tracks_ITSSA = kFALSE;

	delete selec;

	//ITSSA
	Int_t tmp_Ntrk = fESDEvent->GetNumberOfTracks();
	fMult_ITSSA->Fill(tmp_Ntrk);
	Int_t Ntrk_ITSSA = 0;

	for (Int_t i = 0; i < tmp_Ntrk; i++) {
		AliESDtrack *track = fESDEvent->GetTrack(i);
		if(!track) continue;
		if(!fTrackCuts->AcceptTrack(track)) continue;
		Ntrk_ITSSA++;
	}
	fCheck2tracks_ITSSA = (Ntrk_ITSSA == 2) ? kTRUE : kFALSE;
	if (fCheck2tracks_ITSSA && fV0Gap) {
		fRunVs2t_ITSSA->Fill(fESDEvent->GetRunNumber());
	}


	// TAESOO PWA TREE INFO ---------------------------------------------------
	if (fV0Gap || fADGap) {
		//printf("== This events are double gap ==\n");
	}
	else {
		//printf("== This events are not double gap! ==\n");
		PostOutputs();
		return;
	}

	if(fV0Gap) fMult_DG_ITSSA->Fill(tmp_Ntrk);
	if(fV0Gap) fTrackCutsInfo->Fill(Nsel);
	/*
	if(fDGV0SPD) {
		fADATime_V0SPD->Fill(ad->GetADATime());
		fADCTime_V0SPD->Fill(ad->GetADCTime());
	}
	if(fDGV0ADSPD) {
		fADATime_V0ADSPD->Fill(ad->GetADATime());
		fADCTime_V0ADSPD->Fill(ad->GetADCTime());
	}
	*/

	fCheck2tracks = (Nsel == 2) ? kTRUE : kFALSE;
	fCheck4tracks = (Nsel == 4) ? kTRUE : kFALSE;

	if (fCheck2tracks) {
		for (Int_t i = 0; i < 2; i++) {//Two Pion
			fTwoPionMask[i] = 0;
			for (Int_t j = 0; j < 9; j++) {
				fTwoPionTrack[i][j] = 0;
				fTwoPionTPCSigma[i][j] = 0;
				fTwoPionTOFSigma[i][j] = 0;
				fTwoPionITSSigma[i][j] = 0;
				fTwoPionDetMask[i][j] = 0;
				fTwoPionTrack_ITSSA[i][j] = 0;
				if ( j < 5) fTwoPionBayesProb[i][j] = 0;
			}
		}
		for (Int_t i = 0; i < 2; i++) {
			AliESDtrack *track = fESDEvent->GetTrack(indices.At(i));
			if(!track) {
				PostOutputs();
				return;
			}
			fTwoPionTrack[i][0] = track->Px();
			fTwoPionTrack[i][1] = track->Py();
			fTwoPionTrack[i][2] = track->Pz();
			fTwoPionTrack[i][3] = track->E();
			fTwoPionTrack[i][4] = track->GetSign();
			//TPC (by L.Goerlich)
			fTwoPionTPCSigma[i][0] = track->GetTPCsignal();
			fTwoPionTPCSigma[i][1] = (Double_t)(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track) == AliPIDResponse::kDetPidOk);
			fTwoPionTPCSigma[i][2] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
			fTwoPionTPCSigma[i][3] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
			fTwoPionTPCSigma[i][4] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);
			fTwoPionTPCSigma[i][5] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
			//TOF (by L.Goerlich)
			fTwoPionTOFSigma[i][0] = track->GetTOFsignal();
			fTwoPionTOFSigma[i][1] = (Double_t)(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk);
			fTwoPionTOFSigma[i][2] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
			fTwoPionTOFSigma[i][3] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
			fTwoPionTOFSigma[i][4] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);
			//ITS (by L.Goerlich)
			fTwoPionITSSigma[i][0] = track->GetITSsignal();
			fTwoPionITSSigma[i][1] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kPion);
			fTwoPionITSSigma[i][2] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kPion);
			fTwoPionITSSigma[i][3] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kProton);
			fTwoPionITSSigma[i][4] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kKaon);
			fTwoPionITSSigma[i][5] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kElectron);
			fTwoPionITSSigma[i][6] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kKaon);
			fTwoPionITSSigma[i][7] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kProton);
			fTwoPionITSSigma[i][8] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron);
			//Bayesian,combined(UInt_t)
			fTwoPionMask[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb[i]);
			fTwoPionDetMask[i][0] = (fTwoPionMask[i] & AliPIDResponse::kDetITS);
			fTwoPionDetMask[i][1] = (fTwoPionMask[i] & AliPIDResponse::kDetTPC);
			fTwoPionDetMask[i][2] = (fTwoPionMask[i] & AliPIDResponse::kDetTRD);
			fTwoPionDetMask[i][3] = (fTwoPionMask[i] & AliPIDResponse::kDetTOF);

		}

		fCheckTwoPion = kTRUE;
	}

	if (fCheck4tracks) {
		for (Int_t i = 0; i < 4; i++) {//Four Pion
			fFourPionMask[i] = 0;
			for (Int_t j = 0; j < 9; j++) {
				fFourPionTrack[i][j] = 0;
				fFourPionTPCSigma[i][j] = 0;
				fFourPionTOFSigma[i][j] = 0;
				fFourPionITSSigma[i][j] = 0;
				fFourPionDetMask[i][j] = 0;
				fFourPionTrack_ITSSA[i][j] = 0;
				if ( j < 5) fFourPionBayesProb[i][j] = 0;
			}
		}
		for (Int_t i = 0; i < 4; i++) {
			AliESDtrack *track = fESDEvent->GetTrack(indices.At(i));
			if(!track) {
				PostOutputs();
				return;
			}
			fFourPionTrack[i][0] = track->Px();
			fFourPionTrack[i][1] = track->Py();
			fFourPionTrack[i][2] = track->Pz();
			fFourPionTrack[i][3] = track->E();
			fFourPionTrack[i][4] = track->GetSign();
			//TPC (by L.Goerlich)
			fFourPionTPCSigma[i][0] = track->GetTPCsignal();
			fFourPionTPCSigma[i][1] = (Double_t)(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track) == AliPIDResponse::kDetPidOk);
			fFourPionTPCSigma[i][2] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
			fFourPionTPCSigma[i][3] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
			fFourPionTPCSigma[i][4] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);
			fFourPionTPCSigma[i][5] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
			//TOF (by L.Goerlich)
			fFourPionTOFSigma[i][0] = track->GetTOFsignal();
			fFourPionTOFSigma[i][1] = (Double_t)(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk);
			fFourPionTOFSigma[i][2] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
			fFourPionTOFSigma[i][3] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
			fFourPionTOFSigma[i][4] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);
			//ITS (by L.Goerlich)
			fFourPionITSSigma[i][0] = track->GetITSsignal();
			fFourPionITSSigma[i][1] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kPion);
			fFourPionITSSigma[i][2] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kPion);
			fFourPionITSSigma[i][3] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kProton);
			fFourPionITSSigma[i][4] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kKaon);
			fFourPionITSSigma[i][5] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kElectron);
			fFourPionITSSigma[i][6] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kKaon);
			fFourPionITSSigma[i][7] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kProton);
			fFourPionITSSigma[i][8] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron);
			//Bayesian,combined(UInt_t)
			fFourPionMask[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fFourPionBayesProb[i]);
			fFourPionDetMask[i][0] = (fFourPionMask[i] & AliPIDResponse::kDetITS);
			fFourPionDetMask[i][1] = (fFourPionMask[i] & AliPIDResponse::kDetTPC);
			fFourPionDetMask[i][2] = (fFourPionMask[i] & AliPIDResponse::kDetTRD);
			fFourPionDetMask[i][3] = (fFourPionMask[i] & AliPIDResponse::kDetTOF);
		}
		fCheckFourPion = kTRUE;
	}
	//-------------------------------------------------------------------------

	/*
	//For ITSSA track----------------------------------------------------------
	if (fCheck2tracks_ITSSA) {
		for (Int_t j = 0; j < 10; j++) {
			for (Int_t i = 0; i < 2; i++) {
				fTwoPionTrack_ITSSA[i][j] = 0;
			}
		}
		Int_t tmp_i = 0;
		for (Int_t i = 0; i < tmp_Ntrk; i++) {
			AliESDtrack *track = fESDEvent->GetTrack(i);
			if(!track) {
				continue;
			}
			if(!fTrackCuts_ITSSA->AcceptTrack(track)) continue;

			fTwoPionTrack_ITSSA[tmp_i][0] = track->Px();
			fTwoPionTrack_ITSSA[tmp_i][1] = track->Py();
			fTwoPionTrack_ITSSA[tmp_i][2] = track->Pz();
			fTwoPionTrack_ITSSA[tmp_i][3] = track->E();
			fTwoPionTrack_ITSSA[tmp_i][4] = track->GetSign();
			tmp_i++;
		}
		fCheckTwoPion_ITSSA = kTRUE;
	}
	//-------------------------------------------------------------------------
	*/
	
	indices = 0x0;

	if(!fIsMC) fTree->Fill();

	PostOutputs();
	return;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::PostOutputs()
{
	if (fIsMC) fTree->Fill();
	PostData(1,fTree);
	PostData(2,fList);
	return;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::DoSPDCheck(
		const AliESDEvent *ESDEvent,
		TH1D *hCluster,
		TH1D *hHardware,
		TH1D *hwcSPD
		)
{
	AliTriggerAnalysis triggerAnalysis;
	const Int_t fastOr = triggerAnalysis.SPDFiredChips(ESDEvent, 1); //From hardware bits
	if (hHardware) hHardware->Fill(fastOr);

	//Which chip is fired?
	if (hwcSPD) {
		const AliMultiplicity *mult = ESDEvent->GetMultiplicity();

		for (Int_t iChip = 0; iChip < 1200; iChip++) {
			if (mult->TestFastOrFiredChips(iChip)) {
				hwcSPD->Fill((Double_t)iChip);
			}
		}
	}

	const Int_t fastOr2 = triggerAnalysis.SPDFiredChips(ESDEvent, 0); //From cluster
	if (hCluster) hCluster->Fill(fastOr2);

	return;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::CutEvent(
		const AliESDEvent *ESDEvent,
		TH1 *hpriVtxX, TH1 *hpriVtxY, TH1 *hpriVtxZ)
{
	// Primary vertex cut(10cm)------------------------------------------------
	Bool_t kpr0 = kTRUE;
	const AliESDVertex *vertex = ESDEvent->GetPrimaryVertexTracks();
	if (!vertex) return kFALSE;
	if(vertex->GetNContributors() <1) {
		// SPD vertex
		vertex = ESDEvent->GetPrimaryVertexSPD();
		if(vertex->GetNContributors()<1) {
			// NO VERTEX, SKIP EVENT
			kpr0 = kFALSE;
		}
	}
	const Bool_t kpriv = kpr0 && (fabs(vertex->GetZ()) < 15.);
	if(!kpriv) return kFALSE;

	if(hpriVtxX) hpriVtxX->Fill(vertex->GetX());
	if(hpriVtxY) hpriVtxY->Fill(vertex->GetY());
	if(hpriVtxZ) hpriVtxZ->Fill(vertex->GetZ());

	fVertex[0] = vertex->GetX();
	fVertex[1] = vertex->GetY();
	fVertex[2] = vertex->GetZ();
	//-------------------------------------------------------------------------
	
	return kTRUE;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::CheckInput()
{
	//General protection------------------------------------------------------
	if (const AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler)) {
	  fESDEvent = (AliESDEvent*)esdH->GetEvent();
	}

	if(!fESDEvent){
		//printf("AliAnalysisTaskCDPWA No valid event\n");
		return kFALSE;
	}

	fPIDResponse = (AliPIDResponse*)fInputHandler->GetPIDResponse();
	if(!fPIDResponse) {
		//printf("PIDResponse is not working!!\n");
		return kFALSE;
	}
	//-------------------------------------------------------------------------

	//Check Magnetic Field-----------------------------------------------------
	if(TMath::Abs(fESDEvent->GetMagneticField()) < 1) {
		//printf("AliAnalysisTaskCDPWA strange Bfield! %f\n", fESDEvent->GetMagneticField());
		return kFALSE;
	}
	//-------------------------------------------------------------------------

	//Get MC event-------------------------------------------------------------
	fMCEvent = MCEvent();
	if (fMCEvent) fIsMC = kTRUE;
	//-------------------------------------------------------------------------

	return kTRUE;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::CheckOnlineTrigger(
		const AliESDEvent *ESDevent
		)
{
	//Fill histgrams for trigger
	if (ESDevent->IsTriggerClassFired("CINT5-B-NOPF-ALLNOTRD")) fHistTrigger->Fill(kTri_CINT5);
	if (ESDevent->IsTriggerClassFired("CINT7-B-NOPF-ALLNOTRD")) fHistTrigger->Fill(kTri_CINT7);
	if (ESDevent->IsTriggerClassFired("CINT10-B-NOPF-ALLNOTRD")) fHistTrigger->Fill(kTri_CINT10);
	if (ESDevent->IsTriggerClassFired("C0SMB-B-NOPF-ALLNOTRD")) fHistTrigger->Fill(kTri_C0SMB);
	if (ESDevent->IsTriggerClassFired("CADAND-B-NOPF-ALLNOTRD")) fHistTrigger->Fill(kTri_CADAND);
	if (ESDevent->IsTriggerClassFired("CDG6-B-NOPF-CENTNOTRD")) fHistTrigger->Fill(kTri_CDG6);
	if (ESDevent->IsTriggerClassFired("CDG6-B-SPD2-CENTNOTRD")) fHistTrigger->Fill(kTri_CDG6_SPD2);
	if (ESDevent->IsTriggerClassFired("CDG7-B-SPD2-CENTNOTRD")) fHistTrigger->Fill(kTri_CDG7_SPD2);

	//Check that this event is MB_OR online triggered
	if (ESDevent->IsTriggerClassFired("CINT10-B-NOPF-ALLNOTRD") ||
			ESDevent->IsTriggerClassFired("C0SMB-B-NOPF-ALLNOTRD")) return kTRUE;
	else {
		//printf("AliAnalysisTaskCDTest1::CheckOnlineTrigger not passed!");
		return kFALSE;
	}
	
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::DoTimeV0AD(
		const AliESDEvent *ESDevent,
		TH1D *v0ahist, TH1D *v0chist,
		TH1D *adahist, TH1D *adchist
		)
{
	// Fill time measurement before doing offline measurements
	// V0
	AliESDVZERO *esdv0 = ESDevent->GetVZEROData();
	if (!esdv0) return;

	// Fill only BB-triggered event
	Bool_t BB_V0A = kFALSE;
	Bool_t BB_V0C = kFALSE;

	for (Int_t i = 0; i < 32; i++) {
		if (esdv0->GetBBFlag(i)) BB_V0C = kTRUE;
		if (esdv0->GetBBFlag(i+32)) BB_V0A = kTRUE;
	}

	if (BB_V0A) v0ahist->Fill(esdv0->GetV0ATime());
	if (BB_V0C) v0chist->Fill(esdv0->GetV0CTime());

	// AD
	AliESDAD *esdad = (AliESDAD*)ESDevent->GetADData();
	if (!esdad) return;

	// Fill only BB-triggered event
	Bool_t BB_ADA = kFALSE;
	Bool_t BB_ADC = kFALSE;

	for (Int_t i = 0; i < 4; i++) {//Coincidence
		if (esdad->GetBBFlag(i) && esdad->GetBBFlag(i+4)) BB_ADC = kTRUE;
		if (esdad->GetBBFlag(i+8) && esdad->GetBBFlag(i+12)) BB_ADA = kTRUE;
	}

	if (BB_ADA) adahist->Fill(esdad->GetADATime());
	if (BB_ADC) adchist->Fill(esdad->GetADCTime());

	return;
	
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::SPDLoc2Glo(const Int_t id, const Double_t *loc,
		Double_t *glo)
{
	//
	//SPDLoc2Glo, do not touch
	//

	static TGeoHMatrix mat;
	Int_t vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(id);
	if (vid<0) {
		//printf("AliCDMesonUtils Did not find module with such ID %d\n",id);
		return kFALSE;
	}
	AliITSAlignMille2Module::SensVolMatrix(vid,&mat);
	mat.LocalToMaster(loc,glo);
	return kTRUE;
}
//------------------------------------------------------------------------------
Int_t AliAnalysisTaskCDPWA::CheckChipEta(const Int_t chipKey, const TString scut,const Double_t vtxPos[],TH2 *hitMapSPDinner, TH2 *hitMapSPDouter)
{
	//
	//CheckChipEta
	//

	// retrieves the position in eta for a given chip and applies the cut
	// results:
	// 0 <= out of range
	// -1 <= negative pseudo-rapidity position, in range (C-Side)
	// 1 <= positive pseudo-rapidity position, in range (A-Side)
	//
	// scut: "[0.9" or "]0.9", only 3 digits for the value!!


	const Bool_t kincl = (scut[0] == '[');
	const TString cutval = scut(1,3);
	const Double_t etacut = fabs(cutval.Atof());

	//no eta cut, save time
	if(kincl && etacut>=2)
		return kTRUE;

	Int_t etaside = 1;
	//------------------------------- NOT TO TOUCH ------------------------>>
	UInt_t module=999, offchip=999;
	AliSPDUtils::GetOfflineFromOfflineChipKey(chipKey,module,offchip);
	UInt_t hs = AliSPDUtils::GetOnlineHSFromOffline(module);
	if(hs<2) offchip = 4 - offchip; // inversion  in the inner layer...

	const Int_t col[]={
		hs<2? 0 : 31,
		hs<2? 31 : 0,
		hs<2? 31 : 0,
		hs<2? 0 : 31};
	const Int_t aa[]={0, 0, 255, 255};
	const AliITSsegmentationSPD seg;

	for(Int_t ic=0; ic<4; ic++){
		Float_t localchip[3]={0.,0.,0.};
		seg.DetToLocal(aa[ic],col[ic]+32*offchip,localchip[0],localchip[2]);
		// local coordinate of the chip center
		//printf("local coordinates %d %d: %f %f \n",chipKey, ic, localchip[0],localchip[2]);
		const Double_t local[3] = {localchip[0],localchip[1],localchip[2]};
		Double_t glochip[3]={0.,0.,0.};
		if(!SPDLoc2Glo(module,local,glochip)){
			return kFALSE;
		}

		//-------------------------------------------------------------------<<

		const TVector3 pos(glochip[0]-vtxPos[0], glochip[1]-vtxPos[1],
				glochip[2]-vtxPos[2]);
		//pos.Print();

		if (chipKey < 400) { // inner SPD layer
			if (hitMapSPDinner) {
				Double_t phi = pos.Phi(); // output in the range -Pi +Pi
				if (phi < 0.) phi += TMath::TwoPi(); // remap to the interval [0, TwoPi)
				const Double_t eta = pos.Eta();
				hitMapSPDinner->Fill(eta, phi);
			}
		}
		else {
			if (hitMapSPDouter) { // outer SPD layer
				Double_t phi = pos.Phi(); // output in the range -Pi +Pi
				if (phi < 0.) phi += TMath::TwoPi(); // remap to the interval [0, TwoPi)
				const Double_t eta = pos.Eta();
				hitMapSPDouter->Fill(eta, phi);
			}
		}

		if( kincl && fabs(pos.Eta()) > etacut)
			return kFALSE;

		if(!kincl){
			if(fabs(pos.Eta()) < etacut)
				return kFALSE;
			else if(pos.Eta()<0)
				etaside = -1;
			else
				etaside = 1;
		}
	}

	return etaside;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::CheckV0Hit(const AliESDEvent *ESDEvent, TH1D* hitV0A, TH1D* hitV0C)
{
	AliTriggerAnalysis triggerAnalysis;
	const Bool_t khw = kFALSE;
	const Bool_t v0A = (triggerAnalysis.V0Trigger(ESDEvent, AliTriggerAnalysis::kASide, khw) == AliTriggerAnalysis::kV0BB);
	const Bool_t v0C =
		(triggerAnalysis.V0Trigger(ESDEvent, AliTriggerAnalysis::kCSide, khw) ==
		 AliTriggerAnalysis::kV0BB);

	if (v0A || v0C) {
		hitV0A->Fill(1);
		hitV0C->Fill(1);
		return kFALSE;
	}

	hitV0A->Fill(0);
	hitV0C->Fill(0);
	return kTRUE;

}
//______________________________________________________________________________
Int_t AliAnalysisTaskCDPWA::GetFastORmultiplicity(const AliESDEvent* ESDEvent)
{
	// determine the number of fired fastOR chips in both layers within
	// -0.9 < eta < 0.9

	const AliMultiplicity *mult = ESDEvent->GetMultiplicity();

	// position of the primary vertex
	Double_t tmp[3] = { 0., 0., 0. };
	ESDEvent->GetPrimaryVertex()->GetXYZ(tmp);
	Double_t vtxPos[3] = { tmp[0], tmp[1], tmp[2] };
	Int_t multiplicity = 0;
	for (Int_t iChipKey=0; iChipKey < 1200; iChipKey++) {
		if(mult->TestFastOrFiredChips(iChipKey)){
			//here you check if the FastOr bit is 1 or 0
			const Int_t iseta = CheckChipEta(iChipKey, "[0.9]", vtxPos, 0x0, 0x0);
			if(iseta==0)
				continue;
			else
				++multiplicity;
		}
	}
	return multiplicity;
}
//______________________________________________________________________________
Int_t AliAnalysisTaskCDPWA::DetermineGap(const AliESDEvent *esd,const Bool_t wCent, const Int_t type)
{
	// determines the gap configuration for all gap tagging detectors based on the
	// data set which is available

	Bool_t BB_V0A = kFALSE;
	Bool_t BB_V0C = kFALSE;
	Bool_t BB_ADA = kFALSE;
	Bool_t BB_ADC = kFALSE;

	//Read V0
	AliESDVZERO *v0 = NULL;
	v0 = (AliESDVZERO*)(esd->GetVZEROData());
	if(!v0) return -999;
	for (Int_t i = 0; i < 32; i++) {
		if (v0->GetBBFlag(i)) BB_V0C = kTRUE;
		if (v0->GetBBFlag(i+32)) BB_V0A = kTRUE;
	}

	//Read AD
	AliESDAD *ad = NULL;
	ad = (AliESDAD*)(esd->GetADData());
	if(!ad) return -999;
	for (Int_t i = 0; i < 4; i++) {
		if (ad->GetBBFlag(i) && ad->GetBBFlag(i+4)) BB_ADC=kTRUE;
		if (ad->GetBBFlag(i+8) && ad->GetBBFlag(i+12)) BB_ADA=kTRUE;
	}
	const AliMultiplicity *mult = esd->GetMultiplicity();
	Int_t nChips = 0;
	Int_t firstChip = 0;
	Int_t lastChip = 1200;
	Bool_t centAct = kFALSE;
	for (Int_t i=firstChip; i<lastChip; i++) {
		if (mult->TestFastOrFiredChips(i)) centAct=kTRUE;
	}
	if (wCent) {
		if (type ==0) {
			if (!BB_V0A && !BB_V0C && centAct) return kBinDG;
			else if (!BB_V0A && BB_V0C && !centAct) return kBinGA;
			else if (!BB_V0C && BB_V0A && !centAct) return kBinGC;
			else return kBinNG;
		}
		if (type == 1) {//only AD
			if (!BB_ADA && !BB_ADC && centAct) return kBinDG;
			else if (!BB_ADA && BB_ADC && !centAct) return kBinGA;
			else if (!BB_ADC && BB_ADA && !centAct) return kBinGC;
			else return kBinNG;
		}
		if (type == 2) {//V0 and AD
			if (!BB_ADA && !BB_V0A && !BB_V0C && !BB_ADC) return kBinDG;
			else if (!BB_ADA && !BB_V0A && (BB_ADC || BB_V0C)) return kBinGA;
			else if (!BB_ADC && !BB_V0C && (BB_ADA || BB_V0A)) return kBinGA;
			else return kBinNG;
		}
	}
	else {//No Central Activity
		if (type == 0) {//only V0
			if (!BB_V0A && !BB_V0C) return kBinDG;
			else if (!BB_V0A && BB_V0C) return kBinGA;
			else if (!BB_V0C && BB_V0A) return kBinGC;
			else return kBinNG;
		}
		if (type == 1) {//only AD
			if (!BB_ADA && !BB_ADC) return kBinDG;
			else if (!BB_ADA && BB_ADC) return kBinGA;
			else if (!BB_ADC && BB_ADA) return kBinGC;
			else return kBinNG;
		}
		if (type == 2) {//V0 and AD
			if (!BB_ADA && !BB_V0A && !BB_V0C && !BB_ADC) return kBinDG;
			else if (!BB_ADA && !BB_V0A && (BB_ADC || BB_V0C)) return kBinGA;
			else if (!BB_ADC && !BB_V0C && (BB_ADA || BB_V0A)) return kBinGA;
			else return kBinNG;
		}
	}
	return kBinNG;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::DetermineMCprocessType()
{

	// Get MC information------------------------------------------------------
	fMCprocess = -1; //detailed MC sub process information
	fMCprocessType = kBinND; // ND is default, also for data

	if (fMCEvent) {
		AliGenEventHeader* header = fMCEvent->GenEventHeader();
		if (!header) return;
		// Pythia6
		if (TString(header->IsA()->GetName()) == "AliGenPythiaEventHeader") {
			fMCprocess = ((AliGenPythiaEventHeader*)header)->ProcessType();
			switch(fMCprocess) {
				case 92: fMCprocessType = kBinSD1; break;
				case 93: fMCprocessType = kBinSD2; break;
				case 94: fMCprocessType = kBinDD; break;
				default: fMCprocessType = AliAnalysisTaskCDPWA::kBinND; break;
			}
		}
		// Phojet
		else if (TString(header->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
			fMCprocess = ((AliGenDPMjetEventHeader*)header)->ProcessType();
			switch(fMCprocess) {
				case 5: fMCprocessType = kBinSD1; break;
				case 6: fMCprocessType = kBinSD2; break;
				case 7: fMCprocessType = kBinDD; break;
				case 4: fMCprocessType = kBinCD; break;
				default: fMCprocessType = kBinND; break;
			}
		}
	}
	//-------------------------------------------------------------------------
}
//______________________________________________________________________________
Int_t AliAnalysisTaskCDPWA::DoMCTruth()
{
	if (!fMCEvent) return -1;
	AliStack* stack = fMCEvent->Stack();
	if (!stack) return -1;

	// Multiplicity ----------------------------------------------------------
	// determine number of charged physical primaries on the stack
	Int_t nPhysicalPrimaries = 0;
	Int_t nProton = 0;
	Int_t nPrimaries = stack->GetNprimary();
	TParticle *part = NULL;
	for (Int_t iTracks = 0; iTracks < nPrimaries; ++iTracks) {
		part = stack->Particle(iTracks);
		if (!part) continue;
		if (stack->IsPhysicalPrimary(iTracks) && (part->GetPDG()->Charge() != 0.) && part->GetStatusCode() == 1 && TMath::Abs(((Int_t)(part->GetPdgCode()))) == 211) {
			nPhysicalPrimaries++;
		}
		if (stack->IsPhysicalPrimary(iTracks) && (part->GetPDG()->Charge() != 0.) && part->GetStatusCode() == 1 && TMath::Abs(((Int_t)(part->GetPdgCode()))) == 2212) {
			nProton++;
		}
	}

	// Sergey's MC test(2 pion, proton generated)
	if (nPhysicalPrimaries !=2 || nProton !=2) return -1;
	//-------------------------------------------------------------------------
	
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
	part = NULL;
	for (Int_t iTracks = 0; iTracks <  nPrimaries; ++iTracks) {
		part = (TParticle*)stack->Particle(iTracks);
		if (!part) continue;
		if (stack->IsPhysicalPrimary(iTracks) && (part->GetPDG()->Charge() != 0.) && part->GetStatusCode() == 1 && TMath::Abs(((Int_t)(part->GetPdgCode()))) == 211) {
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
	//-------------------------------------------------------------------------

	//Store Generated proton and pion momentum in the TTree--------------------
	TParticle *part1 = NULL;
	Int_t n_proton = 0;
	Int_t n_pion = 0;
	for (Int_t j = 0; j < nPrimaries; j++) {
		part1 = stack->Particle(j);
		if (!part1) continue;
		if (stack->IsPhysicalPrimary(j) && (part1->GetPDG()->Charge() != 0.) && part1->GetStatusCode() == 1 && TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 2212) {
			fMCGenProtonTrack[n_proton][0] = part1->Px();
			fMCGenProtonTrack[n_proton][1] = part1->Py();
			fMCGenProtonTrack[n_proton][2] = part1->Pz();
			fMCGenProtonTrack[n_proton][3] = part1->Energy();
			fMCGenProtonTrack[n_proton][4] = part1->GetPdgCode();
			n_proton++;
			continue;
		}
		if (stack->IsPhysicalPrimary(j) && (part1->GetPDG()->Charge() != 0.) && part1->GetStatusCode() == 1 && TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 211) {
			fMCGenPionTrack[n_pion][0] = part1->Px();
			fMCGenPionTrack[n_pion][1] = part1->Py();
			fMCGenPionTrack[n_pion][2] = part1->Pz();
			fMCGenPionTrack[n_pion][3] = part1->Energy();
			fMCGenPionTrack[n_pion][4] = part1->GetPdgCode();
			n_pion++;
			continue;
		}
	}
	//-------------------------------------------------------------------------
	return nPhysicalPrimaries;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::Terminate(Option_t *)
{
}
//______________________________________________________________________________
