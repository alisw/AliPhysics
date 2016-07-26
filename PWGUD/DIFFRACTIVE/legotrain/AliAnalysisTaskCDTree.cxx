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
#include "TGeoManager.h"
#include "AliPID.h"
#include "AliPIDCombined.h"

#include "AliMultiplicitySelectionCP.h"
#include "AliAnalysisTaskCDTree.h"

//______________________________________________________________________________
AliAnalysisTaskCDTree::AliAnalysisTaskCDTree(const char* name):
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
	, fGapInformation(0)
	, fTree(0x0)
	, fCheckTwoPion(0)
	, fCheckFourPion(0)
	, fCheckV0FMD(0)
	, fCheckTwoPion_ITSSA(0)
	, fIsMC(0)
	, fList(0x0)
	, fHistEvent(0x0)
	, fHistEventProcesses(0x0)
	, fHistPrimVtxX(0x0)
	, fHistPrimVtxY(0x0)
	, fHistPrimVtxZ(0x0)
	, fHitV0A(0x0)
	, fHitV0C(0x0)
	, fMultRegionsMC(0x0)
	, fHistSPDFiredChips(0x0)
	, fRunVsMBOR(0x0)
	, fRunVsMBAND(0x0)
	, fRunVsDG(0x0)
	, fRunVsDG_wFMD(0x0)
	, fRunVs2t(0x0)
	, fRunVs2t_ITSSA(0x0)
	, fRunVs4t(0x0)
	, fRunVs4t_ITSSA(0x0)
	, fMultNG(0x0)
	, fMultNG_MS(0x0)
	, fMultDG(0x0)
	, fMultDG_MS(0x0)
	, fMassNG(0x0)
	, fMassNG_st2t(0x0)
	, fMassNG_MS(0x0)
	, fMassDG(0x0)
	, fMassDG_st2t(0x0)
	, fMassDG_MS(0x0)
	, fTrackCutsInfo(0x0)
	, fSPDTrkvsCls_bf(0x0)
	, fSPDTrkvsCls_af(0x0)
	, fPrimaries(0x0)
	, fRunNumber(0)
	, fPeriod(0)
{
	//
	// standard constructor (the one which should be used)
	//
	// slot in TaskSE must start from 1

	for (Int_t i = 0; i < 3; i++) {
		fVertex[i] = 0.;
	}
	
	// Two-track
	for (Int_t i = 0; i < 2; i++) {
		fTwoPionMask_TPC[i] = 0;
		fTwoPionMask_TOF[i] = 0;
		fTwoPionMask_ITS[i] = 0;
		fTwoPionMask_TRD[i] = 0;
		fTwoPionMask_tot[i] = 0;
		fTwoPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fTwoPionTrack[i][j] = 0.;
			fTwoPionTPCSigma[i][j] = 0.;
			fTwoPionTOFSigma[i][j] = 0.;
			fTwoPionITSSigma[i][j] = 0.;
			fTwoPionTrack_ITSSA[i][j] = 0.;
			if (j < 5) {
				fTwoPionBayesProb_TPC[i][j] = 0.;
				fTwoPionBayesProb_TOF[i][j] = 0.;
				fTwoPionBayesProb_ITS[i][j] = 0.;
				fTwoPionBayesProb_TRD[i][j] = 0.;
				fTwoPionBayesProb_tot[i][j] = 0.;
			}
		}
	}
	// Four-track
	for (Int_t i = 0; i < 4; i++) {
		fFourPionMask_TPC[i] = 0;
		fFourPionMask_TOF[i] = 0;
		fFourPionMask_ITS[i] = 0;
		fFourPionMask_TRD[i] = 0;
		fFourPionMask_tot[i] = 0;
		fFourPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fFourPionTrack[i][j] = 0.;
			fFourPionTPCSigma[i][j] = 0.;
			fFourPionTOFSigma[i][j] = 0.;
			fFourPionITSSigma[i][j] = 0.;
			fFourPionTrack_ITSSA[i][j] = 0.;
			if (j < 5) {
				fFourPionBayesProb_TPC[i][j] = 0.;
				fFourPionBayesProb_TOF[i][j] = 0.;
				fFourPionBayesProb_ITS[i][j] = 0.;
				fFourPionBayesProb_TRD[i][j] = 0.;
				fFourPionBayesProb_tot[i][j] = 0.;
			}
		}
	}

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
AliAnalysisTaskCDTree::AliAnalysisTaskCDTree():
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
	, fGapInformation(0)
	, fTree(0x0)
	, fCheckTwoPion(0)
	, fCheckFourPion(0)
	, fCheckV0FMD(0)
	, fCheckTwoPion_ITSSA(0)
	, fIsMC(0)
	, fList(0x0)
	, fHistEvent(0x0)
	, fHistEventProcesses(0x0)
	, fHistPrimVtxX(0x0)
	, fHistPrimVtxY(0x0)
	, fHistPrimVtxZ(0x0)
	, fHitV0A(0x0)
	, fHitV0C(0x0)
	, fMultRegionsMC(0x0)
	, fHistSPDFiredChips(0x0)
	, fRunVsMBOR(0x0)
	, fRunVsMBAND(0x0)
	, fRunVsDG(0x0)
	, fRunVsDG_wFMD(0x0)
	, fRunVs2t(0x0)
	, fRunVs2t_ITSSA(0x0)
	, fRunVs4t(0x0)
	, fRunVs4t_ITSSA(0x0)
	, fMultNG(0x0)
	, fMultNG_MS(0x0)
	, fMultDG(0x0)
	, fMultDG_MS(0x0)
	, fMassNG(0x0)
	, fMassNG_st2t(0x0)
	, fMassNG_MS(0x0)
	, fMassDG(0x0)
	, fMassDG_st2t(0x0)
	, fMassDG_MS(0x0)
	, fTrackCutsInfo(0x0)
	, fSPDTrkvsCls_bf(0x0)
	, fSPDTrkvsCls_af(0x0)
	, fPrimaries(0x0)
	, fRunNumber(0)
	, fPeriod(0)
{
	for (Int_t i = 0; i < 3; i++) {
		fVertex[i] = 0.;
	}
	// Two-track
	for (Int_t i = 0; i < 2; i++) {
		fTwoPionMask_TPC[i] = 0;
		fTwoPionMask_TOF[i] = 0;
		fTwoPionMask_ITS[i] = 0;
		fTwoPionMask_TRD[i] = 0;
		fTwoPionMask_tot[i] = 0;
		fTwoPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fTwoPionTrack[i][j] = 0.;
			fTwoPionTPCSigma[i][j] = 0.;
			fTwoPionTOFSigma[i][j] = 0.;
			fTwoPionITSSigma[i][j] = 0.;
			fTwoPionTrack_ITSSA[i][j] = 0.;
			if (j < 5) {
				fTwoPionBayesProb_TPC[i][j] = 0.;
				fTwoPionBayesProb_TOF[i][j] = 0.;
				fTwoPionBayesProb_ITS[i][j] = 0.;
				fTwoPionBayesProb_TRD[i][j] = 0.;
				fTwoPionBayesProb_tot[i][j] = 0.;
			}
		}
	}
	// Four-track
	for (Int_t i = 0; i < 4; i++) {
		fFourPionMask_TPC[i] = 0;
		fFourPionMask_TOF[i] = 0;
		fFourPionMask_ITS[i] = 0;
		fFourPionMask_TRD[i] = 0;
		fFourPionMask_tot[i] = 0;
		fFourPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fFourPionTrack[i][j] = 0.;
			fFourPionTPCSigma[i][j] = 0.;
			fFourPionTOFSigma[i][j] = 0.;
			fFourPionITSSigma[i][j] = 0.;
			fFourPionTrack_ITSSA[i][j] = 0.;
			if (j < 5) {
				fFourPionBayesProb_TPC[i][j] = 0.;
				fFourPionBayesProb_TOF[i][j] = 0.;
				fFourPionBayesProb_ITS[i][j] = 0.;
				fFourPionBayesProb_TRD[i][j] = 0.;
				fFourPionBayesProb_tot[i][j] = 0.;
			}
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
AliAnalysisTaskCDTree::~AliAnalysisTaskCDTree()
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
	if (fHistEventProcesses) delete fHistEventProcesses;
	if (fHistPrimVtxX) delete fHistPrimVtxX;
	if (fHistPrimVtxY) delete fHistPrimVtxY;
	if (fHistPrimVtxZ) delete fHistPrimVtxZ;
	if (fHitV0A) delete fHitV0A;
	if (fHitV0C) delete fHitV0C;
	if (fMultRegionsMC) delete fMultRegionsMC;
	if (fHistSPDFiredChips) delete fHistSPDFiredChips;
	if (fRunVsMBOR) delete fRunVsMBOR;
	if (fRunVsMBAND) delete fRunVsMBAND;
	if (fRunVsDG) delete fRunVsDG;
	if (fRunVsDG_wFMD) delete fRunVsDG_wFMD;
	if (fRunVs2t) delete fRunVs2t;
	if (fRunVs2t_ITSSA) delete fRunVs2t_ITSSA;
	if (fRunVs4t) delete fRunVs4t;
	if (fRunVs4t_ITSSA) delete fRunVs4t_ITSSA;
	if (fMultNG) delete fMultNG;
	if (fMultNG_MS) delete fMultNG_MS;
	if (fMultDG) delete fMultDG;
	if (fMultDG_MS) delete fMultDG_MS;
	if (fMassNG) delete fMassNG;
	if (fMassNG_st2t) delete fMassNG_st2t;
	if (fMassNG_MS) delete fMassNG_MS;
	if (fMassDG) delete fMassDG;
	if (fMassDG_st2t) delete fMassDG_st2t;
	if (fMassDG_MS) delete fMassDG_MS;
	if (fTrackCutsInfo) delete fTrackCutsInfo;
	if (fSPDTrkvsCls_bf) delete fSPDTrkvsCls_bf;
	if (fSPDTrkvsCls_af) delete fSPDTrkvsCls_af;
	if (fPrimaries) delete fPrimaries;
}
//______________________________________________________________________________
void AliAnalysisTaskCDTree::UserCreateOutputObjects()
{
	// TTree + TList should be defined-----------------------------------------
	OpenFile(1);
	fTree = new TTree("tree1","PWAtree");
	//For Data + Rec
	fTree->Branch("CheckTwoPion",&fCheckTwoPion);
	fTree->Branch("CheckFourPion",&fCheckFourPion);
	fTree->Branch("CheckV0FMD",&fCheckV0FMD);
	fTree->Branch("CheckTwoPion_ITSSA",&fCheckTwoPion_ITSSA);
	for (Int_t i = 0; i < 3; i++) {
		fTree->Branch(Form("Vertex_%d",i),&fVertex[i]);
	}
	//2 track
	for (Int_t i = 0; i < 2; i++) {
		fTree->Branch(Form("TwoPionMask_TPC_%d",i),&fTwoPionMask_TPC[i]);
		fTree->Branch(Form("TwoPionMask_TOF_%d",i),&fTwoPionMask_TOF[i]);
		fTree->Branch(Form("TwoPionMask_ITS_%d",i),&fTwoPionMask_ITS[i]);
		fTree->Branch(Form("TwoPionMask_TRD_%d",i),&fTwoPionMask_TRD[i]);
		fTree->Branch(Form("TwoPionMask_tot_%d",i),&fTwoPionMask_tot[i]);
		fTree->Branch(Form("TwoPionDetMask_tot_%d",i),&fTwoPionDetMask_tot[i]);
		for (Int_t j = 0; j < 9; j++) {
			fTree->Branch(Form("TwoPionTrack_%d_%d",i,j),&fTwoPionTrack[i][j]);
			fTree->Branch(Form("TwoPionTPCSigma_%d_%d",i,j),&fTwoPionTPCSigma[i][j]);
			fTree->Branch(Form("TwoPionTOFSigma_%d_%d",i,j),&fTwoPionTOFSigma[i][j]);
			fTree->Branch(Form("TwoPionITSSigma_%d_%d",i,j),&fTwoPionITSSigma[i][j]);
			fTree->Branch(Form("TwoPionTrack_ITSSA_%d_%d",i,j),&fTwoPionTrack_ITSSA[i][j]);
			if (j < 5) {
				fTree->Branch(Form("TwoPionBayesProb_TPC_%d_%d",i,j),&fTwoPionBayesProb_TPC[i][j]);
				fTree->Branch(Form("TwoPionBayesProb_TOF_%d_%d",i,j),&fTwoPionBayesProb_TOF[i][j]);
				fTree->Branch(Form("TwoPionBayesProb_ITS_%d_%d",i,j),&fTwoPionBayesProb_ITS[i][j]);
				fTree->Branch(Form("TwoPionBayesProb_TRD_%d_%d",i,j),&fTwoPionBayesProb_TRD[i][j]);
				fTree->Branch(Form("TwoPionBayesProb_tot_%d_%d",i,j),&fTwoPionBayesProb_tot[i][j]);
			}
		}
	}
	//2 track
	for (Int_t i = 0; i < 4; i++) {
		fTree->Branch(Form("FourPionMask_TPC_%d",i),&fFourPionMask_TPC[i]);
		fTree->Branch(Form("FourPionMask_TOF_%d",i),&fFourPionMask_TOF[i]);
		fTree->Branch(Form("FourPionMask_ITS_%d",i),&fFourPionMask_ITS[i]);
		fTree->Branch(Form("FourPionMask_TRD_%d",i),&fFourPionMask_TRD[i]);
		fTree->Branch(Form("FourPionMask_tot_%d",i),&fFourPionMask_tot[i]);
		fTree->Branch(Form("FourPionDetMask_tot_%d",i),&fFourPionDetMask_tot[i]);
		for (Int_t j = 0; j < 9; j++) {
			fTree->Branch(Form("FourPionTrack_%d_%d",i,j),&fFourPionTrack[i][j]);
			fTree->Branch(Form("FourPionTPCSigma_%d_%d",i,j),&fFourPionTPCSigma[i][j]);
			fTree->Branch(Form("FourPionTOFSigma_%d_%d",i,j),&fFourPionTOFSigma[i][j]);
			fTree->Branch(Form("FourPionITSSigma_%d_%d",i,j),&fFourPionITSSigma[i][j]);
			fTree->Branch(Form("FourPionTrack_ITSSA_%d_%d",i,j),&fFourPionTrack_ITSSA[i][j]);
			if (j < 5) {
				fTree->Branch(Form("FourPionBayesProb_TPC_%d_%d",i,j),&fFourPionBayesProb_TPC[i][j]);
				fTree->Branch(Form("FourPionBayesProb_TOF_%d_%d",i,j),&fFourPionBayesProb_TOF[i][j]);
				fTree->Branch(Form("FourPionBayesProb_ITS_%d_%d",i,j),&fFourPionBayesProb_ITS[i][j]);
				fTree->Branch(Form("FourPionBayesProb_TRD_%d_%d",i,j),&fFourPionBayesProb_TRD[i][j]);
				fTree->Branch(Form("FourPionBayesProb_tot_%d_%d",i,j),&fFourPionBayesProb_tot[i][j]);
			}
		}
	}

	//For MC
	for (Int_t i = 0; i < 5; i++) {
		for (Int_t j = 0; j < 2; j++) {
			fTree->Branch(Form("MCGenProtonTrack_%d_%d",j,i),&fMCGenProtonTrack[j][i]);
			fTree->Branch(Form("MCGenPionTrack_%d_%d",j,i),&fMCGenPionTrack[j][i]);
		}
	}
	fTree->Branch("RunNumber",&fRunNumber);
	fTree->Branch("Period",&fPeriod);
	//-------------------------------------------------------------------------

	// TList-------------------------------------------------------------------
	OpenFile(2);
	fList = new TList();
	fList->SetOwner();
	if (!fHistEvent) {
		fHistEvent = new TH1D("fHistEvent","Number of Events for each selection",kAll,0,kAll);
		fHistEvent->GetXaxis()->SetBinLabel(kInput+1,"InputEvent");
		fHistEvent->GetXaxis()->SetBinLabel(kMCCheck+1,"MCEvent");
		fHistEvent->GetXaxis()->SetBinLabel(kOfflineCut+1,"OfflineCut");
		fHistEvent->GetXaxis()->SetBinLabel(kVtxCut+1,"VertexCut");
		fHistEvent->GetXaxis()->SetBinLabel(kPileUpCut+1,"PileupCut");
		fHistEvent->GetXaxis()->SetBinLabel(kClusterCut+1,"ClusterCut");
		fHistEvent->GetXaxis()->SetBinLabel(kMBOR+1,"MBOR");
		fHistEvent->GetXaxis()->SetBinLabel(kMBAND+1,"MBAND");
		fHistEvent->GetXaxis()->SetBinLabel(kNG+1,"NoGap");
		fHistEvent->GetXaxis()->SetBinLabel(kGA+1,"A-side Gap");
		fHistEvent->GetXaxis()->SetBinLabel(kGC+1,"C-side Gap");
		fHistEvent->GetXaxis()->SetBinLabel(kDG+1,"DoubleGap");
		fHistEvent->GetXaxis()->SetBinLabel(kNG_wFMD+1,"NoGap_wFMD");
		fHistEvent->GetXaxis()->SetBinLabel(kGA_wFMD+1,"A-side Gap_wFMD");
		fHistEvent->GetXaxis()->SetBinLabel(kGC_wFMD+1,"C-side Gap_wFMD");
		fHistEvent->GetXaxis()->SetBinLabel(kDG_wFMD+1,"DoubleGap_wFMD");
		fHistEvent->GetXaxis()->SetBinLabel(kCheckDG+1,"CheckDoubleGap");
		fHistEvent->GetXaxis()->SetBinLabel(kCheckV0Hit+1,"CheckV0Hit");
		fHistEvent->GetXaxis()->SetBinLabel(k2Tracks+1,"2 Tracks(DG)");
		fHistEvent->GetXaxis()->SetBinLabel(k4Tracks+1,"4 Tracks(DG)");
		fHistEvent->GetXaxis()->SetBinLabel(k6Tracks+1,"6 Tracks(DG)");
		fList->Add(fHistEvent);
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
		fHistPrimVtxX = new TH1D("fHistPrimVtxX","PrimaryVertex X",1000,-1.5,1.5);
		fList->Add(fHistPrimVtxX);
	}
	if (!fHistPrimVtxY) {
		fHistPrimVtxY = new TH1D("fHistPrimVtxY","PrimaryVertex Y",1000,-1.5,1.5);
		fList->Add(fHistPrimVtxY);
	}
	if (!fHistPrimVtxZ) {
		fHistPrimVtxZ = new TH1D("fHistPrimVtxZ","PrimaryVertex Z",1000,-30,30);
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
	if (!fHistSPDFiredChips) {
		fHistSPDFiredChips = new TH1D("fHistSPDFiredChips","",1200,0,1200);
		fList->Add(fHistSPDFiredChips);
	}
	Int_t firstRun = 114500;
	Int_t lastRun = 178035;
	Int_t totalRun = lastRun-firstRun;

	if (!fRunVsMBOR) {
		fRunVsMBOR = new TH1D("fRunVsMBOR","",totalRun,firstRun,lastRun);
		fList->Add(fRunVsMBOR);
	}
	if (!fRunVsMBAND) {
		fRunVsMBAND = new TH1D("fRunVsMBAND","",totalRun,firstRun,lastRun);
		fList->Add(fRunVsMBAND);
	}
	if (!fRunVsDG) {
		fRunVsDG = new TH1D("fRunVsDG","",totalRun,firstRun,lastRun);
		fList->Add(fRunVsDG);
	}
	if (!fRunVsDG_wFMD) {
		fRunVsDG_wFMD = new TH1D("fRunVsDG_wFMD","",totalRun,firstRun,lastRun);
		fList->Add(fRunVsDG_wFMD);
	}
	if (!fRunVs2t) {
		fRunVs2t = new TH1D("fRunVs2t","",totalRun,firstRun,lastRun);
		fList->Add(fRunVs2t);
	}
	if (!fRunVs2t_ITSSA) {
		fRunVs2t_ITSSA = new TH1D("fRunVs2t_ITSSA","",totalRun,firstRun,lastRun);
		fList->Add(fRunVs2t_ITSSA);
	}
	if (!fRunVs4t) {
		fRunVs4t = new TH1D("fRunVs4t","",totalRun,firstRun,lastRun);
		fList->Add(fRunVs4t);
	}
	if (!fRunVs4t_ITSSA) {
		fRunVs4t_ITSSA = new TH1D("fRunVs4t_ITSSA","",totalRun,firstRun,lastRun);
		fList->Add(fRunVs4t_ITSSA);
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
	if (!fMassNG_st2t) {
		fMassNG_st2t = new TH1D("fMassNG_st2t","",50,0,2);
		fList->Add(fMassNG_st2t);
	}
	if (!fMassNG_MS) {
		fMassNG_MS = new TH1D("fMassNG_MS","",50,0,2);
		fList->Add(fMassNG_MS);
	}
	if (!fMassDG) {
		fMassDG = new TH1D("fMassDG","",50,0,2);
		fList->Add(fMassDG);
	}
	if (!fMassDG_st2t) {
		fMassDG_st2t = new TH1D("fMassDG_st2t","",50,0,2);
		fList->Add(fMassDG_st2t);
	}
	if (!fMassDG_MS) {
		fMassDG_MS = new TH1D("fMassDG_MS","",50,0,2);
		fList->Add(fMassDG_MS);
	}
	if (!fTrackCutsInfo) {
		fTrackCutsInfo = new TH1D("fTrackCutsInfo","",20,-10,10);
		fList->Add(fTrackCutsInfo);
	}
	if (!fSPDTrkvsCls_bf) {
		fSPDTrkvsCls_bf = new TH2D("fSPDTrkvsCls_bf","",200,0,200,1000,0,1000);
		fList->Add(fSPDTrkvsCls_bf);
	}
	if (!fSPDTrkvsCls_af) {
		fSPDTrkvsCls_af = new TH2D("fSPDTrkvsCls_af","",200,0,200,1000,0,1000);
		fList->Add(fSPDTrkvsCls_af);
	}
	if (!fPrimaries) {
		fPrimaries = new TH1D("fPrimaries","",100,0,100);
		fList->Add(fPrimaries);
	}

	// Track Cuts
	fTrackCuts = new AliESDtrackCuts();
	{
		fTrackCuts->GetStandardITSTPCTrackCuts2010(1,1);//0 for d and e
	}
	fTrackCuts_ITSSA = new AliESDtrackCuts();
	{
		fTrackCuts_ITSSA->GetStandardITSSATrackCuts2010(1,0);//0 for no-PID mode
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
void AliAnalysisTaskCDTree::UserExec(Option_t *)
{
	// Initialize all tree varialbles------------------------------------------
	for (Int_t i = 0; i < 3; i++) {
		fVertex[i] = -999.;
	}
	for (Int_t i = 0; i < 2; i++) {
		fTwoPionMask_TPC[i] = 0;
		fTwoPionMask_TOF[i] = 0;
		fTwoPionMask_ITS[i] = 0;
		fTwoPionMask_TRD[i] = 0;
		fTwoPionMask_tot[i] = 0;
		fTwoPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fTwoPionTrack[i][j] = 0.;
			fTwoPionTPCSigma[i][j] = 0.;
			fTwoPionTOFSigma[i][j] = 0.;
			fTwoPionITSSigma[i][j] = 0.;
			fTwoPionTrack_ITSSA[i][j] = 0.;
			if (j < 5) {
				fTwoPionBayesProb_TPC[i][j] = 0.;
				fTwoPionBayesProb_TOF[i][j] = 0.;
				fTwoPionBayesProb_ITS[i][j] = 0.;
				fTwoPionBayesProb_TRD[i][j] = 0.;
				fTwoPionBayesProb_tot[i][j] = 0.;
			}
		}
	}
	// Four-track
	for (Int_t i = 0; i < 4; i++) {
		fFourPionMask_TPC[i] = 0;
		fFourPionMask_TOF[i] = 0;
		fFourPionMask_ITS[i] = 0;
		fFourPionMask_TRD[i] = 0;
		fFourPionMask_tot[i] = 0;
		fFourPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 9; j++) {
			fFourPionTrack[i][j] = 0.;
			fFourPionTPCSigma[i][j] = 0.;
			fFourPionTOFSigma[i][j] = 0.;
			fFourPionITSSigma[i][j] = 0.;
			fFourPionTrack_ITSSA[i][j] = 0.;
			if (j < 5) {
				fFourPionBayesProb_TPC[i][j] = 0.;
				fFourPionBayesProb_TOF[i][j] = 0.;
				fFourPionBayesProb_ITS[i][j] = 0.;
				fFourPionBayesProb_TRD[i][j] = 0.;
				fFourPionBayesProb_tot[i][j] = 0.;
			}
		}
	}
	for (Int_t i = 0; i < 5; i++) {
		for (Int_t j = 0; j < 2; j++) {
			fMCGenProtonTrack[j][i] = 0.;
			fMCGenPionTrack[j][i] = 0.;
		}
	}
	fCheckTwoPion = kFALSE;
	fCheckFourPion = kFALSE;
	fCheckV0FMD = kFALSE;
	fCheckTwoPion_ITSSA = kFALSE;
	fIsMC = kFALSE;
	fGapInformation = 0;
	fRunNumber = -999;
	fPeriod = -999;
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
		fIsMC = kTRUE;
	}
	fHistEvent->Fill(kMCCheck);
	//-------------------------------------------------------------------------

	//Check event selected-----------------------------------------------------
	AliInputEventHandler *inputHandler = (AliInputEventHandler*)
		AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	if (!inputHandler) {
		PostOutputs();
		return;
	}
	if (!inputHandler->IsEventSelected()) {
//		printf("This event is not selected by AliVEvent::kMB");
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kOfflineCut);
	//-------------------------------------------------------------------------

	// EVENT SELECTION---------------------------------------------------------
	Int_t kfo = 0;
	Int_t ninnerp=-999, nouterp=-999;
	Bool_t eventIsValid = CutEvent(fESDEvent, fHistSPDFiredChips,0x0,fHistPrimVtxX, fHistPrimVtxY, fHistPrimVtxZ, fVertex);
	if (!eventIsValid) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kVtxCut); //After Vertex cut
	//-------------------------------------------------------------------------

	// Reject PILE UP ---------------------------------------------------------
	const Bool_t isPileup = fESDEvent->IsPileupFromSPD(2, 0.8, 3., 2., 5.);
	// using only 2 instead of three contributors

	if (isPileup) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kPileUpCut); // After pile up
	fRunNumber = (Int_t)fESDEvent->GetRunNumber();
	if (fRunNumber < 117224) fPeriod = 1;//10b
	else if (fRunNumber < 121041) fPeriod = 2;//10c
	else if (fRunNumber < 126438) fPeriod = 3;//10d
	else if (fRunNumber < 130851) fPeriod = 4;//10e
	else if (fRunNumber < 135032) fPeriod = 5;//10f
	else fPeriod = 6;//12b (8 TeV)
	//-------------------------------------------------------------------------

	//Cluster Cut--------------------------------------------------------------
	const AliMultiplicity *mult = fESDEvent->GetMultiplicity();
	fSPDTrkvsCls_bf->Fill(mult->GetNumberOfTracklets(),fESDEvent->GetNumberOfITSClusters(0)+fESDEvent->GetNumberOfITSClusters(1));
	Double_t cut_slope = 4.;
	Double_t cut_offset = 65.;
	Bool_t IsClusterCut = kFALSE;
	if (fESDEvent->GetNumberOfITSClusters(0) + fESDEvent->GetNumberOfITSClusters(1) <= (cut_offset+cut_slope*mult->GetNumberOfTracklets())) IsClusterCut=kTRUE;
	if (IsClusterCut) fSPDTrkvsCls_af->Fill(mult->GetNumberOfTracklets(),fESDEvent->GetNumberOfITSClusters(0)+fESDEvent->GetNumberOfITSClusters(1));
	if (!IsClusterCut) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kClusterCut);
	//-------------------------------------------------------------------------

	// TRIGGER ANALYSIS -------------------------------------------------------
	AliTriggerAnalysis *fTrigger = new AliTriggerAnalysis;
	Bool_t IsMBOR = (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kMB1)) ? kTRUE : kFALSE;
	Bool_t IsMBAND = (fTrigger->IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0AND)) ? kTRUE : kFALSE;

	if (IsMBOR) {fHistEvent->Fill(kMBOR); fRunVsMBOR->Fill(fESDEvent->GetRunNumber());}
	if (IsMBAND) {fHistEvent->Fill(kMBAND); fRunVsMBAND->Fill(fESDEvent->GetRunNumber());}
	if (!IsMBOR) {
		PostOutputs();
		return;
	}
	//-------------------------------------------------------------------------

	// Determine GAP ----------------------------------------------------------
	Bool_t wCent = kTRUE;
	Bool_t wFMD = kTRUE;
	//FMD First!!
	fGapInformation=DetermineGap(fESDEvent, wCent, wFMD);
	if (fGapInformation==kBinNG) fHistEvent->Fill(kNG_wFMD);
	else if (fGapInformation==kBinGA) fHistEvent->Fill(kGA_wFMD);
	else if (fGapInformation==kBinGC) fHistEvent->Fill(kGC_wFMD);
	else if (fGapInformation==kBinDG) {fHistEvent->Fill(kDG_wFMD); fRunVsDG_wFMD->Fill(fESDEvent->GetRunNumber());}
	if (fGapInformation==kBinDG) fCheckV0FMD = kTRUE;

	//without FMD >> TTree
	wFMD = kFALSE;
	fGapInformation=DetermineGap(fESDEvent, wCent, wFMD);
	if (fGapInformation==kBinNG) fHistEvent->Fill(kNG);
	else if (fGapInformation==kBinGA) fHistEvent->Fill(kGA);
	else if (fGapInformation==kBinGC) fHistEvent->Fill(kGC);
	else if (fGapInformation==kBinDG) {fHistEvent->Fill(kDG); fRunVsDG->Fill(fESDEvent->GetRunNumber());}
	else {
		//		printf("== This events doesn't have gap info! ==\n");
		PostOutputs();
		return;
	}
	//-------------------------------------------------------------------------

	// Martin's selection
	AliMultiplicitySelectionCP *selec = new AliMultiplicitySelectionCP();
	selec->InitDefaultTrackCuts(1);// 1 = b and c, 0 = d and e
	TArrayI indices;
	Int_t Nsel = selec->GetNumberOfITSTPCtracks(fESDEvent,indices);
	Double_t pionmass = 0.139570;
	TLorentzVector lv_track1;
	TLorentzVector lv_track2;
	TLorentzVector lv_sum;

	Int_t nmbTrk_NG = 0;
	Int_t nmbTrk_DG = 0;
	if (fGapInformation==kBinNG) {
		//2tracks combination with Standard ITSTPC cut
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
		//Only 2tracks with standard ITSTPC cut
		if (nmbTrk_NG==2) {
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
					fMassNG_st2t->Fill(lv_sum.M());
				}
			}
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
	else if (fGapInformation==kBinDG) {
		if (Nsel == 2) {
			fHistEvent->Fill(k2Tracks);
			fRunVs2t->Fill(fESDEvent->GetRunNumber());
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
					fMassDG_MS->Fill(lv_sum.M());
				}
			}
		}
		else if (Nsel == 4) {fHistEvent->Fill(k4Tracks); fRunVs4t->Fill(fESDEvent->GetRunNumber());}
		else if (Nsel == 6) fHistEvent->Fill(k6Tracks);
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
		if (nmbTrk_DG==2) {
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
					fMassDG_st2t->Fill(lv_sum.M());
				}
			}
		}
		fMultDG_MS->Fill(Nsel);
		fMultDG->Fill(nmbTrk_DG);
	}

	Bool_t fCheck2tracks = kFALSE;
	Bool_t fCheck4tracks = kFALSE;
	Bool_t fCheck6tracks = kFALSE;

	delete selec;

	// TODO::ITSSA track should be updated!!!
	Int_t tmp_Ntrk = fESDEvent->GetNumberOfTracks();
	Int_t Ntrk_ITSSA = 0;
	for (Int_t i = 0; i < tmp_Ntrk; i++) {
		AliESDtrack *track = fESDEvent->GetTrack(i);
		if(!track) continue;
		if(!fTrackCuts_ITSSA->AcceptTrack(track)) continue;
		Ntrk_ITSSA++;
	}
	Bool_t fCheck2tracks_ITSSA = (Ntrk_ITSSA==2) ? kTRUE : kFALSE;
	Bool_t fCheck4tracks_ITSSA = (Ntrk_ITSSA==4) ? kTRUE : kFALSE;
	if (fCheck2tracks_ITSSA && fGapInformation == kBinDG) fRunVs2t_ITSSA->Fill(fESDEvent->GetRunNumber());
	if (fCheck4tracks_ITSSA && fGapInformation == kBinDG) fRunVs4t_ITSSA->Fill(fESDEvent->GetRunNumber());

	// TAESOO PWA TREE INFO ---------------------------------------------------
	if (fGapInformation != kBinDG) {
//		printf("== This events are not double gap! ==\n");
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kCheckDG);
	fTrackCutsInfo->Fill(Nsel);

	//	if (!CheckV0Hit(fESDEvent,fHitV0A, fHitV0C)) {
	//		printf("== Should be no hit on V0!! ==\n");
	//		PostOutputs();
	//		return;
	//	}
	fHistEvent->Fill(kCheckV0Hit);

	fCheck2tracks = (Nsel == 2) ? kTRUE : kFALSE;
	fCheck4tracks = (Nsel == 4) ? kTRUE : kFALSE;
	fCheck6tracks = (Nsel == 6) ? kTRUE : kFALSE;

	if (fCheck2tracks) {
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
			//TPC 
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
			fTwoPionMask_TPC[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_TPC[i]);
			//TOF
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
			fTwoPionMask_TOF[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_TOF[i]);
			//ITS
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetITS);
			fTwoPionMask_ITS[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_ITS[i]);
			//TRD
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTRD);
			fTwoPionMask_TRD[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_TRD[i]);
			//TPC|TOF|ITS|TRD
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC |
					AliPIDResponse::kDetTOF |
					AliPIDResponse::kDetITS |
					AliPIDResponse::kDetTRD);
			fTwoPionMask_tot[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fTwoPionBayesProb_tot[i]);
			fTwoPionDetMask_tot[i] = (UInt_t)fPIDCombined->GetDetectorMask();
		}
		fCheckTwoPion = kTRUE;
	}

	if (fCheck4tracks) {
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
			//TPC 
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
			fFourPionMask_TPC[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fFourPionBayesProb_TPC[i]);
			//TOF
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
			fFourPionMask_TOF[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fFourPionBayesProb_TOF[i]);
			//ITS
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetITS);
			fFourPionMask_ITS[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fFourPionBayesProb_ITS[i]);
			//TRD
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTRD);
			fFourPionMask_TRD[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fFourPionBayesProb_TRD[i]);
			//TPC|TOF|ITS|TRD
			fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC |
					AliPIDResponse::kDetTOF |
					AliPIDResponse::kDetITS |
					AliPIDResponse::kDetTRD);
			fFourPionMask_tot[i] = fPIDCombined->ComputeProbabilities(track,fPIDResponse,fFourPionBayesProb_tot[i]);
			fFourPionDetMask_tot[i] = (UInt_t)fPIDCombined->GetDetectorMask();

		}
		fCheckFourPion = kTRUE;

	}

	//-------------------------------------------------------------------------

	// TODO :: ITSSA part
	if (fCheck2tracks_ITSSA) {
		fCheckTwoPion_ITSSA = kTRUE;
	}
	
	indices = 0x0;

	if (!fIsMC) {
		if (fCheck2tracks || fCheck4tracks) fTree->Fill();
	}
	PostOutputs();
	return;
}
//______________________________________________________________________________
void AliAnalysisTaskCDTree::PostOutputs()
{
	if (fIsMC) fTree->Fill();
	PostData(1,fTree);
	PostData(2,fList);
	return;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDTree::CutEvent(const AliESDEvent *ESDEvent, TH1 *hspd,
		TH1* hfochans, 
		TH1 *hpriVtxX, TH1 *hpriVtxY, TH1 *hpriVtxZ, Double_t *hVertex)
{
	AliTriggerAnalysis triggerAnalysis;

	// SPD FastOR check-------------------------------------------------------
	const Int_t fastORhw = triggerAnalysis.SPDFiredChips(ESDEvent, 0); //From hardware bits
	if (hspd) hspd->Fill(fastORhw);

	if (hfochans) {
		const AliMultiplicity *mult = ESDEvent->GetMultiplicity();

		for(Int_t iChipKey=0; iChipKey < 1200; iChipKey++){
			if (mult->TestFastOrFiredChips(iChipKey)) {
				hfochans->Fill((Double_t)iChipKey);
			}
		}
	}
	//-------------------------------------------------------------------------

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
	const Bool_t kpriv = kpr0 && (fabs(vertex->GetZ()) < 10.);
	if(!kpriv) return kFALSE;

	hVertex[0] = vertex->GetX();
	hVertex[1] = vertex->GetY();
	hVertex[2] = vertex->GetZ();

	if(hpriVtxX) hpriVtxX->Fill(vertex->GetX());
	if(hpriVtxY) hpriVtxY->Fill(vertex->GetY());
	if(hpriVtxZ) hpriVtxZ->Fill(vertex->GetZ());
	//-------------------------------------------------------------------------
	
	return kTRUE;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDTree::CheckInput()
{
	//General protection------------------------------------------------------
	if (const AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler)) {
		fESDEvent = (AliESDEvent*)esdH->GetEvent();
	}

	if(!fESDEvent){
//		printf("AliAnalysisTaskCDTest No valid event\n");
		return kFALSE;
	}

	fPIDResponse = (AliPIDResponse*)fInputHandler->GetPIDResponse();
	if(!fPIDResponse) {
//		printf("PIDResponse is not working!!\n");
		return kFALSE;
	}
	//-------------------------------------------------------------------------

	//Check Magnetic Field-----------------------------------------------------
	if(TMath::Abs(fESDEvent->GetMagneticField()) < 1) {
//		printf("AliAnalysisTaskCDTree strange Bfield! %f\n", fESDEvent->GetMagneticField());
		return kFALSE;
	}
	//-------------------------------------------------------------------------
	
	//Get MC event-------------------------------------------------------------
	fMCEvent = MCEvent();
//	if (fMCEvent)
//		printf("This is MC Event");
	//-------------------------------------------------------------------------

	return kTRUE;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDTree::SPDLoc2Glo(const Int_t id, const Double_t *loc,
		Double_t *glo)
{
	//
	//SPDLoc2Glo, do not touch
	//

	static TGeoHMatrix mat;
	Int_t vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(id);
	if (vid<0) {
//		printf("AliCDMesonUtils Did not find module with such ID %d\n",id);
		return kFALSE;
	}
	AliITSAlignMille2Module::SensVolMatrix(vid,&mat);
	mat.LocalToMaster(loc,glo);
	return kTRUE;
}
//------------------------------------------------------------------------------
Int_t AliAnalysisTaskCDTree::CheckChipEta(const Int_t chipKey, const TString scut,const Double_t vtxPos[],TH2 *hitMapSPDinner, TH2 *hitMapSPDouter)
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
Bool_t AliAnalysisTaskCDTree::CheckV0Hit(const AliESDEvent *ESDEvent, TH1D* hitV0A, TH1D* hitV0C)
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
Int_t AliAnalysisTaskCDTree::GetFastORmultiplicity(const AliESDEvent* ESDEvent)
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
Int_t AliAnalysisTaskCDTree::DetermineGap(const AliESDEvent *ESDEvent,const Bool_t wCent, const Bool_t wFMD)
{
	// determines the gap configuration for all gap tagging detectors based on the
	// data set which is available
	//

	//Read V0 first
	AliTriggerAnalysis triggerAnalysis;
	const Bool_t khw = kFALSE;
	const Bool_t v0A = triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kV0A);
	const Bool_t v0C = triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kV0C);
	//Check central activity
	Bool_t centAct = kFALSE;
	if (triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kSPDGFO)) {
		centAct = kTRUE;
	}
	//Check FMD activity by offine
	triggerAnalysis.SetFMDThreshold(0.3,0.5);
	const Bool_t FMDA = triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kFMDA);
	const Bool_t FMDC = triggerAnalysis.IsOfflineTriggerFired(ESDEvent,AliTriggerAnalysis::kFMDC);
	if (!wFMD) {
		if (wCent) {
			if (!v0A && !v0C && centAct) return kBinDG;
			else if (!v0A && v0C) return kBinGA;
			else if (v0A && !v0C) return kBinGC;
			else return kBinNG;
		}
		else {
			if (v0A && v0C) return kBinNG;
			else if (!v0A && !v0C) return kBinDG;
			else if (!v0A && v0C) return kBinGA;
			else return kBinGC;
		}
	}
	else {
		if (wCent) {
			if (!v0A && !v0C && !FMDA && !FMDC && centAct) return kBinDG;
			else if (!v0A && !FMDA && (v0C || FMDC)) return kBinGA;
			else if (!v0C && !FMDC && (v0A || FMDA)) return kBinGC;
			else return kBinNG;
		}
		else {
			if (!v0A && !v0C && !FMDA && !FMDC) return kBinDG;
			else if (!v0A && !FMDA && (v0C || FMDC)) return kBinGA;
			else if (!v0C && !FMDC && (v0A || FMDA)) return kBinGC;
			else return kBinNG;
		}
	}

}
//______________________________________________________________________________
void AliAnalysisTaskCDTree::DetermineMCprocessType()
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
				default: fMCprocessType = AliAnalysisTaskCDTree::kBinND; break;
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
Int_t AliAnalysisTaskCDTree::DoMCTruth()
{
	if (!fMCEvent) return -1;
	AliStack* stack = fMCEvent->Stack();
	if (!stack) return -1;

	// Multiplicity ----------------------------------------------------------
	// determine number of charged physical primaries on the stack
	Int_t nPhysicalPrimaries = 0;
	Int_t nProton = 0;
	Int_t nPrimaries = stack->GetNprimary();
	fPrimaries->Fill(nPrimaries);
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
void AliAnalysisTaskCDTree::Terminate(Option_t *)
{
}
//______________________________________________________________________________
