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
#include <TParticle.h>
#include <TParticlePDG.h>

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
	, fIsRun2(0)
	, fSavemode(0)
	, fCombmode(0)
	, fTree(0x0)
	, fTree_Comb(0x0)
	, fCheckTwoPion(0)
	, fCheckFourPion(0)
	, fCheckTwoPion_ITSSA(0)
	, fCheckFourPion_ITSSA(0)
	, fMultiplicity(0)
	, fSPDFired(0)
	, fV0Gap(0)
	, fADGap(0)
	, fFMDGap(0)
	, fZDCGap(0)
	, fIsMC(0)
	, fIsSaveGen(0)
	, fIsPythia8(0)
	, fIsPythia(0)
	, fIsPhojet(0)
	, fIsEPOS(0)
	, fRunNumber(0)
	, fPeriod(0)
	, fComb_IsPileUp(0)
	, fComb_IsPassClusterCut(0)
	, fComb_SPDCluster(0)
	, fComb_SPDTracklets(0)
	, fComb_MC_EventProcess(0)
	, fComb_IsPassMBOR(0)
	, fComb_IsPassVertex(0)
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
	, fRunVsDG_V0FMDSPD(0x0)
	, fRunVsDG_V0FMDZDCSPD(0x0)
	, fRunVs2t(0x0)
	, fRunVs2t_ITSSA(0x0)
	, fRunVs4t(0x0)
	, fRunVs4t_ITSSA(0x0)
	, fMultNG(0x0)
	, fMultNG_MS(0x0)
	, fMultDG(0x0)
	, fMultDG_MS(0x0)
	, fMassNG(0x0)
	, fMassNG_MS(0x0)
	, fMassDG(0x0)
	, fMassDG_MS(0x0)
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
	, fADANmbBB(0)
	, fADCNmbBB(0)
	, fTPCSignal(0x0)
	, fTOFSignal(0x0)
	, fITSSignal(0x0)
	, fTRDSignal(0x0)
	, fMult_Gen(0x0)
	, fMult_Gen_Process(0x0)
	, fMult_Rec_DG_Process(0x0)
	, fMult_Rec_NG_Process(0x0)
{
	//
	// standard constructor (the one which should be used)
	//
	// slot in TaskSE must start from 1

	for (Int_t i = 0; i < 7; i++) {
		fMC_Eta[i] = 0x0;
		fMC_DiffMass[i] = 0x0;
		fMC_DiffMass_PDG[i] = 0x0;
		if (i<6) fRunFiducial[i] = 0x0;
	}
	for (Int_t i = 0; i < 16; i++) {
		fADCharge[i] = 0.;
	}
	
	for (Int_t i = 0; i < 8; i++) {
		fRunOnline[i] = 0x0;
	}
	for (Int_t i = 0; i < 3; i++) {
		fVertex[i] = -999.;
	}
	for (Int_t i = 0; i < 2; i++) {//Two Pion
		fTwoPionMask_TPC[i] = 0;
		fTwoPionMask_TOF[i] = 0;
		fTwoPionMask_ITS[i] = 0;
		fTwoPionMask_TRD[i] = 0;
		fTwoPionMask_tot[i] = 0;
		fTwoPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 10; j++) {
			fTwoPionTrack[i][j] = 0;
			fTwoPionTPCSigma[i][j] = 0;
			fTwoPionTOFSigma[i][j] = 0;
			fTwoPionITSSigma[i][j] = 0;
			fTwoPionTrack_ITSSA[i][j] = 0;
			if ( j < AliPID::kSPECIES) {
				fTwoPionBayesProb_TPC[i][j] = 0;
				fTwoPionBayesProb_TOF[i][j] = 0;
				fTwoPionBayesProb_ITS[i][j] = 0;
				fTwoPionBayesProb_TRD[i][j] = 0;
				fTwoPionBayesProb_tot[i][j] = 0;
				fTwoPionTPCProb[i][j] = 0;
				fTwoPionTOFProb[i][j] = 0;
				fTwoPionITSProb[i][j] = 0;
			}
		}
	}
	for (Int_t i = 0; i < 4; i++) {//Four Pion
		fFourPionMask_TPC[i] = 0;
		fFourPionMask_TOF[i] = 0;
		fFourPionMask_ITS[i] = 0;
		fFourPionMask_TRD[i] = 0;
		fFourPionMask_tot[i] = 0;
		fFourPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 10; j++) {
			fFourPionTrack[i][j] = 0;
			fFourPionTPCSigma[i][j] = 0;
			fFourPionTOFSigma[i][j] = 0;
			fFourPionITSSigma[i][j] = 0;
			fFourPionTrack_ITSSA[i][j] = 0;
			if ( j < AliPID::kSPECIES) {
				fFourPionBayesProb_TPC[i][j] = 0;
				fFourPionBayesProb_TOF[i][j] = 0;
				fFourPionBayesProb_ITS[i][j] = 0;
				fFourPionBayesProb_TRD[i][j] = 0;
				fFourPionBayesProb_tot[i][j] = 0;
				fFourPionTPCProb[i][j] = 0;
				fFourPionTOFProb[i][j] = 0;
				fFourPionITSProb[i][j] = 0;
			}
		}
	}
	for (Int_t i = 0; i < 64; i++) {
		if (i < 2) {
			fComb_V0_Time_Mean[i] = -999.;
			fComb_AD_Time_Mean[i] = -999.;
			fComb_forwardP[i].SetPxPyPzE(0,0,0,0);
			fComb_diffSystem[i].SetPxPyPzE(0,0,0,0);
		}
		if (i < 10) {
			fComb_DetHit[i] = 0;
		}
		if (i < 16) {
			fComb_AD_Time[i] = -999.;
			fComb_AD_ADC[i] = -999.;
		}
		if (i < 64) {
			fComb_V0_Time[i] = -999.;
			fComb_V0_ADC[i] = -999.;
		}
	}


	// MC variable
	for (Int_t i=0; i< 5; i++) {
		for (Int_t j = 0; j < 4; j++) {
			fMCGenProtonTrack[j][i] = 0.;
			fMCGenPionTrack[j][i] = 0.;
		}
	}
	DefineOutput(1, TTree::Class());
	DefineOutput(2, TTree::Class());
	DefineOutput(3, TList::Class());
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
	, fIsRun2(0)
	, fSavemode(0)
	, fCombmode(0)
	, fTree(0x0)
	, fTree_Comb(0x0)
	, fCheckTwoPion(0)
	, fCheckFourPion(0)
	, fCheckTwoPion_ITSSA(0)
	, fCheckFourPion_ITSSA(0)
	, fMultiplicity(0)
	, fSPDFired(0)
	, fV0Gap(0)
	, fADGap(0)
	, fFMDGap(0)
	, fZDCGap(0)
	, fIsMC(0)
	, fIsSaveGen(0)
	, fIsPythia8(0)
	, fIsPythia(0)
	, fIsPhojet(0)
	, fIsEPOS(0)
	, fRunNumber(0)
	, fPeriod(0)
	, fComb_IsPileUp(0)
	, fComb_IsPassClusterCut(0)
	, fComb_SPDCluster(0)
	, fComb_SPDTracklets(0)
	, fComb_MC_EventProcess(0)
	, fComb_IsPassMBOR(0)
	, fComb_IsPassVertex(0)
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
	, fRunVsDG_V0FMDSPD(0x0)
	, fRunVsDG_V0FMDZDCSPD(0x0)
	, fRunVs2t(0x0)
	, fRunVs2t_ITSSA(0x0)
	, fRunVs4t(0x0)
	, fRunVs4t_ITSSA(0x0)
	, fMultNG(0x0)
	, fMultNG_MS(0x0)
	, fMultDG(0x0)
	, fMultDG_MS(0x0)
	, fMassNG(0x0)
	, fMassNG_MS(0x0)
	, fMassDG(0x0)
	, fMassDG_MS(0x0)
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
	, fADANmbBB(0)
	, fADCNmbBB(0)
	, fTPCSignal(0x0)
	, fTOFSignal(0x0)
	, fITSSignal(0x0)
	, fTRDSignal(0x0)
	, fMult_Gen(0x0)
	, fMult_Gen_Process(0x0)
	, fMult_Rec_DG_Process(0x0)
	, fMult_Rec_NG_Process(0x0)
{
	for (Int_t i = 0; i < 7; i++) {
		fMC_Eta[i] = 0x0;
		fMC_DiffMass[i] = 0x0;
		fMC_DiffMass_PDG[i] = 0x0;
		if (i<6) fRunFiducial[i] = 0x0;
	}
	for (Int_t i = 0; i < 16; i++) {
		fADCharge[i] = 0.;
	}
	for (Int_t i = 0; i < 8; i++) {
		fRunOnline[i] = 0x0;
	}
	for (Int_t i = 0; i < 3; i++) {
		fVertex[i] = -999.;
	}
	for (Int_t i = 0; i < 2; i++) {//Two Pion
		fTwoPionMask_TPC[i] = 0;
		fTwoPionMask_TOF[i] = 0;
		fTwoPionMask_ITS[i] = 0;
		fTwoPionMask_TRD[i] = 0;
		fTwoPionMask_tot[i] = 0;
		fTwoPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 10; j++) {
			fTwoPionTrack[i][j] = 0;
			fTwoPionTPCSigma[i][j] = 0;
			fTwoPionTOFSigma[i][j] = 0;
			fTwoPionITSSigma[i][j] = 0;
			fTwoPionTrack_ITSSA[i][j] = 0;
			if ( j < AliPID::kSPECIES) {
				fTwoPionBayesProb_TPC[i][j] = 0;
				fTwoPionBayesProb_TOF[i][j] = 0;
				fTwoPionBayesProb_ITS[i][j] = 0;
				fTwoPionBayesProb_TRD[i][j] = 0;
				fTwoPionBayesProb_tot[i][j] = 0;
				fTwoPionTPCProb[i][j] = 0;
				fTwoPionTOFProb[i][j] = 0;
				fTwoPionITSProb[i][j] = 0;
			}
		}
	}
	for (Int_t i = 0; i < 4; i++) {//Four Pion
		fFourPionMask_TPC[i] = 0;
		fFourPionMask_TOF[i] = 0;
		fFourPionMask_ITS[i] = 0;
		fFourPionMask_TRD[i] = 0;
		fFourPionMask_tot[i] = 0;
		fFourPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 10; j++) {
			fFourPionTrack[i][j] = 0;
			fFourPionTPCSigma[i][j] = 0;
			fFourPionTOFSigma[i][j] = 0;
			fFourPionITSSigma[i][j] = 0;
			fFourPionTrack_ITSSA[i][j] = 0;
			if ( j < AliPID::kSPECIES) {
				fFourPionBayesProb_TPC[i][j] = 0;
				fFourPionBayesProb_TOF[i][j] = 0;
				fFourPionBayesProb_ITS[i][j] = 0;
				fFourPionBayesProb_TRD[i][j] = 0;
				fFourPionBayesProb_tot[i][j] = 0;
				fFourPionTPCProb[i][j] = 0;
				fFourPionTOFProb[i][j] = 0;
				fFourPionITSProb[i][j] = 0;
			}
		}
	}
	for (Int_t i = 0; i < 64; i++) {
		if (i < 2) {
			fComb_V0_Time_Mean[i] = -999.;
			fComb_AD_Time_Mean[i] = -999.;
			fComb_forwardP[i].SetPxPyPzE(0,0,0,0);
			fComb_diffSystem[i].SetPxPyPzE(0,0,0,0);
		}
		if (i < 10) {
			fComb_DetHit[i] = 0;
		}
		if (i < 16) {
			fComb_AD_Time[i] = -999.;
			fComb_AD_ADC[i] = -999.;
		}
		if (i < 64) {
			fComb_V0_Time[i] = -999.;
			fComb_V0_ADC[i] = -999.;
		}
	}

	for (Int_t i=0; i< 5; i++) {
		for (Int_t j = 0; j < 4; j++) {
			fMCGenProtonTrack[j][i] = 0.;
			fMCGenPionTrack[j][i] = 0.;
		}
	}
}
//______________________________________________________________________________
AliAnalysisTaskCDPWA::~AliAnalysisTaskCDPWA()
{
	//Destructor(pointer should be deleted)
	if (fTree) delete fTree;
	if (fTree_Comb) delete fTree_Comb;
	if (fList) delete fList;
	if (fPhysicsSelection) delete fPhysicsSelection;
	if (fTrackCuts) delete fTrackCuts;
	if (fTrackCuts_ITSSA) delete fTrackCuts_ITSSA;
	if (fHistEvent) delete fHistEvent;
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
	for (Int_t i = 0; i < 8; i++) {
		if (fRunOnline[i]) delete fRunOnline[i];
	}
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
	if (fRunVsDG_V0FMDSPD) delete fRunVsDG_V0FMDSPD;
	if (fRunVsDG_V0FMDZDCSPD) delete fRunVsDG_V0FMDZDCSPD;
	if (fRunVs2t) delete fRunVs2t;
	if (fRunVs2t_ITSSA) delete fRunVs2t_ITSSA;
	if (fRunVs4t) delete fRunVs4t;
	if (fRunVs4t_ITSSA) delete fRunVs4t_ITSSA;
	if (fMultNG) delete fMultNG;
	if (fMultNG_MS) delete fMultNG_MS;
	if (fMultDG) delete fMultDG;
	if (fMultDG_MS) delete fMultDG_MS;
	if (fMassNG) delete fMassNG;
	if (fMassNG_MS) delete fMassNG_MS;
	if (fMassDG) delete fMassDG;
	if (fMassDG_MS) delete fMassDG_MS;
	if (fTrackCutsInfo) delete fTrackCutsInfo;
	if (fTrackCutsInfo_ITSSA) delete fTrackCutsInfo_ITSSA;
	if (fhClusterVsTracklets_bf) delete fhClusterVsTracklets_bf;
	if (fhClusterVsTracklets_af) delete fhClusterVsTracklets_af;
	if (fMult_ITSSA) delete fMult_ITSSA;
	if (fMult_DG_ITSSA) delete fMult_DG_ITSSA;
	if (fMult_NG_ITSSA) delete fMult_NG_ITSSA;
	if (fMult_ITSSA_MS) delete fMult_ITSSA_MS;
	if (fMult_DG_ITSSA_MS) delete fMult_DG_ITSSA_MS;
	if (fMult_NG_ITSSA_MS) delete fMult_NG_ITSSA_MS;
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
	if (fTPCSignal) delete fTPCSignal;
	if (fTOFSignal) delete fTOFSignal;
	if (fITSSignal) delete fITSSignal;
	if (fTRDSignal) delete fTRDSignal;
	if (fMult_Gen) delete fMult_Gen;
	if (fMult_Gen_Process) delete fMult_Gen_Process;
	if (fMult_Rec_DG_Process) delete fMult_Rec_DG_Process;
	if (fMult_Rec_NG_Process) delete fMult_Rec_NG_Process;
	for (Int_t i = 0; i < 7; i++) {
		if (fMC_Eta[i]) delete fMC_Eta[i];
		if (fMC_DiffMass[i]) delete fMC_DiffMass[i];
		if (fMC_DiffMass_PDG[i]) delete fMC_DiffMass_PDG[i];
		if (i<6 && fRunFiducial[i]) delete fRunFiducial[i];
	}
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
	fTree->Branch("CheckFourPion_ITSSA",&fCheckFourPion_ITSSA);
	fTree->Branch("SPDFired",&fSPDFired);
	fTree->Branch("V0Gap",&fV0Gap);
	fTree->Branch("ADGap",&fADGap);
	fTree->Branch("FMDGap",&fFMDGap);
	fTree->Branch("ZDCGap",&fZDCGap);
	fTree->Branch("RunNumber",&fRunNumber);
	fTree->Branch("Period",&fPeriod);
	fTree->Branch("Multiplicity",&fMultiplicity);
	//For vertex
	for (Int_t i = 0; i < 3; i++) {
		fTree->Branch(Form("Vertex_%d",i),&fVertex[i]);
	}
	if (fIsRun2) {
		fTree->Branch("ADANmbBB",&fADANmbBB);
		fTree->Branch("ADCNmbBB",&fADCNmbBB);
		for (Int_t i = 0; i < 16; i++) {
			fTree->Branch(Form("ADCharge_%d",i),&fADCharge[i]);
		}
	}
	for (Int_t i = 0; i < 2; i++) {
		fTree->Branch(Form("TwoPionMask_TPC_%d",i),&fTwoPionMask_TPC[i]);
		fTree->Branch(Form("TwoPionMask_TOF_%d",i),&fTwoPionMask_TOF[i]);
		fTree->Branch(Form("TwoPionMask_ITS_%d",i),&fTwoPionMask_ITS[i]);
		fTree->Branch(Form("TwoPionMask_TRD_%d",i),&fTwoPionMask_TRD[i]);
		fTree->Branch(Form("TwoPionMask_tot_%d",i),&fTwoPionMask_tot[i]);
		fTree->Branch(Form("TwoPionDetMask_tot_%d",i),&fTwoPionDetMask_tot[i]);
		for (Int_t j = 0; j < 10; j++) {
			fTree->Branch(Form("TwoPionTrack_%d_%d",i,j),&fTwoPionTrack[i][j]);
			fTree->Branch(Form("TwoPionTPCSigma_%d_%d",i,j),&fTwoPionTPCSigma[i][j]);
			fTree->Branch(Form("TwoPionTOFSigma_%d_%d",i,j),&fTwoPionTOFSigma[i][j]);
			fTree->Branch(Form("TwoPionITSSigma_%d_%d",i,j),&fTwoPionITSSigma[i][j]);
			fTree->Branch(Form("TwoPionTrack_ITSSA_%d_%d",i,j),&fTwoPionTrack_ITSSA[i][j]);
			if ( j < AliPID::kSPECIES) {
				fTree->Branch(Form("TwoPionBayesProb_TPC_%d_%d",i,j),&fTwoPionBayesProb_TPC[i][j]);
				fTree->Branch(Form("TwoPionBayesProb_TOF_%d_%d",i,j),&fTwoPionBayesProb_TOF[i][j]);
				fTree->Branch(Form("TwoPionBayesProb_ITS_%d_%d",i,j),&fTwoPionBayesProb_ITS[i][j]);
				fTree->Branch(Form("TwoPionBayesProb_TRD_%d_%d",i,j),&fTwoPionBayesProb_TRD[i][j]);
				fTree->Branch(Form("TwoPionBayesProb_tot_%d_%d",i,j),&fTwoPionBayesProb_tot[i][j]);
				fTree->Branch(Form("TwoPionTPCProb_%d_%d",i,j),&fTwoPionTPCProb[i][j]);
				fTree->Branch(Form("TwoPionTOFProb_%d_%d",i,j),&fTwoPionTOFProb[i][j]);
				fTree->Branch(Form("TwoPionITSProb_%d_%d",i,j),&fTwoPionITSProb[i][j]);
			}
		}
	}
	for (Int_t i = 0; i < 4; i++) {
		fTree->Branch(Form("FourPionMask_TPC_%d",i),&fFourPionMask_TPC[i]);
		fTree->Branch(Form("FourPionMask_TOF_%d",i),&fFourPionMask_TOF[i]);
		fTree->Branch(Form("FourPionMask_ITS_%d",i),&fFourPionMask_ITS[i]);
		fTree->Branch(Form("FourPionMask_TRD_%d",i),&fFourPionMask_TRD[i]);
		fTree->Branch(Form("FourPionMask_tot_%d",i),&fFourPionMask_tot[i]);
		fTree->Branch(Form("FourPionDetMask_tot_%d",i),&fFourPionDetMask_tot[i]);
		for (Int_t j = 0; j < 10; j++) {
			fTree->Branch(Form("FourPionTrack_%d_%d",i,j),&fFourPionTrack[i][j]);
			fTree->Branch(Form("FourPionTPCSigma_%d_%d",i,j),&fFourPionTPCSigma[i][j]);
			fTree->Branch(Form("FourPionTOFSigma_%d_%d",i,j),&fFourPionTOFSigma[i][j]);
			fTree->Branch(Form("FourPionITSSigma_%d_%d",i,j),&fFourPionITSSigma[i][j]);
			fTree->Branch(Form("FourPionTrack_ITSSA_%d_%d",i,j),&fFourPionTrack_ITSSA[i][j]);
			if ( j < AliPID::kSPECIES) {
				fTree->Branch(Form("FourPionBayesProb_TPC_%d_%d",i,j),&fFourPionBayesProb_TPC[i][j]);
				fTree->Branch(Form("FourPionBayesProb_TOF_%d_%d",i,j),&fFourPionBayesProb_TOF[i][j]);
				fTree->Branch(Form("FourPionBayesProb_ITS_%d_%d",i,j),&fFourPionBayesProb_ITS[i][j]);
				fTree->Branch(Form("FourPionBayesProb_TRD_%d_%d",i,j),&fFourPionBayesProb_TRD[i][j]);
				fTree->Branch(Form("FourPionBayesProb_tot_%d_%d",i,j),&fFourPionBayesProb_tot[i][j]);
				fTree->Branch(Form("FourPionTPCProb_%d_%d",i,j),&fFourPionTPCProb[i][j]);
				fTree->Branch(Form("FourPionTOFProb_%d_%d",i,j),&fFourPionTOFProb[i][j]);
				fTree->Branch(Form("FourPionITSProb_%d_%d",i,j),&fFourPionITSProb[i][j]);
			}
		}
	}
	//For MC
	if (fIsMC) {//Simple parameter
	}
	if (fIsMC && fIsSaveGen) {//For the PWA
		for (Int_t i = 0; i < 5; i++) {
			for (Int_t j = 0; j < 4; j++) {
				fTree->Branch(Form("MCGenProtonTrack_%d_%d",j,i),&fMCGenProtonTrack[j][i]);
				fTree->Branch(Form("MCGenPionTrack_%d_%d",j,i),&fMCGenPionTrack[j][i]);
			}
		}
	}
	//-------------------------------------------------------------------------

	// TTree for combinatorics------------------------------------------------
	OpenFile(2);
	fTree_Comb = new TTree("tree_comb","tree_comb");
	fTree_Comb->Branch("Comb_RunNumber",&fRunNumber);
	fTree_Comb->Branch("Comb_Period",&fPeriod);
	fTree_Comb->Branch("Comb_IsPassMBOR",&fComb_IsPassMBOR);
	fTree_Comb->Branch("Comb_IsPassVertex",&fComb_IsPassVertex);
	fTree_Comb->Branch("Comb_IsPileUp",&fComb_IsPileUp);
	fTree_Comb->Branch("Comb_IsPassClusterCut",&fComb_IsPassClusterCut);
	fTree_Comb->Branch("Comb_SPDCluster",&fComb_SPDCluster);
	fTree_Comb->Branch("Comb_SPDTracklets",&fComb_SPDTracklets);
	for (Int_t i = 0; i < 64; i++) {
		if (i < 2) {
			fTree_Comb->Branch(Form("Comb_V0_Time_Mean_%d",i),&fComb_V0_Time_Mean[i]);
			fTree_Comb->Branch(Form("Comb_AD_TIme_Mean_%d",i),&fComb_AD_Time_Mean[i]);
		}
		if (i < 10) fTree_Comb->Branch(Form("Comb_DetHit_%d",i),&fComb_DetHit[i]);
		if (i < 16) {
			fTree_Comb->Branch(Form("Comb_AD_Time_%d",i),&fComb_AD_Time[i]);
			fTree_Comb->Branch(Form("Comb_AD_ADC_%d",i),&fComb_AD_ADC[i]);
		}
		if (i < 64) {
			fTree_Comb->Branch(Form("Comb_V0_Time_%d",i),&fComb_V0_Time[i]);
			fTree_Comb->Branch(Form("Comb_V0_ADC_%d",i),&fComb_V0_ADC[i]);
		}
		if (i < 3) {
			fTree_Comb->Branch(Form("Comb_Vertex_%d",i),&fVertex[i]);
		}
	}
	if (fIsMC) {
		fTree_Comb->Branch("Comb_MC_EventProcess",&fComb_MC_EventProcess);
		for (Int_t i = 0; i < 2; i++) {
			fTree_Comb->Branch(Form("Comb_forwardP_%d",i),&fComb_forwardP[i]);
			fTree_Comb->Branch(Form("Comb_diffSystem_%d",i),&fComb_diffSystem[i]);
		}
	}
	//-------------------------------------------------------------------------

	// TList-------------------------------------------------------------------
	// For the run by run analysis
	Int_t minRn = 114500;//min run number in Run1
	Int_t maxRn = 179300;//max run number in Run1
	if (fIsRun2) {//Run2
		minRn = 224800;
		maxRn = 260000;
	}
	Int_t diff = maxRn - minRn;

	OpenFile(3);
	fList = new TList();
	fList->SetOwner();
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
		fHistEvent->GetXaxis()->SetBinLabel(kNG_Data+1,"No-gap");
		fHistEvent->GetXaxis()->SetBinLabel(kSGA_Data+1,"Single-gap A");
		fHistEvent->GetXaxis()->SetBinLabel(kSGC_Data+1,"Single-gap C");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0+1,"DG_V0");
		fHistEvent->GetXaxis()->SetBinLabel(kDGAD+1,"DG_AD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0SPD+1,"DG_V0SPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGADSPD+1,"DG_ADSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADSPD+1,"DG_V0ADSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADFMDSPD+1,"DG_V0ADFMDSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0ADFMDZDCSPD+1,"DG_V0ADFMDZDCSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0FMDSPD+1,"DG_V0FMDSPD");
		fHistEvent->GetXaxis()->SetBinLabel(kDGV0FMDZDCSPD+1,"DG_V0FMDZDCSPD");
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
	if (!fRunVsDG_V0FMDSPD) {
		fRunVsDG_V0FMDSPD = new TH1D("fRunVsDG_V0FMDSPD","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_V0FMDSPD);
	}
	if (!fRunVsDG_V0FMDZDCSPD) {
		fRunVsDG_V0FMDZDCSPD = new TH1D("fRunVsDG_V0FMDZDCSPD","",diff,minRn,maxRn);
		fList->Add(fRunVsDG_V0FMDZDCSPD);
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
	if (!fMassDG_MS) {
		fMassDG_MS = new TH1D("fMassDG_MS","",50,0,2);
		fList->Add(fMassDG_MS);
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

	// Track Cuts
	if (!fIsRun2) {
		fTrackCuts = new AliESDtrackCuts();
		fTrackCuts->GetStandardITSTPCTrackCuts2010(1,0);//0 for d and e
		fTrackCuts_ITSSA = new AliESDtrackCuts();
		fTrackCuts_ITSSA->GetStandardITSSATrackCuts2010(1,0);//0 for no-PID mode
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
		fTrackCuts_ITSSA = new AliESDtrackCuts();
		{
			fTrackCuts_ITSSA->GetStandardITSSATrackCuts2010(1,0);//0 for noPID
			fTrackCuts_ITSSA->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
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
	// Initialize all tree varialbles------------------------------------------
	for (Int_t i = 0; i < 2; i++) {//Two Pion
		fTwoPionMask_TPC[i] = 0;
		fTwoPionMask_TOF[i] = 0;
		fTwoPionMask_ITS[i] = 0;
		fTwoPionMask_TRD[i] = 0;
		fTwoPionMask_tot[i] = 0;
		fTwoPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 10; j++) {
			fTwoPionTrack[i][j] = 0;
			fTwoPionTPCSigma[i][j] = 0;
			fTwoPionTOFSigma[i][j] = 0;
			fTwoPionITSSigma[i][j] = 0;
			fTwoPionTrack_ITSSA[i][j] = 0;
			if ( j < AliPID::kSPECIES) {
				fTwoPionBayesProb_TPC[i][j] = 0;
				fTwoPionBayesProb_TOF[i][j] = 0;
				fTwoPionBayesProb_ITS[i][j] = 0;
				fTwoPionBayesProb_TRD[i][j] = 0;
				fTwoPionBayesProb_tot[i][j] = 0;
				fTwoPionTPCProb[i][j] = 0;
				fTwoPionTOFProb[i][j] = 0;
				fTwoPionITSProb[i][j] = 0;
			}
		}
	}
	for (Int_t i = 0; i < 4; i++) {//Four Pion
		fFourPionMask_TPC[i] = 0;
		fFourPionMask_TOF[i] = 0;
		fFourPionMask_ITS[i] = 0;
		fFourPionMask_TRD[i] = 0;
		fFourPionMask_tot[i] = 0;
		fFourPionDetMask_tot[i] = 0;
		for (Int_t j = 0; j < 10; j++) {
			fFourPionTrack[i][j] = 0;
			fFourPionTPCSigma[i][j] = 0;
			fFourPionTOFSigma[i][j] = 0;
			fFourPionITSSigma[i][j] = 0;
			fFourPionTrack_ITSSA[i][j] = 0;
			if ( j < AliPID::kSPECIES) {
				fFourPionBayesProb_TPC[i][j] = 0;
				fFourPionBayesProb_TOF[i][j] = 0;
				fFourPionBayesProb_ITS[i][j] = 0;
				fFourPionBayesProb_TRD[i][j] = 0;
				fFourPionBayesProb_tot[i][j] = 0;
				fFourPionTPCProb[i][j] = 0;
				fFourPionTOFProb[i][j] = 0;
				fFourPionITSProb[i][j] = 0;
			}
		}
	}
	for (Int_t i = 0; i < 5; i++) {
		for (Int_t j = 0; j < 4; j++) {
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
	fCheckFourPion_ITSSA = kFALSE;

	fSPDFired = kFALSE;
	fV0Gap = kFALSE;
	fADGap = kFALSE;
	fFMDGap = kFALSE;
	fZDCGap = kFALSE;
	fRunNumber = -999;
	fPeriod = -999;
	fADANmbBB = 0;
	fADCNmbBB = 0;
	for (Int_t i = 0; i < 16; i++) {
		fADCharge[i] = 0.;
	}
	fMultiplicity = -999;
	fComb_IsPileUp = kFALSE;
	fComb_IsPassClusterCut = kFALSE;
	fComb_SPDCluster = -999;
	fComb_SPDTracklets = -999;
	fComb_IsPassMBOR = kFALSE;
	fComb_IsPassVertex = kFALSE;
	for (Int_t i = 0; i < 64; i++) {
		if (i < 2) {
			fComb_V0_Time_Mean[i] = -999.;
			fComb_AD_Time_Mean[i] = -999.;
		}
		if (i < 10) {
			fComb_DetHit[i] = kFALSE;
		}
		if (i < 16) {
			fComb_AD_Time[i] = -999.;
			fComb_AD_ADC[i] = -999.;
		}
		if (i < 64) {
			fComb_V0_Time[i] = -999.;
			fComb_V0_ADC[i] = -999.;
		}
	}
	if (fIsMC) {
		fComb_MC_EventProcess = -999;
		for (Int_t i = 0; i < 2; i++) {
			fComb_forwardP[i].SetPxPyPzE(0,0,0,0);
			fComb_diffSystem[i].SetPxPyPzE(0,0,0,0);
		}
	}
	//-------------------------------------------------------------------------

	// TODO:: Recommendtation by Christoph Mayer-------------------------------
	/*
	if (fESDEvent->IsIncompleteDAQ()) {
		PostOutputs();
		return;
	}
	*/
	//-------------------------------------------------------------------------

	//Check input file is corrupted or not-------------------------------------
	//Load Input handler and MC
	if (!CheckInput()) {
		PostOutputs();
		return;
	}
	fRunNumber = (Int_t)fESDEvent->GetRunNumber();
	fHistEvent->Fill(kInput);//Check total number
	fRunFiducial[0]->Fill(fRunNumber);//Check run by run

	if (!fIsRun2) {
		if (fRunNumber < 117224) fPeriod = 1;//10b
		else if (fRunNumber < 121041) fPeriod = 2;//10c
		else if (fRunNumber < 126438) fPeriod = 3;//10d
		else if (fRunNumber < 130851) fPeriod = 4;//10e
		else if (fRunNumber < 135032) fPeriod = 5;//10f
		else fPeriod = 6;//12b (8 TeV)
	}
	else {
		if (fRunNumber >= 224891 && fRunNumber < 227750) fPeriod = 1;//15f
		else if (fRunNumber < 232465) fPeriod = 2;//15g
		else if (fRunNumber < 235196) fPeriod = 3;//15h
		else if (fRunNumber < 236892) fPeriod = 4;//15i
		else if (fRunNumber < 238682) fPeriod = 5;//15j
		else if (fRunNumber < 239207) fPeriod = 6;//15k
		else if (fRunNumber < 243374) fPeriod = 7;//15l
		else if (fRunNumber < 244340) fPeriod = 8;//15m
		else if (fRunNumber < 244824) fPeriod = 9;//15n
		else if (fRunNumber <= 246994) fPeriod = 10;//15o
	}
	//-------------------------------------------------------------------------

	//MC TRUTH-----------------------------------------------------------------
	Int_t nMCprimaries = 0;
	DetermineMCprocessType();//Get MC process type in here
	if (fIsMC) fComb_MC_EventProcess = fMCprocessType;
	//Fill Histo for MC processes
	for (Int_t i = 0; i < kBinMCAll; i++) {
		if (fMCprocessType == i) 
			fHistEventProcesses->Fill(i);
	}
	if (fMCEvent) {
		nMCprimaries = DoMCTruth();
		if (fIsSaveGen && nMCprimaries == -1) { // Something is wrong with DIME, DRgen
			PostOutputs();
			return;
		}
	}
	fHistEvent->Fill(kMCCheck);
	//-------------------------------------------------------------------------

	//At least we have to use CINT10 + C0SMB for Run2--------------------------
	if (fIsRun2 && !fIsMC && !CheckOnlineTrigger(fESDEvent)) {
		PostOutputs();
		return;
	}
	fHistEvent->Fill(kOnlineTrigger);
	fRunFiducial[1]->Fill(fRunNumber);
	//-------------------------------------------------------------------------

	//Offline MBOR cut(available only for pass2)-------------------------------
	AliESDAD *ad = NULL;
	if (fIsRun2) ad = (AliESDAD*)fESDEvent->GetADData();//AD is included in the Run2
	AliESDVZERO *vzero = (AliESDVZERO*)fESDEvent->GetVZEROData();
	if (!fIsRun2 && !vzero) {
		PostOutputs();
		return;
	}
	if(fIsRun2 && (!ad || !vzero)) {
		PostOutputs();
		return;
	}
	//Timing check before applying cut
	DoTimeV0AD(fESDEvent, fV0ATime_bf, fV0CTime_bf, fADATime_bf, fADCTime_bf);

	//Apply offline MB_OR cut
	AliInputEventHandler *inputHandler = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
	if (!fIsRun2) {//For Run1
		if (inputHandler->IsEventSelected()) {//passed
			fComb_IsPassMBOR = kTRUE;
			fHistEvent->Fill(kOfflineCut);
			fRunFiducial[2]->Fill(fRunNumber);
		}
	}
	else {//For Run2
		if (inputHandler->IsEventSelected() & AliVEvent::kUserDefined) {
			fComb_IsPassMBOR = kTRUE;
			fHistEvent->Fill(kOfflineCut);
			fRunFiducial[2]->Fill(fRunNumber);
		}
	}
	//Time measurement after cutting out
	if (fComb_IsPassMBOR) {
		DoTimeV0AD(fESDEvent, fV0ATime_af, fV0CTime_af, fADATime_af, fADCTime_af);
		// SPD FastOR chip check to see noisy chips
		DoSPDCheck(fESDEvent, fSPDFiredChipsCluster, fSPDFiredChipsHardware, fSPDwc);
	}
	//-------------------------------------------------------------------------

	// Event selection :: Vertex cut-------------------------------------------
	const Bool_t IsVertexCut = CutEvent(fESDEvent, fHistPrimVtxX, fHistPrimVtxY, fHistPrimVtxZ);
	if (IsVertexCut) {
		if (fComb_IsPassMBOR) {
			fHistEvent->Fill(kVtxCut); //After MBOR/Vertex cut
			fRunFiducial[3]->Fill(fRunNumber);
		}
		fComb_IsPassVertex = kTRUE;
	}
	//-------------------------------------------------------------------------

	// Event selection :: Pile-up cut -----------------------------------------
	Double_t pileUpValue[5] = {2., 0.8, 3., 2., 5.};
	if (fIsRun2) pileUpValue[0] = 5.;//From dN/dpT at 13 TeV analysis
	const Bool_t isPileup = fESDEvent->IsPileupFromSPD(pileUpValue[0], pileUpValue[1], pileUpValue[2], pileUpValue[3], pileUpValue[4]);
	fComb_IsPileUp = isPileup;
	if (!fComb_IsPileUp && fComb_IsPassMBOR && fComb_IsPassVertex) {
		fHistEvent->Fill(kPileUpCut); // After MBOR/Vertex/pile up
		fRunFiducial[4]->Fill(fRunNumber);
	}
	//-------------------------------------------------------------------------

	// Cluster Cut ------------------------------------------------------------
	// Note that this is responsible for USER!!
	const AliMultiplicity *tmp_mult = fESDEvent->GetMultiplicity();
	Int_t nClustersLayer0 = fESDEvent->GetNumberOfITSClusters(0);
	Int_t nClustersLayer1 = fESDEvent->GetNumberOfITSClusters(1);
	Int_t nTracklets      = tmp_mult->GetNumberOfTracklets();
	fComb_SPDCluster = nClustersLayer0 + nClustersLayer1;
	fComb_SPDTracklets = nTracklets;

	if (!fComb_IsPileUp && fComb_IsPassMBOR && fComb_IsPassVertex) fhClusterVsTracklets_bf->Fill(nTracklets,nClustersLayer0+nClustersLayer1);
	Double_t cut_slope = 4.;//4 is default
	Double_t cut_b = 65.;//65 is default
	if (nClustersLayer0 + nClustersLayer1 <= (cut_b+cut_slope*nTracklets)) {
		fComb_IsPassClusterCut = kTRUE;
		if (!fComb_IsPileUp && fComb_IsPassMBOR && fComb_IsPassVertex) {
			fhClusterVsTracklets_af->Fill(nTracklets,nClustersLayer0+nClustersLayer1);
			fHistEvent->Fill(kClusterCut); //Cluster Cut
			fRunFiducial[5]->Fill(fRunNumber);
		}
	}

	//Combinatorics study------------------------------------------------------
	Bool_t IsSPD = kFALSE;
	Bool_t IsV0A = kFALSE;
	Bool_t IsV0C = kFALSE;
	Bool_t IsADA = kFALSE;
	Bool_t IsADC = kFALSE;
	Bool_t IsFMDA = kFALSE;
	Bool_t IsFMDC = kFALSE;
	Bool_t IsZDCA = kFALSE;
	Bool_t IsZDCC = kFALSE;

	if (fCombmode) DoCombStudy(ad,vzero);
	AliTriggerAnalysis fTrigger;
	fTrigger.SetFMDThreshold(0.3,0.5);//FMD Threshold at 7/13TeV is used.

	if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kSPDGFO)) {
		fComb_DetHit[0] = kTRUE;
		fSPDFired = kTRUE;
		IsSPD = kTRUE;
	}
	if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0A)) {
		fComb_DetHit[1] = kTRUE;
		IsV0A = kTRUE;
	}
	if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kV0C)) {
		fComb_DetHit[2] = kTRUE;
		IsV0C = kTRUE;
	}
	if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kFMDA)) {
		fComb_DetHit[3] = kTRUE;
		IsFMDA = kTRUE;
	}
	if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kFMDC)) {
		fComb_DetHit[4] = kTRUE;
		IsFMDC = kTRUE;
	}
	if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kZDCA)) {
		fComb_DetHit[5] = kTRUE;
		IsZDCA = kTRUE;
	}
	if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kZDCC)) {
		fComb_DetHit[6] = kTRUE;
		IsZDCC = kTRUE;
	}
	if (fIsRun2) {
		if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kADA)) {
			fComb_DetHit[7] = kTRUE;
			IsADA = kTRUE;
		}
		if (fTrigger.IsOfflineTriggerFired(fESDEvent,AliTriggerAnalysis::kADC)) {
			fComb_DetHit[8] = kTRUE;
			IsADC = kTRUE;
		}
	}
	//-------------------------------------------------------------------------

	//Fill combinatorics here!!
	if (fCombmode) fTree_Comb->Fill();
	if (!fComb_IsPassMBOR || !fComb_IsPassVertex || fComb_IsPileUp || !fComb_IsPassClusterCut) {//Return to the normal CEP analysis to reduce memory-problem
		PostOutputs();
		return;
	}
	//-------------------------------------------------------------------------
	
	//CEP ANLYSIS START FROM HERE!!

	// QA ANANLYSIS for PID----------------------------------------------------
	Double_t spdc = TMath::C()*1.E-9;// m/ns
	Double_t fTOFlength = 0.;
	Double_t timeTOF = 0.;
	Double_t fTOFTime = 0.;
	Double_t fTOFPIDmom = 0.;
	Double_t fBeta = 0.;
	for (Int_t iTrk = 0; iTrk < fESDEvent->GetNumberOfTracks(); iTrk++) {
		AliESDtrack *trk = (AliESDtrack*)fESDEvent->GetTrack(iTrk);
		if(!trk) continue;
		if(!fTrackCuts->AcceptTrack(trk)) continue;
		fTPCSignal->Fill(trk->P(),trk->GetTPCsignal());
		fITSSignal->Fill(trk->P(),trk->GetITSsignal());
		fTRDSignal->Fill(trk->P(),trk->GetTRDsignal());
		//// *** TOF *** //// for beta extraction
		fTOFlength = trk->GetIntegratedLength();
		if(fTOFlength<=0) continue;
		fTOFPIDmom = trk->P();
		if((fTOFPIDmom*fTOFPIDmom)==0) continue;
		fTOFlength = fTOFlength*0.01;
		timeTOF = trk->GetTOFsignal();
		fTOFTime = timeTOF*1E-3;
		fTOFTime = fTOFTime*spdc;
		fBeta = fTOFlength/fTOFTime;
		fTOFSignal->Fill(trk->P(),fBeta);
	}
	//-------------------------------------------------------------------------

	// TRIGGER ANALYSIS -------------------------------------------------------
	// MBOR and MBAND for Run2
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
	Bool_t fDGV0SPD = (!IsV0A && !IsV0C && IsSPD) ? kTRUE : kFALSE;
	Bool_t fDGADSPD = (!IsADA && !IsADC && IsSPD) ? kTRUE : kFALSE;//Always true for 7TeV

	// Sub-sample for double gap conditions
	// will be stroed in the TTree
	fV0Gap = (!IsV0A && !IsV0C) ? kTRUE : kFALSE;
	fADGap = (!IsADA && !IsADC) ? kTRUE : kFALSE;// Always true for 7TeV
	fFMDGap = (!IsFMDA && !IsFMDC) ? kTRUE : kFALSE;
	fZDCGap = (!IsZDCA && !IsZDCC) ? kTRUE : kFALSE;
	if (fV0Gap) { 
		fHistEvent->Fill(kDGV0);//Total number of events for DG_V0
		fRunVsDG_V0->Fill(fESDEvent->GetRunNumber());//Run-by-Run
	}
	Int_t nmbBB_Online_C = 0;
	Int_t nmbBB_Online_A = 0;
	if (fADGap) {//For the AD gap
		if (fIsRun2) {
			fHistEvent->Fill(kDGAD); 
			fRunVsDG_AD->Fill(fESDEvent->GetRunNumber());
			for (Int_t i = 0; i < 8; i++) {
				if(ad->GetBBFlag(i)) nmbBB_Online_C++;//C-side
				if(ad->GetBBFlag(i+8)) nmbBB_Online_A++;//A-side
			}
			fADANmbBB = nmbBB_Online_A;
			fADCNmbBB = nmbBB_Online_C;
			for (Int_t i = 0; i < 16; i++) {
				fADCharge[i] = ad->GetAdc(i);
			}
		}
	}
	if (fDGV0SPD) {//For V0+SPD gap
		fHistEvent->Fill(kDGV0SPD); 
		fRunVsDG_V0SPD->Fill(fESDEvent->GetRunNumber());
	}
	if (fDGADSPD) {//For AD+SPD gap
		fHistEvent->Fill(kDGADSPD); 
		fRunVsDG_ADSPD->Fill(fESDEvent->GetRunNumber());
	}

	// Local variable for tmp..
	Bool_t fDGV0ADSPD = (fDGV0SPD && fADGap) ? kTRUE : kFALSE;

	// Sub-sample histogram
	if (fDGV0SPD && fADGap) {fHistEvent->Fill(kDGV0ADSPD); fRunVsDG_V0ADSPD->Fill(fESDEvent->GetRunNumber());}
	if (fDGV0SPD && fADGap && fFMDGap) {fHistEvent->Fill(kDGV0ADFMDSPD); fRunVsDG_V0ADFMDSPD->Fill(fESDEvent->GetRunNumber());}
	if (fDGV0SPD && fADGap && fFMDGap && fZDCGap) {fHistEvent->Fill(kDGV0ADFMDZDCSPD); fRunVsDG_V0ADFMDZDCSPD->Fill(fESDEvent->GetRunNumber());}
	if (fDGV0SPD && fFMDGap) {fHistEvent->Fill(kDGV0FMDSPD); fRunVsDG_V0FMDSPD->Fill(fESDEvent->GetRunNumber());}
	if (fDGV0SPD && fFMDGap && fZDCGap) {fHistEvent->Fill(kDGV0FMDZDCSPD); fRunVsDG_V0FMDZDCSPD->Fill(fESDEvent->GetRunNumber());}

	// Define no-gap events to compare(only using V0)
	Bool_t IsNG = kFALSE;
	Bool_t IsGapA = kFALSE;
	Bool_t IsGapC = kFALSE;

	if (IsV0A && !IsV0C) {
		IsGapC = kTRUE;
		fHistEvent->Fill(kSGC_Data);
	}
	if (IsV0C && !IsV0A) {
		IsGapA = kTRUE;
		fHistEvent->Fill(kSGA_Data);
	}
	if (!IsGapC && !IsGapA && !fDGV0SPD) {
		IsNG = kTRUE;
		fHistEvent->Fill(kNG_Data);
	}
	//-------------------------------------------------------------------------

	// Track cuts by Martin's selection----------------------------------------
	Bool_t is10bc = kFALSE;
	if (fPeriod <= 2) is10bc=kTRUE;
	AliMultiplicitySelectionCPPWA *selec = new AliMultiplicitySelectionCPPWA();
	//(clusterCut,useITSSA,isRun2)
	if (!fIsRun2) selec->InitDefaultTrackCuts((Int_t)is10bc,0,kFALSE);// 1 = b,c and 0 = d,e
	else selec->InitDefaultTrackCuts(0,0,kTRUE);
	TArrayI indices;
	Int_t Nsel = selec->GetNumberOfITSTPCtracks(fESDEvent,indices);
	fMultiplicity = Nsel;
	Double_t pionmass = 0.139570;
	TLorentzVector lv_track1;
	TLorentzVector lv_track2;
	TLorentzVector lv_sum;
	//Multiplicity distribution with generated info
	if (fIsMC) {
		if (fDGV0SPD) fMult_Rec_DG_Process->Fill(fMCprocessType,Nsel);
		if (IsNG) fMult_Rec_NG_Process->Fill(fMCprocessType,Nsel);
	}
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
	else if (fDGV0SPD) {//For DG only for V0 with standard trackcuts
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
		else if (Nsel == 4) {
			fHistEvent->Fill(k4Tracks);
			fRunVs4t->Fill(fESDEvent->GetRunNumber());
		}
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
	Bool_t fCheck4tracks_ITSSA = kFALSE;

	delete selec;
	//-------------------------------------------------------------------------

	//ITSSA trackcuts + Martin's selection-------------------------------------
	Int_t tmp_Ntrk = fESDEvent->GetNumberOfTracks();
	Int_t nmbTrack_ITSSA = 0;
	for (Int_t i = 0; i < tmp_Ntrk; i++) {
		AliESDtrack *track = fESDEvent->GetTrack(i);
		if(!track) continue;
		if(!fTrackCuts_ITSSA->AcceptTrack(track)) continue;
		//ITSSA Accepted track
		nmbTrack_ITSSA++;
	}
	//Multiplicity of ITSSA track
	fMult_ITSSA->Fill(nmbTrack_ITSSA);
	if (fDGV0SPD) fMult_DG_ITSSA->Fill(nmbTrack_ITSSA);
	if (IsNG) fMult_NG_ITSSA->Fill(nmbTrack_ITSSA);

	//ITSSA + Martin's trackcuts
	AliMultiplicitySelectionCPPWA *selec_ITSSA = new AliMultiplicitySelectionCPPWA();
	if (!fIsRun2) selec_ITSSA->InitDefaultTrackCuts((Int_t)is10bc,1,kFALSE);// 1 = 10b and 10c, 0 = 10d and 10e
	else selec_ITSSA->InitDefaultTrackCuts(0,1,kTRUE);
	TArrayI indices_ITSSA;
	Int_t Nsel_ITSSA = selec_ITSSA->GetNumberOfITSSAtracks(fESDEvent,indices_ITSSA);
	fCheck2tracks_ITSSA = (Nsel_ITSSA == 2) ? kTRUE : kFALSE;
	fCheck4tracks_ITSSA = (Nsel_ITSSA == 4) ? kTRUE : kFALSE;
	//Multiplicity of ITSSA+Martin's trackcut
	fMult_ITSSA_MS->Fill(Nsel_ITSSA);
	if (fDGV0SPD) fMult_DG_ITSSA_MS->Fill(Nsel_ITSSA);
	if (IsNG) fMult_NG_ITSSA_MS->Fill(Nsel_ITSSA);
	if (fCheck2tracks_ITSSA && fDGV0SPD) {
		fRunVs2t_ITSSA->Fill(fESDEvent->GetRunNumber());
		fHistEvent->Fill(k2Tracks_ITSSA);
	}
	if (fCheck4tracks_ITSSA && fDGV0SPD) {
		fRunVs4t_ITSSA->Fill(fESDEvent->GetRunNumber());
		fHistEvent->Fill(k4Tracks_ITSSA);
	}
	delete selec_ITSSA;
	//-------------------------------------------------------------------------

	// TAESOO PWA TREE INFO ---------------------------------------------------
	if (!fIsRun2 && !fDGV0SPD) {//Run1
		PostOutputs();
		return;
	}
	if (fIsRun2) {
		if ((fV0Gap || fADGap)) {
			//printf("== This events are double gap ==\n");
		}
		else {
			//printf("== This events are not double gap! ==\n");
			PostOutputs();
			return;
		}
	}

	if(fDGV0SPD) {
		fTrackCutsInfo->Fill(Nsel);
		fTrackCutsInfo_ITSSA->Fill(Nsel_ITSSA);
	}
	if(fDGV0SPD) {
		if (fIsRun2) {
			fADATime_V0SPD->Fill(ad->GetADATime());
			fADCTime_V0SPD->Fill(ad->GetADCTime());
		}
	}
	if(fDGV0ADSPD) {
		if (fIsRun2) {
			fADATime_V0ADSPD->Fill(ad->GetADATime());
			fADCTime_V0ADSPD->Fill(ad->GetADCTime());
		}
	}

	fCheck2tracks = (Nsel == 2) ? kTRUE : kFALSE;
	fCheck4tracks = (Nsel == 4) ? kTRUE : kFALSE;

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
			fTwoPionTrack[i][5] = track->GetIntegratedLength();
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
			fTwoPionITSSigma[i][1] = (Double_t)(fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS, track) == AliPIDResponse::kDetPidOk);
			fTwoPionITSSigma[i][2] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kPion);
			fTwoPionITSSigma[i][3] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kKaon);
			fTwoPionITSSigma[i][4] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kProton);
			fTwoPionITSSigma[i][5] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron);
			fTwoPionITSSigma[i][6] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kPion);
			fTwoPionITSSigma[i][7] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kKaon);
			fTwoPionITSSigma[i][8] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kProton);
			fTwoPionITSSigma[i][9] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kElectron);
			//TPC,TOF,ITS Prob (by P.Buehler)
			fPIDResponse->ComputePIDProbability(AliPIDResponse::kTPC,track,AliPID::kSPECIES,fTwoPionTPCProb[i]);
			fPIDResponse->ComputePIDProbability(AliPIDResponse::kTOF,track,AliPID::kSPECIES,fTwoPionTOFProb[i]);
			fPIDResponse->ComputePIDProbability(AliPIDResponse::kITS,track,AliPID::kSPECIES,fTwoPionITSProb[i]);

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
			fFourPionTrack[i][5] = track->GetIntegratedLength();
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
			fFourPionITSSigma[i][1] = (Double_t)(fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS, track) == AliPIDResponse::kDetPidOk);
			fFourPionITSSigma[i][2] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kPion);
			fFourPionITSSigma[i][3] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kKaon);
			fFourPionITSSigma[i][4] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kProton);
			fFourPionITSSigma[i][5] = fPIDResponse->NumberOfSigmasITS(track,AliPID::kElectron);
			fFourPionITSSigma[i][6] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kPion);
			fFourPionITSSigma[i][7] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kKaon);
			fFourPionITSSigma[i][8] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kProton);
			fFourPionITSSigma[i][9] = fPIDResponse->GetSignalDelta(AliPIDResponse::kITS,track,AliPID::kElectron);
			//TPC,TOF,ITS Prob (by P.Buehler)
			fPIDResponse->ComputePIDProbability(AliPIDResponse::kTPC,track,AliPID::kSPECIES,fFourPionTPCProb[i]);
			fPIDResponse->ComputePIDProbability(AliPIDResponse::kTOF,track,AliPID::kSPECIES,fFourPionTOFProb[i]);
			fPIDResponse->ComputePIDProbability(AliPIDResponse::kITS,track,AliPID::kSPECIES,fFourPionITSProb[i]);
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

	//For ITSSA track----------------------------------------------------------
	if (fCheck2tracks_ITSSA) {
		for (Int_t j = 0; j < 10; j++) {
			for (Int_t i = 0; i < 2; i++) {
				fTwoPionTrack_ITSSA[i][j] = 0;
			}
		}
		for (Int_t i = 0; i < 2; i++) {
			AliESDtrack *track = fESDEvent->GetTrack(indices_ITSSA.At(i));
			if(!track) {
				PostOutputs();
				continue;
			}
			fTwoPionTrack_ITSSA[i][0] = track->Px();
			fTwoPionTrack_ITSSA[i][1] = track->Py();
			fTwoPionTrack_ITSSA[i][2] = track->Pz();
			fTwoPionTrack_ITSSA[i][3] = track->E();
			fTwoPionTrack_ITSSA[i][4] = track->GetSign();
		}
		fCheckTwoPion_ITSSA = kTRUE;
	}
	if (fCheck4tracks_ITSSA) {
		for (Int_t j = 0; j < 10; j++) {
			for (Int_t i = 0; i < 4; i++) {
				fFourPionTrack_ITSSA[i][j] = 0;
			}
		}
		for (Int_t i = 0; i < 4; i++) {
			AliESDtrack *track = fESDEvent->GetTrack(indices_ITSSA.At(i));
			if(!track) {
				PostOutputs();
				continue;
			}
			fFourPionTrack_ITSSA[i][0] = track->Px();
			fFourPionTrack_ITSSA[i][1] = track->Py();
			fFourPionTrack_ITSSA[i][2] = track->Pz();
			fFourPionTrack_ITSSA[i][3] = track->E();
			fFourPionTrack_ITSSA[i][4] = track->GetSign();
		}
		fCheckFourPion_ITSSA = kTRUE;
	}
	//-------------------------------------------------------------------------
	
	indices = 0x0;
	indices_ITSSA = 0x0;

	if(!fIsMC) {//Reduce output size (only for 2tracks and 4tracks)
		if(fSavemode) {
		       fTree->Fill();
		}
 		else {
			if(fCheck2tracks || fCheck4tracks || fCheck2tracks_ITSSA || fCheck4tracks_ITSSA) fTree->Fill();
		}
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
	Double_t cutVertex = 10.;
	if (fIsRun2) cutVertex = 15.;//For run2
	const Bool_t kpriv = kpr0 && (fabs(vertex->GetZ()) < cutVertex);
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
	if (fIsMC) fMCEvent = MCEvent();
	//-------------------------------------------------------------------------

	return kTRUE;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskCDPWA::CheckOnlineTrigger(
		const AliESDEvent *ESDevent
		)
{
	//Fill histgrams for trigger
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
	Char_t tname[50];
	for (Int_t i = 0; i < trigger_no; i++) {
		sprintf(tname,"%s",trigger_name[i]);
		if (ESDevent->IsTriggerClassFired(tname)) {
			fHistTrigger->Fill(kTri_CINT5+i);
			fRunOnline[i]->Fill(ESDevent->GetRunNumber());
		}
	}

	//Reject beam-gas trigger in AD

	//Check that this event is MB_OR online triggered
	if (ESDevent->IsTriggerClassFired("CINT10-B-NOPF-ALLNOTRD") ||
			ESDevent->IsTriggerClassFired("C0SMB-B-NOPF-ALLNOTRD")) return kTRUE;
	else {
		//printf("AliAnalysisTaskCDPWA1::CheckOnlineTrigger not passed!");
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
	if (fIsRun2) {
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
	}

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
		TString st_header = header->IsA()->GetName();
		// Pythia
		if (st_header == "AliGenPythiaEventHeader") {//Pythia6
			fIsPythia = kTRUE;
			fMCprocess = ((AliGenPythiaEventHeader*)header)->ProcessType();
			if (!fIsPythia8) {//Pythia6
				switch(fMCprocess) {
					case 91: fMCprocessType = kBinEL; break;
					case 92: fMCprocessType = kBinSD1; break;
					case 93: fMCprocessType = kBinSD2; break;
					case 94: fMCprocessType = kBinDD; break;
					case 95: fMCprocessType = kBinCD; break;//To be checked!!
					default: fMCprocessType = kBinND; break;
				}
			}
			else {//Pythia8
				switch(fMCprocess) {
					case 102: fMCprocessType = kBinEL; break;
					case 103: fMCprocessType = kBinSD1; break;
					case 104: fMCprocessType = kBinSD2; break;
					case 105: fMCprocessType = kBinDD; break;
					case 106: fMCprocessType = kBinCD; break;
					default: fMCprocessType = kBinND; break;
				}
			}
		}
		// Phojet
		else if (st_header == "AliGenDPMjetEventHeader") {
			fIsPhojet = kTRUE;
			fMCprocess = ((AliGenDPMjetEventHeader*)header)->ProcessType();
			switch(fMCprocess) {
				case 2: fMCprocessType = kBinEL; break;
				case 5: fMCprocessType = kBinSD1; break;
				case 6: fMCprocessType = kBinSD2; break;
				case 7: fMCprocessType = kBinDD; break;
				case 4: fMCprocessType = kBinCD; break;
				default: fMCprocessType = kBinND; break;
			}
		}
		else if (st_header == "AliGenEposEventHeader" || st_header == "AliGenEpos3EventHeader" || st_header == "AliGenHepMCEventHeader") {
			fIsEPOS = kTRUE;
			fMCprocessType = kBinND;
		}
		else {// forDIME/DRgen
			fMCprocessType = kBinCD;
		}
	}
	//-------------------------------------------------------------------------
}
//______________________________________________________________________________
Int_t AliAnalysisTaskCDPWA::DoMCTruth()
{
	// Get generated information from stack and save as TTree if necessary---
	if (!fMCEvent) return -1;
	AliStack* stack = fMCEvent->Stack();
	if (!stack) return -1;
	//-------------------------------------------------------------------------

	// Multiplicity ----------------------------------------------------------
	// determine number of charged physical primaries on the stack
	Int_t nPhysicalPrimaries = 0;
	Int_t nProton = 0;
	Int_t nPrimaries = stack->GetNprimary();
	TLorentzVector lv_vect;
	TLorentzVector lv_sum(0,0,0,0);
	TLorentzVector lv_sum_DD1(0,0,0,0);
	TLorentzVector lv_sum_DD2(0,0,0,0);
	TLorentzVector lv_largestEta;
	TLorentzVector lv_smallestEta;
	TLorentzVector lv_minEta_CD;
	TLorentzVector lv_maxEta_CD;
	Double_t etaMax = -999;
	Double_t etaMin = 999;
	Double_t etaMin_CD = 999;
	Double_t etaMax_CD = -999;
	Double_t minEtaMass_CD = 0.;
	Double_t maxEtaMass_CD = 0.;
	Int_t nFinal = 0;

	TLorentzVector checkSum_i(0,0,0,0);

	for (Int_t iTracks = 0; iTracks < nPrimaries; iTracks++) {
		TParticle *part = (TParticle*)stack->Particle(iTracks);
		if (!part) continue;
		//For the DIME, DRgen--------------------------------------
		if (stack->IsPhysicalPrimary(iTracks) && (part->GetPDG()->Charge() != 0.) && part->GetStatusCode() == 1 && (TMath::Abs(((Int_t)(part->GetPdgCode()))) == 211 || TMath::Abs(((Int_t)(part->GetPdgCode()))) == 321)) {
			nPhysicalPrimaries++;
		}
		if (stack->IsPhysicalPrimary(iTracks) && (part->GetPDG()->Charge() != 0.) && part->GetStatusCode() == 1 && TMath::Abs(((Int_t)(part->GetPdgCode()))) == 2212) {
			nProton++;
		}
		//---------------------------------------------------------

		lv_vect.SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
		if (iTracks == 0 || iTracks ==1) {
			checkSum_i = checkSum_i + lv_vect;
		}
		if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD1) {
//			cout << "N=" << iTracks << ",  " << "PDG = " << part->GetPdgCode() << ",  " << "Status = " << part->GetStatusCode() << ",  " << "FM = " << part->GetFirstMother() << ",  " << "Eta=" << part->Eta() <<endl;
		}
		//Check pseudosystem in pythia
		if (fIsPythia) {
			if (part->GetPdgCode() == 9902210) {//excited proton(p_diffr)
				TParticle *part_mo = (TParticle*)stack->Particle(part->GetFirstMother());
				if (!part_mo) continue;
			       	if (part_mo->GetPdgCode()!=2212) continue;//from incoming proton

				if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD1) {
					fMC_DiffMass_PDG[1]->Fill(lv_vect.M());
				}
				else if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD2) fMC_DiffMass_PDG[2]->Fill(lv_vect.M());
				else if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinDD) {
					fMC_DiffMass_PDG[3]->Fill(lv_vect.M());
					if (part_mo->Pz() > 0) fMC_DiffMass_PDG[4]->Fill(lv_vect.M());
					else fMC_DiffMass_PDG[5]->Fill(lv_vect.M());
				}
			}
			else if (part->GetPdgCode() == 9900110) {//??(rho_diffr)
				if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinCD) fMC_DiffMass_PDG[6]->Fill(lv_vect.M());
			}
		}
		if (fIsPhojet) {
		}

		//Reject incoming proton and keep final states
		//AliPWG0??
//		if (!AliPWG0Helper::IsPrimaryCharged(part,nPrimaries,0)) continue;
		if (!stack->IsPhysicalPrimary(iTracks) || part->GetStatusCode() != 1) continue;// || part->GetStatusCode() != 1) continue;
		nFinal++;


		//For the single diffraction
		//Intacted proton should have largest eta
		if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD2) {
			lv_sum = lv_sum + lv_vect;
			if (part->GetPdgCode() == 2212 && (etaMax < part->Eta())) {//Largest Eta for proton
				etaMax = part->Eta();
				lv_largestEta = lv_vect;
			}
			fMC_Eta[2]->Fill(part->Eta());

		}
		else if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD1) {
			lv_sum = lv_sum + lv_vect;
			if (part->GetPdgCode() == 2212 && (etaMin > part->Eta())) {//Largest Eta for proton
				etaMin = part->Eta();
				lv_largestEta = lv_vect;
			}
			fMC_Eta[1]->Fill(part->Eta());

		}
		//For the double diffraction
		else if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinDD) {
			if (fIsPythia) {
				//Largest Eta for proton
				if (part->GetPdgCode() == 2212 && (etaMax < part->Eta())) {
					etaMax = part->Eta();
					lv_largestEta = lv_vect;
				}
				if (part->GetPdgCode() == 2212 && (etaMin > part->Eta())) {
					etaMin = part->Eta();
					lv_smallestEta = lv_vect;
				}
				lv_sum = lv_sum + lv_vect;
				fMC_Eta[3]->Fill(part->Eta());

				//Find pseudo particle
				TParticle *part_mo = (TParticle*)stack->Particle(part->GetFirstMother());
				if (!part_mo) continue;
				while (part_mo->GetPdgCode() != 9902210) {//Find pseudoparticle
					if (part_mo->GetFirstMother() < 0) break;
					part_mo = (TParticle*)stack->Particle(part_mo->GetFirstMother());
					if (!part_mo) break;
				}
				if (part_mo->GetPdgCode() != 9902210 || part_mo->GetFirstMother() == -1) continue;//Protection
				
				// proton -> pseudoparticle decay
				// part_mo == 9902210
				TParticle *part_top = (TParticle*)stack->Particle(part_mo->GetFirstMother());
				if (!part_top) continue;
				if (part_top->GetPdgCode()!=2212) continue;//protection 

				if (part_top->Pz() > 0) {
					lv_sum_DD1 = lv_sum_DD1 + lv_vect;
					fMC_Eta[4]->Fill(part->Eta());
				}
				else {
					lv_sum_DD2 = lv_sum_DD2 + lv_vect;
					fMC_Eta[5]->Fill(part->Eta());
				}
			}
			if (fIsPhojet) {
				//Find mother particle for proton
				fMC_Eta[3]->Fill(part->Eta());//total
			}
		}
		//For the central diffraction
		else if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinCD) {
			lv_sum = lv_sum + lv_vect;
			if (part->GetPdgCode() == 2212 && etaMin_CD > part->Eta()) {//Minimum eta
				etaMin_CD = part->Eta();
				lv_minEta_CD = lv_vect;
			}
			if (part->GetPdgCode() == 2212 && etaMax_CD < part->Eta()) {//Maximum eta
				etaMax_CD = part->Eta();
				lv_maxEta_CD = lv_vect;
			}
			fMC_Eta[6]->Fill(part->Eta());
		}

		//Store eta informations
		else if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinND) fMC_Eta[0]->Fill(part->Eta());
	}
	//Multiplicity distribution for final states particle
	if (fMCprocessType == AliAnalysisTaskCDPWA::kBinSD1 || fMCprocessType == AliAnalysisTaskCDPWA::kBinSD2) {
		fMult_Gen->Fill(nFinal-1);
		fMult_Gen_Process->Fill(fMCprocessType,nFinal-1);
	}
	else if (fMCprocessType == AliAnalysisTaskCDPWA::kBinDD) {
		fMult_Gen->Fill(nFinal);
		fMult_Gen_Process->Fill(fMCprocessType,nFinal);
	}
	else if (fMCprocessType == AliAnalysisTaskCDPWA::kBinCD) {
		fMult_Gen->Fill(nFinal-2);
		fMult_Gen_Process->Fill(fMCprocessType,nFinal-2);
	}
	else {
		fMult_Gen->Fill(nFinal);
		fMult_Gen_Process->Fill(fMCprocessType,nFinal);
	}

	TLorentzVector lv_diffMass;
	//For single-diffraction
	if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD1 || fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD2) {
		lv_diffMass = lv_sum - lv_largestEta;//For single
		if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD1) {
			fMC_DiffMass[1]->Fill(lv_diffMass.M());
			fComb_diffSystem[0] = lv_diffMass;
			fComb_forwardP[0] = lv_largestEta;
		}
		if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinSD2) {
			fMC_DiffMass[2]->Fill(lv_diffMass.M());
			fComb_diffSystem[0] = lv_diffMass;
			fComb_forwardP[0] = lv_largestEta;
		}
	}

	//For double-diffraction
	else if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinDD) {
		/* Should we subtract this???
		lv_sum = lv_sum - lv_largestEta - lv_smallestEta;
		lv_sum_DD1 = lv_sum_DD1 - lv_largestEta;
		lv_sum_DD2 = lv_sum_DD2 - lv_smallestEta;
		*/
		fMC_DiffMass[3]->Fill(lv_sum.M());
		fMC_DiffMass[4]->Fill(lv_sum_DD1.M());
		fMC_DiffMass[5]->Fill(lv_sum_DD2.M());
		fComb_diffSystem[0] = lv_sum_DD1;
		fComb_diffSystem[1] = lv_sum_DD2;
	}
	//For central diffraction
	else if (fComb_MC_EventProcess == AliAnalysisTaskCDPWA::kBinCD) {
		lv_diffMass = lv_sum - lv_minEta_CD - lv_maxEta_CD;
		fMC_DiffMass[6]->Fill(lv_diffMass.M());
		fComb_diffSystem[0] = lv_diffMass;
		fComb_forwardP[0] = lv_maxEta_CD;
		fComb_forwardP[1] = lv_minEta_CD;
	}
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
	//-------------------------------------------------------------------------

	//Store Generated proton and pion momentum in the TTree--------------------
	//This part is for DIME and DRgen
	if (fIsSaveGen) {
		if (nProton !=2) return -1;
		if (nPhysicalPrimaries != 2 && nPhysicalPrimaries != 4) return -1;
		TParticle *part1 = NULL;
		Int_t n_proton = 0;
		Int_t n_pion = 0;
		for (Int_t j = 0; j < nPrimaries; j++) {
			part1 = (TParticle*)stack->Particle(j);
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
			if (stack->IsPhysicalPrimary(j) && (part1->GetPDG()->Charge() != 0.) && part1->GetStatusCode() == 1 && (TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 211 || TMath::Abs(((Int_t)(part1->GetPdgCode()))) == 321)) {
				fMCGenPionTrack[n_pion][0] = part1->Px();
				fMCGenPionTrack[n_pion][1] = part1->Py();
				fMCGenPionTrack[n_pion][2] = part1->Pz();
				fMCGenPionTrack[n_pion][3] = part1->Energy();
				fMCGenPionTrack[n_pion][4] = part1->GetPdgCode();
				n_pion++;
				continue;
			}
		}
	}
	//-------------------------------------------------------------------------
	return nPhysicalPrimaries;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::DoCombStudy(const AliESDAD *ad, const AliESDVZERO *vzero)
{
	//Save AD, V0 charge and adc to the TTree
	//AD
	if (fIsRun2) {
		fComb_AD_Time_Mean[0] = ad->GetADATime();
		fComb_AD_Time_Mean[1] = ad->GetADCTime();
		for (Int_t i = 0; i < 16; i++) {
			fComb_AD_Time[i] = ad->GetTime(i);
			fComb_AD_ADC[i] = ad->GetAdc(i);
		}
	}
	//VZERO
	fComb_V0_Time_Mean[0] = vzero->GetV0ATime();
	fComb_V0_Time_Mean[1] = vzero->GetV0CTime();
	for (Int_t i = 0; i < 64; i++) {
		fComb_V0_Time[i] = vzero->GetTime(i);
		fComb_V0_ADC[i] = vzero->GetAdc(i);
	}
	return;
}
//______________________________________________________________________________
void AliAnalysisTaskCDPWA::Terminate(Option_t *)
{
}
//______________________________________________________________________________
