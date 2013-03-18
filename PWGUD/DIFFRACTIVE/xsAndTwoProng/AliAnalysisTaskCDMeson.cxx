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

#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <THnSparse.h>

#include <TRandom3.h>

#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliESDInputHandler.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliStack.h"

#include "AliCDMesonBase.h"
#include "AliCDMesonUtils.h"
#include "AliCDMesonTracks.h"
#include "AliAnalysisTaskCDMeson.h"


//------------------------------------------------------------------------------
AliAnalysisTaskCDMeson::AliAnalysisTaskCDMeson(const char* name, Long_t state):
	AliAnalysisTaskSE(name)
	, fDoAOD(kFALSE)
	, fAnalysisStatus(state)
	, fAnalysisStatusString(new TObjString(Form("0x%lx", fAnalysisStatus)))
	, fNoGapFraction(0.01)
	, fReducedGapFraction(0.02)
	, fMaxVtxDst(0.5) // value to be checked with the vertex study histograms

	, fESDEvent(0x0)
	, fAODEvent(0x0)
	, fPIDResponse(0x0)
	, fPhysicsSelection(0x0)
	, fTracks(new AliCDMesonTracks())
	, fVtxDst(-1.)
	, fVtxZ(-20)
	, fResidualTracks(0)
	, fResidualTracklets(0)
	, fMCprocessType(0)
	, fMCprocess(-1)

	, fRun(-999)
	, fPIDmode(0)

	, fTheta(0.)
	, fPhi(0.)
	, fMass(0.)
	, fCurrentGapCondition(0)

	, fGapRun(0x0)
	, fHist(0x0)
	, fThnMother(0x0)
	, fThnMotherSoft(0x0)
	, fThnMultiplicity(0x0)
	, fThnMC(0x0)
	, fThnEmptyEvents(0x0)
	, fPWAtree(0x0)

	, fv0ntrk(0x0)
	, fv0fmdntrk(0x0)
	, fv0fmdspdntrk(0x0)
	, fv0fmdspdtpcntrk(0x0)

	, fv0Rmntrk(0x0)
	, fv0fmdRmntrk(0x0)
	, fv0fmdspdRmntrk(0x0)
	, fv0fmdspdtpcRmntrk(0x0)

	, fMultStudy(0x0)
	, fMultStudyV0dg(0x0)
	, fMultStudyV0FMDdg(0x0)
	, fMultStudyV0FMDSPDdg(0x0)
	, fMultStudyV0FMDSPDTPCdg(0x0)

	, fMultResponseMC(0x0)
	, fMultResponseV0dgMC(0x0)
	, fMultResponseV0FMDdgMC(0x0)
	, fMultResponseV0FMDSPDdgMC(0x0)
	, fMultResponseV0FMDSPDTPCdgMC(0x0)

	, fMultRegionsMC(0x0)

	, fhspdV0dg(0x0)
	, fhspdV0FMDdg(0x0)
	, fhspdV0FMDSPDdg(0x0)
	, fhspdV0FMDSPDTPCdg(0x0)
	, fhspdAfterCuts(0x0)

	, fGapResponseMCv0Dg(0x0)
	, fGapResponseMCv0fmdDg(0x0)
	, fGapResponseMCv0fmdspdDg(0x0)
	, fGapResponseMCv0fmdspdtpcDg(0x0)

	, fhspd(0x0)
	, fhfo(0x0)
	, fhpriVtxDist(0x0)
	, fhpriVtxPos(0x0)
	, fhpriv(0x0)
	, fhntrk(0x0)
	, fhStatsFlow(0x0)
	, fhFOchans(0x0)

	, fVZEROhists(0x0)

	, fv0Map(0x0)
	, fv0fmdMap(0x0)
	, fv0fmdspdMap(0x0)
	, fv0fmdspdtpcMap(0x0)

	, fv0MapCutted(0x0)
	, fv0fmdMapCutted(0x0)
	, fv0fmdspdMapCutted(0x0)
	, fv0fmdspdtpcMapCutted(0x0)

	, fHitMapSPDinner(0x0)
	, fHitMapSPDouter(0x0)
	, fHitMapSPDtrklt(0x0)
	, fHitMapFMDa(0x0)
	, fHitMapFMDc(0x0)

	, fFMDsum1I(0x0)
	, fFMDsum2I(0x0)
	, fFMDsum2O(0x0)
	, fFMDsum3I(0x0)
	, fFMDsum3O(0x0)

	, fPriVtxX(0x0)
	, fPriVtxY(0x0)
	, fPriVtxZ(0x0)
	, fPriVtxDst(0x0)
	, fPriVtxDstV0dg(0x0)
	, fPriVtxDstV0FMDdg(0x0)
	, fPriVtxDstV0FMDSPDdg(0x0)
	, fPriVtxDstV0FMDSPDTPCdg(0x0)

	, fTPCGapDCAaSide(0x0)
	, fTPCGapDCAcSide(0x0)

	, fComb2trkPIDuls(0x0)
	, fComb2trkPIDls(0x0)
	, fComb2trkPIDulsDG(0x0)
	, fComb2trkPIDlsDG(0x0)
	, fComb2trkPIDulsMC(0x0)
	, fComb2trkPIDlsMC(0x0)
	, fComb2trkPIDulsDGmc(0x0)
	, fComb2trkPIDlsDGmc(0x0)
	, fMCProcessUls(0x0)
	, fMCProcessLs(0x0)

	, fMAllTrackMass(0x0)
{
	//
	// standard constructor (the one which should be used)
	//
	// slot in TaskSE must start from 1
	Int_t iOutputSlot = 1;
	DefineOutput(iOutputSlot++, TList::Class());
	if (!(fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) {
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitTHnMother)) {
			DefineOutput(iOutputSlot++, THnSparse::Class());
		}
		if (!fAnalysisStatus
		    || ((fAnalysisStatus & AliCDMesonBase::kBitSoftTracks)
		        && (fAnalysisStatus & AliCDMesonBase::kBitTHnMother))) {
			DefineOutput(iOutputSlot++, THnSparse::Class());
		}
		if (fAnalysisStatus & AliCDMesonBase::kBitMultStudy) {
			DefineOutput(iOutputSlot++, THnSparse::Class());
		}
		if (fAnalysisStatus & AliCDMesonBase::kBitTHnMC) {
			DefineOutput(iOutputSlot++, THnSparse::Class());
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitPWAtree)) {
			DefineOutput(iOutputSlot++, TTree::Class());
		}
	}
	else { // create empty event THnSparse
		DefineOutput(iOutputSlot++, THnSparse::Class());
	}

	if (fAnalysisStatus & AliCDMesonBase::kBitEEStudy) { // empty event study
		fPhysicsSelection = new AliPhysicsSelection;
		fPhysicsSelection->SetSkipTriggerClassSelection();
		// accept all (A,C,E and I) triggers of a certain class
	}

	for (Int_t iGap = 0; iGap < kMax; ++iGap) {
		fGapInformation[iGap] = 0;
		fGapInformationWCent[iGap] = 0;
	}

	for (Int_t i = 0; i < 3; ++i) {
		fMomentum[i] = 0.;
	}

	// reset track pair pointers
	fTrkPair[0] = 0x0;
	fTrkPair[1] = 0x0;
}


//------------------------------------------------------------------------------
AliAnalysisTaskCDMeson::AliAnalysisTaskCDMeson():
	AliAnalysisTaskSE()
	, fDoAOD(kFALSE)
	, fAnalysisStatus(0x0)
	, fAnalysisStatusString(0x0)
	, fNoGapFraction(0.01)
	, fReducedGapFraction(0.02)
	, fMaxVtxDst(0.5)

	, fESDEvent(0x0)
	, fAODEvent(0x0)
	, fPIDResponse(0x0)
	, fPhysicsSelection(0x0)
	, fTracks(0x0)
	, fVtxDst(-1.)
	, fVtxZ(-20)
	, fResidualTracks(0)
	, fResidualTracklets(0)
	, fMCprocessType(0)
	, fMCprocess(-1)

	, fRun(-999)
	, fPIDmode(0)

	, fTheta(0.)
	, fPhi(0.)
	, fMass(0.)
	, fCurrentGapCondition(0)

	, fGapRun(0x0)
	, fHist(0x0)
	, fThnMother(0x0)
	, fThnMotherSoft(0x0)
	, fThnMultiplicity(0x0)
	, fThnMC(0x0)
	, fThnEmptyEvents(0x0)
	, fPWAtree(0x0)

	, fv0ntrk(0x0)
	, fv0fmdntrk(0x0)
	, fv0fmdspdntrk(0x0)
	, fv0fmdspdtpcntrk(0x0)

	, fv0Rmntrk(0x0)
	, fv0fmdRmntrk(0x0)
	, fv0fmdspdRmntrk(0x0)
	, fv0fmdspdtpcRmntrk(0x0)

	, fMultStudy(0x0)
	, fMultStudyV0dg(0x0)
	, fMultStudyV0FMDdg(0x0)
	, fMultStudyV0FMDSPDdg(0x0)
	, fMultStudyV0FMDSPDTPCdg(0x0)

	, fMultResponseMC(0x0)
	, fMultResponseV0dgMC(0x0)
	, fMultResponseV0FMDdgMC(0x0)
	, fMultResponseV0FMDSPDdgMC(0x0)
	, fMultResponseV0FMDSPDTPCdgMC(0x0)

	, fMultRegionsMC(0x0)

	, fhspdV0dg(0x0)
	, fhspdV0FMDdg(0x0)
	, fhspdV0FMDSPDdg(0x0)
	, fhspdV0FMDSPDTPCdg(0x0)
	, fhspdAfterCuts(0x0)

	, fGapResponseMCv0Dg(0x0)
	, fGapResponseMCv0fmdDg(0x0)
	, fGapResponseMCv0fmdspdDg(0x0)
	, fGapResponseMCv0fmdspdtpcDg(0x0)

	, fhspd(0x0)
	, fhfo(0x0)
	, fhpriVtxDist(0x0)
	, fhpriVtxPos(0x0)
	, fhpriv(0x0)
	, fhntrk(0x0)
	, fhStatsFlow(0x0)
	, fhFOchans(0x0)

	, fVZEROhists(0x0)

	, fv0Map(0x0)
	, fv0fmdMap(0x0)
	, fv0fmdspdMap(0x0)
	, fv0fmdspdtpcMap(0x0)

	, fv0MapCutted(0x0)
	, fv0fmdMapCutted(0x0)
	, fv0fmdspdMapCutted(0x0)
	, fv0fmdspdtpcMapCutted(0x0)

	, fHitMapSPDinner(0x0)
	, fHitMapSPDouter(0x0)
	, fHitMapSPDtrklt(0x0)
	, fHitMapFMDa(0x0)
	, fHitMapFMDc(0x0)

	, fFMDsum1I(0x0)
	, fFMDsum2I(0x0)
	, fFMDsum2O(0x0)
	, fFMDsum3I(0x0)
	, fFMDsum3O(0x0)

	, fPriVtxX(0x0)
	, fPriVtxY(0x0)
	, fPriVtxZ(0x0)
	, fPriVtxDst(0x0)
	, fPriVtxDstV0dg(0x0)
	, fPriVtxDstV0FMDdg(0x0)
	, fPriVtxDstV0FMDSPDdg(0x0)
	, fPriVtxDstV0FMDSPDTPCdg(0x0)

	, fTPCGapDCAaSide(0x0)
	, fTPCGapDCAcSide(0x0)

	, fComb2trkPIDuls(0x0)
	, fComb2trkPIDls(0x0)
	, fComb2trkPIDulsDG(0x0)
	, fComb2trkPIDlsDG(0x0)
	, fComb2trkPIDulsMC(0x0)
	, fComb2trkPIDlsMC(0x0)
	, fComb2trkPIDulsDGmc(0x0)
	, fComb2trkPIDlsDGmc(0x0)
	, fMCProcessUls(0x0)
	, fMCProcessLs(0x0)

	, fMAllTrackMass(0x0)
{
	//
	// default constructor (should not be used)
	//

	for (Int_t iGap = 0; iGap < kMax; ++iGap) {
		fGapInformation[iGap] = 0;
		fGapInformationWCent[iGap] = 0;
	}

	for (Int_t i = 0; i < 3; ++i) {
		fMomentum[i] = 0.;
	}

	fTrkPair[0] = 0x0;
	fTrkPair[1] = 0x0;
}


//------------------------------------------------------------------------------
AliAnalysisTaskCDMeson::~AliAnalysisTaskCDMeson()
{
	//
	//Destructor
	//
	//delete fESDEvent;
	//delete fESDpid;

	if (!(fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) { // normal operation
		if ((!fAnalysisStatus || (fAnalysisStatus && AliCDMesonBase::kBitTHnMother))
		    && fThnMother
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fThnMother;
			fThnMother = 0x0;
		}
		if ((!fAnalysisStatus
		     || ((fAnalysisStatus && AliCDMesonBase::kBitSoftTracks)
		         && (fAnalysisStatus && AliCDMesonBase::kBitTHnMother)))
		    && fThnMotherSoft
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fThnMotherSoft;
			fThnMotherSoft = 0x0;
		}
		if ((!fAnalysisStatus
		     || (fAnalysisStatus && AliCDMesonBase::kBitMultStudy))
		    && fThnMultiplicity
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fThnMultiplicity;
			fThnMultiplicity = 0x0;
		}
		if ((!fAnalysisStatus || (fAnalysisStatus && AliCDMesonBase::kBitTHnMC))
		    && fThnMC
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fThnMC;
			fThnMC = 0x0;
		}
	}
	else { // empty event study
		if (fPhysicsSelection
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fPhysicsSelection;
			fPhysicsSelection = 0x0;
		}
		if (fThnEmptyEvents
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fThnEmptyEvents;
			fThnEmptyEvents = 0x0;
		}
	}

	/*
	if (fHist
	    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
	        AliAnalysisManager::kProofAnalysis)) {
		fHist->Clear();
		delete fHist;
		fHist = 0x0;
	}
	*/

	if (fTracks
	    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
	        AliAnalysisManager::kProofAnalysis)) {
		delete fTracks;
		fTracks = 0x0;
	}

	if (fVZEROhists) {
		delete fVZEROhists;
		fVZEROhists = 0x0;
	}
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::UserCreateOutputObjects()
{
	//
	//CreateOutputObjects
	//

	//= THnSparse ================================================================
	if (!(fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) { // normal operation
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitTHnMother)) {
			fThnMother = AliCDMesonBase::GetThnMother(/*"CDMeson_Mother"*/);
		}
		if (!fAnalysisStatus
		    || ((fAnalysisStatus & AliCDMesonBase::kBitSoftTracks)
		        && (fAnalysisStatus & AliCDMesonBase::kBitTHnMother))) {
			fThnMotherSoft =  AliCDMesonBase::GetThnMother("CDMeson_MotherSoft");
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitMultStudy)){
			fThnMultiplicity =  AliCDMesonBase::GetThnMultiplicity();
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitTHnMC)) {
			fThnMC = AliCDMesonBase::GetThnMother("CDMeson_MotherMC");
		}
	}
	else { // empty event studies
		fThnEmptyEvents = AliCDMesonBase::GetThnEmptyEvents();
	}

	//= TList for Histograms =====================================================
	fHist = new TList;
	fHist->SetOwner(); // ensures that the histograms are all deleted on exit!

	//fHist->Add(fAnalysisStatusString); // WARNING: CRASHES MERGING

	//= GAP RUN =
	if (!(fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) { // normal operation
		// LHC10 data
		// runa1=114649, runb1=117630, runc0=117631, runc1=121526, rund0=121527,
		// rund1=126460, rune0=126461, rune1 = 130930; runf0=130931
		// all checked on elog
		// LHC11 data
		//const Int_t rune0 = 160676;
		//const Int_t runf1 = 165462;

		// LHC10 + LHC11 data
		const Int_t runb0 = 114650; // first run of LHC10b
		const Int_t runf1 = 165462; // last run of LHC11f

		const Int_t run0 = runb0-1, run1 = runf1 + 1;

		const Int_t nrun = (Int_t)(run1-run0);

		const Int_t bins[2] = {nrun, (Double_t)AliCDMesonBase::kBitGapMax - 1.};
		const Double_t xmin[2] = {(Double_t)run0, 1.};
		const Double_t xmax[2] = {
			(Double_t)run1, (Double_t)AliCDMesonBase::kBitGapMax
		};

		fGapRun = new THnSparseI("GapRun", ";run number;gap condition", 2, bins,
		                         xmin, xmax);
		fHist->Add(fGapRun);
	}

	//= MULTIPLICITY PER GAP CONDITION =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitMultPerGapHists)) {
		fv0ntrk = new TH2D("b00_v0ntrk", ";number of tracks;gap condition",
		                   80, 0., 80., 4, 1., 5.);
		//x: ntrk; y: V0
		fHist->Add(fv0ntrk);

		fv0fmdntrk = new TH2D("b01_v0fmdntrk", ";number of tracks;gap condition",
		                      80, 0., 80., 4, 1., 5.);
		//x: ntrk; y: V0FMD
		fHist->Add(fv0fmdntrk);

		fv0fmdspdntrk =
			new TH2D("b02_v0fmdspdntrk", ";number of tracks;gap condition",
		         80, 0., 80., 4, 1., 5.);
		//x: ntrk; y: V0FMDSPD
		fHist->Add(fv0fmdspdntrk);

		fv0fmdspdtpcntrk =
			new TH2D("b03_v0fmdspdtpcntrk", ";number of tracks;gap condition",
			         80, 0., 80., 4, 1., 5.);
		//x: ntrk; y: V0FMDSPDTPC
		fHist->Add(fv0fmdspdtpcntrk);
	}

	//= MULTIPLICITY REMOVED PER GAP CONDITION =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitRmMultPerGapHists)) {
		fv0Rmntrk = new TH2F("b05_v0Rmntrk",
		                     ";number of removed tracks due to cuts;gap condition",
		                     80, 0., 80., 4, 1., 5.);
		//x: ntrk; y: V0
		fHist->Add(fv0Rmntrk);

		fv0fmdRmntrk =
			new TH2F("b06_v0fmdRmntrk",
			         ";number of removed tracks due to cuts;gap condition",
			         80, 0., 80., 4, 1., 5.);
		//x: ntrk; y: V0FMD
		fHist->Add(fv0fmdRmntrk);

		fv0fmdspdRmntrk =
			new TH2F("b07_v0fmdspdRmntrk",
			         ";number of removed tracks due to cuts;gap condition",
			         80, 0., 80., 4, 1., 5.);
		//x: ntrk; y: V0FMDSPD
		fHist->Add(fv0fmdspdRmntrk);

		fv0fmdspdtpcRmntrk =
			new TH2F("b08_v0fmdspdtpcRmntrk",
			         ";number of removed tracks due to cuts;gap condition",
			         80, 0., 80., 4, 1., 5.);
		//x: ntrk; y: V0FMDSPDTPC
		fHist->Add(fv0fmdspdtpcRmntrk);
	}

	//= SOFT TRACK INFLUENCE =
	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitSoftTracks)) {
		fMultStudy = new TH2F("b10_MultStudySoftInfluence",
		                      ";multiplicity;multiplicity (including ITSsa)",
		                      20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultStudy);

		fMultStudyV0dg =
			new TH2F("b11_MultStudySoftInfluence_V0_DoubleGap",
			         ";multiplicity;multiplicity (including ITSsa)",
			         20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultStudyV0dg);

		fMultStudyV0FMDdg =
			new TH2F("b12_MultStudySoftInfluence_V0FMD_DoubleGap",
			         ";multiplicity;multiplicity (including ITSsa)",
			         20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultStudyV0FMDdg);

		fMultStudyV0FMDSPDdg
			= new TH2F("b13_MultStudySoftInfluence_V0FMDSPD_DoubleGap",
			           ";multiplicity;multiplicity (including ITSsa)",
			           20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultStudyV0FMDSPDdg);

		fMultStudyV0FMDSPDTPCdg
			= new TH2F("b14_MultStudySoftInfluence_V0FMDSPDTPC_DoubleGap",
			           ";multiplicity;multiplicity (including ITSsa)",
			           20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultStudyV0FMDSPDTPCdg);
	}

	//= MULTIPLICITY RESPONSE DATA VS MC =
	if (fAnalysisStatus & AliCDMesonBase::kBitMultResponseMC) {
		fMultResponseMC
			= new TH2F("b15_MultResponseMC",
			           ";multiplicity from MC truth; reconstructed multiplicity",
			           20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultResponseMC);

		fMultResponseV0dgMC
			= new TH2F("b16_MultResponseV0dgMC",
			           ";multiplicity from MC truth; reconstructed multiplicity",
			           20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultResponseV0dgMC);

		fMultResponseV0FMDdgMC
			= new TH2F("b17_MultResponseV0FMDdgMC",
			           ";multiplicity from MC truth; reconstructed multiplicity",
			           20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultResponseV0FMDdgMC);

		fMultResponseV0FMDSPDdgMC
			= new TH2F("b18_MultResponseV0FMDSPDdgMC",
			           ";multiplicity from MC truth; reconstructed multiplicity",
			           20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultResponseV0FMDSPDdgMC);

		fMultResponseV0FMDSPDTPCdgMC
			= new TH2F("b19_MultResponseV0FMDSPDTPCdgMC",
			           ";multiplicity from MC truth; reconstructed multiplicity",
			           20, 0., 20., 20, 0., 20.);
		fHist->Add(fMultResponseV0FMDSPDTPCdgMC);

		fMultRegionsMC
			= new TH2F("b20_fMultRegionsMC", ";multiplicity;",
			           50, 0., 50., 4, 0., 4.);
		TAxis* a = fMultRegionsMC->GetYaxis();
		a->SetBinLabel(1, "A-side");
		a->SetBinLabel(2, "Center");
		a->SetBinLabel(3, "C-side");
		a->SetBinLabel(4, "A and C side");
		fHist->Add(fMultRegionsMC);

		fGapResponseMCv0Dg =
			new TH1I("b26_fGapResponseMCv0Dg", ";;counts", 3, 0., 3.);
		a = fGapResponseMCv0Dg->GetXaxis();
		a->SetBinLabel(1, "DG in MC, no DG in data");
		a->SetBinLabel(2, "DG in both MC and data");
		a->SetBinLabel(3, "DG in data, no DG in MC");
		fHist->Add(fGapResponseMCv0Dg);
		fGapResponseMCv0fmdDg =
			new TH1I("b27_fGapResponseMCv0fmdDg", ";;counts", 3, 0., 3.);
		a = fGapResponseMCv0fmdDg->GetXaxis();
		a->SetBinLabel(1, "DG in MC, no DG in data");
		a->SetBinLabel(2, "DG in both MC and data");
		a->SetBinLabel(3, "DG in data, no DG in MC");
		fHist->Add(fGapResponseMCv0fmdDg);
		fGapResponseMCv0fmdspdDg =
			new TH1I("b28_fGapResponseMCv0fmdspdDg", ";;counts", 3, 0., 3.);
		a = fGapResponseMCv0fmdspdDg->GetXaxis();
		a->SetBinLabel(1, "DG in MC, no DG in data");
		a->SetBinLabel(2, "DG in both MC and data");
		a->SetBinLabel(3, "DG in data, no DG in MC");
		fHist->Add(fGapResponseMCv0fmdspdDg);
		fGapResponseMCv0fmdspdtpcDg
			= new TH1I("b29_fGapResponseMCv0fmdspdtpcDg", ";;counts", 3, 0., 3.);
		a = fGapResponseMCv0fmdspdtpcDg->GetXaxis();
		a->SetBinLabel(1, "DG in MC, no DG in data");
		a->SetBinLabel(2, "DG in both MC and data");
		a->SetBinLabel(3, "DG in data, no DG in MC");
		fHist->Add(fGapResponseMCv0fmdspdtpcDg);
	}

	//= FAST-OR MULTIPLICITY STUDY =
	if (fAnalysisStatus & AliCDMesonBase::kBitFastORmultStudy) {
		fhspdV0dg =
			new TH1I("b21_spd_V0dg",
			         "V0 double gap;fired SPD FastOR chips; counts", 50, 0., 50.);
		fHist->Add(fhspdV0dg);
		fhspdV0FMDdg =
			new TH1I("b22_spd_V0FMDdg",
			         "V0-FMD double gap;fired SPD FastOR chips; counts", 50, 0., 50.);
		fHist->Add(fhspdV0FMDdg);
		fhspdV0FMDSPDdg =
			new TH1I("b23_spd_V0FMDSPDdg",
			         "V0-FMD-SPD double gap;fired SPD FastOR chips; counts",
			         50, 0., 50.);
		fHist->Add(fhspdV0FMDSPDdg);
		fhspdV0FMDSPDTPCdg =
			new TH1I("b24_spd_V0FMDSPDTPCdg",
			         "V0-FMD-SPD-TPC double gap;fired SPD FastOR chips; counts",
			         50, 0., 50.);
		fHist->Add(fhspdV0FMDSPDTPCdg);
		fhspdAfterCuts =
			new TH1I("b25_spd_afterCuts",
			         "after cuts;fired SPD FastOR chips; counts",
			         50, 0., 50.);
		fHist->Add(fhspdAfterCuts);
	}

	//= STATISTICS FLOW =
	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitStatsFlow)) {
		fhspd = new TH1I("a00_spd", ";fired SPD FastOR chips; counts", 50, 0., 50.);
		fHist->Add(fhspd);

		if (fAnalysisStatus & AliCDMesonBase::kBitFastORStudy) {
			fhfo =
				new TH2I("a05_fo", ";fired chips inner layer;fired chips outer layer",
				         50, 0., 50., 50, 0., 50.);
			fHist->Add(fhfo);

			fhFOchans =
				new TH1F("a06_foChans", ";channel key; counts", 1200, 0., 1200.);
			fHist->Add(fhFOchans);
		}

		fhpriVtxDist = new TH1I("a08_priVtxDist", ";vertex position (cm); counts",
		                        601, -30.05, 30.05);
		fHist->Add(fhpriVtxDist);

		fhpriVtxPos = new TH1I("a09_priVtxPos", ";vertex position (boolean);counts",
		                  2, 0., 2.);
		fhpriVtxPos->GetXaxis()->SetBinLabel(1,"inside");
		fhpriVtxPos->GetXaxis()->SetBinLabel(2,"outside");
		fHist->Add(fhpriVtxPos);


		fhpriv = new TH1I("a10_priv", ";vertex quality (boolean);counts",
		                  2, 0., 2.);
		fhpriv->GetXaxis()->SetBinLabel(1,"bad vertex");
		fhpriv->GetXaxis()->SetBinLabel(2,"good vertex");
		fHist->Add(fhpriv);

		fhntrk = new TH1I("a20_ntrk", ";track multiplicity after cuts;counts",
		                  10, 0., 10.);
		fHist->Add(fhntrk);

		fhStatsFlow = AliCDMesonBase::GetHistStatsFlow();
		fHist->Add(fhStatsFlow);
	}

	//= ETA-PHI MAPS FOR TRACKS =
	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitEtaPhiMaps)) {
		fv0Map = new TH2F("d0_v0Map", ";#eta;#phi", 40, -2., 2.,
		                  60, 0., TMath::TwoPi());
		fHist->Add(fv0Map);

		fv0fmdMap = new TH2F("d0_v0fmdMap", ";#eta;#phi", 40, -2., 2.,
		                     60, 0., TMath::TwoPi());
		fHist->Add(fv0fmdMap);

		fv0fmdspdMap = new TH2F("d0_v0fmdspdMap", ";#eta;#phi", 40, -2., 2.,
		                        60, 0., TMath::TwoPi());
		fHist->Add(fv0fmdspdMap);

		fv0fmdspdtpcMap = new TH2F("d0_v0fmdspdtpcMap", ";#eta;#phi", 30, -1.5, 1.5,
		                           60, 0., TMath::TwoPi());
		fHist->Add(fv0fmdspdtpcMap);
	}

	//= ETA-PHI MAPS FOR CUTTED TRACKS =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitEtaPhiMapsWithCuts)) {
		fv0MapCutted = new TH2F("d0_v0MapCutted", ";#eta;#phi", 24, -1.2, 1.2,
		                        60, 0., TMath::TwoPi());
		fHist->Add(fv0MapCutted);

		fv0fmdMapCutted = new TH2F("d0_v0fmdMapCutted", ";#eta;#phi", 24, -1.2, 1.2,
		                           60, 0., TMath::TwoPi());
		fHist->Add(fv0fmdMapCutted);

		fv0fmdspdMapCutted = new TH2F("d0_v0fmdspdMapCutted", ";#eta;#phi",
		                              24, -1.2, 1.2, 60, 0., TMath::TwoPi());
		fHist->Add(fv0fmdspdMapCutted);

		fv0fmdspdtpcMapCutted = new TH2F("d0_v0fmdspdtpcMapCutted", ";#eta;#phi",
		                                 24, -1.2, 1.2, 60, 0., TMath::TwoPi());
		fHist->Add(fv0fmdspdtpcMapCutted);
	}

	//= SPD HIT MAPS =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitHitMapSPD) ||
	    (fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) {
		fHitMapSPDinner = new TH2F("d1_HitMapSPDinner", ";#eta;#phi", 72, -3.6, 3.6,
		                           20, 0., TMath::TwoPi());
		fHist->Add(fHitMapSPDinner);

		fHitMapSPDouter = new TH2F("d2_HitMapSPDouter", ";#eta;#phi", 48, -2.4, 2.4,
		                           40, 0., TMath::TwoPi());
		fHist->Add(fHitMapSPDouter);

		fHitMapSPDtrklt = new TH2F("d3_HitMapSPDtrklt", ";#eta;#phi", 60, -3., 3.,
		                           72, 0., TMath::TwoPi());
		fHist->Add(fHitMapSPDtrklt);
	}

	//= FMD HIT MAPS =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitHitMapFMD) ||
	    (fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) {
		fHitMapFMDa = new TH2F("d4_HitMapFMDa", ";#eta;#phi",
		                       36, 1.5, 5.1, 40, 0., TMath::TwoPi());
		fHist->Add(fHitMapFMDa);
		fHitMapFMDc = new TH2F("d4_HitMapFMDc", ";#eta;#phi",
		                       20, -3.5, -1.5, 40, 0., TMath::TwoPi());
		fHist->Add(fHitMapFMDc);
	}

	//= FMD SUMMATION =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitFMDsum) ||
	    (fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) {
		fFMDsum1I = new TH1F("d5_FMDsum1I", ";sum;counts", 50, 0., 50.);
		fHist->Add(fFMDsum1I);
		fFMDsum2I = new TH1F("d5_FMDsum2I", ";sum;counts", 50, 0., 50.);
		fHist->Add(fFMDsum2I);
		fFMDsum2O = new TH1F("d5_FMDsum2O", ";sum;counts", 50, 0., 50.);
		fHist->Add(fFMDsum2O);
		fFMDsum3I = new TH1F("d5_FMDsum3I", ";sum;counts", 50, 0., 50.);
		fHist->Add(fFMDsum3I);
		fFMDsum3O = new TH1F("d5_FMDsum3O", ";sum;counts", 50, 0., 50.);
		fHist->Add(fFMDsum3O);
	}

	//= VERTEX STUDIES =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitVtxStudies)) {
		fPriVtxX = new TH1I("e0_PriVtxX", ";x (cm);counts", 100, -10., 10.);
		fHist->Add(fPriVtxX);
		fPriVtxY = new TH1I("e0_PriVtxY", ";y (cm);counts", 100, -10., 10.);
		fHist->Add(fPriVtxY);
		fPriVtxZ = new TH1I("e0_PriVtxZ", ";z (cm);counts", 100, -10., 10.);
		fHist->Add(fPriVtxZ);

		fPriVtxDst = new TH1I("e1_PriVtxDst", ";dst (cm);counts", 101, 0., 2.525);
		fHist->Add(fPriVtxDst);

		fPriVtxDstV0dg =
			new TH1I("e2_PriVtxDstV0dg", ";dst (cm);counts", 101, 0., 2.525);
		fHist->Add(fPriVtxDstV0dg);

		fPriVtxDstV0FMDdg =
			new TH1I("e3_PriVtxDstV0FMDdg", ";dst (cm);counts", 101, 0., 2.525);
		fHist->Add(fPriVtxDstV0FMDdg);

		fPriVtxDstV0FMDSPDdg =
			new TH1I("e4_PriVtxDstV0FMDSPDdg", ";dst (cm);counts", 101, 0., 2.525);
		fHist->Add(fPriVtxDstV0FMDSPDdg);

		fPriVtxDstV0FMDSPDTPCdg =
			new TH1I("e5_PriVtxDstV0FMDSPDTPCdg", ";dst (cm);counts", 101, 0., 2.525);
		fHist->Add(fPriVtxDstV0FMDSPDTPCdg);
	}

	//= TPC GAP STUDY =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitTPCGapStudy)) {
		if (!fDoAOD) {
			fTPCGapDCAaSide = new TH2F("i00_TPCGapDCAaSide", ";dca_{xy};dca_{z}",
			                           2000, 0., 1., 1000, 0., 5.);
			fHist->Add(fTPCGapDCAaSide);
			fTPCGapDCAcSide = new TH2F("i01_TPCGapDCAcSide", ";dca_{xy};dca_{z}",
			                           2000, 0., 1., 1000, 0., 5.);
			fHist->Add(fTPCGapDCAcSide);
		}
	}

	//= VZERO TRIGGER STUDY =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitVZEROStudy)) {
		if (!fDoAOD) { // not yet possible for AODs
			fVZEROhists = AliCDMesonBase::GetHistVZEROStudies(fHist);
		}
	}

	//= PID MATRIX FOR DATA AND MC =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitPIDStudy)) {
		fComb2trkPIDuls = AliCDMesonBase::GetHistPIDStudies("f0_Comb2trkPIDuls");
		fHist->Add(fComb2trkPIDuls);

		fComb2trkPIDls = AliCDMesonBase::GetHistPIDStudies("f1_Comb2trkPIDls");
		fHist->Add(fComb2trkPIDls);

		fComb2trkPIDulsDG =
			AliCDMesonBase::GetHistPIDStudies("f2_Comb2trkPIDulsDG");
		fHist->Add(fComb2trkPIDulsDG);

		fComb2trkPIDlsDG = AliCDMesonBase::GetHistPIDStudies("f3_Comb2trkPIDlsDG");
		fHist->Add(fComb2trkPIDlsDG);

		if (fAnalysisStatus & AliCDMesonBase::kBitTHnMC) {
			fComb2trkPIDulsMC =
				AliCDMesonBase::GetHistPIDStudies("f4_Comb2trkPIDulsMC");
			fHist->Add(fComb2trkPIDulsMC);

			fComb2trkPIDlsMC = AliCDMesonBase::GetHistPIDStudies("f5_Comb2trkPIDlsMC");
			fHist->Add(fComb2trkPIDlsMC);

			fComb2trkPIDulsDGmc =
				AliCDMesonBase::GetHistPIDStudies("f6_Comb2trkPIDulsDGmc");
			fHist->Add(fComb2trkPIDulsDGmc);

			fComb2trkPIDlsDGmc =
				AliCDMesonBase::GetHistPIDStudies("f7_Comb2trkPIDlsDGmc");
			fHist->Add(fComb2trkPIDlsDGmc);
		}
	}

	//= MC PROCESS IDs =
	if (!fAnalysisStatus || fAnalysisStatus & AliCDMesonBase::kBitMCProcess) {
		fMCProcessUls = new TH2F("g0_MCProcessUls",
		                         ";MC Process ID;double-gap topology",
		                         105, 1., 106., 12, 0., 12.);
		TAxis* axis = fMCProcessUls->GetYaxis();
		axis->SetBinLabel(1, "V0");
		axis->SetBinLabel(2, "FMD");
		axis->SetBinLabel(3, "SPD");
		axis->SetBinLabel(4, "TPC");
		axis->SetBinLabel(5, "V0-FMD");
		axis->SetBinLabel(6, "V0-FMD-SPD");
		axis->SetBinLabel(7, "V0-FMD-SPD-TPC");
		axis->SetBinLabel(8, "TPC-SPD");
		axis->SetBinLabel(9, "TPC-SPD-FMD");
		axis->SetBinLabel(10, "TPC-SPD-FMD-V0");
		axis->SetBinLabel(11, "SPD-FMD");
		axis->SetBinLabel(12, "SPD-FMD-V0");
		fHist->Add(fMCProcessUls);

		fMCProcessLs = new TH2F("g1_MCProcessLs",
		                        ";MC Process ID;double-gap topology",
		                        105, 1., 106., 12, 0., 12.);
		axis = fMCProcessLs->GetYaxis();
		axis->SetBinLabel(1, "V0");
		axis->SetBinLabel(2, "FMD");
		axis->SetBinLabel(3, "SPD");
		axis->SetBinLabel(4, "TPC");
		axis->SetBinLabel(5, "V0-FMD");
		axis->SetBinLabel(6, "V0-FMD-SPD");
		axis->SetBinLabel(7, "V0-FMD-SPD-TPC");
		axis->SetBinLabel(8, "TPC-SPD");
		axis->SetBinLabel(9, "TPC-SPD-FMD");
		axis->SetBinLabel(10, "TPC-SPD-FMD-V0");
		axis->SetBinLabel(11, "SPD-FMD");
		axis->SetBinLabel(12, "SPD-FMD-V0");
		fHist->Add(fMCProcessLs);
	}

	//= ALL TRACK INVARIANT MASS =================================================
	if (fAnalysisStatus & AliCDMesonBase::kBitAllTrackMass) {
		fMAllTrackMass = new TH1F("j00_AllTrackInvMass",
		                          ";all track invariant mass (GeV);counts",
		                          95, 0.25, 5.);
		fHist->Add(fMAllTrackMass);
	}
	//= PWA TREE =================================================================
	if (!fAnalysisStatus || fAnalysisStatus & AliCDMesonBase::kBitPWAtree) {
		fPWAtree = new TTree("cd_pwa", "pwa");
		fPWAtree->Branch("theta", &fTheta);
		fPWAtree->Branch("phi", &fPhi);
		fPWAtree->Branch("m", &fMass);
		fPWAtree->Branch("theta", &fTheta);
		fPWAtree->Branch("pX", &fMomentum[0]);
		fPWAtree->Branch("pY", &fMomentum[1]);
		fPWAtree->Branch("pZ", &fMomentum[2]);
		fPWAtree->Branch("gapCondition", &fCurrentGapCondition);
	}

	PostOutputs();
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::UserExec(Option_t *)
{
	//
	// Executed for every event which passed the physics selection
	//

	//printf("Entry: %ld\n", (Long_t)Entry()); // print current event number

	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitStatsFlow)) {
		fhStatsFlow->Fill(AliCDMesonBase::kBinTotalInput); // stats flow
	}

	//= INPUT DATA SANITY TESTS ==================================================
	if (!CheckInput()) {
		PostOutputs();
		return;
	}

	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitStatsFlow)) {
		fhStatsFlow->Fill(AliCDMesonBase::kBinGoodInput); // stats flow
	}

	//= ANALYZE ONLY PHOJET CD EVENTS ============================================
	if (fAnalysisStatus & AliCDMesonBase::kBitCDonly) {
		if (fMCEvent) {
			AliGenEventHeader* header = fMCEvent->GenEventHeader();
			if (header) {
				// Phojet
				if (TString(header->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
					fMCprocess = ((AliGenDPMjetEventHeader*)header)->ProcessType();
					if (fMCprocess != 4) { // no central diffraction
						PostOutputs();
						return;
					}
					else{ // CD found
						fhStatsFlow->Fill(AliCDMesonBase::kBinCDonlyEvents); // statsflow
					}
				}
			}
		}
	}


	//= V0-AND/OR ================================================================
	if (!fDoAOD) {
		if (AliCDMesonUtils::V0AND(fESDEvent))
			fhStatsFlow->Fill(AliCDMesonBase::kBinV0OR);
		if (AliCDMesonUtils::V0OR(fESDEvent))
			fhStatsFlow->Fill(AliCDMesonBase::kBinV0AND);
	}

	//= EMPTY EVENT STUDY ========================================================
	// only implemented on ESDs!
	if (fAnalysisStatus & AliCDMesonBase::kBitEEStudy) {
		DoEmptyEventStudy();
		PostOutputs();
		return;
	}

	//= MC TRUTH =================================================================
	// for now only implemented on ESDs
	Int_t nMCprimaries = 0;
	DetermineMCprocessType();
	if (fAnalysisStatus & AliCDMesonBase::kBitTHnMC) {
		nMCprimaries = DoMCTruth();
	}

	//= EVENT SELECTION ==========================================================
	Int_t kfo =
		((fAnalysisStatus & AliCDMesonBase::kBitFastORStudy) && !fDoAOD) ? 1 : 0;
	Int_t ninnerp=-999, nouterp=-999;
	Bool_t eventIsValid = (fDoAOD) ?
		AliCDMesonUtils::CutEvent(fAODEvent, fhpriv, fPriVtxX, fPriVtxY, fPriVtxZ,
		                          fhpriVtxPos, fhpriVtxDist) :
		AliCDMesonUtils::CutEvent(fESDEvent, fhspd, fhpriv, fhpriVtxPos,
		                          fhpriVtxDist, fhfo, fhFOchans, kfo, ninnerp,
		                          nouterp, fPriVtxX, fPriVtxY, fPriVtxZ);
	if (!eventIsValid) {
		PostOutputs();
		return;
	}

	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitStatsFlow)) {
		fhStatsFlow->Fill(AliCDMesonBase::kBinEventsAfterCuts); // stats flow
	}

	//= PILE UP ==================================================================
	const Bool_t isPileup = (fDoAOD) ?
		fAODEvent->IsPileupFromSPD(2, 0.8, 3., 2., 5.) :
		fESDEvent->IsPileupFromSPD(2, 0.8, 3., 2., 5.);
	// using only 2 instead of three contributors

	if (isPileup) {
		PostOutputs();
		return;
	}

	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitStatsFlow)) {
		fhStatsFlow->Fill(AliCDMesonBase::kBinEventsWithOutPileUp); // stats flow
	}

	//= VZERO TRIGGER STUDY ======================================================
	// dtermine the VZERO signal distributions
	if (fVZEROhists) {
		AliCDMesonUtils::DoVZEROStudy(fESDEvent, fVZEROhists, fRun);
	}

	//= GAP ======================================================================
	// determine the complete gap configuration (for all gap-tagging detectors)
	if (!DetermineGap()) {
		PostOutputs();
		return;
	}

	// fill gaprun (determines the gap fraction per run)
	const Double_t x[2] = {(Double_t)fRun, (Double_t)fCurrentGapCondition};
	if (fGapRun) fGapRun->Fill(x);

	//= VERTEX COINCIDENCE AND POSITION ==========================================
	AnalyzeVtx();

	//= TRACK CUTS ===============================================================
	Bool_t doSoft = !fAnalysisStatus ||
		(fAnalysisStatus & AliCDMesonBase::kBitSoftTracks);

	fTracks->ProcessEvent(fAODEvent, fESDEvent, doSoft); // apply cuts
	DoMultiplicityStudy(nMCprimaries); // fill corresponding histograms
	if (fMAllTrackMass &&
	    (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG)) {
		// calculate the invariant mass of all tracks in the event including soft
		// ones if the event is a full double gap event
		fMAllTrackMass->Fill(fTracks->GetInvariantMass());
	}

	// is multiplicity within the desired range of  2 to 3?
	Int_t nch = fTracks->GetTracks();
	Int_t ncombined = fTracks->GetCombinedTracks();
	Bool_t wSoft = ((ncombined >= 2) && (ncombined <=3)); // including soft tracks
	Bool_t woSoft = ((nch >= 2) && (nch <= 3)); // excluding soft tracks

	if (!wSoft && !woSoft) {
		// multiplicity out of range (both with soft and without soft track)
		PostOutputs();
		return;
	}

	//= TRACK PAIRS ==============================================================
	// loop over all track combinations
	for(Int_t ii=0; ii<ncombined; ii++){
		for(Int_t jj=ii+1; jj<ncombined; jj++){
			// assign current tracks
			fTrkPair[0] = fTracks->GetTrack(ii);
			fTrkPair[1] = fTracks->GetTrack(jj);

			// analyze track pairs
			if (wSoft) DoTrackPair(kTRUE);
			if (woSoft && (ii < nch) && (jj < nch)) DoTrackPair(kFALSE);
		}
	}

	//============================================================================
	PostOutputs();
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::PostOutputs()
{
	//
	// PostData
	//

	Int_t iOutputSlot = 1; // dynamic output slot number handling
	PostData(iOutputSlot++, fHist);
	if (!(fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) { // normal operation
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitTHnMother)) {
			// THnSparse for invariant mass distributions
			PostData(iOutputSlot++, fThnMother);
		}
		if (!fAnalysisStatus
		    || ((fAnalysisStatus & AliCDMesonBase::kBitSoftTracks)
		        && (fAnalysisStatus & AliCDMesonBase::kBitTHnMother))) {
			// same including soft tracks
			PostData(iOutputSlot++, fThnMotherSoft);
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitMultStudy)) {
			// multiplicity study
			PostData(iOutputSlot++, fThnMultiplicity);
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitTHnMC)) {
			// MC truth study is active
			PostData(iOutputSlot++, fThnMC);
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitPWAtree)) {
			// PWA tree active
			PostData(iOutputSlot++, fPWAtree);
		}
	}
	else if (fAnalysisStatus & AliCDMesonBase::kBitEEStudy) {
		// empty event study
		PostData(iOutputSlot++, fThnEmptyEvents);
	}
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::FillEtaPhiMaps() const
{
	//
	// Controls which event contributes to which map
	//

	AliVEvent* event = (fDoAOD) ? (AliVEvent*)fAODEvent : (AliVEvent*)fESDEvent;

	if (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG) {
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0fmdspdtpcMap,
		                               fv0fmdspdtpcMapCutted);
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0fmdspdMap,
		                               fv0fmdspdMapCutted);
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0fmdMap, fv0fmdMapCutted);
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0Map, fv0MapCutted);
	}
	else if (fGapInformation[kV0FMDSPD] == AliCDMesonBase::kBinDG) {
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0fmdspdMap,
		                               fv0fmdspdMapCutted);
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0fmdMap, fv0fmdMapCutted);
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0Map, fv0MapCutted);
	}
	else if (fGapInformation[kV0FMD] == AliCDMesonBase::kBinDG) {
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0fmdMap, fv0fmdMapCutted);
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0Map, fv0MapCutted);
	}
	else if (fGapInformation[kV0] == AliCDMesonBase::kBinDG) {
		AliCDMesonUtils::FillEtaPhiMap(event, fTracks, fv0Map, fv0MapCutted);
	}
}


//------------------------------------------------------------------------------
Bool_t AliAnalysisTaskCDMeson::CheckInput()
{
	//
	//general protection
	//
	if (const AliESDInputHandler *esdH =
	    dynamic_cast<AliESDInputHandler*>(fInputHandler)){
		fESDEvent = esdH->GetEvent();
	}
	else if (const AliAODInputHandler *aodH =
	         dynamic_cast<AliAODInputHandler*>(fInputHandler)){
		fAODEvent = aodH->GetEvent();
		fDoAOD = kTRUE;
	}
	fNoGapFraction = (fDoAOD) ? 1. : fNoGapFraction; // process all running on AOD
	fPIDResponse = (AliPIDResponse*)fInputHandler->GetPIDResponse();

	if(!fESDEvent && !fAODEvent){
		printf("AliAnalysisTaskCDMeson No valid event\n");
		return kFALSE;
	}

	if(!fPIDResponse){
		printf("AliAnalysisTaskCDMeson No PIDd\n");
		// PID is fixed to unknown
		//return kFALSE;
	}

	if(fDoAOD && fAODEvent && fabs(fAODEvent->GetMagneticField())<1){
		printf("AliAnalysisTaskCDMeson strange Bfield! %f\n",
		       fAODEvent->GetMagneticField());
		return kFALSE;
	}
	else if((!fDoAOD) && fESDEvent && fabs(fESDEvent->GetMagneticField())<1){
		printf("AliAnalysisTaskCDMeson strange Bfield! %f\n",
		       fESDEvent->GetMagneticField());
		return kFALSE;
	}

	Int_t tmprun = 0;
	if (fDoAOD) {
		tmprun = fAODEvent->GetRunNumber();
	}
	else {
		tmprun = fESDEvent->GetRunNumber();
	}

	if(fRun!=tmprun){
		fRun = tmprun;
		AliCDMesonUtils::SPDLoadGeom(fRun);

		// change PID cuts if necessary
		if ((fRun >= 162933) && (fRun <=  165462)) {
			fPIDmode = 1;
		}
		else {
			fPIDmode = 0;
		}
	}

	// get MC event
	fMCEvent = MCEvent();

	return kTRUE;
}

//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::DoEmptyEventStudy() {
	// Analyses the observables needed for the gap determination in empty events
	// for cross checks same observables are stored for A/C/I(B)-triggers as well
	//


	if (!(fAnalysisStatus & AliCDMesonBase::kBitEEStudy)) {
		// check whether empty event analysis is activated
		return;
	}

	if(!fESDEvent || fDoAOD) { // ensure that we are running on ESDs
		return;
	}

	Int_t eventType = AliCDMesonUtils::GetEventType(fESDEvent);
	TRandom3 rnd(0);
	if ((eventType != AliCDMesonBase::kBinEventUnknown) &&
	    ((rnd.Rndm() > 0.95) ||
	     (AliCDMesonBase::kBinEventE == eventType))) {
		Int_t fmdA = 0, fmdC = 0, spdIA = 0, spdIC = 0, spdOA = 0, spdOC = 0,
			spdTrklSum = 0, spdTrklCentral = 0, spdTrklForwardA = 0,
			spdTrklForwardC = 0;

		Float_t fmdSums[] = { 0., 0., 0., 0., 0. };
		AliCDMesonUtils::GetMultFMD(fESDEvent, fmdA, fmdC, fmdSums);
		// Obtain FMD multiplicity using FMDHitCombinations
		AliCDMesonUtils::GetMultSPD(fESDEvent, spdIA, spdIC, spdOA, spdOC);
		// Obtain SPD multiplicity using FastOR information
		AliCDMesonUtils::GetSPDTrackletMult(fESDEvent, spdTrklSum,
		                                    spdTrklForwardA, spdTrklForwardC,
		                                    spdTrklCentral);
		// obtain SPD multiplicity using tracklets via AliMultiplicity

		AliCDMesonBase::FillThnEmptyEvents(fThnEmptyEvents, eventType, fmdA, fmdC,
		                                   spdIA, spdIC, spdOA, spdOC,
		                                   spdTrklForwardA, spdTrklForwardC,
		                                   (Int_t)fmdSums[0], (Int_t)fmdSums[1],
		                                   (Int_t)fmdSums[2], (Int_t)fmdSums[3],
		                                   (Int_t)fmdSums[4]);
	}
	// the following code is only executed for empty events, this ensures, that
	// all the histograms contain only information on empty events, while the
	// THnSparse contains all information on A/C/I events
	if (AliCDMesonBase::kBinEventE == eventType) {
		TH1F* fmdSums[] = { fFMDsum1I, fFMDsum2I, fFMDsum2O, fFMDsum3I,
		                    fFMDsum3O };
		AliCDMesonUtils::GetGapConfig(fESDEvent, fHitMapSPDinner, fHitMapSPDouter,
		                              fHitMapSPDtrklt, fHitMapFMDa, fHitMapFMDc,
		                              (TH1**)fmdSums, fTPCGapDCAaSide,
		                              fTPCGapDCAcSide);
		// fill hit maps
	}
	return;
}


//------------------------------------------------------------------------------
Bool_t AliAnalysisTaskCDMeson::DetermineGap()
{
	// determines the gap configuration for all gap tagging detectors based on the
	// data set which is available
	//

	fCurrentGapCondition = 0x0; // initialize gap condition

	if (fDoAOD) {
		AliAODHandler* aodHandler =
			(AliAODHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
		TTree *aodTree = aodHandler->GetTree();
		if (aodTree) {
			aodTree->SetBranchAddress("gapCondition", &fCurrentGapCondition);
			aodTree->GetEvent(Entry()); // seems to be needed! (loads current event)
		}
		if (!fCurrentGapCondition) {
			fCurrentGapCondition = 0xfffe;
			puts("AliAnalysisTaskCDMeson - error while gap condition determination using AODs\n");
			return kFALSE;
		}
		// DEBUGGING STUFF
		//aodTree->ls();
		//aodTree->MakeCode();
		//aodTree->MakeClass();
		//aodTree->GetEvent(Entry());
		//aodTree->Show();
		//TBranch* branch = aodTree->GetBranch("gapCondition");
		//printf("AliAnalysisTaskCDMeson - branch=x%x\n",branch);
		//if (!fCurrentGapCondition) {
		//branch->SetAddress(&fCurrentGapCondition);
			//}
		//Int_t entry = Entry();
		//printf("Entry()=%d\n", entry);
		//branch->GetEvent(entry);
		//printf("AliAnalysisTaskCDMeson - gapcondition=0x%x\n", fCurrentGapCondition);
	}
	else {
		if (fAnalysisStatus & AliCDMesonBase::kBitReadPreprocessedGap) {
			// retrieve preprocessed gap information
			AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
			TTree* esdTree = am->GetInputEventHandler()->GetTree(); // get ESD tree
			if (esdTree) {
				esdTree->SetBranchAddress("gapCondition", &fCurrentGapCondition);
				esdTree->GetEvent(Entry()); // seems to be needed! (loads current event)
			}
		}
		else {
			// determine gaps from ESDs
			TH1F* fmdSums[] = {fFMDsum1I, fFMDsum2I, fFMDsum2O, fFMDsum3I, fFMDsum3O};
			fCurrentGapCondition = AliCDMesonUtils::GetGapConfig(fESDEvent,
			                                                     fHitMapSPDinner,
			                                                     fHitMapSPDouter,
			                                                     fHitMapSPDtrklt,
			                                                     fHitMapFMDa,
			                                                     fHitMapFMDc,
			                                                     (TH1**)fmdSums,
			                                                     fTPCGapDCAaSide,
			                                                     fTPCGapDCAcSide);
		}
		if (!fCurrentGapCondition) {
			fCurrentGapCondition = 0xfffe;
			puts("AliAnalysisTaskCDMeson - error while gap condition determination using ESDs\n");
			return kFALSE;
		}
	}
	//printf("AliAnalysisTaskCDMeson - Event: %ld, gapCondition: 0x%x\n", (Long_t)Entry(),
	//       fCurrentGapCondition);


	// disentagle the contributions to the gap conditions of different "tightness"
	fGapInformation[kV0] =
		AliCDMesonBase::GetGapBin("V0", fCurrentGapCondition, kFALSE);
	fGapInformation[kV0FMD] =
		AliCDMesonBase::GetGapBin("V0FMD", fCurrentGapCondition, kFALSE);
	fGapInformation[kV0FMDSPD] =
		AliCDMesonBase::GetGapBin("V0FMDSPD", fCurrentGapCondition, kFALSE);
	fGapInformation[kV0FMDSPDTPC] =
		AliCDMesonBase::GetGapBin("V0FMDSPDTPC", fCurrentGapCondition, kFALSE);
	fGapInformation[kFMD] =
		AliCDMesonBase::GetGapBin("FMD",fCurrentGapCondition, kFALSE);
	fGapInformation[kSPD] =
		AliCDMesonBase::GetGapBin("SPD",fCurrentGapCondition, kFALSE);
	fGapInformation[kTPC] =
		AliCDMesonBase::GetGapBin("TPC",fCurrentGapCondition, kFALSE);
	fGapInformation[kTPCSPD] =
		AliCDMesonBase::GetGapBin("TPCSPD",fCurrentGapCondition, kFALSE);
	fGapInformation[kTPCSPDFMD] =
		AliCDMesonBase::GetGapBin("TPCSPDFMD",fCurrentGapCondition, kFALSE);
	fGapInformation[kTPCSPDFMDV0] =
		AliCDMesonBase::GetGapBin("TPCSPDFMDV0",fCurrentGapCondition, kFALSE);
	fGapInformation[kSPDFMD] =
		AliCDMesonBase::GetGapBin("SPDFMD",fCurrentGapCondition, kFALSE);
	fGapInformation[kSPDFMDV0] =
		AliCDMesonBase::GetGapBin("SPDFMDV0",fCurrentGapCondition, kFALSE);

	fGapInformationWCent[kV0] = AliCDMesonBase::GetGapBin("V0", fCurrentGapCondition);
	fGapInformationWCent[kV0FMD] =
		AliCDMesonBase::GetGapBin("V0FMD", fCurrentGapCondition);
	fGapInformationWCent[kV0FMDSPD] =
		AliCDMesonBase::GetGapBin("V0FMDSPD", fCurrentGapCondition);
	fGapInformationWCent[kV0FMDSPDTPC] =
		AliCDMesonBase::GetGapBin("V0FMDSPDTPC", fCurrentGapCondition);
	fGapInformationWCent[kFMD] =
		AliCDMesonBase::GetGapBin("FMD",fCurrentGapCondition);
	fGapInformationWCent[kSPD] =
		AliCDMesonBase::GetGapBin("SPD",fCurrentGapCondition);
	fGapInformationWCent[kTPC] =
		AliCDMesonBase::GetGapBin("TPC",fCurrentGapCondition);
	fGapInformationWCent[kTPCSPD] =
		AliCDMesonBase::GetGapBin("TPCSPD",fCurrentGapCondition);
	fGapInformationWCent[kTPCSPDFMD] =
		AliCDMesonBase::GetGapBin("TPCSPDFMD",fCurrentGapCondition);
	fGapInformationWCent[kTPCSPDFMDV0] =
		AliCDMesonBase::GetGapBin("TPCSPDFMDV0",fCurrentGapCondition);
	fGapInformationWCent[kSPDFMD] =
		AliCDMesonBase::GetGapBin("SPDFMD",fCurrentGapCondition);
	fGapInformationWCent[kSPDFMDV0] =
		AliCDMesonBase::GetGapBin("SPDFMDV0",fCurrentGapCondition);

	return kTRUE;
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::DoMultiplicityStudy(Int_t nMCprimaries)
{
	// stores the multiplicity distributions for different gap conditions and
	// compares the multiplicity distribution for "normal primary tracks (ITSTPC)"
	// and the multiplicity distribution furthermore taking also soft tracks into
	// account
	//

	// retrieve values from the track object
	Int_t ntrk0 = fTracks->GetTracksBeforeCuts();
	Int_t nch = fTracks->GetTracks();
	Int_t ncombined = fTracks->GetCombinedTracks();
	Int_t nITSpureSA = fTracks->GetITSpureSACount();

	// determine the residual tracks / tracklets
	fResidualTracks = ntrk0 - ncombined - nITSpureSA;
	fResidualTracklets = fTracks->GetRemainingTrackletsCentralBarrel();

	// soft track specific stuff
	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitSoftTracks)) {
		fMultStudy->Fill((Double_t)nch, (Double_t)ncombined);
		if (fGapInformation[kV0] == AliCDMesonBase::kBinDG) {
			fMultStudyV0dg->Fill((Double_t)nch, (Double_t)ncombined);
		}
		if (fGapInformation[kV0FMD] == AliCDMesonBase::kBinDG) {
			fMultStudyV0FMDdg->Fill((Double_t)nch, (Double_t)ncombined);
		}
		if (fGapInformation[kV0FMDSPD] == AliCDMesonBase::kBinDG) {
			fMultStudyV0FMDSPDdg->Fill((Double_t)nch, (Double_t)ncombined);
		}
		if (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG) {
			fMultStudyV0FMDSPDTPCdg->Fill((Double_t)nch, (Double_t)ncombined);
		}
	}


	// multiplicity distributions for different gaps
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitMultPerGapHists)) {
		fv0ntrk->Fill(ncombined, fGapInformation[kV0]);
		fv0fmdntrk->Fill(ncombined, fGapInformation[kV0FMD]);
		fv0fmdspdntrk->Fill(ncombined, fGapInformation[kV0FMDSPD]);
		fv0fmdspdtpcntrk->Fill(ncombined, fGapInformation[kV0FMDSPDTPC]);
	}

	// multiplicity removed by cuts for different gaps
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliCDMesonBase::kBitRmMultPerGapHists)) {
		fv0Rmntrk->Fill(ntrk0-ncombined-nITSpureSA, fGapInformation[kV0]);
		fv0fmdRmntrk->Fill(ntrk0-ncombined-nITSpureSA, fGapInformation[kV0FMD]);
		fv0fmdspdRmntrk->Fill(ntrk0-ncombined-nITSpureSA,
		                      fGapInformation[kV0FMDSPD]);
		fv0fmdspdtpcRmntrk->Fill(ntrk0-ncombined-nITSpureSA,
		                         fGapInformation[kV0FMDSPDTPC]);
	}


	// fill maps with eta phi information for QA purposes
	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitEtaPhiMaps) ||
	    (fAnalysisStatus & AliCDMesonBase::kBitEtaPhiMapsWithCuts)) {
		FillEtaPhiMaps();
	}

	// fill stats and general multiplicity histograms
	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitStatsFlow)) {
		fhntrk->Fill(ncombined); // multiplicity distribution

		if (fGapInformationWCent[kV0] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinv0Gap);
		}
		if (fGapInformationWCent[kV0FMD] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinv0fmdGap);
		}
		if (fGapInformationWCent[kV0FMDSPD] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinv0fmdspdGap);
		}
		if (fGapInformationWCent[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinv0fmdspdtpcGap);
		}
		if (fGapInformationWCent[kTPC] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBintpcGap);
		}
		if (fGapInformationWCent[kTPCSPD] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBintpcspdGap);
		}
		if (fGapInformationWCent[kTPCSPDFMD] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBintpcspdfmdGap);
		}
		if (fGapInformationWCent[kTPCSPDFMDV0] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBintpcspdfmdv0Gap);
		}
		if (fGapInformationWCent[kFMD] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinfmdGap);
		}
		if (fGapInformationWCent[kSPD] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinspdGap);
		}
		if (fGapInformationWCent[kSPDFMD] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinspdfmdGap);
		}
		if (fGapInformationWCent[kSPDFMDV0] == AliCDMesonBase::kBinDG) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinspdfmdv0Gap);
		}
		if (nch == 2) fhStatsFlow->Fill(AliCDMesonBase::kBinTwoTrackEvents);
		else if (nch == 3) fhStatsFlow->Fill(AliCDMesonBase::kBinThreeTrackEvents);
	}

	// fill multiplicity response histograms (data vs. MC)
	if (fAnalysisStatus & AliCDMesonBase::kBitMultResponseMC) {
		fMultResponseMC->Fill(nMCprimaries, ncombined);
		if (fGapInformation[kV0] == AliCDMesonBase::kBinDG) {
			fMultResponseV0dgMC->Fill(nMCprimaries, ncombined);
		}
		if (fGapInformation[kV0FMD] == AliCDMesonBase::kBinDG) {
			fMultResponseV0FMDdgMC->Fill(nMCprimaries, ncombined);
		}
		if (fGapInformation[kV0FMDSPD] == AliCDMesonBase::kBinDG) {
			fMultResponseV0FMDSPDdgMC->Fill(nMCprimaries, ncombined);
		}
		if (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG) {
			fMultResponseV0FMDSPDTPCdgMC->Fill(nMCprimaries, ncombined);
		}
	}

	// event cleanliness
	if (fResidualTracks == 0 && fhStatsFlow) {
		fhStatsFlow->Fill(AliCDMesonBase::kBinResidualTracks);
	}
	if (fResidualTracklets == 0 && fhStatsFlow) {
		fhStatsFlow->Fill(AliCDMesonBase::kBinResidualTracklets);
	}


	// multiplicity THnSparse
	if (!fAnalysisStatus
	    || (fAnalysisStatus & AliCDMesonBase::kBitMultStudy)) {
		AliCDMesonBase::FillThnMultiplicity(fThnMultiplicity, nch,
		                                    ncombined-nch, ncombined,
		                                    fGapInformation[kV0],
		                                    fGapInformation[kFMD],
		                                    fGapInformation[kSPD],
		                                    fGapInformation[kTPC],
		                                    fResidualTracks, fResidualTracklets,
		                                    fVtxZ, fVtxDst, fMCprocessType);
	}

	if (fAnalysisStatus & AliCDMesonBase::kBitFastORmultStudy) {
		if (fhspdV0dg && fhspdV0FMDdg && fhspdV0FMDSPDdg && fhspdV0FMDSPDTPCdg &&
		    !fDoAOD) {
			Int_t mult = AliCDMesonUtils::GetFastORmultiplicity(fESDEvent);

			if (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG) {
				fhspdV0FMDSPDTPCdg->Fill(mult);
				fhspdV0FMDSPDdg->Fill(mult);
				fhspdV0FMDdg->Fill(mult);
				fhspdV0dg->Fill(mult);
				fhspdAfterCuts->Fill(mult);
			}
			else if (fGapInformation[kV0FMDSPD] == AliCDMesonBase::kBinDG) {
				fhspdV0FMDSPDdg->Fill(mult);
				fhspdV0FMDdg->Fill(mult);
				fhspdV0dg->Fill(mult);
				fhspdAfterCuts->Fill(mult);
			}
			else if (fGapInformation[kV0FMD] == AliCDMesonBase::kBinDG) {
				fhspdV0FMDdg->Fill(mult);
				fhspdV0dg->Fill(mult);
				fhspdAfterCuts->Fill(mult);
			}
			else if (fGapInformation[kV0] == AliCDMesonBase::kBinDG) {
				fhspdV0dg->Fill(mult);
				fhspdAfterCuts->Fill(mult);
			}
			else {
				fhspdAfterCuts->Fill(mult);
			}
		}
	}
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::DoTrackPair(Bool_t soft /* = kFALSE */)
{
	//
	// process a two track pair
	//

	AliCDMesonUtils::SwapTrack((const AliVTrack**)fTrkPair); // random order

	const Int_t combCh =
		AliCDMesonUtils::GetCombCh((const AliVTrack**)fTrkPair); // uls or ls?

	const Int_t combPID = DoPID(combCh); // retrieve PID information

	// analyze vertex quality
	Int_t vtxInRng = 0;
	Int_t vtxCoincidence = 0;

	vtxInRng = (fabs(fVtxZ) < 4.) ? 1 : 0;
	// 4cm range, determined from detector geometry

	vtxCoincidence = (fVtxDst < fMaxVtxDst) ? 1 : 0;
	// vertex coincidence of the track and SPD vertex

	// initialize kinematic values
	Double_t mass = -999., pt = -999., cts = -999., oa = -999.;
	const TVector3 p1(fTrkPair[0]->Px(), fTrkPair[0]->Py(), fTrkPair[0]->Pz());
	const TVector3 p2(fTrkPair[1]->Px(), fTrkPair[1]->Py(), fTrkPair[1]->Pz());
	const TVector3* momenta[2] = { &p1, &p2 };
	AliCDMesonUtils::GetMassPtCtsOA(combPID, momenta, mass, pt, cts, oa);

	//= THnMother ================================================================
	if (!fAnalysisStatus ||
	    ((fAnalysisStatus & AliCDMesonBase::kBitTHnMother)
	     && (!soft || fAnalysisStatus & AliCDMesonBase::kBitSoftTracks))) {
		TRandom3 rnd(0);
		if (((fGapInformation[kV0] == AliCDMesonBase::kBinDG)
		     && (!(fAnalysisStatus & AliCDMesonBase::kBitReduceGapEvents)
		         || (fGapInformation[kV0FMD] == AliCDMesonBase::kBinDG)
		         || (rnd.Rndm() > (1. - fReducedGapFraction))))
		    || (rnd.Rndm() > (1. - fNoGapFraction)) || fDoAOD
		    || (vtxInRng && vtxCoincidence && (fResidualTracks == 0)
		        && (fResidualTracklets == 0))) {
			// reduce amount of no-gap events by a factor of 10 in order to spare
			// memory and shrink the final data volume

			// do switching needed to process either events including soft or without
			THnSparseD* thn = (soft) ? fThnMotherSoft : fThnMother; // output variable
			Int_t mult = (soft) ? fTracks->GetCombinedTracks() : fTracks->GetTracks();
			// correct multiplicity value
			Int_t residualTracks = (soft) ?
				fResidualTracks : fResidualTracks + fTracks->GetSoftTracks();
			// if the analysis is done without soft tracks, they will be added to the
			// residualTracks

			// fill MC histograms
			if (fMCEvent && (!fAnalysisStatus ||
			                 fAnalysisStatus & AliCDMesonBase::kBitMCProcess)) {
				if (mult == 2) {
					FillMChists(combCh);
				}
			}

			AliCDMesonBase::FillThnMother(thn, (Double_t)mult, (Double_t)combCh,
			                              (Double_t)combPID,
			                              (Double_t)fGapInformation[kV0],
			                              (Double_t)fGapInformation[kFMD],
			                              (Double_t)fGapInformation[kSPD],
			                              (Double_t)fGapInformation[kTPC],
			                              mass, pt, cts, oa, fTrkPair[0]->Pt(),
			                              (residualTracks == 0) ? 0. : 1.,
			                              (Double_t)vtxInRng,
			                              (Double_t)fMCprocessType,
			                              (Double_t)vtxCoincidence,
			                              (fResidualTracklets == 0) ? 0. : 1.);
		}
	}

	//= PWA ======================================================================
	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitPWAtree)) {
		if ((soft && (fTracks->GetCombinedTracks() == 2))
		    || ((fTracks->GetTracks() == 2)
		        && !(!fAnalysisStatus
		             || (fAnalysisStatus & AliCDMesonBase::kBitSoftTracks)))) {
			// TODO refine this condition
			// TODO add flag for 3 track events and soft tracks and the other quality
			// flags
			// PWA is only done for two track events
			AliCDMesonUtils::GetPWAinfo(combPID, (const AliVTrack**)fTrkPair, fTheta,
			                            fPhi, fMass, fMomentum);
			TRandom3 rnd(0);
			if (((fGapInformation[kV0] == AliCDMesonBase::kBinDG)
			     && (!(fAnalysisStatus & AliCDMesonBase::kBitReduceGapEvents)
			         || (fGapInformation[kV0FMD] == AliCDMesonBase::kBinDG)
			         || (rnd.Rndm() > (1. - fReducedGapFraction)))) ||
			    (rnd.Rndm() > (1. - fNoGapFraction))) {
				fPWAtree->Fill();
			}
		}
	}
}


//------------------------------------------------------------------------------
Int_t AliAnalysisTaskCDMeson::DoPID(Int_t combCh)
{
	// determine the PID of the two tracks for full double gap events, store
	// the results in another histogram than for the remaining events
	//

	Int_t combPID = 0;
	if (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG) {
		// double-gap events
		combPID = (combCh == AliCDMesonBase::kBinPM) ? // assignment to uls/ls histo
			AliCDMesonUtils::GetCombPID(fPIDResponse, (const AliVTrack**)fTrkPair,
			                            fPIDmode, fComb2trkPIDulsDG) :
			AliCDMesonUtils::GetCombPID(fPIDResponse, (const AliVTrack**)fTrkPair,
			                            fPIDmode, fComb2trkPIDlsDG);
	}
	else { // non full double-gap events
		combPID = (combCh == AliCDMesonBase::kBinPM) ? // assignment to uls/ls histo
			AliCDMesonUtils::GetCombPID(fPIDResponse, (const AliVTrack**)fTrkPair,
			                            fPIDmode, fComb2trkPIDuls) :
			AliCDMesonUtils::GetCombPID(fPIDResponse, (const AliVTrack**)fTrkPair,
			                            fPIDmode, fComb2trkPIDls);
	}

	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitStatsFlow)) {
		if (combPID == AliCDMesonBase::kBinPion ||
		    combPID == AliCDMesonBase::kBinPionE) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinPionEvents);
		}
		else if (combPID == AliCDMesonBase::kBinKaon ||
				         combPID == AliCDMesonBase::kBinKaonE) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinKaonEvents);
		}
		else if (combPID == AliCDMesonBase::kBinProton ||
		         combPID == AliCDMesonBase::kBinProtonE) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinProtonEvents);
		}
		else if (combPID == AliCDMesonBase::kBinElectron ||
		         combPID == AliCDMesonBase::kBinElectronE) {
			fhStatsFlow->Fill(AliCDMesonBase::kBinElectronEvents);
		}
		else {
			fhStatsFlow->Fill(AliCDMesonBase::kBinUnknownPIDEvents);
		}
	}

	return combPID;
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::DetermineMCprocessType()
{
	//
	// retrieves the MC process type from the AliGenEventHeader and classifies
	// them
	//

	// get MC information
	fMCprocess = -1; //detailed MC sub process information
	fMCprocessType = AliCDMesonBase::kBinND; // ND is default, also for data

	if (fMCEvent) {
		AliGenEventHeader* header = fMCEvent->GenEventHeader();
		if (header) {
			// Pythia6
			if (TString(header->IsA()->GetName()) == "AliGenPythiaEventHeader") {
				fMCprocess = ((AliGenPythiaEventHeader*)header)->ProcessType();
				switch(fMCprocess) {
				case 92:
				case 93:
				case 103:
				case 104: fMCprocessType = AliCDMesonBase::kBinSD; break;
				case 94:
				case 105: fMCprocessType = AliCDMesonBase::kBinDD; break;
				default: fMCprocessType = AliCDMesonBase::kBinND; break;
				}
			}
			// Phojet
			else if (TString(header->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
				fMCprocess = ((AliGenDPMjetEventHeader*)header)->ProcessType();
				switch(fMCprocess) {
				case 5:
				case 6: fMCprocessType = AliCDMesonBase::kBinSD; break;
				case 7: fMCprocessType = AliCDMesonBase::kBinDD; break;
				case 4: fMCprocessType = AliCDMesonBase::kBinCD; break;
				default: fMCprocessType = AliCDMesonBase::kBinND; break;
				}
			}
		}
	}
}

//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::FillMChists(Int_t combCh)
{
	// fill the MC histograms for different gap conditions. Only real two track
	// events are taken into account
	//

	if (fMCprocess > 0) { // MC process
		for (Int_t iGap = kV0; iGap < kMax; ++iGap) {
			if (fGapInformation[iGap] == AliCDMesonBase::kBinDG) {
				if (combCh == AliCDMesonBase::kBinPM) {
					fMCProcessUls->Fill(fMCprocess, iGap);
				}
				else if (combCh == AliCDMesonBase::kBinPPMM) {
					fMCProcessLs->Fill(fMCprocess, iGap);
				}
			}
		}
	}
}


//------------------------------------------------------------------------------
Int_t AliAnalysisTaskCDMeson::DoMCTruth()
{
	//
	// analyses the MC truth and does a gap selection based on the eta values of
	// the primaries
	//

	if (!fMCEvent) return -1;
	AliStack* stack = fMCEvent->Stack();
	if (!stack) return -1;

	//= Multiplicity =============================================================
	// determine number of charged physical primaries on the stack
	Int_t nPhysicalPrimaries = 0;
	Int_t nPrimaries = stack->GetNprimary();
	for (Int_t iTracks = 0; iTracks < nPrimaries; ++iTracks) {
		TParticle* part = stack->Particle(iTracks);
		if (stack->IsPhysicalPrimary(iTracks) && (part->GetPDG()->Charge() != 0.) &&
		    (part->Eta() < 0.9) && (part->Eta() > -0.9)) { // TODO add parameter
			++nPhysicalPrimaries;
		}
	}


	//= Track Gap && Multiplicity in Region Bins =================================
	Int_t gapA = 0;
	Int_t gapAv0 = 0;
	Int_t gapAv0fmd = 0;
	Int_t gapC = 0;
	Int_t gapCv0 = 0;
	Int_t gapCv0fmd = 0;
	Int_t central = 0;
	for (Int_t iTracks = 0; iTracks <  nPrimaries; ++iTracks) {
		TParticle* part = (TParticle*)stack->Particle(iTracks);
		if (part && stack->IsPhysicalPrimary(iTracks) &&
		    (part->GetPDG()->Charge() != 0.)) {
			if (part->Eta() > -0.9 && part->Eta() < 0.9) central++;
			if (part->Eta() > 0.9 && part->Eta() < 5.1) gapA++;
			if (part->Eta() > 2.8 && part->Eta() < 5.1) gapAv0++;
			if (part->Eta() > 1.7 && part->Eta() < 5.1) gapAv0fmd++;
			if (part->Eta() < -0.9 && part->Eta() > -3.7) gapC++;
			if (part->Eta() < -1.9 && part->Eta() > -3.7) gapCv0++;
			if (part->Eta() < -1.9 && part->Eta() > -3.7) gapCv0fmd++;
		}
	}

	if (fMultRegionsMC) {
		// multiplicity distribution separated in A-side, central barrel and C-side
		fMultRegionsMC->Fill(gapA, 0);
		fMultRegionsMC->Fill(central, 1);
		fMultRegionsMC->Fill(gapC, 2);
		fMultRegionsMC->Fill(gapA+gapC, 3);
	}

	// WARNING: gap response determination based on primaries only, method too
	// simple to obtain reliable results; at least decays should be taken into
	// account properly
	if (fGapResponseMCv0Dg) {
		Bool_t gapFromMC = (gapAv0 == 0) && (gapCv0 == 0);
		Bool_t gapFromData = (fGapInformation[kV0] == AliCDMesonBase::kBinDG);

		if (gapFromMC && !gapFromData) {
			fGapResponseMCv0Dg->Fill(0);
		}
		else if (gapFromMC && gapFromData) {
			fGapResponseMCv0Dg->Fill(1);
		}
		else if (!gapFromMC && gapFromData) {
			fGapResponseMCv0Dg->Fill(2);
		}
	}
	if (fGapResponseMCv0fmdDg) {
		Bool_t gapFromMC = (gapAv0fmd == 0) && (gapCv0fmd == 0);
		Bool_t gapFromData = (fGapInformation[kV0FMD] == AliCDMesonBase::kBinDG);

		if (gapFromMC && !gapFromData) {
			fGapResponseMCv0fmdDg->Fill(0);
		}
		else if (gapFromMC && gapFromData) {
			fGapResponseMCv0fmdDg->Fill(1);
		}
		else if (!gapFromMC && gapFromData) {
			fGapResponseMCv0fmdDg->Fill(2);
		}
	}
	if (fGapResponseMCv0fmdspdDg) {
		Bool_t gapFromMC = (gapA == 0) && (gapC == 0);
		Bool_t gapFromData = (fGapInformation[kV0FMDSPD] == AliCDMesonBase::kBinDG);

		if (gapFromMC && !gapFromData) {
			fGapResponseMCv0fmdspdDg->Fill(0);
		}
		else if (gapFromMC && gapFromData) {
			fGapResponseMCv0fmdspdDg->Fill(1);
		}
		else if (!gapFromMC && gapFromData) {
			fGapResponseMCv0fmdspdDg->Fill(2);
		}
	}
	if (fGapResponseMCv0fmdspdtpcDg) {
		Bool_t gapFromMC = (gapA == 0) && (gapC == 0);
		Bool_t gapFromData =
			(fGapInformation[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG);

		if (gapFromMC && !gapFromData) {
			fGapResponseMCv0fmdspdtpcDg->Fill(0);
		}
		else if (gapFromMC && gapFromData) {
			fGapResponseMCv0fmdspdtpcDg->Fill(1);
		}
		else if (!gapFromMC && gapFromData) {
			fGapResponseMCv0fmdspdtpcDg->Fill(2);
		}
	}

	if ((nPhysicalPrimaries < 2 ) || (nPhysicalPrimaries >= 4 ))
		return nPhysicalPrimaries;

	Int_t gapCond = AliCDMesonBase::kBinNG;
	if ((gapA == 0) && (gapC == 0) && (central > 0))
		gapCond = AliCDMesonBase::kBinDG;
	else if ((gapA == 0) && (gapC > 0)) gapCond = AliCDMesonBase::kBinGA;
	else if ((gapA > 0) && (gapC == 0)) gapCond = AliCDMesonBase::kBinGC;
	else if ((gapA > 0) && (gapC > 0)) gapCond = AliCDMesonBase::kBinNG;

	//= Two Track ================================================================
	for (Int_t iTracks = 0; iTracks < stack->GetNprimary(); ++iTracks) {
		TParticle* part1 =0x0;
		part1 = (TParticle*)stack->Particle(iTracks);
		if (part1) {
			if (stack->IsPhysicalPrimary(iTracks) &&
			    (part1->GetPDG()->Charge() != 0.) && (part1->Eta() < 0.9) &&
			    (part1->Eta() > -0.9)) {
				for (Int_t jTracks = iTracks; jTracks < stack->GetNprimary(); ++jTracks) {
					TParticle* part2 = 0x0;
					part2 = (TParticle*)stack->Particle(jTracks);
					if (part2) {
						if (stack->IsPhysicalPrimary(jTracks) &&
						    (part2->GetPDG()->Charge() != 0.) && (part2->Eta() < 0.9) &&
						    (part2->Eta() > -0.9)) {
							TParticle* particles[2] = { part1, part2 };
							DoMCTrackPair(particles, gapCond, nPhysicalPrimaries);
						}
					}
				}
			}
		}
	}
	return nPhysicalPrimaries;
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::DoMCTrackPair(TParticle* particles[],
                                           Int_t gapCond, Int_t multiplicity)
{
	//
	// analyses a two track pair from MC truth
	//

	// swap tracks randomly
	TRandom3 rnd(0);
	if (rnd.Rndm() >= 0.5) {
		TParticle* t = particles[1];
		particles[1] = particles[0];
		particles[0] = t;
	}

	// determine whether the pair is like sign or unlike sign
	const Int_t combCh =
		(particles[0]->GetPDG()->Charge()*particles[1]->GetPDG()->Charge() > 0.) ?
		AliCDMesonBase::kBinPPMM : AliCDMesonBase::kBinPM;

	// determine PID
	Int_t combPID = 0;
	if (gapCond == AliCDMesonBase::kBinDG) {
		combPID = (combCh == AliCDMesonBase::kBinPM) ? // assignment to uls/ls histo
			AliCDMesonUtils::GetCombPID((const TParticle**)particles,
			                            fComb2trkPIDulsDGmc) :
			AliCDMesonUtils::GetCombPID((const TParticle**)particles,
			                            fComb2trkPIDlsDGmc);
	}
	else {
		combPID = (combCh == AliCDMesonBase::kBinPM) ? // assignment to uls/ls histo
			AliCDMesonUtils::GetCombPID((const TParticle**)particles,
			                            fComb2trkPIDulsMC) :
			AliCDMesonUtils::GetCombPID((const TParticle**)particles,
			                            fComb2trkPIDlsMC);
	}

	// vertex quality is always nice for MC truth ;)
	Int_t vtxInRng = 1;
	Int_t vtxCoincidence = 1;

	// initialize kinematic values
	Double_t mass = -999., pt = -999., cts = -999., oa = -999.;
	const TVector3 p1(particles[0]->Px(), particles[0]->Py(), particles[0]->Pz());
	const TVector3 p2(particles[1]->Px(), particles[1]->Py(), particles[1]->Pz());
	const TVector3* momenta[2] = { &p1, &p2 };
	AliCDMesonUtils::GetMassPtCtsOA(combPID, momenta, mass, pt, cts, oa);

	//= THnMother ================================================================
	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitTHnMother)) {
		if ((gapCond != AliCDMesonBase::kBinNG) ||
		    (rnd.Rndm() > (1. - fNoGapFraction)) || fDoAOD) {

			AliCDMesonBase::FillThnMother(fThnMC, multiplicity, (Double_t)combCh,
			                              (Double_t)combPID, (Double_t)gapCond,
			                              (Double_t)gapCond, (Double_t)gapCond,
			                              (Double_t)gapCond, mass, pt, cts, oa,
			                              particles[0]->Pt(), 0., (Double_t)vtxInRng,
			                              (Double_t)fMCprocessType,
			                              (Double_t)vtxCoincidence, 1.);
		}
	}

	//= PWA ======================================================================
	// TODO ? should this be implemented???
}


//--------------------------------------------------------------------------
void AliAnalysisTaskCDMeson::AnalyzeVtx()
{
	// calculates the distance between the vertex obtain from tracks and the
	// vertex obtain from spd tracklets
	// stores the z position of the primary vertex from tracks

	fVtxDst = 0.; // reset distance

	// retrieve the pointers of the current primary vertices
	AliVVertex* trackVtx = (fDoAOD) ?
		(AliVVertex*)fAODEvent->GetPrimaryVertex() :
		(AliVVertex*)fESDEvent->GetPrimaryVertexTracks();
	AliVVertex* spdVtx = (fDoAOD) ?
		(AliVVertex*)fAODEvent->GetPrimaryVertexSPD() :
		(AliVVertex*)fESDEvent->GetPrimaryVertexSPD();

	fVtxZ = trackVtx->GetZ(); // store the vertex z position

	if (fDoAOD && (trackVtx == spdVtx)) { // no primary track vertex in the AOD
		fVtxDst = -5; // set arbitrary distance (counted in underflow bin!)
	}
	else { // do proper calculation of the geometrical distance
		fVtxDst += (trackVtx->GetX()-spdVtx->GetX())
			* (trackVtx->GetX()-spdVtx->GetX());
		fVtxDst += (trackVtx->GetY()-spdVtx->GetY())
			* (trackVtx->GetY()-spdVtx->GetY());
		fVtxDst += (trackVtx->GetZ()-spdVtx->GetZ())
			* (trackVtx->GetZ()-spdVtx->GetZ());
		fVtxDst = TMath::Sqrt(fVtxDst);
	}

	if (!fAnalysisStatus || (fAnalysisStatus & AliCDMesonBase::kBitVtxStudies)) {
		fPriVtxDst->Fill(fVtxDst);
		if (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBase::kBinDG) {
			fPriVtxDstV0FMDSPDTPCdg->Fill(fVtxDst);
			fPriVtxDstV0FMDSPDdg->Fill(fVtxDst);
			fPriVtxDstV0FMDdg->Fill(fVtxDst);
			fPriVtxDstV0dg->Fill(fVtxDst);
		}
		else if (fGapInformation[kV0FMDSPD] == AliCDMesonBase::kBinDG) {
			fPriVtxDstV0FMDSPDdg->Fill(fVtxDst);
			fPriVtxDstV0FMDdg->Fill(fVtxDst);
			fPriVtxDstV0dg->Fill(fVtxDst);
		}
		else if (fGapInformation[kV0FMD] == AliCDMesonBase::kBinDG) {
			fPriVtxDstV0FMDdg->Fill(fVtxDst);
			fPriVtxDstV0dg->Fill(fVtxDst);
		}
		else if (fGapInformation[kV0] == AliCDMesonBase::kBinDG) {
			fPriVtxDstV0dg->Fill(fVtxDst);
		}
	}
}
