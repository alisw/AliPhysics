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

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <THnSparse.h>
#include <TObject.h>
#include <TKey.h>

#include <TRandom3.h>

#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliESDInputHandler.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliPhysicsSelection.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliStack.h"

#include "AliAnalysisTaskPB.h"
#include "AliPBBase.h"
#include "AliPBUtils.h"
#include "AliPBTracks.h"



//------------------------------------------------------------------------------
AliAnalysisTaskPB::AliAnalysisTaskPB(const char* name, Long_t state):
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
	, fPIDCombined1(0x0)
	, fPIDCombined2(0x0)
	, fPIDCombined3(0x0)
	, fPhysicsSelection(0x0)
	, fTracks(new AliPBTracks())
	, fVtxDst(-1.)
	, fVtxZ(-20)
	, fResidualTracks(0)
	, fResidualTrackletsCB(0)
	, fResidualTrackletsFW(0)
	, fMCprocessType(0)
	, fMCprocess(-1)

	, fRun(-999)
	, fEvent(0)
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
  
  , fmulADA(0x0)
  , fmulADC(0x0)
  , ftimeADA(0x0)
  , ftimeADC(0x0)

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
  
  , fCEPEvent(0x0)
{
	//
	// standard constructor (the one which should be used)
	//
	// slot in TaskSE must start from 1
	Int_t iOutputSlot = 1;
	DefineOutput(iOutputSlot++, TList::Class());
	if (!(fAnalysisStatus & AliPBBase::kBitEEStudy)) {
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitTHnMother)) {
			DefineOutput(iOutputSlot++, THnSparse::Class());
		}
		if (!fAnalysisStatus
		    || ((fAnalysisStatus & AliPBBase::kBitSoftTracks)
		        && (fAnalysisStatus & AliPBBase::kBitTHnMother))) {
			DefineOutput(iOutputSlot++, THnSparse::Class());
		}
		if (fAnalysisStatus & AliPBBase::kBitMultStudy) {
			DefineOutput(iOutputSlot++, THnSparse::Class());
		}
		if (fAnalysisStatus & AliPBBase::kBitTHnMC) {
			DefineOutput(iOutputSlot++, THnSparse::Class());
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitPWAtree)) {
			DefineOutput(iOutputSlot++, TTree::Class());
		}
	}
	else { // create empty event THnSparse
		DefineOutput(iOutputSlot++, THnSparse::Class());
	}

	if (fAnalysisStatus & AliPBBase::kBitEEStudy) { // empty event study
		fPhysicsSelection = new AliPhysicsSelection;
		fPhysicsSelection->SetSkipTriggerClassSelection();
		// accept all (A,C,E and I) triggers of a certain class
	}

	for (Int_t iGap = 0; iGap < kMax; ++iGap) {
		fGapInformation[iGap] = 0;
		fGapInformationWCent[iGap] = 0;
	}

	for (Int_t ii = 0; ii<2; ii++) {
    for (Int_t jj = 0; jj <= AliPID::kSPECIES; jj++) {
      fPIDnSigmaTPC[ii][jj] = 0.;
      fPIDnSigmaProbTPC[ii][jj] = 0.;
      fPIDnSigmaTOF[ii][jj] = 0.;
      fPIDnSigmaProbTOF[ii][jj] = 0.;
      fPIDBayesWP[ii][jj]  = 0.;
      fPIDBayesNP[ii][jj]  = 0.;
    }
    for (Int_t jj = 0; jj < kTrackInfo; jj++) fTrackInfo[ii][jj] = 0.;
	}

	for (Int_t ii = 0; ii < 3; ii++) {
		fMomentum[ii] = 0.;
	}

	// reset track pair pointers
	fTrkPair[0] = 0x0;
	fTrkPair[1] = 0x0;
}


//------------------------------------------------------------------------------
AliAnalysisTaskPB::AliAnalysisTaskPB():
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
	, fPIDCombined1(0x0)
	, fPIDCombined2(0x0)
	, fPIDCombined3(0x0)
	, fPhysicsSelection(0x0)
	, fTracks(0x0)
	, fVtxDst(-1.)
	, fVtxZ(-20)
	, fResidualTracks(0)
	, fResidualTrackletsCB(0)
	, fResidualTrackletsFW(0)
	, fMCprocessType(0)
	, fMCprocess(-1)

	, fRun(-999)
	, fEvent(0)
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
  
  , fmulADA(0x0)
  , fmulADC(0x0)
  , ftimeADA(0x0)
  , ftimeADC(0x0)

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

  , fCEPEvent(0x0)
{
	//
	// default constructor (should not be used)
	//

	for (Int_t iGap = 0; iGap < kMax; ++iGap) {
		fGapInformation[iGap] = 0;
		fGapInformationWCent[iGap] = 0;
	}

	for (Int_t ii = 0; ii<2; ii++) {
    for (Int_t jj = 0; jj <= AliPID::kSPECIES; jj++) {
      fPIDnSigmaTPC[ii][jj] = 0.;
      fPIDnSigmaProbTPC[ii][jj] = 0.;
      fPIDnSigmaTOF[ii][jj] = 0.;
      fPIDnSigmaProbTOF[ii][jj] = 0.;
      fPIDBayesWP[ii][jj]  = 0.;
      fPIDBayesNP[ii][jj]  = 0.;
    }
    for (Int_t jj = 0; jj < kTrackInfo; jj++) fTrackInfo[ii][jj] = 0.;
	}

	for (Int_t ii = 0; ii < 3; ii++) {
		fMomentum[ii] = 0.;
	}

	fTrkPair[0] = 0x0;
	fTrkPair[1] = 0x0;
}


//------------------------------------------------------------------------------
AliAnalysisTaskPB::~AliAnalysisTaskPB()
{
	//
	//Destructor
	//
	//delete fESDEvent;
	//delete fESDpid;

	if (!(fAnalysisStatus & AliPBBase::kBitEEStudy)) { // normal operation
		if ((!fAnalysisStatus || (fAnalysisStatus && AliPBBase::kBitTHnMother))
		    && fThnMother
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fThnMother;
			fThnMother = 0x0;
		}
		if ((!fAnalysisStatus
		     || ((fAnalysisStatus && AliPBBase::kBitSoftTracks)
		         && (fAnalysisStatus && AliPBBase::kBitTHnMother)))
		    && fThnMotherSoft
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fThnMotherSoft;
			fThnMotherSoft = 0x0;
		}
		if ((!fAnalysisStatus
		     || (fAnalysisStatus && AliPBBase::kBitMultStudy))
		    && fThnMultiplicity
		    && (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() !=
		        AliAnalysisManager::kProofAnalysis)) {
			delete fThnMultiplicity;
			fThnMultiplicity = 0x0;
		}
		if ((!fAnalysisStatus || (fAnalysisStatus && AliPBBase::kBitTHnMC))
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
void AliAnalysisTaskPB::UserCreateOutputObjects()
{
	//
	//CreateOutputObjects
	//

	//= THnSparse ================================================================
	if (!(fAnalysisStatus & AliPBBase::kBitEEStudy)) { // normal operation
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitTHnMother)) {
			fThnMother = AliPBBase::GetThnMother("PB_Mother");
		}
		if (!fAnalysisStatus
		    || ((fAnalysisStatus & AliPBBase::kBitSoftTracks)
		        && (fAnalysisStatus & AliPBBase::kBitTHnMother))) {
			fThnMotherSoft =  AliPBBase::GetThnMother("PB_MotherSoft");
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitMultStudy)){
			fThnMultiplicity =  AliPBBase::GetThnMultiplicity();
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitTHnMC)) {
			fThnMC = AliPBBase::GetThnMother("PB_MotherMC");
		}
	}
	else { // empty event studies
		fThnEmptyEvents = AliPBBase::GetThnEmptyEvents();
	}

	//= TList for Histograms =====================================================
	fHist = new TList;
	fHist->SetOwner(); // ensures that the histograms are all deleted on exit!

	//fHist->Add(fAnalysisStatusString); // WARNING: CRASHES MERGING

	//= GAP RUN =
	if (!(fAnalysisStatus & AliPBBase::kBitEEStudy)) { // normal operation
		// LHC10 data
		// runa1=114649, runb1=117630, runc0=117631, runc1=121526, rund0=121527,
		// rund1=126460, rune0=126461, rune1 = 130930; runf0=130931
		// all checked on elog
		const Int_t runb0 = 114900;
		const Int_t runf1 = 131000;
		
    // LHC11 data
		//const Int_t rune0 = 160676;
		//const Int_t runf1 = 165462;

		// LHC10 + LHC11 data
		//const Int_t runb0 = 114650; // first run of LHC10b
		//const Int_t runf1 = 165462; // last run of LHC11f
		
    // LHC15
		//const Int_t runb0 = 224980;
		//const Int_t runf1 = 244630;

		const Int_t run0 = runb0-1, run1 = runf1 + 1;

		const Int_t nrun = (Int_t)(run1-run0);

		const Int_t bins[2] = {nrun, static_cast<Int_t>(AliPBBase::kBitGapMax - 1.)};
		const Double_t xmin[2] = {(Double_t)run0, 1.};
		const Double_t xmax[2] = {
			(Double_t)run1, (Double_t)AliPBBase::kBitGapMax
		};

		fGapRun = new THnSparseI("GapRun", ";run number;gap condition", 2, bins,
		                         xmin, xmax);
		fHist->Add(fGapRun);
	}

	// Different gap conditions are considered
  //
  // V0:
  // V0FMD:
  // V0FMDSPD:
  // V0FMDSPDTPC:
  
  //= MULTIPLICITY PER GAP CONDITION =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliPBBase::kBitMultPerGapHists)) {
		// x: ntrk; y: V0
		fv0ntrk = new TH2D("b00_v0ntrk", ";number of tracks;gap condition",
		                   80, 0., 80., 4, 1., 5.);
		fHist->Add(fv0ntrk);

		// x: ntrk; y: V0FMD
		fv0fmdntrk = new TH2D("b01_v0fmdntrk", ";number of tracks;gap condition",
		                      80, 0., 80., 4, 1., 5.);
		fHist->Add(fv0fmdntrk);

		//x: ntrk; y: V0FMDSPD
		fv0fmdspdntrk =
			new TH2D("b02_v0fmdspdntrk", ";number of tracks;gap condition",
		         80, 0., 80., 4, 1., 5.);
		fHist->Add(fv0fmdspdntrk);

		//x: ntrk; y: V0FMDSPDTPC
		fv0fmdspdtpcntrk =
			new TH2D("b03_v0fmdspdtpcntrk", ";number of tracks;gap condition",
			         80, 0., 80., 4, 1., 5.);
		fHist->Add(fv0fmdspdtpcntrk);
	}

	//= MULTIPLICITY REMOVED PER GAP CONDITION =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliPBBase::kBitRmMultPerGapHists)) {
		//x: ntrk; y: V0
		fv0Rmntrk = new TH2F("b05_v0Rmntrk",
		                     ";number of removed tracks due to cuts;gap condition",
		                     80, 0., 80., 4, 1., 5.);
		fHist->Add(fv0Rmntrk);

		//x: ntrk; y: V0FMD
		fv0fmdRmntrk =
			new TH2F("b06_v0fmdRmntrk",
			         ";number of removed tracks due to cuts;gap condition",
			         80, 0., 80., 4, 1., 5.);
		fHist->Add(fv0fmdRmntrk);

		//x: ntrk; y: V0FMDSPD
		fv0fmdspdRmntrk =
			new TH2F("b07_v0fmdspdRmntrk",
			         ";number of removed tracks due to cuts;gap condition",
			         80, 0., 80., 4, 1., 5.);
		fHist->Add(fv0fmdspdRmntrk);

		//x: ntrk; y: V0FMDSPDTPC
		fv0fmdspdtpcRmntrk =
			new TH2F("b08_v0fmdspdtpcRmntrk",
			         ";number of removed tracks due to cuts;gap condition",
			         80, 0., 80., 4, 1., 5.);
		fHist->Add(fv0fmdspdtpcRmntrk);
	}

	//= SOFT TRACK INFLUENCE =
	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitSoftTracks)) {
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
	if (fAnalysisStatus & AliPBBase::kBitMultResponseMC) {
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
	if (fAnalysisStatus & AliPBBase::kBitFastORmultStudy) {
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
  if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitStatsFlow)) {
		fhspd = new TH1I("a00_spd", ";fired SPD FastOR chips; counts", 50, 0., 50.);
		fHist->Add(fhspd);

		if (fAnalysisStatus & AliPBBase::kBitFastORStudy) {
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

		fhStatsFlow = AliPBBase::GetHistStatsFlow();
		fHist->Add(fhStatsFlow);
	}


	//= ETA-PHI MAPS FOR TRACKS =
	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitEtaPhiMaps)) {
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
	    (fAnalysisStatus & AliPBBase::kBitEtaPhiMapsWithCuts)) {
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
	    (fAnalysisStatus & AliPBBase::kBitHitMapSPD) ||
	    (fAnalysisStatus & AliPBBase::kBitEEStudy)) {
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
	    (fAnalysisStatus & AliPBBase::kBitHitMapFMD) ||
	    (fAnalysisStatus & AliPBBase::kBitEEStudy)) {
		fHitMapFMDa = new TH2F("d4_HitMapFMDa", ";#eta;#phi",
		                       36, 1.5, 5.1, 40, 0., TMath::TwoPi());
		fHist->Add(fHitMapFMDa);
		fHitMapFMDc = new TH2F("d4_HitMapFMDc", ";#eta;#phi",
		                       20, -3.5, -1.5, 40, 0., TMath::TwoPi());
		fHist->Add(fHitMapFMDc);
	}

	//= FMD SUMMATION =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliPBBase::kBitFMDsum) ||
	    (fAnalysisStatus & AliPBBase::kBitEEStudy)) {
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
  
  //= AD study
	if (!fAnalysisStatus ||
    (fAnalysisStatus & AliPBBase::kBitADStudies)) {
    fmulADA = new TH1F("k0_ADAmul",";ADA hit multiplicity;counts",1000,0,1000);
		fHist->Add(fmulADA);
    fmulADC = new TH1F("k0_ADCmul",";ADC hit multiplicity;counts",1000,0,1000);
		fHist->Add(fmulADC);
    ftimeADA = new TH1F("k0_ADActime",";ADA time;counts",1000,0,100);
		fHist->Add(ftimeADA);
    ftimeADC = new TH1F("k0_ADCtime",";ADC time;counts",1000,0,100);
		fHist->Add(ftimeADC);
	}

	//= VERTEX STUDIES =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliPBBase::kBitVtxStudies)) {
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
	    (fAnalysisStatus & AliPBBase::kBitTPCGapStudy)) {
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
	    (fAnalysisStatus & AliPBBase::kBitVZEROStudy)) {
		if (!fDoAOD) { // not yet possible for AODs
			fVZEROhists = AliPBBase::GetHistVZEROStudies(fHist);
		}
	}

	//= PID MATRIX FOR DATA AND MC =
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliPBBase::kBitPIDStudy)) {
		fComb2trkPIDuls = AliPBBase::GetHistPIDStudies("f0_Comb2trkPIDuls");
		fHist->Add(fComb2trkPIDuls);

		fComb2trkPIDls = AliPBBase::GetHistPIDStudies("f1_Comb2trkPIDls");
		fHist->Add(fComb2trkPIDls);

		fComb2trkPIDulsDG =
			AliPBBase::GetHistPIDStudies("f2_Comb2trkPIDulsDG");
		fHist->Add(fComb2trkPIDulsDG);

		fComb2trkPIDlsDG = AliPBBase::GetHistPIDStudies("f3_Comb2trkPIDlsDG");
		fHist->Add(fComb2trkPIDlsDG);

		if (fAnalysisStatus & AliPBBase::kBitTHnMC) {
			fComb2trkPIDulsMC =
				AliPBBase::GetHistPIDStudies("f4_Comb2trkPIDulsMC");
			fHist->Add(fComb2trkPIDulsMC);

			fComb2trkPIDlsMC = AliPBBase::GetHistPIDStudies("f5_Comb2trkPIDlsMC");
			fHist->Add(fComb2trkPIDlsMC);

			fComb2trkPIDulsDGmc =
				AliPBBase::GetHistPIDStudies("f6_Comb2trkPIDulsDGmc");
			fHist->Add(fComb2trkPIDulsDGmc);

			fComb2trkPIDlsDGmc =
				AliPBBase::GetHistPIDStudies("f7_Comb2trkPIDlsDGmc");
			fHist->Add(fComb2trkPIDlsDGmc);
		}
	}

	//= MC PROCESS IDs =
	if (!fAnalysisStatus || fAnalysisStatus & AliPBBase::kBitMCProcess) {
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
	if (fAnalysisStatus & AliPBBase::kBitAllTrackMass) {
		fMAllTrackMass = new TH1F("j00_AllTrackInvMass",
		                          ";all track invariant mass (GeV);counts",
		                          95, 0.25, 5.);
		fHist->Add(fMAllTrackMass);
	}
	//= PWA TREE =================================================================
	if (!fAnalysisStatus || fAnalysisStatus & AliPBBase::kBitPWAtree) {
		fPWAtree = new TTree("cd_pwa", "pwa");

		for (Int_t ii=0; ii<2; ii++) {
      for (Int_t kk=0; kk<=AliPID::kSPECIES; kk++) {
        fPWAtree->Branch(Form("pid_nSigmaTPC_%i_%i",ii,kk),&fPIDnSigmaTPC[ii][kk]);
        fPWAtree->Branch(Form("pid_nSigmaProbTPC_%i_%i",ii,kk),&fPIDnSigmaProbTPC[ii][kk]);
        fPWAtree->Branch(Form("pid_nSigmaTOF_%i_%i",ii,kk),&fPIDnSigmaTOF[ii][kk]);
        fPWAtree->Branch(Form("pid_nSigmaProbTOF_%i_%i",ii,kk),&fPIDnSigmaProbTOF[ii][kk]);
        fPWAtree->Branch(Form("pid_BayesWP_%i_%i",ii,kk),&fPIDBayesWP[ii][kk]);
        fPWAtree->Branch(Form("pid_BayesNP_%i_%i",ii,kk),&fPIDBayesNP[ii][kk]);
      }
      for (Int_t kk=0; kk<kTrackInfo; kk++) {
		    fPWAtree->Branch(Form("trackInfo_%i_%i",ii,kk),&fTrackInfo[ii][kk]);
      }
    }

		fPWAtree->Branch("theta", &fTheta);
		fPWAtree->Branch("phi", &fPhi);
		fPWAtree->Branch("m", &fMass);
		fPWAtree->Branch("theta", &fTheta);
		fPWAtree->Branch("pX", &fMomentum[0]);
		fPWAtree->Branch("pY", &fMomentum[1]);
		fPWAtree->Branch("pZ", &fMomentum[2]);
		fPWAtree->Branch("gapCondition", &fCurrentGapCondition);
	}

	//= FULL TREE ================================================================
	if (!fAnalysisStatus || fAnalysisStatus & AliPBBase::kBitFulltree) {
		fPWAtree = new TTree("cd_full", "full");

		for (Int_t ii=0; ii<2; ii++) {
      for (Int_t kk=0; kk<=AliPID::kSPECIES; kk++) {
        fPWAtree->Branch(Form("pid_nSigmaTPC_%i_%i",ii,kk),     &fPIDnSigmaTPC[ii][kk]);
        fPWAtree->Branch(Form("pid_nSigmaProbTPC_%i_%i",ii,kk), &fPIDnSigmaProbTPC[ii][kk]);
        fPWAtree->Branch(Form("pid_nSigmaTOF_%i_%i",ii,kk),     &fPIDnSigmaTOF[ii][kk]);
        fPWAtree->Branch(Form("pid_nSigmaProbTOF_%i_%i",ii,kk), &fPIDnSigmaProbTOF[ii][kk]);
        fPWAtree->Branch(Form("pid_BayesWP_%i_%i",ii,kk),    &fPIDBayesWP[ii][kk]);
        fPWAtree->Branch(Form("pid_BayesNP_%i_%i",ii,kk),    &fPIDBayesNP[ii][kk]);
      }
      for (Int_t kk=0; kk<kTrackInfo; kk++) {
		    fPWAtree->Branch(Form("trackInfo_%i_%i",ii,kk),&fTrackInfo[ii][kk]);
      }
    }

		fPWAtree->Branch("theta",&fTheta);
		fPWAtree->Branch("phi",  &fPhi);
		fPWAtree->Branch("m",    &fMass);
		fPWAtree->Branch("pX",   &fMomentum[0]);
		fPWAtree->Branch("pY",   &fMomentum[1]);
		fPWAtree->Branch("pZ",   &fMomentum[2]);
		fPWAtree->Branch("gapCondition", &fCurrentGapCondition);
	}

	//PID Combined
	fPIDCombined1 = new AliPIDCombined;
  fPIDCombined1->SetSelectedSpecies(AliPID::kSPECIES);  //This is default
	fPIDCombined1->SetEnablePriors(kTRUE);
	fPIDCombined1->SetDefaultTPCPriors();                 //Need more update..
  
	fPIDCombined2 = new AliPIDCombined;
  fPIDCombined2->SetSelectedSpecies(AliPID::kSPECIES);  //This is default
	fPIDCombined2->SetEnablePriors(kFALSE);               // priors are set to 1
  
  // set CEP specific priors
  //fPIDCombined3 = new AliPIDCombined;
  //fPIDCombined3->SetSelectedSpecies(AliPID::kSPECIES);  //This is default
	//fPIDCombined3->SetEnablePriors(kTRUE);               // priors are set to 1
  
  //TString fnameMyPriors = TString("/home/pbuehler/physics/projects/alice/CEP/working/forpass4/res/20160501/MergedPriors.root");
  //TH1F *priordistr[AliPID::kSPECIES];
  //GetMyPriors(fnameMyPriors,priordistr);
  //for (Int_t ii=0; ii<AliPID::kSPECIES; ii++)
  //  fPIDCombined3->SetPriorDistribution((AliPID::EParticleType)ii,priordistr[ii]);
  
  //UInt_t Maskin =
  UInt_t Maskin = AliPIDResponse::kDetTPC;

  //UInt_t Maskin =
  //  AliPIDResponse::kDetTPC |
  //  AliPIDResponse::kDetTOF;

  //UInt_t Maskin =
  //AliPIDResponse::kDetTPC |
  //  AliPIDResponse::kDetTOF |
  //  AliPIDResponse::kDetITS;

  //UInt_t Maskin =
  //AliPIDResponse::kDetTPC |
  //  AliPIDResponse::kDetTOF |
  //  AliPIDResponse::kDetITS |
  //  AliPIDResponse::kDetTRD;

  fPIDCombined1->SetDetectorMask(Maskin);
  fPIDCombined2->SetDetectorMask(Maskin);
  //fPIDCombined3->SetDetectorMask(Maskin);
  
  fCEPEvent = new CEPEventBuffer();

  Int_t split = 2;                                  // branches are splitted
  Int_t bsize = 16000; 
  fPWAtree->Branch("CEPEvents",&fCEPEvent, bsize, split);
 
	PostOutputs();
}


//------------------------------------------------------------------------------
void AliAnalysisTaskPB::UserExec(Option_t *)
{
	//
	// Executed for every event which passed the physics selection
	//
	// printf("Entry: %ld\n", (Long_t)Entry()); // print current event number

	// increment event number
  fEvent++;
  //printf("\n\nNext Event ***********\n");
  printf("Event number: %i\n",fEvent);
  
  // 
  if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitStatsFlow)) {
		fhStatsFlow->Fill(AliPBBase::kBinTotalInput); // stats flow
	}

	//= INPUT DATA SANITY TESTS ==================================================
	//printf("Checking Input ...\n");
  if (!CheckInput()) {
		PostOutputs();
		return;
	}

	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitStatsFlow)) {
		fhStatsFlow->Fill(AliPBBase::kBinGoodInput); // stats flow
	}

	//= ANALYZE ONLY PHOJET CD EVENTS ============================================
	if (fAnalysisStatus & AliPBBase::kBitCDonly) {
		if (fMCEvent) {
	    //printf("Checking MC event header ...\n");
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
						fhStatsFlow->Fill(AliPBBase::kBinCDonlyEvents); // statsflow
					}
				}
			}
		}
	}


	//= V0-AND/OR ================================================================
	if (!fDoAOD) {
	  // check trigger status
  	if (AliPBUtils::V0AND(fESDEvent)) fhStatsFlow->Fill(AliPBBase::kBinV0OR);
		if (AliPBUtils::V0OR(fESDEvent))  fhStatsFlow->Fill(AliPBBase::kBinV0AND);
	}

	//= EMPTY EVENT STUDY ========================================================
	// only implemented on ESDs!
	if (fAnalysisStatus & AliPBBase::kBitEEStudy) {
		DoEmptyEventStudy();
		PostOutputs();
		return;
	}

	//= MC TRUTH =================================================================
	// for now only implemented on ESDs
	Int_t nMCprimaries = 0;
  //printf("\nDetermining MC process type ...\n");
	DetermineMCprocessType();
	if (fAnalysisStatus & AliPBBase::kBitTHnMC) {
		nMCprimaries = DoMCTruth();
	}

	//= EVENT SELECTION ==========================================================
	// applied cuts
  //  . vertex exists
  //  . vertex_z < 10 cm
  //printf("\nDoing event selection ...\n");
  Int_t kfo =
		((fAnalysisStatus & AliPBBase::kBitFastORStudy) && !fDoAOD) ? 1 : 0;
	Int_t ninnerp=-999, nouterp=-999;
	Bool_t eventIsValid = (fDoAOD) ?
		AliPBUtils::CutEvent(fAODEvent, fhpriv, fPriVtxX, fPriVtxY, fPriVtxZ,
      fhpriVtxPos, fhpriVtxDist) :
		AliPBUtils::CutEvent(fESDEvent, fhspd, fhpriv, fhpriVtxPos,
      fhpriVtxDist, fhfo, fhFOchans, kfo, ninnerp,
      nouterp, fPriVtxX, fPriVtxY, fPriVtxZ);
	
  //printf("This is event is %i valid\n",eventIsValid);
  if (!eventIsValid) {
		//PostOutputs();
		//return;
	}

	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitStatsFlow)) {
		fhStatsFlow->Fill(AliPBBase::kBinEventsAfterCuts); // stats flow
	}

	//= PILE UP ==================================================================
  // using only 2 instead of three contributors
  //printf("\nChecking for pileup ...\n");
	const Bool_t isPileup = (fDoAOD) ?
	  fAODEvent->IsPileupFromSPD
	  (
	  	2,		// minContributors, default = 3
	  	0.8,	// minZdist, default = 0.8
	  	3., 	// nSigmaZdist, default = 3.
	  	2., 	// nSigmaDiamXY, default = 2.
	  	5.		// nSigmaDiamZ, default = 5.
	  ) :
	  fESDEvent->IsPileupFromSPD
	  (
	  	2,		// minContributors, default = 3
	  	0.8,	// minZdist, default = 0.8
	  	3., 	// nSigmaZdist, default = 3.
	  	2., 	// nSigmaDiamXY, default = 2.
	  	5.		// nSigmaDiamZ, default = 5.
	  );
	//printf("This is a pileup %i event\n",isPileup);
  
  if (isPileup) {
		//PostOutputs();
		//return;
	}

	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitStatsFlow)) {
		fhStatsFlow->Fill(AliPBBase::kBinEventsWithOutPileUp); // stats flow
	}

	//= VZERO TRIGGER STUDY ======================================================
	// dtermine the VZERO signal distributions
	if (fVZEROhists) {
		AliPBUtils::DoVZEROStudy(fESDEvent, fVZEROhists, fRun);
	}

	//= GAP ======================================================================
	// determine the complete gap configuration (for all gap-tagging detectors)
	Bool_t gapcond = DetermineGap();
  
  //if (!gapcond) {
	//	PostOutputs();
	//	return;
	//}

	// fill gaprun (determines the gap fraction per run)
	const Double_t x[2] = {(Double_t)fRun, (Double_t)fCurrentGapCondition};
	if (fGapRun) fGapRun->Fill(x);

	//= VERTEX COINCIDENCE AND POSITION ==========================================
	AnalyzeVtx();

	//= TRACK CUTS ===============================================================
	// counts tracks after cuts were applied
  // and fills buffer with selected cuts
  // TObjArray* fTracks contains normal tracks
  // TObjArray* fSoftTracks contains soft tracks
  Bool_t doSoft = !fAnalysisStatus ||
		(fAnalysisStatus & AliPBBase::kBitSoftTracks);
	fTracks->ProcessEvent(fAODEvent, fESDEvent, doSoft); // apply cuts
	
  DoMultiplicityStudy(nMCprimaries); // fill corresponding histograms
  
	if (fMAllTrackMass &&
	    (fGapInformation[kV0FMDSPDTPC] == AliPBBase::kBinDG)) {
		// calculate the invariant mass of all tracks in the event including soft
		// ones if the event is a full double gap event
		fMAllTrackMass->Fill(fTracks->GetInvariantMass());
	}

	// is multiplicity within the desired range of  2 to 3?
  // more than 2 tracks are needed for background studies
	// GetTracks:             ITS + TPC
	// GetCombinedTracks:     ITS + TPC and ITS only
  Int_t nch       = fTracks->GetTracks();
	Int_t ncombined = fTracks->GetCombinedTracks();

	Bool_t wSoft = ((ncombined >= 2) && (ncombined <= 3) ); // including soft tracks
	Bool_t woSoft = ((nch >= 2) && (nch <= 3));             // excluding soft tracks
	//Bool_t wSoft  = ( ncombined == 2 ); // including soft tracks
	//Bool_t woSoft = ( nch == 2 );       // excluding soft tracks

	//if (!wSoft && !woSoft) {
	//	// multiplicity out of range (both with soft and without soft track)
	//	PostOutputs();
	//	return;
	//}

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
  // save tree for 2 and 4 tracks events
	// is multiplicity within the desired range: 2 or 4?
	// GetTracks:             ITS + TPC
	// GetCombinedTracks:     ITS + TPC and ITS only

	Bool_t fCheck2tracks = (ncombined == 2) ? kTRUE : kFALSE;
	Bool_t fCheck4tracks = (ncombined == 4) ? kTRUE : kFALSE;
	
  Bool_t wSoft2tracks  = (ncombined == 2);  // including soft tracks
  Bool_t wSoft4tracks  = (ncombined == 4);  // including soft tracks
	Bool_t woSoft2tracks =(nch == 2);         // excluding soft tracks
	Bool_t woSoft4tracks =(nch == 4);         // excluding soft tracks

	//if (!wSoft2tracks && !wSoft4tracks && !woSoft2tracks && !woSoft4tracks) {
	//	// multiplicity out of range (both with soft and without soft track)
	//	PostOutputs();
	//	return;
	//}
  
	//============================================================================
  // save information into a CEPEventBuffer
  //  
	Int_t nITSpureSA = fTracks->GetITSpureSACount();

  //printf("\nNumber of tracks: %i %i\n",nch,ncombined-nch);
  //printf("Number of residual tracks: %i %i %i\n",fResidualTracks,fResidualTrackletsCB,fResidualTrackletsFW);
  
  fCEPEvent->Reset();
  
  if (ncombined < 7) {
  
    // prepare MC stack
    AliStack *stack = NULL;
    Int_t nPrimaries = 0;
    if (fMCEvent) {
      stack = fMCEvent->Stack();
      nPrimaries = stack->GetNprimary();
    }

    // set event parameters
    fCEPEvent->SetRunNumber(fRun);
    fCEPEvent->SetEventNumber(fEvent);
    fCEPEvent->SetnumResiduals(fResidualTrackletsCB+fResidualTrackletsFW);
    fCEPEvent->SetGapCondition(fCurrentGapCondition);
  
    // add normal tracks to event buffer
    Double_t mom[3];
    Double_t stat,nsig,probs[AliPID::kSPECIES];
    for (Int_t ii=0; ii<nch; ii++) {
    
      AliVTrack *tmptrk = fTracks->GetTrack(ii);
      
      // create new track
      CEPTrackBuffer *trk = new CEPTrackBuffer();
      trk->SetisSoft(0);
      trk->SetChargeSign((Int_t)tmptrk->Charge());
      tmptrk->GetPxPyPz(mom);
      trk->SetMomentum(TVector3(mom));
      
      // set PID information
      // ... TPC
      stat = fPIDResponse->ComputePIDProbability(AliPIDResponse::kTPC,tmptrk,AliPID::kSPECIES,probs);
      trk->SetPIDTPCStatus(stat);
      trk->SetPIDTPCSignal(tmptrk->GetTPCsignal());
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++) {
        stat = fPIDResponse->NumberOfSigmas(
          AliPIDResponse::kTPC,tmptrk,(AliPID::EParticleType)jj,nsig);
        trk->SetPIDTPCnSigma(jj,nsig);
        trk->SetPIDTPCProbability(jj,probs[jj]);
      }
      
      // ... TOF
      stat = fPIDResponse->ComputePIDProbability(AliPIDResponse::kTOF,tmptrk,AliPID::kSPECIES,probs);
      trk->SetPIDTOFStatus(stat);
      trk->SetPIDTOFSignal(tmptrk->GetTOFsignal());
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++) {
        stat = fPIDResponse->NumberOfSigmas(
          AliPIDResponse::kTOF,tmptrk,(AliPID::EParticleType)jj,nsig);
        trk->SetPIDTOFnSigma(jj,nsig);
        trk->SetPIDTOFProbability(jj,probs[jj]);
      }
      
      // ... Bayes
      stat = fPIDCombined1->ComputeProbabilities(tmptrk, fPIDResponse, probs);
      trk->SetPIDBayesStatus(stat);
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++)
        trk->SetPIDBayesProbability(jj,probs[jj]);
      
      // get MC truth
      Int_t MCind = tmptrk->GetLabel();
      if (nPrimaries > 0 && MCind >= 0) {
        TParticle* part = stack->Particle(MCind);
        //printf("MC mass: %f\n", part->GetMass());
        
        // set MC mass and momentum
        TLorentzVector lv;
        part->Momentum(lv);
        
        trk->SetMCPID(part->GetPdgCode());
        trk->SetMCMass(part->GetMass());
        trk->SetMCMomentum(TVector3(lv.Px(),lv.Py(),lv.Pz()));
      }
      
      fCEPEvent->AddTrack(trk);
    }

    // add soft tracks to event buffer
    for (Int_t ii=0; ii<(ncombined-nch); ii++) {
    
      AliVTrack *tmptrk = fTracks->GetTrack(nch+ii);
      
      // create new track
      CEPTrackBuffer *trk = new CEPTrackBuffer();
      trk->SetisSoft(1);
      trk->SetChargeSign((Int_t)tmptrk->Charge());
      tmptrk->GetPxPyPz(mom);
      trk->SetMomentum(TVector3(mom));
      
      // set PID information
      // ... TPC
      stat = fPIDResponse->ComputePIDProbability(AliPIDResponse::kTPC,tmptrk,AliPID::kSPECIES,probs);
      trk->SetPIDTPCStatus(stat);
      trk->SetPIDTPCSignal(tmptrk->GetTPCsignal());
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++) {
        stat = fPIDResponse->NumberOfSigmas(
          AliPIDResponse::kTPC,tmptrk,(AliPID::EParticleType)jj,nsig);
        trk->SetPIDTPCnSigma(jj,nsig);
        trk->SetPIDTPCProbability(jj,probs[jj]);
      }
      
      // ... TOF
      stat = fPIDResponse->ComputePIDProbability(AliPIDResponse::kTOF,tmptrk,AliPID::kSPECIES,probs);
      trk->SetPIDTOFStatus(stat);
      trk->SetPIDTOFSignal(tmptrk->GetTOFsignal());
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++) {
        stat = fPIDResponse->NumberOfSigmas(
          AliPIDResponse::kTOF,tmptrk,(AliPID::EParticleType)jj,nsig);
        trk->SetPIDTOFnSigma(jj,nsig);
        trk->SetPIDTOFProbability(jj,probs[jj]);
      }
      
      // ... Bayes
      stat = fPIDCombined1->ComputeProbabilities(tmptrk, fPIDResponse, probs);
      trk->SetPIDBayesStatus(stat);
      for (Int_t jj=0; jj<AliPID::kSPECIES; jj++)
        trk->SetPIDBayesProbability(jj,probs[jj]);
      
      // get MC truth
      Int_t MCind = tmptrk->GetLabel();
      if (nPrimaries > 0 && MCind >= 0) {
        TParticle* part = stack->Particle(MCind);
        //printf("MC mass: %f\n", part->GetMass());
        
        // set MC mass and momentum
        TLorentzVector lv;
        part->Momentum(lv);
        
        trk->SetMCPID(part->GetPdgCode());
        trk->SetMCMass(part->GetMass());
        trk->SetMCMomentum(TVector3(lv.Px(),lv.Py(),lv.Pz()));
      }

      fCEPEvent->AddTrack(trk);
    
    }

    // save event  
    //printf("RunNumber: %i\n",fCEPEvent->GetRunNumber());
    //printf("EventNumber: %i\n",fCEPEvent->GetEventNumber());
    //printf("NumberTracks: %i\n",fCEPEvent->GetnumTracks());
    //printf("NumberSoftTracks: %i\n",fCEPEvent->GetnumSoftTracks());
    //printf("GapCondition: %i\n",fCEPEvent->GetGapCondition());
    
    // fill tree
    fPWAtree->Fill();
    
    PostOutputs();
  
  }
  
}


//------------------------------------------------------------------------------
void AliAnalysisTaskPB::PostOutputs()
{
	//
	// PostData
	//

	Int_t iOutputSlot = 1; // dynamic output slot number handling
	PostData(iOutputSlot++, fHist);
	
  if (!(fAnalysisStatus & AliPBBase::kBitEEStudy)) { // normal operation
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitTHnMother)) {
			// THnSparse for invariant mass distributions
			PostData(iOutputSlot++, fThnMother);
		}
		if (!fAnalysisStatus
		    || ((fAnalysisStatus & AliPBBase::kBitSoftTracks)
		        && (fAnalysisStatus & AliPBBase::kBitTHnMother))) {
			// same including soft tracks
			PostData(iOutputSlot++, fThnMotherSoft);
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitMultStudy)) {
			// multiplicity study
			PostData(iOutputSlot++, fThnMultiplicity);
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitTHnMC)) {
			// MC truth study is active
			PostData(iOutputSlot++, fThnMC);
		}
		if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitPWAtree)) {
			// PWA tree active
	    //printf("Writing data ...\n");
			PostData(iOutputSlot++, fPWAtree);
		}
	}
	else if (fAnalysisStatus & AliPBBase::kBitEEStudy) {
		// empty event study
		PostData(iOutputSlot++, fThnEmptyEvents);
	}
}


//------------------------------------------------------------------------------
void AliAnalysisTaskPB::FillEtaPhiMaps() const
{
	//
	// Controls which event contributes to which map
	//

	AliVEvent* event = (fDoAOD) ? (AliVEvent*)fAODEvent : (AliVEvent*)fESDEvent;

	if (fGapInformation[kV0FMDSPDTPC] == AliPBBase::kBinDG) {
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0fmdspdtpcMap,
		                               fv0fmdspdtpcMapCutted);
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0fmdspdMap,
		                               fv0fmdspdMapCutted);
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0fmdMap, fv0fmdMapCutted);
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0Map, fv0MapCutted);
	}
	else if (fGapInformation[kV0FMDSPD] == AliPBBase::kBinDG) {
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0fmdspdMap,
		                               fv0fmdspdMapCutted);
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0fmdMap, fv0fmdMapCutted);
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0Map, fv0MapCutted);
	}
	else if (fGapInformation[kV0FMD] == AliPBBase::kBinDG) {
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0fmdMap, fv0fmdMapCutted);
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0Map, fv0MapCutted);
	}
	else if (fGapInformation[kV0] == AliPBBase::kBinDG) {
		AliPBUtils::FillEtaPhiMap(event, fTracks, fv0Map, fv0MapCutted);
	}
}


//------------------------------------------------------------------------------
Bool_t AliAnalysisTaskPB::CheckInput()
{
	//
	//general protection
	//
	if (const AliESDInputHandler *esdH =
	    dynamic_cast<AliESDInputHandler*>(fInputHandler)) {
	  fESDEvent = (AliESDEvent*)esdH->GetEvent();
	}
	else if (const AliAODInputHandler *aodH =
	         dynamic_cast<AliAODInputHandler*>(fInputHandler)) {
		fAODEvent = aodH->GetEvent();
		fDoAOD = kTRUE;
	}
	else {
	  
	}
	fNoGapFraction = (fDoAOD) ? 1. : fNoGapFraction; // process all running on AOD
	fPIDResponse = (AliPIDResponse*)fInputHandler->GetPIDResponse();

	if(!fESDEvent && !fAODEvent){
		printf("AliAnalysisTaskPB No valid event\n");
		return kFALSE;
	}

	if(!fPIDResponse){
		printf("AliAnalysisTaskPB No PIDd\n");
		// PID is fixed to unknown
		//return kFALSE;
	}

	if(fDoAOD && fAODEvent && fabs(fAODEvent->GetMagneticField())<1){
		printf("AliAnalysisTaskPB strange Bfield! %f\n",
		       fAODEvent->GetMagneticField());
		return kFALSE;
	}
	else if((!fDoAOD) && fESDEvent && fabs(fESDEvent->GetMagneticField())<1){
		printf("AliAnalysisTaskPB strange Bfield! %f\n",
		       fESDEvent->GetMagneticField());
		return kFALSE;
	}

	Int_t tmprun = 0;
	if (fDoAOD && fAODEvent) {
		tmprun = fAODEvent->GetRunNumber();
	}
	else if (fESDEvent) {
		tmprun = fESDEvent->GetRunNumber();
	}

	if(fRun != tmprun){
		fRun = tmprun;
		AliPBUtils::SPDLoadGeom(fRun);

		// change PID cuts if necessary
		if ((fRun >= 114000) && (fRun <=  131000)) {
			fPIDmode = 0;
		}
		else {
			fPIDmode = 1;
		}
	}
  // printf("fRun,fPIDmode: %i %i\n",fRun,fPIDmode);

	// get MC event
	fMCEvent = MCEvent();

	return kTRUE;
}

//------------------------------------------------------------------------------
void AliAnalysisTaskPB::DoEmptyEventStudy() {
  // Analyses the observables needed for the gap determination in empty events
  // for cross checks same observables are stored for A/C/I(B)-triggers as well
  //

	if (!(fAnalysisStatus & AliPBBase::kBitEEStudy)) {
		// check whether empty event analysis is activated
		return;
	}

	if(!fESDEvent || fDoAOD) { // ensure that we are running on ESDs
		return;
	}

	// beam-beam interaction (I), beam-gas (A/C), empty (E) are considered
	Int_t eventType = AliPBUtils::GetEventType(fESDEvent);
  
  TRandom3 rnd(0);
	if ((eventType != AliPBBase::kBinEventUnknown) &&
	    ((rnd.Rndm() > 0.95) ||
	     (AliPBBase::kBinEventE == eventType))) {
		Int_t fmdA = 0, fmdC = 0, spdIA = 0, spdIC = 0, spdOA = 0, spdOC = 0,
			spdTrklSum = 0, spdTrklCentral = 0, spdTrklForwardA = 0,
			spdTrklForwardC = 0;

		Float_t fmdSums[] = { 0., 0., 0., 0., 0. };
		AliPBUtils::GetMultFMD(fESDEvent, fmdA, fmdC, fmdSums);
		// Obtain FMD multiplicity using FMDHitCombinations
		AliPBUtils::GetMultSPD(fESDEvent, spdIA, spdIC, spdOA, spdOC);
		// Obtain SPD multiplicity using FastOR information
		AliPBUtils::GetSPDTrackletMult(fESDEvent, spdTrklSum,
		                                    spdTrklForwardA, spdTrklForwardC,
		                                    spdTrklCentral);
		// obtain SPD multiplicity using tracklets via AliMultiplicity

		AliPBBase::FillThnEmptyEvents(fThnEmptyEvents, eventType, fmdA, fmdC,
		                                   spdIA, spdIC, spdOA, spdOC,
		                                   spdTrklForwardA, spdTrklForwardC,
		                                   (Int_t)fmdSums[0], (Int_t)fmdSums[1],
		                                   (Int_t)fmdSums[2], (Int_t)fmdSums[3],
		                                   (Int_t)fmdSums[4]);
	}
	// the following code is only executed for empty events, this ensures, that
	// all the histograms contain only information on empty events, while the
	// THnSparse contains all information on A/C/I events
	if (AliPBBase::kBinEventE == eventType) {
		TH1F* fmdSums[] = { fFMDsum1I, fFMDsum2I, fFMDsum2O, fFMDsum3I,
		                    fFMDsum3O };
		AliPBUtils::GetGapConfig(fESDEvent,
      fHitMapSPDinner, fHitMapSPDouter, fHitMapSPDtrklt,
      fHitMapFMDa, fHitMapFMDc, (TH1**)fmdSums,
      fmulADA, fmulADC, ftimeADA, ftimeADC,
      fTPCGapDCAaSide, fTPCGapDCAcSide);
		// fill hit maps
	}
	return;
}


//------------------------------------------------------------------------------
Bool_t AliAnalysisTaskPB::DetermineGap()
{
  // determines the gap configuration for all gap tagging detectors based on
  // the data set which is available
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
			puts("AliAnalysisTaskPB - error while gap condition determination using AODs\n");
			return kFALSE;
		}
		// DEBUGGING STUFF
		//aodTree->ls();
		//aodTree->MakeCode();
		//aodTree->MakeClass();
		//aodTree->GetEvent(Entry());
		//aodTree->Show();
		//TBranch* branch = aodTree->GetBranch("gapCondition");
		//printf("AliAnalysisTaskPB - branch=x%x\n",branch);
		//if (!fCurrentGapCondition) {
		//branch->SetAddress(&fCurrentGapCondition);
			//}
		//Int_t entry = Entry();
		//printf("Entry()=%d\n", entry);
		//branch->GetEvent(entry);
		//printf("AliAnalysisTaskPB - gapcondition=0x%x\n", fCurrentGapCondition);
	}
	else {
		if (fAnalysisStatus & AliPBBase::kBitReadPreprocessedGap) {
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
      // 
			TH1F* fmdSums[] = {fFMDsum1I, fFMDsum2I, fFMDsum2O, fFMDsum3I, fFMDsum3O};
		  fCurrentGapCondition = AliPBUtils::GetGapConfig(fESDEvent,
        fHitMapSPDinner, fHitMapSPDouter, fHitMapSPDtrklt,
        fHitMapFMDa, fHitMapFMDc, (TH1**)fmdSums,
        fmulADA, fmulADC, ftimeADA, ftimeADC,
        fTPCGapDCAaSide, fTPCGapDCAcSide);
		}
		if (!fCurrentGapCondition) {
			fCurrentGapCondition = 0xfffe;
			puts("AliAnalysisTaskPB - error while gap condition determination using ESDs\n");
			return kFALSE;
		}
	}
  //printf("gapCondition: %i\n",fCurrentGapCondition);

	// disentagle the contributions to the gap conditions of different "tightness"
	fGapInformation[kV0] =
		AliPBBase::GetGapBin("V0", fCurrentGapCondition, kFALSE);
	fGapInformation[kV0FMD] =
		AliPBBase::GetGapBin("V0FMD", fCurrentGapCondition, kFALSE);
	fGapInformation[kV0AD] =
		AliPBBase::GetGapBin("V0AD", fCurrentGapCondition, kFALSE);
	fGapInformation[kV0FMDAD] =
		AliPBBase::GetGapBin("V0FMDAD", fCurrentGapCondition, kFALSE);
	fGapInformation[kV0FMDSPD] =
		AliPBBase::GetGapBin("V0FMDSPD", fCurrentGapCondition, kFALSE);
	fGapInformation[kV0FMDSPDTPC] =
		AliPBBase::GetGapBin("V0FMDSPDTPC", fCurrentGapCondition, kFALSE);
	fGapInformation[kFMD] =
		AliPBBase::GetGapBin("FMD",fCurrentGapCondition, kFALSE);
	fGapInformation[kAD] =
		AliPBBase::GetGapBin("AD",fCurrentGapCondition, kFALSE);
	fGapInformation[kSPD] =
		AliPBBase::GetGapBin("SPD",fCurrentGapCondition, kFALSE);
	fGapInformation[kTPC] =
		AliPBBase::GetGapBin("TPC",fCurrentGapCondition, kFALSE);
	fGapInformation[kTPCSPD] =
		AliPBBase::GetGapBin("TPCSPD",fCurrentGapCondition, kFALSE);
	fGapInformation[kTPCSPDFMD] =
		AliPBBase::GetGapBin("TPCSPDFMD",fCurrentGapCondition, kFALSE);
	fGapInformation[kTPCSPDFMDV0] =
		AliPBBase::GetGapBin("TPCSPDFMDV0",fCurrentGapCondition, kFALSE);
	fGapInformation[kSPDFMD] =
		AliPBBase::GetGapBin("SPDFMD",fCurrentGapCondition, kFALSE);
	fGapInformation[kSPDFMDV0] =
		AliPBBase::GetGapBin("SPDFMDV0",fCurrentGapCondition, kFALSE);

	fGapInformationWCent[kV0] = 
    AliPBBase::GetGapBin("V0", fCurrentGapCondition);
	fGapInformationWCent[kV0FMD] =
		AliPBBase::GetGapBin("V0FMD", fCurrentGapCondition);
	fGapInformationWCent[kV0AD] =
		AliPBBase::GetGapBin("V0AD", fCurrentGapCondition);
	fGapInformationWCent[kV0FMDAD] =
		AliPBBase::GetGapBin("V0FMDAD", fCurrentGapCondition);
	fGapInformationWCent[kV0FMDSPD] =
		AliPBBase::GetGapBin("V0FMDSPD", fCurrentGapCondition);
	fGapInformationWCent[kV0FMDSPDTPC] =
		AliPBBase::GetGapBin("V0FMDSPDTPC", fCurrentGapCondition);
	fGapInformationWCent[kFMD] =
		AliPBBase::GetGapBin("FMD",fCurrentGapCondition);
	fGapInformationWCent[kAD] =
		AliPBBase::GetGapBin("AD",fCurrentGapCondition);
	fGapInformationWCent[kSPD] =
		AliPBBase::GetGapBin("SPD",fCurrentGapCondition);
	fGapInformationWCent[kTPC] =
		AliPBBase::GetGapBin("TPC",fCurrentGapCondition);
	fGapInformationWCent[kTPCSPD] =
		AliPBBase::GetGapBin("TPCSPD",fCurrentGapCondition);
	fGapInformationWCent[kTPCSPDFMD] =
		AliPBBase::GetGapBin("TPCSPDFMD",fCurrentGapCondition);
	fGapInformationWCent[kTPCSPDFMDV0] =
		AliPBBase::GetGapBin("TPCSPDFMDV0",fCurrentGapCondition);
	fGapInformationWCent[kSPDFMD] =
		AliPBBase::GetGapBin("SPDFMD",fCurrentGapCondition);
	fGapInformationWCent[kSPDFMDV0] =
		AliPBBase::GetGapBin("SPDFMDV0",fCurrentGapCondition);

	return kTRUE;
}


//------------------------------------------------------------------------------
void AliAnalysisTaskPB::DoMultiplicityStudy(Int_t nMCprimaries)
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
	//fResidualTracks = ntrk0 - ncombined - nITSpureSA;
	fResidualTracks = ntrk0 - ncombined;
	fResidualTrackletsCB = fTracks->GetRemainingTrackletsCentralBarrel();
	fResidualTrackletsFW = fTracks->GetRemainingTrackletsForward();

	// soft track specific stuff
	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitSoftTracks)) {
		// nch: 
    // ncombined:
    
    fMultStudy->Fill((Double_t)nch, (Double_t)ncombined);
		if (fGapInformation[kV0] == AliPBBase::kBinDG) {
			fMultStudyV0dg->Fill((Double_t)nch, (Double_t)ncombined);
		}
		if (fGapInformation[kV0FMD] == AliPBBase::kBinDG) {
			fMultStudyV0FMDdg->Fill((Double_t)nch, (Double_t)ncombined);
		}
		if (fGapInformation[kV0FMDSPD] == AliPBBase::kBinDG) {
			fMultStudyV0FMDSPDdg->Fill((Double_t)nch, (Double_t)ncombined);
		}
		if (fGapInformation[kV0FMDSPDTPC] == AliPBBase::kBinDG) {
			fMultStudyV0FMDSPDTPCdg->Fill((Double_t)nch, (Double_t)ncombined);
		}
	}


	// multiplicity distributions for different gaps
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliPBBase::kBitMultPerGapHists)) {
		fv0ntrk->Fill(ncombined, fGapInformation[kV0]);
		fv0fmdntrk->Fill(ncombined, fGapInformation[kV0FMD]);
		fv0fmdspdntrk->Fill(ncombined, fGapInformation[kV0FMDSPD]);
		fv0fmdspdtpcntrk->Fill(ncombined, fGapInformation[kV0FMDSPDTPC]);
	}

	// multiplicity removed by cuts for different gaps
	if (!fAnalysisStatus ||
	    (fAnalysisStatus & AliPBBase::kBitRmMultPerGapHists)) {
		fv0Rmntrk->Fill(ntrk0-ncombined-nITSpureSA, fGapInformation[kV0]);
		fv0fmdRmntrk->Fill(ntrk0-ncombined-nITSpureSA, fGapInformation[kV0FMD]);
		fv0fmdspdRmntrk->Fill(ntrk0-ncombined-nITSpureSA,
		                      fGapInformation[kV0FMDSPD]);
		fv0fmdspdtpcRmntrk->Fill(ntrk0-ncombined-nITSpureSA,
		                         fGapInformation[kV0FMDSPDTPC]);
	}


	// fill maps with eta phi information for QA purposes
	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitEtaPhiMaps) ||
	    (fAnalysisStatus & AliPBBase::kBitEtaPhiMapsWithCuts)) {
		FillEtaPhiMaps();
	}

	// fill stats and general multiplicity histograms
	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitStatsFlow)) {
		fhntrk->Fill(ncombined); // multiplicity distribution

		if (fGapInformationWCent[kV0] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBinv0Gap);
		}
		if (fGapInformationWCent[kV0FMD] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBinv0fmdGap);
		}
		if (fGapInformationWCent[kV0FMDSPD] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBinv0fmdspdGap);
		}
		if (fGapInformationWCent[kV0FMDSPDTPC] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBinv0fmdspdtpcGap);
		}
		if (fGapInformationWCent[kTPC] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBintpcGap);
		}
		if (fGapInformationWCent[kTPCSPD] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBintpcspdGap);
		}
		if (fGapInformationWCent[kTPCSPDFMD] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBintpcspdfmdGap);
		}
		if (fGapInformationWCent[kTPCSPDFMDV0] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBintpcspdfmdv0Gap);
		}
		if (fGapInformationWCent[kFMD] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBinfmdGap);
		}
		if (fGapInformationWCent[kSPD] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBinspdGap);
		}
		if (fGapInformationWCent[kSPDFMD] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBinspdfmdGap);
		}
		if (fGapInformationWCent[kSPDFMDV0] == AliPBBase::kBinDG) {
			fhStatsFlow->Fill(AliPBBase::kBinspdfmdv0Gap);
		}
		if (nch == 2) fhStatsFlow->Fill(AliPBBase::kBinTwoTrackEvents);
		else if (nch == 3) fhStatsFlow->Fill(AliPBBase::kBinThreeTrackEvents);
	}

	// fill multiplicity response histograms (data vs. MC)
	if (fAnalysisStatus & AliPBBase::kBitMultResponseMC) {
		fMultResponseMC->Fill(nMCprimaries, ncombined);
		if (fGapInformation[kV0] == AliPBBase::kBinDG) {
			fMultResponseV0dgMC->Fill(nMCprimaries, ncombined);
		}
		if (fGapInformation[kV0FMD] == AliPBBase::kBinDG) {
			fMultResponseV0FMDdgMC->Fill(nMCprimaries, ncombined);
		}
		if (fGapInformation[kV0FMDSPD] == AliPBBase::kBinDG) {
			fMultResponseV0FMDSPDdgMC->Fill(nMCprimaries, ncombined);
		}
		if (fGapInformation[kV0FMDSPDTPC] == AliPBBase::kBinDG) {
			fMultResponseV0FMDSPDTPCdgMC->Fill(nMCprimaries, ncombined);
		}
	}

	// event cleanliness
	if (fResidualTracks == 0 && fhStatsFlow) {
		fhStatsFlow->Fill(AliPBBase::kBinResidualTracks);
	}
	if (fResidualTrackletsCB == 0 && fhStatsFlow) {
		fhStatsFlow->Fill(AliPBBase::kBinResidualTracklets);
	}


	// multiplicity THnSparse
	if (!fAnalysisStatus
	    || (fAnalysisStatus & AliPBBase::kBitMultStudy)) {
		AliPBBase::FillThnMultiplicity(fThnMultiplicity, nch,
		                                    ncombined-nch, ncombined,
		                                    fGapInformation[kV0],
		                                    fGapInformation[kFMD],
		                                    fGapInformation[kSPD],
		                                    fGapInformation[kTPC],
		                                    fResidualTracks, fResidualTrackletsCB,
		                                    fVtxZ, fVtxDst, fMCprocessType);
	}

	if (fAnalysisStatus & AliPBBase::kBitFastORmultStudy) {
		if (fhspdV0dg && fhspdV0FMDdg && fhspdV0FMDSPDdg && fhspdV0FMDSPDTPCdg &&
		    !fDoAOD) {
			Int_t mult = AliPBUtils::GetFastORmultiplicity(fESDEvent);

			if (fGapInformation[kV0FMDSPDTPC] == AliPBBase::kBinDG) {
				fhspdV0FMDSPDTPCdg->Fill(mult);
				fhspdV0FMDSPDdg->Fill(mult);
				fhspdV0FMDdg->Fill(mult);
				fhspdV0dg->Fill(mult);
				fhspdAfterCuts->Fill(mult);
			}
			else if (fGapInformation[kV0FMDSPD] == AliPBBase::kBinDG) {
				fhspdV0FMDSPDdg->Fill(mult);
				fhspdV0FMDdg->Fill(mult);
				fhspdV0dg->Fill(mult);
				fhspdAfterCuts->Fill(mult);
			}
			else if (fGapInformation[kV0FMD] == AliPBBase::kBinDG) {
				fhspdV0FMDdg->Fill(mult);
				fhspdV0dg->Fill(mult);
				fhspdAfterCuts->Fill(mult);
			}
			else if (fGapInformation[kV0] == AliPBBase::kBinDG) {
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
void AliAnalysisTaskPB::DoTrackPair(Bool_t soft /* = kFALSE */)
{
	//
	// process a two track pair (fTrkPair)
	//

	// printf("<I - AliAnalysisTaskPB::DoTrackPair>\n");
  
  // order of tracks is randomized because the properties of only one track is
  // stored
  AliPBUtils::SwapTrack((const AliVTrack**)fTrkPair); // random order

	// check if track pair is like- or unlike-signe
  const Int_t combCh =
		AliPBUtils::GetCombCh((const AliVTrack**)fTrkPair);

	// retrieve PID information of two track system
  // the possible pids are listed in AliPBBase.h
  const Int_t combPID = DoPID(combCh); 

	// analyze vertex quality
	Int_t vtxInRng = 0;
	Int_t vtxCoincidence = 0;

	vtxInRng = (fabs(fVtxZ) < 4.) ? 1 : 0;
	// 4cm range, determined from detector geometry

	// vertex coincidence of the track and SPD vertex
	vtxCoincidence = (fVtxDst < fMaxVtxDst) ? 1 : 0;

	// initialize kinematic values
	Double_t mass = -999., pt = -999., cts = -999., oa = -999.;
	const TVector3 p1(fTrkPair[0]->Px(), fTrkPair[0]->Py(), fTrkPair[0]->Pz());
	const TVector3 p2(fTrkPair[1]->Px(), fTrkPair[1]->Py(), fTrkPair[1]->Pz());
	const TVector3* momenta[2] = { &p1, &p2 };
	AliPBUtils::GetMassPtCtsOA(combPID, momenta, mass, pt, cts, oa);

	//= THnMother ================================================================
	if (!fAnalysisStatus ||
	    ((fAnalysisStatus & AliPBBase::kBitTHnMother)
	     && (!soft || fAnalysisStatus & AliPBBase::kBitSoftTracks))) {
		TRandom3 rnd(0);
		if (((fGapInformation[kV0] == AliPBBase::kBinDG)
		     && (!(fAnalysisStatus & AliPBBase::kBitReduceGapEvents)
		         || (fGapInformation[kV0FMD] == AliPBBase::kBinDG)
		         || (rnd.Rndm() > (1. - fReducedGapFraction))))
		    || (rnd.Rndm() > (1. - fNoGapFraction)) || fDoAOD
		    || (vtxInRng && vtxCoincidence && (fResidualTracks == 0)
		        && (fResidualTrackletsCB == 0))) {
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
			                 fAnalysisStatus & AliPBBase::kBitMCProcess)) {
				if (mult == 2) {
					FillMChists(combCh);
				}
			}

			AliPBBase::FillThnMother(thn, (Double_t)mult, (Double_t)combCh,
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
			                              (fResidualTrackletsCB == 0) ? 0. : 1.);
		}
	}

	//= PWA ======================================================================
	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitPWAtree)) {
		if ((soft && (fTracks->GetCombinedTracks() == 2))
		    || ((fTracks->GetTracks() == 2)
		        && !(!fAnalysisStatus
		             || (fAnalysisStatus & AliPBBase::kBitSoftTracks)))) {
			// TODO refine this condition
			// TODO add flag for 3 track events and soft tracks and the other quality
			// flags
			// PWA is only done for two track events
			AliPBUtils::GetPWAinfo(combPID, (const AliVTrack**)fTrkPair, fTheta,
        fPhi, fMass, fMomentum);
			
      TRandom3 rnd(0);
			/*
      if (((fGapInformation[kV0] == AliPBBase::kBinDG)
			     && (!(fAnalysisStatus & AliPBBase::kBitReduceGapEvents)
			         || (fGapInformation[kV0FMD] == AliPBBase::kBinDG)
			         || (rnd.Rndm() > (1. - fReducedGapFraction)))) ||
			    (rnd.Rndm() > (1. - fNoGapFraction))) {
      */

			if (
        ( (fGapInformation[kV0] == AliPBBase::kBinDG) &&
          (combCh == AliPBBase::kBinPM) &&
          (
            !(fAnalysisStatus & AliPBBase::kBitReduceGapEvents) ||
            (fGapInformation[kV0FMD] == AliPBBase::kBinDG) ||
            (rnd.Rndm() > (1. - fReducedGapFraction))
          )
        ) ||
        (rnd.Rndm() > (1. - fNoGapFraction)) )
      {
        // fPid is filled in DoPID
        fprintf(stderr,"\n");
        fprintf(stderr," * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
        fprintf(stderr,"Track 1\n");
        fprintf(stderr,"fPID @ filling: %e|%e / %e|%e / %e|%e / %e|%e / %e|%e \n",
          fPIDBayesWP[0][1],fPIDBayesNP[0][1],fPIDBayesWP[0][2],fPIDBayesNP[0][2],
          fPIDBayesWP[0][3],fPIDBayesNP[0][3],fPIDBayesWP[0][4],fPIDBayesNP[0][4],
          fPIDBayesWP[0][5],fPIDBayesNP[0][5]);
        fprintf(stderr,"\n");
        fprintf(stderr,"Track 2\n");
        fprintf(stderr,"fPID @ filling: %e|%e / %e|%e / %e|%e / %e|%e / %e|%e \n",
          fPIDBayesWP[1][1],fPIDBayesNP[1][1],fPIDBayesWP[1][2],fPIDBayesNP[1][2],
          fPIDBayesWP[1][3],fPIDBayesNP[1][3],fPIDBayesWP[1][4],fPIDBayesNP[1][4],
          fPIDBayesWP[1][5],fPIDBayesNP[1][5]);
        fprintf(stderr,"\n");
        fprintf(stderr," * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n");
        fprintf(stderr,"\n");
        
        // fill TrackInfo
        for (Int_t ii=0; ii<2; ii++) {
          fTrackInfo[ii][0] = fTrkPair[ii]->Charge();
          fTrackInfo[ii][1] = fTrkPair[ii]->E();
          fTrackInfo[ii][2] = fTrkPair[ii]->Px();
          fTrackInfo[ii][3] = fTrkPair[ii]->Py();
          fTrackInfo[ii][4] = fTrkPair[ii]->Pz();
          fTrackInfo[ii][5] = fTrkPair[ii]->Eta();
        }
        
        // fill tree
        //fPWAtree->Fill();
        
			}
		}
	}
}

//------------------------------------------------------------------------------
Int_t AliAnalysisTaskPB::DoPID(Int_t combCh)
{
	// determine the PID of the two tracks for full double gap events
  // store the results in another histogram than for the remaining events
	//

  // printf("<I - AliAnalysisTaskPB::DoPID>\n");
	
  Int_t combPID = 0;
	if (fGapInformation[kV0FMDSPDTPC] == AliPBBase::kBinDG) {
		// double-gap events
		combPID = (combCh == AliPBBase::kBinPM) ? // assignment to uls/ls histo
			AliPBUtils::GetCombPID(fPIDResponse, (const AliVTrack**)fTrkPair,
			                            fPIDmode, fComb2trkPIDulsDG) :
			AliPBUtils::GetCombPID(fPIDResponse, (const AliVTrack**)fTrkPair,
			                            fPIDmode, fComb2trkPIDlsDG);
	}
	else { // non full double-gap events
		combPID = (combCh == AliPBBase::kBinPM) ? // assignment to uls/ls histo
			AliPBUtils::GetCombPID(fPIDResponse, (const AliVTrack**)fTrkPair,
			                            fPIDmode, fComb2trkPIDuls) :
			AliPBUtils::GetCombPID(fPIDResponse, (const AliVTrack**)fTrkPair,
			                            fPIDmode, fComb2trkPIDls);
	}

	// get PID for the single tracks
	// kpid[AliPID::kSPECIES]:
  // 0: electron
  // 1: muon
  // 2: pion
  // 3: kaon
  // 4: proton
  Double_t kpid[AliPID::kSPECIES] = {0};
  Float_t ptot;
  for (Int_t ii=0; ii<2; ii++) {            // tracks

    // signals
    fPIDnSigmaTPC[ii][0] = fTrkPair[ii]->GetTPCsignal();
    fPIDnSigmaTOF[ii][0] = fTrkPair[ii]->GetTOFsignal();

    // nSigma values
    Double_t pidTPC,probTPC[AliPID::kSPECIES],pidTOF,probTOF[AliPID::kSPECIES];

    // with TPC
    fPIDnSigmaProbTPC[ii][0] = fPIDResponse->ComputePIDProbability(
      AliPIDResponse::kTPC,fTrkPair[ii],AliPID::kSPECIES,probTPC);
    // with TOF
    fPIDnSigmaProbTOF[ii][0] = fPIDResponse->ComputePIDProbability(
      AliPIDResponse::kTOF,fTrkPair[ii],AliPID::kSPECIES,probTOF);

    for (Int_t jj=0; jj<AliPID::kSPECIES; jj++) {

      Int_t pidStat = fPIDResponse->NumberOfSigmas(
        AliPIDResponse::kTPC,fTrkPair[ii],(AliPID::EParticleType)jj,pidTPC);
      fPIDnSigmaTPC[ii][jj+1] = pidTPC;
      fPIDnSigmaProbTPC[ii][jj+1] = probTPC[jj];

      pidStat = fPIDResponse->NumberOfSigmas(
        AliPIDResponse::kTOF,fTrkPair[ii],(AliPID::EParticleType)jj,pidTOF);
      fPIDnSigmaTOF[ii][jj+1] = pidTOF;
      fPIDnSigmaProbTOF[ii][jj+1] = probTOF[jj];

    }
      
    // Bayes with prior
    fPIDBayesWP[ii][0] = fPIDCombined1->ComputeProbabilities(fTrkPair[ii], fPIDResponse, kpid);

    ptot = 0.;
    for (Int_t kk=0; kk<AliPID::kSPECIES; kk++) ptot += kpid[kk];
    if (abs(ptot-1.) < 1.E-3) {
      for (Int_t kk=0; kk<AliPID::kSPECIES; kk++)
        fPIDBayesWP[ii][kk+1] = kpid[kk];
    } else {
      for (Int_t kk=0; kk<AliPID::kSPECIES; kk++)
        fPIDBayesWP[ii][kk+1] = -999.9;
    }

    // Bayes without prior
    fPIDBayesNP[ii][0] = fPIDCombined2->ComputeProbabilities(fTrkPair[ii], fPIDResponse, kpid);
        
    ptot = 0.;
    for (Int_t kk=0; kk<AliPID::kSPECIES; kk++) ptot += kpid[kk];
    if (abs(ptot-1.) < 1.E-3) {
      for (Int_t kk=0; kk<AliPID::kSPECIES; kk++)
        fPIDBayesNP[ii][kk+1] = kpid[kk];
    } else {
      for (Int_t kk=0; kk<AliPID::kSPECIES; kk++)
        fPIDBayesNP[ii][kk+1] = -999.9;
    }

    // combine the pid information of the two tracks
    // is using Bayes without priors
    // find largest probability for each track
    Int_t imax[2] = {0,0};
    for (Int_t ii=1; ii<AliPID::kSPECIES; ii++) {
      if (kpid[ii]>kpid[imax[0]]) imax[0] = ii;
      if (kpid[ii]>kpid[imax[1]]) imax[1] = ii;
    }
    //printf("pids %i (%e) / %i (%e)\n",
    //  imax[0],kpid[imax[0]],imax[1],kpid[imax[1]]);
      
    combPID = AliPBBase::kBinPIDUnknown;
    if (imax[0] == imax[1]) {
      if (imax[0] == 0) combPID = AliPBBase::kBinElectron;
      if (imax[0] == 1) combPID = AliPBBase::kBinPIDUnknown;
      if (imax[0] == 2) combPID = AliPBBase::kBinPion;
      if (imax[0] == 3) combPID = AliPBBase::kBinKaon;
      if (imax[0] == 4) combPID = AliPBBase::kBinProton;
    }
    
	}
	
  // fill fhStatsFlow with PID information
  if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitStatsFlow)) {
		if (combPID == AliPBBase::kBinPion ||
		    combPID == AliPBBase::kBinPionE) {
			fhStatsFlow->Fill(AliPBBase::kBinPionEvents);
		}
		else if (combPID == AliPBBase::kBinKaon ||
				         combPID == AliPBBase::kBinKaonE) {
			fhStatsFlow->Fill(AliPBBase::kBinKaonEvents);
		}
		else if (combPID == AliPBBase::kBinProton ||
		         combPID == AliPBBase::kBinProtonE) {
			fhStatsFlow->Fill(AliPBBase::kBinProtonEvents);
		}
		else if (combPID == AliPBBase::kBinElectron ||
		         combPID == AliPBBase::kBinElectronE) {
			fhStatsFlow->Fill(AliPBBase::kBinElectronEvents);
		}
		else {
			fhStatsFlow->Fill(AliPBBase::kBinUnknownPIDEvents);
		}
	}

	return combPID;
}

//------------------------------------------------------------------------------
void AliAnalysisTaskPB::DetermineMCprocessType()
{
	//
	// retrieves the MC process type from the AliGenEventHeader and classifies
	// them
	//

	// get MC information
	fMCprocess = -1; //detailed MC sub process information
	fMCprocessType = AliPBBase::kBinND; // ND is default, also for data

	if (fMCEvent) {
		AliGenEventHeader* header = fMCEvent->GenEventHeader();
		if (header) {
      printf("MC generator name: %s\n",TString(header->GetName()).Data());

			// Pythia6
			if (TString(header->IsA()->GetName()) == "AliGenPythiaEventHeader") {
				fMCprocess = ((AliGenPythiaEventHeader*)header)->ProcessType();
				switch(fMCprocess) {
				case 92:
				case 93:
				case 103:
				case 104: fMCprocessType = AliPBBase::kBinSD; break;
				case 94:
				case 105: fMCprocessType = AliPBBase::kBinDD; break;
				default:  fMCprocessType = AliPBBase::kBinND; break;
				}
			}
			// Phojet
			else if (TString(header->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
				fMCprocess = ((AliGenDPMjetEventHeader*)header)->ProcessType();
				switch(fMCprocess) {
				case 5:
				case 6:  fMCprocessType = AliPBBase::kBinSD; break;
				case 7:  fMCprocessType = AliPBBase::kBinDD; break;
				case 4:  fMCprocessType = AliPBBase::kBinCD; break;
				default: fMCprocessType = AliPBBase::kBinND; break;
				}
			}
		}
	}
}

//------------------------------------------------------------------------------
void AliAnalysisTaskPB::FillMChists(Int_t combCh)
{
	// fill the MC histograms for different gap conditions. Only real two track
	// events are taken into account
	//

	if (fMCprocess > 0) { // MC process
		for (Int_t iGap = kV0; iGap < kMax; ++iGap) {
			if (fGapInformation[iGap] == AliPBBase::kBinDG) {
				if (combCh == AliPBBase::kBinPM) {
					fMCProcessUls->Fill(fMCprocess, iGap);
				}
				else if (combCh == AliPBBase::kBinPPMM) {
					fMCProcessLs->Fill(fMCprocess, iGap);
				}
			}
		}
	}
}


//------------------------------------------------------------------------------
Int_t AliAnalysisTaskPB::DoMCTruth()
{
	//
	// analyses the MC truth and do a gap selection based on the eta values of
	// the primaries
	//

	if (!fMCEvent) return -1;
	AliStack* stack = fMCEvent->Stack();
	if (!stack) return -1;
  //printf("We got a MC stack ...\n");

	//= Multiplicity =============================================================
	// determine number of charged physical primaries on the stack
	Int_t nPhysicalPrimaries = 0;
	Int_t nPrimaries = stack->GetNprimary();
  //printf("There are %i/%i/%i particles\n",
  //  stack->GetNtrack(),nPrimaries,stack->GetNtransported());
	for (Int_t iTracks = 0; iTracks < nPrimaries; ++iTracks) {
    
    TParticle* part = stack->Particle(iTracks);
    //printf("PDG code of primary %i: %i\n",iTracks,part->GetPdgCode());
		if (stack->IsPhysicalPrimary(iTracks) && (part->GetPDG()->Charge() != 0.) &&
		    (part->Eta() < 0.9) && (part->Eta() > -0.9)) { // TODO add parameter
			++nPhysicalPrimaries;
		}
	}
  
  
  // select the two pions and compute invariant mass
	Int_t inds[2] = {0,0};
  for (Int_t iTracks = 0; iTracks < nPrimaries; ++iTracks) {
    
    TParticle* part = stack->Particle(iTracks);
    if (part->GetPdgCode() ==  211) inds[0] = iTracks;
    if (part->GetPdgCode() == -211) inds[1] = iTracks;

  }
  TLorentzVector v1, v2, vsum;
  stack->Particle(inds[0])->Momentum(v1);
  stack->Particle(inds[1])->Momentum(v2);
  vsum = v1+v2;
  //printf("MC particle mass: %f\n",vsum.Mag());

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
		Bool_t gapFromData = (fGapInformation[kV0] == AliPBBase::kBinDG);

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
		Bool_t gapFromData = (fGapInformation[kV0FMD] == AliPBBase::kBinDG);

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
		Bool_t gapFromData = (fGapInformation[kV0FMDSPD] == AliPBBase::kBinDG);

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
			(fGapInformation[kV0FMDSPDTPC] == AliPBBase::kBinDG);

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

	//printf("number of physical primaries: %i\n",nPhysicalPrimaries);
  if ((nPhysicalPrimaries < 2 ) || (nPhysicalPrimaries >= 4 ))
		return nPhysicalPrimaries;

	Int_t gapCond = AliPBBase::kBinNG;
	if ((gapA == 0) && (gapC == 0) && (central > 0))
		gapCond = AliPBBase::kBinDG;
	else if ((gapA == 0) && (gapC > 0)) gapCond = AliPBBase::kBinGA;
	else if ((gapA > 0) && (gapC == 0)) gapCond = AliPBBase::kBinGC;
	else if ((gapA > 0) && (gapC > 0)) gapCond = AliPBBase::kBinNG;

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
void AliAnalysisTaskPB::DoMCTrackPair(TParticle* particles[],
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
		AliPBBase::kBinPPMM : AliPBBase::kBinPM;

	// determine PID
	Int_t combPID = 0;
	if (gapCond == AliPBBase::kBinDG) {
		combPID = (combCh == AliPBBase::kBinPM) ? // assignment to uls/ls histo
			AliPBUtils::GetCombPID((const TParticle**)particles,
			                            fComb2trkPIDulsDGmc) :
			AliPBUtils::GetCombPID((const TParticle**)particles,
			                            fComb2trkPIDlsDGmc);
	}
	else {
		combPID = (combCh == AliPBBase::kBinPM) ? // assignment to uls/ls histo
			AliPBUtils::GetCombPID((const TParticle**)particles,
			                            fComb2trkPIDulsMC) :
			AliPBUtils::GetCombPID((const TParticle**)particles,
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
	
  //printf("\nCombined PID: %i\n",combPID);
  AliPBUtils::GetMassPtCtsOA(combPID, momenta, mass, pt, cts, oa);

	//= THnMother ================================================================
  if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitTHnMother)) {
    if ((gapCond != AliPBBase::kBinNG) ||
		    (rnd.Rndm() > (1. - fNoGapFraction)) || fDoAOD) {

      AliPBBase::FillThnMother(fThnMC, multiplicity, (Double_t)combCh,
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
void AliAnalysisTaskPB::AnalyzeVtx()
{
	// calculates the distance between the vertex obtain from tracks and the
	// vertex obtain from spd tracklets
	// stores the z position of the primary vertex from tracks

	fVtxDst = 0.;     // reset distance

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
  //printf("Vertex distance: %f\n",fVtxDst);

	if (!fAnalysisStatus || (fAnalysisStatus & AliPBBase::kBitVtxStudies)) {
		fPriVtxDst->Fill(fVtxDst);
		if (fGapInformation[kV0FMDSPDTPC] == AliPBBase::kBinDG) {
			fPriVtxDstV0FMDSPDTPCdg->Fill(fVtxDst);
			fPriVtxDstV0FMDSPDdg->Fill(fVtxDst);
			fPriVtxDstV0FMDdg->Fill(fVtxDst);
			fPriVtxDstV0dg->Fill(fVtxDst);
		}
		else if (fGapInformation[kV0FMDSPD] == AliPBBase::kBinDG) {
			fPriVtxDstV0FMDSPDdg->Fill(fVtxDst);
			fPriVtxDstV0FMDdg->Fill(fVtxDst);
			fPriVtxDstV0dg->Fill(fVtxDst);
		}
		else if (fGapInformation[kV0FMD] == AliPBBase::kBinDG) {
			fPriVtxDstV0FMDdg->Fill(fVtxDst);
			fPriVtxDstV0dg->Fill(fVtxDst);
		}
		else if (fGapInformation[kV0] == AliPBBase::kBinDG) {
			fPriVtxDstV0dg->Fill(fVtxDst);
		}
	}
}

//--------------------------------------------------------------------------
void AliAnalysisTaskPB::GetMyPriors(TString fnameMyPriors, TH1F** mypriors) {

  // open file
  TFile *priorff = TFile::Open(fnameMyPriors.Data(),"READ");
  if (priorff) {
    printf("File %s opened!\n",fnameMyPriors.Data());
  
    for (Int_t ii=0;ii<AliPID::kSPECIES;ii++) {
      mypriors[ii] = (TH1F*)priorff->Get(Form("pppriors%i",ii));
      mypriors[ii]->SetLineStyle(kSolid);
    }
  }
  
}

// -----------------------------------------------------------------------------
