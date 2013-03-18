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
// Example analysis for diffractive studies
//
// Author:
//  Felix Reidt <Felix.Reidt@cern.ch>


#include <TH1.h>
#include <TH2.h>
#include <TList.h>


#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliESDInputHandler.h"
#include "AliPIDResponse.h"
#include "AliPhysicsSelection.h"

#include "AliCDMesonBaseStripped.h"
#include "AliCDMesonUtilsStripped.h"
#include "AliCDMesonTracks.h"
#include "AliAnalysisTaskCDex.h"


//------------------------------------------------------------------------------
AliAnalysisTaskCDex::AliAnalysisTaskCDex(const char* name):
	AliAnalysisTaskSE(name)
	, fDoAOD(kFALSE)
	, fMaxVtxDst(0.5) // value to be checked with the vertex study histograms

	, fESDEvent(0x0)
	, fAODEvent(0x0)
	, fPIDResponse(0x0)
	, fTracks(new AliCDMesonTracks())
	, fVtxDst(-1.)
	, fVtxZ(-20)
	, fResidualTracks(0)
	, fResidualTracklets(0)
	, fMCprocessType(0)
	, fMCprocess(-1)

	, fRun(-999)
	, fCurrentGapCondition(0)

	, fHist(0x0)

	, fv0ntrk(0x0)
	, fv0fmdntrk(0x0)
	, fv0fmdspdntrk(0x0)
	, fv0fmdspdtpcntrk(0x0)

	, fhStatsFlow(0x0)
{
	//
	// standard constructor
	//
	// slot in TaskSE must start from 1
	DefineOutput(1, TList::Class());

	// initialize gap information
	for (Int_t iGap = 0; iGap < kMax; ++iGap) {
		fGapInformation[iGap] = 0;
	}
}


//------------------------------------------------------------------------------
AliAnalysisTaskCDex::~AliAnalysisTaskCDex()
{
	//
	// destructor
	//

	if (fHist) {
		fHist->Clear();
		delete fHist;
		fHist = 0x0;
	}

	if (fTracks) {
		delete fTracks;
		fTracks = 0x0;
	}
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDex::UserCreateOutputObjects()
{
	//
	// createOutputObjects
	//

	//= TList for Histograms =====================================================
	fHist = new TList;
	fHist->SetOwner(); // ensures that the histograms are all deleted on exit!

	//= MULTIPLICITY PER GAP CONDITION =
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


	//= STATISTICS FLOW =
	fhStatsFlow = AliCDMesonBaseStripped::GetHistStatsFlow();
	fHist->Add(fhStatsFlow);
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDex::UserExec(Option_t *)
{
	//
	// executed for every event which passed the physics selection
	//
	// in order to select only correct minimum bias events,
	// SetCollisionCandidates(AliVEvent::kMB) should be used

	fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinTotalInput); // stats flow

	//= INPUT DATA SANITY TESTS ==================================================
	if(!CheckInput()) {
		PostOutputs();
		return;
	}

	fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinGoodInput); // stats flow

	//= EVENT SELECTION ==========================================================
	Bool_t eventIsValid = (fDoAOD) ?
		AliCDMesonUtilsStripped::CutEvent(fAODEvent) :
		AliCDMesonUtilsStripped::CutEvent(fESDEvent);
	if (!eventIsValid) {
		PostOutputs();
		return;
	}

	fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinEventsAfterCuts); // stats flow

	//= PILE UP ==================================================================
	const Bool_t isPileup = (fDoAOD) ?
		fAODEvent->IsPileupFromSPD(2, 0.8, 3., 2., 5.) :
		fESDEvent->IsPileupFromSPD(2, 0.8, 3., 2., 5.);
	// using only 2 instead of three contributors

	if (isPileup) {
		PostOutputs();
		return;
	}

	fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinEventsWithOutPileUp);

	//= GAP ======================================================================
	// determine the complete gap configuration (for all gap-tagging detectors)
	if (!DetermineGap()) {
		PostOutputs();
		return;
	}

	//= VERTEX COINCIDENCE AND POSITION ==========================================
	AnalyzeVtx();
	if (!(abs(fVtxZ) < 4.)) { // vertex from tracks within +/-4cm
		//PostOutputs();
		//return;
	}

	//= TRACK CUTS ===============================================================
	fTracks->ProcessEvent(fAODEvent, fESDEvent, kTRUE);
	// apply cuts (including soft)
	DoMultiplicityStudy(); // fill corresponding histograms

	// is multiplicity within the desired range of  2 to 3?
	Int_t nch = fTracks->GetTracks();
	Int_t ncombined = fTracks->GetCombinedTracks();

	//============================================================================
	//=== USER ANALYSIS CODE =====================================================
	//============================================================================
	for (Int_t iTrack = 0; iTrack < ncombined; ++iTrack) { // including soft
		AliVTrack* trk = fTracks->GetTrack(iTrack);
		trk->GetID(); // prevent warning...
	}
	for (Int_t iTrack = 0; iTrack < nch; ++iTrack) { // excluding soft tracks
		AliVTrack* trk = fTracks->GetTrack(iTrack);
		trk->GetID(); // prevent warning...
	}

	if (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBaseStripped::kBinDG) {
		// event is a full double gap event
	}

	//============================================================================
	PostOutputs();
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDex::PostOutputs()
{
	//
	// PostData
	//
	// this function is main of use with multiple output containers

	PostData(1, fHist);
}


//------------------------------------------------------------------------------
Bool_t AliAnalysisTaskCDex::CheckInput()
{
	//
	// general protection of the task against malicious input data
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
	fPIDResponse = (AliPIDResponse*)fInputHandler->GetPIDResponse();

	if(!fESDEvent && !fAODEvent){
		printf("AliAnalysisTaskex - No valid event\n");
		return kFALSE;
	}

	if(!fPIDResponse){
		printf("AliAnalysisTaskex -  No PIDd\n");
		// PID is fixed to unknown
		//return kFALSE;
	}

	if(fDoAOD && fAODEvent && fabs(fAODEvent->GetMagneticField())<1){
		printf("AliAnalysisTaskex - strange Bfield! %f\n",
		       fAODEvent->GetMagneticField());
		return kFALSE;
	}
	else if((!fDoAOD) && fESDEvent && fabs(fESDEvent->GetMagneticField())<1){
		printf("AliAnalysisTaskex - strange Bfield! %f\n",
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

	if(fRun!=tmprun){
		fRun = tmprun;
		AliCDMesonUtilsStripped::SPDLoadGeom(fRun);
	}

	return kTRUE;
}


//------------------------------------------------------------------------------
Bool_t AliAnalysisTaskCDex::DetermineGap()
{
	// determines the gap configuration for all gap tagging detectors based on the
	// data set which is available
	//

	if (fDoAOD) {
		AliAODHandler* aodHandler =
			(AliAODHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
		TTree *aodTree = aodHandler->GetTree();
		aodTree->SetBranchAddress("gapCondition", &fCurrentGapCondition);
		aodTree->GetEvent(Entry()); // seems to be needed! (loads current event)
		if (!fCurrentGapCondition) {
			fCurrentGapCondition = 0xfffe;
			puts("AliAnalysisTaskCDex - ");
			puts("error while gap condition determination using AODs\n");
			return kFALSE;
		}
	}
	else {
		// gap determination from detector information
		fCurrentGapCondition = AliCDMesonUtilsStripped::GetGapConfig(fESDEvent);

		// gap determination from preprocessed detector information
		/*
		AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
		TTree* esdTree = am->GetInputEventHandler()->GetTree(); // get ESD tree
		if (esdTree) {
			esdTree->SetBranchAddress("gapCondition", &fCurrentGapCondition);
			esdTree->GetEvent(Entry()); // seems to be needed! (loads current event)
		}
		*/

		if (!fCurrentGapCondition) {
			fCurrentGapCondition = 0xfffe;
			puts("AliAnalysisTaskCDex - ");
			puts("error while gap condition determination using ESDs\n");
			return kFALSE;
		}
	}

	// disentagle the contributions to the gap conditions of different "tightness"
	fGapInformation[kV0] =
		AliCDMesonBaseStripped::GetGapBin("V0", fCurrentGapCondition);
	fGapInformation[kV0FMD] =
		AliCDMesonBaseStripped::GetGapBin("V0FMD", fCurrentGapCondition);
	fGapInformation[kV0FMDSPD] =
		AliCDMesonBaseStripped::GetGapBin("V0FMDSPD", fCurrentGapCondition);
	fGapInformation[kV0FMDSPDTPC] =
		AliCDMesonBaseStripped::GetGapBin("V0FMDSPDTPC", fCurrentGapCondition);
	fGapInformation[kFMD] =
		AliCDMesonBaseStripped::GetGapBin("FMD",fCurrentGapCondition);
	fGapInformation[kSPD] =
		AliCDMesonBaseStripped::GetGapBin("SPD",fCurrentGapCondition);
	fGapInformation[kTPC] =
		AliCDMesonBaseStripped::GetGapBin("TPC",fCurrentGapCondition);

	return kTRUE;
}


//------------------------------------------------------------------------------
void AliAnalysisTaskCDex::DoMultiplicityStudy()
{
	// stores the multiplicity distributions for different gap conditions and
	// adds some information to the statsFlow histogram
	//

	// retrieve values from the track object
	Int_t ntrk0 = fTracks->GetTracksBeforeCuts(); // number of tracks before cuts
	//Int_t nch = fTracks->GetTracks(); // number of good ITS-TPC primaries
	Int_t ncombined = fTracks->GetCombinedTracks(); // number ITSTPC and ITS only
	Int_t nITSpureSA = fTracks->GetITSpureSACount(); // number ITS standalone

	// determine the residual tracks / tracklets
	fResidualTracks = ntrk0 - ncombined - nITSpureSA;
	fResidualTracklets = fTracks->GetRemainingTrackletsCentralBarrel();

	// multiplicity distributions for different gaps
	fv0ntrk->Fill(ncombined, fGapInformation[kV0]);
	fv0fmdntrk->Fill(ncombined, fGapInformation[kV0FMD]);
	fv0fmdspdntrk->Fill(ncombined, fGapInformation[kV0FMDSPD]);
	fv0fmdspdtpcntrk->Fill(ncombined, fGapInformation[kV0FMDSPDTPC]);

	if (fGapInformation[kV0] == AliCDMesonBaseStripped::kBinDG) {
		fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinv0Gap);
	}
	if (fGapInformation[kV0FMD] == AliCDMesonBaseStripped::kBinDG) {
		fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinv0fmdGap);
	}
	if (fGapInformation[kV0FMDSPD] == AliCDMesonBaseStripped::kBinDG) {
		fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinv0fmdspdGap);
	}
	if (fGapInformation[kV0FMDSPDTPC] == AliCDMesonBaseStripped::kBinDG) {
		fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinv0fmdspdtpcGap);
	}

	// event cleanliness
	if (fResidualTracks == 0) {
		fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinResidualTracks);
	}
	if (fResidualTracklets == 0) {
		fhStatsFlow->Fill(AliCDMesonBaseStripped::kBinResidualTracklets);
	}
}


//--------------------------------------------------------------------------
void AliAnalysisTaskCDex::AnalyzeVtx()
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
}
