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

#ifndef ALIANALYSISTASKCDEX_H
#define ALIANALYSISTASKCDEX_H

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

class AliESDEvent;
class AliVTrack;
class AliPIDResponse;
class AliPhysicsSelection;

class TH1I;
class TH1F;
class TH2I;
class TH2F;
class TH2D;
class TList;
class THnSparse;
class TObjString;

class AliCDMesonTracks;

class AliAnalysisTaskCDex : public AliAnalysisTaskSE
{
public:

	AliAnalysisTaskCDex(const char* name);
	virtual ~AliAnalysisTaskCDex();

	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t *);

private:
	enum { //
		kV0 = 0,
		kFMD,
		kSPD,
		kTPC,
		kV0FMD,
		kV0FMDSPD,
		kV0FMDSPDTPC,
		kMax
	};

	AliAnalysisTaskCDex(const AliAnalysisTaskCDex  &p);
	AliAnalysisTaskCDex& operator=(const AliAnalysisTaskCDex  &p);

	void FillEtaPhiMaps() const; // controls which event contributes to which map


	// functions called by the UserExec(...), not to be called elsewhere!
	//-------------------------------------------------------------------
	Bool_t CheckInput();
	void PostOutputs(); // cares about posting the output before exiting UserExec,
	// WARNING: PostOutputs should only be used directly in the UserExec!!
	Bool_t DetermineGap(); // determines the gap of all available detectors
	void AnalyzeVtx(); // calcs the distance of the pri. vertex from tracks and SPD
	void DoMultiplicityStudy(); // multiplicity distributions for different gaps

	// analysis task status
	//---------------------
	Bool_t fDoAOD; // true for running on AODs
	Double_t fMaxVtxDst; // maximum distance of the track and SPD vertex

	// event information
	//------------------
	AliESDEvent *fESDEvent; // esd event object
	AliAODEvent *fAODEvent; // esd event object
	AliPIDResponse *fPIDResponse; // pid object (for ESDs and AODs)
	AliCDMesonTracks *fTracks; // object taking care about the track cuts
	Double_t fVtxDst; // distance of the primary vertex from tracks and from SPD
	Double_t fVtxZ; // z-position of the primary vertex from tracks
	Int_t fResidualTracks; // tracks rejected by cuts within the event
	Int_t fResidualTracklets; // SPD tracklets not assigned to tracks
	Int_t fMCprocessType; // MC process type, 0 for data
	Int_t fMCprocess; // detailed MC sub process information

	Int_t fRun; // number of the run which is about to be processed
	Int_t fCurrentGapCondition; // gap condition of the current event
	Int_t fGapInformation[kMax]; // gap condition for different detectors
	                             // individually and their combinations

	// output objects
	//---------------

	TList *fHist; // output list (contains all histograms)

	// Multiplicity distributions for the different gap conditions
	TH2D *fv0ntrk; //v0bit vs. nch
	TH2D *fv0fmdntrk; //v0fmdbit vs. nch
	TH2D *fv0fmdspdntrk; //v0fmdspdbit vs. nch
	TH2D *fv0fmdspdtpcntrk; //v0fmdspdtpcbit vs. nch

	// Statistics flow diagrams
	TH1F *fhStatsFlow; // stepwise statistics flow

	ClassDef(AliAnalysisTaskCDex, 1);
};


#endif
