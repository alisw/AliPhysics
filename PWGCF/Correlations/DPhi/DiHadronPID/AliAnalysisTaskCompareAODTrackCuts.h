#ifndef ALIANALYSISTASKCOMPAREAODTRACKCUTS_H
#define ALIANALYSISTASKCOMPAREAODTRACKCUTS_H

#include "AliAODEventCutsDiHadronPID.h"
#include "AliAODTrackCutsDiHadronPID.h"
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
* See cxx source for full Copyright notice */ 
/* $Id$ */

class AliAnalysisTaskCompareAODTrackCuts : public AliAnalysisTaskSE {

// Basic AnalysisTask Functions.
public:
	AliAnalysisTaskCompareAODTrackCuts();
	AliAnalysisTaskCompareAODTrackCuts(const char* name);
	virtual ~AliAnalysisTaskCompareAODTrackCuts();
	
private:
	AliAnalysisTaskCompareAODTrackCuts(const AliAnalysisTaskCompareAODTrackCuts&);
	AliAnalysisTaskCompareAODTrackCuts& operator=(const AliAnalysisTaskCompareAODTrackCuts&);		

public:
	virtual void UserCreateOutputObjects(); 
	virtual void UserExec(Option_t*);
	virtual void Terminate(Option_t*);	

// Mismatch related functions.
	Bool_t LoadExternalMismatchHistos();			// For each mismatch method, external histo's are needed.
	Double_t GenerateRandomHit(Double_t eta);		// Generates a random time for a certain eta.

// Task Settings.
	void SetMC(Bool_t isMC = kTRUE) {fIsMC = isMC;} 
	void SetVerbose(Bool_t verbose = kTRUE) {fVerbose = verbose;}
	void SetCalculateTOFMismatch(Bool_t calculatetofmismatch = kTRUE/* const Int_t method*/) {fCalculateTOFMismatch = calculatetofmismatch;}
	void SetUseMismatchFileFromGridHomeDir(Bool_t usefromhomedir = kTRUE) {fUseMismatchFileFromHomeDir = usefromhomedir;}
	void SetUseNSigmaOnPIDAxes(Bool_t UseNSigma = kTRUE);

// Managing Cuts.
    void SetEventCuts(AliAODEventCutsDiHadronPID* eventcuts);
    void AddTrackCuts(AliAODTrackCutsDiHadronPID* trackcuts);

// Override from AliAnalysisTaskSE.
	void SetDebugLevel(Int_t debuglvl);

private:
	void FillGlobalTracksArray();
	AliAODTrack* GetGlobalTrack(AliAODTrack* track);

private:

	// PID Response Object.
	AliPIDResponse*					fPIDResponse;				//! PID Response.

	// Output List.
	TList* 							fOutputList;				//! Output List.

	// Settings (streamed!).
	Bool_t							fIsMC;						// ran over MC or not.
	Bool_t							fVerbose;					// Verbose mode.
	Bool_t 							fCalculateTOFMismatch;		// Compute mismatch or not. (Needs input histograms!)
	Bool_t							fUseMismatchFileFromHomeDir;// Take TOFmistmachHistos.root from the home dir, or take the one copied when submitting jobs
	Bool_t							fUseNSigmaOnPIDAxes;		// Uses NSigma on the PID axes of all histograms.

	// Event Cut Object (streamed!).
	AliAODEventCutsDiHadronPID*		fEventCuts;					// Event Cuts.
	
	// Array of Track Cut Objects (streamed!).
	TObjArray*						fTrackCuts;					// TObjArray with all Track Cut Objects.

	// Inclusive track times.
	TH2F* 							fInclusiveTimes;			//!

	// TOF mismatch stuff.
	TH1F*							fT0Fill;					//
	TH2F*							fLvsEta;					//
	TH2F*							fGlobalPtvsTPCPt;			//
	TObjArray*						fLvsEtaProjections;			//	


	// Event and Track related objects.
	AliAODEvent*					fCurrentAODEvent;			//!
	AliAODTrack*					fCurrentAODTrack;			//!
	AliTrackDiHadronPID*			fCurrentDiHadronPIDTrack;	//!
	TObjArray*						fGlobalTracksArray;			//! Array of Global Tracks.

	ClassDef(AliAnalysisTaskCompareAODTrackCuts,1);

};

#endif
